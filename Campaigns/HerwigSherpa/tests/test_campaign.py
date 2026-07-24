#!/usr/bin/env python3

import argparse
from concurrent.futures import Future, ThreadPoolExecutor
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import threading
import time
import types
import unittest
from unittest import mock


CAMPAIGN_DIR = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "run_campaign", CAMPAIGN_DIR / "run_campaign.py"
)
CAMPAIGN = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(CAMPAIGN)


class FakeDbn:

    def __init__(self, sum_w, entries, value=0.0):
        self._sum_w = sum_w
        self._entries = entries
        self._value = value

    def sumW(self):
        return self._sum_w

    def numEntries(self):
        return self._entries

    def val(self):
        return self._value


class FakeHistogram:

    def __init__(self, bins, underflow, overflow):
        self._bins = bins
        self._flow_bins = [underflow, *bins, overflow]

    def bins(self):
        return self._bins

    def bin(self, index):
        return self._flow_bins[index]


class FakeEstimate:

    def __init__(self, edges):
        self.xEdges = edges
        self._bins = [FakeDbn(0.0, 0.0, 0.0) for _ in edges[:-1]]

    def bins(self):
        return self._bins


class CampaignConfigurationTest(unittest.TestCase):

    def setUp(self):
        self.config = CAMPAIGN.load_config()

    def test_static_validation(self):
        CAMPAIGN.static_checks(self.config)

    def test_sample_matrix_and_event_budget(self):
        self.assertEqual(
            CAMPAIGN.profile_names(self.config),
            ["gg-low", "gg-high", "all-low", "all-high"],
        )
        self.assertEqual(
            CAMPAIGN.sample_ids(self.config),
            [
                "herwig",
                "herwig-qcdinsplanar",
                "herwig-random3-dipole",
                "sherpa",
            ],
        )
        total = (
            sum(region["events"] for region in self.config["regions"].values())
            * len(self.config["processes"])
            * len(self.config["samples"])
        )
        self.assertEqual(total, 560000)

    def test_sample_physics_settings(self):
        baseline = CAMPAIGN.sample_definition(self.config, "herwig")
        planar = CAMPAIGN.sample_definition(
            self.config, "herwig-qcdinsplanar"
        )
        dipole = CAMPAIGN.sample_definition(
            self.config, "herwig-random3-dipole"
        )
        self.assertEqual(
            (baseline["colour"], baseline["shower"]),
            ("Random3", "angular"),
        )
        self.assertEqual(
            (planar["colour"], planar["shower"]),
            ("QCDINSPlanar", "angular"),
        )
        self.assertEqual(
            (dipole["colour"], dipole["shower"]),
            ("Random3", "dipole"),
        )

    def test_sherpa_process_families(self):
        self.assertEqual(CAMPAIGN.sherpa_processes("gg"), [(21, 21)])
        all_processes = CAMPAIGN.sherpa_processes("all")
        self.assertEqual(len(all_processes), 56)
        self.assertEqual(len(set(all_processes)), 56)
        self.assertNotIn((1, 1), all_processes)
        self.assertNotIn((-1, -1), all_processes)
        self.assertIn((1, -1), all_processes)
        self.assertIn((21, 5), all_processes)
        self.assertIn((21, -5), all_processes)

    def test_sherpa_native_thresholds(self):
        sherpa = self.config["sherpa"]
        self.assertEqual(sherpa["charm_threshold_gev"], 20.0)
        self.assertEqual(sherpa["bottom_threshold_gev"], 100.0)
        self.assertEqual(sherpa["include_quarks"], 5)

    def test_sherpa_incoming_threshold_patch_is_tracked(self):
        patch = (
            CAMPAIGN_DIR
            / "patches/sherpa-instanton-incoming-thresholds.patch"
        ).read_text(encoding="utf-8")
        self.assertIn("IncomingFlavoursActive() const", patch)
        self.assertIn("code==kf_b && m_Ehat<m_bthreshold", patch)
        self.assertIn("code==kf_c && m_Ehat<m_cthreshold", patch)
        self.assertIn("cols[0].size()!=cols[1].size()", patch)

    def test_sherpa_card_uses_matched_pdf_slots(self):
        card = CAMPAIGN.render_sherpa_profile(self.config, "gg", "low")
        self.assertIn("PDF_SET: [NNPDF31_nnlo_as_0118]", card)
        self.assertIn("MPI_PDF_SET: [NNPDF31_nnlo_as_0118]", card)
        self.assertIn("MI_HANDLER: None", card)
        self.assertNotIn("ME_SIGNAL_GENERATOR", card)

    def test_generated_cards_are_current(self):
        cards = CAMPAIGN.rendered_cards(self.config)
        self.assertEqual(len(cards), 16)
        for path, content in cards.items():
            self.assertTrue(path.is_file(), path)
            self.assertEqual(path.read_text(encoding="utf-8"), content)

    def test_baseline_card_remains_the_legacy_random3_card(self):
        sample = CAMPAIGN.sample_definition(self.config, "herwig")
        card = CAMPAIGN.render_herwig_profile(
            self.config, sample, "gg", "high"
        )
        self.assertIn(
            "read Campaigns/HerwigSherpa/cards/herwig/common.in", card
        )
        self.assertIn(
            "set /Herwig/Analysis/Rivet:Filename herwig-gg-high.yoda", card
        )
        self.assertIn(
            "saverun Campaign-Herwig-gg-high EventGenerator", card
        )
        self.assertNotIn("DipoleShowerHandler", card)
        self.assertNotIn("QCDINSPlanar", card)

    def test_qcdinsplanar_cards_only_override_colour(self):
        sample = CAMPAIGN.sample_definition(
            self.config, "herwig-qcdinsplanar"
        )
        for process in self.config["processes"]:
            for region in self.config["regions"]:
                card = CAMPAIGN.render_herwig_profile(
                    self.config, sample, process, region
                )
                self.assertIn(
                    "MEInstanton:ColourConnections QCDINSPlanar", card
                )
                self.assertNotIn("Dipole_AutoTunes_gss", card)
                self.assertNotIn("DipoleShowerHandler", card)

    def test_dipole_cards_restore_matched_settings_after_tune(self):
        sample = CAMPAIGN.sample_definition(
            self.config, "herwig-random3-dipole"
        )
        card = CAMPAIGN.render_herwig_profile(
            self.config, sample, "all", "high"
        )
        tune = card.index("read snippets/Dipole_AutoTunes_gss.in")
        restored_masses = card.index(
            "set /Herwig/Particles/b:HardProcessMass 4.18*GeV"
        )
        self.assertLess(tune, restored_masses)
        self.assertIn(
            "EventHandler:CascadeHandler "
            "/Herwig/DipoleShower/DipoleShowerHandler",
            card,
        )
        self.assertIn(
            "DipoleShowerHandler:PDFA /Herwig/Partons/InstantonPDF", card
        )
        self.assertIn(
            "DipoleShowerHandler:PDFB /Herwig/Partons/InstantonPDF", card
        )
        self.assertIn("DipoleShowerHandler:MPIHandler NULL", card)
        self.assertIn(
            "#set /Herwig/DipoleShower/DipoleShowerHandler:MPIHandler "
            "/Herwig/UnderlyingEvent/MPIHandler",
            card,
        )
        masses = {
            "u": "0.0023",
            "d": "0.005",
            "s": "0.093",
            "c": "1.27",
            "b": "4.18",
        }
        for flavour, mass in masses.items():
            self.assertIn(
                f"set /Herwig/Particles/{flavour}:NominalMass 0*GeV",
                card,
            )
            self.assertIn(
                f"set /Herwig/Particles/{flavour}bar:NominalMass 0*GeV",
                card,
            )
            self.assertIn(
                f"set /Herwig/Particles/{flavour}:HardProcessMass "
                f"{mass}*GeV",
                card,
            )
            self.assertIn(
                f"set /Herwig/Particles/{flavour}bar:HardProcessMass "
                f"{mass}*GeV",
                card,
            )
        self.assertNotIn("OffShellInShower", card)
        for stock_name in (
            "FIMdx2dgxDipoleKernel",
            "FIMux2ugxDipoleKernel",
            "FIMcx2cgxDipoleKernel",
            "FIMsx2sgxDipoleKernel",
            "FIMbx2bgxDipoleKernel",
            "FFMdx2dgxDipoleKernel",
            "FFMux2ugxDipoleKernel",
            "FFMcx2cgxDipoleKernel",
            "FFMsx2sgxDipoleKernel",
            "FFMbx2bgxDipoleKernel",
        ):
            self.assertIn(f"{stock_name}:UseKernel No", card)
        for adapter, kinematics in (
            ("InstantonFFMqx2qgxDipoleKernel", "FFMassiveKinematics"),
            ("InstantonFIMqx2qgxDipoleKernel", "FIMassiveKinematics"),
        ):
            self.assertIn(f"create Herwig::{adapter}", card)
            self.assertIn(
                f"{adapter}:PDFRatio "
                "/Herwig/DipoleShower/Kernels/PDFRatio",
                card,
            )
            self.assertIn(
                f"{adapter}:SplittingKinematics "
                f"/Herwig/DipoleShower/Kinematics/{kinematics}",
                card,
            )
            self.assertIn(f"{adapter}:CMWScheme Factor", card)
            self.assertIn(
                "DipoleShowerHandler:Kernels 0 "
                f"/Herwig/DipoleShower/Kernels/{adapter}",
                card,
            )
        for obsolete_adapter in (
            "InstantonIFqx2qgxDipoleKernel",
            "InstantonIFMqx2qgxDipoleKernel",
            "InstantonFIqx2qgxDipoleKernel",
            "InstantonIIqx2qgxDipoleKernel",
        ):
            self.assertNotIn(obsolete_adapter, card)
        self.assertEqual(card.count(":CMWScheme Factor"), 2)
        self.assertIn("MEInstanton:ColourConnections Random3", card)

    def test_dipole_adapter_preserves_actual_final_state_emitters(self):
        source = (CAMPAIGN.ROOT / "InstantonDipoleKernels.cc").read_text(
            encoding="utf-8"
        )
        self.assertIn("isFinalLightQuark", source)
        self.assertIn(
            "!(index.emitterData()->mass() == ZERO\n"
            "           && index.spectatorData()->mass() == ZERO)",
            source,
        )
        self.assertIn("index.emitterData()->mass() > ZERO", source)
        self.assertEqual(source.count("return index.emitterData();"), 2)
        self.assertNotIn("flavour()", source)

    def test_dipole_mass_split_is_persistent_in_me(self):
        source = (CAMPAIGN.ROOT / "MEInstanton.cc").read_text(
            encoding="utf-8"
        )
        self.assertIn("canonical->hardProcessMass()", source)
        self.assertIn("setupHardOutgoingQuarks();", source)
        self.assertIn("HardOutgoingQuarks", source)
        self.assertIn("hideFromParticleLookup", source)
        self.assertIn(
            'describeHerwigMEInstanton("Herwig::MEInstanton", '
            '"Instantons.so", 3);',
            source,
        )

    def test_plotting_environment_ignores_missing_texmf_override(self):
        with tempfile.TemporaryDirectory() as temporary:
            work_dir = Path(temporary)
            missing = work_dir / "missing-texmf"
            with (
                mock.patch.object(CAMPAIGN, "WORK_DIR", work_dir),
                mock.patch.dict(
                    CAMPAIGN.os.environ,
                    {
                        "TEXMFCNF": str(missing),
                        "TEXMFHOME": str(missing),
                        "TEXINPUTS": str(missing),
                    },
                    clear=False,
                ),
            ):
                env = CAMPAIGN.plotting_environment()

            self.assertEqual(
                env["MPLCONFIGDIR"], str(work_dir / "matplotlib")
            )
            self.assertEqual(
                env["RIVET_PLOT_PATH"], str(CAMPAIGN.ROOT / "Rivet")
            )
            self.assertEqual(
                env["RIVET_INFO_PATH"], str(CAMPAIGN.ROOT / "Rivet")
            )
            self.assertNotIn("TEXMFCNF", env)
            self.assertNotIn("TEXMFHOME", env)
            self.assertNotIn("TEXINPUTS", env)

    def test_process_dependent_gluon_caps(self):
        cap = self.config["max_final_partons"]
        gluon_cap = self.config["gluon_cap"]
        expected = {
            (4, 2): 21,
            (4, 1): 22,
            (4, 0): 23,
            (5, 2): 19,
            (5, 1): 20,
            (5, 0): 21,
        }
        for key, maximum in expected.items():
            flavours, incoming_gluons = key
            outgoing_quarks = 2 * flavours - (2 - incoming_gluons)
            self.assertEqual(
                min(gluon_cap, cap - outgoing_quarks),
                maximum,
            )

    def test_shard_seeds_are_unique_and_legacy_seeds_are_stable(self):
        seeds = {
            CAMPAIGN.stable_seed(sample, profile, shard)
            for sample in CAMPAIGN.sample_ids(self.config)
            for profile in CAMPAIGN.profile_names(self.config)
            for shard in range(self.config["shards"])
        }
        self.assertEqual(len(seeds), 160)
        self.assertEqual(
            CAMPAIGN.stable_seed("herwig", "gg-high", 0), 668213947
        )
        self.assertEqual(
            CAMPAIGN.stable_seed("sherpa", "gg-high", 0), 329223182
        )

    def test_mass_and_migration_edges(self):
        edges = CAMPAIGN.MASS_EDGES
        self.assertEqual(len(edges), 125)
        self.assertEqual(edges[:81], tuple(float(x) for x in range(0, 801, 10)))
        self.assertEqual(
            edges[81:], tuple(float(x) for x in range(850, 3001, 50))
        )
        self.assertEqual(
            CAMPAIGN.MIGRATION_RATIO_EDGES,
            tuple(index / 10.0 for index in range(101)),
        )

    def test_sample_selection_preserves_requested_order(self):
        selected = CAMPAIGN.select_samples(
            self.config, "sherpa,herwig-qcdinsplanar"
        )
        self.assertEqual(
            [sample["id"] for sample in selected],
            ["sherpa", "herwig-qcdinsplanar"],
        )
        with self.assertRaises(RuntimeError):
            CAMPAIGN.select_samples(self.config, "herwig,herwig")
        with self.assertRaises(RuntimeError):
            CAMPAIGN.select_samples(self.config, "unknown")

    def test_cross_section_parsers(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            herwig = directory / "herwig.log"
            herwig.write_text(
                "estimated total cross section is "
                "( 1.25e-3 +/- 2.5e-5 ) nb\n",
                encoding="utf-8",
            )
            self.assertEqual(
                CAMPAIGN.parse_herwig_cross_section(herwig),
                (1.25, 0.025),
            )

            sherpa = directory / "sherpa.log"
            sherpa.write_text(
                "\x1b[1m2_1__G__G__Instanton\x1b[0m : "
                "2.0 pb +- ( 0.3 pb = 15 % )\n"
                "2_1__G__u__Instanton : "
                "3.0 pb +- ( 0.4 pb = 13 % )\n",
                encoding="utf-8",
            )
            self.assertEqual(
                CAMPAIGN.parse_sherpa_cross_section(sherpa),
                (5.0, 0.5),
            )

    def test_herwig_cross_section_compatibility(self):
        cross_sections = {
            "herwig": {
                "gg-low": {"value": 10.0, "uncertainty": 1.0},
            },
            "herwig-qcdinsplanar": {
                "gg-low": {"value": 11.0, "uncertainty": 1.0},
            },
            "herwig-random3-dipole": {
                "gg-low": {"value": 20.0, "uncertainty": 1.0},
            },
        }
        result = CAMPAIGN.herwig_cross_section_comparison(
            self.config, cross_sections
        )
        self.assertTrue(
            result["herwig-qcdinsplanar"]["gg-low"]["compatible_3sigma"]
        )
        self.assertFalse(
            result["herwig-random3-dipole"]["gg-low"]["compatible_3sigma"]
        )

    def test_herwig_integration_is_cached_per_sample(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            campaign_dir = root / "campaign"
            build_dir = root / "build"
            work_dir = root / "work"
            plugin_dir = build_dir / "herwig"
            rivet_dir = build_dir / "rivet"
            (campaign_dir / "cards/herwig").mkdir(parents=True)
            plugin_dir.mkdir(parents=True)
            rivet_dir.mkdir(parents=True)

            card = campaign_dir / "cards/herwig/gg-low.in"
            common = campaign_dir / "cards/herwig/common.in"
            herwig = root / "Herwig"
            card.write_text("card\n", encoding="utf-8")
            common.write_text("common\n", encoding="utf-8")
            herwig.write_text("binary\n", encoding="utf-8")
            (plugin_dir / "CampaignInstantons.so").write_text(
                "plugin\n", encoding="utf-8"
            )
            (rivet_dir / "RivetQCD_INSTANTON_KKS.so").write_text(
                "rivet\n", encoding="utf-8"
            )

            calls = []

            def fake_run(command, *, cwd=None, env=None, log_path=None):
                calls.append(command)
                (cwd / "Campaign-Herwig-gg-low.run").write_text(
                    "run\n", encoding="utf-8"
                )

            with (
                mock.patch.object(CAMPAIGN, "CAMPAIGN_DIR", campaign_dir),
                mock.patch.object(CAMPAIGN, "BUILD_DIR", build_dir),
                mock.patch.object(CAMPAIGN, "WORK_DIR", work_dir),
                mock.patch.object(
                    CAMPAIGN, "capture", return_value="Herwig 7.3.0"
                ),
                mock.patch.object(CAMPAIGN, "run_checked", fake_run),
            ):
                first = CAMPAIGN.prepare_herwig_run(
                    "herwig", "gg-low", card, herwig, plugin_dir, {}
                )
                second = CAMPAIGN.prepare_herwig_run(
                    "herwig", "gg-low", card, herwig, plugin_dir, {}
                )

            self.assertEqual(first, second)
            self.assertEqual(len(calls), 1)


class YodaAccountingTest(unittest.TestCase):

    def fake_yoda(self, objects):
        return types.SimpleNamespace(read=lambda _: objects)

    def test_overflow_fractions_include_flow_bins(self):
        histogram = FakeHistogram(
            [FakeDbn(2.0, 2.0), FakeDbn(6.0, 6.0)],
            FakeDbn(0.0, 0.0),
            FakeDbn(2.0, 2.0),
        )
        fake_yoda = self.fake_yoda(
            {
                f"/RAW/QCD_INSTANTON_KKS/{name}": histogram
                for name in CAMPAIGN.JET_MASS_HISTOGRAMS
            }
        )
        with mock.patch.dict(sys.modules, {"yoda": fake_yoda}):
            result = CAMPAIGN.jets_mreco_overflow_totals(
                [Path("one.yoda")]
            )
        for histogram_name in CAMPAIGN.JET_MASS_HISTOGRAMS:
            self.assertEqual(
                result[histogram_name]["weighted_fraction"], 0.2
            )
            self.assertEqual(
                result[histogram_name]["unweighted_fraction"], 0.2
            )

    def test_validate_yoda_checks_edges_and_event_count(self):
        regular = [
            FakeDbn(1.0, 1.0)
            for _ in CAMPAIGN.MASS_EDGES[:-1]
        ]
        histogram = FakeHistogram(
            regular, FakeDbn(0.0, 0.0), FakeDbn(0.0, 0.0)
        )
        objects = {
            "/QCD_INSTANTON_KKS/_all": FakeDbn(1.0, 2.0, 1.0),
            "/RAW/QCD_INSTANTON_KKS/_all": FakeDbn(2.0, 2.0),
            "/RAW/QCD_INSTANTON_KKS/_truth_mass_valid": FakeDbn(2.0, 2.0),
        }
        for histogram_name in CAMPAIGN.JET_MASS_HISTOGRAMS:
            objects[
                f"/QCD_INSTANTON_KKS/{histogram_name}"
            ] = FakeEstimate(CAMPAIGN.MASS_EDGES)
            objects[
                f"/RAW/QCD_INSTANTON_KKS/{histogram_name}"
            ] = histogram
        objects[
            f"/QCD_INSTANTON_KKS/{CAMPAIGN.TRUTH_MASS_HISTOGRAM}"
        ] = FakeEstimate(CAMPAIGN.MASS_EDGES)
        for histogram_name in CAMPAIGN.MIGRATION_RATIO_HISTOGRAMS:
            objects[
                f"/QCD_INSTANTON_KKS/{histogram_name}"
            ] = FakeEstimate(CAMPAIGN.MIGRATION_RATIO_EDGES)
        with mock.patch.dict(
            sys.modules, {"yoda": self.fake_yoda(objects)}
        ):
            result = CAMPAIGN.validate_yoda(
                Path("sample.yoda"), expected_events=2
            )
        self.assertEqual(
            result["jets_mreco_inclusive_eta45"]["overflow_entries"],
            0.0,
        )
        self.assertEqual(
            result["jets_mreco_central"]["overflow_entries"],
            0.0,
        )


class CampaignTrackerTest(unittest.TestCase):

    def setUp(self):
        self.config = CAMPAIGN.load_config()

    def test_generator_progress_parsers(self):
        self.assertEqual(
            CAMPAIGN.parse_herwig_event_progress(
                "event> init 20\nevent> 7 20\n", 20
            ),
            (7, 20),
        )
        self.assertEqual(
            CAMPAIGN.parse_sherpa_event_progress(
                "Event 1\nEvent 13\n", 25
            ),
            (13, 25),
        )
        self.assertEqual(
            CAMPAIGN.parse_herwig_integration_progress(
                "Integrate 1 of 5\nIntegrate 4 of 5\n"
            ),
            (4, 5),
        )

    def test_monitor_payload_and_atomic_files(self):
        with tempfile.TemporaryDirectory() as temporary:
            output_root = Path(temporary) / "results/shards"
            complete = CAMPAIGN.CampaignTask(
                "herwig", "gg-low", 0, 10, output_root, False
            )
            running = CAMPAIGN.CampaignTask(
                "sherpa", "gg-low", 1, 10, output_root, False
            )
            running.directory.mkdir(parents=True)
            (running.directory / "run.log").write_text(
                "Event 1\nEvent 3\n", encoding="utf-8"
            )
            started_at = time.time() - 5.0
            payload = CAMPAIGN.build_campaign_monitor_payload(
                config=self.config,
                tasks=[complete, running],
                outcomes={complete.label: "DONE"},
                active={
                    running.label: CAMPAIGN.TaskActivity(
                        running, "sherpa", started_at
                    )
                },
                run_id="tracker-test",
                phase="running-production",
                started_at=started_at,
                workers=2,
            )

            self.assertEqual(payload["shards"]["completed"], 1)
            self.assertEqual(payload["shards"]["running"], 1)
            self.assertEqual(payload["events"]["completed"], 10)
            self.assertEqual(payload["events"]["active_generated"], 3)
            self.assertEqual(payload["active"][0]["progress"], "3/10 ( 30.0%)")

            rendered = CAMPAIGN.write_campaign_monitor(
                output_root, payload, max_listed=12
            )
            monitor = output_root.parent / "monitor"
            stored = json.loads(
                (monitor / "status.json").read_text(encoding="utf-8")
            )
            self.assertEqual(stored["run_id"], "tracker-test")
            self.assertIn("Active shards:", rendered)
            self.assertEqual(
                (monitor / "status.txt").read_text(encoding="utf-8"),
                rendered,
            )

    def test_tracker_cli_options(self):
        parser = CAMPAIGN.make_parser()
        defaults = parser.parse_args(["run"])
        self.assertEqual(defaults.progress_interval, 5.0)
        self.assertEqual(defaults.max_listed, 12)

        selected = parser.parse_args(
            [
                "smoke",
                "--progress-interval",
                "-1",
                "--max-listed",
                "3",
            ]
        )
        self.assertEqual(selected.progress_interval, -1.0)
        self.assertEqual(selected.max_listed, 3)

        status = parser.parse_args(["status", "--smoke", "--json"])
        self.assertTrue(status.smoke)
        self.assertTrue(status.json)


class CampaignExecutionTest(unittest.TestCase):

    def setUp(self):
        self.config = CAMPAIGN.load_config()

    def task(self, root, sample_id="herwig"):
        return CAMPAIGN.CampaignTask(
            sample_id, "gg-low", 0, 2, root / "shards", True
        )

    def test_complete_shard_is_skipped(self):
        with tempfile.TemporaryDirectory() as temporary:
            task = self.task(Path(temporary))
            with (
                mock.patch.object(
                    CAMPAIGN,
                    "inspect_shard",
                    return_value=("complete", ""),
                ),
                mock.patch.object(CAMPAIGN, "run_herwig") as run,
            ):
                result = CAMPAIGN.execute_task(
                    self.config, task, force=False, run_id="test"
                )
            self.assertEqual(result, "SKIP")
            run.assert_not_called()

    def test_stale_shard_requires_force(self):
        with tempfile.TemporaryDirectory() as temporary:
            task = self.task(Path(temporary))
            task.directory.mkdir(parents=True)
            (task.directory / "old").write_text("old", encoding="utf-8")
            with mock.patch.object(
                CAMPAIGN,
                "inspect_shard",
                return_value=("stale", "card_sha256"),
            ):
                with self.assertRaisesRegex(RuntimeError, "--force"):
                    CAMPAIGN.execute_task(
                        self.config, task, force=False, run_id="test"
                    )

    def test_incomplete_shard_is_archived_and_retried(self):
        with tempfile.TemporaryDirectory() as temporary:
            task = self.task(Path(temporary))
            task.directory.mkdir(parents=True)
            (task.directory / "partial").write_text(
                "partial", encoding="utf-8"
            )

            def fake_run(config, sample, selected_task):
                selected_task.directory.mkdir(parents=True)
                (selected_task.directory / "new").write_text(
                    "new", encoding="utf-8"
                )

            with (
                mock.patch.object(
                    CAMPAIGN,
                    "inspect_shard",
                    return_value=("incomplete", "missing YODA"),
                ),
                mock.patch.object(CAMPAIGN, "run_herwig", fake_run),
            ):
                result = CAMPAIGN.execute_task(
                    self.config, task, force=False, run_id="test"
                )
            self.assertEqual(result, "DONE")
            archived = (
                task.output_root.parent
                / "attempts/test/herwig/gg-low/shard-00/partial"
            )
            self.assertTrue(archived.is_file())
            self.assertTrue((task.directory / "new").is_file())

    def test_force_archives_a_complete_shard_before_replacement(self):
        with tempfile.TemporaryDirectory() as temporary:
            task = self.task(Path(temporary))
            task.directory.mkdir(parents=True)
            (task.directory / "complete").write_text(
                "old", encoding="utf-8"
            )

            def fake_run(config, sample, selected_task):
                selected_task.directory.mkdir(parents=True)
                (selected_task.directory / "complete").write_text(
                    "new", encoding="utf-8"
                )

            with (
                mock.patch.object(
                    CAMPAIGN,
                    "inspect_shard",
                    return_value=("complete", ""),
                ),
                mock.patch.object(CAMPAIGN, "run_herwig", fake_run),
            ):
                result = CAMPAIGN.execute_task(
                    self.config, task, force=True, run_id="test"
                )
            self.assertEqual(result, "DONE")
            archived = (
                task.output_root.parent
                / "superseded/test/herwig/gg-low/shard-00/complete"
            )
            self.assertEqual(archived.read_text(encoding="utf-8"), "old")
            self.assertEqual(
                (task.directory / "complete").read_text(encoding="utf-8"),
                "new",
            )

    def test_overlapping_invocations_execute_shard_once(self):
        with tempfile.TemporaryDirectory() as temporary:
            task = self.task(Path(temporary))
            calls = 0
            calls_lock = threading.Lock()

            def fake_inspect(config, sample, selected_task):
                if (selected_task.directory / "complete").is_file():
                    return "complete", ""
                return "missing", ""

            def fake_run(config, sample, selected_task):
                nonlocal calls
                with calls_lock:
                    calls += 1
                time.sleep(0.05)
                selected_task.directory.mkdir(parents=True)
                (selected_task.directory / "complete").write_text(
                    "complete", encoding="utf-8"
                )

            with (
                mock.patch.object(
                    CAMPAIGN, "inspect_shard", fake_inspect
                ),
                mock.patch.object(CAMPAIGN, "run_herwig", fake_run),
                ThreadPoolExecutor(max_workers=2) as executor,
            ):
                futures = [
                    executor.submit(
                        CAMPAIGN.execute_task,
                        self.config,
                        task,
                        force=False,
                        run_id="test",
                    )
                    for _ in range(2)
                ]
                results = sorted(future.result() for future in futures)
            self.assertEqual(results, ["DONE", "SKIP"])
            self.assertEqual(calls, 1)

    def test_run_selection_respects_worker_limit(self):
        sample = [CAMPAIGN.sample_definition(self.config, "herwig")]
        active = 0
        maximum = 0
        lock = threading.Lock()

        def fake_execute(config, task, *, force, run_id):
            nonlocal active, maximum
            with lock:
                active += 1
                maximum = max(maximum, active)
            time.sleep(0.03)
            with lock:
                active -= 1
            return "DONE"

        with (
            mock.patch.object(CAMPAIGN, "static_checks"),
            mock.patch.object(CAMPAIGN, "ensure_runtime_artifacts"),
            mock.patch.object(CAMPAIGN, "execute_task", fake_execute),
            tempfile.TemporaryDirectory() as temporary,
        ):
            CAMPAIGN.run_selection(
                config=self.config,
                samples=sample,
                profile_selection="gg-low",
                shards=range(6),
                event_override=1,
                output_root=Path(temporary),
                jobs=2,
                force=False,
            )
        self.assertEqual(maximum, 2)

    def test_run_selection_collects_failures_after_other_tasks(self):
        sample = [CAMPAIGN.sample_definition(self.config, "herwig")]
        visited = []
        lock = threading.Lock()

        def fake_execute(config, task, *, force, run_id):
            with lock:
                visited.append(task.shard)
            if task.shard == 1:
                raise RuntimeError("deliberate")
            return "DONE"

        with (
            mock.patch.object(CAMPAIGN, "static_checks"),
            mock.patch.object(CAMPAIGN, "ensure_runtime_artifacts"),
            mock.patch.object(CAMPAIGN, "execute_task", fake_execute),
            tempfile.TemporaryDirectory() as temporary,
        ):
            with self.assertRaisesRegex(RuntimeError, "1 campaign task"):
                CAMPAIGN.run_selection(
                    config=self.config,
                    samples=sample,
                    profile_selection="gg-low",
                    shards=range(4),
                    event_override=1,
                    output_root=Path(temporary),
                    jobs=2,
                    force=False,
                )
        self.assertEqual(sorted(visited), [0, 1, 2, 3])

    def test_run_selection_writes_completed_tracker(self):
        sample = [CAMPAIGN.sample_definition(self.config, "herwig")]
        with (
            mock.patch.object(CAMPAIGN, "static_checks"),
            mock.patch.object(CAMPAIGN, "ensure_runtime_artifacts"),
            mock.patch.object(
                CAMPAIGN, "execute_task", return_value="DONE"
            ),
            tempfile.TemporaryDirectory() as temporary,
        ):
            output_root = Path(temporary) / "results/shards"
            CAMPAIGN.run_selection(
                config=self.config,
                samples=sample,
                profile_selection="gg-low",
                shards=[0, 1],
                event_override=None,
                output_root=output_root,
                jobs=2,
                force=False,
                progress_interval=-1,
            )
            payload = json.loads(
                (
                    output_root.parent / "monitor/status.json"
                ).read_text(encoding="utf-8")
            )
        self.assertEqual(payload["phase"], "campaign-complete")
        self.assertEqual(payload["shards"]["completed"], 2)
        self.assertEqual(payload["shards"]["pending"], 0)

    def test_run_selection_interleaves_integrations(self):
        submitted = []

        class ImmediateExecutor:

            def __init__(self, max_workers):
                self.max_workers = max_workers

            def submit(self, function, *args, **kwargs):
                submitted.append(args[0])
                future = Future()
                future.set_result(function(*args, **kwargs))
                return future

            def shutdown(self, wait=True, cancel_futures=False):
                pass

        with (
            mock.patch.object(CAMPAIGN, "static_checks"),
            mock.patch.object(CAMPAIGN, "ensure_runtime_artifacts"),
            mock.patch.object(
                CAMPAIGN, "execute_task", return_value="DONE"
            ),
            mock.patch.object(
                CAMPAIGN, "ThreadPoolExecutor", ImmediateExecutor
            ),
            tempfile.TemporaryDirectory() as temporary,
        ):
            CAMPAIGN.run_selection(
                config=self.config,
                samples=CAMPAIGN.sample_definitions(self.config),
                profile_selection="all",
                shards=[0],
                event_override=1,
                output_root=Path(temporary),
                jobs=4,
                force=False,
            )
        self.assertEqual(
            [(task.sample_id, task.profile) for task in submitted[:4]],
            [
                ("herwig", "gg-low"),
                ("herwig", "gg-high"),
                ("herwig", "all-low"),
                ("herwig", "all-high"),
            ],
        )

    def test_run_selection_terminates_children_on_interrupt(self):
        futures = []
        shutdown_calls = []

        class PendingExecutor:

            def __init__(self, max_workers):
                self.max_workers = max_workers

            def submit(self, function, *args, **kwargs):
                future = Future()
                futures.append(future)
                return future

            def shutdown(self, wait=True, cancel_futures=False):
                shutdown_calls.append((wait, cancel_futures))

        sample = [CAMPAIGN.sample_definition(self.config, "herwig")]
        self.addCleanup(CAMPAIGN._CANCEL_REQUESTED.clear)
        with (
            mock.patch.object(CAMPAIGN, "static_checks"),
            mock.patch.object(CAMPAIGN, "ensure_runtime_artifacts"),
            mock.patch.object(
                CAMPAIGN, "ThreadPoolExecutor", PendingExecutor
            ),
            mock.patch.object(
                CAMPAIGN, "wait", side_effect=KeyboardInterrupt
            ),
            mock.patch.object(
                CAMPAIGN, "terminate_active_processes"
            ) as terminate,
            tempfile.TemporaryDirectory() as temporary,
        ):
            with self.assertRaises(KeyboardInterrupt):
                CAMPAIGN.run_selection(
                    config=self.config,
                    samples=sample,
                    profile_selection="gg-low",
                    shards=[0, 1],
                    event_override=1,
                    output_root=Path(temporary),
                    jobs=2,
                    force=False,
                )
            interrupted = json.loads(
                (
                    Path(temporary) / "monitor/status.json"
                ).read_text(encoding="utf-8")
            )
        terminate.assert_called_once_with()
        self.assertTrue(all(future.cancelled() for future in futures))
        self.assertEqual(shutdown_calls, [(True, True)])
        self.assertEqual(interrupted["phase"], "campaign-interrupted")

    def test_sample_and_generator_arguments_are_mutually_exclusive(self):
        args = argparse.Namespace(
            sample="herwig-qcdinsplanar", generator="herwig"
        )
        with self.assertRaisesRegex(RuntimeError, "either"):
            CAMPAIGN.selected_samples(self.config, args)


if __name__ == "__main__":
    unittest.main()
