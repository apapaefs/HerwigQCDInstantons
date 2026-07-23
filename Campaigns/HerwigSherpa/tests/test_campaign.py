#!/usr/bin/env python3

import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest import mock


CAMPAIGN_DIR = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "run_campaign", CAMPAIGN_DIR / "run_campaign.py"
)
CAMPAIGN = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(CAMPAIGN)


class CampaignConfigurationTest(unittest.TestCase):

    def setUp(self):
        self.config = CAMPAIGN.load_config()

    def test_static_validation(self):
        CAMPAIGN.static_checks(self.config)

    def test_profile_event_budget(self):
        self.assertEqual(
            CAMPAIGN.profile_names(self.config),
            ["gg-low", "gg-high", "all-low", "all-high"],
        )
        total = (
            sum(region["events"] for region in self.config["regions"].values())
            * len(self.config["processes"])
            * len(self.config["generators"])
        )
        self.assertEqual(total, 280000)

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

    def test_generated_cards_are_current(self):
        for path, content in CAMPAIGN.rendered_cards(self.config).items():
            self.assertTrue(path.is_file(), path)
            self.assertEqual(path.read_text(encoding="utf-8"), content)

    def test_shard_seeds_are_unique(self):
        seeds = {
            CAMPAIGN.stable_seed(generator, profile, shard)
            for generator in self.config["generators"]
            for profile in CAMPAIGN.profile_names(self.config)
            for shard in range(self.config["shards"])
        }
        self.assertEqual(len(seeds), 80)

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

    def test_herwig_integration_is_cached(self):
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
                    "gg-low", card, herwig, plugin_dir, {}
                )
                second = CAMPAIGN.prepare_herwig_run(
                    "gg-low", card, herwig, plugin_dir, {}
                )

            self.assertEqual(first, second)
            self.assertEqual(len(calls), 1)


if __name__ == "__main__":
    unittest.main()
