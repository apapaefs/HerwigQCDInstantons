#!/usr/bin/env python3
"""Reproducible Herwig/Sherpa QCD-instanton comparison campaign."""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shlex
import shutil
import subprocess
import sys
from typing import Dict, Iterable, List, Sequence, Tuple


CAMPAIGN_DIR = Path(__file__).resolve().parent
ROOT = CAMPAIGN_DIR.parents[1]
CONFIG_PATH = CAMPAIGN_DIR / "campaign.json"
BUILD_DIR = CAMPAIGN_DIR / ".build"
LOCAL_DIR = CAMPAIGN_DIR / ".local"
WORK_DIR = CAMPAIGN_DIR / ".work"
RESULTS_DIR = CAMPAIGN_DIR / "results"

KKS_TABLE: Tuple[Tuple[float, float, float, float, float], ...] = (
    (10.7, 0.99, 0.416, 4.59, 4.922e9),
    (11.4, 1.04, 0.405, 4.68, 3.652e9),
    (13.4, 1.16, 0.382, 4.90, 1.671e9),
    (15.7, 1.31, 0.360, 5.13, 728.9e6),
    (22.9, 1.76, 0.315, 5.44, 85.94e6),
    (29.7, 2.12, 0.293, 6.02, 17.25e6),
    (40.8, 2.72, 0.267, 6.47, 2.121e6),
    (56.1, 3.50, 0.245, 6.92, 229.0e3),
    (61.8, 3.64, 0.223, 7.28, 72.97e3),
    (89.6, 4.98, 0.206, 7.67, 2.733e3),
    (118.0, 6.21, 0.195, 8.25, 235.4),
    (174.4, 8.72, 0.180, 8.60, 6.720),
    (246.9, 11.76, 0.169, 9.04, 0.284),
    (349.9, 15.90, 0.159, 9.49, 0.012),
    (496.3, 21.58, 0.150, 9.93, 5.112e-4),
    (704.8, 29.37, 0.142, 10.37, 21.65e-6),
    (1001.8, 40.07, 0.135, 10.81, 0.9017e-6),
    (1425.6, 54.83, 0.128, 11.26, 36.45e-9),
    (2030.6, 75.21, 0.122, 11.70, 1.419e-9),
    (2895.5, 103.4, 0.117, 12.14, 52.07e-12),
)


def load_config() -> Dict[str, object]:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def profile_names(config: Dict[str, object]) -> List[str]:
    return [
        f"{process}-{region}"
        for process in config["processes"]
        for region in config["regions"]
    ]


def split_profile(profile: str) -> Tuple[str, str]:
    process, region = profile.split("-", 1)
    return process, region


def stable_seed(generator: str, profile: str, shard: int) -> int:
    key = f"{generator}:{profile}:{shard}".encode("ascii")
    return int(hashlib.sha256(key).hexdigest()[:8], 16) % 900000000 + 1


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def shell_join(command: Sequence[object]) -> str:
    return " ".join(shlex.quote(str(item)) for item in command)


def run_checked(
    command: Sequence[object],
    *,
    cwd: Path | None = None,
    env: Dict[str, str] | None = None,
    log_path: Path | None = None,
) -> None:
    command_text = shell_join(command)
    print(f"+ {command_text}", flush=True)
    if log_path is None:
        subprocess.run(
            [str(item) for item in command],
            cwd=cwd,
            env=env,
            check=True,
        )
        return

    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        log.write(f"$ {command_text}\n")
        log.flush()
        subprocess.run(
            [str(item) for item in command],
            cwd=cwd,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )


def capture(command: Sequence[object], *, allow_nonzero: bool = False) -> str:
    result = subprocess.run(
        [str(item) for item in command],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if result.returncode != 0 and not allow_nonzero:
        raise subprocess.CalledProcessError(
            result.returncode, [str(item) for item in command], result.stdout
        )
    return result.stdout.strip()


def find_executable(name: str, env_name: str | None = None) -> Path:
    if env_name and os.environ.get(env_name):
        candidate = Path(os.environ[env_name]).expanduser().resolve()
        if candidate.is_file():
            return candidate
        raise RuntimeError(f"{env_name} does not name an executable: {candidate}")
    located = shutil.which(name)
    if located:
        return Path(located).resolve()
    raise RuntimeError(f"Cannot find {name}; activate the generator environment.")


def find_sherpa_source(config: Dict[str, object]) -> Path:
    if os.environ.get("SHERPA_SOURCE"):
        source = Path(os.environ["SHERPA_SOURCE"]).expanduser().resolve()
    else:
        version = config["sherpa"]["version"]
        source = (
            Path.home()
            / "Projects/Pull/summer2026_new/Sherpa"
            / f"sherpa-v{version}"
        )
    if not (source / "CMakeLists.txt").is_file():
        raise RuntimeError(
            f"Cannot find Sherpa source at {source}; set SHERPA_SOURCE."
        )
    return source


def find_sherpa_binary(config: Dict[str, object]) -> Path:
    if os.environ.get("SHERPA_BIN"):
        return find_executable("Sherpa", "SHERPA_BIN")
    local = LOCAL_DIR / f"sherpa-{config['sherpa']['version']}" / "bin/Sherpa"
    if local.is_file():
        return local.resolve()
    source = find_sherpa_source(config)
    existing = source / "build-gcc-openloops/install/bin/Sherpa"
    if existing.is_file():
        return existing.resolve()
    return find_executable("Sherpa")


def find_campaign_sherpa_binary(config: Dict[str, object]) -> Path:
    """Locate the patched campaign build, unless explicitly overridden."""
    if os.environ.get("SHERPA_BIN"):
        return find_executable("Sherpa", "SHERPA_BIN")
    local = LOCAL_DIR / f"sherpa-{config['sherpa']['version']}" / "bin/Sherpa"
    if local.is_file():
        return local.resolve()
    raise RuntimeError(
        "Build the campaign-local Sherpa before running, or set SHERPA_BIN "
        "to an explicitly patched compatible executable."
    )


def prepare_sherpa_source(config: Dict[str, object]) -> Path:
    source = find_sherpa_source(config)
    mirror = BUILD_DIR / f"sherpa-source-{config['sherpa']['version']}"
    if not mirror.exists():
        shutil.copytree(
            source,
            mirror,
            ignore=shutil.ignore_patterns(
                ".git",
                "build",
                "build-*",
                "__pycache__",
                "*.pyc",
            ),
        )

    find_rivet = mirror / "cmake/Modules/FindRivet.cmake"
    text = find_rivet.read_text(encoding="utf-8")
    old = 'set(RIVET_HEPMC3_MIN_VERSION "3.2.6")'
    compatible = 'set(RIVET_HEPMC3_MIN_VERSION "3.2.5")'
    override = f"    {compatible}\n"
    marker = "    set(RIVET_MKHTML_ARGS )"
    if compatible not in text:
        if old not in text or marker not in text:
            raise RuntimeError(
                "Cannot apply the Sherpa HepMC3 compatibility patch."
            )
        find_rivet.write_text(
            text.replace(marker, override + marker, 1),
            encoding="utf-8",
        )

    rambo = mirror / "PHASIC++/Channels/Rambo.C"
    text = rambo.read_text(encoding="utf-8")
    old = """  if (E<p_ms[0]+p_ms[1]) THROW(fatal_error, "sqrt(s) smaller than particle masses");
  double x=1.0/2.0+(p_ms[0]*p_ms[0]-p_ms[1]*p_ms[1])/(2.0*E*E);
  p[0]=ATOOLS::Vec4D(x*E,0.0,0.0,sqrt(ATOOLS::sqr(x*E)-p_ms[0]*p_ms[0]));"""
    new = """  if (E<sqrt(p_ms[0])+sqrt(p_ms[1]))
    THROW(fatal_error, "sqrt(s) smaller than particle masses");
  double x=1.0/2.0+(p_ms[0]-p_ms[1])/(2.0*E*E);
  p[0]=ATOOLS::Vec4D(x*E,0.0,0.0,sqrt(ATOOLS::sqr(x*E)-p_ms[0]));"""
    if old in text:
        rambo.write_text(text.replace(old, new, 1), encoding="utf-8")
    elif new not in text:
        raise RuntimeError("Cannot apply the Sherpa RAMBO mass patch.")

    instanton = mirror / "EXTRA_XS/Special/Instanton.C"
    text = instanton.read_text(encoding="utf-8")
    old = """    if (E>m_Ehatmax || E<m_Ehatmin) {
      msg_Debugging()<<METHOD<<" yields false for E = "<<E<<".\\n";
      return false;
    }
    std::map<double, xsec_data*>::iterator dit;"""
    new = """    if (E>m_Ehatmax || E<m_Ehatmin) {
      msg_Debugging()<<METHOD<<" yields false for E = "<<E<<".\\n";
      return false;
    }
    if (E==m_Ehatmax) {
      m_lower = m_upper = m_data.rbegin()->second;
      m_rho      = m_upper->m_rho;
      m_Ngluons  = m_upper->m_Ngluons;
      m_sigmahat = m_upper->m_sigmahat;
      return true;
    }
    std::map<double, xsec_data*>::iterator dit;"""
    if old in text:
        instanton.write_text(text.replace(old, new, 1), encoding="utf-8")
    elif new not in text:
        raise RuntimeError("Cannot apply the Sherpa instanton endpoint patch.")
    return mirror


def sherpa_processes(process: str) -> List[Tuple[int, int]]:
    if process == "gg":
        return [(21, 21)]
    if process != "all":
        raise ValueError(f"Unknown process profile: {process}")

    result: List[Tuple[int, int]] = [(21, 21)]
    fermions = list(range(1, 6)) + list(range(-1, -6, -1))
    result.extend((21, flavour) for flavour in fermions)
    for first in range(1, 6):
        for second in range(first + 1, 6):
            result.append((first, second))
            result.append((-first, -second))
    result.extend((quark, -antiquark) for quark in range(1, 6)
                  for antiquark in range(1, 6))
    return result


def render_herwig_profile(
    config: Dict[str, object],
    process: str,
    region: str,
    *,
    smoke: bool = False,
) -> str:
    values = config["regions"][region]
    profile = f"{process}-{region}"
    process_option = "GG" if process == "gg" else "All"
    smoke_sampler = ""
    if smoke:
        smoke_sampler = """
# Keep complete process coverage while making smoke initialization inexpensive.
set /Herwig/Samplers/FlatBinSampler:NIterations 1
set /Herwig/Samplers/FlatBinSampler:InitialPoints 500
"""
    return f"""# -*- ThePEG-repository -*-
# Generated by run_campaign.py render; edit campaign.json or common.in.

read Campaigns/HerwigSherpa/cards/herwig/common.in

set /Herwig/Generators/EventGenerator:NumberOfEvents {values['events']}
set /Herwig/Generators/EventGenerator:RandomNumberGenerator:Seed {stable_seed('herwig', profile, 0)}
set /Herwig/MatrixElements/MEInstanton:Processes {process_option}
set /Herwig/Cuts/InstantonCuts:MHatMin {values['min_mass_gev']}*GeV
set /Herwig/Cuts/InstantonCuts:MHatMax {values['max_mass_gev']}*GeV
set /Herwig/Analysis/Rivet:Filename herwig-{profile}.yoda
{smoke_sampler}
cd /Herwig/Generators
saverun Campaign-Herwig-{profile} EventGenerator
"""


def format_kks_table() -> str:
    rows = []
    for energy, rho, alpha_s, mean, cross_section in KKS_TABLE:
        rows.append(
            "  "
            + f"[{energy:.6g}, {rho:.6g}, {alpha_s:.6g}, "
            + f"{mean:.6g}, {cross_section:.12g}]"
        )
    return ",\n".join(rows)


def render_sherpa_profile(
    config: Dict[str, object], process: str, region: str
) -> str:
    values = config["regions"][region]
    masses = config["masses_gev"]
    sherpa = config["sherpa"]
    profile = f"{process}-{region}"
    process_lines = "\n".join(
        f'  - "{first} {second} -> 999":\n'
        "      Order: {EW: 2, QCD: Any}"
        for first, second in sherpa_processes(process)
    )
    return f"""# Generated by run_campaign.py render; edit campaign.json instead.

BEAMS: 2212
BEAM_ENERGIES: {config['beam_energy_gev']}
EVENTS: {values['events']}
RANDOM_SEED: {stable_seed('sherpa', profile, 0)}

PDF_LIBRARY: [LHAPDFSherpa]
PDF_SET: [{config['pdf']['set']}]
PDF_SET_VERSIONS: [{config['pdf']['member']}]
MPI_PDF_LIBRARY: [LHAPDFSherpa]
MPI_PDF_SET: [{config['pdf']['set']}]
MPI_PDF_SET_VERSIONS: [{config['pdf']['member']}]
ALPHAS: {{USE_PDF: 1}}
FREEZE_PDF_FOR_LOW_Q: 1

MI_HANDLER: None
SHOWER_GENERATOR: CSS
FRAGMENTATION: Ahadic
HADRON_DECAYS:
  Max_Proper_Lifetime: 10.0
HARD_DECAYS:
  Enabled: false

PARTICLE_DATA:
  1: {{Mass: {masses['d']}}}
  2: {{Mass: {masses['u']}}}
  3: {{Mass: {masses['s']}}}
  4: {{Mass: {masses['c']}}}
  5: {{Mass: {masses['b']}}}

SCALES: Democratic
INSTANTON_XSECS: [
{format_kks_table()}
]
INSTANTON_MIN_MASS: {values['min_mass_gev']}
INSTANTON_MAX_MASS: {values['max_mass_gev']}
INSTANTON_NGLUONS_MODIFIER: 1.0
INSTANTON_SIGMAHAT_MODIFIER: {sherpa['sigmahat_modifier']}
INSTANTON_ALPHAS_FACTOR: 1.0
INSTANTON_SCALE_CHOICE: shat/N
INSTANTON_SCALE_FACTOR: 1.0
INSTANTON_INCLUDE_QUARKS: {sherpa['include_quarks']}
INSTANTON_C_PRODUCTION_THRESHOLD: {sherpa['charm_threshold_gev']}
INSTANTON_B_PRODUCTION_THRESHOLD: {sherpa['bottom_threshold_gev']}

PROCESSES:
{process_lines}

ANALYSIS: Rivet
ANALYSIS_OUTPUT: sherpa-{profile}
RIVET:
  --analyses: [{config['analysis']}]
"""


def rendered_cards(config: Dict[str, object]) -> Dict[Path, str]:
    result: Dict[Path, str] = {}
    for process in config["processes"]:
        for region in config["regions"]:
            profile = f"{process}-{region}"
            result[
                CAMPAIGN_DIR / f"cards/herwig/{profile}.in"
            ] = render_herwig_profile(config, process, region)
            result[
                CAMPAIGN_DIR / f"cards/sherpa/{profile}.yaml"
            ] = render_sherpa_profile(config, process, region)
    return result


def command_render(args: argparse.Namespace) -> None:
    config = load_config()
    changed: List[Path] = []
    for path, content in rendered_cards(config).items():
        current = path.read_text(encoding="utf-8") if path.exists() else None
        if current == content:
            continue
        changed.append(path)
        if not args.check:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content, encoding="utf-8")
    if args.check and changed:
        relative = ", ".join(str(path.relative_to(ROOT)) for path in changed)
        raise RuntimeError(f"Rendered cards are stale: {relative}")
    if not args.check:
        print(f"Rendered {len(rendered_cards(config))} cards; "
              f"{len(changed)} changed.")


def static_checks(config: Dict[str, object], require_cards: bool = True) -> None:
    all_processes = sherpa_processes("all")
    if len(all_processes) != 56 or len(set(all_processes)) != 56:
        raise RuntimeError("Sherpa All must contain 56 unique incoming states.")
    for first, second in all_processes:
        if first == second and abs(first) <= 5:
            raise RuntimeError("Equal-flavour qq or qbar-qbar process found.")

    total_events = (
        sum(values["events"] for values in config["regions"].values())
        * len(config["processes"])
        * len(config["generators"])
    )
    if total_events != 280000:
        raise RuntimeError(f"Unexpected campaign event total: {total_events}")
    if len(KKS_TABLE) != 20 or KKS_TABLE[0][-1] != 4.922e9:
        raise RuntimeError("The Sherpa KKS table is incomplete or inconsistent.")

    cap = config["max_final_partons"]
    expected = {
        (4, 2): 21,
        (4, 1): 22,
        (4, 0): 23,
        (5, 2): 19,
        (5, 1): 20,
        (5, 0): 21,
    }
    for (flavours, incoming_gluons), maximum in expected.items():
        outgoing_quarks = 2 * flavours - (2 - incoming_gluons)
        calculated = min(config["gluon_cap"], cap - outgoing_quarks)
        if calculated != maximum:
            raise RuntimeError("Process-dependent gluon cap is inconsistent.")

    if require_cards:
        for path, content in rendered_cards(config).items():
            if not path.is_file() or path.read_text(encoding="utf-8") != content:
                raise RuntimeError(f"Card is missing or stale: {path}")


def command_doctor(args: argparse.Namespace) -> None:
    config = load_config()
    static_checks(config)
    source = find_sherpa_source(config)
    tools = {
        "Herwig": [find_executable("Herwig", "HERWIG_BIN"), "--version"],
        "Sherpa": [find_sherpa_binary(config), "--version"],
        "Rivet": [find_executable("rivet-config"), "--version"],
        "LHAPDF": [find_executable("lhapdf-config"), "--version"],
        "CMake": [find_executable("cmake"), "--version"],
    }
    report = {"sherpa_source": str(source), "tools": {}}
    for name, command in tools.items():
        output = capture(command)
        report["tools"][name] = {
            "path": str(command[0]),
            "version": output.splitlines()[0],
        }
    lhapdf = find_executable("lhapdf")
    available_sets = capture([lhapdf, "ls"]).splitlines()
    pdf_name = config["pdf"]["set"]
    if pdf_name not in available_sets:
        raise RuntimeError(f"Required PDF set is unavailable: {pdf_name}")
    report["pdf"] = {"set": pdf_name, "member": config["pdf"]["member"]}
    print(json.dumps(report, indent=2, sort_keys=True))


def detect_cxx() -> List[str]:
    if os.environ.get("CXX"):
        configured = shlex.split(os.environ["CXX"])
        located = shutil.which(configured[0])
        if located:
            return [str(Path(located).resolve()), *configured[1:]]
    for name in ("g++-16", "g++-15", "g++-14", "g++", "c++"):
        located = shutil.which(name)
        if located:
            return [str(Path(located).resolve())]
    raise RuntimeError("Cannot find a C++ compiler.")


def build_herwig_plugin() -> None:
    output_dir = BUILD_DIR / "herwig"
    output_dir.mkdir(parents=True, exist_ok=True)
    cxx = detect_cxx()
    cppflags = shlex.split(
        capture(
            [find_executable("Herwig-config"), "--cppflags"],
            allow_nonzero=True,
        )
    )
    gsl_libs = shlex.split(capture([find_executable("gsl-config"), "--libs"]))
    if platform.system() == "Darwin":
        link_flags = ["-bundle", "-Wl,-undefined,dynamic_lookup"]
    else:
        link_flags = ["-shared"]
    command = [
        *cxx,
        "-std=gnu++17",
        "-O2",
        "-fPIC",
        *link_flags,
        *cppflags,
        ROOT / "MEInstanton.cc",
        "-o",
        output_dir / "CampaignInstantons.so",
        *gsl_libs,
    ]
    run_checked(command, cwd=ROOT)


def build_rivet_plugin() -> None:
    output_dir = BUILD_DIR / "rivet"
    output_dir.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["CXX"] = detect_cxx()[0]
    run_checked(
        [
            find_executable("rivet-build"),
            output_dir / "RivetQCD_INSTANTON_KKS.so",
            ROOT / "Rivet/QCD_INSTANTON_KKS.cc",
        ],
        cwd=ROOT,
        env=env,
    )


def build_sherpa(config: Dict[str, object], jobs: int) -> None:
    source = prepare_sherpa_source(config)
    build = BUILD_DIR / f"sherpa-build-{config['sherpa']['version']}"
    install = LOCAL_DIR / f"sherpa-{config['sherpa']['version']}"
    prefix = capture(
        [find_executable("Herwig-config"), "--prefix"],
        allow_nonzero=True,
    )
    cxx = detect_cxx()[0]
    cc = next(
        (shutil.which(name) for name in ("gcc-16", "gcc-15", "gcc", "cc")
         if shutil.which(name)),
        None,
    )
    fc = next(
        (shutil.which(name)
         for name in ("gfortran-16", "gfortran-15", "gfortran")
         if shutil.which(name)),
        None,
    )
    configure = [
        find_executable("cmake"),
        "-S",
        source,
        "-B",
        build,
        f"-DCMAKE_INSTALL_PREFIX={install}",
        f"-DCMAKE_PREFIX_PATH={prefix}",
        f"-DCMAKE_CXX_COMPILER={cxx}",
        "-DSHERPA_ENABLE_RIVET=ON",
        "-DSHERPA_ENABLE_HEPMC3=ON",
        "-DSHERPA_ENABLE_LHAPDF=ON",
        "-DSHERPA_ENABLE_OPENLOOPS=OFF",
        "-DSHERPA_ENABLE_EXAMPLES=OFF",
        "-DSHERPA_ENABLE_TESTING=OFF",
        "-DSHERPA_ENABLE_MPI=OFF",
    ]
    if cc:
        configure.append(f"-DCMAKE_C_COMPILER={cc}")
    if fc:
        configure.append(f"-DCMAKE_Fortran_COMPILER={fc}")
    run_checked(configure, cwd=ROOT)
    run_checked(
        [find_executable("cmake"), "--build", build, "--parallel", jobs],
        cwd=ROOT,
    )
    run_checked(
        [find_executable("cmake"), "--install", build],
        cwd=ROOT,
    )


def command_build(args: argparse.Namespace) -> None:
    config = load_config()
    components = (
        ("herwig", "rivet", "sherpa")
        if args.component == "all"
        else (args.component,)
    )
    if "herwig" in components:
        build_herwig_plugin()
    if "rivet" in components:
        build_rivet_plugin()
    if "sherpa" in components:
        build_sherpa(config, args.jobs)


def campaign_environment() -> Dict[str, str]:
    env = os.environ.copy()
    rivet_paths = [str(BUILD_DIR / "rivet")]
    if env.get("RIVET_ANALYSIS_PATH"):
        rivet_paths.append(env["RIVET_ANALYSIS_PATH"])
    env["RIVET_ANALYSIS_PATH"] = os.pathsep.join(rivet_paths)
    return env


def write_provenance(
    directory: Path,
    *,
    generator: str,
    profile: str,
    shard: int,
    events: int,
    seed: int,
    card: Path,
    command: Sequence[object],
    run_file: Path | None = None,
) -> None:
    git_commit = capture(["git", "-C", ROOT, "rev-parse", "HEAD"])
    provenance = {
        "generator": generator,
        "profile": profile,
        "shard": shard,
        "events": events,
        "seed": seed,
        "card": str(card.relative_to(ROOT)),
        "card_sha256": sha256(card),
        "command": [str(item) for item in command],
        "git_commit": git_commit,
    }
    executable = Path(command[0])
    if executable.is_file():
        provenance["generator_executable"] = str(executable.resolve())
        provenance["generator_executable_sha256"] = sha256(executable)
        version = capture(
            [executable, "--version"], allow_nonzero=True
        ).splitlines()
        if version:
            provenance["generator_version"] = version[0]
    for name, path in (
        ("instanton_plugin", BUILD_DIR / "herwig/CampaignInstantons.so"),
        ("rivet_plugin", BUILD_DIR / "rivet/RivetQCD_INSTANTON_KKS.so"),
        ("herwig_run", run_file),
    ):
        if path is not None and path.is_file():
            provenance[f"{name}_sha256"] = sha256(path)
    (directory / "provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def prepare_herwig_run(
    profile: str,
    card: Path,
    herwig: Path,
    plugin_path: Path,
    env: Dict[str, str],
) -> Path:
    """Create or reuse the production integration for one Herwig profile."""
    rivet_plugin = BUILD_DIR / "rivet/RivetQCD_INSTANTON_KKS.so"
    common_card = CAMPAIGN_DIR / "cards/herwig/common.in"
    if not rivet_plugin.is_file():
        raise RuntimeError("Build the Rivet plugin before running.")

    inputs = {
        "card": sha256(card),
        "common_card": sha256(common_card),
        "herwig": str(herwig),
        "herwig_version": capture([herwig, "--version"]).splitlines()[0],
        "instanton_plugin": sha256(plugin_path / "CampaignInstantons.so"),
        "rivet_plugin": sha256(rivet_plugin),
    }
    fingerprint = hashlib.sha256(
        json.dumps(inputs, sort_keys=True).encode("utf-8")
    ).hexdigest()
    directory = WORK_DIR / "herwig-runs" / profile / fingerprint[:16]
    directory.mkdir(parents=True, exist_ok=True)
    run_file = directory / f"Campaign-Herwig-{profile}.run"

    # Separate campaign invocations can request different shards concurrently.
    # Serialize the one-time integration without constraining event generation.
    with (directory / ".integration.lock").open("w", encoding="utf-8") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if not run_file.is_file():
            read_command = [
                herwig,
                "read",
                "-I",
                ROOT,
                "-L",
                plugin_path,
                card,
            ]
            run_checked(
                read_command,
                cwd=directory,
                env=env,
                log_path=directory / "read.log",
            )
            (directory / "integration.json").write_text(
                json.dumps(
                    {
                        "fingerprint": fingerprint,
                        "inputs": inputs,
                        "command": [str(item) for item in read_command],
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n",
                encoding="utf-8",
            )
        else:
            print(f"Reusing Herwig integration {run_file}", flush=True)
        fcntl.flock(lock, fcntl.LOCK_UN)
    return run_file


def run_herwig(
    config: Dict[str, object],
    profile: str,
    shard: int,
    events: int,
    directory: Path,
    *,
    smoke: bool,
) -> None:
    canonical_card = CAMPAIGN_DIR / f"cards/herwig/{profile}.in"
    herwig = find_executable("Herwig", "HERWIG_BIN")
    plugin_path = BUILD_DIR / "herwig"
    if not (plugin_path / "CampaignInstantons.so").is_file():
        raise RuntimeError("Build the Herwig plugin before running.")
    env = campaign_environment()
    directory.mkdir(parents=True, exist_ok=True)
    card = canonical_card
    if smoke:
        process, region = split_profile(profile)
        card = directory / f"{profile}-smoke.in"
        card.write_text(
            render_herwig_profile(
                config,
                process,
                region,
                smoke=True,
            ),
            encoding="utf-8",
        )
        read_command = [
            herwig,
            "read",
            "-I",
            ROOT,
            "-L",
            plugin_path,
            card,
        ]
        run_checked(
            read_command,
            cwd=directory,
            env=env,
            log_path=directory / "read.log",
        )
        run_file = directory / f"Campaign-Herwig-{profile}.run"
    else:
        run_file = prepare_herwig_run(
            profile,
            card,
            herwig,
            plugin_path,
            env,
        )
    seed = stable_seed("herwig", profile, shard)
    run_command = [
        herwig,
        "run",
        "-L",
        plugin_path,
        "-N",
        events,
        "-s",
        seed,
        run_file,
    ]
    write_provenance(
        directory,
        generator="herwig",
        profile=profile,
        shard=shard,
        events=events,
        seed=seed,
        card=card,
        command=run_command,
        run_file=run_file,
    )
    run_checked(
        run_command,
        cwd=directory,
        env=env,
        log_path=directory / "run.log",
    )
    validate_yoda(locate_yoda(directory))


def run_sherpa(
    config: Dict[str, object],
    profile: str,
    shard: int,
    events: int,
    directory: Path,
) -> None:
    card = CAMPAIGN_DIR / f"cards/sherpa/{profile}.yaml"
    sherpa = find_campaign_sherpa_binary(config)
    directory.mkdir(parents=True, exist_ok=True)
    seed = stable_seed("sherpa", profile, shard)
    run_command = [
        sherpa,
        "-e",
        events,
        "-R",
        seed,
        "-r",
        "Results",
        "-A",
        f"sherpa-{profile}",
        card,
    ]
    write_provenance(
        directory,
        generator="sherpa",
        profile=profile,
        shard=shard,
        events=events,
        seed=seed,
        card=card,
        command=run_command,
    )
    run_checked(
        run_command,
        cwd=directory,
        env=campaign_environment(),
        log_path=directory / "run.log",
    )
    validate_yoda(locate_yoda(directory))


def select_profiles(config: Dict[str, object], selection: str) -> List[str]:
    available = profile_names(config)
    if selection == "all":
        return available
    if selection not in available:
        raise RuntimeError(f"Unknown profile {selection}; choose {available}.")
    return [selection]


def select_generators(config: Dict[str, object], selection: str) -> List[str]:
    if selection == "all":
        return list(config["generators"])
    if selection not in config["generators"]:
        raise RuntimeError(f"Unknown generator: {selection}")
    return [selection]


def run_selection(
    *,
    config: Dict[str, object],
    generator_selection: str,
    profile_selection: str,
    shards: Iterable[int],
    event_override: int | None,
    output_root: Path,
) -> None:
    for generator in select_generators(config, generator_selection):
        for profile in select_profiles(config, profile_selection):
            _, region = split_profile(profile)
            production_events = config["regions"][region]["events"]
            per_shard = production_events // config["shards"]
            events = event_override if event_override is not None else per_shard
            for shard in shards:
                directory = (
                    output_root
                    / generator
                    / profile
                    / f"shard-{shard:02d}"
                )
                if generator == "herwig":
                    run_herwig(
                        config,
                        profile,
                        shard,
                        events,
                        directory,
                        smoke=event_override is not None,
                    )
                else:
                    run_sherpa(config, profile, shard, events, directory)


def command_smoke(args: argparse.Namespace) -> None:
    config = load_config()
    run_selection(
        config=config,
        generator_selection=args.generator,
        profile_selection=args.profile,
        shards=[0],
        event_override=args.events,
        output_root=WORK_DIR / "smoke",
    )


def parse_shards(value: str, total: int) -> List[int]:
    if value == "all":
        return list(range(total))
    result = sorted({int(item) for item in value.split(",")})
    if not result or result[0] < 0 or result[-1] >= total:
        raise RuntimeError(f"Shard selection must lie in 0..{total - 1}.")
    return result


def command_run(args: argparse.Namespace) -> None:
    config = load_config()
    run_selection(
        config=config,
        generator_selection=args.generator,
        profile_selection=args.profile,
        shards=parse_shards(args.shards, config["shards"]),
        event_override=None,
        output_root=RESULTS_DIR / "shards",
    )


def locate_yoda(directory: Path) -> Path:
    candidates = sorted(directory.glob("*.yoda")) + sorted(
        directory.glob("*.yoda.gz")
    )
    if len(candidates) != 1:
        raise RuntimeError(
            f"Expected one YODA file in {directory}, found {len(candidates)}."
        )
    return candidates[0]


def command_merge(args: argparse.Namespace) -> None:
    config = load_config()
    output_dir = RESULTS_DIR / "merged"
    output_dir.mkdir(parents=True, exist_ok=True)
    yodamerge = find_executable("yodamerge")
    for generator in select_generators(config, args.generator):
        for profile in select_profiles(config, args.profile):
            inputs = [
                locate_yoda(
                    RESULTS_DIR
                    / "shards"
                    / generator
                    / profile
                    / f"shard-{shard:02d}"
                )
                for shard in range(config["shards"])
            ]
            output = output_dir / f"{generator}-{profile}.yoda"
            run_checked([yodamerge, "-o", output, *inputs], cwd=ROOT)


def command_plot(args: argparse.Namespace) -> None:
    config = load_config()
    rivet_mkhtml = find_executable("rivet-mkhtml")
    env = campaign_environment()
    for profile in select_profiles(config, args.profile):
        herwig = RESULTS_DIR / f"merged/herwig-{profile}.yoda"
        sherpa = RESULTS_DIR / f"merged/sherpa-{profile}.yoda"
        if not herwig.is_file() or not sherpa.is_file():
            raise RuntimeError(f"Missing merged YODA files for {profile}.")
        for output_format in ("PDF", "PNG"):
            output = RESULTS_DIR / "plots" / output_format.lower() / profile
            run_checked(
                [
                    rivet_mkhtml,
                    "-o",
                    output,
                    "-f",
                    output_format,
                    "--pwd",
                    f"{herwig}:Title=Herwig VariableKKS",
                    f"{sherpa}:Title=Sherpa threshold variable",
                ],
                cwd=ROOT,
                env=env,
            )


FLOAT_PATTERN = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
ANSI_ESCAPE = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")
HERWIG_XSEC = re.compile(
    rf"estimated total cross section is\s*\(\s*({FLOAT_PATTERN})"
    rf"\s*\+/-\s*({FLOAT_PATTERN})\s*\)\s*(pb|nb|fb)",
    re.IGNORECASE,
)
SHERPA_PROCESS_XSEC = re.compile(
    rf"^(2_1__\S+)\s*:\s*({FLOAT_PATTERN})\s*pb"
    rf"\s*\+-\s*\(\s*({FLOAT_PATTERN})\s*pb",
    re.MULTILINE,
)


def parse_herwig_cross_section(path: Path) -> Tuple[float, float]:
    matches = HERWIG_XSEC.findall(path.read_text(encoding="utf-8"))
    if not matches:
        raise RuntimeError(f"Cannot find a Herwig cross section in {path}.")
    value, uncertainty, unit = matches[-1]
    to_pb = {"pb": 1.0, "nb": 1000.0, "fb": 0.001}[unit.lower()]
    return float(value) * to_pb, float(uncertainty) * to_pb


def parse_sherpa_cross_section(path: Path) -> Tuple[float, float]:
    text = ANSI_ESCAPE.sub("", path.read_text(encoding="utf-8"))
    processes: Dict[str, Tuple[float, float]] = {}
    for name, value, uncertainty in SHERPA_PROCESS_XSEC.findall(text):
        processes[name] = (float(value), float(uncertainty))
    if not processes:
        raise RuntimeError(f"Cannot find Sherpa process cross sections in {path}.")
    value = sum(result[0] for result in processes.values())
    uncertainty = math.sqrt(sum(result[1] ** 2 for result in processes.values()))
    return value, uncertainty


def combine_cross_sections(
    estimates: Sequence[Tuple[float, float]],
) -> Tuple[float, float]:
    if not estimates:
        raise RuntimeError("No cross-section estimates were supplied.")
    weighted = [
        (value, uncertainty)
        for value, uncertainty in estimates
        if uncertainty > 0.0 and math.isfinite(uncertainty)
    ]
    if len(weighted) == len(estimates):
        denominator = sum(1.0 / uncertainty**2 for _, uncertainty in weighted)
        numerator = sum(
            value / uncertainty**2 for value, uncertainty in weighted
        )
        return numerator / denominator, math.sqrt(1.0 / denominator)
    return sum(value for value, _ in estimates) / len(estimates), 0.0


def herwig_integration_log(shard: Path, smoke: bool) -> Path:
    if smoke:
        return shard / "read.log"
    provenance = json.loads(
        (shard / "provenance.json").read_text(encoding="utf-8")
    )
    run_file = Path(provenance["command"][-1])
    return run_file.parent / "read.log"


def command_summarize(args: argparse.Namespace) -> None:
    config = load_config()
    base = WORK_DIR / "smoke" if args.smoke else RESULTS_DIR / "shards"
    report = {
        "mode": "smoke" if args.smoke else "production",
        "units": "pb",
        "cross_sections": {},
    }
    for generator in select_generators(config, args.generator):
        generator_report = {}
        for profile in select_profiles(config, args.profile):
            shards = sorted((base / generator / profile).glob("shard-*"))
            if not shards:
                raise RuntimeError(
                    f"No {generator} {profile} shards found under {base}."
                )
            if generator == "herwig":
                logs = {
                    herwig_integration_log(shard, args.smoke)
                    for shard in shards
                }
                estimates = [
                    parse_herwig_cross_section(path) for path in sorted(logs)
                ]
            else:
                estimates = [
                    parse_sherpa_cross_section(shard / "run.log")
                    for shard in shards
                ]
            value, uncertainty = combine_cross_sections(estimates)
            generator_report[profile] = {
                "value": value,
                "uncertainty": uncertainty,
                "independent_estimates": len(estimates),
            }
        report["cross_sections"][generator] = generator_report

    output = (
        WORK_DIR / "smoke/cross-sections.json"
        if args.smoke
        else RESULTS_DIR / "cross-sections.json"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2, sort_keys=True))


def validate_yoda(path: Path) -> None:
    try:
        import yoda
    except ImportError as error:
        raise RuntimeError(
            "The YODA Python module is required for output validation."
        ) from error

    objects = yoda.read(str(path))
    analysis_path = "/QCD_INSTANTON_KKS/"
    if f"{analysis_path}_all" not in objects:
        raise RuntimeError(f"Rivet all-events counter is missing from {path}.")
    for name, analysis_object in objects.items():
        if not name.startswith(analysis_path):
            continue
        values = []
        if hasattr(analysis_object, "bins"):
            values.extend(bin_.val() for bin_ in analysis_object.bins())
        elif hasattr(analysis_object, "val"):
            values.append(analysis_object.val())
        if any(not math.isfinite(value) for value in values):
            raise RuntimeError(
                f"Non-finite analysis value found in {path}: {name}"
            )


def command_validate(args: argparse.Namespace) -> None:
    config = load_config()
    static_checks(config)
    report = {
        "campaign": config["name"],
        "total_events": (
            sum(region["events"] for region in config["regions"].values())
            * len(config["processes"])
            * len(config["generators"])
        ),
        "sherpa_all_processes": len(sherpa_processes("all")),
        "cards": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in rendered_cards(config)
        },
    }
    if args.smoke:
        checked = []
        for generator in config["generators"]:
            for profile in profile_names(config):
                directory = (
                    WORK_DIR
                    / "smoke"
                    / generator
                    / profile
                    / "shard-00"
                )
                path = locate_yoda(directory)
                validate_yoda(path)
                checked.append(str(path.relative_to(ROOT)))
        report["validated_smoke_yoda"] = checked
    merged = RESULTS_DIR / "merged"
    if args.results:
        checked = []
        for generator in config["generators"]:
            for profile in profile_names(config):
                path = merged / f"{generator}-{profile}.yoda"
                validate_yoda(path)
                checked.append(str(path.relative_to(ROOT)))
        report["validated_yoda"] = checked
    output = RESULTS_DIR / "validation.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2, sort_keys=True))


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    render = subparsers.add_parser("render", help="render canonical cards")
    render.add_argument("--check", action="store_true")
    render.set_defaults(function=command_render)

    doctor = subparsers.add_parser("doctor", help="check local dependencies")
    doctor.set_defaults(function=command_doctor)

    build = subparsers.add_parser("build", help="build campaign plugins")
    build.add_argument(
        "--component",
        choices=("all", "herwig", "rivet", "sherpa"),
        default="all",
    )
    build.add_argument("--jobs", type=int, default=max(1, os.cpu_count() or 1))
    build.set_defaults(function=command_build)

    smoke = subparsers.add_parser("smoke", help="run small validation samples")
    smoke.add_argument("--generator", choices=("all", "herwig", "sherpa"),
                       default="all")
    smoke.add_argument("--profile", default="all")
    smoke.add_argument("--events", type=int, default=20)
    smoke.set_defaults(function=command_smoke)

    run = subparsers.add_parser("run", help="run production shards")
    run.add_argument("--generator", choices=("all", "herwig", "sherpa"),
                     default="all")
    run.add_argument("--profile", default="all")
    run.add_argument("--shards", default="all",
                     help="all or a comma-separated zero-based list")
    run.set_defaults(function=command_run)

    merge = subparsers.add_parser("merge", help="merge production YODA files")
    merge.add_argument("--generator", choices=("all", "herwig", "sherpa"),
                       default="all")
    merge.add_argument("--profile", default="all")
    merge.set_defaults(function=command_merge)

    plot = subparsers.add_parser("plot", help="plot Herwig/Sherpa overlays")
    plot.add_argument("--profile", default="all")
    plot.set_defaults(function=command_plot)

    summarize = subparsers.add_parser(
        "summarize", help="extract generator cross sections"
    )
    summarize.add_argument(
        "--generator",
        choices=("all", "herwig", "sherpa"),
        default="all",
    )
    summarize.add_argument("--profile", default="all")
    summarize.add_argument("--smoke", action="store_true")
    summarize.set_defaults(function=command_summarize)

    validate = subparsers.add_parser("validate", help="validate configuration")
    validate.add_argument("--smoke", action="store_true")
    validate.add_argument("--results", action="store_true")
    validate.set_defaults(function=command_validate)
    return parser


def main() -> int:
    parser = make_parser()
    args = parser.parse_args()
    try:
        args.function(args)
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
