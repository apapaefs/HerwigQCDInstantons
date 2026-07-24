#!/usr/bin/env python3
"""Reproducible Herwig/Sherpa QCD-instanton comparison campaign."""

from __future__ import annotations

import argparse
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
import fcntl
from functools import lru_cache
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import signal
import shlex
import shutil
import subprocess
import sys
import threading
import time
from datetime import datetime, timezone
from typing import Dict, Iterable, List, NamedTuple, Sequence, Tuple


CAMPAIGN_DIR = Path(__file__).resolve().parent
ROOT = CAMPAIGN_DIR.parents[1]
CONFIG_PATH = CAMPAIGN_DIR / "campaign.json"
BUILD_DIR = CAMPAIGN_DIR / ".build"
LOCAL_DIR = CAMPAIGN_DIR / ".local"
WORK_DIR = CAMPAIGN_DIR / ".work"
RESULTS_DIR = CAMPAIGN_DIR / "results"

PROVENANCE_SCHEMA = 2
MONITOR_SCHEMA = 1
PROGRESS_TAIL_BYTES = 131072
MASS_EDGES = tuple(
    [float(edge) for edge in range(0, 801, 10)]
    + [float(edge) for edge in range(850, 3001, 50)]
)
MIGRATION_RATIO_EDGES = tuple(index / 10.0 for index in range(101))
JET_MASS_HISTOGRAMS = (
    "jets_mreco_inclusive_eta45",
    "jets_mreco_central",
)
TRUTH_MASS_HISTOGRAM = "instanton_mass_truth"
MIGRATION_RATIO_HISTOGRAMS = (
    "jets_mreco_inclusive_over_truth",
    "jets_mreco_central_over_truth",
)

_ACTIVE_PROCESSES: set[subprocess.Popen[object]] = set()
_ACTIVE_PROCESSES_LOCK = threading.Lock()
_CANCEL_REQUESTED = threading.Event()

HERWIG_EVENT_PROGRESS = re.compile(
    r"event>\s+(?P<current>init|\d+)"
    r"(?:\s+(?P<total>\d+)|/(?P<total_alt>\d+))?"
)
HERWIG_INTEGRATION_PROGRESS = re.compile(
    r"Integrate\s+(?P<current>\d+)\s+of\s+(?P<total>\d+)"
)
SHERPA_EVENT_PROGRESS = re.compile(
    r"(?:^|\n)\s*Event\s+(?P<current>\d+)\b"
)

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


def sample_definitions(config: Dict[str, object]) -> List[Dict[str, str]]:
    return list(config["samples"])


def sample_ids(config: Dict[str, object]) -> List[str]:
    return [sample["id"] for sample in sample_definitions(config)]


def sample_definition(
    config: Dict[str, object], sample_id: str
) -> Dict[str, str]:
    for sample in sample_definitions(config):
        if sample["id"] == sample_id:
            return sample
    raise RuntimeError(f"Unknown sample: {sample_id}")


def split_profile(profile: str) -> Tuple[str, str]:
    process, region = profile.split("-", 1)
    return process, region


def stable_seed(generator: str, profile: str, shard: int) -> int:
    key = f"{generator}:{profile}:{shard}".encode("ascii")
    return int(hashlib.sha256(key).hexdigest()[:8], 16) % 900000000 + 1


def content_sha256(content: str) -> str:
    return hashlib.sha256(content.encode("utf-8")).hexdigest()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


@lru_cache(maxsize=None)
def _current_sha256(path: str, modified_ns: int, size: int) -> str:
    del modified_ns, size
    return sha256(Path(path))


def current_sha256(path: Path) -> str:
    resolved = path.resolve()
    stat = resolved.stat()
    return _current_sha256(
        str(resolved), stat.st_mtime_ns, stat.st_size
    )


@lru_cache(maxsize=1)
def source_commit() -> str:
    return capture(["git", "-C", ROOT, "rev-parse", "HEAD"])


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
    command_args = [str(item) for item in command]
    log = None
    try:
        if log_path is not None:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            log = log_path.open("w", encoding="utf-8")
            log.write(f"$ {command_text}\n")
            log.flush()
        process = subprocess.Popen(
            command_args,
            cwd=cwd,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT if log is not None else None,
            start_new_session=True,
        )
        with _ACTIVE_PROCESSES_LOCK:
            _ACTIVE_PROCESSES.add(process)
        returncode = process.wait()
        if returncode != 0:
            raise subprocess.CalledProcessError(returncode, command_args)
    finally:
        process = locals().get("process")
        if process is not None:
            with _ACTIVE_PROCESSES_LOCK:
                _ACTIVE_PROCESSES.discard(process)
        if log is not None:
            log.close()


def terminate_active_processes() -> None:
    with _ACTIVE_PROCESSES_LOCK:
        processes = list(_ACTIVE_PROCESSES)
    for process in processes:
        if process.poll() is not None:
            continue
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            continue
    deadline = time.monotonic() + 5.0
    for process in processes:
        remaining = max(0.0, deadline - time.monotonic())
        try:
            process.wait(timeout=remaining)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass


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

    text = instanton.read_text(encoding="utf-8")
    old = """    std::vector<double>     m_masses;
    void   Initialise();
    bool   DefineFlavours();"""
    new = """    std::vector<double>     m_masses;
    void   Initialise();
    bool   IncomingFlavoursActive() const;
    bool   DefineFlavours();"""
    if old in text:
        text = text.replace(old, new, 1)
    elif new not in text:
        raise RuntimeError(
            "Cannot declare the Sherpa incoming-flavour threshold guard."
        )

    old = """}

double XS_instanton::operator()(const Vec4D_Vector& momenta) {"""
    new = """}

bool XS_instanton::IncomingFlavoursActive() const {
  for (size_t i=0;i<2;i++) {
    const Flavour & flav = m_flavs[i];
    if (!flav.IsQuark()) continue;
    const kf_code code = flav.Kfcode();
    if (code>m_includeQ) return false;
    if (code==kf_b && m_Ehat<m_bthreshold) return false;
    if (code==kf_c && m_Ehat<m_cthreshold) return false;
  }
  return true;
}

double XS_instanton::operator()(const Vec4D_Vector& momenta) {"""
    if old in text:
        text = text.replace(old, new, 1)
    elif new not in text:
        raise RuntimeError(
            "Cannot define the Sherpa incoming-flavour threshold guard."
        )

    old = """  if (m_Ehat<m_Ehatmin || m_Ehat>m_Ehatmax ||
      !m_data.Interpolate(m_Ehat)) return 0.;"""
    new = """  if (m_Ehat<m_Ehatmin || m_Ehat>m_Ehatmax ||
      !IncomingFlavoursActive() ||
      !m_data.Interpolate(m_Ehat)) return 0.;"""
    if text.count(old) == 1:
        text = text.replace(old, new, 1)
    elif text.count(new) != 1:
        raise RuntimeError(
            "Cannot guard the Sherpa instanton matrix element by flavour."
        )

    old = """  if (m_Ehat<m_Ehatmin || m_Ehat>m_Ehatmax ||
      !m_data.Interpolate(m_Ehat)) return false;"""
    new = """  if (m_Ehat<m_Ehatmin || m_Ehat>m_Ehatmax ||
      !IncomingFlavoursActive() ||
      !m_data.Interpolate(m_Ehat)) return false;"""
    if text.count(old) == 1:
        text = text.replace(old, new, 1)
    elif text.count(new) != 1:
        raise RuntimeError(
            "Cannot guard the Sherpa instanton final state by flavour."
        )

    old = """    m_colours[i].resize(2);
  }
  size_t pos[2], parts[2], colindex = 500;"""
    new = """    m_colours[i].resize(2);
  }
  if (cols[0].empty() || cols[0].size()!=cols[1].size()) return false;
  size_t pos[2], parts[2], colindex = 500;"""
    if old in text:
        text = text.replace(old, new, 1)
    elif new not in text:
        raise RuntimeError(
            "Cannot apply the Sherpa instanton colour-balance guard."
        )
    instanton.write_text(text, encoding="utf-8")
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


def herwig_card_path(sample_id: str, profile: str) -> Path:
    directory = "herwig" if sample_id == "herwig" else sample_id
    return CAMPAIGN_DIR / f"cards/{directory}/{profile}.in"


def herwig_run_name(sample_id: str, profile: str) -> str:
    if sample_id == "herwig":
        return f"Campaign-Herwig-{profile}"
    return f"Campaign-{sample_id}-{profile}"


def sample_yoda_name(sample_id: str, profile: str) -> str:
    return f"{sample_id}-{profile}.yoda"


def render_dipole_kernel_adapters() -> str:
    """Keep physical zero-mode emitters in the two massive final kernels."""
    kernel_directory = "/Herwig/DipoleShower/Kernels"
    handler = "/Herwig/DipoleShower/DipoleShowerHandler"
    lines = [
        "# The hard outgoing zero modes use private physical-mass ParticleData.",
        "# These two adapters preserve that emitter record after radiation;",
        "# their splitting functions and kinematics are otherwise unchanged.",
    ]

    for stock_name in (
        "FFMdx2dgxDipoleKernel",
        "FFMux2ugxDipoleKernel",
        "FFMcx2cgxDipoleKernel",
        "FFMsx2sgxDipoleKernel",
        "FFMbx2bgxDipoleKernel",
        "FIMdx2dgxDipoleKernel",
        "FIMux2ugxDipoleKernel",
        "FIMcx2cgxDipoleKernel",
        "FIMsx2sgxDipoleKernel",
        "FIMbx2bgxDipoleKernel",
    ):
        lines.append(
            f"set {kernel_directory}/{stock_name}:UseKernel No"
        )

    for class_name, kinematics in (
        (
            "InstantonFFMqx2qgxDipoleKernel",
            "FFMassiveKinematics",
        ),
        (
            "InstantonFIMqx2qgxDipoleKernel",
            "FIMassiveKinematics",
        ),
    ):
        path = f"{kernel_directory}/{class_name}"
        lines.extend(
            [
                f"create Herwig::{class_name} {path}",
                f"set {path}:PDFRatio {kernel_directory}/PDFRatio",
                (
                    f"set {path}:SplittingKinematics "
                    f"/Herwig/DipoleShower/Kinematics/{kinematics}"
                ),
                f"set {path}:CMWScheme Factor",
                f"insert {handler}:Kernels 0 {path}",
            ]
        )
    return "\n".join(lines)


def render_herwig_variant(
    config: Dict[str, object], sample: Dict[str, str]
) -> str:
    colour = sample["colour"]
    shower = sample["shower"]
    if shower == "angular":
        if colour == "Random3":
            return ""
        return f"""
# Colour-flow systematic; all other generator settings match the baseline.
set /Herwig/MatrixElements/MEInstanton:ColourConnections {colour}
"""

    if shower != "dipole" or colour != "Random3":
        raise RuntimeError(
            f"Unsupported Herwig sample configuration: {sample['id']}"
        )

    masses = config["masses_gev"]
    mass_lines = []
    for flavour in ("u", "d", "s", "c", "b"):
        mass = masses[flavour]
        mass_lines.extend(
            [
                f"set /Herwig/Particles/{flavour}:NominalMass 0*GeV",
                f"set /Herwig/Particles/{flavour}bar:NominalMass 0*GeV",
                (
                    f"set /Herwig/Particles/{flavour}:HardProcessMass "
                    f"{mass}*GeV"
                ),
                (
                    f"set /Herwig/Particles/{flavour}bar:HardProcessMass "
                    f"{mass}*GeV"
                ),
            ]
        )
    kernel_adapters = render_dipole_kernel_adapters()

    return f"""
# Full Herwig dipole-shower variation. The native tune also changes shower
# alpha_s, infrared cutoffs, constituent masses, and hadronization parameters.
set /Herwig/EventHandlers/EventHandler:CascadeHandler /Herwig/DipoleShower/DipoleShowerHandler
read snippets/Dipole_AutoTunes_gss.in

# Restore the settings held fixed in this Herwig/Sherpa comparison.
set /Herwig/DipoleShower/DipoleShowerHandler:PDFA /Herwig/Partons/InstantonPDF
set /Herwig/DipoleShower/DipoleShowerHandler:PDFB /Herwig/Partons/InstantonPDF
set /Herwig/DipoleShower/DipoleShowerHandler:MPIHandler NULL
# MPI-on alternative:
#set /Herwig/DipoleShower/DipoleShowerHandler:MPIHandler /Herwig/UnderlyingEvent/MPIHandler
{chr(10).join(mass_lines)}
# Five-flavour ISR uses massless canonical quark data. MEInstanton detects
# the distinct HardProcessMass values and makes private physical-mass data
# for outgoing zero modes, so their MAMBO masses remain matched to Sherpa.
{kernel_adapters}
set /Herwig/MatrixElements/MEInstanton:ColourConnections Random3
"""


def render_herwig_profile(
    config: Dict[str, object],
    sample: Dict[str, str],
    process: str,
    region: str,
    *,
    smoke: bool = False,
) -> str:
    values = config["regions"][region]
    profile = f"{process}-{region}"
    sample_id = sample["id"]
    process_option = "GG" if process == "gg" else "All"
    variant = render_herwig_variant(config, sample)
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
{variant}
set /Herwig/Generators/EventGenerator:NumberOfEvents {values['events']}
set /Herwig/Generators/EventGenerator:RandomNumberGenerator:Seed {stable_seed(sample_id, profile, 0)}
set /Herwig/MatrixElements/MEInstanton:Processes {process_option}
set /Herwig/Cuts/InstantonCuts:MHatMin {values['min_mass_gev']}*GeV
set /Herwig/Cuts/InstantonCuts:MHatMax {values['max_mass_gev']}*GeV
set /Herwig/Analysis/Rivet:Filename {sample_yoda_name(sample_id, profile)}
{smoke_sampler}
cd /Herwig/Generators
saverun {herwig_run_name(sample_id, profile)} EventGenerator
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
    for sample in sample_definitions(config):
        for process in config["processes"]:
            for region in config["regions"]:
                profile = f"{process}-{region}"
                if sample["engine"] == "herwig":
                    result[
                        herwig_card_path(sample["id"], profile)
                    ] = render_herwig_profile(
                        config, sample, process, region
                    )
                elif sample["engine"] == "sherpa":
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
    samples = sample_definitions(config)
    expected_samples = [
        "herwig",
        "herwig-qcdinsplanar",
        "herwig-random3-dipole",
        "sherpa",
    ]
    if sample_ids(config) != expected_samples:
        raise RuntimeError(
            f"Campaign samples must be ordered as {expected_samples}."
        )
    if len(set(sample_ids(config))) != len(samples):
        raise RuntimeError("Campaign sample IDs must be unique.")
    for sample in samples:
        if sample["engine"] not in ("herwig", "sherpa"):
            raise RuntimeError(f"Unknown sample engine: {sample['engine']}")
    if sample_definition(config, "herwig-qcdinsplanar") != {
        "id": "herwig-qcdinsplanar",
        "engine": "herwig",
        "colour": "QCDINSPlanar",
        "shower": "angular",
        "title": "Herwig QCDINSPlanar (angular)",
    }:
        raise RuntimeError("The QCDINSPlanar sample definition is inconsistent.")
    dipole = sample_definition(config, "herwig-random3-dipole")
    if dipole["colour"] != "Random3" or dipole["shower"] != "dipole":
        raise RuntimeError("The dipole sample must use Random3.")

    all_processes = sherpa_processes("all")
    if len(all_processes) != 56 or len(set(all_processes)) != 56:
        raise RuntimeError("Sherpa All must contain 56 unique incoming states.")
    for first, second in all_processes:
        if first == second and abs(first) <= 5:
            raise RuntimeError("Equal-flavour qq or qbar-qbar process found.")

    total_events = (
        sum(values["events"] for values in config["regions"].values())
        * len(config["processes"])
        * len(samples)
    )
    if total_events != 560000:
        raise RuntimeError(f"Unexpected campaign event total: {total_events}")
    if (
        MASS_EDGES[:81] != tuple(float(x) for x in range(0, 801, 10))
        or MASS_EDGES[81:]
        != tuple(float(x) for x in range(850, 3001, 50))
    ):
        raise RuntimeError("The instanton-mass binning is inconsistent.")
    if MIGRATION_RATIO_EDGES != tuple(
        index / 10.0 for index in range(101)
    ):
        raise RuntimeError("The mass-migration binning is inconsistent.")
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
        ROOT / "InstantonDipoleKernels.cc",
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


def plotting_environment() -> Dict[str, str]:
    env = campaign_environment()
    matplotlib_dir = WORK_DIR / "matplotlib"
    matplotlib_dir.mkdir(parents=True, exist_ok=True)
    env["MPLCONFIGDIR"] = str(matplotlib_dir)
    for variable in ("RIVET_PLOT_PATH", "RIVET_INFO_PATH"):
        paths = [str(ROOT / "Rivet")]
        if env.get(variable):
            paths.append(env[variable])
        env[variable] = os.pathsep.join(paths)

    # Some local activation scripts export a Rivet TEXMFCNF directory that
    # does not exist. Let the system TeX installation discover its own files.
    texmfcnf = env.get("TEXMFCNF")
    if texmfcnf and not Path(texmfcnf).exists():
        for name in ("TEXMFCNF", "TEXMFHOME", "TEXINPUTS"):
            env.pop(name, None)
    return env


class CampaignTask(NamedTuple):
    sample_id: str
    profile: str
    shard: int
    events: int
    output_root: Path
    smoke: bool

    @property
    def label(self) -> str:
        return f"{self.sample_id}/{self.profile}/shard-{self.shard:02d}"

    @property
    def directory(self) -> Path:
        return (
            self.output_root
            / self.sample_id
            / self.profile
            / f"shard-{self.shard:02d}"
        )


class TaskActivity(NamedTuple):
    task: CampaignTask
    engine: str
    started_at: float


def positive_integer(value: str) -> int:
    result = int(value)
    if result < 1:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return result


def nonnegative_integer(value: str) -> int:
    result = int(value)
    if result < 0:
        raise argparse.ArgumentTypeError("value must be a non-negative integer")
    return result


def finite_seconds(value: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise argparse.ArgumentTypeError("value must be finite")
    return result


def relative_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except ValueError:
        return str(path.resolve())


def write_json_atomic(path: Path, data: Dict[str, object]) -> None:
    write_text_atomic(
        path, json.dumps(data, indent=2, sort_keys=True) + "\n"
    )


def write_text_atomic(path: Path, content: str) -> None:
    temporary = path.with_name(
        f".{path.name}.{os.getpid()}.{threading.get_ident()}.tmp"
    )
    try:
        temporary.write_text(content, encoding="utf-8")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def campaign_monitor_dir(output_root: Path) -> Path:
    if output_root.name in ("shards", "smoke"):
        return output_root.parent / "monitor"
    return output_root / "monitor"


def read_log_tail(path: Path, limit: int = PROGRESS_TAIL_BYTES) -> str:
    try:
        with path.open("rb") as stream:
            stream.seek(0, os.SEEK_END)
            size = stream.tell()
            stream.seek(max(0, size - limit))
            return stream.read().decode("utf-8", errors="replace")
    except OSError:
        return ""


def log_modified_since(path: Path, started_at: float) -> bool:
    try:
        return path.is_file() and path.stat().st_mtime >= started_at
    except OSError:
        return False


def parse_herwig_event_progress(
    text: str, expected_events: int
) -> Tuple[int, int] | None:
    matches = list(HERWIG_EVENT_PROGRESS.finditer(text))
    if not matches:
        return None
    match = matches[-1]
    current_text = match.group("current")
    current = 0 if current_text == "init" else int(current_text)
    total_text = match.group("total") or match.group("total_alt")
    total = int(total_text) if total_text else expected_events
    return min(max(current, 0), total), total


def parse_sherpa_event_progress(
    text: str, expected_events: int
) -> Tuple[int, int] | None:
    matches = list(SHERPA_EVENT_PROGRESS.finditer(text))
    if not matches:
        return None
    current = int(matches[-1].group("current"))
    return min(max(current, 0), expected_events), expected_events


def parse_herwig_integration_progress(
    text: str,
) -> Tuple[int, int] | None:
    matches = list(HERWIG_INTEGRATION_PROGRESS.finditer(text))
    if not matches:
        return None
    match = matches[-1]
    return int(match.group("current")), int(match.group("total"))


def active_herwig_integration_log(
    activity: TaskActivity,
) -> Path | None:
    task = activity.task
    local = task.directory / "read.log"
    if log_modified_since(local, activity.started_at):
        return local

    if task.sample_id == "herwig":
        root = WORK_DIR / "herwig-runs" / task.profile
    else:
        root = WORK_DIR / "herwig-runs" / task.sample_id / task.profile
    if not root.is_dir():
        return None
    candidates = []
    for path in root.glob("*/read.log"):
        try:
            modified = path.stat().st_mtime
        except OSError:
            continue
        if modified >= activity.started_at:
            candidates.append((modified, path))
    if not candidates:
        return None
    return max(candidates, key=lambda item: item[0])[1]


def task_progress(activity: TaskActivity) -> Tuple[str, int]:
    task = activity.task
    run_log = task.directory / "run.log"
    run_text = (
        read_log_tail(run_log)
        if log_modified_since(run_log, activity.started_at)
        else ""
    )
    if activity.engine == "herwig":
        progress = parse_herwig_event_progress(run_text, task.events)
    else:
        progress = parse_sherpa_event_progress(run_text, task.events)
    if progress is not None:
        current, total = progress
        percent = 100.0 * current / total if total else 0.0
        return f"{current}/{total} ({percent:5.1f}%)", current

    if activity.engine == "herwig":
        integration_log = active_herwig_integration_log(activity)
        integration = parse_herwig_integration_progress(
            read_log_tail(integration_log)
            if integration_log is not None
            else ""
        )
        if integration is not None:
            current, total = integration
            return f"integrate {current}/{total}", 0
    return "setup", 0


def format_duration(seconds: float) -> str:
    total = max(0, int(seconds))
    hours, remainder = divmod(total, 3600)
    minutes, seconds = divmod(remainder, 60)
    if hours:
        return f"{hours:d}:{minutes:02d}:{seconds:02d}"
    return f"{minutes:02d}:{seconds:02d}"


def format_count(value: int) -> str:
    return f"{value:,}"


def render_table(
    headers: Sequence[str],
    rows: Sequence[Sequence[object]],
    *,
    right_aligned: Iterable[int] = (),
) -> str:
    rendered = [[str(value) for value in row] for row in rows]
    widths = [
        max(
            len(headers[index]),
            *(len(row[index]) for row in rendered),
        )
        for index in range(len(headers))
    ]
    right = set(right_aligned)

    def render_row(row: Sequence[str]) -> str:
        cells = []
        for index, value in enumerate(row):
            if index in right:
                cells.append(value.rjust(widths[index]))
            else:
                cells.append(value.ljust(widths[index]))
        return "  ".join(cells).rstrip()

    separator = ["-" * width for width in widths]
    return "\n".join(
        [
            render_row(list(headers)),
            render_row(separator),
            *(render_row(row) for row in rendered),
        ]
    )


def build_campaign_monitor_payload(
    *,
    config: Dict[str, object],
    tasks: Sequence[CampaignTask],
    outcomes: Dict[str, str],
    active: Dict[str, TaskActivity],
    run_id: str,
    phase: str,
    started_at: float,
    workers: int,
    message: str = "",
) -> Dict[str, object]:
    now = time.time()
    groups: Dict[Tuple[str, str], Dict[str, object]] = {}
    counts = {
        "completed": 0,
        "running": 0,
        "pending": 0,
        "failed": 0,
        "skipped": 0,
        "total": len(tasks),
    }
    completed_events = 0

    for task in tasks:
        key = (task.sample_id, task.profile)
        group = groups.setdefault(
            key,
            {
                "sample": task.sample_id,
                "profile": task.profile,
                "completed": 0,
                "running": 0,
                "pending": 0,
                "failed": 0,
                "skipped": 0,
                "total": 0,
                "events": 0,
            },
        )
        group["total"] += 1
        group["events"] += task.events
        outcome = outcomes.get(task.label)
        if outcome == "DONE":
            state = "completed"
            completed_events += task.events
        elif outcome == "SKIP":
            state = "skipped"
            completed_events += task.events
        elif outcome == "FAIL":
            state = "failed"
        elif task.label in active:
            state = "running"
        else:
            state = "pending"
        counts[state] += 1
        group[state] += 1

    active_rows = []
    active_generated = 0
    for label in sorted(active):
        activity = active[label]
        progress, generated = task_progress(activity)
        active_generated += generated
        task = activity.task
        active_rows.append(
            {
                "sample": task.sample_id,
                "profile": task.profile,
                "shard": task.shard,
                "engine": activity.engine,
                "progress": progress,
                "current_events": generated,
                "events": task.events,
                "seed": stable_seed(
                    task.sample_id, task.profile, task.shard
                ),
                "runtime_s": round(now - activity.started_at, 1),
            }
        )

    mode = "smoke" if tasks and tasks[0].smoke else "production"
    return {
        "schema_version": MONITOR_SCHEMA,
        "campaign": config["name"],
        "mode": mode,
        "run_id": run_id,
        "phase": phase,
        "started_at": datetime.fromtimestamp(
            started_at, timezone.utc
        ).isoformat(),
        "updated_at": datetime.fromtimestamp(now, timezone.utc).isoformat(),
        "elapsed_s": round(now - started_at, 1),
        "workers": workers,
        "selection": {
            "samples": list(dict.fromkeys(task.sample_id for task in tasks)),
            "profiles": list(dict.fromkeys(task.profile for task in tasks)),
            "shards": sorted({task.shard for task in tasks}),
            "output_root": relative_path(tasks[0].output_root)
            if tasks
            else "",
        },
        "shards": counts,
        "events": {
            "completed": completed_events,
            "active_generated": active_generated,
            "target": sum(task.events for task in tasks),
        },
        "groups": list(groups.values()),
        "active": active_rows,
        "message": message,
    }


def render_campaign_monitor(
    payload: Dict[str, object], max_listed: int = 12
) -> str:
    shards = payload["shards"]
    events = payload["events"]
    selection = payload["selection"]
    lines = [
        f"Campaign: {payload['campaign']}",
        (
            f"Mode: {payload['mode']}  Run: {payload['run_id']}  "
            f"Phase: {payload['phase']}"
        ),
        (
            f"Elapsed: {format_duration(payload['elapsed_s'])}  "
            f"Workers: {payload['workers']}"
        ),
        (
            f"Selection: samples={','.join(selection['samples'])}  "
            f"profiles={','.join(selection['profiles'])}  "
            f"shards={len(selection['shards'])}"
        ),
        "",
        (
            "Shards: "
            f"{shards['completed']} completed, "
            f"{shards['skipped']} skipped, "
            f"{shards['running']} running, "
            f"{shards['pending']} pending, "
            f"{shards['failed']} failed / {shards['total']} total"
        ),
        (
            "Events: "
            f"{format_count(events['completed'])} completed"
            f" + {format_count(events['active_generated'])} active"
            f" / {format_count(events['target'])} target"
        ),
    ]
    if payload.get("message"):
        lines.extend(["", f"Message: {payload['message']}"])

    incomplete_groups = [
        group
        for group in payload["groups"]
        if (
            group["completed"]
            + group["skipped"]
            < group["total"]
            or group["failed"]
        )
    ]
    if incomplete_groups:
        group_rows = [
            (
                f"{group['sample']}/{group['profile']}",
                group["completed"],
                group["skipped"],
                group["running"],
                group["pending"],
                group["failed"],
                group["total"],
            )
            for group in incomplete_groups[:max_listed]
        ]
        lines.extend(
            [
                "",
                "Incomplete groups:",
                render_table(
                    (
                        "Run",
                        "Done",
                        "Skip",
                        "Run",
                        "Wait",
                        "Fail",
                        "Total",
                    ),
                    group_rows,
                    right_aligned=range(1, 7),
                ),
            ]
        )
        if len(incomplete_groups) > max_listed:
            lines.append(
                f"... {len(incomplete_groups) - max_listed} more groups"
            )

    active_rows = payload["active"]
    if active_rows:
        rows = [
            (
                f"{row['sample']}/{row['profile']}",
                row["shard"],
                row["engine"],
                row["progress"],
                row["seed"],
                format_duration(row["runtime_s"]),
            )
            for row in active_rows[:max_listed]
        ]
        lines.extend(
            [
                "",
                "Active shards:",
                render_table(
                    ("Run", "Shard", "Engine", "Progress", "Seed", "Runtime"),
                    rows,
                    right_aligned=(1, 4, 5),
                ),
            ]
        )
        if len(active_rows) > max_listed:
            lines.append(f"... {len(active_rows) - max_listed} more shards")
    return "\n".join(lines) + "\n"


def write_campaign_monitor(
    output_root: Path,
    payload: Dict[str, object],
    *,
    max_listed: int,
) -> str:
    directory = campaign_monitor_dir(output_root)
    directory.mkdir(parents=True, exist_ok=True)
    rendered = render_campaign_monitor(payload, max_listed=max_listed)
    write_json_atomic(directory / "status.json", payload)
    write_text_atomic(directory / "status.txt", rendered)
    return rendered


def emit_campaign_monitor(rendered: str) -> None:
    if not sys.stdout.isatty():
        return
    print("\033[2J\033[H", end="")
    print(rendered, end="", flush=True)


def require_fresh_artifact(
    output: Path, sources: Sequence[Path], description: str
) -> None:
    if not output.is_file():
        raise RuntimeError(f"Build the {description} before running.")
    newest_source = max(source.stat().st_mtime_ns for source in sources)
    if output.stat().st_mtime_ns < newest_source:
        raise RuntimeError(
            f"The {description} is older than its source; rebuild it."
        )


def ensure_runtime_artifacts(
    config: Dict[str, object], sample: Dict[str, str]
) -> None:
    require_fresh_artifact(
        BUILD_DIR / "rivet/RivetQCD_INSTANTON_KKS.so",
        [ROOT / "Rivet/QCD_INSTANTON_KKS.cc"],
        "campaign Rivet plugin",
    )
    if sample["engine"] == "herwig":
        require_fresh_artifact(
            BUILD_DIR / "herwig/CampaignInstantons.so",
            [
                ROOT / "MEInstanton.cc",
                ROOT / "MEInstanton.h",
                ROOT / "InstantonDipoleKernels.cc",
                ROOT / "InstantonDipoleKernels.h",
            ],
            "campaign Herwig plugin",
        )
        find_executable("Herwig", "HERWIG_BIN")
    else:
        find_campaign_sherpa_binary(config)


@lru_cache(maxsize=None)
def executable_metadata(executable: Path) -> Dict[str, str]:
    resolved = executable.resolve()
    version = capture([resolved, "--version"], allow_nonzero=True).splitlines()
    return {
        "generator_executable": str(resolved),
        "generator_executable_sha256": current_sha256(resolved),
        "generator_version": version[0] if version else "unknown",
    }


def provenance_record(
    sample: Dict[str, str],
    task: CampaignTask,
    card: Path,
    command: Sequence[object],
    *,
    run_file: Path | None = None,
    status: str,
) -> Dict[str, object]:
    executable = Path(command[0]).resolve()
    record: Dict[str, object] = {
        "schema_version": PROVENANCE_SCHEMA,
        "status": status,
        "sample": sample["id"],
        "engine": sample["engine"],
        "generator": sample["engine"],
        "profile": task.profile,
        "shard": task.shard,
        "events": task.events,
        "seed": stable_seed(sample["id"], task.profile, task.shard),
        "card": relative_path(card),
        "card_sha256": current_sha256(card),
        "command": [str(item) for item in command],
        "git_commit": source_commit(),
        "rivet_source_sha256": current_sha256(
            ROOT / "Rivet/QCD_INSTANTON_KKS.cc"
        ),
        "rivet_plugin_sha256": current_sha256(
            BUILD_DIR / "rivet/RivetQCD_INSTANTON_KKS.so"
        ),
        "yoda": sample_yoda_name(sample["id"], task.profile),
        "started_at": datetime.now(timezone.utc).isoformat(),
    }
    record.update(executable_metadata(executable))
    if sample["engine"] == "herwig":
        record["instanton_source_sha256"] = current_sha256(
            ROOT / "MEInstanton.cc"
        )
        record["instanton_header_sha256"] = current_sha256(
            ROOT / "MEInstanton.h"
        )
        record["dipole_adapter_source_sha256"] = current_sha256(
            ROOT / "InstantonDipoleKernels.cc"
        )
        record["dipole_adapter_header_sha256"] = current_sha256(
            ROOT / "InstantonDipoleKernels.h"
        )
        record["instanton_plugin_sha256"] = current_sha256(
            BUILD_DIR / "herwig/CampaignInstantons.so"
        )
    if run_file is not None:
        record["herwig_run"] = str(run_file.resolve())
        record["herwig_run_sha256"] = current_sha256(run_file)
    return record


PROVENANCE_MATCH_FIELDS = (
    "schema_version",
    "sample",
    "engine",
    "profile",
    "shard",
    "events",
    "seed",
    "card",
    "card_sha256",
    "command",
    "generator_executable",
    "generator_executable_sha256",
    "generator_version",
    "rivet_source_sha256",
    "rivet_plugin_sha256",
    "yoda",
    "instanton_source_sha256",
    "instanton_header_sha256",
    "dipole_adapter_source_sha256",
    "dipole_adapter_header_sha256",
    "instanton_plugin_sha256",
    "herwig_run",
    "herwig_run_sha256",
)


def provenance_mismatches(
    stored: Dict[str, object], expected: Dict[str, object]
) -> List[str]:
    return [
        field
        for field in PROVENANCE_MATCH_FIELDS
        if stored.get(field) != expected.get(field)
    ]


def begin_provenance(
    directory: Path, record: Dict[str, object]
) -> Path:
    pending = directory / "provenance.pending.json"
    write_json_atomic(pending, record)
    return pending


def complete_provenance(
    directory: Path, record: Dict[str, object]
) -> None:
    completed = dict(record)
    completed["status"] = "complete"
    completed["completed_at"] = datetime.now(timezone.utc).isoformat()
    write_json_atomic(directory / "provenance.json", completed)
    (directory / "provenance.pending.json").unlink(missing_ok=True)


def prepare_herwig_run(
    sample_id: str,
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
    if sample_id == "herwig":
        directory = WORK_DIR / "herwig-runs" / profile / fingerprint[:16]
    else:
        directory = (
            WORK_DIR
            / "herwig-runs"
            / sample_id
            / profile
            / fingerprint[:16]
        )
    directory.mkdir(parents=True, exist_ok=True)
    run_file = directory / f"{herwig_run_name(sample_id, profile)}.run"

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
    sample: Dict[str, str],
    task: CampaignTask,
) -> None:
    profile = task.profile
    directory = task.directory
    canonical_card = herwig_card_path(sample["id"], profile)
    herwig = find_executable("Herwig", "HERWIG_BIN")
    plugin_path = BUILD_DIR / "herwig"
    env = campaign_environment()
    directory.mkdir(parents=True, exist_ok=True)
    card = canonical_card
    if task.smoke:
        process, region = split_profile(profile)
        card = directory / f"{sample['id']}-{profile}-smoke.in"
        card.write_text(
            render_herwig_profile(
                config,
                sample,
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
        run_file = directory / f"{herwig_run_name(sample['id'], profile)}.run"
    else:
        run_file = prepare_herwig_run(
            sample["id"],
            profile,
            card,
            herwig,
            plugin_path,
            env,
        )
    seed = stable_seed(sample["id"], profile, task.shard)
    run_command = [
        herwig,
        "run",
        "-L",
        plugin_path,
        "-N",
        task.events,
        "-s",
        seed,
        run_file,
    ]
    record = provenance_record(
        sample,
        task,
        card,
        run_command,
        run_file=run_file,
        status="running",
    )
    begin_provenance(directory, record)
    run_checked(
        run_command,
        cwd=directory,
        env=env,
        log_path=directory / "run.log",
    )
    validate_yoda(locate_yoda(directory), expected_events=task.events)
    complete_provenance(directory, record)


def run_sherpa(
    config: Dict[str, object],
    sample: Dict[str, str],
    task: CampaignTask,
) -> None:
    profile = task.profile
    directory = task.directory
    card = CAMPAIGN_DIR / f"cards/sherpa/{profile}.yaml"
    sherpa = find_campaign_sherpa_binary(config)
    directory.mkdir(parents=True, exist_ok=True)
    seed = stable_seed(sample["id"], profile, task.shard)
    run_command = [
        sherpa,
        "-e",
        task.events,
        "-R",
        seed,
        "-r",
        "Results",
        "-A",
        f"{sample['id']}-{profile}",
        card,
    ]
    record = provenance_record(
        sample,
        task,
        card,
        run_command,
        status="running",
    )
    begin_provenance(directory, record)
    run_checked(
        run_command,
        cwd=directory,
        env=campaign_environment(),
        log_path=directory / "run.log",
    )
    validate_yoda(locate_yoda(directory), expected_events=task.events)
    complete_provenance(directory, record)


def select_profiles(config: Dict[str, object], selection: str) -> List[str]:
    available = profile_names(config)
    if selection == "all":
        return available
    if selection not in available:
        raise RuntimeError(f"Unknown profile {selection}; choose {available}.")
    return [selection]


def select_samples(
    config: Dict[str, object], selection: str
) -> List[Dict[str, str]]:
    if selection == "all":
        return sample_definitions(config)
    requested = [item.strip() for item in selection.split(",") if item.strip()]
    if not requested:
        raise RuntimeError("The sample selection is empty.")
    if "all" in requested:
        raise RuntimeError("Use 'all' by itself in a sample selection.")
    if len(requested) != len(set(requested)):
        raise RuntimeError("A sample may only be selected once.")
    available = sample_ids(config)
    unknown = [sample for sample in requested if sample not in available]
    if unknown:
        raise RuntimeError(
            f"Unknown sample(s) {unknown}; choose from {available}."
        )
    return [sample_definition(config, sample) for sample in requested]


def selected_samples(
    config: Dict[str, object], args: argparse.Namespace
) -> List[Dict[str, str]]:
    sample_selection = getattr(args, "sample", None)
    generator_selection = getattr(args, "generator", None)
    if sample_selection is not None and generator_selection is not None:
        raise RuntimeError("Use either --sample or --generator, not both.")
    selection = sample_selection or generator_selection or "all"
    return select_samples(config, selection)


def card_for_task(
    config: Dict[str, object],
    sample: Dict[str, str],
    task: CampaignTask,
) -> Path:
    if sample["engine"] == "sherpa":
        return CAMPAIGN_DIR / f"cards/sherpa/{task.profile}.yaml"
    if not task.smoke:
        return herwig_card_path(sample["id"], task.profile)
    return task.directory / f"{sample['id']}-{task.profile}-smoke.in"


def desired_smoke_card(
    config: Dict[str, object],
    sample: Dict[str, str],
    task: CampaignTask,
) -> str:
    process, region = split_profile(task.profile)
    return render_herwig_profile(
        config, sample, process, region, smoke=True
    )


def expected_existing_provenance(
    config: Dict[str, object],
    sample: Dict[str, str],
    task: CampaignTask,
    stored: Dict[str, object],
) -> Dict[str, object]:
    card = card_for_task(config, sample, task)
    if task.smoke and sample["engine"] == "herwig":
        desired = desired_smoke_card(config, sample, task)
        if not card.is_file() or card.read_text(encoding="utf-8") != desired:
            raise RuntimeError("smoke card")

    seed = stable_seed(sample["id"], task.profile, task.shard)
    if sample["engine"] == "herwig":
        run_path = stored.get("herwig_run")
        if run_path is None:
            command = stored.get("command", [])
            run_path = command[-1] if command else None
        if run_path is None:
            raise RuntimeError("Herwig run file")
        run_file = Path(str(run_path))
        if not run_file.is_file():
            raise RuntimeError("Herwig run file")
        command = [
            find_executable("Herwig", "HERWIG_BIN"),
            "run",
            "-L",
            BUILD_DIR / "herwig",
            "-N",
            task.events,
            "-s",
            seed,
            run_file,
        ]
        return provenance_record(
            sample,
            task,
            card,
            command,
            run_file=run_file,
            status="running",
        )

    command = [
        find_campaign_sherpa_binary(config),
        "-e",
        task.events,
        "-R",
        seed,
        "-r",
        "Results",
        "-A",
        f"{sample['id']}-{task.profile}",
        card,
    ]
    return provenance_record(
        sample,
        task,
        card,
        command,
        status="running",
    )


def inspect_shard(
    config: Dict[str, object],
    sample: Dict[str, str],
    task: CampaignTask,
) -> Tuple[str, str]:
    directory = task.directory
    if not directory.exists():
        return "missing", ""
    if not any(directory.iterdir()):
        return "missing", ""

    complete_path = directory / "provenance.json"
    pending_path = directory / "provenance.pending.json"
    provenance_path = complete_path if complete_path.is_file() else pending_path
    if not provenance_path.is_file():
        return "stale", "missing provenance"
    try:
        stored = json.loads(provenance_path.read_text(encoding="utf-8"))
        expected = expected_existing_provenance(
            config, sample, task, stored
        )
    except (json.JSONDecodeError, OSError, RuntimeError) as error:
        return "stale", str(error)

    mismatches = provenance_mismatches(stored, expected)
    if mismatches:
        return "stale", ", ".join(mismatches)

    try:
        yoda_path = locate_yoda(directory)
        expected_name = sample_yoda_name(sample["id"], task.profile)
        if yoda_path.name not in (expected_name, expected_name + ".gz"):
            raise RuntimeError(f"unexpected YODA filename {yoda_path.name}")
        validate_yoda(yoda_path, expected_events=task.events)
    except RuntimeError as error:
        return "incomplete", str(error)

    if provenance_path == pending_path:
        complete_provenance(directory, stored)
    elif stored.get("status") != "complete":
        return "incomplete", "provenance status is not complete"
    return "complete", ""


def archive_shard(
    task: CampaignTask, category: str, run_id: str
) -> Path:
    relative = task.directory.relative_to(task.output_root)
    target = task.output_root.parent / category / run_id / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    suffix = 1
    original = target
    while target.exists():
        target = original.with_name(f"{original.name}-{suffix}")
        suffix += 1
    shutil.move(str(task.directory), str(target))
    return target


def shard_lock_path(task: CampaignTask) -> Path:
    return (
        task.output_root.parent
        / ".locks"
        / task.output_root.name
        / task.sample_id
        / task.profile
        / f"shard-{task.shard:02d}.lock"
    )


def execute_task(
    config: Dict[str, object],
    task: CampaignTask,
    *,
    force: bool,
    run_id: str,
) -> str:
    sample = sample_definition(config, task.sample_id)
    lock_path = shard_lock_path(task)
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("w", encoding="utf-8") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if _CANCEL_REQUESTED.is_set():
            raise RuntimeError("cancelled")

        state, detail = inspect_shard(config, sample, task)
        if state == "complete" and not force:
            return "SKIP"
        if state == "stale" and not force:
            raise RuntimeError(
                f"stale output ({detail}); use --force to archive and rerun"
            )
        if task.directory.exists() and any(task.directory.iterdir()):
            category = "superseded" if force or state == "stale" else "attempts"
            archived = archive_shard(task, category, run_id)
            print(f"ARCHIVE {task.label} -> {archived}", flush=True)

        if _CANCEL_REQUESTED.is_set():
            raise RuntimeError("cancelled")
        print(f"START {task.label}", flush=True)
        if sample["engine"] == "herwig":
            run_herwig(config, sample, task)
        else:
            run_sherpa(config, sample, task)
        return "DONE"


def run_selection(
    *,
    config: Dict[str, object],
    samples: Sequence[Dict[str, str]],
    profile_selection: str,
    shards: Iterable[int],
    event_override: int | None,
    output_root: Path,
    jobs: int,
    force: bool,
    progress_interval: float = 5.0,
    max_listed: int = 12,
) -> None:
    static_checks(config)
    for sample in samples:
        ensure_runtime_artifacts(config, sample)

    tasks: List[CampaignTask] = []
    profiles = select_profiles(config, profile_selection)
    selected_shards = list(shards)
    # Put distinct sample/profile integrations at the front of the queue.
    # Otherwise the first workers can all wait behind one integration lock.
    for shard in selected_shards:
        for sample in samples:
            for profile in profiles:
                _, region = split_profile(profile)
                production_events = config["regions"][region]["events"]
                per_shard, remainder = divmod(
                    production_events, config["shards"]
                )
                events = (
                    event_override
                    if event_override is not None
                    else per_shard + int(shard < remainder)
                )
                tasks.append(
                    CampaignTask(
                        sample["id"],
                        profile,
                        shard,
                        events,
                        output_root,
                        event_override is not None,
                    )
                )

    run_id = (
        datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        + f"-{os.getpid()}"
    )
    started_at = time.time()
    failures: List[Tuple[str, BaseException]] = []
    outcomes: Dict[str, str] = {}
    active: Dict[str, TaskActivity] = {}
    active_lock = threading.Lock()
    _CANCEL_REQUESTED.clear()

    def tracked_execute(task: CampaignTask) -> str:
        sample = sample_definition(config, task.sample_id)
        with active_lock:
            active[task.label] = TaskActivity(
                task, sample["engine"], time.time()
            )
        try:
            return execute_task(
                config,
                task,
                force=force,
                run_id=run_id,
            )
        finally:
            with active_lock:
                active.pop(task.label, None)

    def update_monitor(phase: str, message: str = "") -> None:
        with active_lock:
            active_snapshot = dict(active)
        payload = build_campaign_monitor_payload(
            config=config,
            tasks=tasks,
            outcomes=dict(outcomes),
            active=active_snapshot,
            run_id=run_id,
            phase=phase,
            started_at=started_at,
            workers=jobs,
            message=message,
        )
        rendered = write_campaign_monitor(
            output_root, payload, max_listed=max_listed
        )
        emit_campaign_monitor(rendered)

    running_phase = (
        "running-smoke" if event_override is not None else "running-production"
    )
    update_monitor(running_phase)
    executor = ThreadPoolExecutor(max_workers=jobs)
    futures: Dict[Future[str], CampaignTask] = {
        executor.submit(tracked_execute, task): task
        for task in tasks
    }
    completed = 0
    remaining = set(futures)
    try:
        while remaining:
            timeout = (
                None
                if progress_interval < 0
                else max(0.1, progress_interval)
            )
            finished, remaining = wait(
                remaining,
                timeout=timeout,
                return_when=FIRST_COMPLETED,
            )
            for future in finished:
                task = futures[future]
                completed += 1
                try:
                    status = future.result()
                    outcomes[task.label] = status
                    print(
                        f"[{completed}/{len(tasks)}] {status} {task.label}",
                        flush=True,
                    )
                except BaseException as error:
                    outcomes[task.label] = "FAIL"
                    failures.append((task.label, error))
                    print(
                        f"[{completed}/{len(tasks)}] "
                        f"FAIL {task.label}: {error}",
                        file=sys.stderr,
                        flush=True,
                    )
            if progress_interval >= 0:
                update_monitor(running_phase)
    except KeyboardInterrupt:
        _CANCEL_REQUESTED.set()
        for future in futures:
            future.cancel()
        terminate_active_processes()
        executor.shutdown(wait=True, cancel_futures=True)
        update_monitor("campaign-interrupted", "Interrupted by user.")
        raise
    else:
        executor.shutdown(wait=True)

    if failures:
        details = "; ".join(
            f"{label}: {error}" for label, error in failures[:10]
        )
        if len(failures) > 10:
            details += f"; and {len(failures) - 10} more"
        update_monitor("campaign-failed", details)
        raise RuntimeError(f"{len(failures)} campaign task(s) failed: {details}")
    update_monitor("campaign-complete")


def command_smoke(args: argparse.Namespace) -> None:
    config = load_config()
    run_selection(
        config=config,
        samples=selected_samples(config, args),
        profile_selection=args.profile,
        shards=[0],
        event_override=args.events,
        output_root=WORK_DIR / "smoke",
        jobs=args.jobs,
        force=args.force,
        progress_interval=args.progress_interval,
        max_listed=args.max_listed,
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
        samples=selected_samples(config, args),
        profile_selection=args.profile,
        shards=parse_shards(args.shards, config["shards"]),
        event_override=None,
        output_root=RESULTS_DIR / "shards",
        jobs=args.jobs,
        force=args.force,
        progress_interval=args.progress_interval,
        max_listed=args.max_listed,
    )


def command_status(args: argparse.Namespace) -> None:
    output_root = (
        WORK_DIR / "smoke" if args.smoke else RESULTS_DIR / "shards"
    )
    directory = campaign_monitor_dir(output_root)
    path = directory / ("status.json" if args.json else "status.txt")
    if not path.is_file():
        mode = "smoke" if args.smoke else "production"
        raise RuntimeError(f"No {mode} campaign tracker exists at {path}.")
    print(path.read_text(encoding="utf-8"), end="")


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
    input_root = (
        WORK_DIR / "smoke" if args.smoke else RESULTS_DIR / "shards"
    )
    output_dir = (
        WORK_DIR / "smoke-merged" if args.smoke else RESULTS_DIR / "merged"
    )
    shards = [0] if args.smoke else range(config["shards"])
    output_dir.mkdir(parents=True, exist_ok=True)
    yodamerge = find_executable("yodamerge")
    for sample in selected_samples(config, args):
        sample_id = sample["id"]
        for profile in select_profiles(config, args.profile):
            inputs = [
                locate_yoda(
                    input_root
                    / sample_id
                    / profile
                    / f"shard-{shard:02d}"
                )
                for shard in shards
            ]
            output = output_dir / f"{sample_id}-{profile}.yoda"
            run_checked([yodamerge, "-o", output, *inputs], cwd=ROOT)


def command_plot(args: argparse.Namespace) -> None:
    config = load_config()
    rivet_mkhtml = find_executable("rivet-mkhtml")
    env = plotting_environment()
    samples = selected_samples(config, args)
    merged_root = (
        WORK_DIR / "smoke-merged"
        if args.smoke
        else RESULTS_DIR / "merged"
    )
    plot_root = (
        WORK_DIR / "smoke-plots"
        if args.smoke
        else RESULTS_DIR / "plots"
    )
    all_ids = sample_ids(config)
    selected_ids = [sample["id"] for sample in samples]
    selection_directory = (
        None
        if selected_ids == all_ids
        else "-vs-".join(selected_ids)
    )
    for profile in select_profiles(config, args.profile):
        inputs = []
        for sample in samples:
            path = merged_root / f"{sample['id']}-{profile}.yoda"
            if not path.is_file():
                raise RuntimeError(f"Missing merged YODA file: {path}")
            inputs.append(f"{path}:Title={sample['title']}")
        for output_format in ("PDF", "PNG"):
            output_root = plot_root / output_format.lower()
            if selection_directory is not None:
                output_root /= selection_directory
            output = output_root / profile
            run_checked(
                [
                    rivet_mkhtml,
                    "-o",
                    output,
                    "-f",
                    output_format,
                    "--pwd",
                    *inputs,
                ],
                cwd=ROOT,
                env=env,
            )
            extension = output_format.lower()
            if not any(output.rglob(f"*.{extension}")):
                raise RuntimeError(
                    f"rivet-mkhtml produced no {output_format} plots in {output}."
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


def yoda_number(obj: object, attribute: str) -> float:
    value = getattr(obj, attribute)
    return float(value() if callable(value) else value)


def histogram_overflow_totals(
    paths: Sequence[Path], histogram_name: str
) -> Dict[str, float]:
    try:
        import yoda
    except ImportError as error:
        raise RuntimeError(
            "The YODA Python module is required for output validation."
        ) from error

    totals = {
        "overflow_weight": 0.0,
        "total_weight": 0.0,
        "overflow_entries": 0.0,
        "total_entries": 0.0,
    }
    raw_path = f"/RAW/QCD_INSTANTON_KKS/{histogram_name}"
    for path in paths:
        objects = yoda.read(str(path))
        if raw_path not in objects:
            raise RuntimeError(
                f"Raw {histogram_name} is missing from {path}."
            )
        histogram = objects[raw_path]
        bins = list(histogram.bins())
        underflow = histogram.bin(0)
        overflow = histogram.bin(len(bins) + 1)
        regular_weight = sum(yoda_number(bin_, "sumW") for bin_ in bins)
        regular_entries = sum(
            yoda_number(bin_, "numEntries") for bin_ in bins
        )
        totals["overflow_weight"] += yoda_number(overflow, "sumW")
        totals["overflow_entries"] += yoda_number(overflow, "numEntries")
        totals["total_weight"] += (
            yoda_number(underflow, "sumW")
            + regular_weight
            + yoda_number(overflow, "sumW")
        )
        totals["total_entries"] += (
            yoda_number(underflow, "numEntries")
            + regular_entries
            + yoda_number(overflow, "numEntries")
        )
    if any(not math.isfinite(value) for value in totals.values()):
        raise RuntimeError(
            f"Non-finite {histogram_name} overflow accounting."
        )
    totals["weighted_fraction"] = (
        totals["overflow_weight"] / totals["total_weight"]
        if totals["total_weight"] != 0.0
        else 0.0
    )
    totals["unweighted_fraction"] = (
        totals["overflow_entries"] / totals["total_entries"]
        if totals["total_entries"] != 0.0
        else 0.0
    )
    return totals


def jets_mreco_overflow_totals(
    paths: Sequence[Path],
) -> Dict[str, Dict[str, float]]:
    return {
        histogram: histogram_overflow_totals(paths, histogram)
        for histogram in JET_MASS_HISTOGRAMS
    }


def herwig_cross_section_comparison(
    config: Dict[str, object],
    cross_sections: Dict[str, Dict[str, Dict[str, float]]],
) -> Dict[str, Dict[str, Dict[str, object]]]:
    baseline = cross_sections.get("herwig")
    if baseline is None:
        return {}
    comparison: Dict[str, Dict[str, Dict[str, object]]] = {}
    for sample in sample_definitions(config):
        sample_id = sample["id"]
        if (
            sample["engine"] != "herwig"
            or sample_id == "herwig"
            or sample_id not in cross_sections
        ):
            continue
        profile_report = {}
        for profile, reference in baseline.items():
            if profile not in cross_sections[sample_id]:
                continue
            variation = cross_sections[sample_id][profile]
            combined_error = math.hypot(
                reference["uncertainty"], variation["uncertainty"]
            )
            difference = variation["value"] - reference["value"]
            pull = (
                difference / combined_error
                if combined_error > 0.0
                else (0.0 if difference == 0.0 else math.inf)
            )
            profile_report[profile] = {
                "difference_pb": difference,
                "combined_uncertainty_pb": combined_error,
                "pull": pull,
                "compatible_3sigma": math.isfinite(pull)
                and abs(pull) <= 3.0,
            }
        comparison[sample_id] = profile_report
    return comparison


def command_summarize(args: argparse.Namespace) -> None:
    config = load_config()
    base = WORK_DIR / "smoke" if args.smoke else RESULTS_DIR / "shards"
    samples = selected_samples(config, args)
    report = {
        "mode": "smoke" if args.smoke else "production",
        "units": "pb",
        "cross_sections": {},
        "jets_mreco_overflow": {},
    }
    for sample in samples:
        sample_id = sample["id"]
        sample_report = {}
        overflow_report = {}
        for profile in select_profiles(config, args.profile):
            shards = sorted((base / sample_id / profile).glob("shard-*"))
            if not shards:
                raise RuntimeError(
                    f"No {sample_id} {profile} shards found under {base}."
                )
            if sample["engine"] == "herwig":
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
            sample_report[profile] = {
                "value": value,
                "uncertainty": uncertainty,
                "independent_estimates": len(estimates),
            }
            overflow_report[profile] = jets_mreco_overflow_totals(
                [locate_yoda(shard) for shard in shards]
            )
        report["cross_sections"][sample_id] = sample_report
        report["jets_mreco_overflow"][sample_id] = overflow_report

    comparison = herwig_cross_section_comparison(
        config, report["cross_sections"]
    )
    if comparison:
        report["herwig_cross_section_compatibility"] = comparison
    incompatible = (
        [
            f"{sample_id}/{profile}"
            for sample_id, profiles in comparison.items()
            for profile, result in profiles.items()
            if not result["compatible_3sigma"]
        ]
        if not args.smoke
        else []
    )

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
    if incompatible:
        raise RuntimeError(
            "Herwig hard-process cross sections are incompatible at "
            f"3 sigma: {', '.join(incompatible)}"
        )


def validate_yoda(
    path: Path, *, expected_events: int | None = None
) -> Dict[str, Dict[str, float]]:
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
    expected_histograms = {
        **{name: MASS_EDGES for name in JET_MASS_HISTOGRAMS},
        TRUTH_MASS_HISTOGRAM: MASS_EDGES,
        **{
            name: MIGRATION_RATIO_EDGES
            for name in MIGRATION_RATIO_HISTOGRAMS
        },
    }
    for histogram_name, expected_edges in expected_histograms.items():
        histogram_path = f"{analysis_path}{histogram_name}"
        if histogram_path not in objects:
            raise RuntimeError(
                f"Rivet {histogram_name} is missing from {path}."
            )
        edge_values = objects[histogram_path].xEdges
        if callable(edge_values):
            edge_values = edge_values()
        edges = [float(edge) for edge in edge_values]
        if len(edges) != len(expected_edges) or any(
            not math.isclose(
                found, expected, rel_tol=0.0, abs_tol=1.0e-9
            )
            for found, expected in zip(edges, expected_edges)
        ):
            raise RuntimeError(
                f"Unexpected {histogram_name} binning in {path}: "
                f"{edges[0] if edges else 'empty'}.."
                f"{edges[-1] if edges else 'empty'} with {len(edges)} edges."
            )
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
    raw_all_path = f"/RAW{analysis_path}_all"
    raw_truth_path = f"/RAW{analysis_path}_truth_mass_valid"
    if raw_all_path not in objects:
        raise RuntimeError(f"Raw all-events counter is missing from {path}.")
    if raw_truth_path not in objects:
        raise RuntimeError(
            f"Raw truth-mass counter is missing from {path}."
        )
    entries = yoda_number(objects[raw_all_path], "numEntries")
    truth_entries = yoda_number(objects[raw_truth_path], "numEntries")
    if not math.isclose(
        truth_entries, entries, rel_tol=0.0, abs_tol=0.5
    ):
        raise RuntimeError(
            f"Expected a valid truth mass for all events in {path}; "
            f"found {truth_entries} of {entries}."
        )
    if expected_events is not None:
        if not math.isclose(
            entries, float(expected_events), rel_tol=0.0, abs_tol=0.5
        ):
            raise RuntimeError(
                f"Expected {expected_events} events in {path}, found {entries}."
            )
    return jets_mreco_overflow_totals([path])


def command_validate(args: argparse.Namespace) -> None:
    config = load_config()
    static_checks(config)
    samples = selected_samples(config, args)
    profiles = select_profiles(config, args.profile)
    report = {
        "campaign": config["name"],
        "total_events": (
            sum(region["events"] for region in config["regions"].values())
            * len(config["processes"])
            * len(sample_definitions(config))
        ),
        "sherpa_all_processes": len(sherpa_processes("all")),
        "cards": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in rendered_cards(config)
        },
    }
    if args.smoke:
        checked = []
        overflow = {}
        for sample in samples:
            sample_id = sample["id"]
            overflow[sample_id] = {}
            for profile in profiles:
                directory = (
                    WORK_DIR
                    / "smoke"
                    / sample_id
                    / profile
                    / "shard-00"
                )
                path = locate_yoda(directory)
                provenance = json.loads(
                    (directory / "provenance.json").read_text(
                        encoding="utf-8"
                    )
                )
                overflow[sample_id][profile] = validate_yoda(
                    path, expected_events=int(provenance["events"])
                )
                checked.append(str(path.relative_to(ROOT)))
        report["validated_smoke_yoda"] = checked
        report["smoke_jets_mreco_overflow"] = overflow
    merged = RESULTS_DIR / "merged"
    if args.results:
        checked = []
        overflow = {}
        for sample in samples:
            sample_id = sample["id"]
            overflow[sample_id] = {}
            for profile in profiles:
                path = merged / f"{sample_id}-{profile}.yoda"
                overflow[sample_id][profile] = validate_yoda(path)
                checked.append(str(path.relative_to(ROOT)))
        report["validated_yoda"] = checked
        report["jets_mreco_overflow"] = overflow
    output = RESULTS_DIR / "validation.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2, sort_keys=True))


def add_sample_arguments(parser: argparse.ArgumentParser) -> None:
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument(
        "--sample",
        help="all or a comma-separated list of campaign sample IDs",
    )
    selection.add_argument(
        "--generator",
        choices=("all", "herwig", "sherpa"),
        help="backward-compatible alias for an original sample ID",
    )


def add_progress_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--progress-interval",
        type=finite_seconds,
        default=5.0,
        metavar="SECONDS",
        help=(
            "tracker refresh interval (default: 5; negative disables "
            "intermediate refreshes)"
        ),
    )
    parser.add_argument(
        "--max-listed",
        type=nonnegative_integer,
        default=12,
        metavar="N",
        help="maximum active shards and incomplete groups shown (default: 12)",
    )


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
    build.add_argument(
        "--jobs",
        type=positive_integer,
        default=max(1, os.cpu_count() or 1),
    )
    build.set_defaults(function=command_build)

    smoke = subparsers.add_parser("smoke", help="run small validation samples")
    add_sample_arguments(smoke)
    smoke.add_argument("--profile", default="all")
    smoke.add_argument("--events", type=positive_integer, default=20)
    smoke.add_argument("--jobs", type=positive_integer, default=1)
    smoke.add_argument("--force", action="store_true")
    add_progress_arguments(smoke)
    smoke.set_defaults(function=command_smoke)

    run = subparsers.add_parser("run", help="run production shards")
    add_sample_arguments(run)
    run.add_argument("--profile", default="all")
    run.add_argument("--shards", default="all",
                     help="all or a comma-separated zero-based list")
    run.add_argument("--jobs", type=positive_integer, default=1)
    run.add_argument("--force", action="store_true")
    add_progress_arguments(run)
    run.set_defaults(function=command_run)

    status = subparsers.add_parser(
        "status", help="show the latest campaign tracker snapshot"
    )
    status.add_argument(
        "--smoke",
        action="store_true",
        help="show the smoke tracker instead of production",
    )
    status.add_argument(
        "--json",
        action="store_true",
        help="print the machine-readable tracker snapshot",
    )
    status.set_defaults(function=command_status)

    merge = subparsers.add_parser("merge", help="merge YODA shard files")
    add_sample_arguments(merge)
    merge.add_argument("--profile", default="all")
    merge.add_argument(
        "--smoke",
        action="store_true",
        help="merge the one-shard smoke outputs instead of production",
    )
    merge.set_defaults(function=command_merge)

    plot = subparsers.add_parser("plot", help="plot sample overlays")
    add_sample_arguments(plot)
    plot.add_argument("--profile", default="all")
    plot.add_argument(
        "--smoke",
        action="store_true",
        help="plot merged smoke outputs instead of production",
    )
    plot.set_defaults(function=command_plot)

    summarize = subparsers.add_parser(
        "summarize", help="extract generator cross sections"
    )
    add_sample_arguments(summarize)
    summarize.add_argument("--profile", default="all")
    summarize.add_argument("--smoke", action="store_true")
    summarize.set_defaults(function=command_summarize)

    validate = subparsers.add_parser("validate", help="validate configuration")
    add_sample_arguments(validate)
    validate.add_argument("--profile", default="all")
    validate.add_argument("--smoke", action="store_true")
    validate.add_argument("--results", action="store_true")
    validate.set_defaults(function=command_validate)
    return parser


def main() -> int:
    parser = make_parser()
    args = parser.parse_args()
    try:
        args.function(args)
    except KeyboardInterrupt:
        _CANCEL_REQUESTED.set()
        terminate_active_processes()
        print("error: interrupted", file=sys.stderr)
        return 130
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
