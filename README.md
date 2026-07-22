# Herwig QCD Instantons

This repository provides a Herwig 7.3 model for phenomenological
QCD-instanton events,

```text
g g -> one q qbar pair per active flavour + additional gluons
```

together with a repaired MAMBO many-body phase-space generator, reference
Herwig input cards, and a particle-level Rivet analysis. The supplied cards use
the Khoze-Krauss-Schott (KKS) calculation in
[arXiv:1911.09726](https://arxiv.org/abs/1911.09726).

The implementation is not a first-principles exclusive instanton amplitude. It
enumerates the requested high-multiplicity final states, generates their phase
space, and assigns either simple multiplicity weights or interpolated inclusive
KKS partonic rates. Flavour selection and shower colour flow are explicit
phenomenological models.

## Contents

- [`MEInstanton.cc`](MEInstanton.cc) and
  [`MEInstanton.h`](MEInstanton.h): process registration, matrix-element
  weights, flavour selection, scales, and colour connections.
- [`Phasespace/MamboPhasespace.cc`](Phasespace/MamboPhasespace.cc) and
  [`Phasespace/MamboPhasespace.h`](Phasespace/MamboPhasespace.h): stochastic,
  non-invertible many-body phase space with internal accept/reject unweighting.
- [`LHC-Instanton*.in`](LHC-Instanton.in): ready-to-read Herwig cards.
- [`Rivet/`](Rivet/README.md): analysis source, metadata, plot configuration,
  and build instructions for `QCD_INSTANTON_KKS`.

## Requirements

- Herwig 7.3 with ThePEG and LHAPDF support.
- GSL, including `gsl-config`. `VariableKKS` uses GSL's hypergeometric
  implementation while constructing its interpolation table.
- The `NNPDF31_nnlo_as_0118` LHAPDF set for the supplied cards.
- Rivet for the analysis plugin and Rivet-enabled cards.

## Installation

MAMBO is compiled into the Herwig core, while `MEInstanton` is a loadable
contrib plugin. Let `HERWIG_SRC` be the configured Herwig 7.3 source tree and
activate the corresponding installation.

### 1. Install MAMBO

```sh
export HERWIG_SRC=/path/to/Herwig-7.3-source
source /path/to/Herwig-install/bin/activate

cp Phasespace/MamboPhasespace.cc \
   Phasespace/MamboPhasespace.h \
  "$HERWIG_SRC/MatrixElement/Matchbox/Phasespace/"
```

Add `MamboPhasespace.h` to `ALL_H_FILES` and `MamboPhasespace.cc` to
`ALL_CC_FILES` in
`$HERWIG_SRC/MatrixElement/Matchbox/Phasespace/Makefile.am`. The
[`Phasespace/Makefile.am`](Phasespace/Makefile.am) file shows the required
entries; do not replace a differently versioned Herwig makefile wholesale.

Regenerate the configured build files when required, then rebuild and install
Herwig. For an in-source build this is typically:

```sh
make -C "$HERWIG_SRC/MatrixElement/Matchbox/Phasespace"
make -C "$HERWIG_SRC" install
```

Use the corresponding build-directory paths for an out-of-source build.

### 2. Install `MEInstanton`

Copy this repository into Herwig's `Contrib` directory and use Herwig's contrib
build machinery:

```sh
cp -R /path/to/HerwigQCDInstantons \
  "$HERWIG_SRC/Contrib/HerwigQCDInstantons"
cd "$HERWIG_SRC/Contrib"
bash make_makefiles.sh
make -C HerwigQCDInstantons
make -C HerwigQCDInstantons install
```

This installs `Instantons.so` in the active Herwig plugin directory.

### 3. Build the Rivet plugin

From the repository or contrib checkout:

```sh
make -C Rivet
```

This creates `Rivet/RivetQCD_INSTANTON_KKS.so`.

## Quick Start

Run from the repository root so that the Rivet cards can resolve their local
`Rivet/` plugin path.

Without Rivet:

```sh
Herwig read LHC-Instanton.in
Herwig run LHC-Instanton.run -N 1000
```

With Rivet:

```sh
make -C Rivet
Herwig read LHC-Instanton-Rivet-Table3-Low.in
Herwig run LHC-Instanton-Rivet-Table3-Low.run -N 1000
```

The `-N` value overrides the event count stored in the card.

### Reference cards

| Card | Rivet | `MHatMin` | Intended use |
| --- | --- | ---: | --- |
| `LHC-Instanton.in` | No | `20 GeV` | Generator-level reference |
| `LHC-Instanton-Rivet.in` | Yes | `20 GeV` | General particle-level sample |
| `LHC-Instanton-Rivet-Table3-Low.in` | Yes | `30 GeV` | KKS low-mass region |
| `LHC-Instanton-Rivet-Table3-High.in` | Yes | `500 GeV` | KKS high-mass region |

Every card sets `MHatMax = 2895.5 GeV`, the upper edge of the KKS table. Their
shared reference configuration is:

| Component | Supplied setting |
| --- | --- |
| Matrix-element model | `KKS` |
| Flavour model | `VariableKKS`, `KKSBottomMass = 4.18 GeV` |
| Gluon cap | `NAdditional = 25` |
| Common hard scale | `sHatOverN` |
| Phase space | repaired MAMBO |
| Hard-process quark masses | nominal Herwig `ParticleData` masses |
| Colour flow | `Random3` |
| Shower reconstruction | `General` |

These are card choices, not the C++ interface defaults documented below.

## Process and Multiplicity

For `n_q` selected quark pairs and `n_g` additional gluons, the registered
subprocess is

```text
g g -> d dbar [u ubar] [s sbar] [c cbar] [b bbar] + n_g gluons
```

with flavours included in PDG order. `NAdditional` is the maximum `n_g`.
Every multiplicity from `0` through `NAdditional` is registered, so

```text
set MEInstanton:NAdditional 10
```

creates eleven gluon-multiplicity channels. There is no hard-coded ten-gluon
limit. Larger values increase the number of subprocesses and the phase-space
cost.

The interface default is `NAdditional = 0`. The supplied cards use `25`.
At the largest tabulated KKS mean, `<N_g> = 12.14`, the range `0..25`
contains about `99.9635%` of the ordinary Poisson probability.

In KKS mode the Poisson distribution is normalized over the retained
`0..NAdditional` channels. Changing the cap therefore repartitions the
inclusive KKS rate rather than changing it. `PureMultiplicity` retains its
ordinary, untruncated multiplicity factors.

## Matrix-Element Models

### `MEModeling`

Interface default: `PureMultiplicity`.

| Value | Meaning |
| --- | --- |
| `PureMultiplicity` | Flat base matrix element times the selected toy multiplicity factor |
| `KKS` | Interpolated inclusive KKS partonic rate, mean gluon multiplicity, hard scale, and optional KKS flavour partition |

The supplied cards select:

```text
set MEInstanton:MEModeling KKS
```

### Pure-multiplicity parameters

`MultiplicityParametrisation` defaults to `Poisson`.

| Value | Multiplicity factor |
| --- | --- |
| `Poisson` | `mean^n * exp(-mean) / n!` |
| `Gaussian` | `exp(-(n-A)^2/B) / sqrt(pi*B)` |
| `Flat` | No multiplicity-dependent factor |
| `UserDefined` | Reserved; currently equivalent to `Flat` |

| Parameter | Default | Constraint |
| --- | ---: | --- |
| `PoissonMean` | `3.0` | `>= 0` |
| `GaussianParamA` | `5.0` | `>= 0` |
| `GaussianParamB` | `200.0` | `> 0` |

These settings do not affect `KKS`.

## Flavour Models

`QuarkPairs` controls the number of quark pairs.

| Value | Meaning |
| --- | --- |
| `Fixed` | Use `NQuarkPair`; this is the interface default |
| `Variable` | Legacy 50/50 selection between four and five pairs |
| `VariableKKS` | Scale- and mass-dependent four/five-pair probabilities |

`NQuarkPair` defaults to `4`, is restricted to `1..5`, and applies only to
`Fixed`. `VariableKKS` is accepted only with `MEModeling KKS`.

`KKSBottomMass` defaults to `4.18 GeV` and must be non-negative. It enters
only the KKS flavour-selection condition; it does not set the kinematic bottom
mass.

### `VariableKKS`

At each of the 20 KKS table nodes the implementation computes

```text
u         = sqrt(shat) * rho
rho_tilde = alpha_s(1/rho) * u / (4*pi)
S'(chi)   = rho_tilde
```

using the tabulated `alpha_s(1/rho)`, a bounded root solve, and the full KKS
valley action. With

```text
z = (2 + chi^2 + chi*sqrt(4 + chi^2))/2
```

the fermion overlap is

```text
omega = (3*pi/8) * z^(-3/2)
        * 2F1(3/2, 3/2; 4; 1 - 1/z^2).
```

The resulting `omega` values are interpolated at runtime. For
`m_b*rho > 1`,

```text
W4 = kappa4^2 * omega^8
W5 = 0
```

and otherwise

```text
W4 = kappa4^2 * (m_b*rho)^2 * omega^8
W5 = kappa5^2 * omega^10

kappa4 = 0.008
kappa5 = 0.01
```

`W4` and `W5` are normalized to probabilities, so they partition rather than
rescale the inclusive KKS cross section. This preserves the published table
instead of claiming separately calculated exclusive four- and five-flavour
rates.

This differs from stock Sherpa, which uses configurable deterministic
`sqrt(shat)` thresholds for charm and bottom, `20` and `100 GeV` by default.
Herwig `VariableKKS` instead closes the bottom channel when `m_b*rho > 1`
and otherwise samples both channels according to `W4/(W4+W5)` and
`W5/(W4+W5)`. The legacy `Variable` option remains the independent 50/50
model.

## KKS Table and Scales

The linearly interpolated KKS tables cover

```text
10.7 GeV <= sqrt(shat) <= 2895.5 GeV
```

and both bounds are enforced. The corrected first cross-section node is
`4.922e9 pb`; the table contains the inclusive four-plus-five-flavour result.

The published partonic cross sections already contain their coupling
dependence and initial-state averaging, including the KKS `1/36` factor. Do
not apply another `1/36` or an additional perturbative
`alpha_s(mu_R)` factor. The tabulated `alpha_s(1/rho)` values are used to
derive `chi` and `omega` for `VariableKKS`; they are not multiplied into
`sigmahat` again.

`FactorizationScale` is the historical interface name for the common KKS hard
scale. Its interface default is `InvRho`.

| Value | Scale |
| --- | --- |
| `InvRho` | `mu^2 = (1/rho)^2` |
| `sHat` | `mu^2 = shat` |
| `sHatOverN` | `mu^2 = shat/<N_g>`, using the interpolated KKS mean |

`MEInstanton::scale()` returns this value, so ThePEG evaluates the incoming
PDFs at the selected scale and records it on the hard-process particles used to
initialize the shower. There is no second manual PDF ratio.

The supplied cards select `sHatOverN`, matching Sherpa's `Democratic`
`shat/N` setup:

```text
mu_F^2 = mu_R^2 = mu_Q^2 = shat/<N_g>
```

Non-finite or non-positive scales and phase-space Jacobians, together with
invalid interpolated factors, are rejected before division.

## Colour Models

`ColourConnections` defaults to `Simple`; the supplied cards select
`Random3`. All choices are phenomenological leading-colour assignments for
the shower and hadronization, not exact instanton colour amplitudes.

| Value | Construction |
| --- | --- |
| `Simple` | Fixed incoming singlet; direct quark pairs and two-gluon singlets, with one gluon inserted in the last quark chain when `n_g` is odd |
| `Random` | Random outgoing endpoint map with the incoming singlet fixed |
| `Random2` | Random outgoing map containing one crossed incoming source/sink pair; the other incoming line remains fixed |
| `Random3` | Both crossed incoming gluons participate in a fully mixed endpoint map |
| `QCDINSPlanar` | Random antiquark permutation with every incoming and outgoing gluon placed on an open `q-g-...-g-qbar` string |

All modes support `NQuarkPair=1..5` and arbitrary retained gluon multiplicity.
The randomized constructions use one-to-one endpoint maps and forbid a gluon
from connecting directly to itself.

### `Random3` and Sherpa

`Random3` uses the same broad endpoint pool as Sherpa's instanton model: both
ends of both incoming gluons and all outgoing quark, antiquark, and gluon
endpoints participate. The sampling algorithms differ. Herwig shuffles a
complete anticolour permutation and rejects a map containing any gluon
self-connection. Sherpa selects pairs sequentially and repairs the last pair
through an outgoing leg when necessary. They therefore allow the same broad
topology class but do not assign identical probabilities to individual flows.

### `QCDINSPlanar`

`QCDINSPlanar` adapts the prescription in
[QCDINS 2.0](https://arxiv.org/abs/hep-ph/9911516). It randomly pairs quarks
with antiquarks, distributes all gluons among the resulting strings, and joins
adjacent partons. Every gluon belongs to a component containing a quark and an
antiquark, although an interior gluon may have only gluons as immediate
neighbours. Independent pure-gluon loops are excluded.

The motivation is a leading-`N_c` planar representation of the
colour-averaged instanton final state. It is not an exclusive colour
probability derived from the KKS amplitude. QCDINS was also formulated for a
DIS `q' g` subprocess; placing both crossed incoming gluons on the strings is
the explicit `g g` adaptation used here. Compare `QCDINSPlanar` with
`Random3` as a colour-model systematic.

Each input card explains all five choices and leaves only this line active:

```text
set MEInstanton:ColourConnections Random3
```

## Phase Space and Masses

The supplied cards create repaired MAMBO and select it with:

```text
set MEInstanton:Phasespace /Herwig/MatrixElements/Matchbox/Phasespace/MamboPS
```

MAMBO is stochastic and non-invertible. It consumes Herwig's random stream
internally and presents one dummy integration coordinate to the sampler. It
preserves the requested final-state masses, checks thresholds and convergence,
and bounds all rejection loops.

| MAMBO parameter | Default | Range | Meaning |
| --- | ---: | --- | --- |
| `MaxWeight` | `10.0` | `1e-12..1e12` | Strict accept/reject weight bound; a violation throws |
| `MaxTrials` | `100000` | `1..100000000` | Maximum attempts in each rejection loop |

For deterministic, invertible phase space, comment the MAMBO line and uncomment
the adjacent alternative:

```text
set MEInstanton:Phasespace /Herwig/MatrixElements/Matchbox/Phasespace/InvertiblePhasespace
```

The cards call `UnsetHardProcessMass` for `u,d,s,c,b` and their antiparticles.
Both phase-space choices therefore receive the same nominal Herwig
`ParticleData` masses. `KKSBottomMass` remains a separate parameter used
only for flavour selection.

## Shower Reconstruction

Every supplied card uses:

```text
set /Herwig/Shower/KinematicsReconstructor:ReconstructionOption General
```

In a controlled Herwig 7.3 high-multiplicity instanton sample, the default
`Colour3` reconstruction rejected about `88%` of shower attempts.
`General` reduced the rejection rate to about `3%`. The default would
therefore select a small, strongly biased survivor population before Rivet or
any other event analysis sees the sample.

## Rivet Analysis

`QCD_INSTANTON_KKS` is a particle-level interpretation of the KKS Section 4
observables, not a Delphes or detector-level reproduction.

Object definitions:

- tracks: charged stable particles with `pT > 0.5 GeV` and `|eta| < 2.5`;
- jets: anti-`kT`, `R=0.4` particle jets with `pT > 20 GeV` and
  `|eta| < 4.5`;
- reconstructed-mass proxies: invariant masses of the selected track and jet
  systems;
- low track window: `25 < Mreco < 35 GeV`, with a separate
  `20 < Mreco < 30 GeV` track-`ST` histogram;
- high jet window: `320 < Mreco < 480 GeV`.

The analysis fills multiplicity, scalar `ST`, average pairwise
`Delta phi`, and sphericity. Sphericity is calculated after boosting the
selected tracks or jets into their combined reconstructed-system rest frame.
Only histograms with finite positive integrals are normalized to unit area.

Counters record accepted-event denominators and selection numerators. Rivet
cannot observe hard-process, phase-space, or shower attempts rejected before an
event reaches the analysis, so these counters are accepted-event efficiencies,
not generator-attempt efficiencies.

### Plotting

After producing a YODA file, generate HTML plus the default PDF and PNG plots
from the repository root:

```sh
unset TEXMFCNF TEXMFHOME TEXINPUTS
MPLCONFIGDIR=/private/tmp/mplconfig-rivet \
  rivet-mkhtml --no-rivet-refs \
  -c Rivet/QCD_INSTANTON_KKS.plot \
  LHC-Instanton-Rivet-Table3-Low.yoda \
  -o plots-low
```

Replace the YODA filename and output directory for another sample. Omitting
`--format`, as above, avoids installations that pass a literal `pdf,png`
string to Matplotlib instead of splitting it into two formats.

See [`Rivet/README.md`](Rivet/README.md) for standalone HepMC usage.

## Validation and Reproducibility

Read every card before starting a production run:

```sh
for card in LHC-Instanton*.in; do
  Herwig read "$card"
done
```

For a MAMBO versus invertible comparison, change only the adjacent
`Phasespace` line and keep the random seed, masses, cuts, multiplicity cap,
model settings, and integration statistics fixed. Compare integrated rates and
principal parton-level distributions within Monte Carlo uncertainty.

Generated `.run`, `.log`, `.out`, `.tex`, `.yoda`, shared-library, and
plot products are excluded from version control.
