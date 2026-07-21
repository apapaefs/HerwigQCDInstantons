# Herwig QCD Instantons

This repository provides a Herwig 7.3 model for phenomenological QCD-instanton
events,

```text
g g -> one q qbar pair per active flavour + additional gluons
```

together with a repaired MAMBO many-body phase-space generator, reference
Herwig cards, and a particle-level Rivet analysis. The supplied cards use the
Khoze-Krauss-Schott (KKS) model of
[arXiv:1911.09726](https://arxiv.org/abs/1911.09726).

This is not a first-principles event-by-event instanton amplitude. It generates
the allowed high-multiplicity final states and assigns either simple
multiplicity weights or interpolated KKS partonic rates.

## Repository map

- `MEInstanton.{h,cc}`: process registration, model weights, flavour selection,
  scales, and colour connections.
- `Phasespace/MamboPhasespace.{h,cc}`: stochastic, non-invertible many-body
  phase space with internal accept/reject unweighting.
- `LHC-Instanton*.in`: ready-to-read Herwig cards.
- `Rivet/`: source, metadata, plotting configuration, and build file for
  `QCD_INSTANTON_KKS`.

## Requirements

- Herwig 7.3 with ThePEG and LHAPDF support.
- GSL, including `gsl-config`; `VariableKKS` uses GSL's hypergeometric
  implementation while constructing its interpolation table.
- The `NNPDF31_nnlo_as_0118` LHAPDF set for the supplied KKS cards.
- Rivet for the analysis cards and plugin.

## Installation

MAMBO is a Herwig-core phase-space class, while `MEInstanton` is a loadable
contrib plugin. Let `HERWIG_SRC` denote the configured Herwig source tree and
activate the corresponding installation first.

```sh
source /path/to/Herwig/bin/activate
cp Phasespace/MamboPhasespace.{cc,h} \
  "$HERWIG_SRC/MatrixElement/Matchbox/Phasespace/"
```

Add `MamboPhasespace.h` to `ALL_H_FILES` and `MamboPhasespace.cc` to
`ALL_CC_FILES` in
`$HERWIG_SRC/MatrixElement/Matchbox/Phasespace/Makefile.am`. Regenerate the
configured build files if necessary, then rebuild and install Herwig:

```sh
make -C "$HERWIG_SRC/MatrixElement/Matchbox/Phasespace"
make -C "$HERWIG_SRC" install
```

Install the matrix-element plugin through Herwig's contrib machinery:

```sh
cp -R /path/to/HerwigQCDInstantons \
  "$HERWIG_SRC/Contrib/HerwigQCDInstantons"
cd "$HERWIG_SRC/Contrib"
bash make_makefiles.sh
make -C HerwigQCDInstantons
make -C HerwigQCDInstantons install
```

Build the Rivet plugin when it is needed:

```sh
make -C Rivet
```

## Running the examples

From the repository root:

```sh
Herwig read LHC-Instanton.in
Herwig run LHC-Instanton.run -N 1000
```

The Rivet variants are:

- `LHC-Instanton-Rivet.in`: general example with direct YODA output.
- `LHC-Instanton-Rivet-Table3-Low.in`: 30 GeV lower mass cut.
- `LHC-Instanton-Rivet-Table3-High.in`: 500 GeV lower mass cut.

All KKS cards cap the tabulated range at `2895.5*GeV`.

## Matrix-element options

### `MEModeling`

Default: `PureMultiplicity`.

- `PureMultiplicity`: flat base matrix element multiplied by the selected toy
  gluon-multiplicity factor.
- `KKS`: interpolated KKS inclusive partonic cross section, mean gluon
  multiplicity, instanton scale, PDF-scale reweighting, and optional KKS
  flavour partition.

```text
set MEInstanton:MEModeling KKS
```

### `NAdditional`

Default: `0`. Valid range: any non-negative integer.

This inherited `BlobME` parameter is the maximum number of additional gluons.
The code registers every multiplicity from `0` through `NAdditional`, so `10`
means eleven gluon channels. There is no hard-coded ten-gluon limit; larger
values are supported, but process count and phase-space cost grow quickly.

In KKS mode, the Poisson probabilities are normalized over the retained range
`0..NAdditional`. Changing the cap therefore redistributes the truncated tail
without changing the tabulated inclusive rate. `PureMultiplicity` retains its
original untruncated Poisson or Gaussian factors.

### `QuarkPairs`

Default: `Fixed`.

- `Fixed`: use `NQuarkPair` pairs.
- `Variable`: legacy equal KKS weighting of the four- and five-pair channels.
- `VariableKKS`: scale- and mass-dependent four/five-pair probabilities. This
  option is valid only with `MEModeling KKS`; other combinations are rejected
  during initialization.

`NQuarkPair` defaults to `4` and is restricted to `1..5`. It applies only to
`Fixed`. Flavours follow PDG order: `d`, `u`, `s`, `c`, `b`.

`KKSBottomMass` defaults to `4.18*GeV` and must be non-negative. It controls the
KKS active-bottom condition; it is independent of the kinematic mass assigned
to the hard-process bottom legs.

```text
set MEInstanton:QuarkPairs VariableKKS
set MEInstanton:KKSBottomMass 4.18*GeV
```

### `VariableKKS` definition

At each of the 20 KKS table nodes the implementation forms

```text
u         = sqrt(shat) * rho
rho_tilde = alpha_s(1/rho) * u / (4*pi)
S'(chi)   = rho_tilde
```

The last equation is solved by bounded bisection using the full valley action
`S(chi)` from KKS. With

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
W4 = kappa4^2 * omega^8,    W5 = 0.
```

Otherwise,

```text
W4 = kappa4^2 * (m_b*rho)^2 * omega^8,
W5 = kappa5^2 * omega^10,
kappa4 = 0.008,             kappa5 = 0.01.
```

`W4` and `W5` are normalized to probabilities. Their sum therefore partitions,
rather than rescales, the inclusive KKS cross section. This is an explicit
modelling choice: it uses the common tabulated saddle point and preserves the
published inclusive table instead of recomputing separate four- and
five-flavour cross sections.

The tabulated `alpha_s(1/rho)` values now determine `chi` and `omega`; they are
not multiplied into the already-inclusive `sigmahat` a second time.

### KKS interpolation and scales

The linear tables cover

```text
10.7 GeV <= sqrt(shat) <= 2895.5 GeV.
```

Both bounds are enforced. The first cross-section node is
`4.922e9 pb`. The interpolated table stores the inclusive four-plus-five
flavour result.

`FactorizationScale` defaults to `InvRho`:

- `InvRho`: use `(1/rho)^2` from the table.
- `sHat`: use the partonic invariant mass squared.

KKS mode reweights the incoming gluon PDFs from `sHat` to this scale. Zero,
negative, or non-finite PDFs, scales, Jacobians, and interpolated factors reject
the point before any division.

### Pure-multiplicity parameters

`MultiplicityParametrisation` defaults to `Poisson`:

- `Poisson`: `mean^n * exp(-mean) / n!`.
- `Gaussian`: `exp(-(n-A)^2/B) / sqrt(pi*B)`.
- `Flat`: no multiplicity factor.
- `UserDefined`: reserved; currently equivalent to `Flat`.

Parameters and constraints:

| Parameter | Default | Valid values |
| --- | ---: | --- |
| `PoissonMean` | `3.0` | `>= 0` |
| `GaussianParamA` | `5.0` | `>= 0` |
| `GaussianParamB` | `200.0` | `> 0` |

These settings are ignored by `KKS`, which uses its interpolated mean and a
numerically stable Poisson distribution normalized over `0..NAdditional`.

### `ColourConnections`

Default: `Simple`.

- `Simple`: deterministic incoming singlet; quark and gluon lines are paired,
  with one gluon joining the last quark pair for odd gluon multiplicity.
- `Random`: randomized final-state map with singlet incoming gluons.
- `Random2`: randomized map containing one incoming colour line.
- `Random3`: randomized map containing both incoming colour lines; used by the
  supplied cards.

All modes support `NQuarkPair=1..5` and arbitrary retained gluon multiplicity.
Random modes generate a one-to-one map, use integer random indices, and reject
gluon self-connections. These remain phenomenological shower colour models,
not exact instanton colour amplitudes.

## Phase-space options

The cards create repaired MAMBO and select it by default:

```text
set MEInstanton:Phasespace /Herwig/MatrixElements/Matchbox/Phasespace/MamboPS
#set MEInstanton:Phasespace /Herwig/MatrixElements/Matchbox/Phasespace/InvertiblePhasespace
```

MAMBO is stochastic and non-invertible. It consumes Herwig's random stream
internally and declares one dummy integration coordinate. It preserves
requested masses, checks thresholds and convergence, and returns unweighted
accepted configurations.

MAMBO parameters:

| Parameter | Default | Valid values | Meaning |
| --- | ---: | --- | --- |
| `MaxWeight` | `10.0` | `1e-12..1e12` | Strict accept/reject upper bound; a violation throws. |
| `MaxTrials` | `100000` | `1..100000000` | Safety limit for stochastic rejection loops. |

Herwig's `InvertiblePhasespace` is deterministic and invertible. Uncomment its
line and comment the MAMBO line for sampler studies or an independent
phase-space comparison.

The reference cards explicitly set `u,d,s,c,b` and their antiparticles to zero
hard-process mass. This gives MAMBO and invertible phase space the same
massless convention and leaves `KKSBottomMass` as the physical mass used only
for flavour selection.

## Shower reconstruction

Every supplied card contains:

```text
set /Herwig/Shower/KinematicsReconstructor:ReconstructionOption General
```

In a controlled Herwig 7.3 high-multiplicity instanton sample, the default
`Colour3` reconstruction rejected about 88% of shower attempts. `General`
reduced that rate to about 3%. Keeping `Colour3` would therefore select a small,
strongly biased survivor population before analysis.

## Rivet analysis

`QCD_INSTANTON_KKS` defines:

- tracks: charged stable particles with `pT > 0.5 GeV`, `|eta| < 2.5`;
- jets: anti-kT `R=0.4` particle jets with `pT > 20 GeV`, `|eta| < 4.5`;
- reconstructed-mass proxies: invariant masses of the selected track and jet
  systems;
- low track windows: `25 < Mreco < 35 GeV`, plus the
  `20 < Mreco < 30 GeV` `ST` variant;
- high jet window: `320 < Mreco < 480 GeV`;
- observables: multiplicity, scalar `ST`, average pairwise `Delta phi`, and
  sphericity.

Sphericity is calculated after boosting the selected tracks or jets into their
combined reconstructed-system rest frame. Histograms with finite positive
integrals are normalized to unit area.

The counters provide the accepted-event denominator and selection numerators
needed to form efficiencies among events delivered to Rivet. Rivet cannot
observe hard-process, phase-space, or shower attempts rejected before an event
reaches the analysis, so it cannot diagnose shower rejection efficiency by
itself.

This is a particle-level interpretation of the KKS Section 4 observables, not a
Delphes or detector-level reproduction.

## Practical checks

Read each card before launching a long run:

```sh
for card in LHC-Instanton*.in; do Herwig read "$card"; done
```

For phase-space comparisons, read once with the active MAMBO line and once with
the adjacent `InvertiblePhasespace` line. Use the same seed, hard-process mass
convention, cuts, multiplicity cap, and model settings. Compare integrated
rates within Monte Carlo uncertainty and inspect principal parton-level mass,
multiplicity, and momentum distributions.

Generated `.run`, `.log`, `.out`, `.tex`, `.yoda`, shared-library, and plot
artifacts are intentionally excluded from version control.
