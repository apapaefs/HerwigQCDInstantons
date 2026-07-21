# Rivet analysis: QCD_INSTANTON_KKS

Particle-level Rivet implementation of the Section 4 observables in Khoze,
Krauss and Schott, arXiv:1911.09726.

## What it measures

- Tracks: charged final-state particles with `pT > 0.5 GeV` and `|eta| < 2.5`.
- Jets: anti-kT `R = 0.4` particle jets with `pT > 20 GeV` and `|eta| < 4.5`.
- Mass proxies: invariant masses of the selected track and jet systems.
- Low-mass track window: `25 < M_I^reco < 35 GeV`, plus the `20 < M_I^reco < 30 GeV` track-`ST` variant used in the paper's Fig. 9 caption.
- High-mass jet window: `320 < M_I^reco < 480 GeV`.
- Shapes: multiplicity, scalar `ST`, average pairwise `Delta phi`, and
  sphericity. Sphericity is computed after boosting the selected track or jet
  system into its combined reconstructed-system rest frame.

Histograms with a finite positive integral are normalized to unit area in
`finalize()`, matching the shape comparisons in the paper. The counters provide
the accepted-event denominator and selection numerators for efficiencies:
Rivet sees only events delivered by Herwig and cannot count hard-process or
shower attempts rejected before analysis.

## Build

```sh
source /path/to/Herwig/bin/activate
make
```

This creates `RivetQCD_INSTANTON_KKS.so`.

## Run

```sh
rivet --pwd -a QCD_INSTANTON_KKS events.hepmc -o qcd-instanton.yoda
rivet-mkhtml --pwd qcd-instanton.yoda -o plots
```

For the direct Herwig example from the repository root:

```sh
unset TEXMFCNF TEXMFHOME TEXINPUTS
export RIVET_PLOT_PATH="$PWD/Rivet${RIVET_PLOT_PATH:+:$RIVET_PLOT_PATH}"
export RIVET_INFO_PATH="$PWD/Rivet${RIVET_INFO_PATH:+:$RIVET_INFO_PATH}"
MPLCONFIGDIR=/private/tmp/mplconfig-rivet rivet-mkhtml LHC-Instanton-Rivet.yoda -o plots
```

This is not a Delphes reproduction of the original plots; it is a Rivet
particle-level analysis designed for fast comparisons of QCD-instanton and
background samples.

From another directory, replace `--pwd` with `--analysis-path /path/to/HerwigQCDInstantons/Rivet`.
