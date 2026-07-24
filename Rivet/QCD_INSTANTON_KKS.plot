BEGIN PLOT /QCD_INSTANTON_KKS/tracks_mreco
Title=Track-based reconstructed instanton mass proxy
XLabel=$M^{\mathrm{reco}}_{I,\mathrm{tracks}}$ [GeV]
YLabel=Normalized events
END PLOT

BEGIN PLOT /QCD_INSTANTON_KKS/jets_mreco_inclusive_eta45
Title=Inclusive jet-system mass, $|\eta_j| < 4.5$
XLabel=$M^{\mathrm{reco}}_{I,\mathrm{jets}}$ [GeV]
YLabel=Normalized events
END PLOT

BEGIN PLOT /QCD_INSTANTON_KKS/jets_mreco_central
Title=Central jet-system mass, $|\eta_j| < 2.5$
XLabel=$M^{\mathrm{reco,central}}_{I,\mathrm{jets}}$ [GeV]
YLabel=Normalized events
END PLOT

BEGIN PLOT /QCD_INSTANTON_KKS/instanton_mass_truth
Title=Pre-shower hard-process instanton mass
XLabel=$\sqrt{\hat{s}}$ [GeV]
YLabel=Normalized events
END PLOT

BEGIN PLOT /QCD_INSTANTON_KKS/jets_mreco_inclusive_over_truth
Title=Inclusive jet-mass migration, $|\eta_j| < 4.5$
XLabel=$M^{\mathrm{reco}}_{I,\mathrm{jets}}/\sqrt{\hat{s}}$
YLabel=Normalized events
END PLOT

BEGIN PLOT /QCD_INSTANTON_KKS/jets_mreco_central_over_truth
Title=Central jet-mass migration, $|\eta_j| < 2.5$
XLabel=$M^{\mathrm{reco,central}}_{I,\mathrm{jets}}/\sqrt{\hat{s}}$
YLabel=Normalized events
END PLOT

BEGIN PLOT /QCD_INSTANTON_KKS/low_tracks_n
Title=Number of tracks, $25 < M^{\mathrm{reco}}_{I,\mathrm{tracks}} < 35$ GeV
XLabel=Number of tracks
YLabel=Normalized events
END PLOT

BEGIN PLOT /QCD_INSTANTON_KKS/low_tracks_st
Title=Track $S_T$, $25 < M^{\mathrm{reco}}_{I,\mathrm{tracks}} < 35$ GeV
XLabel=$S_T$ [GeV]
YLabel=Normalized events
END PLOT

BEGIN PLOT /QCD_INSTANTON_KKS/low_tracks_st_20_30
Title=Track $S_T$, $20 < M^{\mathrm{reco}}_{I,\mathrm{tracks}} < 30$ GeV
XLabel=$S_T$ [GeV]
YLabel=Normalized events
END PLOT

BEGIN PLOT /QCD_INSTANTON_KKS/low_tracks_avg_dphi
Title=Average pairwise track $\Delta\phi$, $25 < M^{\mathrm{reco}}_{I,\mathrm{tracks}} < 35$ GeV
XLabel=$\langle\Delta\phi\rangle$ [rad]
YLabel=Normalized events
END PLOT

BEGIN PLOT /QCD_INSTANTON_KKS/low_tracks_sphericity
Title=Rest-frame track sphericity, $25 < M^{\mathrm{reco}}_{I,\mathrm{tracks}} < 35$ GeV
XLabel=Sphericity
YLabel=Normalized events
END PLOT

BEGIN PLOT /QCD_INSTANTON_KKS/high_jets_n
Title=Number of $|\eta_j| < 4.5$ jets, $320 < M^{\mathrm{reco}}_{I,\mathrm{jets}} < 480$ GeV
XLabel=Number of jets
YLabel=Normalized events
END PLOT

BEGIN PLOT /QCD_INSTANTON_KKS/high_jets_st
Title=$|\eta_j| < 4.5$ jet $S_T$, $320 < M^{\mathrm{reco}}_{I,\mathrm{jets}} < 480$ GeV
XLabel=$S_T$ [GeV]
YLabel=Normalized events
END PLOT

BEGIN PLOT /QCD_INSTANTON_KKS/high_jets_avg_dphi
Title=Average $|\eta_j| < 4.5$ jet $\Delta\phi$, $320 < M^{\mathrm{reco}}_{I,\mathrm{jets}} < 480$ GeV
XLabel=$\langle\Delta\phi\rangle$ [rad]
YLabel=Normalized events
END PLOT

BEGIN PLOT /QCD_INSTANTON_KKS/high_jets_sphericity
Title=Rest-frame $|\eta_j| < 4.5$ jet sphericity, $320 < M^{\mathrm{reco}}_{I,\mathrm{jets}} < 480$ GeV
XLabel=Sphericity
YLabel=Normalized events
END PLOT
