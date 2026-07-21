BEGIN PLOT /QCD_INSTANTON_KKS/tracks_mreco
Title=Track-based reconstructed instanton mass proxy
XLabel=$M^{\mathrm{reco}}_{I,\mathrm{tracks}}$ [GeV]
YLabel=Normalized events
END PLOT

BEGIN PLOT /QCD_INSTANTON_KKS/jets_mreco
Title=Jet-based reconstructed instanton mass proxy
XLabel=$M^{\mathrm{reco}}_{I,\mathrm{jets}}$ [GeV]
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
Title=Number of jets, $320 < M^{\mathrm{reco}}_{I,\mathrm{jets}} < 480$ GeV
XLabel=Number of jets
YLabel=Normalized events
END PLOT

BEGIN PLOT /QCD_INSTANTON_KKS/high_jets_st
Title=Jet $S_T$, $320 < M^{\mathrm{reco}}_{I,\mathrm{jets}} < 480$ GeV
XLabel=$S_T$ [GeV]
YLabel=Normalized events
END PLOT

BEGIN PLOT /QCD_INSTANTON_KKS/high_jets_avg_dphi
Title=Average pairwise jet $\Delta\phi$, $320 < M^{\mathrm{reco}}_{I,\mathrm{jets}} < 480$ GeV
XLabel=$\langle\Delta\phi\rangle$ [rad]
YLabel=Normalized events
END PLOT

BEGIN PLOT /QCD_INSTANTON_KKS/high_jets_sphericity
Title=Rest-frame jet sphericity, $320 < M^{\mathrm{reco}}_{I,\mathrm{jets}} < 480$ GeV
XLabel=Sphericity
YLabel=Normalized events
END PLOT
