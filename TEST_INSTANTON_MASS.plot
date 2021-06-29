BEGIN PLOT /TEST_INSTANTON_MASS/MassSumParton
Title= C.O.M energy of the partons from instantons
XLabel=$\sqrt{s_{instanton}}(GeV)$
YLabel=Cross section(pb)
# + any additional plot settings you might like, see make-plots documentation
END PLOT
BEGIN PLOT /TEST_INSTANTON_MASS/RapidityInstanton
XLabel= $Rapidity_{instanton}$
YLabel= Cross section(pb)
Title= Rapidity of instantons
END PLOT

#BEGIN PLOT /TEST_INSTANTON_MASS/EtaInstanton
#XLabel= $\eta_{instanton}$
#YLabel= Cross section(pb)
#Title= Pseudorapidity of instantons
#END PLOT

#BEGIN PLOT /TEST_INSTANTON_MASS/ThetaInstanton
#XLabel= $\theta_{instanton}$
#YLabel= Cross section(pb)
#Title= Polar angle of instantons
#END PLOT
# ... add more histograms as you need them ...
