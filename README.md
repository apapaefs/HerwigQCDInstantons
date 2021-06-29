#Plot Kinematic observables
In `LHC-Instanton.in`, `LHC-HardQCD.in`, `LHC-MB.in`, the `hepmc` generation and plot methods with `rivet` are added. 
To get the histograms of some selected observables after showering, first compile the rivet routine (weight the events to luminosity 1 pb^-1)
```
source sourceRivet.sh
rivet-build CMS_2021_Instanton.cc
```
Then do the Herwig read and run to produce `*.yoda` files containing the histograms and `*.hepmc` files containing the events.
To run the signal-background comparison
```
rivet-mkhtml -o rivet-compare-Instanton-HardQCD-MinBias-mass20 --errs LHC-Instanton-nf4.yoda:'LineColor=orange':'Title=Random2 with UE nf=4' LHC-HardQCD-default.yoda:'LineColor=blue':'Title=HardQCD with UE' LHC-HardQCD-default-noUE.yoda:'LineColor=blue!60!white':'LineStyle=dashed':'Title=HardQCD no UE' LHC-MB-default.yoda:'LineColor=green':'Title=MinBias default' PLOT:'PlotSize=16,10':'RightMargin=6':'LegendXPos=1.01':'LegendYPos=0.9' -t 'Comparison of Instanton, HardQCD and MinBias'
```
The `TEST_INSTANTON_MASS.cc` script can plot the instanton mass spectrum and the rapidity distribution. First comment out line 16 , 17, 18, 20 to turn off the showering, hadronization and multiparton interaction, comment out line 126-129 and instead run line 131-134. Compile and run
```
source sourceRivet.sh
rivet-build TEST_INSTANTON_MASS.cc
```
After Herwig read and run, you can plot the instanton `.yoda` file and get the mass and rapidity distribution.
