// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Math/Units.hh"
#include "Rivet/Projections/UnstableParticles.hh"
namespace Rivet {


  /// @brief Add a short analysis description here
  class TEST_INSTANTON_MASS : public Analysis {
  public:

    /// Constructor
//    DEFAULT_RIVET_ANALYSIS_CTOR(TEST_INSTANTON_MASS);
    TEST_INSTANTON_MASS()
      : Analysis("TEST_INSTANTON_MASS") {  }

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs;
      declare(fs,"FS");

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
//      FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
//      declare(jetfs, "jets");

      // FinalState of prompt photons and bare muons and electrons in the event
//      PromptFinalState photons(Cuts::abspid == PID::PHOTON);
//      PromptFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);
      FinalState partons(Cuts::abspid <= 5 || Cuts::abspid == 21);
      declare(partons,"PARTON");
      // Dress the prompt bare leptons with prompt photons within dR < 0.1,
      // and apply some fiducial cuts on the dressed leptons
//      Cut lepton_cuts = Cuts::abseta < 2.5 && Cuts::pT > 20*GeV;
//      DressedLeptons dressed_leps(photons, bare_leps, 0.1, lepton_cuts);
//      declare(dressed_leps, "leptons");

      // Missing momentum
//      declare(MissingMomentum(fs), "MET");

      // Book histograms
      // specify custom binning
      book(_h["MassSumParton"], "MassSumParton", 50, 0.0, 500.0);
      book(_h["RapidityInstanton"],"RapidityInstanton",50, -5.0, 5.0);
//      book(_h["EtaInstanton"],"EtaInstanton",50, -1000.0, 1000.0);
//      book(_h["ThetaInstanton"],"ThetaInstanton",50,0,pi);
//      book(_h["YYYY"], "myh2", logspace(20, 1e-2, 1e3));
//      book(_h["ZZZZ"], "myh3", {0.0, 1.0, 2.0, 4.0, 8.0, 16.0});
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
//      book(_h["AAAA"], 1, 1, 1);
//      book(_p["BBBB"], 2, 1, 1);
//      book(_c["CCCC"], 3, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Retrieve dressed leptons, sorted by pT
//      vector<DressedLepton> leptons = apply<DressedLeptons>(event, "leptons").dressedLeptons();

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
//      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30*GeV);

      // Remove all jets within dR < 0.2 of a dressed lepton
//      idiscardIfAnyDeltaRLess(jets, leptons, 0.2);

      // Select jets ghost-associated to B-hadrons with a certain fiducial selection
//      Jets bjets = filter_select(jets, [](const Jet& jet) {
//        return  jet.bTagged(Cuts::pT > 5*GeV && Cuts::abseta < 2.5);
//      });

      // Veto event if there are no b-jets
//      if (bjets.empty())  vetoEvent;

      // Apply a missing-momentum cut
//      if (apply<MissingMomentum>(event, "MET").missingPt() < 30*GeV)  vetoEvent;

      const double weight = 1.0;
      FourMomentum InstantonMom(0,0,0,0);
      double InstantonRapidity = 0;
//      double InstantonEta = 0;
      double InstantonMass = 0;
//      double InstantonTheta = 0;
      const Particles& instantonpartons = apply<FinalState>(event, "PARTON").particlesByPt();
      for(const Particle& p : instantonpartons){
        InstantonMom = InstantonMom + p.momentum();
      }
      InstantonMass = InstantonMom.mass()/GeV;
      InstantonRapidity = InstantonMom.rapidity();
//      InstantonEta = InstantonMom.pseudorapidity();
//      InstantonTheta = InstantonMom.theta();
      // Fill histogram with leading b-jet pT
//      _h["XXXX"]->fill(bjets[0].pT()/GeV);

      _h["MassSumParton"]->fill(InstantonMass,weight);
      _h["RapidityInstanton"]->fill(InstantonRapidity,weight);
//      _h["EtaInstanton"]->fill(InstantonEta);
//      _h["ThetaInstanton"]->fill(InstantonTheta);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
cout<<"Cross section "<<crossSection()<<", picobarn "<<picobarn<<", sumW "<<sumW()<<endl;
      normalize(_h["MassSumParton"],crossSection()/picobarn); // normalize to unity
      normalize(_h["RapidityInstanton"],crossSection()/picobarn);
//      normalize(_h["EtaInstanton"],crossSection()/picobarn/sumW());
//      normalize(_h["ThetaInstanton"],crossSection()/picobarn/sumW());
//      normalize(_h["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in fb (no cuts)
//      scale(_h["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)

    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
//    map<string, Profile1DPtr> _p;
//    map<string, CounterPtr> _c;
    //@}


  };


  DECLARE_RIVET_PLUGIN(TEST_INSTANTON_MASS);

}
