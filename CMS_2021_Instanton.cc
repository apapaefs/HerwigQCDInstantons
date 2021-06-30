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
  class CMS_2021_Instanton : public Analysis {
  public:

    /// Constructor
//    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2021_Instanton);
    CMS_2021_Instanton()
      : Analysis("CMS_2021_Instanton") {  }

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fsfull;
      const FinalState fs(Cuts::abseta < 2.4);
      const FinalState fs05(Cuts::abseta < 2.4 && Cuts::pT > 0.5*GeV);
      const FinalState fs09(Cuts::abseta < 2.4 && Cuts::pT > 0.9*GeV);
      declare(fs, "FS");
      declare(fsfull,"FSFULL");
//      ChargedFinalState cfs(Cuts::abseta < 2.4 && Cuts::pT > 0.5*GeV);
      const ChargedFinalState cfs(Cuts::abseta < 2.4);
      const ChargedFinalState cfs05(Cuts::abseta < 2.4 && Cuts::pT > 0.5*GeV);
      const ChargedFinalState cfs09(Cuts::abseta < 2.4 && Cuts::pT > 0.9*GeV);
      declare(cfs, "CFS");
      declare(UnstableParticles(), "UFS");

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs, "jets");

      FastJets jetfs05(fs05, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs05, "jets05");

      FastJets jetfs09(fs09, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs09, "jets09");

      FastJets jetcfs(cfs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetcfs, "chargedjets");

      FastJets jetcfs05(cfs05, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetcfs05, "chargedjets05");

      FastJets jetcfs09(cfs09, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetcfs09, "chargedjets09");

      // FinalState of prompt photons and bare muons and electrons in the event
//      PromptFinalState photons(Cuts::abspid == PID::PHOTON);
//      PromptFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);

      // Dress the prompt bare leptons with prompt photons within dR < 0.1,
      // and apply some fiducial cuts on the dressed leptons
//      Cut lepton_cuts = Cuts::abseta < 2.5 && Cuts::pT > 20*GeV;
//      DressedLeptons dressed_leps(photons, bare_leps, 0.1, lepton_cuts);
//      declare(dressed_leps, "leptons");

      // Missing momentum
//      declare(MissingMomentum(fs), "MET");

      // Book histograms
      // specify custom binning
//      book(_h["XXXX"], "myh1", 20, 0.0, 100.0);
//      book(_h["YYYY"], "myh2", logspace(20, 1e-2, 1e3));
//      book(_h["ZZZZ"], "myh3", {0.0, 1.0, 2.0, 4.0, 8.0, 16.0});
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
//      book(_h["AAAA"], 1, 1, 1);
//      book(_p["BBBB"], 2, 1, 1);
//      book(_c["CCCC"], 3, 1, 1);
//      book(_hist_T_05    ,2,1,1);
//      book(_hist_TM_05    ,4,1,1);
//      book(_hist_S_05    ,6,1,1);
      book(_h["Thrust"], "Thrust", 20, 0.0, 0.5);
      book(_h["ThrustMajor"], "ThrustMajor", 20, 0.0, 1.0);
      book(_h["Spherocity"], "Spherocity", 20, 0.0, 1.0);
      book(_h["NCharged"],"NCharged",50,0.0,500.0);
      book(_h["PtSumCharged"],"PtSumCharged",50,0.0,1000.0);
      book(_h["PtCharged"],"PtCharged",20,0.0,200.0);
      book(_h["MassSumCharged"],"MassSumCharged",50,0.0,1500.0);
      book(_h["IPCharged"],"IPCharged",50,0.0,500.0);
      book(_h["NDisplacedCharged"],"NDisplacedCharged",40,0.0,80.0);
      book(_h["MassMeanCharged"],"MassMeanCharged",30,0.0,30.0);

      book(_h["Thrust05"], "Thrust05", 20, 0.0, 0.5);
      book(_h["ThrustMajor05"], "ThrustMajor05", 20, 0.0, 1.0);
      book(_h["Spherocity05"], "Spherocity05", 20, 0.0, 1.0);
      book(_h["NCharged05"],"NCharged05",30,0.0,300.0);
      book(_h["PtSumCharged05"],"PtSumCharged05",50,0.0,1000.0);
      book(_h["PtCharged05"],"PtCharged05",20,0.0,200.0);
      book(_h["MassSumCharged05"],"MassSumCharged05",50,0.0,1500.0);
      book(_h["IPCharged05"],"IPCharged05",50,0.0,500.0);
      book(_h["NDisplacedCharged05"],"NDisplacedCharged05",40,0.0,80.0);
      book(_h["MassMeanCharged05"],"MassMeanCharged05",30,0.0,30.0);

      book(_h["Thrust09"], "Thrust09", 20, 0.0, 0.5);
      book(_h["ThrustMajor09"], "ThrustMajor09", 20, 0.0, 1.0);
      book(_h["Spherocity09"], "Spherocity09", 20, 0.0, 1.0);
      book(_h["NCharged09"],"NCharged09",20,0.0,200.0);
      book(_h["PtSumCharged09"],"PtSumCharged09",50,0.0,1000.0);
      book(_h["PtCharged09"],"PtCharged09",35,0.0,350.0);
      book(_h["MassSumCharged09"],"MassSumCharged09",50,0.0,1500.0);
      book(_h["IPCharged09"],"IPCharged09",50,0.0,500.0);
      book(_h["NDisplacedCharged09"],"NDisplacedCharged09",40,0.0,80.0);
      book(_h["NChargedJets"],"NChargedJets",40,0.0,40.0);
      book(_h["NChargedJets05"],"NChargedJets05",40,0.0,40.0);
      book(_h["NChargedJets09"],"NChargedJets09",40,0.0,40.0);
      book(_h["MassMeanCharged09"],"MassMeanCharged09",30,0.0,30.0);

      book(_h["ThrustAll"], "ThrustAll", 20, 0.0, 0.5);
      book(_h["ThrustMajorAll"], "ThrustMajorAll", 20, 0.0, 1.0);
      book(_h["SpherocityAll"], "SpherocityAll", 20, 0.0, 1.0);
      book(_h["NParticle"],"NParticle",50,0.0,1000.0);
      book(_h["PtSumParticle"],"PtSumParticle",30,0.0,1500.0);
      book(_h["PtParticle"],"PtParticle",20,0.0,200.0);
      book(_h["MassSumParticle"],"MassSumParticle",50,0.0,2500.0);
      book(_h["MassSumParticleFulleta"],"MassSumParticleFulleta",50,0.0,2500.0);
      book(_h["IPParticle"],"IPParticle",50,0.0,500.0);
      book(_h["NDisplacedParticle"],"NDisplacedParticle",30,0.0,150.0);
      book(_h["MassMeanParticle"],"MassMeanParticle",30,0.0,30.0);

      book(_h["ThrustAll05"], "ThrustAll05", 20, 0.0, 0.5);
      book(_h["ThrustMajorAll05"], "ThrustMajorAll05", 20, 0.0, 1.0);
      book(_h["SpherocityAll05"], "SpherocityAll05", 20, 0.0, 1.0);
      book(_h["NParticle05"],"NParticle05",40,0.0,400.0);
      book(_h["PtSumParticle05"],"PtSumParticle05",30,0.0,1500.0);
      book(_h["PtParticle05"],"PtParticle05",20,0.0,200.0);
      book(_h["MassSumParticle05"],"MassSumParticle05",50,0.0,2500.0);
      book(_h["IPParticle05"],"IPParticle05",50,0.0,500.0);
      book(_h["NDisplacedParticle05"],"NDisplacedParticle05",30,0.0,150.0);
      book(_h["MassMeanParticle05"],"MassMeanParticle05",30,0.0,30.0);

      book(_h["ThrustAll09"], "ThrustAll09", 20, 0.0, 0.5);
      book(_h["ThrustMajorAll09"], "ThrustMajorAll09", 20, 0.0, 1.0);
      book(_h["SpherocityAll09"], "SpherocityAll09", 20, 0.0, 1.0);
      book(_h["NParticle09"],"NParticle09",30,0.0,300.0);
      book(_h["PtSumParticle09"],"PtSumParticle09",30,0.0,1500.0);
      book(_h["PtParticle09"],"PtParticle09",20,0.0,200.0);
      book(_h["MassSumParticle09"],"MassSumParticle09",50,0.0,2500.0);
      book(_h["IPParticle09"],"IPParticle09",50,0.0,500.0);
      book(_h["NDisplacedParticle09"],"NDisplacedParticle09",30,0.0,150.0);
      book(_h["MassMeanParticle09"],"MassMeanParticle09",30,0.0,30.0);

      book(_h["NChargedJets"],"NChargedJets",40,0.0,40.0);
      book(_h["NChargedJets05"],"NChargedJets05",40,0.0,40.0);
      book(_h["NChargedJets09"],"NChargedJets09",40,0.0,40.0);
      book(_h["NJets"],"NJets",50,0.0,50.0);
      book(_h["NJets05"],"NJets05",50,0.0,50.0);
      book(_h["NJets09"],"NJets09",50,0.0,50.0);
      book(_h["PtJets"],"PtJets",40,0.0,500.0);
      book(_h["PtJets05"],"PtJets05",40,0.0,500.0);
      book(_h["PtJets09"],"PtJets09",40,0.0,500.0);
      book(_h["PtChargedJets"],"PtChargedJets",40,0.0,500.0);
      book(_h["PtChargedJets05"],"PtChargedJets05",40,0.0,500.0);
      book(_h["PtChargedJets09"],"PtChargedJets09",40,0.0,500.0);

      book(_h["NBHadron"],"NBHadron",20,0.0,20.0);
      book(_h["NCHadron"],"NCHadron",20,0.0,20.0);
      book(_h["NBCHadron"],"NBCHadron",40,0.0,40.0);
      book(_h["NCHadronPrompt"],"NCHadronPrompt",20,0.0,20.0); // C Hadron not from b-quark
      book(_h["DeltaPhiBB"],"DeltaPhiBB",20,0.0,pi);
      book(_h["DeltaPhiCC"],"DeltaPhiCC",20,0.0,pi);
      book(_h["DeltaRBB"],"DeltaRBB",20,0.0,10);
      book(_h["DeltaRCC"],"DeltaRCC",20,0.0,10);
      book(_h["DeltaEtaBB"],"DeltaEtaBB",20,0.0,5);
      book(_h["DeltaEtaCC"],"DeltaEtaCC",20,0.0,5);

      book(_h["NBChargeddecayFS"],"NBChargeddecayFS",40,0.0,40.0);
      book(_h["NCChargeddecayFS"],"NCChargeddecayFS",40,0.0,40.0);
      book(_h["NBdecayDisplacedCharged"],"NBdecayDisplacedCharged",40,0.0,40.0);
      book(_h["NCdecayDisplacedCharged"],"NCdecayDisplacedCharged",40,0.0,40.0);
      book(_h["IPBdecayCharged"],"IPBdecayCharged",50,0.0,500.0);
      book(_h["IPCdecayCharged"],"IPCdecayCharged",50,0.0,500.0);

      book(_h["NBdecayFS"],"NBdecayFS",50,0.0,100.0);
      book(_h["NCdecayFS"],"NCdecayFS",50,0.0,100.0);

      book(_h["Spherocity09_verylow"], "Spherocity09_verylow", 20, 0.0, 1.0);
      book(_h["Spherocity09_low"], "Spherocity09_low", 20, 0.0, 1.0);
      book(_h["Spherocity09_medium"], "Spherocity09_medium", 20, 0.0, 1.0);
      book(_h["Spherocity09_high"], "Spherocity09_high", 20, 0.0, 1.0);
      book(_h["Spherocity09_veryhigh"], "Spherocity09_veryhigh", 20, 0.0, 1.0);

      book(_h["Thrust09_verylow"], "Thrust09_verylow", 20, 0.0, 0.5);
      book(_h["Thrust09_low"], "Thrust09_low", 20, 0.0, 0.5);
      book(_h["Thrust09_medium"], "Thrust09_medium", 20, 0.0, 0.5);
      book(_h["Thrust09_high"], "Thrust09_high", 20, 0.0, 0.5);
      book(_h["Thrust09_veryhigh"], "Thrust09_veryhigh", 20, 0.0, 0.5);

      book(_h["NCharged09_verylow"],"NCharged09_verylow",30,0.0,150.0);
      book(_h["NCharged09_low"],"NCharged09_low",30,0.0,150.0);
      book(_h["NCharged09_medium"],"NCharged09_medium",30,0.0,150.0);
      book(_h["NCharged09_high"],"NCharged09_high",30,0.0,150.0);
      book(_h["NCharged09_veryhigh"],"NCharged09_veryhigh",30,0.0,150.0);

      book(_h["NDisplacedCharged09_verylow"],"NDisplacedCharged09_verylow",40,0.0,40.0);
      book(_h["NDisplacedCharged09_low"],"NDisplacedCharged09_low",40,0.0,40.0);
      book(_h["NDisplacedCharged09_medium"],"NDisplacedCharged09_medium",40,0.0,40.0);
      book(_h["NDisplacedCharged09_high"],"NDisplacedCharged09_high",40,0.0,40.0);
      book(_h["NDisplacedCharged09_veryhigh"],"NDisplacedCharged09_veryhigh",40,0.0,40.0);

      book(_h["MassMeanCharged09_verylow"],"MassMeanCharged09_verylow",30,0.0,30.0);
      book(_h["MassMeanCharged09_low"],"MassMeanCharged09_low",30,0.0,30.0);
      book(_h["MassMeanCharged09_medium"],"MassMeanCharged09_medium",30,0.0,30.0);
      book(_h["MassMeanCharged09_high"],"MassMeanCharged09_high",30,0.0,30.0);
      book(_h["MassMeanCharged09_veryhigh"],"MassMeanCharged09_veryhigh",30,0.0,30.0);

      book(_h["NChargedJets09_verylow"],"NChargedJets09_verylow",40,0.0,40.0);
      book(_h["NChargedJets09_low"],"NChargedJets09_low",40,0.0,40.0);
      book(_h["NChargedJets09_medium"],"NChargedJets09_medium",40,0.0,40.0);
      book(_h["NChargedJets09_high"],"NChargedJets09_high",40,0.0,40.0);
      book(_h["NChargedJets09_veryhigh"],"NChargedJets09_veryhigh",40,0.0,40.0);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = 1.0;

      // Retrieve dressed leptons, sorted by pT
//      vector<DressedLepton> leptons = apply<DressedLeptons>(event, "leptons").dressedLeptons();

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 3.0*GeV);
      Jets jets05 = apply<FastJets>(event, "jets05").jetsByPt(Cuts::pT > 3.0*GeV);
      Jets jets09 = apply<FastJets>(event, "jets09").jetsByPt(Cuts::pT > 3.0*GeV);
      Jets chargedjets = apply<FastJets>(event, "chargedjets").jetsByPt(Cuts::pT > 3.0*GeV);
      Jets chargedjets05 = apply<FastJets>(event, "chargedjets05").jetsByPt(Cuts::pT > 3.0*GeV);
      Jets chargedjets09 = apply<FastJets>(event, "chargedjets09").jetsByPt(Cuts::pT > 3.0*GeV);
      for (const Jet& j : jets){
        _h["PtJets"]->fill(j.pT(),weight);
      }
      for (const Jet& j : jets05){
        _h["PtJets05"]->fill(j.pT(),weight);
      }
      for (const Jet& j : jets09){
        _h["PtJets09"]->fill(j.pT(),weight);
      }
      for (const Jet& cj : chargedjets){
        _h["PtChargedJets"]->fill(cj.pT(),weight);
      }
      for (const Jet& cj : chargedjets05){
        _h["PtChargedJets05"]->fill(cj.pT(),weight);
      }
      for (const Jet& cj : chargedjets09){
        _h["PtChargedJets09"]->fill(cj.pT(),weight);
      }

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

      // Fill histogram with leading b-jet pT
//      _h["XXXX"]->fill(bjets[0].pT()/GeV);

      const Particles& chargedparticles = apply<ChargedFinalState>(event, "CFS").particlesByPt();
//      if (particles.size() < 2) vetoEvent;
      vector<Vector3> momenta;
      double num = 0;
      double ptSum = 0;
      double IP;
      FourMomentum Charged4mom(0,0,0,0);
      double massSum = 0;
      int numDisplace = 0;

      vector<Vector3> momenta05;
      double num05 = 0;
      double ptSum05 = 0;
      FourMomentum Charged4mom05(0,0,0,0);
      double massSum05 = 0;
      int numDisplace05 = 0;

      vector<Vector3> momenta09;
      double num09 = 0;
      double ptSum09 = 0;
      FourMomentum Charged4mom09(0,0,0,0);
      double massSum09 = 0;
      int numDisplace09 = 0;

      const Particles& particles = apply<FinalState>(event, "FS").particlesByPt();
      const Particles& particlesfull = apply<FinalState>(event, "FSFULL").particlesByPt();
//      if (particles.size() < 2) vetoEvent;
      vector<Vector3> momentaAll;
      double numAll = 0;
      double ptSumAll = 0;
      double IPAll;
      FourMomentum All4mom(0,0,0,0);
      double massSumAll = 0;
      int numDisplaceAll = 0;

      vector<Vector3> momentaAll05;
      double numAll05 = 0;
      double ptSumAll05 = 0;
      FourMomentum All4mom05(0,0,0,0);
      double massSumAll05 = 0;
      int numDisplaceAll05 = 0;

      vector<Vector3> momentaAll09;
      double numAll09 = 0;
      double ptSumAll09 = 0;
      FourMomentum All4mom09(0,0,0,0);
      double massSumAll09 = 0;
      int numDisplaceAll09 = 0;

      FourMomentum All4momfulleta(0,0,0,0);
      double massSumAllfulleta = 0;

      int numBHadron = 0;
      int numCHadron = 0;
      int numCHadronPrompt = 0;

      int numBdecayFS = 0;
      int numCdecayFS = 0;

      int numBdecayChargedFS = 0;
      int numCdecayChargedFS = 0;

      int numBdecayDisplacedCharged = 0;
      int numCdecayDisplacedCharged = 0;

      Particles bhadrons;

      const UnstableParticles &ufs = apply<UnstableFinalState> (event, "UFS");
      for(const Particle& p: ufs.particles()) {

        if (!( p.isHadron() && p.hasBottom()) ) continue;

        bool hasBdaughter = false;
        for(const Particle& pp : p.children()){
          if (pp.isHadron() && pp.hasBottom()) {
            hasBdaughter = true;
          }
        }
        if (hasBdaughter) continue;

        bhadrons += p ;
      }

      numBHadron = bhadrons.size();
      if(numBHadron>=2){

        std::sort(bhadrons.begin(), bhadrons.end(),
            [](const Particle& lhs, const Particle& rhs){
              return lhs.pT() > rhs.pT();
        });
        _h["DeltaPhiBB"]->fill(deltaPhi(bhadrons[0].momentum(),bhadrons[1].momentum()));
        _h["DeltaRBB"]->fill(deltaR(bhadrons[0].momentum(),bhadrons[1].momentum()));
        _h["DeltaEtaBB"]->fill(deltaEta(bhadrons[0].momentum(),bhadrons[1].momentum()));
      }
      Particles chadrons;
      Particles chadronsprompt;
/*      for(ConstGenParticlePtr p: HepMCUtils::particles(event.genEvent())) {

        if (!( PID::isHadron( p->pdg_id() ) && PID::hasCharm( p->pdg_id() )) ) continue;

        ConstGenVertexPtr dv = p->end_vertex();
        bool hasCdaughter = false;
        if ( PID::isHadron( p->pdg_id() ) && PID::hasCharm( p->pdg_id() )) { // c-hadron selection
          if (dv) {            
              for(ConstGenParticlePtr pp: HepMCUtils::particles(dv, Relatives::CHILDREN)){
              if ( PID::isHadron( pp->pdg_id() ) && PID::hasCharm( pp->pdg_id()) ) {
                hasCdaughter = true;
              }
            }
          }
        }
        if (hasCdaughter) continue;

        chadrons += Particle(*p);
      }
*/

      for(const Particle& p: ufs.particles()) {
        
        if (!( p.isHadron() && p.hasCharm()) ) continue;

        bool hasCdaughter = false;
        for(const Particle& pp : p.children()){
          if (pp.isHadron() && pp.hasCharm()) {
            hasCdaughter = true;
          }
        }
        if (hasCdaughter) continue;

        chadrons += p ;
      }
      numCHadron = chadrons.size();

      for(const Particle& p: ufs.particles()) {

        if (!( p.isHadron() && p.hasCharm() && !p.fromBottom()) ) continue;

        bool hasCdaughter = false;
        for(const Particle& pp : p.children()){
          if (pp.isHadron() && pp.hasCharm()) {
            hasCdaughter = true;
          }
        }
        if (hasCdaughter) continue;

        chadronsprompt += p ;
      }
      numCHadronPrompt = chadronsprompt.size();

      if(numCHadronPrompt>=2){

        std::sort(chadronsprompt.begin(), chadronsprompt.end(),
           [](const Particle& lhs, const Particle& rhs){
              return lhs.pT() > rhs.pT();
        });
        _h["DeltaPhiCC"]->fill(deltaPhi(chadronsprompt[0].momentum(),chadronsprompt[1].momentum()));
        _h["DeltaRCC"]->fill(deltaR(chadronsprompt[0].momentum(),chadronsprompt[1].momentum()));
        _h["DeltaEtaCC"]->fill(deltaEta(chadronsprompt[0].momentum(),chadronsprompt[1].momentum()));
      }
      for (const Particle& p : chargedparticles) {
        num += 1;
        ptSum += p.pT()/GeV;
        Vector3 mom = p.p3();
//        Vector3 IPvector = p.closestApproach();  NOTE BY WJ: the closestApproach() function is class Particle is wrong !
//        IP = IPvector.polarRadius();
        IP = getPVDCA(p);

        if(p.fromBottom()){
          numBdecayChargedFS++;
          _h["IPBdecayCharged"]->fill(IP/mm,weight);
          if(IP/mm>0.2) numBdecayDisplacedCharged++;
        }
        if(p.fromCharm()&&!p.fromBottom()){
          numCdecayChargedFS++;
          _h["IPCdecayCharged"]->fill(IP/mm,weight);
          if(IP/mm>0.2) numCdecayDisplacedCharged++;
        }
        if(IP/mm>0.2) numDisplace++;
        Charged4mom = Charged4mom + p.momentum();
        momenta.push_back(mom);
        _h["PtCharged"]->fill(p.pT()/GeV,weight);
        _h["IPCharged"]->fill(IP/mm,weight);

        if(p.pT()/GeV <= 0.5) continue;
        num05 += 1;
        ptSum05 += p.pT()/GeV;
        if(IP/mm>0.2) numDisplace05++;
        Charged4mom05 = Charged4mom05 + p.momentum();
        momenta05.push_back(mom);
        _h["PtCharged05"]->fill(p.pT()/GeV,weight);
        _h["IPCharged05"]->fill(IP/mm,weight);

        if(p.pT()/GeV <= 0.9) continue;
        num09 += 1;
        ptSum09 += p.pT()/GeV;
        if(IP/mm>0.2) numDisplace09++;
        Charged4mom09 = Charged4mom09 + p.momentum();
        momenta09.push_back(mom);
        _h["PtCharged09"]->fill(p.pT()/GeV,weight);
        _h["IPCharged09"]->fill(IP/mm,weight);
      }
      massSum = Charged4mom.mass()/GeV;
      massSum05 = Charged4mom05.mass()/GeV;
      massSum09 = Charged4mom09.mass()/GeV;
//      if (momenta.size() == 2) {
//        momenta.push_back(Vector3(1e-10*MeV, 0., 0.));
//      }

      Thrust thrust;
      thrust.calc(momenta);

      const double T  = 1.0 - thrust.thrust();
      const double TM = thrust.thrustMajor();

      Sphericity sphericity;
      sphericity.calc(momenta);

      double S = sphericity.sphericity();
      if ( std::isnan(S) )  S = -1.0; 

      Thrust thrust05;
      thrust05.calc(momenta05);

      const double T05  = 1.0 - thrust05.thrust();
      const double TM05 = thrust05.thrustMajor();

      Sphericity sphericity05;
      sphericity05.calc(momenta05);

      double S05 = sphericity05.sphericity();
      if ( std::isnan(S05) )  S05 = -1.0;

      Thrust thrust09;
      thrust09.calc(momenta09);

      const double T09  = 1.0 - thrust09.thrust();
      const double TM09 = thrust09.thrustMajor();

      Sphericity sphericity09;
      sphericity09.calc(momenta09);

      double S09 = sphericity09.sphericity();
      if ( std::isnan(S09) )  S09 = -1.0;


      for (const Particle& p : particles) {
        numAll += 1;
        ptSumAll += p.pT()/GeV;
        Vector3 mom = p.p3();
//        Vector3 IPvector = p.closestApproach();
//        IPAll = IPvector.polarRadius();
        IPAll = getPVDCA(p);
        if(IPAll/mm>0.2) numDisplaceAll++;

        if(p.fromBottom()) numBdecayFS++;
        if(p.fromCharm()&&!p.fromBottom()) numCdecayFS++;

        All4mom = All4mom + p.momentum();
        momentaAll.push_back(mom);
        _h["PtParticle"]->fill(p.pT()/GeV,weight);
        _h["IPParticle"]->fill(IPAll/mm,weight);

        if(p.pT()/GeV <= 0.5) continue;
        numAll05 += 1;
        ptSumAll05 += p.pT()/GeV;
        if(IPAll/mm>0.2) numDisplaceAll05++;
        All4mom05 = All4mom05 + p.momentum();
        momentaAll05.push_back(mom);
        _h["PtParticle05"]->fill(p.pT()/GeV,weight);
        _h["IPParticle05"]->fill(IPAll/mm,weight);

        if(p.pT()/GeV <= 0.9) continue;
        numAll09 += 1;
        ptSumAll09 += p.pT()/GeV;
        if(IPAll/mm>0.2) numDisplaceAll09++;
        All4mom09 = All4mom09 + p.momentum();
        momentaAll09.push_back(mom);
        _h["PtParticle09"]->fill(p.pT()/GeV,weight);
        _h["IPParticle09"]->fill(IPAll/mm,weight);
      }
      massSumAll = All4mom.mass()/GeV;
      massSumAll05 = All4mom05.mass()/GeV;
      massSumAll09 = All4mom09.mass()/GeV;

      Thrust thrustAll;
      thrustAll.calc(momentaAll);

      const double TAll  = 1.0 - thrustAll.thrust();
      const double TMAll = thrustAll.thrustMajor();

      Sphericity sphericityAll;
      sphericityAll.calc(momentaAll);

      double SAll = sphericityAll.sphericity();
      if ( std::isnan(SAll) )  SAll = -1.0;

      Thrust thrustAll05;
      thrustAll05.calc(momentaAll05);

      const double TAll05  = 1.0 - thrustAll05.thrust();
      const double TMAll05 = thrustAll05.thrustMajor();

      Sphericity sphericityAll05;
      sphericityAll05.calc(momentaAll05);

      double SAll05 = sphericityAll05.sphericity();
      if ( std::isnan(SAll05) )  SAll05 = -1.0;

      Thrust thrustAll09;
      thrustAll09.calc(momentaAll09);

      const double TAll09  = 1.0 - thrustAll09.thrust();
      const double TMAll09 = thrustAll09.thrustMajor();

      Sphericity sphericityAll09;
      sphericityAll09.calc(momentaAll09);

      double SAll09 = sphericityAll09.sphericity();
      if ( std::isnan(SAll09) )  SAll09 = -1.0;

      for (const Particle& p : particlesfull) {
        All4momfulleta = All4momfulleta + p.momentum();
      }

      massSumAllfulleta = All4momfulleta.mass()/GeV;

      _h["NBdecayFS"]->fill(numBdecayFS,weight);
      _h["NCdecayFS"]->fill(numCdecayFS,weight);
      _h["NBChargeddecayFS"]->fill(numBdecayChargedFS,weight);
      _h["NCChargeddecayFS"]->fill(numCdecayChargedFS,weight);

      _h["NBHadron"]->fill(numBHadron,weight);
      _h["NCHadron"]->fill(numCHadron,weight);
      _h["NCHadronPrompt"]->fill(numCHadronPrompt,weight);
      _h["NBCHadron"]->fill(numBHadron+numCHadronPrompt,weight);
      _h["NBdecayDisplacedCharged"]->fill(numBdecayDisplacedCharged,weight);
      _h["NCdecayDisplacedCharged"]->fill(numCdecayDisplacedCharged,weight);

      _h["NChargedJets"]->fill(chargedjets.size(),weight);
      _h["NJets"]->fill(jets.size(),weight);
      _h["NChargedJets05"]->fill(chargedjets05.size(),weight);
      _h["NJets05"]->fill(jets05.size(),weight);
      _h["NChargedJets09"]->fill(chargedjets09.size(),weight);
      _h["NJets09"]->fill(jets09.size(),weight);

      _h["Thrust"]->fill(T , weight);
      _h["ThrustMajor"]->fill(TM, weight);
      _h["Spherocity"]->fill(S , weight);
      _h["NCharged"]->fill(num,weight);
      _h["PtSumCharged"]->fill(ptSum,weight);
      _h["MassSumCharged"]->fill(massSum,weight);
      _h["NDisplacedCharged"]->fill(numDisplace,weight);
      _h["MassMeanCharged"]->fill(massSum/num,weight);

      _h["Thrust05"]->fill(T05 , weight);
      _h["ThrustMajor05"]->fill(TM05, weight);
      _h["Spherocity05"]->fill(S05 , weight);
      _h["NCharged05"]->fill(num05,weight);
      _h["PtSumCharged05"]->fill(ptSum05,weight);
      _h["MassSumCharged05"]->fill(massSum05,weight);
      _h["NDisplacedCharged05"]->fill(numDisplace05,weight);
      _h["MassMeanCharged05"]->fill(massSum05/num05,weight);

      _h["Thrust09"]->fill(T09 , weight);
      _h["ThrustMajor09"]->fill(TM09, weight);
      _h["Spherocity09"]->fill(S09 , weight);
      _h["NCharged09"]->fill(num09,weight);
      _h["PtSumCharged09"]->fill(ptSum09,weight);
      _h["MassSumCharged09"]->fill(massSum09,weight);
      _h["NDisplacedCharged09"]->fill(numDisplace09,weight);
      _h["MassMeanCharged09"]->fill(massSum09/num09,weight);

      if (massSum09>=20&&massSum09<=40){
        _h["Spherocity09_verylow"]->fill(S09 , weight);
        _h["Thrust09_verylow"]->fill(T09 , weight);
        _h["NCharged09_verylow"]->fill(num09,weight);
        _h["NDisplacedCharged09_verylow"]->fill(numDisplace09,weight);
        _h["MassMeanCharged09_verylow"]->fill(massSum09/num09,weight);
        _h["NChargedJets09_verylow"]->fill(jets09.size(),weight);
      }
      else if (massSum09>40&&massSum09<=80){
        _h["Spherocity09_low"]->fill(S09 , weight);
        _h["Thrust09_low"]->fill(T09 , weight);
        _h["NCharged09_low"]->fill(num09,weight);
        _h["NDisplacedCharged09_low"]->fill(numDisplace09,weight);
        _h["MassMeanCharged09_low"]->fill(massSum09/num09,weight);
        _h["NChargedJets09_low"]->fill(jets09.size(),weight);
      }
      else if (massSum09>80&&massSum09<=200){
        _h["Spherocity09_medium"]->fill(S09 , weight);
        _h["Thrust09_medium"]->fill(T09 , weight);
        _h["NCharged09_medium"]->fill(num09,weight);
        _h["NDisplacedCharged09_medium"]->fill(numDisplace09,weight);
        _h["MassMeanCharged09_medium"]->fill(massSum09/num09,weight);
        _h["NChargedJets09_medium"]->fill(jets09.size(),weight);
      }
      else if (massSum09>200&&massSum09<=300){
        _h["Spherocity09_high"]->fill(S09 , weight);
        _h["Thrust09_high"]->fill(T09 , weight);
        _h["NCharged09_high"]->fill(num09,weight);
        _h["NDisplacedCharged09_high"]->fill(numDisplace09,weight);
        _h["MassMeanCharged09_high"]->fill(massSum09/num09,weight);
        _h["NChargedJets09_high"]->fill(jets09.size(),weight);
      }
      else if (massSum09>300&&massSum09<=500){
        _h["Spherocity09_veryhigh"]->fill(S09 , weight);
        _h["Thrust09_veryhigh"]->fill(T09 , weight);
        _h["NCharged09_veryhigh"]->fill(num09,weight);
        _h["NDisplacedCharged09_veryhigh"]->fill(numDisplace09,weight);
        _h["MassMeanCharged09_veryhigh"]->fill(massSum09/num09,weight);
        _h["NChargedJets09_veryhigh"]->fill(jets09.size(),weight);
      }


      _h["ThrustAll"]->fill(TAll , weight);
      _h["ThrustMajorAll"]->fill(TMAll, weight);
      _h["SpherocityAll"]->fill(SAll , weight);
      _h["NParticle"]->fill(numAll,weight);
      _h["PtSumParticle"]->fill(ptSumAll,weight);
      _h["MassSumParticle"]->fill(massSumAll,weight);
      _h["MassSumParticleFulleta"]->fill(massSumAllfulleta,weight);
      _h["NDisplacedParticle"]->fill(numDisplaceAll,weight);
      _h["MassMeanParticle"]->fill(massSumAll/numAll,weight);

      _h["ThrustAll05"]->fill(TAll05 , weight);
      _h["ThrustMajorAll05"]->fill(TMAll05, weight);
      _h["SpherocityAll05"]->fill(SAll05 , weight);
      _h["NParticle05"]->fill(numAll05,weight);
      _h["PtSumParticle05"]->fill(ptSumAll05,weight);
      _h["MassSumParticle05"]->fill(massSumAll05,weight);
      _h["NDisplacedParticle05"]->fill(numDisplaceAll05,weight);
      _h["MassMeanParticle05"]->fill(massSumAll05/numAll05,weight);

      _h["ThrustAll09"]->fill(TAll09 , weight);
      _h["ThrustMajorAll09"]->fill(TMAll09, weight);
      _h["SpherocityAll09"]->fill(SAll09 , weight);
      _h["NParticle09"]->fill(numAll09,weight);
      _h["PtSumParticle09"]->fill(ptSumAll09,weight);
      _h["MassSumParticle09"]->fill(massSumAll09,weight);
      _h["NDisplacedParticle09"]->fill(numDisplaceAll09,weight);
      _h["MassMeanParticle09"]->fill(massSumAll09/numAll09,weight);
   }


    /// Normalise histograms etc., after the run
    void finalize() {

 //     normalize(_h["XXXX"]); // scale to unity
//      normalize(_h["YYYY"], crossSection()/picobarn); // scale to generated cross-section in fb (no cuts)
//      scale(_h["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)
      scale(_h["Thrust"],crossSection()/picobarn/sumW());
      scale(_h["ThrustMajor"],crossSection()/picobarn/sumW());
      scale(_h["Spherocity"],crossSection()/picobarn/sumW());
      scale(_h["NCharged"],crossSection()/picobarn/sumW());
      scale(_h["PtCharged"],crossSection()/picobarn/sumW());
      scale(_h["IPCharged"],crossSection()/picobarn/sumW());
      scale(_h["PtSumCharged"],crossSection()/picobarn/sumW());
      scale(_h["MassSumCharged"],crossSection()/picobarn/sumW());
      scale(_h["NDisplacedCharged"],crossSection()/picobarn/sumW());
      scale(_h["MassMeanCharged"],crossSection()/picobarn/sumW());

      scale(_h["Thrust05"],crossSection()/picobarn/sumW());
      scale(_h["ThrustMajor05"],crossSection()/picobarn/sumW());
      scale(_h["Spherocity05"],crossSection()/picobarn/sumW());
      scale(_h["NCharged05"],crossSection()/picobarn/sumW());
      scale(_h["PtCharged05"],crossSection()/picobarn/sumW());
      scale(_h["IPCharged05"],crossSection()/picobarn/sumW());
      scale(_h["PtSumCharged05"],crossSection()/picobarn/sumW());
      scale(_h["MassSumCharged05"],crossSection()/picobarn/sumW());
      scale(_h["NDisplacedCharged05"],crossSection()/picobarn/sumW());
      scale(_h["MassMeanCharged05"],crossSection()/picobarn/sumW());

      scale(_h["Thrust09"],crossSection()/picobarn/sumW());
      scale(_h["ThrustMajor09"],crossSection()/picobarn/sumW());
      scale(_h["Spherocity09"],crossSection()/picobarn/sumW());
      scale(_h["NCharged09"],crossSection()/picobarn/sumW());
      scale(_h["PtCharged09"],crossSection()/picobarn/sumW());
      scale(_h["IPCharged09"],crossSection()/picobarn/sumW());
      scale(_h["PtSumCharged09"],crossSection()/picobarn/sumW());
      scale(_h["MassSumCharged09"],crossSection()/picobarn/sumW());
      scale(_h["NDisplacedCharged09"],crossSection()/picobarn/sumW());
      scale(_h["MassMeanCharged09"],crossSection()/picobarn/sumW());

      scale(_h["Spherocity09_verylow"],crossSection()/picobarn/sumW());
      scale(_h["Thrust09_verylow"],crossSection()/picobarn/sumW());
      scale(_h["NCharged09_verylow"],crossSection()/picobarn/sumW());
      scale(_h["NDisplacedCharged09_verylow"],crossSection()/picobarn/sumW());
      scale(_h["MassMeanCharged09_verylow"],crossSection()/picobarn/sumW());
      scale(_h["NChargedJets09_verylow"],crossSection()/picobarn/sumW());

      scale(_h["Spherocity09_low"],crossSection()/picobarn/sumW());
      scale(_h["Thrust09_low"],crossSection()/picobarn/sumW());
      scale(_h["NCharged09_low"],crossSection()/picobarn/sumW());
      scale(_h["NDisplacedCharged09_low"],crossSection()/picobarn/sumW());
      scale(_h["MassMeanCharged09_low"],crossSection()/picobarn/sumW());
      scale(_h["NChargedJets09_low"],crossSection()/picobarn/sumW());

      scale(_h["Spherocity09_medium"],crossSection()/picobarn/sumW());
      scale(_h["Thrust09_medium"],crossSection()/picobarn/sumW());
      scale(_h["NCharged09_medium"],crossSection()/picobarn/sumW());
      scale(_h["NDisplacedCharged09_medium"],crossSection()/picobarn/sumW());
      scale(_h["MassMeanCharged09_medium"],crossSection()/picobarn/sumW());
      scale(_h["NChargedJets09_medium"],crossSection()/picobarn/sumW());
 
      scale(_h["Spherocity09_high"],crossSection()/picobarn/sumW());
      scale(_h["Thrust09_high"],crossSection()/picobarn/sumW());
      scale(_h["NCharged09_high"],crossSection()/picobarn/sumW());
      scale(_h["NDisplacedCharged09_high"],crossSection()/picobarn/sumW());
      scale(_h["MassMeanCharged09_high"],crossSection()/picobarn/sumW());
      scale(_h["NChargedJets09_high"],crossSection()/picobarn/sumW());

      scale(_h["Spherocity09_veryhigh"],crossSection()/picobarn/sumW());
      scale(_h["Thrust09_veryhigh"],crossSection()/picobarn/sumW());
      scale(_h["NCharged09_veryhigh"],crossSection()/picobarn/sumW());
      scale(_h["NDisplacedCharged09_veryhigh"],crossSection()/picobarn/sumW());
      scale(_h["MassMeanCharged09_veryhigh"],crossSection()/picobarn/sumW());
      scale(_h["NChargedJets09_veryhigh"],crossSection()/picobarn/sumW());

      scale(_h["ThrustAll"],crossSection()/picobarn/sumW());
      scale(_h["ThrustMajorAll"],crossSection()/picobarn/sumW());
      scale(_h["SpherocityAll"],crossSection()/picobarn/sumW());
      scale(_h["NParticle"],crossSection()/picobarn/sumW());
      scale(_h["PtParticle"],crossSection()/picobarn/sumW());
      scale(_h["IPParticle"],crossSection()/picobarn/sumW());
      scale(_h["PtSumParticle"],crossSection()/picobarn/sumW());
      scale(_h["MassSumParticle"],crossSection()/picobarn/sumW());
      scale(_h["MassSumParticleFulleta"],crossSection()/picobarn/sumW());
      scale(_h["NDisplacedParticle"],crossSection()/picobarn/sumW());
      scale(_h["MassMeanParticle"],crossSection()/picobarn/sumW());

      scale(_h["ThrustAll05"],crossSection()/picobarn/sumW());
      scale(_h["ThrustMajorAll05"],crossSection()/picobarn/sumW());
      scale(_h["SpherocityAll05"],crossSection()/picobarn/sumW());
      scale(_h["NParticle05"],crossSection()/picobarn/sumW());
      scale(_h["PtParticle05"],crossSection()/picobarn/sumW());
      scale(_h["IPParticle05"],crossSection()/picobarn/sumW());
      scale(_h["PtSumParticle05"],crossSection()/picobarn/sumW());
      scale(_h["MassSumParticle05"],crossSection()/picobarn/sumW());
      scale(_h["NDisplacedParticle05"],crossSection()/picobarn/sumW());
      scale(_h["MassMeanParticle05"],crossSection()/picobarn/sumW());

      scale(_h["ThrustAll09"],crossSection()/picobarn/sumW());
      scale(_h["ThrustMajorAll09"],crossSection()/picobarn/sumW());
      scale(_h["SpherocityAll09"],crossSection()/picobarn/sumW());
      scale(_h["NParticle09"],crossSection()/picobarn/sumW());
      scale(_h["PtParticle09"],crossSection()/picobarn/sumW());
      scale(_h["IPParticle09"],crossSection()/picobarn/sumW());
      scale(_h["PtSumParticle09"],crossSection()/picobarn/sumW());
      scale(_h["MassSumParticle09"],crossSection()/picobarn/sumW());
      scale(_h["NDisplacedParticle09"],crossSection()/picobarn/sumW());
      scale(_h["MassMeanParticle09"],crossSection()/picobarn/sumW());

      scale(_h["NChargedJets"],crossSection()/picobarn/sumW());
      scale(_h["NJets"],crossSection()/picobarn/sumW());
      scale(_h["NChargedJets05"],crossSection()/picobarn/sumW());
      scale(_h["NJets05"],crossSection()/picobarn/sumW());
      scale(_h["NChargedJets09"],crossSection()/picobarn/sumW());
      scale(_h["NJets09"],crossSection()/picobarn/sumW());
      scale(_h["PtChargedJets"],crossSection()/picobarn/sumW());
      scale(_h["PtChargedJets05"],crossSection()/picobarn/sumW()); 
      scale(_h["PtChargedJets09"],crossSection()/picobarn/sumW());
      scale(_h["PtJets"],crossSection()/picobarn/sumW());
      scale(_h["PtJets05"],crossSection()/picobarn/sumW());
      scale(_h["PtJets09"],crossSection()/picobarn/sumW());

      scale(_h["NBdecayFS"],crossSection()/picobarn/sumW());
      scale(_h["NCdecayFS"],crossSection()/picobarn/sumW());
      scale(_h["NBChargeddecayFS"],crossSection()/picobarn/sumW());
      scale(_h["NCChargeddecayFS"],crossSection()/picobarn/sumW());
      scale(_h["NBdecayDisplacedCharged"],crossSection()/picobarn/sumW());
      scale(_h["NCdecayDisplacedCharged"],crossSection()/picobarn/sumW());
      scale(_h["IPBdecayCharged"],crossSection()/picobarn/sumW());
      scale(_h["IPCdecayCharged"],crossSection()/picobarn/sumW());

      scale(_h["NBHadron"],crossSection()/picobarn/sumW());
      scale(_h["NCHadron"],crossSection()/picobarn/sumW());
      scale(_h["NBCHadron"],crossSection()/picobarn/sumW());
      scale(_h["NCHadronPrompt"],crossSection()/picobarn/sumW());
      scale(_h["DeltaPhiBB"],crossSection()/picobarn/sumW());
      scale(_h["DeltaPhiCC"],crossSection()/picobarn/sumW());
      scale(_h["DeltaRBB"],crossSection()/picobarn/sumW());
      scale(_h["DeltaRCC"],crossSection()/picobarn/sumW());
      scale(_h["DeltaEtaBB"],crossSection()/picobarn/sumW());
      scale(_h["DeltaEtaCC"],crossSection()/picobarn/sumW());
    }

    //@}


    /// @name Histograms
    //@{
 //   map<string, Histo1DPtr> _h;
//    map<string, Profile1DPtr> _p;
//    map<string, CounterPtr> _c;
    //@}

    private:

    double getPVDCA(const Particle& p) {
      ConstGenVertexPtr vtx = p.genParticle()->production_vertex();
      if ( 0 == vtx ) return -1.;
      const Vector3 u = p.momentum().p3().unit();
      Vector3 d(vtx->position().x(), vtx->position().y(), vtx->position().z());
      double proj = d.dot(u);
      d -= (u * proj);
      return d.mod();

    }
    map<string, Histo1DPtr> _h;
   
  };

  DECLARE_RIVET_PLUGIN(CMS_2021_Instanton);

}
