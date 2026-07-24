// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Math/LorentzTrans.hh"

#include <cmath>

namespace Rivet {


  /// Particle-level version of the QCD instanton observables in arXiv:1911.09726, Sec. 4.
  class QCD_INSTANTON_KKS : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(QCD_INSTANTON_KKS);


    void init() {
      const ChargedFinalState tracks(Cuts::pT > 0.5*GeV && Cuts::abseta < 2.5);
      declare(tracks, "Tracks");

      const FinalState fs(Cuts::abseta < 5.0);
      const FastJets jets(fs, JetAlg::ANTIKT, 0.4, JetMuons::NONE, JetInvisibles::NONE);
      declare(jets, "AntiKt04Jets");

      book(_h["tracks_mreco"], "tracks_mreco", 60, 0.0, 120.0);

      std::vector<double> jetMassEdges;
      for (double edge = 0.0; edge <= 800.0; edge += 10.0) {
        jetMassEdges.push_back(edge);
      }
      for (double edge = 850.0; edge <= 3000.0; edge += 50.0) {
        jetMassEdges.push_back(edge);
      }
      book(
        _h["jets_mreco_inclusive_eta45"],
        "jets_mreco_inclusive_eta45",
        jetMassEdges
      );
      book(
        _h["jets_mreco_central"],
        "jets_mreco_central",
        jetMassEdges
      );
      book(
        _h["instanton_mass_truth"],
        "instanton_mass_truth",
        jetMassEdges
      );
      book(
        _h["jets_mreco_inclusive_over_truth"],
        "jets_mreco_inclusive_over_truth",
        100, 0.0, 10.0
      );
      book(
        _h["jets_mreco_central_over_truth"],
        "jets_mreco_central_over_truth",
        100, 0.0, 10.0
      );

      book(_h["low_tracks_n"], "low_tracks_n", 40, -0.5, 39.5);
      book(_h["low_tracks_st"], "low_tracks_st", 45, 0.0, 45.0);
      book(_h["low_tracks_st_20_30"], "low_tracks_st_20_30", 45, 0.0, 45.0);
      book(_h["low_tracks_avg_dphi"], "low_tracks_avg_dphi", 50, 0.0, M_PI);
      book(_h["low_tracks_sphericity"], "low_tracks_sphericity", 50, 0.0, 1.0);

      book(_h["high_jets_n"], "high_jets_n", 25, -0.5, 24.5);
      book(_h["high_jets_st"], "high_jets_st", 60, 0.0, 600.0);
      book(_h["high_jets_avg_dphi"], "high_jets_avg_dphi", 50, 0.0, M_PI);
      book(_h["high_jets_sphericity"], "high_jets_sphericity", 50, 0.0, 1.0);

      book(_c["all"], "_all");
      book(_c["low_tracks_25_35"], "_low_tracks_25_35");
      book(_c["low_tracks_20_30"], "_low_tracks_20_30");
      book(_c["high_jets_320_480"], "_high_jets_320_480");
      book(_c["truth_mass_valid"], "_truth_mass_valid");
    }


    void analyze(const Event& event) {
      _c["all"]->fill();

      const Particles tracks = apply<ChargedFinalState>(event, "Tracks").particlesByPt();
      const FastJets& jetProjection = apply<FastJets>(event, "AntiKt04Jets");
      const Jets jetsInclusive = jetProjection.jetsByPt(
        Cuts::pT > 20*GeV && Cuts::abseta < 4.5
      );
      const Jets jetsCentral = jetProjection.jetsByPt(
        Cuts::pT > 20*GeV && Cuts::abseta < 2.5
      );

      const double tracksMreco = invariantMass(tracks);
      const double jetsMrecoInclusive = invariantMass(jetsInclusive);
      const double jetsMrecoCentral = invariantMass(jetsCentral);
      const double truthMass = instantonMassTruth(event);
      if (std::isfinite(tracksMreco) && tracksMreco >= 0.0) {
        _h["tracks_mreco"]->fill(tracksMreco);
      }
      if (std::isfinite(jetsMrecoInclusive) && jetsMrecoInclusive >= 0.0) {
        _h["jets_mreco_inclusive_eta45"]->fill(jetsMrecoInclusive);
      }
      if (std::isfinite(jetsMrecoCentral) && jetsMrecoCentral >= 0.0) {
        _h["jets_mreco_central"]->fill(jetsMrecoCentral);
      }
      if (std::isfinite(truthMass) && truthMass > 0.0) {
        _c["truth_mass_valid"]->fill();
        _h["instanton_mass_truth"]->fill(truthMass);
        fillFiniteRatio(
          "jets_mreco_inclusive_over_truth",
          jetsMrecoInclusive,
          truthMass
        );
        fillFiniteRatio(
          "jets_mreco_central_over_truth",
          jetsMrecoCentral,
          truthMass
        );
      }

      if (tracksMreco > 25.0 && tracksMreco < 35.0) {
        _c["low_tracks_25_35"]->fill();
        _h["low_tracks_n"]->fill(tracks.size());
        _h["low_tracks_st"]->fill(scalarPtSum(tracks));
        fillPairAndShapeObservables(tracks, "low_tracks");
      }

      if (tracksMreco > 20.0 && tracksMreco < 30.0) {
        _c["low_tracks_20_30"]->fill();
        _h["low_tracks_st_20_30"]->fill(scalarPtSum(tracks));
      }

      // Retain the published high-mass selection on the original eta < 4.5
      // jet collection. The central mass is an additional shower diagnostic.
      if (jetsMrecoInclusive > 320.0 && jetsMrecoInclusive < 480.0) {
        _c["high_jets_320_480"]->fill();
        _h["high_jets_n"]->fill(jetsInclusive.size());
        _h["high_jets_st"]->fill(scalarPtSum(jetsInclusive));
        fillPairAndShapeObservables(jetsInclusive, "high_jets");
      }
    }


    void finalize() {
      for (auto& hist : _h) {
        const double integral = hist.second->integral();
        if (std::isfinite(integral) && integral > 0.0) normalize(hist.second);
      }
    }


  private:

    double instantonMassTruth(const Event& event) const {
      const GenEvent* genEvent = event.genEvent();
      if (genEvent == nullptr || genEvent->pdf_info() == nullptr) return -1.0;

      // Both generators store the hard-process momentum fractions in HepMC:
      // sqrt(x1*x2*s) equals the Herwig 2->N mass and Sherpa's PDG-999 mass.
      const PdfInfo pdfInfo = *(genEvent->pdf_info());
      const double x1 = pdfInfo.x[0];
      const double x2 = pdfInfo.x[1];
      const double colliderEnergy = event.sqrtS();
      if (!std::isfinite(x1) || !std::isfinite(x2)
          || !std::isfinite(colliderEnergy)
          || x1 <= 0.0 || x2 <= 0.0 || colliderEnergy <= 0.0) {
        return -1.0;
      }

      const double mass = colliderEnergy*std::sqrt(x1*x2)/GeV;
      return std::isfinite(mass) && mass > 0.0 ? mass : -1.0;
    }


    void fillFiniteRatio(
      const string& histogram,
      double numerator,
      double denominator
    ) {
      if (!std::isfinite(numerator) || numerator < 0.0
          || !std::isfinite(denominator) || denominator <= 0.0) {
        return;
      }
      const double ratio = numerator/denominator;
      if (std::isfinite(ratio) && ratio >= 0.0) {
        _h[histogram]->fill(ratio);
      }
    }


    template <typename Objects>
    double invariantMass(const Objects& objects) const {
      FourMomentum total(0.0, 0.0, 0.0, 0.0);
      for (const auto& obj : objects) {
        total += obj.momentum();
      }
      const double mass = total.mass()/GeV;
      return std::isfinite(mass) ? mass : -1.0;
    }


    template <typename Objects>
    double scalarPtSum(const Objects& objects) const {
      double st = 0.0;
      for (const auto& obj : objects) {
        st += obj.pT();
      }
      return st/GeV;
    }


    template <typename Objects>
    double averageDeltaPhi(const Objects& objects) const {
      if (objects.size() < 2) return -1.0;

      double dphiSum = 0.0;
      size_t npairs = 0;
      for (size_t i = 0; i < objects.size(); ++i) {
        for (size_t j = i + 1; j < objects.size(); ++j) {
          dphiSum += std::fabs(deltaPhi(objects[i].phi(), objects[j].phi()));
          ++npairs;
        }
      }
      return npairs > 0 ? dphiSum/npairs : -1.0;
    }


    template <typename Objects>
    double sphericity(const Objects& objects) const {
      if (objects.size() < 2) return -1.0;

      FourMomentum total(0.0, 0.0, 0.0, 0.0);
      for (const auto& obj : objects) total += obj.momentum();
      const double beta2 = total.betaVec().mod2();
      if (!std::isfinite(total.E()) || total.E() <= 0.0
          || !std::isfinite(beta2) || beta2 >= 1.0) return -1.0;
      const LorentzTransform toRestFrame =
        LorentzTransform::mkFrameTransformFromBeta(total.betaVec());

      vector<Vector3> momenta;
      momenta.reserve(objects.size());
      for (const auto& obj : objects) {
        momenta.push_back(toRestFrame.transform(obj.momentum()).p3());
      }

      Sphericity sph;
      sph.calc(momenta);
      const double s = sph.sphericity();
      return std::isfinite(s) ? s : -1.0;
    }


    template <typename Objects>
    void fillPairAndShapeObservables(const Objects& objects, const string& prefix) {
      const double dphi = averageDeltaPhi(objects);
      if (dphi >= 0.0) _h[prefix + "_avg_dphi"]->fill(dphi);

      const double sph = sphericity(objects);
      if (sph >= 0.0) _h[prefix + "_sphericity"]->fill(sph);
    }


    map<string, Histo1DPtr> _h;
    map<string, CounterPtr> _c;

  };


  RIVET_DECLARE_PLUGIN(QCD_INSTANTON_KKS);

}
