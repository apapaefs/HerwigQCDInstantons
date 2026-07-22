// -*- C++ -*-
#ifndef Herwig_MEInstanton_H
#define Herwig_MEInstanton_H

#include "Herwig/MatrixElement/BlobME.h"
#include "Herwig/Utilities/Interpolator.h"

#include <limits>

namespace Herwig {

using namespace ThePEG;

/**
 * Phenomenological matrix element for QCD-instanton event generation.
 *
 * The class registers selected gluon- and quark-initiated final states built
 * from one instanton zero-mode pair per active flavour and a configurable
 * number of gluons. It can apply either a simple multiplicity model or the
 * tabulated KKS model of arXiv:1911.09726.
 *
 * @see \ref MEInstantonInterfaces "The interfaces"
 */
class MEInstanton : public BlobME {

public:
  MEInstanton();
  MEInstanton(const MEInstanton &) = default;
  ~MEInstanton() override;

  /** Validate options and construct the KKS interpolation tables. */
  void doinit() override;

  /** Synchronize the cached gluon cap at the beginning of a run. */
  void doinitrun() override;

  /** The KKS table is not assigned a fixed perturbative alpha_s order. */
  unsigned int orderInAlphaS() const override {
    return std::numeric_limits<unsigned int>::max();
  }

  unsigned int orderInAlphaEW() const override { return 0; }

  /** Return the dimensionless weight expected by BlobMEBase. */
  double me2() const override;

  /** Register every retained gluon and quark-flavour channel. */
  multimap<tcPDPair, tcPDVector> processes() const override;

  /** Construct the selected large-Nc shower colour flow. */
  list<BlobMEBase::ColourConnection> colourConnections() const override;

  /**
   * Return the maximum base multiplicity before BlobME adds NAdditional.
   * Variable flavour modes must reserve space for five quark pairs.
   */
  size_t nOutgoing() const override;

  /** Return the fixed NQuarkPair setting. */
  size_t nQuarkPair() const { return theNQuarkPair; }

  /**
   * Return the common hard-process scale. In KKS mode this is the selected
   * FactorizationScale option; PureMultiplicity retains shat.
   */
  Energy2 scale() const override;

  void persistentOutput(PersistentOStream & os) const;
  void persistentInput(PersistentIStream & is, int version);

  /** Register the user-facing ThePEG interfaces. */
  static void Init();

protected:
  IBPtr clone() const override;
  IBPtr fullclone() const override;

private:
  MEInstanton & operator=(const MEInstanton &) = delete;

  /** Count active zero-mode flavour pairs in the selected subprocess. */
  size_t currentNQuarkPairs() const;

  /** Count incoming gluons in the selected subprocess. */
  size_t currentNIncomingGluons() const;

  /** Count outgoing gluons in the selected subprocess. */
  size_t currentNFinalGluons() const;

  /** Convert incoming content into the selected outgoing-gluon shift. */
  size_t gluonMultiplicityShift(size_t nIncomingGluons) const;

  /** Return the unshifted multiplicity to which the model weight applies. */
  size_t currentGluonMultiplicity() const;

  /** Evaluate the selected KKS scale for the current partonic energy. */
  Energy2 selectedKKSScale() const;

  /** Build interpolators for the KKS table and derived fermion overlap. */
  void setupInterpolators();

  /** Fixed number of quark pairs used by QuarkPairs=Fixed. */
  size_t theNQuarkPair;

  /** Cached maximum number of additional gluons. */
  size_t theNgluonMax;

  /** PureMultiplicity gluon-distribution option. */
  unsigned int theMultiplicityOption;

  /** Matrix-element model option. */
  unsigned int theModelOption;

  /** Enabled family of incoming partonic subprocesses. */
  unsigned int theProcessOption;

  /** Convention relating sampled and outgoing gluon multiplicities. */
  unsigned int theGluonCountingOption;

  /** Parameters of the PureMultiplicity Gaussian distribution. */
  double theGaussianParamA;
  double theGaussianParamB;

  /** Mean of the PureMultiplicity Poisson distribution. */
  double thePoissonMean;

  /** Shower colour-flow option. */
  unsigned int theColourOption;

  /** Common KKS hard-scale option. */
  unsigned int theScaleOption;

  /** Fixed, legacy variable, or KKS-variable flavour selection. */
  unsigned int theQuarkPairOption;

  /** Bottom mass entering the KKS active-flavour condition. */
  Energy theKKSBottomMass;

  /** Interpolators for the published and derived KKS table columns. */
  Interpolator<double, double>::Ptr theInverseRhoInterpolator;
  Interpolator<double, double>::Ptr theAlphaSInterpolator;
  Interpolator<double, double>::Ptr theMeanGluonsInterpolator;
  Interpolator<double, double>::Ptr theCrossSectionInterpolator;
  Interpolator<double, double>::Ptr theFermionOverlapInterpolator;
};

} // namespace Herwig

#endif // Herwig_MEInstanton_H
