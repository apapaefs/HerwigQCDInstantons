// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEInstanton class.
//

#include "MEInstanton.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include <gsl/gsl_sf_hyperg.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

using namespace Herwig;

namespace {

constexpr double minKKSEnergy = 10.7;
constexpr double maxKKSEnergy = 2895.5;
constexpr double kappa4 = 0.008;
constexpr double kappa5 = 0.01;
constexpr double picobarnToInverseGeV2 = 2.568e-9;
constexpr size_t kKKSNodeCount = 20;

enum ModelOption : unsigned int {
  PureMultiplicityModel = 0,
  KKSModel = 1
};

enum MultiplicityOption : unsigned int {
  PoissonMultiplicity = 0,
  GaussianMultiplicity = 1,
  FlatMultiplicity = 2,
  UserDefinedMultiplicity = 3
};

enum ColourOption : unsigned int {
  SimpleColour = 0,
  RandomColour = 1,
  Random2Colour = 2,
  Random3Colour = 3,
  QCDINSPlanarColour = 4
};

enum ScaleOption : unsigned int {
  InverseRhoScale = 0,
  SHatScale = 1,
  SHatOverMeanGluonsScale = 2
};

enum QuarkPairOption : unsigned int {
  FixedQuarkPairs = 0,
  LegacyVariableQuarkPairs = 1,
  KKSVariableQuarkPairs = 2
};

// Valley action used in the KKS saddle-point equation, arXiv:1911.09726.
double instantonAction(double chi) {
  const double z =
    0.5*(2.0 + chi*chi + chi*std::sqrt(4.0 + chi*chi));
  const double zminus = z - 1.0/z;
  const double zplus = z + 1.0/z;
  return 3.0*((6.0*z*z - 14.0)/(zminus*zminus) - 17.0/3.0 -
              std::log(z)*
                ((z - 5.0/z)*zplus*zplus/(zminus*zminus*zminus) - 1.0));
}

// Fourth-order central derivative. It is evaluated only while building the
// 20-node interpolation table, not for every generated phase-space point.
double instantonActionDerivative(double chi) {
  const double h = 1.e-5*std::max(1.0, std::abs(chi));
  return (instantonAction(chi - 2.0*h) -
          8.0*instantonAction(chi - h) +
          8.0*instantonAction(chi + h) -
          instantonAction(chi + 2.0*h))/(12.0*h);
}

// Solve S'(chi)=target in the physical KKS interval. Bisection is slower than
// Newton iteration but guarantees that the precomputation remains bounded.
double solveKKSChi(double target) {
  double lower = 1.0;
  double upper = 3.0;
  double fLower = instantonActionDerivative(lower) - target;
  const double fUpper = instantonActionDerivative(upper) - target;
  if (!std::isfinite(fLower) || !std::isfinite(fUpper)
      || fLower*fUpper > 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  for (unsigned int iteration = 0; iteration < 100; ++iteration) {
    const double middle = 0.5*(lower + upper);
    const double fMiddle = instantonActionDerivative(middle) - target;
    if (!std::isfinite(fMiddle)) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    if (std::abs(fMiddle) < 1.e-11 || upper - lower < 1.e-11) {
      return middle;
    }
    if (fLower*fMiddle > 0.0) {
      lower = middle;
      fLower = fMiddle;
    } else {
      upper = middle;
    }
  }
  return 0.5*(lower + upper);
}

// Fermion zero-mode overlap in the instanton-anti-instanton background.
double omegaFermion(double chi) {
  const double z =
    0.5*(2.0 + chi*chi + chi*std::sqrt(4.0 + chi*chi));
  const double argument = 1.0 - 1.0/(z*z);
  return (3.0*Constants::pi/8.0)*std::pow(z, -1.5)*
         gsl_sf_hyperg_2F1(1.5, 1.5, 4.0, argument);
}

double poissonProbability(size_t multiplicity, double mean) {
  if (!std::isfinite(mean) || mean < 0.0) return 0.0;
  if (mean == 0.0) return multiplicity == 0 ? 1.0 : 0.0;
  const double logProbability =
    -mean + multiplicity*std::log(mean) -
    std::lgamma(static_cast<double>(multiplicity) + 1.0);
  return std::exp(logProbability);
}

// Normalize the retained 0..maximum channels in log space. This keeps the
// tabulated inclusive KKS cross section independent of the configured cap.
double truncatedPoissonProbability(size_t multiplicity, size_t maximum,
                                   double mean) {
  if (!std::isfinite(mean) || mean < 0.0 || multiplicity > maximum) return 0.0;
  if (mean == 0.0) return multiplicity == 0 ? 1.0 : 0.0;

  const double logMean = std::log(mean);
  double largestLogWeight = -std::numeric_limits<double>::infinity();
  for (size_t n = 0; n <= maximum; ++n) {
    const double logWeight =
      n*logMean - std::lgamma(static_cast<double>(n) + 1.0);
    largestLogWeight = std::max(largestLogWeight, logWeight);
  }

  double normalization = 0.0;
  for (size_t n = 0; n <= maximum; ++n) {
    const double logWeight =
      n*logMean - std::lgamma(static_cast<double>(n) + 1.0);
    normalization += std::exp(logWeight - largestLogWeight);
  }

  const double selectedLogWeight =
    multiplicity*logMean -
    std::lgamma(static_cast<double>(multiplicity) + 1.0);
  if (!std::isfinite(normalization) || normalization <= 0.0) return 0.0;
  return std::exp(selectedLogWeight - largestLogWeight)/normalization;
}

// Convert the KKS four- and five-flavour fermion factors into probabilities.
// Normalizing W4+W5 partitions, rather than rescales, the inclusive table.
double kksFlavourProbability(size_t nQuarkPairs, double omega,
                             double bottomMassRho) {
  if (!std::isfinite(omega) || omega <= 0.0
      || !std::isfinite(bottomMassRho) || bottomMassRho < 0.0) {
    return 0.0;
  }

  double weight4 = kappa4*kappa4*std::pow(omega, 8);
  double weight5 = 0.0;
  if (bottomMassRho <= 1.0) {
    weight4 *= bottomMassRho*bottomMassRho;
    weight5 = kappa5*kappa5*std::pow(omega, 10);
  }

  const double normalization = weight4 + weight5;
  if (!std::isfinite(normalization) || normalization <= 0.0) return 0.0;
  if (nQuarkPairs == 4) return weight4/normalization;
  if (nQuarkPairs == 5) return weight5/normalization;
  return 0.0;
}

bool hasNoSelfConnections(const std::vector<int> & colours,
                          const std::vector<int> & anticolours) {
  if (colours.size() != anticolours.size()) return false;
  for (size_t i = 0; i < colours.size(); ++i) {
    if (colours[i] == anticolours[i]) return false;
  }
  return true;
}

// Fisher-Yates shuffle using ThePEG's event-generator random stream.
void randomShuffle(std::vector<int> & values) {
  for (size_t i = values.size(); i > 1; --i) {
    const size_t selected = static_cast<size_t>(UseRandom::irnd(i));
    std::swap(values[i - 1], values[selected]);
  }
}

// Draw a bijection from colour sources to anticolour sinks. Equal external-leg
// indices identify the two ends of one gluon and are therefore forbidden.
std::vector<int> randomColourMap(const std::vector<int> & colours,
                                 const std::vector<int> & anticolours) {
  if (colours.empty() || colours.size() != anticolours.size()) {
    throw Exception() << "MEInstanton: invalid colour-map inputs."
                      << Exception::runerror;
  }

  for (unsigned int attempt = 0; attempt < 1000; ++attempt) {
    std::vector<int> result = anticolours;
    randomShuffle(result);
    if (hasNoSelfConnections(colours, result)) return result;
  }

  // This should be unreachable in practice, but a cyclic search gives a
  // deterministic final attempt before reporting an invalid endpoint set.
  for (size_t shift = 0; shift < anticolours.size(); ++shift) {
    std::vector<int> result(anticolours.size());
    for (size_t i = 0; i < anticolours.size(); ++i) {
      result[i] = anticolours[(i + shift) % anticolours.size()];
    }
    if (hasNoSelfConnections(colours, result)) return result;
  }

  throw Exception() << "MEInstanton: no one-to-one colour map without "
                    << "a forbidden self-connection exists."
                    << Exception::runerror;
}

// BlobMEBase describes incoming legs in the crossed, all-outgoing convention.
// Consequently an incoming colour source is recorded as an anticolour end,
// while an incoming anticolour sink is recorded as a colour end.
BlobMEBase::ColourConnection connectColourEndpoints(int colourSource,
                                                     int anticolourSink) {
  BlobMEBase::ColourConnection line;
  if (colourSource < 2) {
    line.addAntiColour(colourSource);
  } else {
    line.addColour(colourSource);
  }

  if (anticolourSink < 2) {
    line.addColour(anticolourSink);
  } else {
    line.addAntiColour(anticolourSink);
  }
  return line;
}

} // namespace

MEInstanton::MEInstanton()
  : theNQuarkPair(4), theNgluonMax(1),
    theMultiplicityOption(PoissonMultiplicity),
    theModelOption(PureMultiplicityModel), theGaussianParamA(5.0),
    theGaussianParamB(200.0), thePoissonMean(3.0),
    theColourOption(SimpleColour), theScaleOption(InverseRhoScale),
    theQuarkPairOption(FixedQuarkPairs), theKKSBottomMass(4.18*GeV) {}

MEInstanton::~MEInstanton() = default;

IBPtr MEInstanton::clone() const {
  return new_ptr(*this);
}

IBPtr MEInstanton::fullclone() const {
  return new_ptr(*this);
}

void MEInstanton::setupInterpolators() {
  // Table 1 of arXiv:1911.09726. Energies and inverse radii are in GeV;
  // partonic cross sections are in pb.
  static const std::array<double, kKKSNodeCount> energies = {{
    10.7, 11.4, 13.4, 15.7, 22.9, 29.7, 40.8, 56.1, 61.8, 89.6,
    118.0, 174.4, 246.9, 349.9, 496.3, 704.8, 1001.8, 1425.6,
    2030.6, 2895.5
  }};
  static const std::array<double, kKKSNodeCount> inverseRho = {{
    0.99, 1.04, 1.16, 1.31, 1.76, 2.12, 2.72, 3.50, 3.64, 4.98,
    6.21, 8.72, 11.76, 15.90, 21.58, 29.37, 40.07, 54.83, 75.21,
    103.4
  }};
  static const std::array<double, kKKSNodeCount> alphaSAtInverseRho = {{
    0.416, 0.405, 0.382, 0.360, 0.315, 0.293, 0.267, 0.245, 0.223,
    0.206, 0.195, 0.180, 0.169, 0.159, 0.150, 0.142, 0.135, 0.128,
    0.122, 0.117
  }};
  static const std::array<double, kKKSNodeCount> meanGluons = {{
    4.59, 4.68, 4.90, 5.13, 5.44, 6.02, 6.47, 6.92, 7.28, 7.67,
    8.25, 8.60, 9.04, 9.49, 9.93, 10.37, 10.81, 11.26, 11.70,
    12.14
  }};
  static const std::array<double, kKKSNodeCount> crossSections = {{
    4.922e9, 3.652e9, 1.671e9, 728.9e6, 85.94e6, 17.25e6,
    2.121e6, 229.0e3, 72.97e3, 2.733e3, 235.4, 6.720, 0.284,
    0.012, 5.112e-4, 21.65e-6, 0.9017e-6, 36.45e-9, 1.419e-9,
    52.07e-12
  }};

  // VariableKKS needs the fermion overlap, which is not tabulated directly.
  // Derive it once at each node from the published alpha_s and rho columns.
  std::array<double, kKKSNodeCount> fermionOverlaps;
  for (size_t i = 0; i < energies.size(); ++i) {
    const double u = energies[i]/inverseRho[i];
    const double rhoTilde =
      alphaSAtInverseRho[i]*u/(4.0*Constants::pi);
    const double chi = solveKKSChi(rhoTilde);
    if (!std::isfinite(chi)) {
      throw InitException()
        << "MEInstanton: failed to solve the KKS saddle point at table node "
        << i << "." << Exception::abortnow;
    }

    fermionOverlaps[i] = omegaFermion(chi);
    if (!std::isfinite(fermionOverlaps[i]) || fermionOverlaps[i] <= 0.0) {
      throw InitException()
        << "MEInstanton: invalid KKS fermion overlap at table node " << i
        << "." << Exception::abortnow;
    }
  }

  theInverseRhoInterpolator = make_InterpolatorPtr(inverseRho, energies, 1);
  theAlphaSInterpolator =
    make_InterpolatorPtr(alphaSAtInverseRho, energies, 1);
  theMeanGluonsInterpolator = make_InterpolatorPtr(meanGluons, energies, 1);
  theCrossSectionInterpolator =
    make_InterpolatorPtr(crossSections, energies, 1);
  theFermionOverlapInterpolator =
    make_InterpolatorPtr(fermionOverlaps, energies, 1);
}

Energy2 MEInstanton::scale() const {
  if (theModelOption == KKSModel) return selectedKKSScale();
  return sHat();
}

Energy2 MEInstanton::selectedKKSScale() const {
  const Energy2 shat = sHat();
  const double shatValue = shat/GeV/GeV;
  if (!std::isfinite(shatValue) || shatValue <= 0.0) return shat;

  const double partonicEnergy = std::sqrt(shatValue);
  // me2() rejects points outside the table. Returning shat here avoids asking
  // an interpolator to extrapolate while ThePEG is preparing that point.
  if (partonicEnergy < minKKSEnergy || partonicEnergy > maxKKSEnergy) {
    return shat;
  }

  if (theScaleOption == SHatScale) return shat;

  if (theScaleOption == SHatOverMeanGluonsScale) {
    const double meanGluons = (*theMeanGluonsInterpolator)(partonicEnergy);
    if (!std::isfinite(meanGluons) || meanGluons <= 0.0) return shat;
    return shat/meanGluons;
  }

  const double inverseRho = (*theInverseRhoInterpolator)(partonicEnergy);
  if (!std::isfinite(inverseRho) || inverseRho <= 0.0) return shat;
  return sqr(inverseRho)*GeV*GeV;
}

double MEInstanton::me2() const {
  const size_t nQuarkPairs = currentNQuarkPairs();
  const size_t baseMultiplicity = 2 + 2*nQuarkPairs;
  if (meMomenta().size() < baseMultiplicity) return 0.0;

  const size_t nGluons = meMomenta().size() - baseMultiplicity;
  if (nGluons > theNgluonMax) return 0.0;

  const double shat = sHat()/GeV/GeV;
  if (!std::isfinite(shat) || shat <= 0.0) return 0.0;
  const double partonicEnergy = std::sqrt(shat);
  double weight = 1.0;

  if (theModelOption == PureMultiplicityModel) {
    if (theMultiplicityOption == PoissonMultiplicity) {
      weight *= poissonProbability(nGluons, thePoissonMean);
    } else if (theMultiplicityOption == GaussianMultiplicity) {
      if (!std::isfinite(theGaussianParamB) || theGaussianParamB <= 0.0) {
        return 0.0;
      }
      weight *=
        std::exp(-sqr(static_cast<double>(nGluons) - theGaussianParamA)/
                 theGaussianParamB)/
        std::sqrt(Constants::pi*theGaussianParamB);
    }
    return std::isfinite(weight) && weight >= 0.0 ? weight : 0.0;
  }

  if (theModelOption != KKSModel || partonicEnergy < minKKSEnergy
      || partonicEnergy > maxKKSEnergy) {
    return 0.0;
  }

  // ThePEG evaluates the incoming PDFs at scale(). No manual PDF ratio is
  // needed here; this check only protects the remaining calculation.
  const double scaleValue = selectedKKSScale()/GeV/GeV;
  if (!std::isfinite(scaleValue) || scaleValue <= 0.0) return 0.0;

  const double meanGluons = (*theMeanGluonsInterpolator)(partonicEnergy);
  const double partonicCrossSection =
    (*theCrossSectionInterpolator)(partonicEnergy);
  if (!std::isfinite(meanGluons) || meanGluons < 0.0
      || !std::isfinite(partonicCrossSection) || partonicCrossSection < 0.0) {
    return 0.0;
  }

  weight *=
    truncatedPoissonProbability(nGluons, theNgluonMax, meanGluons);
  weight *= partonicCrossSection*picobarnToInverseGeV2;

  const double phaseSpaceJacobian = jacobian();
  if (!std::isfinite(phaseSpaceJacobian) || phaseSpaceJacobian <= 0.0) {
    return 0.0;
  }

  // BlobMEBase later multiplies me2() by the generated phase-space Jacobian
  // and divides by the 2*shat flux. Undo those factors so the integration
  // reproduces the tabulated partonic cross section.
  weight *= 2.0*shat/phaseSpaceJacobian;

  if (theQuarkPairOption == LegacyVariableQuarkPairs) {
    weight *= 0.5;
  } else if (theQuarkPairOption == KKSVariableQuarkPairs) {
    const double inverseRho =
      (*theInverseRhoInterpolator)(partonicEnergy);
    const double omega =
      (*theFermionOverlapInterpolator)(partonicEnergy);
    if (!std::isfinite(inverseRho) || inverseRho <= 0.0) return 0.0;

    const double bottomMassRho = (theKKSBottomMass/GeV)/inverseRho;
    const double flavourProbability = kksFlavourProbability(
      nQuarkPairs, omega, bottomMassRho);
    if (!std::isfinite(flavourProbability) || flavourProbability < 0.0) {
      return 0.0;
    }
    weight *= flavourProbability;
  }

  return std::isfinite(weight) && weight >= 0.0 ? weight : 0.0;
}

void MEInstanton::doinit() {
  if (theNQuarkPair < 1 || theNQuarkPair > 5) {
    throw InitException() << "MEInstanton: NQuarkPair must be between 1 and 5."
                          << Exception::abortnow;
  }
  if (theQuarkPairOption == KKSVariableQuarkPairs
      && theModelOption != KKSModel) {
    throw InitException()
      << "MEInstanton: QuarkPairs VariableKKS requires MEModeling KKS."
      << Exception::abortnow;
  }
  if (!std::isfinite(thePoissonMean) || thePoissonMean < 0.0) {
    throw InitException() << "MEInstanton: PoissonMean must be non-negative."
                          << Exception::abortnow;
  }
  if (!std::isfinite(theGaussianParamA) || theGaussianParamA < 0.0) {
    throw InitException()
      << "MEInstanton: GaussianParamA must be non-negative."
      << Exception::abortnow;
  }
  if (!std::isfinite(theGaussianParamB) || theGaussianParamB <= 0.0) {
    throw InitException() << "MEInstanton: GaussianParamB must be positive."
                          << Exception::abortnow;
  }
  if (!std::isfinite(theKKSBottomMass/GeV) || theKKSBottomMass < ZERO) {
    throw InitException() << "MEInstanton: KKSBottomMass must be non-negative."
                          << Exception::abortnow;
  }

  theNgluonMax = nAdditional();
  setupInterpolators();
  BlobME::doinit();
}

void MEInstanton::doinitrun() {
  theNgluonMax = nAdditional();
  BlobME::doinitrun();
}

multimap<tcPDPair, tcPDVector> MEInstanton::processes() const {
  multimap<tcPDPair, tcPDVector> processMap;
  tcPDPtr gluon = getParticleData(ParticleID::g);
  std::vector<tcPDPtr> quarks;
  std::vector<tcPDPtr> antiquarks;
  quarks.reserve(5);
  antiquarks.reserve(5);
  for (int id = 1; id <= 5; ++id) {
    tcPDPtr quark = getParticleData(id);
    quarks.push_back(quark);
    antiquarks.push_back(quark->CC());
  }
  const tcPDPair incoming = std::make_pair(gluon, gluon);

  // Each multiplicity is a separate BlobME subprocess. Variable flavour modes
  // register both bases; me2() supplies the relative channel probability.
  const auto addProcesses = [&](size_t nQuarkPairs) {
    tcPDVector outgoing;
    outgoing.reserve(2*nQuarkPairs + theNgluonMax);
    for (size_t i = 0; i < nQuarkPairs; ++i) {
      outgoing.push_back(quarks[i]);
      outgoing.push_back(antiquarks[i]);
    }

    processMap.insert(std::make_pair(incoming, outgoing));
    for (size_t n = 0; n < theNgluonMax; ++n) {
      outgoing.push_back(gluon);
      processMap.insert(std::make_pair(incoming, outgoing));
    }
  };

  if (theQuarkPairOption == FixedQuarkPairs) {
    addProcesses(nQuarkPair());
  } else {
    addProcesses(4);
    addProcesses(5);
  }
  return processMap;
}

size_t MEInstanton::currentNQuarkPairs() const {
  if (theQuarkPairOption == FixedQuarkPairs) return theNQuarkPair;

  // In a variable mode, the selected subprocess is the authoritative source
  // for whether this point contains four or five pairs.
  size_t countedQuarks = 0;
  for (size_t i = 2; i < mePartonData().size(); ++i) {
    if (!mePartonData()[i]) continue;
    const long id = std::abs(mePartonData()[i]->id());
    if (id >= 1 && id <= 5) ++countedQuarks;
  }
  return countedQuarks/2;
}

list<BlobMEBase::ColourConnection> MEInstanton::colourConnections() const {
  list<BlobMEBase::ColourConnection> result;
  const size_t nQuarkPairs = currentNQuarkPairs();
  const size_t baseMultiplicity = 2 + 2*nQuarkPairs;
  if (nQuarkPairs == 0 || meMomenta().size() < baseMultiplicity
      || meMomenta().size() >
           static_cast<size_t>(std::numeric_limits<int>::max())) {
    throw Exception()
      << "MEInstanton: inconsistent multiplicity in colourConnections()."
      << Exception::runerror;
  }

  const size_t nGluons = meMomenta().size() - baseMultiplicity;
  const int firstGluon = static_cast<int>(2 + 2*nQuarkPairs);

  if (theColourOption == SimpleColour) {
    // Keep the incoming gluons in a two-line singlet disconnected from the
    // final state.
    BlobMEBase::ColourConnection incoming1;
    BlobMEBase::ColourConnection incoming2;
    incoming1.addColour(0);
    incoming1.addAntiColour(1);
    incoming2.addColour(1);
    incoming2.addAntiColour(0);
    result.push_back(incoming1);
    result.push_back(incoming2);

    // Directly connect all but the last quark-antiquark pair. The last pair
    // absorbs one gluon when the gluon multiplicity is odd.
    for (size_t q = 0; q + 1 < nQuarkPairs; ++q) {
      BlobMEBase::ColourConnection line;
      line.addColour(static_cast<int>(2 + 2*q));
      line.addAntiColour(static_cast<int>(3 + 2*q));
      result.push_back(line);
    }

    const int lastQuark = static_cast<int>(2 + 2*(nQuarkPairs - 1));
    const int lastAntiquark = lastQuark + 1;
    size_t firstPairedGluon = 0;
    if (nGluons % 2 == 0) {
      BlobMEBase::ColourConnection lastPair;
      lastPair.addColour(lastQuark);
      lastPair.addAntiColour(lastAntiquark);
      result.push_back(lastPair);
    } else {
      BlobMEBase::ColourConnection quarkToGluon;
      BlobMEBase::ColourConnection gluonToAntiquark;
      quarkToGluon.addColour(lastQuark);
      quarkToGluon.addAntiColour(firstGluon);
      gluonToAntiquark.addColour(firstGluon);
      gluonToAntiquark.addAntiColour(lastAntiquark);
      result.push_back(quarkToGluon);
      result.push_back(gluonToAntiquark);
      firstPairedGluon = 1;
    }

    // Every remaining pair of gluons forms a closed two-line singlet.
    for (size_t g = firstPairedGluon; g + 1 < nGluons; g += 2) {
      const int first = firstGluon + static_cast<int>(g);
      const int second = first + 1;
      BlobMEBase::ColourConnection line1;
      BlobMEBase::ColourConnection line2;
      line1.addColour(first);
      line1.addAntiColour(second);
      line2.addColour(second);
      line2.addAntiColour(first);
      result.push_back(line1);
      result.push_back(line2);
    }
    return result;
  }

  if (theColourOption == QCDINSPlanarColour) {
    // QCDINS 2.0 models its leading-colour final state as N_f open
    // q-g-...-g-qbar strings. For this gg-initiated adaptation, both crossed
    // incoming gluons are inserted alongside every outgoing gluon. Thus each
    // gluon belongs to a component anchored by a quark and an antiquark.
    std::vector<int> antiquarks;
    antiquarks.reserve(nQuarkPairs);
    for (size_t q = 0; q < nQuarkPairs; ++q) {
      antiquarks.push_back(static_cast<int>(3 + 2*q));
    }
    randomShuffle(antiquarks);

    std::vector<int> gluons;
    gluons.reserve(nGluons + 2);
    gluons.push_back(0);
    gluons.push_back(1);
    for (size_t g = 0; g < nGluons; ++g) {
      gluons.push_back(firstGluon + static_cast<int>(g));
    }
    randomShuffle(gluons);

    std::vector<std::vector<int>> strings(nQuarkPairs);
    for (const int gluon : gluons) {
      const size_t selected =
        static_cast<size_t>(UseRandom::irnd(nQuarkPairs));
      strings[selected].push_back(gluon);
    }

    for (size_t q = 0; q < nQuarkPairs; ++q) {
      int colourSource = static_cast<int>(2 + 2*q);
      for (const int gluon : strings[q]) {
        result.push_back(connectColourEndpoints(colourSource, gluon));
        colourSource = gluon;
      }
      result.push_back(connectColourEndpoints(colourSource, antiquarks[q]));
    }
    return result;
  }

  if (theColourOption < RandomColour || theColourOption > Random3Colour) {
    throw Exception() << "MEInstanton: unknown colour-connection option."
                      << Exception::runerror;
  }

  std::vector<int> colourSources;
  std::vector<int> anticolourSinks;
  colourSources.reserve(nQuarkPairs + nGluons + 2);
  anticolourSinks.reserve(nQuarkPairs + nGluons + 2);

  // Random leaves both incoming lines fixed. Random2 exposes one crossed
  // source/sink pair to the random map. Random3 exposes both incoming gluons,
  // matching Sherpa's endpoint pool. Herwig draws a complete valid
  // permutation, whereas Sherpa selects pairs sequentially and repairs a
  // forbidden self-match if it is the only final pairing left.
  if (theColourOption == Random2Colour) {
    colourSources.push_back(0);
    anticolourSinks.push_back(1);
  } else if (theColourOption == Random3Colour) {
    colourSources.push_back(0);
    colourSources.push_back(1);
    anticolourSinks.push_back(0);
    anticolourSinks.push_back(1);
  }

  for (size_t q = 0; q < nQuarkPairs; ++q) {
    colourSources.push_back(static_cast<int>(2 + 2*q));
    anticolourSinks.push_back(static_cast<int>(3 + 2*q));
  }
  for (size_t g = 0; g < nGluons; ++g) {
    const int index = firstGluon + static_cast<int>(g);
    colourSources.push_back(index);
    anticolourSinks.push_back(index);
  }

  const std::vector<int> colourMap =
    randomColourMap(colourSources, anticolourSinks);
  for (size_t i = 0; i < colourSources.size(); ++i) {
    result.push_back(connectColourEndpoints(colourSources[i], colourMap[i]));
  }

  // Add the incoming lines that were deliberately excluded from the random
  // map. Random3 has none: all four incoming colour ends already participate.
  if (theColourOption == RandomColour) {
    BlobMEBase::ColourConnection incoming1;
    BlobMEBase::ColourConnection incoming2;
    incoming1.addColour(0);
    incoming1.addAntiColour(1);
    incoming2.addColour(1);
    incoming2.addAntiColour(0);
    result.push_back(incoming1);
    result.push_back(incoming2);
  } else if (theColourOption == Random2Colour) {
    BlobMEBase::ColourConnection incoming;
    incoming.addColour(0);
    incoming.addAntiColour(1);
    result.push_back(incoming);
  }
  return result;
}

size_t MEInstanton::nOutgoing() const {
  return theQuarkPairOption == FixedQuarkPairs ? 2*nQuarkPair() : 10;
}

void MEInstanton::persistentOutput(PersistentOStream & os) const {
  // Keep this order synchronized with persistentInput for run-file
  // compatibility. Renaming the C++ members does not change the stream layout.
  os << theNQuarkPair << theNgluonMax << theMultiplicityOption
     << theModelOption << theGaussianParamA << theGaussianParamB
     << thePoissonMean << theColourOption << theInverseRhoInterpolator
     << theAlphaSInterpolator << theMeanGluonsInterpolator
     << theCrossSectionInterpolator << theScaleOption << theQuarkPairOption
     << ounit(theKKSBottomMass, GeV) << theFermionOverlapInterpolator;
}

void MEInstanton::persistentInput(PersistentIStream & is, int) {
  is >> theNQuarkPair >> theNgluonMax >> theMultiplicityOption
     >> theModelOption >> theGaussianParamA >> theGaussianParamB
     >> thePoissonMean >> theColourOption >> theInverseRhoInterpolator
     >> theAlphaSInterpolator >> theMeanGluonsInterpolator
     >> theCrossSectionInterpolator >> theScaleOption >> theQuarkPairOption
     >> iunit(theKKSBottomMass, GeV) >> theFermionOverlapInterpolator;
}

DescribeClass<MEInstanton, Herwig::BlobME>
  describeHerwigMEInstanton("Herwig::MEInstanton", "Instantons.so");

void MEInstanton::Init() {
  static ClassDocumentation<MEInstanton> documentation(
    "Phenomenological QCD-instanton matrix element with configurable "
    "multiplicity, flavour, scale and colour models."
  );

  static Parameter<MEInstanton, size_t> interfaceNQuarkPair(
    "NQuarkPair",
    "The fixed number of quark pairs.",
    &MEInstanton::theNQuarkPair,
    4, 1, 5,
    false, false, Interface::limited
  );

  static Switch<MEInstanton, unsigned int> interfaceColourConnections(
    "ColourConnections",
    "How to connect the large-Nc shower colour lines.",
    &MEInstanton::theColourOption,
    SimpleColour,
    false, false
  );
  static SwitchOption interfaceColourConnectionsSimple(
    interfaceColourConnections,
    "Simple",
    "Use a deterministic singlet-oriented assignment.",
    SimpleColour
  );
  static SwitchOption interfaceColourConnectionsRandom(
    interfaceColourConnections,
    "Random",
    "Randomize final-state lines and keep the incoming singlet fixed.",
    RandomColour
  );
  static SwitchOption interfaceColourConnectionsRandom2(
    interfaceColourConnections,
    "Random2",
    "Randomize the final state with one crossed incoming colour line.",
    Random2Colour
  );
  static SwitchOption interfaceColourConnectionsRandom3(
    interfaceColourConnections,
    "Random3",
    "Randomize both crossed incoming colour lines with the final state.",
    Random3Colour
  );
  static SwitchOption interfaceColourConnectionsQCDINSPlanar(
    interfaceColourConnections,
    "QCDINSPlanar",
    "Place every incoming and outgoing gluon on a randomly constructed "
    "planar quark-antiquark string inspired by QCDINS 2.0.",
    QCDINSPlanarColour
  );

  static Switch<MEInstanton, unsigned int> interfaceFactorizationScale(
    "FactorizationScale",
    "The common hard-process scale used in KKS mode.",
    &MEInstanton::theScaleOption,
    InverseRhoScale,
    false, false
  );
  static SwitchOption interfaceFactorizationScaleInvRho(
    interfaceFactorizationScale,
    "InvRho",
    "Use inverse-rho squared as the common scale.",
    InverseRhoScale
  );
  static SwitchOption interfaceFactorizationScaleSHat(
    interfaceFactorizationScale,
    "sHat",
    "Use the partonic centre-of-mass energy squared.",
    SHatScale
  );
  static SwitchOption interfaceFactorizationScaleSHatOverN(
    interfaceFactorizationScale,
    "sHatOverN",
    "Use shat divided by the interpolated mean gluon multiplicity.",
    SHatOverMeanGluonsScale
  );

  static Switch<MEInstanton, unsigned int> interfaceQuarkPairs(
    "QuarkPairs",
    "How to select the number of quark pairs.",
    &MEInstanton::theQuarkPairOption,
    FixedQuarkPairs,
    false, false
  );
  static SwitchOption interfaceNQuarkPairsFixed(
    interfaceQuarkPairs,
    "Fixed",
    "Use NQuarkPair.",
    FixedQuarkPairs
  );
  static SwitchOption interfaceNQuarkPairsVariable(
    interfaceQuarkPairs,
    "Variable",
    "Use the legacy equal mixture of four and five pairs.",
    LegacyVariableQuarkPairs
  );
  static SwitchOption interfaceNQuarkPairsVariableKKS(
    interfaceQuarkPairs,
    "VariableKKS",
    "Use KKS scale-, mass- and overlap-dependent four/five-pair weights.",
    KKSVariableQuarkPairs
  );

  static Parameter<MEInstanton, Energy> interfaceKKSBottomMass(
    "KKSBottomMass",
    "Bottom-quark mass entering the KKS active-flavour condition.",
    &MEInstanton::theKKSBottomMass,
    GeV, 4.18*GeV, ZERO, ZERO,
    false, false, Interface::lowerlim
  );

  static Switch<MEInstanton, unsigned int>
    interfaceMultiplicityParametrisation(
      "MultiplicityParametrisation",
      "How to weight gluon multiplicities in PureMultiplicity mode.",
      &MEInstanton::theMultiplicityOption,
      PoissonMultiplicity,
      false, false
    );
  static SwitchOption interfaceMultiplicityParametrisationPoisson(
    interfaceMultiplicityParametrisation,
    "Poisson",
    "Use a Poisson distribution with mean PoissonMean.",
    PoissonMultiplicity
  );
  static SwitchOption interfaceMultiplicityParametrisationGaussian(
    interfaceMultiplicityParametrisation,
    "Gaussian",
    "Use a Gaussian controlled by GaussianParamA and GaussianParamB.",
    GaussianMultiplicity
  );
  static SwitchOption interfaceMultiplicityParametrisationFlat(
    interfaceMultiplicityParametrisation,
    "Flat",
    "Apply no multiplicity-dependent factor.",
    FlatMultiplicity
  );
  static SwitchOption interfaceMultiplicityParametrisationUserDefined(
    interfaceMultiplicityParametrisation,
    "UserDefined",
    "Reserved for a user-defined factor; currently equivalent to Flat.",
    UserDefinedMultiplicity
  );

  static Switch<MEInstanton, unsigned int> interfaceMEModeling(
    "MEModeling",
    "How to model the matrix element.",
    &MEInstanton::theModelOption,
    PureMultiplicityModel,
    false, false
  );
  static SwitchOption interfaceMEModelingPureMultiplicity(
    interfaceMEModeling,
    "PureMultiplicity",
    "Use a flat matrix element with a configurable multiplicity factor.",
    PureMultiplicityModel
  );
  static SwitchOption interfaceMEModelingKKS(
    interfaceMEModeling,
    "KKS",
    "Use the tabulated Khoze-Krauss-Schott model (arXiv:1911.09726).",
    KKSModel
  );

  static Parameter<MEInstanton, double> interfaceGaussianParamA(
    "GaussianParamA",
    "Centre of the PureMultiplicity Gaussian.",
    &MEInstanton::theGaussianParamA,
    5.0, 0.0, 0.0,
    false, false, Interface::lowerlim
  );

  static Parameter<MEInstanton, double> interfaceGaussianParamB(
    "GaussianParamB",
    "Positive denominator in the PureMultiplicity Gaussian exponent.",
    &MEInstanton::theGaussianParamB,
    200.0, 0.0, 0.0,
    false, false, Interface::lowerlim
  );

  static Parameter<MEInstanton, double> interfacePoissonMean(
    "PoissonMean",
    "Non-negative mean of the PureMultiplicity Poisson distribution.",
    &MEInstanton::thePoissonMean,
    3.0, 0.0, 0.0,
    false, false, Interface::lowerlim
  );
}
