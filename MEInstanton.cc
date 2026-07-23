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

enum ProcessOption : unsigned int {
  GluonGluonProcesses = 0,
  QuarkGluonProcesses = 1,
  QuarkQuarkProcesses = 2,
  QuarkAntiquarkProcesses = 3,
  AllProcesses = 4
};

enum GluonCountingOption : unsigned int {
  FinalStateGluonCounting = 0,
  FixedTotalGluonCounting = 1
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

// Find a deterministic perfect matching in the complete endpoint graph with
// the diagonal removed. The diagonal corresponds to joining a gluon directly
// to itself; quark legs occur in only one of the two endpoint vectors.
bool augmentColourMap(size_t sourcePosition,
                      const std::vector<int> & colours,
                      const std::vector<int> & anticolours,
                      std::vector<int> & sourceForSink,
                      std::vector<bool> & visitedSinks) {
  for (size_t sinkPosition = 0; sinkPosition < anticolours.size();
       ++sinkPosition) {
    if (visitedSinks[sinkPosition]
        || colours[sourcePosition] == anticolours[sinkPosition]) {
      continue;
    }
    visitedSinks[sinkPosition] = true;
    if (sourceForSink[sinkPosition] < 0
        || augmentColourMap(
             static_cast<size_t>(sourceForSink[sinkPosition]), colours,
             anticolours, sourceForSink, visitedSinks)) {
      sourceForSink[sinkPosition] = static_cast<int>(sourcePosition);
      return true;
    }
  }
  return false;
}

std::vector<int> deterministicColourMap(
    const std::vector<int> & colours,
    const std::vector<int> & anticolours) {
  if (colours.size() != anticolours.size()) {
    throw Exception() << "MEInstanton: unbalanced colour endpoints."
                      << Exception::runerror;
  }
  if (colours.empty()) return {};

  std::vector<int> sourceForSink(anticolours.size(), -1);
  for (size_t sourcePosition = 0; sourcePosition < colours.size();
       ++sourcePosition) {
    std::vector<bool> visitedSinks(anticolours.size(), false);
    if (!augmentColourMap(sourcePosition, colours, anticolours,
                          sourceForSink, visitedSinks)) {
      throw Exception() << "MEInstanton: no one-to-one colour map without "
                        << "a forbidden self-connection exists."
                        << Exception::runerror;
    }
  }

  std::vector<int> result(colours.size(), -1);
  for (size_t sinkPosition = 0; sinkPosition < anticolours.size();
       ++sinkPosition) {
    result[static_cast<size_t>(sourceForSink[sinkPosition])] =
      anticolours[sinkPosition];
  }
  return result;
}

// Draw a bijection from colour sources to anticolour sinks. Equal external-leg
// indices identify the two ends of one gluon and are therefore forbidden.
std::vector<int> randomColourMap(const std::vector<int> & colours,
                                 const std::vector<int> & anticolours) {
  if (colours.size() != anticolours.size()) {
    throw Exception() << "MEInstanton: invalid colour-map inputs."
                      << Exception::runerror;
  }
  if (colours.empty()) return {};

  for (unsigned int attempt = 0; attempt < 1000; ++attempt) {
    std::vector<int> result = anticolours;
    randomShuffle(result);
    if (hasNoSelfConnections(colours, result)) return result;
  }

  return deterministicColourMap(colours, anticolours);
}

// BlobMEBase describes incoming legs in the crossed, all-outgoing convention.
// Consequently an incoming colour source is recorded as an anticolour end,
// while an incoming anticolour sink is recorded as a colour end.
BlobMEBase::ColourConnection connectColourEndpoints(int colourSource,
                                                     int anticolourSink) {
  if (colourSource == anticolourSink) {
    throw Exception() << "MEInstanton: attempted to connect a gluon to itself."
                      << Exception::runerror;
  }
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

template <typename ParticlePointer>
bool isInstantonQuark(const ParticlePointer & particle) {
  if (!particle) return false;
  const long id = std::abs(particle->id());
  return id >= 1 && id <= 5;
}

struct ColourEndpointSets {
  std::vector<int> sources;
  std::vector<int> sinks;
  std::vector<int> incomingSources;
  std::vector<int> incomingSinks;
  std::vector<int> outgoingSources;
  std::vector<int> outgoingSinks;
  std::vector<int> fermionSources;
  std::vector<int> fermionSinks;
  std::vector<int> gluons;
};

// Classify colour ends in the all-outgoing convention used by BlobMEBase.
// This is the same crossing rule used by Sherpa's instanton MakeColours(): an
// incoming quark is an anticolour sink and an incoming antiquark is a colour
// source, while a gluon contributes one end of each type.
ColourEndpointSets collectColourEndpoints(const cPDVector & partons) {
  if (partons.size() < 2
      || partons.size() >
           static_cast<size_t>(std::numeric_limits<int>::max())) {
    throw Exception() << "MEInstanton: invalid external-particle list."
                      << Exception::runerror;
  }

  ColourEndpointSets endpoints;
  for (size_t position = 0; position < partons.size(); ++position) {
    if (!partons[position]) {
      throw Exception() << "MEInstanton: null external ParticleData pointer."
                        << Exception::runerror;
    }

    const int index = static_cast<int>(position);
    const bool incoming = position < 2;
    const long id = partons[position]->id();
    if (id == ParticleID::g) {
      endpoints.sources.push_back(index);
      endpoints.sinks.push_back(index);
      endpoints.gluons.push_back(index);
      if (incoming) {
        endpoints.incomingSources.push_back(index);
        endpoints.incomingSinks.push_back(index);
      } else {
        endpoints.outgoingSources.push_back(index);
        endpoints.outgoingSinks.push_back(index);
      }
      continue;
    }

    if (std::abs(id) < 1 || std::abs(id) > 5) {
      throw Exception() << "MEInstanton: unsupported coloured external leg "
                        << id << "." << Exception::runerror;
    }

    const bool isSource = incoming ? id < 0 : id > 0;
    if (isSource) {
      endpoints.sources.push_back(index);
      endpoints.fermionSources.push_back(index);
      (incoming ? endpoints.incomingSources : endpoints.outgoingSources)
        .push_back(index);
    } else {
      endpoints.sinks.push_back(index);
      endpoints.fermionSinks.push_back(index);
      (incoming ? endpoints.incomingSinks : endpoints.outgoingSinks)
        .push_back(index);
    }
  }

  if (endpoints.sources.size() != endpoints.sinks.size()
      || endpoints.fermionSources.size() != endpoints.fermionSinks.size()) {
    throw Exception() << "MEInstanton: external state does not conserve colour."
                      << Exception::runerror;
  }
  return endpoints;
}

int sourceMember(int index) {
  return index < 2 ? -index - 1 : index + 1;
}

int sinkMember(int index) {
  return index < 2 ? index + 1 : -index - 1;
}

void validateColourConnections(
    const list<BlobMEBase::ColourConnection> & connections,
    const ColourEndpointSets & endpoints) {
  std::vector<int> expected;
  expected.reserve(endpoints.sources.size() + endpoints.sinks.size());
  for (const int source : endpoints.sources) {
    expected.push_back(sourceMember(source));
  }
  for (const int sink : endpoints.sinks) {
    expected.push_back(sinkMember(sink));
  }

  std::vector<int> actual;
  actual.reserve(expected.size());
  for (const auto & connection : connections) {
    if (connection.members.size() != 2) {
      throw Exception() << "MEInstanton: colour flow contains a non-pair line."
                        << Exception::runerror;
    }
    actual.insert(actual.end(), connection.members.begin(),
                  connection.members.end());
  }

  std::sort(expected.begin(), expected.end());
  std::sort(actual.begin(), actual.end());
  if (actual != expected) {
    throw Exception() << "MEInstanton: colour flow does not use every external "
                      << "colour end exactly once." << Exception::runerror;
  }
}

bool eraseEndpoint(std::vector<int> & endpoints, int endpoint) {
  const auto found = std::find(endpoints.begin(), endpoints.end(), endpoint);
  if (found == endpoints.end()) return false;
  endpoints.erase(found);
  return true;
}

int takeOutgoingEndpoint(std::vector<int> & endpoints) {
  const auto found = std::find_if(
    endpoints.begin(), endpoints.end(), [](int endpoint) {
      return endpoint >= 2;
    });
  if (found == endpoints.end()) {
    throw Exception() << "MEInstanton: cannot balance fixed incoming colour "
                      << "ends with the outgoing state."
                      << Exception::runerror;
  }
  const int endpoint = *found;
  endpoints.erase(found);
  return endpoint;
}

// Cross incoming quarks through one zero-mode pair per active flavour.
// Reusing the same zero mode twice (for example u u in the initial state) is
// forbidden; u ubar is valid because it crosses the two distinct legs.
bool crossedFermionFinalState(
    const tcPDPair & incoming, size_t nFlavours,
    const std::vector<tcPDPtr> & quarks,
    const std::vector<tcPDPtr> & antiquarks,
    tcPDVector & outgoing) {
  if (nFlavours < 1 || nFlavours > quarks.size()
      || nFlavours > antiquarks.size()) {
    return false;
  }

  std::array<bool, 5> includeQuark = {{true, true, true, true, true}};
  std::array<bool, 5> includeAntiquark = {{true, true, true, true, true}};
  const std::array<tcPDPtr, 2> incomingPartons = {{
    incoming.first, incoming.second
  }};

  for (tcPDPtr particle : incomingPartons) {
    if (!particle) return false;
    const long id = particle->id();
    if (id == ParticleID::g) continue;
    const size_t flavour = static_cast<size_t>(std::abs(id));
    if (flavour < 1 || flavour > nFlavours) return false;
    const size_t index = flavour - 1;
    bool & crossedLeg = id > 0 ? includeAntiquark[index]
                               : includeQuark[index];
    if (!crossedLeg) return false;
    crossedLeg = false;
  }

  outgoing.clear();
  outgoing.reserve(2*nFlavours);
  for (size_t flavour = 0; flavour < nFlavours; ++flavour) {
    if (includeQuark[flavour]) outgoing.push_back(quarks[flavour]);
    if (includeAntiquark[flavour]) {
      outgoing.push_back(antiquarks[flavour]);
    }
  }
  return true;
}

} // namespace

MEInstanton::MEInstanton()
  : theNQuarkPair(4), theNgluonMax(1), theMaxFinalPartons(0),
    theMultiplicityOption(PoissonMultiplicity),
    theModelOption(PureMultiplicityModel),
    theProcessOption(GluonGluonProcesses),
    theGluonCountingOption(FinalStateGluonCounting),
    theGaussianParamA(5.0), theGaussianParamB(200.0), thePoissonMean(3.0),
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
  if ((theQuarkPairOption == FixedQuarkPairs
       && nQuarkPairs != theNQuarkPair)
      || (theQuarkPairOption != FixedQuarkPairs
          && nQuarkPairs != 4 && nQuarkPairs != 5)) {
    return 0.0;
  }

  // The model probability is attached to the unshifted KKS/PureMultiplicity
  // draw. FixedTotal can add one or two literal final-state gluons after that
  // draw for qg or two-quark initial states.
  const size_t nGluons = currentGluonMultiplicity();
  const size_t maximumGluons = currentGluonMultiplicityMaximum();
  if (nGluons == std::numeric_limits<size_t>::max()
      || maximumGluons == std::numeric_limits<size_t>::max()
      || nGluons > maximumGluons) {
    return 0.0;
  }

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
    truncatedPoissonProbability(nGluons, maximumGluons, meanGluons);
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
  if (theProcessOption > AllProcesses) {
    throw InitException() << "MEInstanton: unknown Processes option."
                          << Exception::abortnow;
  }
  if (theGluonCountingOption > FixedTotalGluonCounting) {
    throw InitException() << "MEInstanton: unknown GluonCounting option."
                          << Exception::abortnow;
  }
  if (theProcessOption == QuarkQuarkProcesses
      && theQuarkPairOption == FixedQuarkPairs && theNQuarkPair < 2) {
    throw InitException()
      << "MEInstanton: Processes QQ needs at least two active flavours."
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
  if (theMaxFinalPartons > 0) {
    const size_t maximumFlavours =
      theQuarkPairOption == FixedQuarkPairs ? theNQuarkPair : 5;
    size_t minimumFinalQuarks = 2*maximumFlavours;
    if (theProcessOption == QuarkGluonProcesses) {
      --minimumFinalQuarks;
    } else if (theProcessOption == QuarkQuarkProcesses
               || theProcessOption == QuarkAntiquarkProcesses) {
      minimumFinalQuarks -= 2;
    }
    if (theMaxFinalPartons < minimumFinalQuarks) {
      throw InitException()
        << "MEInstanton: MaxFinalPartons is too small for every enabled "
        << "zero-mode flavour channel." << Exception::abortnow;
    }
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

  // Each retained multiplicity is a separate BlobME subprocess. For incoming
  // quarks, crossedFermionFinalState removes the conjugate zero-mode leg from
  // the final state. Variable flavour modes register N_f=4 and N_f=5; me2()
  // supplies their relative KKS or legacy probabilities.
  const auto addIncoming = [&](size_t nFlavours,
                               const tcPDPair & incoming) {
    tcPDVector outgoing;
    if (!crossedFermionFinalState(
          incoming, nFlavours, quarks, antiquarks, outgoing)) {
      return;
    }

    const size_t nIncomingGluons =
      static_cast<size_t>(incoming.first->id() == ParticleID::g)
      + static_cast<size_t>(incoming.second->id() == ParticleID::g);
    const size_t shift = gluonMultiplicityShift(nIncomingGluons);
    if (shift == std::numeric_limits<size_t>::max()) return;
    const size_t maximumGluons =
      gluonMultiplicityMaximum(outgoing.size(), nIncomingGluons);
    if (maximumGluons == std::numeric_limits<size_t>::max()) return;

    outgoing.reserve(outgoing.size() + shift + maximumGluons);
    for (size_t n = 0; n < shift; ++n) outgoing.push_back(gluon);
    processMap.insert(std::make_pair(incoming, outgoing));
    for (size_t n = 0; n < maximumGluons; ++n) {
      outgoing.push_back(gluon);
      processMap.insert(std::make_pair(incoming, outgoing));
    }
  };

  const auto addProcessFamilies = [&](size_t nFlavours) {
    // Matchbox constructs the beam-reversed XComb for each non-identical
    // incoming pair, so only one canonical ordering is registered here.
    const bool includeGG =
      theProcessOption == GluonGluonProcesses
      || theProcessOption == AllProcesses;
    const bool includeQG =
      theProcessOption == QuarkGluonProcesses
      || theProcessOption == AllProcesses;
    const bool includeQQ =
      theProcessOption == QuarkQuarkProcesses
      || theProcessOption == AllProcesses;
    const bool includeQQbar =
      theProcessOption == QuarkAntiquarkProcesses
      || theProcessOption == AllProcesses;

    if (includeGG) addIncoming(nFlavours, std::make_pair(gluon, gluon));

    if (includeQG) {
      for (size_t flavour = 0; flavour < nFlavours; ++flavour) {
        addIncoming(nFlavours, std::make_pair(gluon, quarks[flavour]));
        addIncoming(nFlavours, std::make_pair(gluon, antiquarks[flavour]));
      }
    }

    if (includeQQ) {
      // A one-instanton zero-mode set contains only one crossed leg of each
      // sign per flavour, so same-sign equal-flavour initial states are absent.
      for (size_t first = 0; first < nFlavours; ++first) {
        for (size_t second = first + 1; second < nFlavours; ++second) {
          addIncoming(nFlavours,
                      std::make_pair(quarks[first], quarks[second]));
          addIncoming(nFlavours,
                      std::make_pair(antiquarks[first],
                                     antiquarks[second]));
        }
      }
    }

    if (includeQQbar) {
      for (size_t quark = 0; quark < nFlavours; ++quark) {
        for (size_t antiquark = 0; antiquark < nFlavours; ++antiquark) {
          addIncoming(nFlavours,
                      std::make_pair(quarks[quark],
                                     antiquarks[antiquark]));
        }
      }
    }
  };

  if (theQuarkPairOption == FixedQuarkPairs) {
    addProcessFamilies(nQuarkPair());
  } else {
    addProcessFamilies(4);
    addProcessFamilies(5);
  }
  return processMap;
}

size_t MEInstanton::currentNQuarkPairs() const {
  if (mePartonData().size() < 2) return 0;

  // An incoming quark replaces one outgoing conjugate zero-mode leg. Counting
  // fermions on both sides therefore reconstructs the active N_f for gg, qg,
  // and every valid two-quark/antiquark channel.
  size_t zeroModeFermions = 0;
  for (const auto & particle : mePartonData()) {
    if (!particle) return 0;
    if (isInstantonQuark(particle)) {
      ++zeroModeFermions;
    } else if (particle->id() != ParticleID::g) {
      return 0;
    }
  }
  return zeroModeFermions % 2 == 0 ? zeroModeFermions/2 : 0;
}

size_t MEInstanton::currentNIncomingGluons() const {
  if (mePartonData().size() < 2) {
    return std::numeric_limits<size_t>::max();
  }
  size_t result = 0;
  for (size_t position = 0; position < 2; ++position) {
    const auto & particle = mePartonData()[position];
    if (!particle
        || (!isInstantonQuark(particle)
            && particle->id() != ParticleID::g)) {
      return std::numeric_limits<size_t>::max();
    }
    if (particle->id() == ParticleID::g) ++result;
  }
  return result;
}

size_t MEInstanton::currentNFinalGluons() const {
  if (mePartonData().size() < 2) {
    return std::numeric_limits<size_t>::max();
  }
  size_t result = 0;
  for (size_t position = 2; position < mePartonData().size(); ++position) {
    const auto & particle = mePartonData()[position];
    if (!particle
        || (!isInstantonQuark(particle)
            && particle->id() != ParticleID::g)) {
      return std::numeric_limits<size_t>::max();
    }
    if (particle->id() == ParticleID::g) ++result;
  }
  return result;
}

size_t MEInstanton::gluonMultiplicityShift(size_t nIncomingGluons) const {
  if (nIncomingGluons > 2) return std::numeric_limits<size_t>::max();
  return theGluonCountingOption == FixedTotalGluonCounting
           ? 2 - nIncomingGluons
           : 0;
}

size_t MEInstanton::gluonMultiplicityMaximum(
    size_t nOutgoingQuarks, size_t nIncomingGluons) const {
  const size_t shift = gluonMultiplicityShift(nIncomingGluons);
  if (shift == std::numeric_limits<size_t>::max()) return shift;
  if (theMaxFinalPartons == 0) return theNgluonMax;
  if (nOutgoingQuarks > theMaxFinalPartons
      || shift > theMaxFinalPartons - nOutgoingQuarks) {
    return std::numeric_limits<size_t>::max();
  }
  return std::min(
    theNgluonMax, theMaxFinalPartons - nOutgoingQuarks - shift);
}

size_t MEInstanton::currentGluonMultiplicity() const {
  const size_t nIncomingGluons = currentNIncomingGluons();
  const size_t nFinalGluons = currentNFinalGluons();
  if (nIncomingGluons == std::numeric_limits<size_t>::max()
      || nFinalGluons == std::numeric_limits<size_t>::max()) {
    return std::numeric_limits<size_t>::max();
  }
  const size_t shift = gluonMultiplicityShift(nIncomingGluons);
  if (shift == std::numeric_limits<size_t>::max()
      || nFinalGluons < shift) {
    return std::numeric_limits<size_t>::max();
  }
  return nFinalGluons - shift;
}

size_t MEInstanton::currentGluonMultiplicityMaximum() const {
  if (mePartonData().size() < 2) {
    return std::numeric_limits<size_t>::max();
  }
  const size_t nIncomingGluons = currentNIncomingGluons();
  const size_t nFinalGluons = currentNFinalGluons();
  const size_t nFinalParticles = mePartonData().size() - 2;
  if (nIncomingGluons == std::numeric_limits<size_t>::max()
      || nFinalGluons == std::numeric_limits<size_t>::max()
      || nFinalGluons > nFinalParticles) {
    return std::numeric_limits<size_t>::max();
  }
  return gluonMultiplicityMaximum(
    nFinalParticles - nFinalGluons, nIncomingGluons);
}

list<BlobMEBase::ColourConnection> MEInstanton::colourConnections() const {
  list<BlobMEBase::ColourConnection> result;
  if (currentNQuarkPairs() == 0
      || mePartonData().size() != meMomenta().size()) {
    throw Exception()
      << "MEInstanton: inconsistent multiplicity in colourConnections()."
      << Exception::runerror;
  }
  const ColourEndpointSets endpoints =
    collectColourEndpoints(mePartonData());

  if (theColourOption == SimpleColour) {
    // Pair each zero-mode fermion with the endpoint of the same flavour. This
    // remains well-defined after crossing one or two fermions to the initial
    // state. Gluons form two-line singlets, except that an odd leftover gluon
    // is inserted into the final fermion string.
    std::vector<int> unusedFermionSinks = endpoints.fermionSinks;
    std::vector<std::pair<int, int>> fermionPairs;
    fermionPairs.reserve(endpoints.fermionSources.size());
    for (const int source : endpoints.fermionSources) {
      const long flavour =
        std::abs(mePartonData()[static_cast<size_t>(source)]->id());
      const auto sink = std::find_if(
        unusedFermionSinks.begin(), unusedFermionSinks.end(),
        [&](int candidate) {
          return std::abs(
                   mePartonData()[static_cast<size_t>(candidate)]->id())
                 == flavour;
        });
      if (sink == unusedFermionSinks.end()) {
        throw Exception()
          << "MEInstanton: cannot pair the zero-mode fermion colours."
          << Exception::runerror;
      }
      fermionPairs.push_back(std::make_pair(source, *sink));
      unusedFermionSinks.erase(sink);
    }

    std::vector<int> gluons = endpoints.gluons;
    int insertedGluon = -1;
    if (gluons.size() % 2 != 0) {
      const auto outgoing = std::find_if(
        gluons.begin(), gluons.end(), [](int index) { return index >= 2; });
      const auto selected = outgoing != gluons.end() ? outgoing
                                                     : gluons.begin();
      insertedGluon = *selected;
      gluons.erase(selected);
    }

    for (size_t pair = 0; pair < fermionPairs.size(); ++pair) {
      const int source = fermionPairs[pair].first;
      const int sink = fermionPairs[pair].second;
      if (insertedGluon >= 0 && pair + 1 == fermionPairs.size()) {
        result.push_back(connectColourEndpoints(source, insertedGluon));
        result.push_back(connectColourEndpoints(insertedGluon, sink));
      } else {
        result.push_back(connectColourEndpoints(source, sink));
      }
    }
    for (size_t gluon = 0; gluon < gluons.size(); gluon += 2) {
      result.push_back(
        connectColourEndpoints(gluons[gluon], gluons[gluon + 1]));
      result.push_back(
        connectColourEndpoints(gluons[gluon + 1], gluons[gluon]));
    }
    validateColourConnections(result, endpoints);
    return result;
  }

  if (theColourOption == QCDINSPlanarColour) {
    // QCDINS 2.0 models the leading-colour state as N_f open fermion strings.
    // Crossing can move either end of a string into the initial state, but all
    // incoming and outgoing gluons still lie on one of those strings.
    std::vector<int> fermionSinks = endpoints.fermionSinks;
    std::vector<int> gluons = endpoints.gluons;
    randomShuffle(fermionSinks);
    randomShuffle(gluons);

    std::vector<std::vector<int>> strings(endpoints.fermionSources.size());
    for (const int gluon : gluons) {
      const size_t selected =
        static_cast<size_t>(UseRandom::irnd(strings.size()));
      strings[selected].push_back(gluon);
    }

    for (size_t string = 0; string < strings.size(); ++string) {
      int colourSource = endpoints.fermionSources[string];
      for (const int gluon : strings[string]) {
        result.push_back(connectColourEndpoints(colourSource, gluon));
        colourSource = gluon;
      }
      result.push_back(
        connectColourEndpoints(colourSource, fermionSinks[string]));
    }
    validateColourConnections(result, endpoints);
    return result;
  }

  if (theColourOption < RandomColour || theColourOption > Random3Colour) {
    throw Exception() << "MEInstanton: unknown colour-connection option."
                      << Exception::runerror;
  }

  std::vector<int> fixedSources = endpoints.incomingSources;
  std::vector<int> fixedSinks = endpoints.incomingSinks;
  std::vector<int> randomSources = endpoints.outgoingSources;
  std::vector<int> randomSinks = endpoints.outgoingSinks;

  // Random2 moves at most one incoming source and one distinct incoming sink
  // into the random pool. For gg this is exactly the historical source-0,
  // sink-1 choice. With a one-sided qq or qbar-qbar initial state it exposes
  // the one endpoint type that exists.
  if (theColourOption == Random2Colour) {
    int exposedSource = -1;
    if (!fixedSources.empty()) {
      exposedSource = fixedSources.front();
      eraseEndpoint(fixedSources, exposedSource);
      randomSources.push_back(exposedSource);
    }

    const auto exposedSink = std::find_if(
      fixedSinks.begin(), fixedSinks.end(),
      [&](int candidate) { return candidate != exposedSource; });
    if (exposedSink != fixedSinks.end()) {
      randomSinks.push_back(*exposedSink);
      fixedSinks.erase(exposedSink);
    }
  } else if (theColourOption == Random3Colour) {
    fixedSources.clear();
    fixedSinks.clear();
    randomSources = endpoints.sources;
    randomSinks = endpoints.sinks;
  }

  // A quark initial state has only one colour end. Complete the deterministic
  // incoming subgraph with the minimum number of outgoing endpoints, leaving
  // balanced pools for the random map.
  while (fixedSources.size() < fixedSinks.size()) {
    fixedSources.push_back(takeOutgoingEndpoint(randomSources));
  }
  while (fixedSinks.size() < fixedSources.size()) {
    fixedSinks.push_back(takeOutgoingEndpoint(randomSinks));
  }

  const std::vector<int> fixedMap =
    deterministicColourMap(fixedSources, fixedSinks);
  for (size_t source = 0; source < fixedSources.size(); ++source) {
    result.push_back(
      connectColourEndpoints(fixedSources[source], fixedMap[source]));
  }

  const std::vector<int> randomMap =
    randomColourMap(randomSources, randomSinks);
  for (size_t source = 0; source < randomSources.size(); ++source) {
    result.push_back(
      connectColourEndpoints(randomSources[source], randomMap[source]));
  }

  validateColourConnections(result, endpoints);
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
     << ounit(theKKSBottomMass, GeV) << theFermionOverlapInterpolator
     << theProcessOption << theGluonCountingOption << theMaxFinalPartons;
}

void MEInstanton::persistentInput(PersistentIStream & is, int version) {
  is >> theNQuarkPair >> theNgluonMax >> theMultiplicityOption
     >> theModelOption >> theGaussianParamA >> theGaussianParamB
     >> thePoissonMean >> theColourOption >> theInverseRhoInterpolator
     >> theAlphaSInterpolator >> theMeanGluonsInterpolator
     >> theCrossSectionInterpolator >> theScaleOption >> theQuarkPairOption
     >> iunit(theKKSBottomMass, GeV) >> theFermionOverlapInterpolator;
  if (version > 0) {
    is >> theProcessOption >> theGluonCountingOption;
  } else {
    theProcessOption = GluonGluonProcesses;
    theGluonCountingOption = FinalStateGluonCounting;
  }
  if (version > 1) {
    is >> theMaxFinalPartons;
  } else {
    theMaxFinalPartons = 0;
  }
}

DescribeClass<MEInstanton, Herwig::BlobME>
  describeHerwigMEInstanton("Herwig::MEInstanton", "Instantons.so", 2);

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

  static Parameter<MEInstanton, size_t> interfaceMaxFinalPartons(
    "MaxFinalPartons",
    "Maximum number of literal outgoing hard partons. Zero leaves only "
    "NAdditional as the gluon-multiplicity cap.",
    &MEInstanton::theMaxFinalPartons,
    0, 0, 0,
    false, false, Interface::lowerlim
  );

  static Switch<MEInstanton, unsigned int> interfaceProcesses(
    "Processes",
    "The family of incoming partonic states to register.",
    &MEInstanton::theProcessOption,
    GluonGluonProcesses,
    false, false
  );
  static SwitchOption interfaceProcessesGG(
    interfaceProcesses,
    "GG",
    "Register only gluon-gluon initial states.",
    GluonGluonProcesses
  );
  static SwitchOption interfaceProcessesQG(
    interfaceProcesses,
    "QG",
    "Register qg, gq, qbar-g and g-qbar initial states.",
    QuarkGluonProcesses
  );
  static SwitchOption interfaceProcessesQQ(
    interfaceProcesses,
    "QQ",
    "Register distinct-flavour qq and qbar-qbar initial states.",
    QuarkQuarkProcesses
  );
  static SwitchOption interfaceProcessesQQbar(
    interfaceProcesses,
    "QQbar",
    "Register q-qbar and qbar-q initial states.",
    QuarkAntiquarkProcesses
  );
  static SwitchOption interfaceProcessesAll(
    interfaceProcesses,
    "All",
    "Register every supported incoming family.",
    AllProcesses
  );

  static Switch<MEInstanton, unsigned int> interfaceGluonCounting(
    "GluonCounting",
    "How the sampled gluon multiplicity is translated after crossing.",
    &MEInstanton::theGluonCountingOption,
    FinalStateGluonCounting,
    false, false
  );
  static SwitchOption interfaceGluonCountingFinalState(
    interfaceGluonCounting,
    "FinalState",
    "Use the sampled multiplicity as the number of outgoing gluons, as in "
    "the Sherpa implementation.",
    FinalStateGluonCounting
  );
  static SwitchOption interfaceGluonCountingFixedTotal(
    interfaceGluonCounting,
    "FixedTotal",
    "Keep n_g+2 total gluon legs: add one outgoing gluon for qg and two "
    "for two-quark initial states.",
    FixedTotalGluonCounting
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
    "Randomize outgoing lines and complete the incoming subgraph "
    "deterministically.",
    RandomColour
  );
  static SwitchOption interfaceColourConnectionsRandom2(
    interfaceColourConnections,
    "Random2",
    "Expose one incoming source and one distinct sink to the random map.",
    Random2Colour
  );
  static SwitchOption interfaceColourConnectionsRandom3(
    interfaceColourConnections,
    "Random3",
    "Randomize every crossed incoming and outgoing colour endpoint.",
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
