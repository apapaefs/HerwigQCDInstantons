// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEInstanton class.
//

#include "MEInstanton.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "Herwig/Utilities/Interpolator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include <gsl/gsl_sf_hyperg.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

using namespace Herwig;

namespace {

const double minKKSShat = 10.7;
const double maxKKSShat = 2895.5;

double instantonAction(double chi) {
  const double z = 0.5*(2.0 + chi*chi + chi*sqrt(4.0 + chi*chi));
  const double zminus = z - 1.0/z;
  const double zplus = z + 1.0/z;
  return 3.0*((6.0*z*z - 14.0)/(zminus*zminus) - 17.0/3.0
              - log(z)*((z - 5.0/z)*zplus*zplus/(zminus*zminus*zminus)
                        - 1.0));
}

double instantonActionDerivative(double chi) {
  const double h = 1.e-5*std::max(1.0, fabs(chi));
  return (instantonAction(chi - 2.0*h)
          - 8.0*instantonAction(chi - h)
          + 8.0*instantonAction(chi + h)
          - instantonAction(chi + 2.0*h))/(12.0*h);
}

double solveKKSChi(double target) {
  double lower = 1.0;
  double upper = 3.0;
  double fLower = instantonActionDerivative(lower) - target;
  double fUpper = instantonActionDerivative(upper) - target;
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
    if (fabs(fMiddle) < 1.e-11 || upper - lower < 1.e-11) {
      return middle;
    }
    if (fLower*fMiddle > 0.0) {
      lower = middle;
      fLower = fMiddle;
    } else {
      upper = middle;
      fUpper = fMiddle;
    }
  }
  return 0.5*(lower + upper);
}

double omegaFermion(double chi) {
  const double z = 0.5*(2.0 + chi*chi + chi*sqrt(4.0 + chi*chi));
  const double argument = 1.0 - 1.0/(z*z);
  return (3.0*Constants::pi/8.0)*pow(z, -1.5)
    *gsl_sf_hyperg_2F1(1.5, 1.5, 4.0, argument);
}

double poissonProbability(size_t multiplicity, double mean) {
  if (!std::isfinite(mean) || mean < 0.0) return 0.0;
  if (mean == 0.0) return multiplicity == 0 ? 1.0 : 0.0;
  const double logProbability = -mean + multiplicity*log(mean)
    - lgamma(static_cast<double>(multiplicity) + 1.0);
  return exp(logProbability);
}

double truncatedPoissonProbability(size_t multiplicity, size_t maximum,
                                   double mean) {
  if (!std::isfinite(mean) || mean < 0.0 || multiplicity > maximum) return 0.0;
  if (mean == 0.0) return multiplicity == 0 ? 1.0 : 0.0;

  const double logMean = log(mean);
  double largestLogWeight = -std::numeric_limits<double>::infinity();
  for (size_t n = 0; n <= maximum; ++n) {
    const double logWeight = n*logMean - lgamma(static_cast<double>(n) + 1.0);
    largestLogWeight = std::max(largestLogWeight, logWeight);
  }

  double normalization = 0.0;
  for (size_t n = 0; n <= maximum; ++n) {
    const double logWeight = n*logMean - lgamma(static_cast<double>(n) + 1.0);
    normalization += exp(logWeight - largestLogWeight);
  }
  const double selectedLogWeight = multiplicity*logMean
    - lgamma(static_cast<double>(multiplicity) + 1.0);
  if (!std::isfinite(normalization) || normalization <= 0.0) return 0.0;
  return exp(selectedLogWeight - largestLogWeight)/normalization;
}

double kksFlavourProbability(size_t nQuarkPairs, double omega,
                             double bottomMassRho) {
  if (!std::isfinite(omega) || omega <= 0.0
      || !std::isfinite(bottomMassRho) || bottomMassRho < 0.0) {
    return 0.0;
  }

  const double kappa4 = 0.008;
  const double kappa5 = 0.01;
  double weight4 = kappa4*kappa4*pow(omega, 8);
  double weight5 = 0.0;
  if (bottomMassRho <= 1.0) {
    weight4 *= bottomMassRho*bottomMassRho;
    weight5 = kappa5*kappa5*pow(omega, 10);
  }
  const double normalization = weight4 + weight5;
  if (!std::isfinite(normalization) || normalization <= 0.0) return 0.0;
  if (nQuarkPairs == 4) return weight4/normalization;
  if (nQuarkPairs == 5) return weight5/normalization;
  return 0.0;
}

bool hasNoSelfConnections(const vector<int>& colours,
                          const vector<int>& anticolours) {
  if (colours.size() != anticolours.size()) return false;
  for (size_t i = 0; i < colours.size(); ++i) {
    if (colours[i] == anticolours[i]) return false;
  }
  return true;
}

vector<int> randomColourMap(const vector<int>& colours,
                            const vector<int>& anticolours) {
  if (colours.empty() || colours.size() != anticolours.size()) {
    throw Exception() << "MEInstanton: invalid colour-map inputs."
                      << Exception::runerror;
  }

  for (unsigned int attempt = 0; attempt < 1000; ++attempt) {
    vector<int> result = anticolours;
    for (size_t i = result.size(); i > 1; --i) {
      const size_t selected = static_cast<size_t>(UseRandom::irnd(i));
      std::swap(result[i - 1], result[selected]);
    }
    if (hasNoSelfConnections(colours, result)) return result;
  }

  for (size_t shift = 0; shift < anticolours.size(); ++shift) {
    vector<int> result(anticolours.size());
    for (size_t i = 0; i < anticolours.size(); ++i) {
      result[i] = anticolours[(i + shift) % anticolours.size()];
    }
    if (hasNoSelfConnections(colours, result)) return result;
  }

  throw Exception() << "MEInstanton: no one-to-one colour map without "
                    << "a forbidden self-connection exists."
                    << Exception::runerror;
}

}

MEInstanton::MEInstanton()
  : theNQuarkPair(4), MultiplicityParametrisation(0), MEModeling(0),
    GaussianParamA(5.0), GaussianParamB(200.0), PoissonMean(3.0),
    theColourConnections(0), facscale_option(0), quarkpair_option(0),
    KKSBottomMass(4.18*GeV) {}

MEInstanton::~MEInstanton() {}

IBPtr MEInstanton::clone() const {
  return new_ptr(*this);
}

IBPtr MEInstanton::fullclone() const {
  return new_ptr(*this);
}

void MEInstanton::setup_interpolator() {
  static const array<double,20> hats = {{
      10.7, 11.4, 13.4, 15.7, 22.9, 29.7, 40.8, 56.1, 61.8, 89.6,
      118.0, 174.4, 246.9, 349.9, 496.3, 704.8, 1001.8, 1425.6,
      2030.6, 2895.5 }};
  static const array<double,20> invrho = {{
      0.99, 1.04, 1.16, 1.31, 1.76, 2.12, 2.72, 3.50, 3.64, 4.98,
      6.21, 8.72, 11.76, 15.90, 21.58, 29.37, 40.07, 54.83, 75.21,
      103.4 }};
  static const array<double,20> alphasrho = {{
      0.416, 0.405, 0.382, 0.360, 0.315, 0.293, 0.267, 0.245, 0.223,
      0.206, 0.195, 0.180, 0.169, 0.159, 0.150, 0.142, 0.135, 0.128,
      0.122, 0.117 }};
  static const array<double,20> meangluons = {{
      4.59, 4.68, 4.90, 5.13, 5.44, 6.02, 6.47, 6.92, 7.28, 7.67,
      8.25, 8.60, 9.04, 9.49, 9.93, 10.37, 10.81, 11.26, 11.70,
      12.14 }};
  static const array<double,20> sigmahat = {{
      4.922E9, 3.652E9, 1.671E9, 728.9E6, 85.94E6, 17.25E6,
      2.121E6, 229.0E3, 72.97E3, 2.733E3, 235.4, 6.720, 0.284,
      0.012, 5.112E-4, 21.65E-6, 0.9017E-6, 36.45E-9, 1.419E-9,
      52.07E-12 }};

  array<double,20> omegaferm;
  for (size_t i = 0; i < hats.size(); ++i) {
    const double u = hats[i]/invrho[i];
    const double scaledRho = alphasrho[i]*u/(4.0*Constants::pi);
    const double chi = solveKKSChi(scaledRho);
    omegaferm[i] = omegaFermion(chi);
    if (!std::isfinite(chi) || !std::isfinite(omegaferm[i])
        || omegaferm[i] <= 0.0) {
      throw InitException()
        << "MEInstanton: failed to solve the KKS saddle point at table node "
        << i << "." << Exception::abortnow;
    }
  }

  interpol_invrho = make_InterpolatorPtr(invrho, hats, 1);
  interpol_alphasrho = make_InterpolatorPtr(alphasrho, hats, 1);
  interpol_meangluons = make_InterpolatorPtr(meangluons, hats, 1);
  interpol_sigmahat = make_InterpolatorPtr(sigmahat, hats, 1);
  interpol_omegaferm = make_InterpolatorPtr(omegaferm, hats, 1);
}

Energy2 MEInstanton::scale() const {
  return sHat();
}

Energy2 MEInstanton::FactorizationScale() const {
  if (facscale_option == 1) return sHat();
  return sqr((*interpol_invrho)(sqrt(sHat()/GeV/GeV)))*GeV*GeV;
}

double MEInstanton::me2() const {
  const size_t nQuarkPairs = GetnQuarkPair();
  const size_t baseMultiplicity = 2 + 2*nQuarkPairs;
  if (meMomenta().size() < baseMultiplicity) return 0.0;
  const size_t nGluons = meMomenta().size() - baseMultiplicity;
  if (nGluons > ngluon_max) return 0.0;

  const double shat = sHat()/GeV/GeV;
  if (!std::isfinite(shat) || shat <= 0.0) return 0.0;
  const double sqrtShat = sqrt(shat);
  double result = 1.0;

  if (MEModeling == 0) {
    if (MultiplicityParametrisation == 0) {
      result *= poissonProbability(nGluons, PoissonMean);
    } else if (MultiplicityParametrisation == 1) {
      if (!std::isfinite(GaussianParamB) || GaussianParamB <= 0.0) return 0.0;
      result *= exp(-sqr(static_cast<double>(nGluons) - GaussianParamA)
                    /GaussianParamB)
        /sqrt(Constants::pi*GaussianParamB);
    }
    return std::isfinite(result) && result >= 0.0 ? result : 0.0;
  }

  if (MEModeling != 1 || sqrtShat < minKKSShat || sqrtShat > maxKKSShat) {
    return 0.0;
  }

  if (!lastParticles().first || !lastParticles().second) return 0.0;
  hadron1 = dynamic_ptr_cast<tcBeamPtr>(lastParticles().first->dataPtr());
  hadron2 = dynamic_ptr_cast<tcBeamPtr>(lastParticles().second->dataPtr());
  if (!hadron1 || !hadron2 || !hadron1->pdf() || !hadron2->pdf()) return 0.0;

  x1 = lastX1();
  x2 = lastX2();
  if (!std::isfinite(x1) || !std::isfinite(x2)
      || x1 <= 0.0 || x1 > 1.0 || x2 <= 0.0 || x2 > 1.0) {
    return 0.0;
  }

  tcPDPtr gluon = getParticleData(ParticleID::g);
  if (!gluon) return 0.0;
  const Energy2 factorizationScale = FactorizationScale();
  const double factorizationScaleValue = factorizationScale/GeV/GeV;
  if (!std::isfinite(factorizationScaleValue)
      || factorizationScaleValue <= 0.0) return 0.0;

  const double pdf1Original = hadron1->pdf()->xfx(hadron1, gluon, sHat(), x1)/x1;
  const double pdf2Original = hadron2->pdf()->xfx(hadron2, gluon, sHat(), x2)/x2;
  const double pdf1Reweighted = hadron1->pdf()->xfx(
    hadron1, gluon, factorizationScale, x1)/x1;
  const double pdf2Reweighted = hadron2->pdf()->xfx(
    hadron2, gluon, factorizationScale, x2)/x2;
  if (!std::isfinite(pdf1Original) || pdf1Original <= 0.0
      || !std::isfinite(pdf2Original) || pdf2Original <= 0.0
      || !std::isfinite(pdf1Reweighted) || pdf1Reweighted <= 0.0
      || !std::isfinite(pdf2Reweighted) || pdf2Reweighted <= 0.0) {
    return 0.0;
  }
  const double pdfRatio1 = pdf1Reweighted/pdf1Original;
  const double pdfRatio2 = pdf2Reweighted/pdf2Original;
  if (!std::isfinite(pdfRatio1) || pdfRatio1 <= 0.0
      || !std::isfinite(pdfRatio2) || pdfRatio2 <= 0.0) {
    return 0.0;
  }
  const double pdfRatio = pdfRatio1*pdfRatio2;
  if (!std::isfinite(pdfRatio) || pdfRatio <= 0.0) return 0.0;
  result *= pdfRatio;

  const double meanGluons = (*interpol_meangluons)(sqrtShat);
  const double partonicCrossSection = (*interpol_sigmahat)(sqrtShat);
  if (!std::isfinite(meanGluons) || meanGluons < 0.0
      || !std::isfinite(partonicCrossSection) || partonicCrossSection < 0.0) {
    return 0.0;
  }
  result *= truncatedPoissonProbability(nGluons, ngluon_max, meanGluons);
  result *= partonicCrossSection*2.568E-9;

  const double phaseSpaceJacobian = jacobian();
  if (!std::isfinite(phaseSpaceJacobian) || phaseSpaceJacobian <= 0.0) {
    return 0.0;
  }
  result /= phaseSpaceJacobian;
  result *= 2.0*shat;

  if (quarkpair_option == 1) {
    result *= 0.5;
  } else if (quarkpair_option == 2) {
    const double inverseRho = (*interpol_invrho)(sqrtShat);
    const double omega = (*interpol_omegaferm)(sqrtShat);
    if (!std::isfinite(inverseRho) || inverseRho <= 0.0) return 0.0;
    const double bottomMassRho = (KKSBottomMass/GeV)/inverseRho;
    const double flavourProbability = kksFlavourProbability(
      nQuarkPairs, omega, bottomMassRho);
    if (!std::isfinite(flavourProbability) || flavourProbability < 0.0) {
      return 0.0;
    }
    result *= flavourProbability;
  }

  return std::isfinite(result) && result >= 0.0 ? result : 0.0;
}

void MEInstanton::doinit() {
  if (theNQuarkPair < 1 || theNQuarkPair > 5) {
    throw InitException() << "MEInstanton: NQuarkPair must be between 1 and 5."
                          << Exception::abortnow;
  }
  if (quarkpair_option == 2 && MEModeling != 1) {
    throw InitException()
      << "MEInstanton: QuarkPairs VariableKKS requires MEModeling KKS."
      << Exception::abortnow;
  }
  if (!std::isfinite(PoissonMean) || PoissonMean < 0.0
      || !std::isfinite(GaussianParamA) || GaussianParamA < 0.0
      || !std::isfinite(GaussianParamB) || GaussianParamB <= 0.0) {
    throw InitException() << "MEInstanton: invalid multiplicity parameter."
                          << Exception::abortnow;
  }
  if (!std::isfinite(KKSBottomMass/GeV) || KKSBottomMass < ZERO) {
    throw InitException() << "MEInstanton: KKSBottomMass must be non-negative."
                          << Exception::abortnow;
  }
  ngluonmax(nAdditional());
  setup_interpolator();
}

void MEInstanton::doinitrun() {
  ngluonmax(nAdditional());
}

multimap<tcPDPair,tcPDVector> MEInstanton::processes() const {
  multimap<tcPDPair,tcPDVector> processmap;
  tcPDPtr gluon = getParticleData(ParticleID::g);
  vector<tcPDPtr> quarks;
  vector<tcPDPtr> antiquarks;
  for (int id = 1; id <= 5; ++id) {
    tcPDPtr quark = getParticleData(id);
    quarks.push_back(quark);
    antiquarks.push_back(quark->CC());
  }
  const tcPDPair incoming = make_pair(gluon, gluon);

  const auto addProcesses = [&](size_t nQuarkPairs) {
    tcPDVector outgoing;
    for (size_t i = 0; i < nQuarkPairs; ++i) {
      outgoing.push_back(quarks[i]);
      outgoing.push_back(antiquarks[i]);
    }
    processmap.insert(make_pair(incoming, outgoing));
    for (size_t n = 0; n < ngluon_max; ++n) {
      outgoing.push_back(gluon);
      processmap.insert(make_pair(incoming, outgoing));
    }
  };

  if (quarkpair_option == 0) {
    addProcesses(nQuarkPair());
  } else if (quarkpair_option == 1 || quarkpair_option == 2) {
    addProcesses(4);
    addProcesses(5);
  }
  return processmap;
}

size_t MEInstanton::GetnQuarkPair() const {
  if (quarkpair_option == 0) return theNQuarkPair;

  size_t countedQuarks = 0;
  for (size_t i = 0; i < mePartonData().size(); ++i) {
    if (!mePartonData()[i]) continue;
    const long id = std::abs(mePartonData()[i]->id());
    if (id >= 1 && id <= 5) ++countedQuarks;
  }
  return countedQuarks/2;
}

list<BlobMEBase::ColourConnection> MEInstanton::colourConnections() const {
  list<BlobMEBase::ColourConnection> result;
  const size_t nQuarkPairs = GetnQuarkPair();
  const size_t baseMultiplicity = 2 + 2*nQuarkPairs;
  if (nQuarkPairs == 0 || meMomenta().size() < baseMultiplicity
      || meMomenta().size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
    throw Exception() << "MEInstanton: inconsistent multiplicity in colourConnections()."
                      << Exception::runerror;
  }
  const size_t nGluons = meMomenta().size() - baseMultiplicity;
  const int firstGluon = static_cast<int>(2 + 2*nQuarkPairs);

  if (theColourConnections == 0) {
    BlobMEBase::ColourConnection incoming1;
    BlobMEBase::ColourConnection incoming2;
    incoming1.addColour(0);
    incoming1.addAntiColour(1);
    incoming2.addColour(1);
    incoming2.addAntiColour(0);
    result.push_back(incoming1);
    result.push_back(incoming2);

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

  vector<int> colours;
  vector<int> anticolours;
  if (theColourConnections == 2) {
    colours.push_back(0);
    anticolours.push_back(1);
  } else if (theColourConnections == 3) {
    colours.push_back(0);
    colours.push_back(1);
    anticolours.push_back(0);
    anticolours.push_back(1);
  }
  for (size_t q = 0; q < nQuarkPairs; ++q) {
    colours.push_back(static_cast<int>(2 + 2*q));
    anticolours.push_back(static_cast<int>(3 + 2*q));
  }
  for (size_t g = 0; g < nGluons; ++g) {
    const int index = firstGluon + static_cast<int>(g);
    colours.push_back(index);
    anticolours.push_back(index);
  }
  const vector<int> colourMap = randomColourMap(colours, anticolours);

  for (size_t i = 0; i < colours.size(); ++i) {
    const int colour = colours[i];
    const int anticolour = colourMap[i];
    BlobMEBase::ColourConnection line;
    if (theColourConnections == 1) {
      line.addColour(colour);
      line.addAntiColour(anticolour);
    } else if (theColourConnections == 2) {
      if (anticolour == 1 && colour != 0) {
        line.addColour(anticolour);
        line.addColour(colour);
      } else if (colour == 0 && anticolour != 1) {
        line.addAntiColour(colour);
        line.addAntiColour(anticolour);
      } else if (anticolour != 1 && colour != 0) {
        line.addColour(colour);
        line.addAntiColour(anticolour);
      } else {
        line.addAntiColour(colour);
        line.addColour(anticolour);
      }
    } else if (theColourConnections == 3) {
      if ((anticolour == 1 && colour != 0)
          || (anticolour == 0 && colour != 1)) {
        line.addColour(anticolour);
        line.addColour(colour);
      } else if ((colour == 0 && anticolour != 1)
                 || (colour == 1 && anticolour != 0)) {
        line.addAntiColour(colour);
        line.addAntiColour(anticolour);
      } else if (anticolour != 0 && anticolour != 1
                 && colour != 0 && colour != 1) {
        line.addColour(colour);
        line.addAntiColour(anticolour);
      } else {
        line.addAntiColour(colour);
        line.addColour(anticolour);
      }
    } else {
      throw Exception() << "MEInstanton: unknown colour-connection option."
                        << Exception::runerror;
    }
    result.push_back(line);
  }

  if (theColourConnections == 1) {
    BlobMEBase::ColourConnection incoming1;
    BlobMEBase::ColourConnection incoming2;
    incoming1.addColour(0);
    incoming1.addAntiColour(1);
    incoming2.addColour(1);
    incoming2.addAntiColour(0);
    result.push_back(incoming1);
    result.push_back(incoming2);
  } else if (theColourConnections == 2) {
    BlobMEBase::ColourConnection incoming;
    incoming.addColour(0);
    incoming.addAntiColour(1);
    result.push_back(incoming);
  }
  return result;
}

size_t MEInstanton::nOutgoing() const {
  return quarkpair_option == 0 ? 2*nQuarkPair() : 10;
}

void MEInstanton::persistentOutput(PersistentOStream& os) const {
  os << theNQuarkPair << ngluon_max << MultiplicityParametrisation
     << MEModeling << GaussianParamA << GaussianParamB << PoissonMean
     << theColourConnections << interpol_invrho << interpol_alphasrho
     << interpol_meangluons << interpol_sigmahat << facscale_option
     << quarkpair_option << ounit(KKSBottomMass, GeV) << interpol_omegaferm;
}

void MEInstanton::persistentInput(PersistentIStream& is, int) {
  is >> theNQuarkPair >> ngluon_max >> MultiplicityParametrisation
     >> MEModeling >> GaussianParamA >> GaussianParamB >> PoissonMean
     >> theColourConnections >> interpol_invrho >> interpol_alphasrho
     >> interpol_meangluons >> interpol_sigmahat >> facscale_option
     >> quarkpair_option >> iunit(KKSBottomMass, GeV) >> interpol_omegaferm;
}

DescribeClass<MEInstanton,Herwig::BlobME>
  describeHerwigMEInstanton("Herwig::MEInstanton", "Instantons.so");

void MEInstanton::Init() {
  static ClassDocumentation<MEInstanton> documentation
    ("Phenomenological QCD-instanton matrix element with configurable "
     "multiplicity, flavour, scale and colour models.");

  static Parameter<MEInstanton,size_t> interfaceNQuarkPair
    ("NQuarkPair",
     "The fixed number of quark pairs.",
     &MEInstanton::theNQuarkPair, 4, 1, 5,
     false, false, Interface::limited);

  static Switch<MEInstanton,unsigned int> interfaceColourConnections
    ("ColourConnections",
     "How to connect the colour lines.",
     &MEInstanton::theColourConnections, 0, false, false);
  static SwitchOption interfaceColourConnectionsSimple
    (interfaceColourConnections, "Simple",
     "A deterministic singlet-oriented assignment.", 0);
  static SwitchOption interfaceColourConnectionsRandom
    (interfaceColourConnections, "Random",
     "Random final-state connections with singlet incoming gluons.", 1);
  static SwitchOption interfaceColourConnectionsRandom2
    (interfaceColourConnections, "Random2",
     "Random final-state and one initial-state gluon connection.", 2);
  static SwitchOption interfaceColourConnectionsRandom3
    (interfaceColourConnections, "Random3",
     "Random connections including both initial-state gluons.", 3);

  static Switch<MEInstanton,unsigned int> interfaceFactorizationScale
    ("FactorizationScale",
     "The factorization scale for KKS modelling.",
     &MEInstanton::facscale_option, 0, false, false);
  static SwitchOption interfaceFactorizationScaleInvRho
    (interfaceFactorizationScale, "InvRho",
     "Use inverse-rho squared as the factorization scale.", 0);
  static SwitchOption interfaceFactorizationScaleSHat
    (interfaceFactorizationScale, "sHat",
     "Use the partonic centre-of-mass energy squared.", 1);

  static Switch<MEInstanton,unsigned int> interfaceQuarkPairs
    ("QuarkPairs",
     "How to select the number of quark pairs.",
     &MEInstanton::quarkpair_option, 0, false, false);
  static SwitchOption interfaceNQuarkPairsFixed
    (interfaceQuarkPairs, "Fixed",
     "Use NQuarkPair.", 0);
  static SwitchOption interfaceNQuarkPairsVariable
    (interfaceQuarkPairs, "Variable",
     "Use the legacy equal mixture of four and five pairs.", 1);
  static SwitchOption interfaceNQuarkPairsVariableKKS
    (interfaceQuarkPairs, "VariableKKS",
     "Use KKS scale-, mass- and fermion-factor-dependent four/five-pair weights.",
     2);

  static Parameter<MEInstanton,Energy> interfaceKKSBottomMass
    ("KKSBottomMass",
     "Bottom-quark mass entering the KKS active-flavour condition.",
     &MEInstanton::KKSBottomMass, GeV, 4.18*GeV, ZERO, ZERO,
     false, false, Interface::lowerlim);

  static Switch<MEInstanton,unsigned int> interfaceMultiplicityParametrisation
    ("MultiplicityParametrisation",
     "How to weight gluon multiplicities in PureMultiplicity mode.",
     &MEInstanton::MultiplicityParametrisation, 0, false, false);
  static SwitchOption interfaceMultiplicityParametrisationPoisson
    (interfaceMultiplicityParametrisation, "Poisson",
     "A Poisson distribution with mean PoissonMean.", 0);
  static SwitchOption interfaceMultiplicityParametrisationGaussian
    (interfaceMultiplicityParametrisation, "Gaussian",
     "A Gaussian controlled by GaussianParamA and GaussianParamB.", 1);
  static SwitchOption interfaceMultiplicityParametrisationFlat
    (interfaceMultiplicityParametrisation, "Flat",
     "No multiplicity-dependent factor.", 2);
  static SwitchOption interfaceMultiplicityParametrisationUserDefined
    (interfaceMultiplicityParametrisation, "UserDefined",
     "Reserved for a user-defined factor; currently equivalent to Flat.", 3);

  static Switch<MEInstanton,unsigned int> interfaceMEModeling
    ("MEModeling",
     "How to model the matrix element.",
     &MEInstanton::MEModeling, 0, false, false);
  static SwitchOption interfaceMEModelingPureMultiplicity
    (interfaceMEModeling, "PureMultiplicity",
     "Flat matrix element with a configurable gluon multiplicity factor.", 0);
  static SwitchOption interfaceMEModelingKKS
    (interfaceMEModeling, "KKS",
     "Tabulated model of Khoze, Krauss and Schott (arXiv:1911.09726).", 1);

  static Parameter<MEInstanton,double> interfaceGaussianParamA
    ("GaussianParamA",
     "Centre of the PureMultiplicity Gaussian.",
     &MEInstanton::GaussianParamA, 5.0, 0.0, 0.0,
     false, false, Interface::lowerlim);

  static Parameter<MEInstanton,double> interfaceGaussianParamB
    ("GaussianParamB",
     "Positive denominator in the PureMultiplicity Gaussian exponent.",
     &MEInstanton::GaussianParamB, 200.0, 0.0, 0.0,
     false, false, Interface::lowerlim);

  static Parameter<MEInstanton,double> interfacePoissonMean
    ("PoissonMean",
     "Non-negative mean of the PureMultiplicity Poisson distribution.",
     &MEInstanton::PoissonMean, 3.0, 0.0, 0.0,
     false, false, Interface::lowerlim);
}
