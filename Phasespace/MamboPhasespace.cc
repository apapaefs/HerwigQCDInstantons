// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MamboPhasespace class.
//

#include "MamboPhasespace.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Parameter.h"
#include "Herwig/Utilities/Kinematics.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace Herwig;

MamboPhasespace::MamboPhasespace()
  : _maxweight(10.0), _maxTrials(100000), _a0(10, 0.0), _a1(10, 0.0) {
  _a0[0] = 0.5;
  _a0[1] = 0.375;
  _a0[2] = 0.375;
  _a0[3] = 0.4921875;
  _a0[4] = 0.84375;
  _a0[5] = 1.854492188;
  _a0[6] = 5.0625;
  _a0[7] = 16.58578491;
  _a0[8] = 63.33398438;
  _a0[9] = 275.6161079;

  _a1[0] = 0.5;
  _a1[1] = 0.75;
  _a1[2] = 1.125;
  _a1[3] = 1.96875;
  _a1[4] = 4.21875;
  _a1[5] = 11.12695313;
  _a1[6] = 35.4375;
  _a1[7] = 132.6862793;
  _a1[8] = 570.0058594;
  _a1[9] = 2756.161079;
}

MamboPhasespace::~MamboPhasespace() {}

IBPtr MamboPhasespace::clone() const {
  return new_ptr(*this);
}

IBPtr MamboPhasespace::fullclone() const {
  return new_ptr(*this);
}

double MamboPhasespace::generateKinematics(vector<Lorentz5Momentum>& momenta,
                                           Energy centreOfMassEnergy,
                                           const double* randomNumbers) const {
  (void)randomNumbers;
  if (!std::isfinite(_maxweight) || _maxweight <= 0.0 || _maxTrials == 0) {
    throw Exception() << "MamboPhasespace: MaxWeight and MaxTrials must be "
                      << "strictly positive." << Exception::runerror;
  }

  for (unsigned int trial = 0; trial < _maxTrials; ++trial) {
    const double weight = calculateMomentum(momenta, centreOfMassEnergy);
    if (!std::isfinite(weight) || weight <= 0.0) continue;

    const double tolerance = 64.0*std::numeric_limits<double>::epsilon()
      *std::max(1.0, _maxweight);
    if (weight > _maxweight + tolerance) {
      throw Exception()
        << "MamboPhasespace: generated weight " << weight
        << " exceeds MaxWeight=" << _maxweight
        << ". Increase MaxWeight before generating events."
        << Exception::runerror;
    }
    if (weight >= _maxweight*UseRandom::rnd()) return 1.0;
  }

  throw Exception() << "MamboPhasespace: failed to accept a phase-space point "
                    << "within MaxTrials=" << _maxTrials << "."
                    << Exception::runerror;
}

double MamboPhasespace::calculateMomentum(
    vector<Lorentz5Momentum>& partMomenta, Energy centreOfMassEnergy) const {
  if (partMomenta.size() < 4) {
    throw Exception() << "MamboPhasespace requires at least two final-state "
                      << "particles." << Exception::runerror;
  }
  const size_t nFinal = partMomenta.size() - 2;
  const double centreOfMassValue = centreOfMassEnergy/GeV;
  if (!std::isfinite(centreOfMassValue) || centreOfMassValue <= 0.0) {
    throw Exception() << "MamboPhasespace: non-positive or non-finite "
                      << "centre-of-mass energy." << Exception::runerror;
  }

  vector<Lorentz5Momentum> momenta(nFinal);
  vector<Energy> requestedMasses(nFinal, ZERO);
  Energy totalMass = ZERO;
  Energy2 totalMassSquared = ZERO;
  for (size_t i = 0; i < nFinal; ++i) {
    requestedMasses[i] = partMomenta[i + 2].mass();
    if (!std::isfinite(requestedMasses[i]/GeV) || requestedMasses[i] < ZERO) {
      throw Exception() << "MamboPhasespace: invalid requested final-state mass."
                        << Exception::runerror;
    }
    momenta[i].setMass(requestedMasses[i]);
    totalMass += requestedMasses[i];
    totalMassSquared += requestedMasses[i]*requestedMasses[i];
  }
  if (!(centreOfMassEnergy > totalMass)) {
    throw Exception() << "MamboPhasespace: final-state masses reach or exceed "
                      << "the available centre-of-mass energy."
                      << Exception::runerror;
  }

  const double n = static_cast<double>(nFinal);
  const Energy2 initialRadicand =
    (n*sqr(centreOfMassEnergy) - sqr(totalMass))/(n*n*(n - 1.0));
  if (!std::isfinite(initialRadicand/GeV/GeV) || initialRadicand <= ZERO) {
    throw Exception() << "MamboPhasespace: invalid initial MAMBO scale."
                      << Exception::runerror;
  }
  const Energy initialScale = sqrt(initialRadicand) - totalMass/n;
  Energy oldScale = (2.0/3.0)*initialScale;
  if (!std::isfinite(oldScale/GeV) || oldScale <= ZERO) {
    throw Exception() << "MamboPhasespace: failed to initialize the MAMBO scale."
                      << Exception::runerror;
  }

  bool scaleConverged = false;
  Energy scale = oldScale;
  for (unsigned int iteration = 0; iteration < 100; ++iteration) {
    Energy sumF = ZERO;
    long double sumFPrime = 0.0;
    Energy sumFFPrime = ZERO;
    Energy2 minusSumF2 = ZERO;
    for (size_t i = 0; i < nFinal; ++i) {
      const long double ratio = fabs(requestedMasses[i]/oldScale);
      Energy f = ZERO;
      long double fPrime = 0.0;
      if (ratio == 0.0) {
        f = 2.0*oldScale;
        fPrime = 2.0;
      } else {
        long double besselRatio = 0.0;
        long double besselDerivative = 0.0;
        BesselFns(ratio, besselRatio, besselDerivative);
        if (!std::isfinite(static_cast<double>(besselRatio))
            || !std::isfinite(static_cast<double>(besselDerivative))) {
          throw Exception() << "MamboPhasespace: non-finite Bessel expansion."
                            << Exception::runerror;
        }
        f = oldScale*(2.0 + ratio*besselRatio);
        fPrime = 2.0 - ratio*ratio*besselDerivative;
      }
      sumF += f;
      sumFPrime += fPrime;
      sumFFPrime += f*fPrime;
      minusSumF2 -= f*f;
    }

    const Energy derivative = 2.0*(sumF*sumFPrime - sumFFPrime);
    const Energy2 residual = sumF*sumF + minusSumF2 + totalMassSquared
      - sqr(centreOfMassEnergy);
    if (!std::isfinite(derivative/GeV)
        || fabs(derivative/GeV) <= std::numeric_limits<double>::min()) {
      throw Exception() << "MamboPhasespace: singular scale iteration."
                        << Exception::runerror;
    }
    scale = oldScale - residual/derivative;
    if (!std::isfinite(scale/GeV) || scale <= ZERO) {
      throw Exception() << "MamboPhasespace: invalid scale iteration."
                        << Exception::runerror;
    }
    const double relativeScaleChange = fabs((scale - oldScale)/GeV)
      /std::max(1.0, fabs(scale/GeV));
    if (relativeScaleChange <= 1.e-12) {
      scaleConverged = true;
      break;
    }
    oldScale = scale;
  }
  if (!scaleConverged) {
    throw Exception() << "MamboPhasespace: scale iteration did not converge."
                      << Exception::runerror;
  }

  vector<long double> alpha(nFinal);
  vector<long double> maximumU(nFinal);
  vector<long double> maximumV(nFinal);
  for (size_t i = 0; i < nFinal; ++i) {
    alpha[i] = 2.0*(requestedMasses[i]/scale);
    const long double xu = (1.0 - alpha[i]
      + sqrt(1.0 + alpha[i]*alpha[i]))/2.0;
    const long double xv = (3.0 - alpha[i]
      + sqrt(9.0 + 4.0*alpha[i] + alpha[i]*alpha[i]))/2.0;
    maximumU[i] = exp(-xu/2.0)*pow(xu*(xu + alpha[i]), 0.25L);
    maximumV[i] = xv*exp(-xv/2.0)*pow(xv*(xv + alpha[i]), 0.25L);
    if (!std::isfinite(static_cast<double>(alpha[i]))
        || !std::isfinite(static_cast<double>(maximumU[i]))
        || !std::isfinite(static_cast<double>(maximumV[i]))
        || maximumU[i] <= 0.0 || maximumV[i] <= 0.0) {
      throw Exception() << "MamboPhasespace: invalid rejection envelope."
                        << Exception::runerror;
    }
  }

  vector<Lorentz5Momentum> trialMomenta(nFinal);
  Lorentz5Momentum totalTrialMomentum;
  long double rescaling = 0.0;
  bool totalAccepted = false;
  for (unsigned int totalTrial = 0; totalTrial < _maxTrials; ++totalTrial) {
    totalTrialMomentum = Lorentz5Momentum();
    bool particlesAccepted = true;
    for (size_t i = 0; i < nFinal; ++i) {
      long double sampledX = 0.0;
      bool particleAccepted = false;
      for (unsigned int particleTrial = 0;
           particleTrial < _maxTrials; ++particleTrial) {
        const long double u = UseRandom::rnd()*maximumU[i];
        const long double v = UseRandom::rnd()*maximumV[i];
        if (u <= 0.0) continue;
        sampledX = v/u;
        const long double bound = exp(-sampledX)
          *sqrt(sampledX*(sampledX + alpha[i]));
        if (std::isfinite(static_cast<double>(sampledX))
            && std::isfinite(static_cast<double>(bound))
            && sampledX >= 0.0 && bound >= 0.0 && u*u <= bound) {
          particleAccepted = true;
          break;
        }
      }
      if (!particleAccepted) {
        particlesAccepted = false;
        break;
      }

      double cosine = 0.0;
      double phi = 0.0;
      Kinematics::generateAngles(cosine, phi);
      const double sine = sqrt(std::max(0.0, 1.0 - cosine*cosine));
      const Energy momentumMagnitude = scale
        *sqrt(sampledX*(sampledX + alpha[i]));
      trialMomenta[i] = Lorentz5Momentum(
        momentumMagnitude*sine*sin(phi),
        momentumMagnitude*sine*cos(phi),
        momentumMagnitude*cosine,
        requestedMasses[i] + scale*sampledX,
        requestedMasses[i]);
      totalTrialMomentum += trialMomenta[i];
    }
    if (!particlesAccepted) continue;

    totalTrialMomentum.rescaleMass();
    const double trialMassSquared = totalTrialMomentum.mass2()/GeV/GeV;
    if (!std::isfinite(trialMassSquared) || trialMassSquared <= 0.0) continue;
    rescaling = sqrt(sqr(centreOfMassEnergy)/totalTrialMomentum.mass2());
    if (std::isfinite(static_cast<double>(rescaling))
        && rescaling > 0.0 && rescaling <= 1.0) {
      totalAccepted = true;
      break;
    }
  }
  if (!totalAccepted) {
    throw Exception() << "MamboPhasespace: total-momentum rejection exceeded "
                      << "MaxTrials=" << _maxTrials << "."
                      << Exception::runerror;
  }

  const Energy trialMass = totalTrialMomentum.mass();
  const Energy boostDenominator = totalTrialMomentum.e() + trialMass;
  if (!std::isfinite(trialMass/GeV) || trialMass <= ZERO
      || !std::isfinite(boostDenominator/GeV) || boostDenominator <= ZERO) {
    return 0.0;
  }

  vector<Lorentz5Momentum> massScaledMomenta(nFinal);
  vector<Energy2> spatialMomentumSquared(nFinal, ZERO);
  for (size_t i = 0; i < nFinal; ++i) {
    const Energy projectedEnergy = (trialMomenta[i]*totalTrialMomentum)/trialMass;
    const long double boostFactor = (projectedEnergy + trialMomenta[i].e())
      /boostDenominator;
    const Energy px = trialMomenta[i].x()
      - totalTrialMomentum.x()*boostFactor;
    const Energy py = trialMomenta[i].y()
      - totalTrialMomentum.y()*boostFactor;
    const Energy pz = trialMomenta[i].z()
      - totalTrialMomentum.z()*boostFactor;
    massScaledMomenta[i] = Lorentz5Momentum(
      rescaling*px, rescaling*py, rescaling*pz,
      rescaling*projectedEnergy);
    spatialMomentumSquared[i] = sqr(massScaledMomenta[i].e())
      - rescaling*rescaling*sqr(requestedMasses[i]);
    const double spatialValue = spatialMomentumSquared[i]/GeV/GeV;
    if (!std::isfinite(spatialValue)) return 0.0;
    if (spatialValue < 0.0) {
      if (spatialValue < -1.e-12) return 0.0;
      spatialMomentumSquared[i] = ZERO;
    }
  }

  long double xi = 1.0;
  bool xiConverged = false;
  for (unsigned int iteration = 0; iteration < 100; ++iteration) {
    Energy residual = -centreOfMassEnergy;
    Energy derivative = ZERO;
    for (size_t i = 0; i < nFinal; ++i) {
      const Energy energy = sqrt(xi*xi*spatialMomentumSquared[i]
                                 + sqr(requestedMasses[i]));
      if (!std::isfinite(energy/GeV) || energy <= ZERO) return 0.0;
      residual += energy;
      derivative += spatialMomentumSquared[i]/energy;
    }
    const Energy denominator = xi*derivative;
    if (!std::isfinite(denominator/GeV) || denominator <= ZERO) return 0.0;
    const long double nextXi = xi - residual/denominator;
    if (!std::isfinite(static_cast<double>(nextXi)) || nextXi <= 0.0) {
      return 0.0;
    }
    if (fabs(nextXi - xi) <= 1.e-12) {
      xi = nextXi;
      xiConverged = true;
      break;
    }
    xi = nextXi;
  }
  if (!xiConverged) return 0.0;

  for (size_t i = 0; i < nFinal; ++i) {
    const Energy energy = sqrt(xi*xi*spatialMomentumSquared[i]
                               + sqr(requestedMasses[i]));
    momenta[i] = Lorentz5Momentum(
      xi*massScaledMomenta[i].x(), xi*massScaledMomenta[i].y(),
      xi*massScaledMomenta[i].z(), energy, requestedMasses[i]);
  }

  Lorentz5Momentum totalMomentum;
  for (size_t i = 0; i < nFinal; ++i) totalMomentum += momenta[i];
  const double conservationTolerance = 1.e-8*std::max(1.0, centreOfMassValue);
  if (!std::isfinite(totalMomentum.e()/GeV)
      || fabs(totalMomentum.x()/GeV) > conservationTolerance
      || fabs(totalMomentum.y()/GeV) > conservationTolerance
      || fabs(totalMomentum.z()/GeV) > conservationTolerance
      || fabs((totalMomentum.e() - centreOfMassEnergy)/GeV)
         > conservationTolerance) {
    return 0.0;
  }

  double energyRatioProduct = 1.0;
  Energy massTermBefore = ZERO;
  Energy massTermAfter = ZERO;
  for (size_t i = 0; i < nFinal; ++i) {
    if (massScaledMomenta[i].e() <= ZERO || momenta[i].e() <= ZERO) return 0.0;
    energyRatioProduct *= massScaledMomenta[i].e()/momenta[i].e();
    massTermBefore += sqr(requestedMasses[i])/massScaledMomenta[i].e();
    massTermAfter += sqr(requestedMasses[i])/momenta[i].e();
  }
  const Energy denominator = centreOfMassEnergy - massTermAfter;
  if (!std::isfinite(denominator/GeV) || denominator <= ZERO) return 0.0;
  const double weight = pow(xi, 3.0*n - 3.0)*energyRatioProduct
    *(centreOfMassEnergy - rescaling*rescaling*massTermBefore)/denominator;
  if (!std::isfinite(weight) || weight <= 0.0) return 0.0;

  for (size_t i = 0; i < nFinal; ++i) partMomenta[i + 2] = momenta[i];
  return weight;
}

double MamboPhasespace::generateTwoToNKinematics(
    const double* randomNumbers, vector<Lorentz5Momentum>& momenta) {
  if (!lastXCombPtr()) {
    throw Exception() << "MamboPhasespace: no XComb is available."
                      << Exception::runerror;
  }
  const Energy2 shat = lastXCombPtr()->lastSHat();
  if (!std::isfinite(shat/GeV/GeV) || shat <= ZERO) return 0.0;
  return generateKinematics(momenta, sqrt(shat), randomNumbers);
}

void MamboPhasespace::persistentOutput(PersistentOStream& os) const {
  os << _maxweight << _maxTrials << _a0 << _a1;
}

void MamboPhasespace::persistentInput(PersistentIStream& is, int) {
  is >> _maxweight >> _maxTrials >> _a0 >> _a1;
}

DescribeClass<MamboPhasespace,MatchboxPhasespace>
  describeHerwigMamboPhasespace("Herwig::MamboPhasespace", "Herwig.so");

void MamboPhasespace::Init() {
  static ClassDocumentation<MamboPhasespace> documentation
    ("Stochastic, non-invertible MAMBO phase-space generation with internal "
     "accept/reject unweighting. The supplied integration coordinate is a "
     "dummy; use InvertiblePhasespace when deterministic inversion is needed.");

  static Parameter<MamboPhasespace,double> interfaceMaximumWeight
    ("MaxWeight",
     "Strict upper bound for MAMBO accept/reject unweighting.",
     &MamboPhasespace::_maxweight, 10.0, 1.e-12, 1.e12,
     false, false, Interface::limited);

  static Parameter<MamboPhasespace,unsigned int> interfaceMaximumTrials
    ("MaxTrials",
     "Maximum attempts allowed in each MAMBO rejection loop.",
     &MamboPhasespace::_maxTrials, 100000, 1, 100000000,
     false, false, Interface::limited);
}
