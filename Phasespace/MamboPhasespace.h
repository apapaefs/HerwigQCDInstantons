// -*- C++ -*-
//
// MamboPhasespace.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MamboPhasespace_H
#define Herwig_MamboPhasespace_H
//
// This is the declaration of the MamboPhasespace class.
//

#include "Herwig/MatrixElement/Matchbox/Phasespace/MatchboxPhasespace.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief MamboPhasespace implements stochastic, non-invertible MAMBO
 * phase-space generation.
 *
 */
class MamboPhasespace: public MatchboxPhasespace {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  //@{
  /**
   * The default constructor.
   */
  MamboPhasespace();

  /**
   * The destructor.
   */
  virtual ~MamboPhasespace();
  //@}


public:

  /**
   * Generate a phase space point and return its weight.
   */
  virtual double generateTwoToNKinematics(const double*,
					  vector<Lorentz5Momentum>& momenta);

  /**
   * Return the number of random numbers required to produce a given
   * multiplicity final state.
   */
  virtual int nDimPhasespace(int nFinal) const {
    (void)nFinal;
    return 1;
  }


public:

  /**
   * Return true, if this phase space generator is invertible
   */
  virtual bool isInvertible() const { return false; }

 

  /** falling factorial */
  long double falling_factorial(unsigned int n, unsigned int m) const {

    if (n>m) { 
      return n * falling_factorial(n-1, m);
    }
    else {
      return 1.0;
    }
    
    
    //return factorial(n)/factorial(m);
    //return (n < 2) ? 1.0 : n * factorial(n - 1);
  }

private:

 
  /**
   * Return rho
   */
  double rho(Energy M, Energy N, Energy m) const {
    return sqrt((sqr(M)-sqr(N+m))*(sqr(M)-sqr(N-m)))/(8.*sqr(M));
  }


  /**
   * Generate momenta in the CMS
   */
  double generateKinematics(vector<Lorentz5Momentum>& P,
			    Energy Ecm,
			    const double* r) const;

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

    /** @name Bessel Functions.*/
  //@{
  /**
   * Compute the values \f$K_0(x)/K_1(x)\f$ and it's derivative using
   * asymptotic expansion for large x values.
   * @param x The argument
   * @param f The value of the ratio
   * @param fp The value of the derivative ratio
   */
  void BesselFns(const long double x,
		 long double & f, long double & fp) const {
    assert(x>=0.);
    if( x < 10. ) {
      f = BesselK0(x)/BesselK1(x);
      fp = ( sqr(f)*x + f - x )/x;
    }
    else
      BesselIExpand(-x, f, fp);
  } 
  
  /**
   * Compute the values \f$I_0(x)/I_1(x)\f$ and it's derivative using
   * asymptotic expansion.
   * @param x The argument
   * @param f The value of the ratio
   * @param fp The value of the derivative ratio
   */
  void BesselIExpand(const long double x,
		     long double & f, long double & fp) const {
    long double y = 1./x;
    f = 1.+ y*(_a0[0] + y*(_a0[1] + y*(_a0[2] + y*(_a0[3] 
        + y*(_a0[4] + y*(_a0[5] + y*(_a0[6] + y*(_a0[7] 
        + y*(_a0[8] + y*_a0[9] )))))))));
    fp = -y*y*(_a1[0] + y*(_a1[1] + y*(_a1[2] + y*(_a1[3] 
        + y*(_a1[4] + y*(_a1[5] + y*(_a1[6] + y*(_a1[7] 
        + y*(_a1[8] + y*_a1[9] )))))))));
  }

  /**
   * Modified Bessel function of first kind \f$I_0(x)\f$.
   *@param x Argument of Bessel Function 
   **/
  long double BesselI0(const long double x) const {
    long double y,ans;
    if(x < 3.75) {
      y = sqr(x/3.75);
      ans = 1. + y*(3.5156229 + y*(3.0899424 + y*(1.2067492 
          + y*(0.2659732 + y*(0.0360768+y*0.0045813)))));
    }
    else {
      y = (3.75/x);
      ans = (exp(x)/sqrt(x))*(0.39894228 + y*(0.01328592 
          + y*(0.00225319 + y*(-0.00157565 + y*(0.00916281 
          + y*(-0.02057706+y*(0.02635537+y*(-0.01647633+y*0.00392377))))))));
    }
    return ans;
  }
  
  /**
   *  Modified Bessel function of first kind \f$I_1(x)\f$.
   *@param x Argument of Bessel Function 
   **/
  long double BesselI1(const long double x) const {
    long double y,ans;
    if(x < 3.75) {
      y = sqr(x/3.75);
      ans = x*(0.5 + y*(0.87890594 + y*(0.51498869 + y*(0.15084934 
          + y*(0.02658733 + y*(0.00301532 + y*0.00032411))))));
    }
    else {
      y = 3.75/x;
      ans = (0.39894228 + y*(-0.03988024 + y*(-0.00362018 
          + y*(0.00163801 + y*(-0.01031555 + y*(0.02282967 
	  + y*(-0.02895312 + y*(0.01787654-y*0.00420059))))))))*(exp(x)/sqrt(x));
    }
    return ans;
  }
  
  /**
   * Modified Bessel function of second kind \f$K_0(x)\f$.
   * @param x Argument of Bessel Function 
   **/
  long double BesselK0(const long double x) const {
    long double y,ans;
    if(x <= 2.0) {
      y = x*x/4.0;
      ans = -log(x/2.0)*BesselI0(x) - 0.57721566 
          + y*(0.42278420 + y*(0.23069756 
          + y*(0.03488590 + y*(0.00262698 + y*(0.00010750+y*0.00000740)))));
    }
    else {
      y = 2.0/x;
      ans = (1.25331414 + y*(-0.07832358 + y*(+0.02189568 
          + y*(-0.01062446 + y*(0.00587872 
          + y*(-0.00251540 + y*0.00053208))))))*(exp(-x)/sqrt(x));
    }
    return ans;
  }
  
  /**
   * Modified Bessel function of second kind \f$K_1(x)\f$.
   * @param x Argument of Bessel Function 
   **/
  long double BesselK1(const long double x) const  {
    long double y,ans;
    if(x <= 2.0) {
      y = x*x/4.;
      ans = log(x/2.)*BesselI1(x) + (1./x)*(1. + y*(0.15443144 
          + y*(-0.67278579 + y*(-0.18156897 
          + y*(-0.01919402+y*(-0.00110404-(y*0.00004686)))))));
    }
    else {
      y = 2./x;
      ans = (exp(-x)/sqrt(x))*(1.25331414 + y*(0.23498619 
          + y*(-0.03655620 + y*(0.01504268 + y*(-0.00780353 
          + y*(0.00325614+y*(-0.00068245)))))));
    }
    return ans;
  }
  //@}
  
   /**
    *Set array of mometum to particles
     *@param mom   Momentum set to be distributed over phase-space
     *@param comEn The mass of the decaying particle
     *@return The weight of the configuration
     **/
  double calculateMomentum(vector<Lorentz5Momentum> & mom,
			   Energy comEn) const;

  

  /**
   * Upper bound used to unweight MAMBO configurations.
   */
  double _maxweight;

  /**
   * Maximum number of attempts in stochastic rejection loops.
   */
  unsigned int _maxTrials;

  /**
   * Store coefficents for aysymptotic expansion of \f$\frac{I_0}{I_1}\f$
   */
  vector<double> _a0;

  /**
   * Store data for aysymptotic expansion of the first derivative
   * \f$\frac{I_0}{I_1}\f$.
   */
  vector<double> _a1;

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MamboPhasespace & operator=(const MamboPhasespace &);

};

}

#endif /* Herwig_MamboPhasespace_H */
