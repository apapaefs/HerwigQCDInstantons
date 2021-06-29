// -*- C++ -*-
#ifndef Herwig_MEInstanton_H
#define Herwig_MEInstanton_H
//
// This is the declaration of the MEInstanton class.
//
#include "Herwig/MatrixElement/HwMEBase.h"
#include "Herwig/MatrixElement/BlobME.h"
#include "Herwig/Utilities/Interpolator.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the MEInstanton class.
 *
 * @see \ref MEInstantonInterfaces "The interfaces"
 * defined for MEInstanton.
 */
class MEInstanton: public Herwig::BlobME {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MEInstanton();

  /**
   * The destructor.
   */
  virtual ~MEInstanton();
  //@}

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

public:

  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix element
   * is given.
   */
  virtual unsigned int orderInAlphaS() const {
    return UINT_MAX;
  }

  /**
   * Return the order in \f$\alpha_{EM}\f$ in which this matrix
   * element is given. Returns 0.
   */
  virtual unsigned int orderInAlphaEW() const {
    return 0;
  }

  /**
   * Return the matrix element for the kinematical configuation
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   */
  virtual double me2() const;

  /**
   * Return the possible processes this matrix element will be able to handle,
   * as a map incoming to outgoing; it is assumed that the number of outgoing
   * partons does not vary.
   */
  virtual multimap<tcPDPair,tcPDVector> processes() const;

  /**
   * Return the colour connections for the process as pairs of id's of
   * external legs connecting colour to anticolour; id's of incoming partons
   * (0 and 1) have the meaning of colour and anti-colour eversed (crossed to
   * the final state).
   */
  virtual list<BlobMEBase::ColourConnection> colourConnections() const;

  /**
   * Return the number of final state particles for the multiplicity set
   * through nAdditional
   */
  virtual size_t nOutgoing() const;

   /**
   * Return the FIXED number of quark pairs
   */

  size_t nQuarkPair() const { return theNQuarkPair; }

   /**
   * Return the CURRENT number of quark pairs
   */

  size_t GetnQuarkPair() const;
  
  /**
   * set the maximum number of gluons
   */

  void ngluonmax(size_t ng) const { ngluon_max = ng; }


  /**
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const;

  /**
   * Return the factorization scale associated with the last set phase space point.
   */
  virtual Energy2 FactorizationScale() const;

 
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
  MEInstanton & operator=(const MEInstanton &);
  
  /**
   * the number of qqbar pairs
   */
  size_t theNQuarkPair;

  /**
   * the maximum number of additional gluons:
   */
  mutable size_t ngluon_max = 1;

  /**
   * whether to use a Gaussian or Poisson parametrizations for the ME2 of the number of gluons:
   * i.e. Gaussian: multiply the ME2 by: exp( -pow((ngluon-GaussianParamA),2)/GaussianParamB)/sqrt(M_PI * GaussianParamB);
   * Poisson: by PoissonMean^ngluon * exp(-PoissonMean)/ngluon!
   */
  unsigned int MultiplicityParametrisation;

  /** 
   * How to model the Matrix Element 
   * 0 = flat with MultiplicityParametrisation giving the kind of distribution for the final state gluons or
   * 1 = according to 1911.09726
   */
  
  unsigned int MEModeling;

  /**
   * the parameters A and B in the GaussianParametrisation of gluon multiplicity:
   */
  double GaussianParamA;
  double GaussianParamB;

  /*
   * The mean of the Poisson distribution used for the gluons if the MultiplicityParametrisation is "Poisson"
   */
  double PoissonMean;

  /**
   * How to perform the colour connections 
   */
  
  unsigned int theColourConnections;

  /**
   * The choice of factorization scale (KKS modeling only)
   */
  unsigned int facscale_option;

  /**
   * How to treat the quark pairs
   */
  unsigned int quarkpair_option;

  /**
   * Setup the interpolators
   */
  void setup_interpolator();

    /**
   * the Interpolators
   */ 
  Interpolator<double, double>::Ptr interpol_invrho;
  Interpolator<double, double>::Ptr interpol_alphasrho;
  Interpolator<double, double>::Ptr interpol_meangluons;
  Interpolator<double, double>::Ptr interpol_sigmahat;

  /**
   *  The hadrons
   */
  mutable tcBeamPtr hadron1;
  mutable tcBeamPtr hadron2;

  /**
   *  momentum fraction
   */
  mutable double x1;
  mutable double x2;

};

}

#endif /* Herwig_MEInstanton_H */
