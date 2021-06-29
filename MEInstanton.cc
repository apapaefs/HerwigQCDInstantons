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

using namespace Herwig;

MEInstanton::MEInstanton() : theNQuarkPair(4), theColourConnections(0), MultiplicityParametrisation(0), MEModeling(0), GaussianParamA(5), GaussianParamB(200.), PoissonMean(3), facscale_option(0), quarkpair_option(0) {}

MEInstanton::~MEInstanton() {}

IBPtr MEInstanton::clone() const {
  return new_ptr(*this);
}

IBPtr MEInstanton::fullclone() const {
  return new_ptr(*this);
}

void MEInstanton::setup_interpolator()  {
  static const array<double,20> hats = {{ 10.7, 11.4, 13.4, 15.7, 22.9, 29.7, 40.8, 56.1, 61.8, 89.6, 118.0, 174.4, 246.9, 349.9, 496.3, 704.8, 1001.8, 1425.6, 2030.6, 2895.5 }}; // the square root of s-hat
  static const array<double,20> invrho = {{ 0.99, 1.04, 1.16, 1.31, 1.76, 2.12, 2.72, 3.50, 3.64, 4.98, 6.21, 8.72, 11.76, 15.90, 21.58, 29.37, 40.07, 54.83, 75.21, 103.4 }}; // the 1/rho
  static const array<double,20> alphasrho = {{ 0.416, 0.405, 0.382, 0.360, 0.315, 0.293, 0.267, 0.245, 0.223, 0.206, 0.195, 0.180, 0.169, 0.159, 0.150, 0.142, 0.135, 0.128, 0.122, 0.117}}; //the alphaS(1/rho)
  static const array<double,20> meangluons = {{ 4.59, 4.68, 4.90, 5.13, 5.44, 6.02, 6.47, 6.92, 7.28, 7.67, 8.25, 8.60, 9.04, 9.49, 9.93, 10.37, 10.81, 11.26, 11.70, 12.14}}; //the mean number of gluons
  static const array<double,20> sigmahat = {{4.992E9, 3.652E9, 1.671E9, 728.9E6, 85.94E6, 17.25E6, 2.121E6, 229.0E3, 72.97E3, 2.733E3, 235.4, 6.720, 0.284, 0.012, 5.112E-4, 21.65E-6, 0.9017E-6, 36.45E-9, 1.419E-9, 52.07E-12}};
  // static const array<double,20> sigmahat = {{1.45813e10, 1.05266e10, 4.51405e9, 1.85274e9, 1.76977e8, 3.55261e7, 3.99487e6, 397757., 113207., 3876.77, 333.886, 8.68739,0.348676, 0.0140647, 0.000571738, 0.0000232145, 9.29353e-7, 3.61946e-8, 1.36042e-9, 4.83086e-11}};


  
  //create the interpolators:
  interpol_invrho = make_InterpolatorPtr(invrho, hats, 1);
  interpol_alphasrho = make_InterpolatorPtr(alphasrho, hats, 1);
  interpol_meangluons = make_InterpolatorPtr(meangluons, hats, 1);
  interpol_sigmahat = make_InterpolatorPtr(sigmahat, hats, 1);

  /* test inteprolators */
  bool test_interpolators = true;
  if(test_interpolators == true) {
    double sqrt_hats = 10.7;
    cout << "sqrt_hats, invrho, alphasrho, meangluons, sigmahat=" << sqrt_hats << "\t" << (*interpol_invrho)(sqrt_hats) << "\t" << (*interpol_alphasrho)(sqrt_hats) << "\t" << (*interpol_meangluons)(sqrt_hats) << "\t" << (*interpol_sigmahat)(sqrt_hats) << endl;
    sqrt_hats = 29.7;
    cout << "sqrt_hats, invrho, alphasrho, meangluons, sigmahat=" << sqrt_hats << "\t" << (*interpol_invrho)(sqrt_hats) << "\t" << (*interpol_alphasrho)(sqrt_hats) << "\t" << (*interpol_meangluons)(sqrt_hats) << "\t" << (*interpol_sigmahat)(sqrt_hats) << endl;
    sqrt_hats = 704.8;
    cout << "sqrt_hats, invrho, alphasrho, meangluons, sigmahat=" << sqrt_hats << "\t" << (*interpol_invrho)(sqrt_hats) << "\t" << (*interpol_alphasrho)(sqrt_hats) << "\t" << (*interpol_meangluons)(sqrt_hats) << "\t" << (*interpol_sigmahat)(sqrt_hats) << endl;
    sqrt_hats = 900.;
    cout << "sqrt_hats, invrho, alphasrho, meangluons, sigmahat=" << sqrt_hats << "\t" << (*interpol_invrho)(sqrt_hats) << "\t" << (*interpol_alphasrho)(sqrt_hats) << "\t" << (*interpol_meangluons)(sqrt_hats) << "\t" << (*interpol_sigmahat)(sqrt_hats) << endl;
  }
}

Energy2 MEInstanton::scale() const {
  //return sqr((*interpol_invrho)(sHat()/GeV/GeV))*GeV*GeV;
  return sHat();
}

Energy2 MEInstanton::FactorizationScale() const {
  if(facscale_option == 0) {
    //cout << "sqrt(sHat), InvRho = " << sqrt(sHat()/GeV/GeV) << "\t" << (*interpol_invrho)(sqrt(sHat()/GeV/GeV)) << endl;
    return sqr((*interpol_invrho)(sqrt(sHat()/GeV/GeV)))*GeV*GeV;
  } else if(facscale_option == 1) {
    return sHat();
  }
  return sqr((*interpol_invrho)(sqrt(sHat()/GeV/GeV)))*GeV*GeV;
}

double MEInstanton::me2() const {

 
  // the square of the matrix element
  double mesq = 1.;

  // get the number of gluons
  int ngluon = (meMomenta().size()-2-GetnQuarkPair()*2);
  //cout << "number of additional gluons = " << ngluon << endl;

  // get the total number of outgoing particles
  int noutgoing = ngluon + GetnQuarkPair()*2;

  double interpolated_meangluons;
  double interpolated_sigmahat;

  //a double of the sqrt(sHat)
  double sqrt_hats = sqrt(sHat()/GeV/GeV);

  
  //cout << "scale() = " << scale()/GeV/GeV << endl;
  // cout << "sqr((*interpol_invrho)(sHat()/GeV/GeV))*GeV*GeV = " << sqr((*interpol_invrho)(sHat()/GeV/GeV)) << endl;
  /* 
   * multiply by an appropriate factor 
   * to take into account the multiplicity parametrisation
   */
  if(MEModeling == 0) { 
    if(MultiplicityParametrisation==0) {
      mesq *= pow(PoissonMean, ngluon) * exp(-PoissonMean)/factorial(ngluon);
    } else if(MultiplicityParametrisation==1) {
      mesq *= exp( -pow((ngluon-GaussianParamA),2)/GaussianParamB)/sqrt(M_PI * GaussianParamB);
    } else if (MultiplicityParametrisation==3) {
      /* 
     * UserDefined multiplicity parametrisation goes here 
     */
    }
  }
  else if(MEModeling == 1) {
     if(sqrt_hats > 2895.5) { return 0.; } // return 0. if the ME is larger than 2895.5 GeV, the maximum value in the interpolator

    //get the hadrons:
    hadron1 = dynamic_ptr_cast<tcBeamPtr>(lastParticles().first->dataPtr());
    hadron2 = dynamic_ptr_cast<tcBeamPtr>(lastParticles().second->dataPtr());
    x1 = lastX1();
    x2 = lastX2();
    
    //cout << "x1, x2= " << x1 << ", " << x2 << endl;
    //cout << "hadron1, hadron2 = " << hadron1->id() << " " << hadron2->id() << endl;
    tcPDPtr gluon = getParticleData(ParticleID::g);
    /* reweigh the ME2 to change the factorization scale:
     * first get the PDF weights used by default,
     * divide by those and multiply by the ones for the new factorisation scale
     */
    double gPDF1_orig   = hadron1->pdf()->xfx(hadron1,gluon,sHat(),x1)/x1;
    double gPDF2_orig   = hadron2->pdf()->xfx(hadron2,gluon,sHat(),x2)/x2;
    //cout << "gPDF1_orig, gPDF2_orig = " << gPDF1_orig << " " << gPDF2_orig << endl;
    double gPDF1_rw   = hadron1->pdf()->xfx(hadron1,gluon,FactorizationScale(),x1)/x1;
    double gPDF2_rw   = hadron2->pdf()->xfx(hadron2,gluon,FactorizationScale(),x2)/x2;
    //cout << "gPDF1_rw, gPDF2_rw = " << gPDF1_rw << " " << gPDF2_rw << endl;
    //cout << "RW factor = " << gPDF1_rw * gPDF2_rw / gPDF2_orig / gPDF1_orig << endl;
    
    //reweigh the ME to change the factorization scale:
    mesq *= gPDF1_rw * gPDF2_rw / gPDF2_orig / gPDF1_orig;

    /**
     * write out the hat-s
     */
    /*    ofstream hatsfile;
    hatsfile.open ("hats.txt", ios::app);
    hatsfile << sqrt_hats << endl;*/

    // get the mean number of gluons from the interpolation
    interpolated_meangluons = (*interpol_meangluons)(sqrt_hats);

    //Poisson norm:
    /*double poisson_norm = 0;
    for(int pp = 0; pp < ngluon_max; pp++) {
      poisson_norm += pow(interpolated_meangluons, pp) * exp(-interpolated_meangluons)/factorial(pp);
      cout << "pow(interpolated_meangluons, pp) * exp(-interpolated_meangluons)/factorial(pp)=" << pow(interpolated_meangluons, pp) * exp(-interpolated_meangluons)/factorial(pp) << endl;
    }
    cout << "poisson_norm = " << poisson_norm << endl;*/
    
    
    // multiply the ME with the appropriate Poisson factor
    double poissonfac = pow(interpolated_meangluons, ngluon) * exp(-interpolated_meangluons)/factorial(ngluon);
    mesq *= poissonfac;
    // get the sigma_hat from the interpolation:
    interpolated_sigmahat = (*interpol_sigmahat)(sqrt_hats)*2.568E-9; //convert pb to GeV^-2 [1 pb = 2.567E-9 GeV^-2]
    mesq *= interpolated_sigmahat;
    //divide by phase space volume:
    //double psvolume = pow(Constants::pi/2, noutgoing-1) * pow(sHat()/GeV/GeV, noutgoing-2) / factorial(noutgoing-1) / factorial(noutgoing-2);
    //mesq /= psvolume;

    //reset jacobian to 1 (remove phase-space weights):
    mesq /= jacobian();

    //from Sherpa (???)
    mesq *= 2.*sHat()/GeV/GeV;

    if(quarkpair_option == 1) mesq *= 0.5; //divide by two if there's 4+5 flavours "active"
    
    //cout << sqrt_hats << "\t" << mesq << endl;
    //cout << "sqrt_hats, interpolated_meangluons, interpolated_sigmahat, mesq, noutgoing = " << sqrt_hats << "\t" << interpolated_meangluons << "\t" << interpolated_sigmahat << "\t" << mesq << "\t" << noutgoing << endl;
    //cout << "sqrt_hats, interpolated_meangluons, interpolated_sigmahat, poissonfac, psvolume, jacobian, mesq, noutgoing = " << sqrt_hats << "\t" << interpolated_meangluons << "\t" << interpolated_sigmahat << "\t" << poissonfac << "\t" << psvolume << "\t" << jacobian() << "\t" << mesq << "\t" << noutgoing << endl;
  }
  
  
  
  return mesq; 
}

void MEInstanton::doinit() {
  //the number of maximum gluons is given by the nAdditional() number of extra partons 
  ngluonmax(this->nAdditional());

  
  setup_interpolator();

}

void MEInstanton::doinitrun() {
  //the number of maximum gluons is given by the nAdditional() number of extra partons 
  ngluonmax(this->nAdditional());
  cout << "MEModeling is set to " << MEModeling << endl;
  cout << "The factorization scale choice is set to " << facscale_option << endl;
}

multimap<tcPDPair,tcPDVector> MEInstanton::processes() const {
  //the processmap to return 
  multimap<tcPDPair,tcPDVector> processmap;
  
  //define the particles to be used in this process
  tcPDPtr g = getParticleData(ParticleID::g);
  vector<tcPDPtr> q, qb;
  for(int i = 1; i <= 5; ++i) {
    tcPDPtr quark = getParticleData(i);
    q.push_back(quark);
    qb.push_back(quark->CC());
  }

  //make the incoming partons
  tcPDPair incoming = make_pair(g, g);



  //fixed number of quark pairs:
  if(quarkpair_option == 0) {
    //the vectors that hold the outgoing partons
    tcPDVector outgoing;
    //push all the quarks and anti-quarks (one of each for the instanton processes). 
    for(int i = 0; i < nQuarkPair(); ++i) {
      outgoing.push_back(q[i]);
      outgoing.push_back(qb[i]);
    }
    // make the process using the incoming and outgoing particles so far (i.e. no additional gluons)
    processmap.insert(make_pair(incoming,outgoing));
    
    //insert (1 to ngluon_max) gluons into the process
    for(unsigned int jj=0; jj < ngluon_max; jj++) { 
      outgoing.push_back(g);
      processmap.insert(make_pair(incoming,outgoing));
    }
  } else if(quarkpair_option == 1) { //variable number of quark pairs
    //push all the quarks and anti-quarks (one of each for the instanton processes).
    for(int maxq = 4; maxq < 6; maxq++) {
      tcPDVector outgoing;	  
      for(int i = 0; i < maxq; ++i) {
	outgoing.push_back(q[i]);
	outgoing.push_back(qb[i]);
      }
      // make the process using the incoming and outgoing particles so far (i.e. no additional gluons)
      processmap.insert(make_pair(incoming,outgoing));
      
      //insert (1 to ngluon_max) gluons into the process
      for(unsigned int jj=0; jj < ngluon_max; jj++) { 
	outgoing.push_back(g);
	processmap.insert(make_pair(incoming,outgoing));
      }
    }
  }

  return processmap;
}

size_t MEInstanton::GetnQuarkPair() const {
  //cout << "theNQuarkPair = " << theNQuarkPair << endl;
  if(quarkpair_option==0) { return theNQuarkPair; }
  else if(quarkpair_option==1) {
    //count the number of quark pairs
    size_t counted_quarks(0);
    for(int pp = 0; pp < mePartonData().size(); pp++) {
      if(fabs(mePartonData()[pp]->id()) < 7 && fabs(mePartonData()[pp]->id()) > 0) counted_quarks++; 
    }
    return counted_quarks/2;

  }
  return theNQuarkPair;
}

list<BlobMEBase::ColourConnection> MEInstanton::colourConnections() const {
  list<BlobMEBase::ColourConnection> res;
  
  //count the number of gluons in the given event 
  int ngluon = (meMomenta().size()-2-GetnQuarkPair()*2);

  //cout << "number of additional gluons = " << ngluon << endl;
  
  if(theColourConnections==0) {
    /*  a simple choice:qqbar pairs apart from the last one in the case of odd number of gluons.
        Then the last pair is connected to the first gluon. The rest of the gluons are 
        paired together. 
    */
    //make the colour connections for the incoming gluons
    //initial-state gluons are connected to each other (i.e. colour singlet)
    BlobMEBase::ColourConnection first; BlobMEBase::ColourConnection second;
    first.addColour(0); first.addAntiColour(1);
    second.addColour(1); second.addAntiColour(0);
    res.push_back(first);
    res.push_back(second);
  

    //make the colour connections for some of the quark lines
    vector<BlobMEBase::ColourConnection> quark_lines;
    quark_lines.resize(4);
    for(unsigned int cc = 0; cc < GetnQuarkPair()-1; cc++) {
      quark_lines[cc].addColour(2*cc+2);
      quark_lines[cc].addAntiColour(2*cc+3);
      res.push_back(quark_lines[cc]);
    }
    //if even number of gluons then connect the remaining quarks as well
    //and connect the gluons in pairs
    if(ngluon%2 == 0) {
      BlobMEBase::ColourConnection quark_lastpair;
      quark_lastpair.addColour(2*(GetnQuarkPair()-1)+2);
      quark_lastpair.addAntiColour(2*(GetnQuarkPair()-1)+3);
      res.push_back(quark_lastpair);
      //loop over gluons and pair adjacent ones 
      for(unsigned int gg = 1; gg <= ngluon/2; gg++) {
        BlobMEBase::ColourConnection firstgl;
        BlobMEBase::ColourConnection secondgl;
        firstgl.addColour(2*GetnQuarkPair()+1+2*gg-1);
        firstgl.addAntiColour(2*GetnQuarkPair()+2+2*gg-1);
        secondgl.addColour(2*GetnQuarkPair()+2+2*gg-1);
        secondgl.addAntiColour(2*GetnQuarkPair()+1+2*gg-1);
        res.push_back(firstgl);
        res.push_back(secondgl);      
      }  
    } else {
      //if odd number of gluons, connect the last quark pair to one of them
      // and the rest in pairs as before. 
      BlobMEBase::ColourConnection quark_lastpair1;
      BlobMEBase::ColourConnection quark_lastpair2;
    
      quark_lastpair1.addColour(2*(GetnQuarkPair()-1)+2);
      quark_lastpair1.addAntiColour(2*(GetnQuarkPair()-1)+2+2);
      res.push_back(quark_lastpair1);
      quark_lastpair2.addAntiColour(2*(GetnQuarkPair()-1)+3);
      quark_lastpair2.addColour(2*(GetnQuarkPair()-1)+2+2);
      res.push_back(quark_lastpair2);
      for(unsigned int gg = 1; gg <= (ngluon-1)/2; gg++) {
        BlobMEBase::ColourConnection firstg;
        BlobMEBase::ColourConnection secondg;
        firstg.addColour(GetnQuarkPair()*2+2+2*gg-1);
        firstg.addAntiColour(GetnQuarkPair()*2+3+2*gg-1);
        secondg.addColour(GetnQuarkPair()*2+3+2*gg-1);
        secondg.addAntiColour(GetnQuarkPair()*2+2+2*gg-1);
        res.push_back(firstg);
        res.push_back(secondg);      
      }
    } 
  } else if(theColourConnections==1) { //purely random selection (colour singlet gg only)
    //construct the array of colour and anticolour numbers
    vector<int> colours; vector<int> anticolours; vector<int> colourmap;
    // cout << "number of particles = " << meMomenta().size() << " number of gluons = " << ngluon << endl;
    //loop over the qqbar
    for(unsigned int cc = 0; cc < GetnQuarkPair(); cc++) {
      colours.push_back(2*cc+2);
      anticolours.push_back(2*cc+3);
    }
    //loop over the gluons
    for(int gg = 0; gg < ngluon; gg++) {
      colours.push_back(2+GetnQuarkPair()*2+gg);
      anticolours.push_back(2+GetnQuarkPair()*2+gg);
    }
    //pick a random element of the anticolour array that has not already been used. 
    int col = UseRandom::rnd(0, int(colours.size()));
    for(int pp = 0; pp < colours.size(); pp++) {
      col = UseRandom::rnd(0, int(anticolours.size()));
      if(anticolours.size()>1) { 
        while(colours[pp]==anticolours[col]) {
          col = UseRandom::rnd(0, int(anticolours.size()));
          /*cout << "anticolours.size() = " << anticolours.size() << endl;
          cout << "col chosen = " << col << endl;
          cout << "trying to connect " << colours[pp] << " to -" << anticolours[col] << endl;*/
        }
        colourmap.push_back(anticolours[col]);
        anticolours.erase(anticolours.begin()+col);
        continue; 
      } else if(anticolours.size()==1) { 
        if(colours[pp]==anticolours[col]) {
          // cout << "only one colour/anti-colour left and they belong to the same gluon!" << endl;
          //pick a random element of the colourmap 
          int switchcol = UseRandom::rnd(0, int(colourmap.size())-1);
          //save the previous anticolour that was placed there
          int acolourold = colourmap[switchcol];
          /*cout << "switching a randomly-chosen colour, element " << switchcol << " corresponding to colour " <<  colours[switchcol] << " to connect to anticolour " << anticolours[col] << endl;
            cout << "this used to correspond to anticolour " << acolourold << endl;*/
          //push back the new connection 
          colourmap.push_back(acolourold);
          //and change the old one to the last anticolour
          colourmap[switchcol] = anticolours[col];
          continue;
        } else {
          colourmap.push_back(anticolours[col]);
          continue; 
        }
      }
    }
    for(int ii = 0; ii < colourmap.size(); ii++) {
      BlobMEBase::ColourConnection conline;
      conline.addColour(colours[ii]);
      conline.addAntiColour(colourmap[ii]);
      //  cout << "connecting " << colours[ii] << " to -" << colourmap[ii] << endl;
      res.push_back(conline);
    }
    
    //make the colour connections for the incoming gluons
    BlobMEBase::ColourConnection first; BlobMEBase::ColourConnection second;
    first.addColour(0); first.addAntiColour(1);
    second.addColour(1); second.addAntiColour(0);
    res.push_back(first);
    res.push_back(second);
    // cout << "done connecting randomly" << endl;
  } else if(theColourConnections==2) { //purely random selection (singlet + octet gg)
    //construct the array of colour and anticolour numbers
    vector<int> colours; vector<int> anticolours; vector<int> colourmap;
    //   cout << "number of particles = " << meMomenta().size() << " number of gluons = " << ngluon << endl;
    //loop over the qqbar
    colours.push_back(0);
    anticolours.push_back(1);
    for(unsigned int cc = 0; cc < GetnQuarkPair(); cc++) {
      colours.push_back(2*cc+2);
      anticolours.push_back(2*cc+3);
    }
    //loop over the gluons
    for(int gg = 0; gg < ngluon; gg++) {
      colours.push_back(2+GetnQuarkPair()*2+gg);
      anticolours.push_back(2+GetnQuarkPair()*2+gg);
    }
    //pick a random element of the anticolour array that has not already been used. 
    int col = UseRandom::rnd(0, int(colours.size()));
    for(int pp = 0; pp < colours.size(); pp++) {
      col = UseRandom::rnd(0, int(anticolours.size()));
      if(anticolours.size()>1) { 
        while(colours[pp]==anticolours[col]) {
          col = UseRandom::rnd(0, int(anticolours.size()));
          /*cout << "anticolours.size() = " << anticolours.size() << endl;
          cout << "col chosen = " << col << endl;
          cout << "trying to connect " << colours[pp] << " to -" << anticolours[col] << endl;*/
        }
        colourmap.push_back(anticolours[col]);
        anticolours.erase(anticolours.begin()+col);
        continue; 
      } else if(anticolours.size()==1) { 
        if(colours[pp]==anticolours[col]) {
          // cout << "only one colour/anti-colour left and they belong to the same gluon!" << endl;
          //pick a random element of the colourmap 
          int switchcol = UseRandom::rnd(0, int(colourmap.size())-1);
          //save the previous anticolour that was placed there
          int acolourold = colourmap[switchcol];
          /*cout << "switching a randomly-chosen colour, element " << switchcol << " corresponding to colour " <<  colours[switchcol] << " to connect to anticolour " << anticolours[col] << endl;
            cout << "this used to correspond to anticolour " << acolourold << endl;*/
          //push back the new connection 
          colourmap.push_back(acolourold);
          //and change the old one to the last anticolour
          colourmap[switchcol] = anticolours[col];
          continue;
        } else {
          colourmap.push_back(anticolours[col]);
          continue; 
        }
      }
    }
    for(int ii = 0; ii < colourmap.size(); ii++) {
      BlobMEBase::ColourConnection conline;
      if(colourmap[ii] == 1 && colours[ii]!= 0) {
        conline.addColour(colourmap[ii]);
        conline.addColour(colours[ii]);
	//        cout << "connecting " << colours[ii] << " to " << colourmap[ii] << endl;

      }
      if(colours[ii] == 0 && colourmap[ii]!=1) {
        conline.addAntiColour(colours[ii]);
        conline.addAntiColour(colourmap[ii]);
        //cout << "connecting -" << colours[ii] << " to -" << colourmap[ii] << endl;
                
      }
      if (colourmap[ii] != 1 && colours[ii] != 0){ 
        conline.addColour(colours[ii]);
        conline.addAntiColour(colourmap[ii]);
        //cout << "connecting " << colours[ii] << " to -" << colourmap[ii] << endl;

      }
      if (colourmap[ii] == 1 && colours[ii] == 0){ 
        conline.addAntiColour(colours[ii]);
        conline.addColour(colourmap[ii]);
	//        cout << "connecting -" << colours[ii] << " to " << colourmap[ii] << endl;
        
      }
      res.push_back(conline);
    }
    
    //make the colour connections for the incoming gluons
    BlobMEBase::ColourConnection first; 
    first.addColour(0); first.addAntiColour(1);
    //    cout << "connecting 0 to -1" << endl;

    res.push_back(first);
    // cout << "done connecting randomly" << endl;
  }

  return res;
}

size_t MEInstanton::nOutgoing() const {
  return (2*nQuarkPair());
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MEInstanton::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << theNQuarkPair << ngluon_max << MultiplicityParametrisation << MEModeling << GaussianParamA << GaussianParamB << PoissonMean << theColourConnections << interpol_invrho << interpol_alphasrho << interpol_meangluons << interpol_sigmahat << facscale_option << quarkpair_option;
}

void MEInstanton::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> theNQuarkPair >> ngluon_max >> MultiplicityParametrisation >> MEModeling >> GaussianParamA >> GaussianParamB >> PoissonMean >> theColourConnections >> interpol_invrho >> interpol_alphasrho >> interpol_meangluons >> interpol_sigmahat >> facscale_option >> quarkpair_option;
}
    


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MEInstanton,Herwig::BlobME>
  describeHerwigMEInstanton("Herwig::MEInstanton", "Instantons.so");

void MEInstanton::Init() {
  
  static ClassDocumentation<MEInstanton> documentation
    ("There is no documentation for the MEInstanton class");

  static Parameter<MEInstanton,size_t> interfaceNQuarkPair
    ("NQuarkPair",
     "The number of quark pairs to consider.",
     &MEInstanton::theNQuarkPair, 4, 1, 6,
     false, false, Interface::limited);

  static Switch<MEInstanton,unsigned int> interfaceColourConnections
    ("ColourConnections",
     "How to connect the colour lines",
     &MEInstanton::theColourConnections, 0, false, false);
  static SwitchOption interfaceColourConnectionsSimple
    (interfaceColourConnections,
     "Simple",
     "A very simple arbitrary choice.",
     0);
  static SwitchOption interfaceColourConnectionsRandom
    (interfaceColourConnections,
     "Random",
     "Completely randomized colour connections. Singlet gg only.",
     1);
    static SwitchOption interfaceColourConnectionsRandom2
    (interfaceColourConnections,
     "Random2",
     "Completely randomized colour connections. Singlet gg and gg with final state.",
     2);


    static Switch<MEInstanton,unsigned int> interfaceFactorizationScale
    ("FactorizationScale",
     "The choice of factorization scale if KKS modeling is chosen",
     &MEInstanton::facscale_option, 0, false, false);
  static SwitchOption interfaceFactorizationScaleInvRho
    (interfaceFactorizationScale,
     "InvRho",
     "Use InvRho**2 as the factorization scale",
     0);
  static SwitchOption interfaceFactorizationScaleInvsHat
    (interfaceFactorizationScale,
     "sHat",
     "Use sHat() as the factorization scale",
     1);

  static Switch<MEInstanton,unsigned int> interfaceQuarkPairs
    ("QuarkPairs",
     "Whether to fix or vary the quark pairs",
     &MEInstanton::quarkpair_option, 0, false, false);
  static SwitchOption interfaceNQuarkPairsFixed
    (interfaceQuarkPairs,
     "Fixed",
     "Fix the quark pairs to the NQuarkPair number",
     0);
  static SwitchOption interfaceNQuarkPairsVariable
    (interfaceQuarkPairs,
     "Variable",
     "Vary the quark pairs (4 or 5)",
     1);
 
    static Switch<MEInstanton,unsigned int> interfaceMultiplicityParametrisation
    ("MultiplicityParametrisation",
     "How to weigh the different gluon multiplicities",
     &MEInstanton::MultiplicityParametrisation, 0, false, false);
  static SwitchOption interfaceMultiplicityParametrisationPoisson
    (interfaceMultiplicityParametrisation,
     "Poisson",
     "The multiplicity is parametrised as a Poisson distribution with mean PoissonMean.",
     0);
  static SwitchOption interfaceMultiplicityParametrisationGaussian
    (interfaceMultiplicityParametrisation,
     "Gaussian",
     "The multiplicity is parametrised as a Gaussian distribution with parameters GaussianParamA (the mean) and GaussianParamB (the denominator in the exponent).",
     1);
  static SwitchOption interfaceMultiplicityParametrisationFlat
    (interfaceMultiplicityParametrisation,
     "Flat",
     "The multiplicity is kept flat (no preference).",
     2);
    static SwitchOption interfaceMultiplicityParametrisationUserDefined
    (interfaceMultiplicityParametrisation,
     "UserDefined",
     "The multiplicity is parametrised via a user-defined function. If no function is provided, then this is identical to flat.",
     3);

    static Switch<MEInstanton,unsigned int> interfaceMEModeling
    ("MEModeling",
     "How to model the matrix element",
     &MEInstanton::MEModeling, 0, false, false);
  static SwitchOption interfaceMEModelingPureMultiplicity
    (interfaceMEModeling,
     "PureMultiplicity",
     "Flat ME with MultiplicityParametrisation giving the kind of distribution for the final state gluons.",
     0);
  static SwitchOption interfaceMEModelingKKS
    (interfaceMEModeling,
     "KKS",
     "Matrix element according to Khoze, Krauss, Schott (1911.09726).",
     1);

    static Parameter<MEInstanton, double> interfaceGaussianParamA
    ("GaussianParamA",
     "GaussianParamA",
     &MEInstanton::GaussianParamA,5. ,-1.E99, 1.E99,
     false, false, Interface::limited);

    static Parameter<MEInstanton, double> interfaceGaussianParamB
    ("GaussianParamB",
     "GaussianParamB",
     &MEInstanton::GaussianParamB,200. ,-1.E99, 1.E99,
     false, false, Interface::limited);
    
    static Parameter<MEInstanton, double> interfacePoissonMean
      ("PoissonMean",
     "PoissonMean",
     &MEInstanton::PoissonMean,3. ,-1.E99, 1.E99,
     false, false, Interface::limited);
    
    
}

