
#include "StandardModel_two_scale_model.hpp"
#include "ew_input.hpp"
#include "logger.hpp"
#include "wrappers.hpp"

using namespace flexiblesusy;

void setup(StandardModel<Two_scale>& mssm)
{
   Eigen::Matrix<double,3,3> Yu;
   Eigen::Matrix<double,3,3> Yd;
   Eigen::Matrix<double,3,3> Ye;
   double Mu;
   double g1;
   double g2;
   double g3;
   double v;
   double lam;
   // susy parameters
   // these are madeup values for now
   Yu << 1.26136e-05, 0, 0,
                   0, 0.00667469, 0,
                   0, 0, 0.857849;

   Yd << 0.000242026, 0, 0,
                   0, 0.00529911, 0,
                   0, 0, 0.193602;

   Ye << 2.84161e-05, 0, 0,
                   0, 0.00587557, 0,
                   0, 0, 0.10199;

   Mu = 627.164;
   g1 = 0.468171;
   g2 = 0.642353;
   g3 = 1.06459;
   v = 247.2;
   //   lam what value?
   lam = 0.5;
   // set parameters
   mssm.set_scale(Electroweak_constants::MZ);
   mssm.set_Yu(Yu);
   mssm.set_Yd(Yd);
   mssm.set_Ye(Ye);
   //this appears in soft pars as dim2 parameter
   //check this carefully
   mssm.set_mu2(Mu*Mu);
   mssm.set_g1(g1);
   mssm.set_g2(g2);
   mssm.set_g3(g3);
   mssm.set_v(v);
   mssm.set_Lambdax(lam);
}

double rho2(double r) {/// checked 
  if (r <= 1.9)
    return 19.0 - 16.5 * r + 43.0 * sqr(r) / 12.0 + 7.0 / 120.0 * sqr(r) * r -
      PI * sqrt(r) * (4.0 - 1.5 * r + 3.0 / 32.0 * sqr(r) + sqr(r) * r /
		      256.0) - sqr(PI) * (2.0 - 2.0 * r + 0.5 * sqr(r)) -
      log(r) * (3.0 * r - 0.5 * sqr(r));
  else {
    double rm1 = 1.0 / r, rm2 = sqr(rm1), rm3 = rm2 * rm1, rm4 = rm3 * rm1,
      rm5 = rm4 * rm1;
    return sqr(log(r)) * (1.5 - 9.0 * rm1 - 15.0 * rm2 - 48.0 * rm3 - 168.0
			  * rm4 - 612.0 * rm5) -
      log(r) * (13.5 + 4.0 * rm1 - 125.0 / 4.0 * rm2 - 558.0 / 5.0 * rm3 -
		8307.0 / 20.0 * rm4 - 109321.0 / 70.0 * rm5)
      + sqr(PI) * (1.0 - 4.0 * rm1 - 5.0 * rm2 - 16.0 * rm3 -
		   56.0 * rm4 - 204.0 * rm5)
      + 49.0 / 4.0 + 2.0 / 3.0 * rm1 + 1613.0 / 48.0 * rm2 + 87.57 * rm3 +
      341959.0 / 1200.0 * rm4 + 9737663.0 / 9800.0 * rm5;
  }
}

///PA:  I don't think we need any of this but keeping it bhere for now
///     just in case.

/*  double fEff(double x) {
  double arg = 1.0 / (1.0 + x);
  return 2.0 / x + 3.5 - (3.0 + 2.0 / x) * log(x) +
    sqr(1.0 + 1.0 / x) * 
    (2.0 * dilog(arg) - sqr(PI) / 3.0 + sqr(log(1.0 + x)));
}

double gEff(double x) {
  double y = sqrt(x / (4.0 - x));
  return (1.0 / x + 0.5) * (atan(y) / y - 1.0) + 9.0 / 8.0 + 0.5 / x -
    (1.0 + 0.5 / x) * 4.0 / x * sqr(atan(y));
}


/// calculates sinsqThetaEff as given in BPMZ
/// should make sure drbar/msbar spectrum has been calculated before 
/// brfore we call this. 

double sinSqThetaEff(StandardModel<Two_scale> sm) {
  /// PA: get renormalisation scale
  double Q = sm.get_scale();
  double mz = sm.get_physical().MVZ;
  double mw = sm.get_physical().MVWm;
  /// PA MZ is the pole mass input, replace with EW::MZ for FS, but here
  /// we want recalculated value (I think).
  if (Q != mz) {
    throw("Should call Softsusy<SoftPars>::sinSqThetaEff() at MZ only\n");
  }
  double kl = 0.;
 
  //  calcDrBarPars();
  double alphaMsbar = dataSet.displayAlpha(ALPHA);
  double alphaDrbar = qedSusythresh(alphaMsbar, Q);
  double sinthDrbar = calcSinthdrbar();
  double costhDrbar = sqrt(1.0 - sqr(sinthDrbar));
  double costhW     = mw / mz;

  double vLmzSq = 0.5 * fEff(1.0 / sqr(costhW)) + 
    4.0 * sqr(costhDrbar) * gEff(1.0 / sqr(costhW)) -
    (1.0 - 6.0 * sqr(sinthDrbar) + 8.0 * sqr(sinthDrbar) * sqr(sinthDrbar)) /
    (4.0 * sqr(costhDrbar)) * fEff(1.0);

  kl = 1.0 + costhDrbar / sinthDrbar * 
    (piZGT(mz, Q) - piZGT(0., Q)) / sqr(mz) +
    alphaDrbar / PI * sqr(costhDrbar) / sqr(sinthDrbar) * 2.0 * log(costhW)
    - alphaDrbar / (4.0 * PI * sqr(sinthDrbar)) * vLmzSq;

  return sqr(sinthDrbar) * kl;
}
*/

/// outrho, outsin represent the DRbar values
double deltaVb(double outrho, double outsin, double alphaDRbar, 
	       StandardModel<Two_scale> sm)  {

  double mz = sm.get_physical().MVZ;
  double mw = sm.get_physical().MVWp;
  double costh   = (mw / mz);
  double cw2     = sqr(costh) ;
  double sw2     = (1.0 - cw2);
  double outcos  = sqrt(1.0 - sqr(outsin));

  double deltaVbSm = outrho * alphaDRbar / (4.0 * PI * sqr(outsin)) *
    (6.0 + log(cw2) / sw2 * 
     (3.5 - 2.5 * sw2 - sqr(outsin) * (5.0 - 1.5 * cw2 / sqr(outcos))));
  

  return deltaVbSm;
}


double dRho(double outrho, double outgmu, double outsin, double alphaDRbar,  
	    StandardModel<Two_scale> sm, double pizztMZ, double piwwtMW) {
  //PA: must decide what to do with mz here
  // double mz = displayMz();
  double mz = sm.get_physical().MVZ;
  double mw = sm.get_physical().MVWp;

  /// 2 loop SM contribution
  double mt   = sm.get_physical().MFu(3); 
  
  double xt = 3.0 * outgmu * sqr(mt) / (8.0 * sqr(PI) * root2);
  
  double deltaRho2LoopSm = alphaDRbar * sqr(sm.get_g3()) / 
    (16.0 * PI * sqr(PI) * sqr(outsin)) * /// bug-fixed 24.08.2002
    (-2.145 * sqr(mt) / sqr(mw) + 1.262 * log(mt / mz) - 2.24 
     - 0.85 * sqr(mz)
     / sqr(mt)) + sqr(xt) * rho2(sm.get_Mhh() / mt) / 3.0;

  double deltaRhoOneLoop = pizztMZ / (outrho * sqr(mz))
    - piwwtMW / sqr(mw);
  
  double deltaRho = deltaRhoOneLoop + deltaRho2LoopSm;

  return deltaRho;
}


double dR(double outrho, double outgmu, double outsin, double alphaDRbar,
	  StandardModel<Two_scale> sm, double pizztMZ, double piwwt0) {
 
  double outcos = cos(asin(outsin));
  /// 2 loop SM contribution
  double mt   = sm.get_physical().MFu(3);
  
  double xt = 3.0 * outgmu * sqr(mt) / (8.0 * sqr(PI) * root2);
  
  double dvb = deltaVb(outrho, outsin, alphaDRbar, sm);

  //  double mz = displayMz(); // MZ
  double mz = sm.get_physical().MVZ;
  double mw = sm.get_physical().MVWp;
  double deltaR =  outrho * piwwt0 / sqr(mw) - 
    pizztMZ / sqr(mz) + dvb;

  /// Dominant two-loop SM term
  double deltaR2LoopSm = alphaDRbar * sqr(sm.get_g3()) / 
    (16.0 * sqr(PI) * PI * sqr(outsin) * sqr(outcos)) *
    (2.145 * sqr(mt) / sqr(mz) + 0.575 * log(mt / mz) - 0.224 
     - 0.144 * sqr(mz) / sqr(mt)) - 
    sqr(xt) * 
    rho2(sm.get_Mhh() / mt) * (1.0 - deltaR) * outrho / 3.0;

  deltaR = deltaR + deltaR2LoopSm; 

  return deltaR;
}

/// Checked 20.11.00
/// Flags noconvergence if there's trouble...then don't believe outrho and
/// outsin produced - they are fudged!
/// PA: This is from softsusy and calculates rhohat in the case where one gets
/// sinthetaW from input Gmu
// void rhohatold(double & outrho, double & outsin, double alphaDRbar,
// 	    StandardModel<Two_scale> sm, double pizztMZ, double piwwt0, 
// 	    double piwwtMW, double tol, int maxTries, int & err) {

//   static double oldrho = 0.23, oldsin = 0.8;

//   //  double mz = displayMz();  //displayMz() returns MZ in the most complicated and unsafe way ever!
//   double mz = get_physical().MVZ //we need to calculate MVZ pole first this way
//   if (sm.get_scale() != mz) {
//     cerr << "rhohat c alled at scale which is not mz" << endl;
//     err = 1;
//     return;
//   }
  
//   static int numTries = 0;
  
//   if ((outrho < TOLERANCE || outsin < TOLERANCE) || fabs(outsin) > 1. ||
//       (numTries - 1 > maxTries)) {  
//     oldrho = 0.23; oldsin = 0.8;
//     numTries = 0;
//     outrho = 0.23; outsin = 0.8; 
//     flagNoRhoConvergence(true);
//     cerr  << "rhohat reached maxtries\n";
//     err = 2;
//     return;
//   }
  
//   /// Difference to last iteration
//   double sumTol;
//   sumTol = fabs(oldrho / outrho - 1.0) + fabs(oldsin / outsin - 1.0);
  
//   if (numTries != 0 && sumTol < tol) {
//     numTries = 0;
//     oldrho = 0.23, oldsin = 0.8;
//     // if (PRINTOUT > 2) cout << "sin rho converged\n";
//     err = 0;
//     return;
//   }

//   numTries = numTries + 1;
  
//   oldrho = outrho; oldsin = outsin; 
    
//   double deltaR = dR(outrho, outsin, alphaDRbar, pizztMZ, piwwt0);

//   // if (PRINTOUT > 2) cout << " st2=" << sumTol << " dr=" << deltaR 
//   // 			 << " outrho=" << outrho << " outsin=" << outsin 
//   // 			 << " aDRbar=" << alphaDRbar << " piZ=" << pizztMZ 
//   // 			 << " piwwt0=" << piwwt0 << endl;
  
//   double sin2thetasqO4 = PI * alphaDRbar / 
//     (root2 * sqr(mz) * GMU * (1.0 - deltaR)); 

//   if (sin2thetasqO4 >= 0.25) sin2thetasqO4 = 0.25;
//   if (sin2thetasqO4 < 0.0) sin2thetasqO4 = 0.0;

//   double sin2theta = sqrt(4.0 * sin2thetasqO4);

//   double theta = 0.5 * asin(sin2theta);
  
//   outsin = sin(theta); 

//   double deltaRho = dRho(outrho, outsin, alphaDRbar, pizztMZ, piwwtMW);

//   if (fabs(deltaRho) < 1.0) outrho = 1.0 / (1.0 - deltaRho);
//   else outrho = 1.0;

//   // if (PRINTOUT > 2) cout << " drho=" << deltaRho << " sw=" << outsin << endl; 

//   rhohatold(outrho, outsin, alphaDRbar, pizztMZ, piwwt0, piwwtMW, tol, maxTries);
// }



/// here we assume MZ and MW already calculated and we are calculatin this
/// to get MW.
/// Initial guess for gmu should be: outgmu = PI * alphaDRbar /(root2 * sqr(mz) * sinthW * costhW );
void rhohat(double & outrho, double & outgmu, double alphaDRbar,
	    StandardModel<Two_scale> sm, double pizztMZ, double piwwt0, 
	    double piwwtMW, double tol, int maxTries, int & err) {

  static double oldrho = 0.23, oldgmu = 1e-15;

  //  double mz = displayMz();  //displayMz() returns MZ in the most complicated and unsafe way ever!
  double mz = sm.get_physical().MVZ; //we need to calculate MVZ pole first this way
  if (sm.get_scale() != mz) {
    cerr << "rhohat c alled at scale which is not mz" << endl;
    err = 1;
    return;
  }
  
  static int numTries = 0;
  
  if ((outrho < TOLERANCE) || (numTries - 1 > maxTries)) {  
    oldrho = 0.23; 
    numTries = 0;
    outrho = 0.23; 
    cerr  << "rhohat reached maxtries\n";
    err = 2;
    return;
  }
  
  /// Difference to last iteration
  double sumTol;
  sumTol = fabs(oldrho / outrho - 1.0) + fabs(oldgmu / outgmu - 1.0);
  
  if (numTries != 0 && sumTol < tol) {
    numTries = 0;
    oldrho = 0.23, oldgmu = 0;
    // if (PRINTOUT > 2) cout << "sin rho converged\n";
    err = 0;
    return;
  }

  numTries = numTries + 1;
  
  oldrho = outrho;  oldgmu = outgmu;

  ///PA: assumes previously calculated
  double mzdrbar = sm.get_MVZ();
  double mwdrbar = sm.get_MVWp();
  //DRbar sinthW
  double sinthW = sqrt(1 - mzdrbar*mzdrbar / (mwdrbar* mwdrbar)  );  
  double costhW =  mzdrbar / mwdrbar;

  double deltaR = dR(outrho, outgmu, sinthW, alphaDRbar, sm, pizztMZ, piwwt0);      
    
  outgmu = PI * alphaDRbar / 
    (root2 * sqr(mz) * sinthW * costhW * (1.0 - deltaR));  

  
  double deltaRho = dRho(outrho, outgmu, sinthW, alphaDRbar, sm, pizztMZ, 
			 piwwtMW);

  if (fabs(deltaRho) < 1.0) outrho = 1.0 / (1.0 - deltaRho);
  else outrho = 1.0;

  // if (PRINTOUT > 2) cout << " drho=" << deltaRho << " sw=" << outsin << endl; 

  rhohat(outrho, outgmu, alphaDRbar, sm, pizztMZ, piwwt0, piwwtMW, tol, maxTries, err);
}



int main()
{
   INFO("=============================");
   INFO("running sm_example()");
   INFO("=============================");

   StandardModel<Two_scale> sm;
   //setup  sm  from  predifined parameters.
   //still to edit this to take sm values
   setup(sm);
   //get the DRbar masses immediately - just for fun.
   sm.calculate_DRbar_parameters();
   //maybe we can take msbar StandardModel parameters from literature
   
   // Alternativelo set up properly ourselves:
   //1. Input sm masses
   //2. calcualte running masses at MZ
   //3.  get msbar parameters

   //once we have have sm running Yukawas, gauge couplings and vevs
   // we do the following
   
   //1 change vev to input (?) value
   double newvev = 280;
   // run to VEV scale
   sm.run_to(newvev);
   //reset vev  (alternatively we can let mu and lambda float)
   sm.set_v(newvev);
   //do EWSB 

   //test print out.
   std::cout << "vev = " << sm.get_v() << std::endl;
   std::cout << "yu = " << sm.get_Yu() << std::endl;
   std::cout << "yd = " << sm.get_Yd() << std::endl;
   std::cout << "ye = " << sm.get_Ye() << std::endl;
   std::cout << "g1= " << sm.get_g1() << std::endl;
   std::cout << "g2 = " << sm.get_g2() << std::endl;
   std::cout << "g3 = " << sm.get_g3() << std::endl;

   std::cout << "mu2= " << sm.get_mu2() << std::endl;
   std::cout << "lam= " << sm.get_Lambdax() << std::endl;
  
   
   //think I must calculate DRbar masses first
   sm.calculate_DRbar_parameters();
   sm.print(std::cout);
   //calcualte pole masses
   
   //currently not calculating right masses
   // simple to hack or create alternative version
   sm.calculate_pole_masses();
   // is it this simple?
   std::cout << "after calculating pole mases MZ =" << sm.get_physical().MVZ << std::endl;  
   sm.print(std::cout);
   //may need problem point test at this point.

   //Now I must calculate GMU here.
    
   //do decays

   return 0;
}
