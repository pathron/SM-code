
#include "SM_two_scale_model.hpp"
#include "ew_input.hpp"
#include "logger.hpp"
#include "wrappers.hpp"

using namespace flexiblesusy;

void setup(SM<Two_scale>& mssm)
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


int main()
{
   INFO("=============================");
   INFO("running sm_example()");
   INFO("=============================");

   SM<Two_scale> sm;
   //setup  sm  from  predifined parameters.
   //still to edit this to take sm values
   setup(sm);
   //get the DRbar masses immediately - just for fun.
   sm.calculate_DRbar_parameters();
   //maybe we can take msbar SM parameters from literature
   
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
   std::cout << "after calculating pole mases MZ =" << sm.get_physical().MVZ << std::endl;  sm.print(std::cout);
   //may need problem point test at this point.

   //do decays




   return 0;
}
