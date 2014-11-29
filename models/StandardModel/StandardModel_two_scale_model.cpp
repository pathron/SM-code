// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

// File generated at Sat 29 Nov 2014 17:55:13

/**
 * @file StandardModel_two_scale_model.cpp
 * @brief implementation of the StandardModel model class
 *
 * Contains the definition of the StandardModel model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Sat 29 Nov 2014 17:55:13 with FlexibleSUSY
 * 1.0.3 and SARAH 4.2.2 .
 */

#include "StandardModel_two_scale_model.hpp"
#include "wrappers.hpp"
#include "linalg2.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "root_finder.hpp"
#include "gsl_utils.hpp"
#include "config.h"
#include "pv.hpp"



#include <cmath>
#include <iostream>

#ifdef ENABLE_THREADS
#include <thread>
#endif

namespace flexiblesusy {

using namespace StandardModel_info;

#define CLASSNAME StandardModel<Two_scale>

#define PHYSICAL(parameter) physical.parameter
#define INPUT(parameter) model->get_input().parameter
#define LOCALINPUT(parameter) input.parameter

#define HIGGS_2LOOP_CORRECTION_AT_AS     higgs_2loop_corrections.at_as
#define HIGGS_2LOOP_CORRECTION_AB_AS     higgs_2loop_corrections.ab_as
#define HIGGS_2LOOP_CORRECTION_AT_AT     higgs_2loop_corrections.at_at
#define HIGGS_2LOOP_CORRECTION_ATAU_ATAU higgs_2loop_corrections.atau_atau

#ifdef ENABLE_THREADS
   std::mutex CLASSNAME::mtx_fortran;
   #define LOCK_MUTEX() mtx_fortran.lock()
   #define UNLOCK_MUTEX() mtx_fortran.unlock()
#else
   #define LOCK_MUTEX()
   #define UNLOCK_MUTEX()
#endif

CLASSNAME::StandardModel(const StandardModel_input_parameters& input_)
   : Two_scale_model()
   , StandardModel_soft_parameters(input_)
   , number_of_ewsb_iterations(100)
   , number_of_mass_iterations(20)
   , ewsb_loop_order(2)
   , pole_mass_loop_order(2)
   , calculate_sm_pole_masses(false)
   , precision(1.0e-3)
   , ewsb_iteration_precision(1.0e-5)
   , physical()
   , problems(StandardModel_info::particle_names)
   , higgs_2loop_corrections()
#ifdef ENABLE_THREADS
   , thread_exception()
#endif
   , MVG(0), MHp(0), MFv(Eigen::Array<double,3,1>::Zero()), MAh(0), Mhh(0), MVP
      (0), MVZ(0), MFd(Eigen::Array<double,3,1>::Zero()), MFu(Eigen::Array<double,
      3,1>::Zero()), MFe(Eigen::Array<double,3,1>::Zero()), MVWp(0)

   , Vd(Eigen::Matrix<std::complex<double>,3,3>::Zero()), Ud(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), Vu(Eigen::Matrix<std::complex<double>,3,
      3>::Zero()), Uu(Eigen::Matrix<std::complex<double>,3,3>::Zero()), Ve(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), Ue(Eigen::Matrix<
      std::complex<double>,3,3>::Zero())


{
}

CLASSNAME::~StandardModel()
{
}

void CLASSNAME::do_calculate_sm_pole_masses(bool flag)
{
   calculate_sm_pole_masses = flag;
}

bool CLASSNAME::do_calculate_sm_pole_masses() const
{
   return calculate_sm_pole_masses;
}

void CLASSNAME::set_ewsb_loop_order(unsigned loop_order)
{
   ewsb_loop_order = loop_order;
}

void CLASSNAME::set_higgs_2loop_corrections(const Higgs_2loop_corrections& higgs_2loop_corrections_)
{
   higgs_2loop_corrections = higgs_2loop_corrections_;
}

void CLASSNAME::set_number_of_ewsb_iterations(std::size_t iterations)
{
   number_of_ewsb_iterations = iterations;
}

void CLASSNAME::set_number_of_mass_iterations(std::size_t iterations)
{
   number_of_mass_iterations = iterations;
}

void CLASSNAME::set_precision(double precision_)
{
   precision = precision_;
   ewsb_iteration_precision = precision_;
}

void CLASSNAME::set_pole_mass_loop_order(unsigned loop_order)
{
   pole_mass_loop_order = loop_order;
}

void CLASSNAME::set_ewsb_iteration_precision(double precision)
{
   ewsb_iteration_precision = precision;
}

double CLASSNAME::get_ewsb_iteration_precision() const
{
   return ewsb_iteration_precision;
}

double CLASSNAME::get_ewsb_loop_order() const
{
   return ewsb_loop_order;
}

const StandardModel_physical& CLASSNAME::get_physical() const
{
   return physical;
}

StandardModel_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const StandardModel_physical& physical_)
{
   physical = physical_;
}

const Problems<StandardModel_info::NUMBER_OF_PARTICLES>& CLASSNAME::get_problems() const
{
   return problems;
}

Problems<StandardModel_info::NUMBER_OF_PARTICLES>& CLASSNAME::get_problems()
{
   return problems;
}

/**
 * Method which calculates the tadpoles at loop order specified in the
 * pointer to the CLASSNAME::Ewsb_parameters struct.
 *
 * @param x GSL vector of EWSB output parameters
 * @param params pointer to CLASSNAME::Ewsb_parameters struct
 * @param f GSL vector with tadpoles
 *
 * @return GSL_EDOM if x contains Nans, GSL_SUCCESS otherwise.
 */
int CLASSNAME::tadpole_equations(const gsl_vector* x, void* params, gsl_vector* f)
{
   if (contains_nan(x, number_of_ewsb_equations)) {
      for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
         gsl_vector_set(f, i, std::numeric_limits<double>::max());
      return GSL_EDOM;
   }

   const CLASSNAME::Ewsb_parameters* ewsb_parameters
      = static_cast<CLASSNAME::Ewsb_parameters*>(params);
   StandardModel* model = ewsb_parameters->model;
   const unsigned ewsb_loop_order = ewsb_parameters->ewsb_loop_order;

   double tadpole[number_of_ewsb_equations];

   model->set_mu2(gsl_vector_get(x, 0));

   tadpole[0] = model->get_ewsb_eq_hh_1();

   if (ewsb_loop_order > 0) {
      model->calculate_DRbar_masses();
      tadpole[0] -= Re(model->tadpole_hh());

      if (ewsb_loop_order > 1) {

      }
   }

   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      gsl_vector_set(f, i, tadpole[i]);

   return GSL_SUCCESS;
}

/**
 * method which solves the EWSB conditions iteratively, trying GSL
 * root finding methods
 *       gsl_multiroot_fsolver_hybrid, gsl_multiroot_fsolver_hybrids, gsl_multiroot_fsolver_broyden
 * in that order until a solution is found.
 */
int CLASSNAME::solve_ewsb_iteratively()
{
   const gsl_multiroot_fsolver_type* solvers[] = {
      gsl_multiroot_fsolver_hybrid, gsl_multiroot_fsolver_hybrids, gsl_multiroot_fsolver_broyden
   };

   double x_init[number_of_ewsb_equations];
   ewsb_initial_guess(x_init);

#ifdef ENABLE_VERBOSE
   std::cout << "Solving EWSB equations ...\n"
      "\tInitial guess: x_init =";
   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      std::cout << " " << x_init[i];
   std::cout << '\n';
#endif

   int status;
   for (std::size_t i = 0; i < sizeof(solvers)/sizeof(*solvers); ++i) {
      VERBOSE_MSG("\tStarting EWSB iteration using solver " << i);
      status = solve_ewsb_iteratively_with(solvers[i], x_init);
      if (status == GSL_SUCCESS) {
         VERBOSE_MSG("\tSolver " << i << " finished successfully!");
         break;
      }
#ifdef ENABLE_VERBOSE
      else {
         WARNING("\tSolver " << i << " could not find a solution!"
                 " (requested precision: " << ewsb_iteration_precision << ")");
      }
#endif
   }

   if (status != GSL_SUCCESS) {
      problems.flag_no_ewsb();
#ifdef ENABLE_VERBOSE
      WARNING("\tCould not find a solution to the EWSB equations!"
              " (requested precision: " << ewsb_iteration_precision << ")");
#endif
   } else {
      problems.unflag_no_ewsb();
   }

   return status;
}

int CLASSNAME::solve_ewsb_iteratively(unsigned loop_order)
{
   // temporarily set `ewsb_loop_order' to `loop_order' and do
   // iteration
   const unsigned old_loop_order = ewsb_loop_order;
   ewsb_loop_order = loop_order;
   const int status = solve_ewsb_iteratively();
   ewsb_loop_order = old_loop_order;
   return status;
}


int CLASSNAME::solve_ewsb_tree_level()
{
   int error = 0;

   const double old_mu2 = mu2;

   mu2 = 0.5*Lambdax*Sqr(v);

   const bool is_finite = std::isfinite(mu2);

   if (!is_finite) {
      mu2 = old_mu2;
      error = 1;
   }


   return error;
}

int CLASSNAME::solve_ewsb_tree_level_via_soft_higgs_masses()
{
   int error = 0;



   return error;
}

int CLASSNAME::solve_ewsb_one_loop()
{
   return solve_ewsb_iteratively(1);
}

int CLASSNAME::solve_ewsb()
{
   VERBOSE_MSG("\tSolving EWSB at " << ewsb_loop_order << "-loop order");

   if (ewsb_loop_order == 0)
      return solve_ewsb_tree_level();

   return solve_ewsb_iteratively(ewsb_loop_order);
}

void CLASSNAME::ewsb_initial_guess(double x_init[number_of_ewsb_equations])
{
   x_init[0] = mu2;

}

int CLASSNAME::solve_ewsb_iteratively_with(const gsl_multiroot_fsolver_type* solver,
                                           const double x_init[number_of_ewsb_equations])
{
   Ewsb_parameters params = {this, ewsb_loop_order};
   Root_finder<number_of_ewsb_equations> root_finder(CLASSNAME::tadpole_equations,
                              &params,
                              number_of_ewsb_iterations,
                              ewsb_iteration_precision);
   root_finder.set_solver_type(solver);
   const int status = root_finder.find_root(x_init);

   return status;
}

void CLASSNAME::print(std::ostream& ostr) const
{
   ostr << "========================================\n"
           "StandardModel\n"
           "========================================\n";
   StandardModel_soft_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "tree-level DRbar masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MHp = " << MHp << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MAh = " << MAh << '\n';
   ostr << "Mhh = " << Mhh << '\n';
   ostr << "MVP = " << MVP << '\n';
   ostr << "MVZ = " << MVZ << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MVWp = " << MVWp << '\n';

   ostr << "----------------------------------------\n"
           "tree-level DRbar mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "Vd = " << Vd << '\n';
   ostr << "Ud = " << Ud << '\n';
   ostr << "Vu = " << Vu << '\n';
   ostr << "Uu = " << Uu << '\n';
   ostr << "Ve = " << Ve << '\n';
   ostr << "Ue = " << Ue << '\n';

   physical.print(ostr);
}

/**
 * wrapper routines for passarino Veltman functions
 */

double CLASSNAME::A0(double m) const
{
   return passarino_veltman::ReA0(m*m, Sqr(get_scale()));
}

double CLASSNAME::B0(double p, double m1, double m2) const
{
   return passarino_veltman::ReB0(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

double CLASSNAME::B1(double p, double m1, double m2) const
{
   return passarino_veltman::ReB1(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

double CLASSNAME::B00(double p, double m1, double m2) const
{
   return passarino_veltman::ReB00(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

double CLASSNAME::B22(double p, double m1, double m2) const
{
   return passarino_veltman::ReB22(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

double CLASSNAME::H0(double p, double m1, double m2) const
{
   return passarino_veltman::ReH0(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

double CLASSNAME::F0(double p, double m1, double m2) const
{
   return passarino_veltman::ReF0(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

double CLASSNAME::G0(double p, double m1, double m2) const
{
   return passarino_veltman::ReG0(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

/**
 * routine which finds the DRbar mass eigenstates and mixings.
 */
void CLASSNAME::calculate_DRbar_masses()
{

   solve_ewsb_tree_level_via_soft_higgs_masses();

   calculate_MVG();
   calculate_MVP();
   calculate_MVZ();
   calculate_MVWp();
   calculate_MHp();
   calculate_MFv();
   calculate_MAh();
   calculate_Mhh();
   calculate_MFd();
   calculate_MFu();
   calculate_MFe();


}

/**
 * Backward compatibility routine which finds the DRbar mass
 * eigenstates and mixings.
 */
void CLASSNAME::calculate_DRbar_parameters()
{
   calculate_DRbar_masses();
}

/**
 * routine which finds the pole mass eigenstates and mixings.
 */
void CLASSNAME::calculate_pole_masses()
{
#ifdef ENABLE_THREADS
   thread_exception = 0;

   std::thread thread_MHp(Thread(this, &CLASSNAME::calculate_MHp_pole));
   std::thread thread_MAh(Thread(this, &CLASSNAME::calculate_MAh_pole));
   std::thread thread_Mhh(Thread(this, &CLASSNAME::calculate_Mhh_pole));

   if (calculate_sm_pole_masses) {
      std::thread thread_MFd(Thread(this, &CLASSNAME::calculate_MFd_pole));
      std::thread thread_MFe(Thread(this, &CLASSNAME::calculate_MFe_pole));
      std::thread thread_MFu(Thread(this, &CLASSNAME::calculate_MFu_pole));
      std::thread thread_MFv(Thread(this, &CLASSNAME::calculate_MFv_pole));
      std::thread thread_MVG(Thread(this, &CLASSNAME::calculate_MVG_pole));
      std::thread thread_MVP(Thread(this, &CLASSNAME::calculate_MVP_pole));
      std::thread thread_MVWp(Thread(this, &CLASSNAME::calculate_MVWp_pole));
      std::thread thread_MVZ(Thread(this, &CLASSNAME::calculate_MVZ_pole));
      thread_MFd.join();
      thread_MFe.join();
      thread_MFu.join();
      thread_MFv.join();
      thread_MVG.join();
      thread_MVP.join();
      thread_MVWp.join();
      thread_MVZ.join();
   }

   thread_MHp.join();
   thread_MAh.join();
   thread_Mhh.join();


   if (thread_exception != 0)
      std::rethrow_exception(thread_exception);
#else
   calculate_MHp_pole();
   calculate_MAh_pole();
   calculate_Mhh_pole();

   if (calculate_sm_pole_masses) {
      calculate_MFd_pole();
      calculate_MFe_pole();
      calculate_MFu_pole();
      calculate_MFv_pole();
      calculate_MVG_pole();
      calculate_MVP_pole();
      calculate_MVWp_pole();
      calculate_MVZ_pole();
   }

#endif
}

void CLASSNAME::copy_DRbar_masses_to_pole_masses()
{
   PHYSICAL(MVG) = MVG;
   PHYSICAL(MHp) = MHp;
   PHYSICAL(MFv) = MFv;
   PHYSICAL(MAh) = MAh;
   PHYSICAL(Mhh) = Mhh;
   PHYSICAL(MVP) = MVP;
   PHYSICAL(MVZ) = MVZ;
   PHYSICAL(MFd) = MFd;
   PHYSICAL(Vd) = Vd;
   PHYSICAL(Ud) = Ud;
   PHYSICAL(MFu) = MFu;
   PHYSICAL(Vu) = Vu;
   PHYSICAL(Uu) = Uu;
   PHYSICAL(MFe) = MFe;
   PHYSICAL(Ve) = Ve;
   PHYSICAL(Ue) = Ue;
   PHYSICAL(MVWp) = MVWp;

}

/**
 * reorders DRbar masses so that golstones are placed at the index
 * specified in the model files definition of the associuated
 * gauge boson (see Z-boson definition in default particles.m file
 * in the Models directory of your SARAH distribution for example)
 */
void CLASSNAME::reorder_DRbar_masses()
{

}

/**
 * reorders pole masses so that golstones are placed at the index
 * specified in the model files definition of the associuated
 * gauge boson (see Z-boson definition in default particles.m file
 * in the Models directory of your SARAH distribution for example)
 */
void CLASSNAME::reorder_pole_masses()
{

}
/**
 * calculates spectrum for model once the DRbar parameters at
 * at low energies are known
 */
void CLASSNAME::calculate_spectrum()
{
   calculate_DRbar_masses();
   if (pole_mass_loop_order > 0)
      calculate_pole_masses();

   // move goldstone bosons to the front
   reorder_DRbar_masses();
   if (pole_mass_loop_order == 0)
      copy_DRbar_masses_to_pole_masses();
   else
      reorder_pole_masses();

   if (problems.have_problem()) {
      clear_DRbar_parameters();
      physical.clear();
   }
}

void CLASSNAME::clear_DRbar_parameters()
{
   MVG = 0.;
   MHp = 0.;
   MFv = Eigen::Matrix<double,3,1>::Zero();
   MAh = 0.;
   Mhh = 0.;
   MVP = 0.;
   MVZ = 0.;
   MFd = Eigen::Matrix<double,3,1>::Zero();
   Vd = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ud = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFu = Eigen::Matrix<double,3,1>::Zero();
   Vu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Uu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFe = Eigen::Matrix<double,3,1>::Zero();
   Ve = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ue = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MVWp = 0.;


}

void CLASSNAME::clear()
{
   StandardModel_soft_parameters::clear();
   clear_DRbar_parameters();
   physical.clear();
   problems.clear();
}

std::string CLASSNAME::name() const
{
   return "StandardModel";
}

void CLASSNAME::run_to(double scale, double eps)
{
   if (eps < 0.0)
      eps = precision;
   StandardModel_soft_parameters::run_to(scale, eps);
}



void CLASSNAME::calculate_MVG()
{
   MVG = 0;
}

void CLASSNAME::calculate_MHp()
{
   MHp = 0.25*(4*mu2 - 2*Lambdax*Sqr(v) + Sqr(g2)*Sqr(v));

   problems.flag_tachyon(StandardModel_info::Hp, MHp < 0.);

   MHp = AbsSqrt(MHp);
}

void CLASSNAME::calculate_MFv()
{
   MFv.setConstant(0);
}

void CLASSNAME::calculate_MAh()
{
   MAh = 0.25*(4*mu2 - 2*Lambdax*Sqr(v) + Sqr(v)*Sqr(g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW())));

   problems.flag_tachyon(StandardModel_info::Ah, MAh < 0.);

   MAh = AbsSqrt(MAh);
}

void CLASSNAME::calculate_Mhh()
{
   Mhh = mu2 - 1.5*Lambdax*Sqr(v);

   problems.flag_tachyon(StandardModel_info::hh, Mhh < 0.);

   Mhh = AbsSqrt(Mhh);
}

void CLASSNAME::calculate_MVP()
{
   MVP = 0;
}

void CLASSNAME::calculate_MVZ()
{
   MVZ = 0.25*Sqr(v)*Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(
      ThetaW()));

   problems.flag_tachyon(StandardModel_info::VZ, MVZ < 0.);

   MVZ = AbsSqrt(MVZ);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fd() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fd;
   mass_matrix_Fd(0,0) = -0.7071067811865475*v*Yd(0,0);
   mass_matrix_Fd(0,1) = -0.7071067811865475*v*Yd(1,0);
   mass_matrix_Fd(0,2) = -0.7071067811865475*v*Yd(2,0);
   mass_matrix_Fd(1,0) = -0.7071067811865475*v*Yd(0,1);
   mass_matrix_Fd(1,1) = -0.7071067811865475*v*Yd(1,1);
   mass_matrix_Fd(1,2) = -0.7071067811865475*v*Yd(2,1);
   mass_matrix_Fd(2,0) = -0.7071067811865475*v*Yd(0,2);
   mass_matrix_Fd(2,1) = -0.7071067811865475*v*Yd(1,2);
   mass_matrix_Fd(2,2) = -0.7071067811865475*v*Yd(2,2);

   return mass_matrix_Fd;
}

void CLASSNAME::calculate_MFd()
{
   const auto mass_matrix_Fd(get_mass_matrix_Fd());

#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fd, MFd, Vd, Ud, eigenvalue_error);
   problems.flag_bad_mass(StandardModel_info::Fd, eigenvalue_error >
      precision * Abs(MFd(0)));
#else
   fs_svd(mass_matrix_Fd, MFd, Vd, Ud);
#endif
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fu() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fu;
   mass_matrix_Fu(0,0) = 0.7071067811865475*v*Yu(0,0);
   mass_matrix_Fu(0,1) = 0.7071067811865475*v*Yu(1,0);
   mass_matrix_Fu(0,2) = 0.7071067811865475*v*Yu(2,0);
   mass_matrix_Fu(1,0) = 0.7071067811865475*v*Yu(0,1);
   mass_matrix_Fu(1,1) = 0.7071067811865475*v*Yu(1,1);
   mass_matrix_Fu(1,2) = 0.7071067811865475*v*Yu(2,1);
   mass_matrix_Fu(2,0) = 0.7071067811865475*v*Yu(0,2);
   mass_matrix_Fu(2,1) = 0.7071067811865475*v*Yu(1,2);
   mass_matrix_Fu(2,2) = 0.7071067811865475*v*Yu(2,2);

   return mass_matrix_Fu;
}

void CLASSNAME::calculate_MFu()
{
   const auto mass_matrix_Fu(get_mass_matrix_Fu());

#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fu, MFu, Vu, Uu, eigenvalue_error);
   problems.flag_bad_mass(StandardModel_info::Fu, eigenvalue_error >
      precision * Abs(MFu(0)));
#else
   fs_svd(mass_matrix_Fu, MFu, Vu, Uu);
#endif
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fe() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fe;
   mass_matrix_Fe(0,0) = -0.7071067811865475*v*Ye(0,0);
   mass_matrix_Fe(0,1) = -0.7071067811865475*v*Ye(1,0);
   mass_matrix_Fe(0,2) = -0.7071067811865475*v*Ye(2,0);
   mass_matrix_Fe(1,0) = -0.7071067811865475*v*Ye(0,1);
   mass_matrix_Fe(1,1) = -0.7071067811865475*v*Ye(1,1);
   mass_matrix_Fe(1,2) = -0.7071067811865475*v*Ye(2,1);
   mass_matrix_Fe(2,0) = -0.7071067811865475*v*Ye(0,2);
   mass_matrix_Fe(2,1) = -0.7071067811865475*v*Ye(1,2);
   mass_matrix_Fe(2,2) = -0.7071067811865475*v*Ye(2,2);

   return mass_matrix_Fe;
}

void CLASSNAME::calculate_MFe()
{
   const auto mass_matrix_Fe(get_mass_matrix_Fe());

#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fe, MFe, Ve, Ue, eigenvalue_error);
   problems.flag_bad_mass(StandardModel_info::Fe, eigenvalue_error >
      precision * Abs(MFe(0)));
#else
   fs_svd(mass_matrix_Fe, MFe, Ve, Ue);
#endif
}

void CLASSNAME::calculate_MVWp()
{
   MVWp = 0.25*Sqr(g2)*Sqr(v);

   problems.flag_tachyon(StandardModel_info::VWp, MVWp < 0.);

   MVWp = AbsSqrt(MVWp);
}


double CLASSNAME::get_ewsb_eq_hh_1() const
{
   double result = mu2*v - 0.5*Power(v,3)*Lambdax;

   return result;
}



double CLASSNAME::CpconjHpHphh() const
{
   double result = 0.0;

   result = v*Lambdax;

   return result;
}

double CLASSNAME::CpconjHpVWpVP() const
{
   double result = 0.0;

   result = 0.3872983346207417*g1*g2*v*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpconjHpVZVWp() const
{
   double result = 0.0;

   result = -0.3872983346207417*g1*g2*v*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpHpgWpCbargZ() const
{
   double result = 0.0;

   result = 0.25*g2*v*(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpconjHpbargWpCgZ() const
{
   double result = 0.0;

   result = 0.05*g2*v*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW())
      );

   return result;
}

double CLASSNAME::CpHpgZbargWp() const
{
   double result = 0.0;

   result = 0.05*g2*v*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW())
      );

   return result;
}

double CLASSNAME::CpconjHpbargZgWp() const
{
   double result = 0.0;

   result = 0.25*g2*v*(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpHpconjHpAhAh() const
{
   double result = 0.0;

   result = Lambdax;

   return result;
}

double CLASSNAME::CpHpconjHphhhh() const
{
   double result = 0.0;

   result = Lambdax;

   return result;
}

std::complex<double> CLASSNAME::CpHpconjHpVZVZ() const
{
   std::complex<double> result;

   result = 0.1*(-7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW()) + 5*Sqr(
      g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW())));

   return result;
}

double CLASSNAME::CpHpconjHpconjHpHp() const
{
   double result = 0.0;

   result = 2*Lambdax;

   return result;
}

double CLASSNAME::CpHpconjHpconjVWpVWp() const
{
   double result = 0.0;

   result = 0.5*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpconjHpVWpAh() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2;

   return result;
}

double CLASSNAME::CpconjHpVWphh() const
{
   double result = 0.0;

   result = -0.5*g2;

   return result;
}

double CLASSNAME::CpconjHpVPHp() const
{
   double result = 0.0;

   result = 0.1*(-3.872983346207417*g1*Cos(ThetaW()) - 5*g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpconjHpVZHp() const
{
   double result = 0.0;

   result = 0.1*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpconjHpbarFdFuPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_0;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1;
      std::complex<double> tmp_2;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2 += Conj(Yu(j1,j2))*Uu(gI2,j1);
      }
      tmp_1 += tmp_2;
      tmp_0 += (Vd(gI1,j2)) * tmp_1;
   }
   result += tmp_0;

   return result;
}

std::complex<double> CLASSNAME::CpconjHpbarFdFuPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_4;
      std::complex<double> tmp_5;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_5 += Conj(Ud(gI1,j1))*Yd(j1,j2);
      }
      tmp_4 += tmp_5;
      tmp_3 += (Conj(Vu(gI2,j2))) * tmp_4;
   }
   result += tmp_3;

   return result;
}

double CLASSNAME::CpconjHpbarFeFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjHpbarFeFvPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_6;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_6 += Conj(Ue(gI1,j1))*Ye(j1,gI2);
   }
   result += tmp_6;

   return result;
}

double CLASSNAME::CpAhhhAh() const
{
   double result = 0.0;

   result = v*Lambdax;

   return result;
}

std::complex<double> CLASSNAME::CpAhbargWpgWp() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.25)*v*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpAhbargWpCgWpC() const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.25)*v*Sqr(g2);

   return result;
}

double CLASSNAME::CpAhAhAhAh() const
{
   double result = 0.0;

   result = 3*Lambdax;

   return result;
}

double CLASSNAME::CpAhAhhhhh() const
{
   double result = 0.0;

   result = Lambdax;

   return result;
}

std::complex<double> CLASSNAME::CpAhAhVZVZ() const
{
   std::complex<double> result;

   result = 0.1*(g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*
      Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

double CLASSNAME::CpAhAhconjHpHp() const
{
   double result = 0.0;

   result = Lambdax;

   return result;
}

double CLASSNAME::CpAhAhconjVWpVWp() const
{
   double result = 0.0;

   result = 0.5*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpAhconjVWpHp() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2;

   return result;
}

std::complex<double> CLASSNAME::CpAhbarFdFdPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_7;
   std::complex<double> tmp_8;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_9;
      std::complex<double> tmp_10;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_10 += Conj(Yd(j1,j2))*Ud(gI2,j1);
      }
      tmp_9 += tmp_10;
      tmp_8 += (Vd(gI1,j2)) * tmp_9;
   }
   tmp_7 += tmp_8;
   result += (std::complex<double>(0.,0.7071067811865475)) * tmp_7;

   return result;
}

std::complex<double> CLASSNAME::CpAhbarFdFdPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_11;
   std::complex<double> tmp_12;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_13;
      std::complex<double> tmp_14;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_14 += Conj(Ud(gI1,j1))*Yd(j1,j2);
      }
      tmp_13 += tmp_14;
      tmp_12 += (Conj(Vd(gI2,j2))) * tmp_13;
   }
   tmp_11 += tmp_12;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_11;

   return result;
}

std::complex<double> CLASSNAME::CpAhbarFeFePR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_15;
   std::complex<double> tmp_16;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_17;
      std::complex<double> tmp_18;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_18 += Conj(Ye(j1,j2))*Ue(gI2,j1);
      }
      tmp_17 += tmp_18;
      tmp_16 += (Ve(gI1,j2)) * tmp_17;
   }
   tmp_15 += tmp_16;
   result += (std::complex<double>(0.,0.7071067811865475)) * tmp_15;

   return result;
}

std::complex<double> CLASSNAME::CpAhbarFeFePL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_19;
   std::complex<double> tmp_20;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_21;
      std::complex<double> tmp_22;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_22 += Conj(Ue(gI1,j1))*Ye(j1,j2);
      }
      tmp_21 += tmp_22;
      tmp_20 += (Conj(Ve(gI2,j2))) * tmp_21;
   }
   tmp_19 += tmp_20;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_19;

   return result;
}

std::complex<double> CLASSNAME::CpAhbarFuFuPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_23;
   std::complex<double> tmp_24;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_25;
      std::complex<double> tmp_26;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_26 += Conj(Yu(j1,j2))*Uu(gI2,j1);
      }
      tmp_25 += tmp_26;
      tmp_24 += (Vu(gI1,j2)) * tmp_25;
   }
   tmp_23 += tmp_24;
   result += (std::complex<double>(0.,0.7071067811865475)) * tmp_23;

   return result;
}

std::complex<double> CLASSNAME::CpAhbarFuFuPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_27;
   std::complex<double> tmp_28;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_29;
      std::complex<double> tmp_30;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_30 += Conj(Uu(gI1,j1))*Yu(j1,j2);
      }
      tmp_29 += tmp_30;
      tmp_28 += (Conj(Vu(gI2,j2))) * tmp_29;
   }
   tmp_27 += tmp_28;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_27;

   return result;
}

double CLASSNAME::CphhAhAh() const
{
   double result = 0.0;

   result = v*Lambdax;

   return result;
}

double CLASSNAME::Cphhhhhh() const
{
   double result = 0.0;

   result = 3*v*Lambdax;

   return result;
}

double CLASSNAME::CphhVZVZ() const
{
   double result = 0.0;

   result = 0.5*v*Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CphhbargWpgWp() const
{
   double result = 0.0;

   result = -0.25*v*Sqr(g2);

   return result;
}

double CLASSNAME::CphhbargWpCgWpC() const
{
   double result = 0.0;

   result = -0.25*v*Sqr(g2);

   return result;
}

double CLASSNAME::CphhbargZgZ() const
{
   double result = 0.0;

   result = -0.25*v*Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))
      ;

   return result;
}

double CLASSNAME::CphhconjHpHp() const
{
   double result = 0.0;

   result = v*Lambdax;

   return result;
}

double CLASSNAME::CphhconjVWpVWp() const
{
   double result = 0.0;

   result = 0.5*v*Sqr(g2);

   return result;
}

double CLASSNAME::CphhhhAhAh() const
{
   double result = 0.0;

   result = Lambdax;

   return result;
}

double CLASSNAME::Cphhhhhhhh() const
{
   double result = 0.0;

   result = 3*Lambdax;

   return result;
}

std::complex<double> CLASSNAME::CphhhhVZVZ() const
{
   std::complex<double> result;

   result = 0.1*(g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*
      Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

double CLASSNAME::CphhhhconjHpHp() const
{
   double result = 0.0;

   result = Lambdax;

   return result;
}

double CLASSNAME::CphhhhconjVWpVWp() const
{
   double result = 0.0;

   result = 0.5*Sqr(g2);

   return result;
}

double CLASSNAME::CphhconjVWpHp() const
{
   double result = 0.0;

   result = 0.5*g2;

   return result;
}

std::complex<double> CLASSNAME::CphhbarFdFdPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_31;
   std::complex<double> tmp_32;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_33;
      std::complex<double> tmp_34;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_34 += Conj(Yd(j1,j2))*Ud(gI2,j1);
      }
      tmp_33 += tmp_34;
      tmp_32 += (Vd(gI1,j2)) * tmp_33;
   }
   tmp_31 += tmp_32;
   result += (0.7071067811865475) * tmp_31;

   return result;
}

std::complex<double> CLASSNAME::CphhbarFdFdPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_35;
   std::complex<double> tmp_36;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_37;
      std::complex<double> tmp_38;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_38 += Conj(Ud(gI1,j1))*Yd(j1,j2);
      }
      tmp_37 += tmp_38;
      tmp_36 += (Conj(Vd(gI2,j2))) * tmp_37;
   }
   tmp_35 += tmp_36;
   result += (0.7071067811865475) * tmp_35;

   return result;
}

std::complex<double> CLASSNAME::CphhbarFeFePR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_39;
   std::complex<double> tmp_40;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_41;
      std::complex<double> tmp_42;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_42 += Conj(Ye(j1,j2))*Ue(gI2,j1);
      }
      tmp_41 += tmp_42;
      tmp_40 += (Ve(gI1,j2)) * tmp_41;
   }
   tmp_39 += tmp_40;
   result += (0.7071067811865475) * tmp_39;

   return result;
}

std::complex<double> CLASSNAME::CphhbarFeFePL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_43;
   std::complex<double> tmp_44;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_45;
      std::complex<double> tmp_46;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_46 += Conj(Ue(gI1,j1))*Ye(j1,j2);
      }
      tmp_45 += tmp_46;
      tmp_44 += (Conj(Ve(gI2,j2))) * tmp_45;
   }
   tmp_43 += tmp_44;
   result += (0.7071067811865475) * tmp_43;

   return result;
}

std::complex<double> CLASSNAME::CphhbarFuFuPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_47;
   std::complex<double> tmp_48;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_49;
      std::complex<double> tmp_50;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_50 += Conj(Yu(j1,j2))*Uu(gI2,j1);
      }
      tmp_49 += tmp_50;
      tmp_48 += (Vu(gI1,j2)) * tmp_49;
   }
   tmp_47 += tmp_48;
   result += (-0.7071067811865475) * tmp_47;

   return result;
}

std::complex<double> CLASSNAME::CphhbarFuFuPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_51;
   std::complex<double> tmp_52;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_53;
      std::complex<double> tmp_54;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_54 += Conj(Uu(gI1,j1))*Yu(j1,j2);
      }
      tmp_53 += tmp_54;
      tmp_52 += (Conj(Vu(gI2,j2))) * tmp_53;
   }
   tmp_51 += tmp_52;
   result += (-0.7071067811865475) * tmp_51;

   return result;
}

double CLASSNAME::CpVZVZhh() const
{
   double result = 0.0;

   result = 0.5*v*Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbargWpgWp() const
{
   double result = 0.0;

   result = g2*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpVZbargWpCgWpC() const
{
   double result = 0.0;

   result = -(g2*Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVZconjHpHp() const
{
   double result = 0.0;

   result = 0.1*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZconjVWpHp() const
{
   double result = 0.0;

   result = -0.3872983346207417*g1*g2*v*Sin(ThetaW());

   return result;
}

std::complex<double> CLASSNAME::CpVZVZAhAh() const
{
   std::complex<double> result;

   result = 0.1*(g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*
      Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZhhhh() const
{
   std::complex<double> result;

   result = 0.1*(g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*
      Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjHpHp() const
{
   std::complex<double> result;

   result = 0.1*(-7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW()) + 5*Sqr(
      g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW())));

   return result;
}

double CLASSNAME::CpVZconjVWpVWp() const
{
   double result = 0.0;

   result = g2*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpVZbarFdFdPL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.16666666666666666*KroneckerDelta(gI1,gI2)*(3*g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbarFdFdPR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.2581988897471611*g1*KroneckerDelta(gI1,gI2)*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpVZbarFeFePL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW()) - 0.7745966692414834*
      g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbarFeFePR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.7745966692414834*g1*KroneckerDelta(gI1,gI2)*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpVZbarFuFuPL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.03333333333333333*KroneckerDelta(gI1,gI2)*(-15*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbarFuFuPR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5163977794943222*g1*KroneckerDelta(gI1,gI2)*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpVZbarFvFvPL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW()) + 0.7745966692414834
      *g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbarFvFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpVZVZconjVWpVWp1() const
{
   double result = 0.0;

   result = -2*Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVZVZconjVWpVWp2() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVZVZconjVWpVWp3() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWpHpAh() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2;

   return result;
}

double CLASSNAME::CpconjVWpHphh() const
{
   double result = 0.0;

   result = 0.5*g2;

   return result;
}

double CLASSNAME::CpconjVWpVPHp() const
{
   double result = 0.0;

   result = 0.3872983346207417*g1*g2*v*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpconjVWpVWphh() const
{
   double result = 0.0;

   result = 0.5*v*Sqr(g2);

   return result;
}

double CLASSNAME::CpconjVWpVZHp() const
{
   double result = 0.0;

   result = -0.3872983346207417*g1*g2*v*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpconjVWpbargPgWp() const
{
   double result = 0.0;

   result = -(g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWpbargWpCgP() const
{
   double result = 0.0;

   result = g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpconjVWpbargWpCgZ() const
{
   double result = 0.0;

   result = g2*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpconjVWpbargZgWp() const
{
   double result = 0.0;

   result = -(g2*Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVWpconjVWpAhAh() const
{
   double result = 0.0;

   result = 0.5*Sqr(g2);

   return result;
}

double CLASSNAME::CpVWpconjVWphhhh() const
{
   double result = 0.0;

   result = 0.5*Sqr(g2);

   return result;
}

double CLASSNAME::CpVWpconjVWpconjHpHp() const
{
   double result = 0.0;

   result = 0.5*Sqr(g2);

   return result;
}

double CLASSNAME::CpconjVWpVWpVP() const
{
   double result = 0.0;

   result = g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpconjVWpVZVWp() const
{
   double result = 0.0;

   result = -(g2*Cos(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWpbarFdFuPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_55;
   std::complex<double> tmp_56;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_56 += Conj(Vu(gI2,j1))*Vd(gI1,j1);
   }
   tmp_55 += tmp_56;
   result += (-0.7071067811865475*g2) * tmp_55;

   return result;
}

double CLASSNAME::CpconjVWpbarFdFuPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjVWpbarFeFvPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.7071067811865475*g2*Ve(gI1,gI2);
   }

   return result;
}

double CLASSNAME::CpconjVWpbarFeFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpVWpconjVWpVPVP1() const
{
   double result = 0.0;

   result = -2*Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVWpconjVWpVPVP2() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVWpconjVWpVPVP3() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVWpconjVWpVZVZ1() const
{
   double result = 0.0;

   result = -2*Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVWpconjVWpVZVZ2() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVWpconjVWpVZVZ3() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVWpconjVWpconjVWpVWp1() const
{
   double result = 0.0;

   result = -Sqr(g2);

   return result;
}

double CLASSNAME::CpVWpconjVWpconjVWpVWp2() const
{
   double result = 0.0;

   result = -Sqr(g2);

   return result;
}

double CLASSNAME::CpVWpconjVWpconjVWpVWp3() const
{
   double result = 0.0;

   result = 2*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPL(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_57;
      std::complex<double> tmp_58;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_58 += Conj(Vd(gI1,j2))*Yd(gO2,j2);
      }
      tmp_57 += tmp_58;
      result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_57;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPR(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_59;
      std::complex<double> tmp_60;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_60 += Conj(Yd(j1,gO1))*Ud(gI1,j1);
      }
      tmp_59 += tmp_60;
      result += (std::complex<double>(0.,0.7071067811865475)) * tmp_59;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdhhFdPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_61;
      std::complex<double> tmp_62;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_62 += Conj(Vd(gI2,j2))*Yd(gO2,j2);
      }
      tmp_61 += tmp_62;
      result += (0.7071067811865475) * tmp_61;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdhhFdPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_63;
      std::complex<double> tmp_64;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_64 += Conj(Yd(j1,gO1))*Ud(gI2,j1);
      }
      tmp_63 += tmp_64;
      result += (0.7071067811865475) * tmp_63;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVGFdPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -(g3*Ud(gI2,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVGFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -(g3*Conj(Vd(gI2,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVPFdPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.2581988897471611*g1*Cos(ThetaW())*Ud(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVPFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.12909944487358055*g1*Conj(Vd(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += 0.5*g2*Conj(Vd(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVZFdPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.2581988897471611*g1*Sin(ThetaW())*Ud(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVZFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.5*g2*Conj(Vd(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += 0.12909944487358055*g1*Conj(Vd(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdconjHpFuPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_65;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_65 += Conj(Vu(gI2,j2))*Yd(gO2,j2);
      }
      result += tmp_65;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdconjHpFuPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_66;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_66 += Conj(Yu(j1,gO1))*Uu(gI2,j1);
      }
      result += tmp_66;
   }

   return result;
}

double CLASSNAME::CpbarUFdconjVWpFuPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdconjVWpFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(Vu(gI2,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPL(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_67;
      std::complex<double> tmp_68;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_68 += Conj(Vu(gI1,j2))*Yu(gO2,j2);
      }
      tmp_67 += tmp_68;
      result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_67;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPR(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_69;
      std::complex<double> tmp_70;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_70 += Conj(Yu(j1,gO1))*Uu(gI1,j1);
      }
      tmp_69 += tmp_70;
      result += (std::complex<double>(0.,0.7071067811865475)) * tmp_69;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuhhFuPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_71;
      std::complex<double> tmp_72;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_72 += Conj(Vu(gI2,j2))*Yu(gO2,j2);
      }
      tmp_71 += tmp_72;
      result += (-0.7071067811865475) * tmp_71;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuhhFuPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_73;
      std::complex<double> tmp_74;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_74 += Conj(Yu(j1,gO1))*Uu(gI2,j1);
      }
      tmp_73 += tmp_74;
      result += (-0.7071067811865475) * tmp_73;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuHpFdPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_75;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_75 += Conj(Vd(gI2,j2))*Yu(gO2,j2);
      }
      result += tmp_75;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuHpFdPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_76;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_76 += Conj(Yd(j1,gO1))*Ud(gI2,j1);
      }
      result += tmp_76;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVGFuPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -(g3*Uu(gI2,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVGFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -(g3*Conj(Vu(gI2,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVPFuPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.5163977794943222*g1*Cos(ThetaW())*Uu(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVPFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.12909944487358055*g1*Conj(Vu(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += -0.5*g2*Conj(Vu(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

double CLASSNAME::CpbarUFuVWpFdPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVWpFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(Vd(gI2,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVZFuPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.5163977794943222*g1*Sin(ThetaW())*Uu(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVZFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.5*g2*Conj(Vu(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += 0.12909944487358055*g1*Conj(Vu(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPL(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_77;
      std::complex<double> tmp_78;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_78 += Conj(Ve(gI1,j2))*Ye(gO2,j2);
      }
      tmp_77 += tmp_78;
      result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_77;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPR(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_79;
      std::complex<double> tmp_80;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_80 += Conj(Ye(j1,gO1))*Ue(gI1,j1);
      }
      tmp_79 += tmp_80;
      result += (std::complex<double>(0.,0.7071067811865475)) * tmp_79;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFehhFePL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_81;
      std::complex<double> tmp_82;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_82 += Conj(Ve(gI2,j2))*Ye(gO2,j2);
      }
      tmp_81 += tmp_82;
      result += (0.7071067811865475) * tmp_81;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFehhFePR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_83;
      std::complex<double> tmp_84;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_84 += Conj(Ye(j1,gO1))*Ue(gI2,j1);
      }
      tmp_83 += tmp_84;
      result += (0.7071067811865475) * tmp_83;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeVPFePR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.7745966692414834*g1*Cos(ThetaW())*Ue(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeVPFePL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.3872983346207417*g1*Conj(Ve(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += 0.5*g2*Conj(Ve(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeVZFePR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.7745966692414834*g1*Sin(ThetaW())*Ue(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeVZFePL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.5*g2*Conj(Ve(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += -0.3872983346207417*g1*Conj(Ve(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeconjHpFvPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += Ye(gO2,gI2);
   }

   return result;
}

double CLASSNAME::CpbarUFeconjHpFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpbarUFeconjVWpFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpbarUFeconjVWpFvPL(unsigned gO1, unsigned gI2) const
{
   double result = 0.0;

   if (gI2 < 3) {
      result += -0.7071067811865475*g2*KroneckerDelta(gI2,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPL(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   std::complex<double> tmp_85;
   std::complex<double> tmp_86;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_87;
      std::complex<double> tmp_88;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_88 += Conj(Ud(gO2,j1))*Yd(j1,j2);
      }
      tmp_87 += tmp_88;
      tmp_86 += (Conj(Vd(gI1,j2))) * tmp_87;
   }
   tmp_85 += tmp_86;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_85;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPR(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   std::complex<double> tmp_89;
   std::complex<double> tmp_90;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_91;
      std::complex<double> tmp_92;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_92 += Conj(Yd(j1,j2))*Ud(gI1,j1);
      }
      tmp_91 += tmp_92;
      tmp_90 += (Vd(gO1,j2)) * tmp_91;
   }
   tmp_89 += tmp_90;
   result += (std::complex<double>(0.,0.7071067811865475)) * tmp_89;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdhhFdPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_93;
   std::complex<double> tmp_94;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_95;
      std::complex<double> tmp_96;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_96 += Conj(Ud(gO2,j1))*Yd(j1,j2);
      }
      tmp_95 += tmp_96;
      tmp_94 += (Conj(Vd(gI2,j2))) * tmp_95;
   }
   tmp_93 += tmp_94;
   result += (0.7071067811865475) * tmp_93;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdhhFdPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_97;
   std::complex<double> tmp_98;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_99;
      std::complex<double> tmp_100;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_100 += Conj(Yd(j1,j2))*Ud(gI2,j1);
      }
      tmp_99 += tmp_100;
      tmp_98 += (Vd(gO1,j2)) * tmp_99;
   }
   tmp_97 += tmp_98;
   result += (0.7071067811865475) * tmp_97;

   return result;
}

double CLASSNAME::CpbarFdVZFdPR(unsigned gO2, unsigned gI2) const
{
   double result = 0.0;

   result = -0.2581988897471611*g1*KroneckerDelta(gI2,gO2)*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpbarFdVZFdPL(unsigned gO1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.16666666666666666*KroneckerDelta(gI2,gO1)*(3*g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdconjHpFuPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_101;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_102;
      std::complex<double> tmp_103;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_103 += Conj(Ud(gO2,j1))*Yd(j1,j2);
      }
      tmp_102 += tmp_103;
      tmp_101 += (Conj(Vu(gI2,j2))) * tmp_102;
   }
   result += tmp_101;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdconjHpFuPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_104;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_105;
      std::complex<double> tmp_106;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_106 += Conj(Yu(j1,j2))*Uu(gI2,j1);
      }
      tmp_105 += tmp_106;
      tmp_104 += (Vd(gO1,j2)) * tmp_105;
   }
   result += tmp_104;

   return result;
}

double CLASSNAME::CpbarFdconjVWpFuPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdconjVWpFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_107;
   std::complex<double> tmp_108;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_108 += Conj(Vu(gI2,j1))*Vd(gO1,j1);
   }
   tmp_107 += tmp_108;
   result += (-0.7071067811865475*g2) * tmp_107;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPL(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   std::complex<double> tmp_109;
   std::complex<double> tmp_110;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_111;
      std::complex<double> tmp_112;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_112 += Conj(Ue(gO2,j1))*Ye(j1,j2);
      }
      tmp_111 += tmp_112;
      tmp_110 += (Conj(Ve(gI1,j2))) * tmp_111;
   }
   tmp_109 += tmp_110;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_109;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPR(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   std::complex<double> tmp_113;
   std::complex<double> tmp_114;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_115;
      std::complex<double> tmp_116;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_116 += Conj(Ye(j1,j2))*Ue(gI1,j1);
      }
      tmp_115 += tmp_116;
      tmp_114 += (Ve(gO1,j2)) * tmp_115;
   }
   tmp_113 += tmp_114;
   result += (std::complex<double>(0.,0.7071067811865475)) * tmp_113;

   return result;
}

std::complex<double> CLASSNAME::CpbarFehhFePL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_117;
   std::complex<double> tmp_118;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_119;
      std::complex<double> tmp_120;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_120 += Conj(Ue(gO2,j1))*Ye(j1,j2);
      }
      tmp_119 += tmp_120;
      tmp_118 += (Conj(Ve(gI2,j2))) * tmp_119;
   }
   tmp_117 += tmp_118;
   result += (0.7071067811865475) * tmp_117;

   return result;
}

std::complex<double> CLASSNAME::CpbarFehhFePR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_121;
   std::complex<double> tmp_122;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_123;
      std::complex<double> tmp_124;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_124 += Conj(Ye(j1,j2))*Ue(gI2,j1);
      }
      tmp_123 += tmp_124;
      tmp_122 += (Ve(gO1,j2)) * tmp_123;
   }
   tmp_121 += tmp_122;
   result += (0.7071067811865475) * tmp_121;

   return result;
}

double CLASSNAME::CpbarFeVZFePR(unsigned gO2, unsigned gI2) const
{
   double result = 0.0;

   result = -0.7745966692414834*g1*KroneckerDelta(gI2,gO2)*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpbarFeVZFePL(unsigned gO1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5*KroneckerDelta(gI2,gO1)*(g2*Cos(ThetaW()) - 0.7745966692414834*
      g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeconjHpFvPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_125;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_125 += Conj(Ue(gO2,j1))*Ye(j1,gI2);
   }
   result += tmp_125;

   return result;
}

double CLASSNAME::CpbarFeconjHpFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpbarFeconjVWpFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeconjVWpFvPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.7071067811865475*g2*Ve(gO1,gI2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPL(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   std::complex<double> tmp_126;
   std::complex<double> tmp_127;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_128;
      std::complex<double> tmp_129;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_129 += Conj(Uu(gO2,j1))*Yu(j1,j2);
      }
      tmp_128 += tmp_129;
      tmp_127 += (Conj(Vu(gI1,j2))) * tmp_128;
   }
   tmp_126 += tmp_127;
   result += (std::complex<double>(0.,-0.7071067811865475)) * tmp_126;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPR(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   std::complex<double> tmp_130;
   std::complex<double> tmp_131;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_132;
      std::complex<double> tmp_133;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_133 += Conj(Yu(j1,j2))*Uu(gI1,j1);
      }
      tmp_132 += tmp_133;
      tmp_131 += (Vu(gO1,j2)) * tmp_132;
   }
   tmp_130 += tmp_131;
   result += (std::complex<double>(0.,0.7071067811865475)) * tmp_130;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuhhFuPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_134;
   std::complex<double> tmp_135;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_136;
      std::complex<double> tmp_137;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_137 += Conj(Uu(gO2,j1))*Yu(j1,j2);
      }
      tmp_136 += tmp_137;
      tmp_135 += (Conj(Vu(gI2,j2))) * tmp_136;
   }
   tmp_134 += tmp_135;
   result += (-0.7071067811865475) * tmp_134;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuhhFuPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_138;
   std::complex<double> tmp_139;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_140;
      std::complex<double> tmp_141;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_141 += Conj(Yu(j1,j2))*Uu(gI2,j1);
      }
      tmp_140 += tmp_141;
      tmp_139 += (Vu(gO1,j2)) * tmp_140;
   }
   tmp_138 += tmp_139;
   result += (-0.7071067811865475) * tmp_138;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuHpFdPL(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_142;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_143;
      std::complex<double> tmp_144;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_144 += Conj(Uu(gO2,j1))*Yu(j1,j2);
      }
      tmp_143 += tmp_144;
      tmp_142 += (Conj(Vd(gI2,j2))) * tmp_143;
   }
   result += tmp_142;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuHpFdPR(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_145;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_146;
      std::complex<double> tmp_147;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_147 += Conj(Yd(j1,j2))*Ud(gI2,j1);
      }
      tmp_146 += tmp_147;
      tmp_145 += (Vu(gO1,j2)) * tmp_146;
   }
   result += tmp_145;

   return result;
}

double CLASSNAME::CpbarFuVPFuPR(unsigned gO2, unsigned gI2) const
{
   double result = 0.0;

   result = -0.5163977794943222*g1*Cos(ThetaW())*KroneckerDelta(gI2,gO2);

   return result;
}

double CLASSNAME::CpbarFuVPFuPL(unsigned gO1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.16666666666666666*KroneckerDelta(gI2,gO1)*(0.7745966692414834*g1
      *Cos(ThetaW()) + 3*g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFuVWpFdPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuVWpFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_148;
   std::complex<double> tmp_149;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_149 += Conj(Vd(gI2,j1))*Vu(gO1,j1);
   }
   tmp_148 += tmp_149;
   result += (-0.7071067811865475*g2) * tmp_148;

   return result;
}

double CLASSNAME::CpbarFuVZFuPR(unsigned gO2, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5163977794943222*g1*KroneckerDelta(gI2,gO2)*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpbarFuVZFuPL(unsigned gO1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.03333333333333333*KroneckerDelta(gI2,gO1)*(-15*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()));

   return result;
}


std::complex<double> CLASSNAME::self_energy_Hp(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpconjHpHphh())*B0(p,MHp,Mhh);
   result += 4*AbsSqr(CpconjHpVWpVP())*B0(p,0,MVWp);
   result += 4*AbsSqr(CpconjHpVZVWp())*B0(p,MVWp,MVZ);
   result += -0.5*A0(MAh)*CpHpconjHpAhAh();
   result += -(A0(MHp)*CpHpconjHpconjHpHp());
   result += 4*A0(MVWp)*CpHpconjHpconjVWpVWp();
   result += -0.5*A0(Mhh)*CpHpconjHphhhh();
   result += 2*A0(MVZ)*CpHpconjHpVZVZ();
   result += -(B0(p,MVZ,MVWp)*CpconjHpbargWpCgZ()*CpHpgWpCbargZ());
   result += -(B0(p,MVWp,MVZ)*CpconjHpbargZgWp()*CpHpgZbargWp());
   result += AbsSqr(CpconjHpVWpAh())*F0(p,MAh,MVWp);
   result += AbsSqr(CpconjHpVWphh())*F0(p,Mhh,MVWp);
   result += AbsSqr(CpconjHpVPHp())*F0(p,MHp,0);
   result += AbsSqr(CpconjHpVZHp())*F0(p,MHp,MVZ);
   std::complex<double> tmp_150;
   std::complex<double> tmp_151;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_152;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_152 += (AbsSqr(CpconjHpbarFdFuPL(gI1,gI2)) + AbsSqr(
            CpconjHpbarFdFuPR(gI1,gI2)))*G0(p,MFd(gI1),MFu(gI2));
      }
      tmp_151 += tmp_152;
   }
   tmp_150 += tmp_151;
   result += (3) * tmp_150;
   std::complex<double> tmp_153;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_154;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_154 += (AbsSqr(CpconjHpbarFeFvPL(gI1,gI2)) + AbsSqr(
            CpconjHpbarFeFvPR(gI1,gI2)))*G0(p,MFe(gI1),MFv(gI2));
      }
      tmp_153 += tmp_154;
   }
   result += tmp_153;
   std::complex<double> tmp_155;
   std::complex<double> tmp_156;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_157;
      std::complex<double> tmp_158;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_158 += B0(p,MFd(gI1),MFu(gI2))*(Conj(CpconjHpbarFdFuPR(gI1,
            gI2))*CpconjHpbarFdFuPL(gI1,gI2) + Conj(CpconjHpbarFdFuPL(gI1,gI2))*
            CpconjHpbarFdFuPR(gI1,gI2))*MFu(gI2);
      }
      tmp_157 += tmp_158;
      tmp_156 += (MFd(gI1)) * tmp_157;
   }
   tmp_155 += tmp_156;
   result += (-6) * tmp_155;
   std::complex<double> tmp_159;
   std::complex<double> tmp_160;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_161;
      std::complex<double> tmp_162;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_162 += B0(p,MFe(gI1),MFv(gI2))*(Conj(CpconjHpbarFeFvPR(gI1,
            gI2))*CpconjHpbarFeFvPL(gI1,gI2) + Conj(CpconjHpbarFeFvPL(gI1,gI2))*
            CpconjHpbarFeFvPR(gI1,gI2))*MFv(gI2);
      }
      tmp_161 += tmp_162;
      tmp_160 += (MFe(gI1)) * tmp_161;
   }
   tmp_159 += tmp_160;
   result += (-2) * tmp_159;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Ah(double p ) const
{
   std::complex<double> result;

   result += -0.5*A0(MAh)*CpAhAhAhAh();
   result += -(A0(MHp)*CpAhAhconjHpHp());
   result += 4*A0(MVWp)*CpAhAhconjVWpVWp();
   result += -0.5*A0(Mhh)*CpAhAhhhhh();
   result += 2*A0(MVZ)*CpAhAhVZVZ();
   result += -(B0(p,MVWp,MVWp)*Sqr(CpAhbargWpCgWpC()));
   result += -(B0(p,MVWp,MVWp)*Sqr(CpAhbargWpgWp()));
   result += AbsSqr(CpAhhhAh())*B0(p,Mhh,MAh);
   result += 2*AbsSqr(CpAhconjVWpHp())*F0(p,MHp,MVWp);
   std::complex<double> tmp_163;
   std::complex<double> tmp_164;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_165;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_165 += (AbsSqr(CpAhbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpAhbarFdFdPR(gI1,gI2)))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_164 += tmp_165;
   }
   tmp_163 += tmp_164;
   result += (3) * tmp_163;
   std::complex<double> tmp_166;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_167;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_167 += (AbsSqr(CpAhbarFeFePL(gI1,gI2)) + AbsSqr(
            CpAhbarFeFePR(gI1,gI2)))*G0(p,MFe(gI1),MFe(gI2));
      }
      tmp_166 += tmp_167;
   }
   result += tmp_166;
   std::complex<double> tmp_168;
   std::complex<double> tmp_169;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_170;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_170 += (AbsSqr(CpAhbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpAhbarFuFuPR(gI1,gI2)))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_169 += tmp_170;
   }
   tmp_168 += tmp_169;
   result += (3) * tmp_168;
   std::complex<double> tmp_171;
   std::complex<double> tmp_172;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_173;
      std::complex<double> tmp_174;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_174 += B0(p,MFd(gI1),MFd(gI2))*(Conj(CpAhbarFdFdPR(gI1,gI2))
            *CpAhbarFdFdPL(gI1,gI2) + Conj(CpAhbarFdFdPL(gI1,gI2))*CpAhbarFdFdPR(
            gI1,gI2))*MFd(gI2);
      }
      tmp_173 += tmp_174;
      tmp_172 += (MFd(gI1)) * tmp_173;
   }
   tmp_171 += tmp_172;
   result += (-6) * tmp_171;
   std::complex<double> tmp_175;
   std::complex<double> tmp_176;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_177;
      std::complex<double> tmp_178;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_178 += B0(p,MFe(gI1),MFe(gI2))*(Conj(CpAhbarFeFePR(gI1,gI2))
            *CpAhbarFeFePL(gI1,gI2) + Conj(CpAhbarFeFePL(gI1,gI2))*CpAhbarFeFePR(
            gI1,gI2))*MFe(gI2);
      }
      tmp_177 += tmp_178;
      tmp_176 += (MFe(gI1)) * tmp_177;
   }
   tmp_175 += tmp_176;
   result += (-2) * tmp_175;
   std::complex<double> tmp_179;
   std::complex<double> tmp_180;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_181;
      std::complex<double> tmp_182;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_182 += B0(p,MFu(gI1),MFu(gI2))*(Conj(CpAhbarFuFuPR(gI1,gI2))
            *CpAhbarFuFuPL(gI1,gI2) + Conj(CpAhbarFuFuPL(gI1,gI2))*CpAhbarFuFuPR(
            gI1,gI2))*MFu(gI2);
      }
      tmp_181 += tmp_182;
      tmp_180 += (MFu(gI1)) * tmp_181;
   }
   tmp_179 += tmp_180;
   result += (-6) * tmp_179;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_hh(double p ) const
{
   std::complex<double> result;

   result += 0.5*AbsSqr(CphhAhAh())*B0(p,MAh,MAh);
   result += -(B0(p,MVWp,MVWp)*Sqr(CphhbargWpCgWpC()));
   result += -(B0(p,MVWp,MVWp)*Sqr(CphhbargWpgWp()));
   result += -(B0(p,MVZ,MVZ)*Sqr(CphhbargZgZ()));
   result += AbsSqr(CphhconjHpHp())*B0(p,MHp,MHp);
   result += 4*AbsSqr(CphhconjVWpVWp())*B0(p,MVWp,MVWp);
   result += -0.5*A0(MAh)*CphhhhAhAh();
   result += -(A0(MHp)*CphhhhconjHpHp());
   result += 4*A0(MVWp)*CphhhhconjVWpVWp();
   result += 0.5*AbsSqr(Cphhhhhh())*B0(p,Mhh,Mhh);
   result += -0.5*A0(Mhh)*Cphhhhhhhh();
   result += 2*A0(MVZ)*CphhhhVZVZ();
   result += 2*AbsSqr(CphhVZVZ())*B0(p,MVZ,MVZ);
   result += 2*AbsSqr(CphhconjVWpHp())*F0(p,MHp,MVWp);
   std::complex<double> tmp_183;
   std::complex<double> tmp_184;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_185;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_185 += (AbsSqr(CphhbarFdFdPL(gI1,gI2)) + AbsSqr(
            CphhbarFdFdPR(gI1,gI2)))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_184 += tmp_185;
   }
   tmp_183 += tmp_184;
   result += (3) * tmp_183;
   std::complex<double> tmp_186;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_187;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_187 += (AbsSqr(CphhbarFeFePL(gI1,gI2)) + AbsSqr(
            CphhbarFeFePR(gI1,gI2)))*G0(p,MFe(gI1),MFe(gI2));
      }
      tmp_186 += tmp_187;
   }
   result += tmp_186;
   std::complex<double> tmp_188;
   std::complex<double> tmp_189;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_190;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_190 += (AbsSqr(CphhbarFuFuPL(gI1,gI2)) + AbsSqr(
            CphhbarFuFuPR(gI1,gI2)))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_189 += tmp_190;
   }
   tmp_188 += tmp_189;
   result += (3) * tmp_188;
   std::complex<double> tmp_191;
   std::complex<double> tmp_192;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_193;
      std::complex<double> tmp_194;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_194 += B0(p,MFd(gI1),MFd(gI2))*(Conj(CphhbarFdFdPR(gI1,gI2))
            *CphhbarFdFdPL(gI1,gI2) + Conj(CphhbarFdFdPL(gI1,gI2))*CphhbarFdFdPR(
            gI1,gI2))*MFd(gI2);
      }
      tmp_193 += tmp_194;
      tmp_192 += (MFd(gI1)) * tmp_193;
   }
   tmp_191 += tmp_192;
   result += (-6) * tmp_191;
   std::complex<double> tmp_195;
   std::complex<double> tmp_196;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_197;
      std::complex<double> tmp_198;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_198 += B0(p,MFe(gI1),MFe(gI2))*(Conj(CphhbarFeFePR(gI1,gI2))
            *CphhbarFeFePL(gI1,gI2) + Conj(CphhbarFeFePL(gI1,gI2))*CphhbarFeFePR(
            gI1,gI2))*MFe(gI2);
      }
      tmp_197 += tmp_198;
      tmp_196 += (MFe(gI1)) * tmp_197;
   }
   tmp_195 += tmp_196;
   result += (-2) * tmp_195;
   std::complex<double> tmp_199;
   std::complex<double> tmp_200;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_201;
      std::complex<double> tmp_202;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_202 += B0(p,MFu(gI1),MFu(gI2))*(Conj(CphhbarFuFuPR(gI1,gI2))
            *CphhbarFuFuPL(gI1,gI2) + Conj(CphhbarFuFuPL(gI1,gI2))*CphhbarFuFuPR(
            gI1,gI2))*MFu(gI2);
      }
      tmp_201 += tmp_202;
      tmp_200 += (MFu(gI1)) * tmp_201;
   }
   tmp_199 += tmp_200;
   result += (-6) * tmp_199;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VZ(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpVZbargWpCgWpC())*B00(p,MVWp,MVWp);
   result += AbsSqr(CpVZbargWpgWp())*B00(p,MVWp,MVWp);
   result += -4*AbsSqr(CpVZconjHpHp())*B00(p,MHp,MHp);
   result += 2*AbsSqr(CpVZconjVWpHp())*B0(p,MVWp,MHp);
   result += 0.5*A0(MAh)*CpVZVZAhAh();
   result += A0(MHp)*CpVZVZconjHpHp();
   result += -(A0(MVWp)*(4*CpVZVZconjVWpVWp1() + CpVZVZconjVWpVWp2() +
      CpVZVZconjVWpVWp3()));
   result += AbsSqr(CpVZVZhh())*B0(p,MVZ,Mhh);
   result += 0.5*A0(Mhh)*CpVZVZhhhh();
   std::complex<double> tmp_203;
   std::complex<double> tmp_204;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_205;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_205 += (AbsSqr(CpVZbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpVZbarFdFdPR(gI1,gI2)))*H0(p,MFd(gI1),MFd(gI2));
         tmp_205 += 4*B0(p,MFd(gI1),MFd(gI2))*MFd(gI1)*MFd(gI2)*Re(Conj(
            CpVZbarFdFdPL(gI1,gI2))*CpVZbarFdFdPR(gI1,gI2));
      }
      tmp_204 += tmp_205;
   }
   tmp_203 += tmp_204;
   result += (3) * tmp_203;
   std::complex<double> tmp_206;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_207;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_207 += (AbsSqr(CpVZbarFeFePL(gI1,gI2)) + AbsSqr(
            CpVZbarFeFePR(gI1,gI2)))*H0(p,MFe(gI1),MFe(gI2));
         tmp_207 += 4*B0(p,MFe(gI1),MFe(gI2))*MFe(gI1)*MFe(gI2)*Re(Conj(
            CpVZbarFeFePL(gI1,gI2))*CpVZbarFeFePR(gI1,gI2));
      }
      tmp_206 += tmp_207;
   }
   result += tmp_206;
   std::complex<double> tmp_208;
   std::complex<double> tmp_209;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_210;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_210 += (AbsSqr(CpVZbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpVZbarFuFuPR(gI1,gI2)))*H0(p,MFu(gI1),MFu(gI2));
         tmp_210 += 4*B0(p,MFu(gI1),MFu(gI2))*MFu(gI1)*MFu(gI2)*Re(Conj(
            CpVZbarFuFuPL(gI1,gI2))*CpVZbarFuFuPR(gI1,gI2));
      }
      tmp_209 += tmp_210;
   }
   tmp_208 += tmp_209;
   result += (3) * tmp_208;
   std::complex<double> tmp_211;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_212;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_212 += (AbsSqr(CpVZbarFvFvPL(gI1,gI2)) + AbsSqr(
            CpVZbarFvFvPR(gI1,gI2)))*H0(p,MFv(gI1),MFv(gI2));
         tmp_212 += 4*B0(p,MFv(gI1),MFv(gI2))*MFv(gI1)*MFv(gI2)*Re(Conj(
            CpVZbarFvFvPL(gI1,gI2))*CpVZbarFvFvPR(gI1,gI2));
      }
      tmp_211 += tmp_212;
   }
   result += tmp_211;
   result += -(AbsSqr(CpVZconjVWpVWp())*(2*A0(MVWp) + 10*B00(p,MVWp,MVWp) + B0(
      p,MVWp,MVWp)*(2*Sqr(MVWp) + 4*Sqr(p))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VWp(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpconjVWpbargPgWp())*B00(p,MVWp,MVP);
   result += AbsSqr(CpconjVWpbargWpCgP())*B00(p,MVP,MVWp);
   result += AbsSqr(CpconjVWpbargWpCgZ())*B00(p,MVZ,MVWp);
   result += AbsSqr(CpconjVWpbargZgWp())*B00(p,MVWp,MVZ);
   result += -4*AbsSqr(CpconjVWpHpAh())*B00(p,MAh,MHp);
   result += -4*AbsSqr(CpconjVWpHphh())*B00(p,Mhh,MHp);
   result += AbsSqr(CpconjVWpVPHp())*B0(p,0,MHp);
   result += AbsSqr(CpconjVWpVWphh())*B0(p,MVWp,Mhh);
   result += AbsSqr(CpconjVWpVZHp())*B0(p,MVZ,MHp);
   result += 0.5*A0(MAh)*CpVWpconjVWpAhAh();
   result += A0(MHp)*CpVWpconjVWpconjHpHp();
   result += -(A0(MVWp)*(4*CpVWpconjVWpconjVWpVWp1() + CpVWpconjVWpconjVWpVWp2(
      ) + CpVWpconjVWpconjVWpVWp3()));
   result += 0.5*A0(Mhh)*CpVWpconjVWphhhh();
   result += 0;
   result += -0.5*A0(MVZ)*(4*CpVWpconjVWpVZVZ1() + CpVWpconjVWpVZVZ2() +
      CpVWpconjVWpVZVZ3());
   std::complex<double> tmp_213;
   std::complex<double> tmp_214;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_215;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_215 += (AbsSqr(CpconjVWpbarFdFuPL(gI1,gI2)) + AbsSqr(
            CpconjVWpbarFdFuPR(gI1,gI2)))*H0(p,MFd(gI1),MFu(gI2));
         tmp_215 += 4*B0(p,MFd(gI1),MFu(gI2))*MFd(gI1)*MFu(gI2)*Re(Conj(
            CpconjVWpbarFdFuPL(gI1,gI2))*CpconjVWpbarFdFuPR(gI1,gI2));
      }
      tmp_214 += tmp_215;
   }
   tmp_213 += tmp_214;
   result += (3) * tmp_213;
   std::complex<double> tmp_216;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_217;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_217 += (AbsSqr(CpconjVWpbarFeFvPL(gI1,gI2)) + AbsSqr(
            CpconjVWpbarFeFvPR(gI1,gI2)))*H0(p,MFe(gI1),MFv(gI2));
         tmp_217 += 4*B0(p,MFe(gI1),MFv(gI2))*MFe(gI1)*MFv(gI2)*Re(Conj(
            CpconjVWpbarFeFvPL(gI1,gI2))*CpconjVWpbarFeFvPR(gI1,gI2));
      }
      tmp_216 += tmp_217;
   }
   result += tmp_216;
   result += -(AbsSqr(CpconjVWpVWpVP())*(A0(MVWp) + 10*B00(p,MVWp,0) + B0(p,
      MVWp,0)*(Sqr(MVWp) + 4*Sqr(p))));
   result += -(AbsSqr(CpconjVWpVZVWp())*(A0(MVWp) + A0(MVZ) + 10*B00(p,MVZ,MVWp
      ) + B0(p,MVZ,MVWp)*(Sqr(MVWp) + Sqr(MVZ) + 4*Sqr(p))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_218;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_218 += B0(p,MFd(gI1),MAh)*Conj(CpbarUFdFdAhPL(gO2,gI1))*
         CpbarUFdFdAhPR(gO1,gI1)*MFd(gI1);
   }
   result += tmp_218;
   std::complex<double> tmp_219;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_219 += B0(p,MFd(gI2),Mhh)*Conj(CpbarUFdhhFdPL(gO2,gI2))*
         CpbarUFdhhFdPR(gO1,gI2)*MFd(gI2);
   }
   result += tmp_219;
   std::complex<double> tmp_220;
   std::complex<double> tmp_221;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_221 += B0(p,MFd(gI2),0)*Conj(CpbarUFdVGFdPR(gO2,gI2))*
         CpbarUFdVGFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_220 += tmp_221;
   result += (-5.333333333333333) * tmp_220;
   std::complex<double> tmp_222;
   std::complex<double> tmp_223;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_223 += B0(p,MFd(gI2),0)*Conj(CpbarUFdVPFdPR(gO2,gI2))*
         CpbarUFdVPFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_222 += tmp_223;
   result += (-4) * tmp_222;
   std::complex<double> tmp_224;
   std::complex<double> tmp_225;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_225 += B0(p,MFd(gI2),MVZ)*Conj(CpbarUFdVZFdPR(gO2,gI2))*
         CpbarUFdVZFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_224 += tmp_225;
   result += (-4) * tmp_224;
   std::complex<double> tmp_226;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_226 += B0(p,MFu(gI2),MHp)*Conj(CpbarUFdconjHpFuPL(gO2,gI2))*
         CpbarUFdconjHpFuPR(gO1,gI2)*MFu(gI2);
   }
   result += tmp_226;
   std::complex<double> tmp_227;
   std::complex<double> tmp_228;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_228 += B0(p,MFu(gI2),MVWp)*Conj(CpbarUFdconjVWpFuPR(gO2,gI2))*
         CpbarUFdconjVWpFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_227 += tmp_228;
   result += (-4) * tmp_227;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_229;
   std::complex<double> tmp_230;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_230 += B1(p,MFd(gI1),MAh)*Conj(CpbarUFdFdAhPR(gO2,gI1))*
         CpbarUFdFdAhPR(gO1,gI1);
   }
   tmp_229 += tmp_230;
   result += (-0.5) * tmp_229;
   std::complex<double> tmp_231;
   std::complex<double> tmp_232;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_232 += B1(p,MFu(gI2),MHp)*Conj(CpbarUFdconjHpFuPR(gO2,gI2))*
         CpbarUFdconjHpFuPR(gO1,gI2);
   }
   tmp_231 += tmp_232;
   result += (-0.5) * tmp_231;
   std::complex<double> tmp_233;
   std::complex<double> tmp_234;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_234 += B1(p,MFu(gI2),MVWp)*Conj(CpbarUFdconjVWpFuPL(gO2,gI2))*
         CpbarUFdconjVWpFuPL(gO1,gI2);
   }
   tmp_233 += tmp_234;
   result += (-1) * tmp_233;
   std::complex<double> tmp_235;
   std::complex<double> tmp_236;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_236 += B1(p,MFd(gI2),Mhh)*Conj(CpbarUFdhhFdPR(gO2,gI2))*
         CpbarUFdhhFdPR(gO1,gI2);
   }
   tmp_235 += tmp_236;
   result += (-0.5) * tmp_235;
   std::complex<double> tmp_237;
   std::complex<double> tmp_238;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_238 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVGFdPL(gO2,gI2))*
         CpbarUFdVGFdPL(gO1,gI2);
   }
   tmp_237 += tmp_238;
   result += (-1.3333333333333333) * tmp_237;
   std::complex<double> tmp_239;
   std::complex<double> tmp_240;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_240 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVPFdPL(gO2,gI2))*
         CpbarUFdVPFdPL(gO1,gI2);
   }
   tmp_239 += tmp_240;
   result += (-1) * tmp_239;
   std::complex<double> tmp_241;
   std::complex<double> tmp_242;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_242 += B1(p,MFd(gI2),MVZ)*Conj(CpbarUFdVZFdPL(gO2,gI2))*
         CpbarUFdVZFdPL(gO1,gI2);
   }
   tmp_241 += tmp_242;
   result += (-1) * tmp_241;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_243;
   std::complex<double> tmp_244;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_244 += B1(p,MFd(gI1),MAh)*Conj(CpbarUFdFdAhPL(gO2,gI1))*
         CpbarUFdFdAhPL(gO1,gI1);
   }
   tmp_243 += tmp_244;
   result += (-0.5) * tmp_243;
   std::complex<double> tmp_245;
   std::complex<double> tmp_246;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_246 += B1(p,MFu(gI2),MHp)*Conj(CpbarUFdconjHpFuPL(gO2,gI2))*
         CpbarUFdconjHpFuPL(gO1,gI2);
   }
   tmp_245 += tmp_246;
   result += (-0.5) * tmp_245;
   std::complex<double> tmp_247;
   std::complex<double> tmp_248;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_248 += B1(p,MFu(gI2),MVWp)*Conj(CpbarUFdconjVWpFuPR(gO2,gI2))*
         CpbarUFdconjVWpFuPR(gO1,gI2);
   }
   tmp_247 += tmp_248;
   result += (-1) * tmp_247;
   std::complex<double> tmp_249;
   std::complex<double> tmp_250;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_250 += B1(p,MFd(gI2),Mhh)*Conj(CpbarUFdhhFdPL(gO2,gI2))*
         CpbarUFdhhFdPL(gO1,gI2);
   }
   tmp_249 += tmp_250;
   result += (-0.5) * tmp_249;
   std::complex<double> tmp_251;
   std::complex<double> tmp_252;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_252 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVGFdPR(gO2,gI2))*
         CpbarUFdVGFdPR(gO1,gI2);
   }
   tmp_251 += tmp_252;
   result += (-1.3333333333333333) * tmp_251;
   std::complex<double> tmp_253;
   std::complex<double> tmp_254;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_254 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVPFdPR(gO2,gI2))*
         CpbarUFdVPFdPR(gO1,gI2);
   }
   tmp_253 += tmp_254;
   result += (-1) * tmp_253;
   std::complex<double> tmp_255;
   std::complex<double> tmp_256;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_256 += B1(p,MFd(gI2),MVZ)*Conj(CpbarUFdVZFdPR(gO2,gI2))*
         CpbarUFdVZFdPR(gO1,gI2);
   }
   tmp_255 += tmp_256;
   result += (-1) * tmp_255;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_257;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_257 += B0(p,MFu(gI1),MAh)*Conj(CpbarUFuFuAhPL(gO2,gI1))*
         CpbarUFuFuAhPR(gO1,gI1)*MFu(gI1);
   }
   result += tmp_257;
   std::complex<double> tmp_258;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_258 += B0(p,MFd(gI2),MHp)*Conj(CpbarUFuHpFdPL(gO2,gI2))*
         CpbarUFuHpFdPR(gO1,gI2)*MFd(gI2);
   }
   result += tmp_258;
   std::complex<double> tmp_259;
   std::complex<double> tmp_260;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_260 += B0(p,MFd(gI2),MVWp)*Conj(CpbarUFuVWpFdPR(gO2,gI2))*
         CpbarUFuVWpFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_259 += tmp_260;
   result += (-4) * tmp_259;
   std::complex<double> tmp_261;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_261 += B0(p,MFu(gI2),Mhh)*Conj(CpbarUFuhhFuPL(gO2,gI2))*
         CpbarUFuhhFuPR(gO1,gI2)*MFu(gI2);
   }
   result += tmp_261;
   std::complex<double> tmp_262;
   std::complex<double> tmp_263;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_263 += B0(p,MFu(gI2),0)*Conj(CpbarUFuVGFuPR(gO2,gI2))*
         CpbarUFuVGFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_262 += tmp_263;
   result += (-5.333333333333333) * tmp_262;
   std::complex<double> tmp_264;
   std::complex<double> tmp_265;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_265 += B0(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_264 += tmp_265;
   result += (-4) * tmp_264;
   std::complex<double> tmp_266;
   std::complex<double> tmp_267;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_267 += B0(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_266 += tmp_267;
   result += (-4) * tmp_266;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_268;
   std::complex<double> tmp_269;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_269 += B1(p,MFu(gI1),MAh)*Conj(CpbarUFuFuAhPR(gO2,gI1))*
         CpbarUFuFuAhPR(gO1,gI1);
   }
   tmp_268 += tmp_269;
   result += (-0.5) * tmp_268;
   std::complex<double> tmp_270;
   std::complex<double> tmp_271;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_271 += B1(p,MFu(gI2),Mhh)*Conj(CpbarUFuhhFuPR(gO2,gI2))*
         CpbarUFuhhFuPR(gO1,gI2);
   }
   tmp_270 += tmp_271;
   result += (-0.5) * tmp_270;
   std::complex<double> tmp_272;
   std::complex<double> tmp_273;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_273 += B1(p,MFd(gI2),MHp)*Conj(CpbarUFuHpFdPR(gO2,gI2))*
         CpbarUFuHpFdPR(gO1,gI2);
   }
   tmp_272 += tmp_273;
   result += (-0.5) * tmp_272;
   std::complex<double> tmp_274;
   std::complex<double> tmp_275;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_275 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVGFuPL(gO2,gI2))*
         CpbarUFuVGFuPL(gO1,gI2);
   }
   tmp_274 += tmp_275;
   result += (-1.3333333333333333) * tmp_274;
   std::complex<double> tmp_276;
   std::complex<double> tmp_277;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_277 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPL(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2);
   }
   tmp_276 += tmp_277;
   result += (-1) * tmp_276;
   std::complex<double> tmp_278;
   std::complex<double> tmp_279;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_279 += B1(p,MFd(gI2),MVWp)*Conj(CpbarUFuVWpFdPL(gO2,gI2))*
         CpbarUFuVWpFdPL(gO1,gI2);
   }
   tmp_278 += tmp_279;
   result += (-1) * tmp_278;
   std::complex<double> tmp_280;
   std::complex<double> tmp_281;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_281 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPL(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2);
   }
   tmp_280 += tmp_281;
   result += (-1) * tmp_280;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_282;
   std::complex<double> tmp_283;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_283 += B1(p,MFu(gI1),MAh)*Conj(CpbarUFuFuAhPL(gO2,gI1))*
         CpbarUFuFuAhPL(gO1,gI1);
   }
   tmp_282 += tmp_283;
   result += (-0.5) * tmp_282;
   std::complex<double> tmp_284;
   std::complex<double> tmp_285;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_285 += B1(p,MFu(gI2),Mhh)*Conj(CpbarUFuhhFuPL(gO2,gI2))*
         CpbarUFuhhFuPL(gO1,gI2);
   }
   tmp_284 += tmp_285;
   result += (-0.5) * tmp_284;
   std::complex<double> tmp_286;
   std::complex<double> tmp_287;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_287 += B1(p,MFd(gI2),MHp)*Conj(CpbarUFuHpFdPL(gO2,gI2))*
         CpbarUFuHpFdPL(gO1,gI2);
   }
   tmp_286 += tmp_287;
   result += (-0.5) * tmp_286;
   std::complex<double> tmp_288;
   std::complex<double> tmp_289;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_289 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVGFuPR(gO2,gI2))*
         CpbarUFuVGFuPR(gO1,gI2);
   }
   tmp_288 += tmp_289;
   result += (-1.3333333333333333) * tmp_288;
   std::complex<double> tmp_290;
   std::complex<double> tmp_291;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_291 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPR(gO1,gI2);
   }
   tmp_290 += tmp_291;
   result += (-1) * tmp_290;
   std::complex<double> tmp_292;
   std::complex<double> tmp_293;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_293 += B1(p,MFd(gI2),MVWp)*Conj(CpbarUFuVWpFdPR(gO2,gI2))*
         CpbarUFuVWpFdPR(gO1,gI2);
   }
   tmp_292 += tmp_293;
   result += (-1) * tmp_292;
   std::complex<double> tmp_294;
   std::complex<double> tmp_295;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_295 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPR(gO1,gI2);
   }
   tmp_294 += tmp_295;
   result += (-1) * tmp_294;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_296;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_296 += B0(p,MFe(gI1),MAh)*Conj(CpbarUFeFeAhPL(gO2,gI1))*
         CpbarUFeFeAhPR(gO1,gI1)*MFe(gI1);
   }
   result += tmp_296;
   std::complex<double> tmp_297;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_297 += B0(p,MFe(gI2),Mhh)*Conj(CpbarUFehhFePL(gO2,gI2))*
         CpbarUFehhFePR(gO1,gI2)*MFe(gI2);
   }
   result += tmp_297;
   std::complex<double> tmp_298;
   std::complex<double> tmp_299;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_299 += B0(p,MFe(gI2),0)*Conj(CpbarUFeVPFePR(gO2,gI2))*
         CpbarUFeVPFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_298 += tmp_299;
   result += (-4) * tmp_298;
   std::complex<double> tmp_300;
   std::complex<double> tmp_301;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_301 += B0(p,MFe(gI2),MVZ)*Conj(CpbarUFeVZFePR(gO2,gI2))*
         CpbarUFeVZFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_300 += tmp_301;
   result += (-4) * tmp_300;
   std::complex<double> tmp_302;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_302 += B0(p,MFv(gI2),MHp)*Conj(CpbarUFeconjHpFvPL(gO2,gI2))*
         CpbarUFeconjHpFvPR(gO1,gI2)*MFv(gI2);
   }
   result += tmp_302;
   std::complex<double> tmp_303;
   std::complex<double> tmp_304;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_304 += B0(p,MFv(gI2),MVWp)*Conj(CpbarUFeconjVWpFvPR(gO2,gI2))*
         CpbarUFeconjVWpFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_303 += tmp_304;
   result += (-4) * tmp_303;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_305;
   std::complex<double> tmp_306;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_306 += B1(p,MFe(gI1),MAh)*Conj(CpbarUFeFeAhPR(gO2,gI1))*
         CpbarUFeFeAhPR(gO1,gI1);
   }
   tmp_305 += tmp_306;
   result += (-0.5) * tmp_305;
   std::complex<double> tmp_307;
   std::complex<double> tmp_308;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_308 += B1(p,MFv(gI2),MHp)*Conj(CpbarUFeconjHpFvPR(gO2,gI2))*
         CpbarUFeconjHpFvPR(gO1,gI2);
   }
   tmp_307 += tmp_308;
   result += (-0.5) * tmp_307;
   std::complex<double> tmp_309;
   std::complex<double> tmp_310;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_310 += B1(p,MFv(gI2),MVWp)*Conj(CpbarUFeconjVWpFvPL(gO2,gI2))*
         CpbarUFeconjVWpFvPL(gO1,gI2);
   }
   tmp_309 += tmp_310;
   result += (-1) * tmp_309;
   std::complex<double> tmp_311;
   std::complex<double> tmp_312;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_312 += B1(p,MFe(gI2),Mhh)*Conj(CpbarUFehhFePR(gO2,gI2))*
         CpbarUFehhFePR(gO1,gI2);
   }
   tmp_311 += tmp_312;
   result += (-0.5) * tmp_311;
   std::complex<double> tmp_313;
   std::complex<double> tmp_314;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_314 += B1(p,MFe(gI2),0)*Conj(CpbarUFeVPFePL(gO2,gI2))*
         CpbarUFeVPFePL(gO1,gI2);
   }
   tmp_313 += tmp_314;
   result += (-1) * tmp_313;
   std::complex<double> tmp_315;
   std::complex<double> tmp_316;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_316 += B1(p,MFe(gI2),MVZ)*Conj(CpbarUFeVZFePL(gO2,gI2))*
         CpbarUFeVZFePL(gO1,gI2);
   }
   tmp_315 += tmp_316;
   result += (-1) * tmp_315;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_317;
   std::complex<double> tmp_318;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_318 += B1(p,MFe(gI1),MAh)*Conj(CpbarUFeFeAhPL(gO2,gI1))*
         CpbarUFeFeAhPL(gO1,gI1);
   }
   tmp_317 += tmp_318;
   result += (-0.5) * tmp_317;
   std::complex<double> tmp_319;
   std::complex<double> tmp_320;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_320 += B1(p,MFv(gI2),MHp)*Conj(CpbarUFeconjHpFvPL(gO2,gI2))*
         CpbarUFeconjHpFvPL(gO1,gI2);
   }
   tmp_319 += tmp_320;
   result += (-0.5) * tmp_319;
   std::complex<double> tmp_321;
   std::complex<double> tmp_322;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_322 += B1(p,MFv(gI2),MVWp)*Conj(CpbarUFeconjVWpFvPR(gO2,gI2))*
         CpbarUFeconjVWpFvPR(gO1,gI2);
   }
   tmp_321 += tmp_322;
   result += (-1) * tmp_321;
   std::complex<double> tmp_323;
   std::complex<double> tmp_324;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_324 += B1(p,MFe(gI2),Mhh)*Conj(CpbarUFehhFePL(gO2,gI2))*
         CpbarUFehhFePL(gO1,gI2);
   }
   tmp_323 += tmp_324;
   result += (-0.5) * tmp_323;
   std::complex<double> tmp_325;
   std::complex<double> tmp_326;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_326 += B1(p,MFe(gI2),0)*Conj(CpbarUFeVPFePR(gO2,gI2))*
         CpbarUFeVPFePR(gO1,gI2);
   }
   tmp_325 += tmp_326;
   result += (-1) * tmp_325;
   std::complex<double> tmp_327;
   std::complex<double> tmp_328;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_328 += B1(p,MFe(gI2),MVZ)*Conj(CpbarUFeVZFePR(gO2,gI2))*
         CpbarUFeVZFePR(gO1,gI2);
   }
   tmp_327 += tmp_328;
   result += (-1) * tmp_327;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VZ_heavy(double p ) const
{
   std::complex<double> result;

   result += -4*AbsSqr(CpVZconjHpHp())*B00(p,MHp,MHp);
   result += 2*AbsSqr(CpVZconjVWpHp())*B0(p,MVWp,MHp);
   result += 0.5*A0(MAh)*CpVZVZAhAh();
   result += A0(MHp)*CpVZVZconjHpHp();
   result += AbsSqr(CpVZVZhh())*B0(p,MVZ,Mhh);
   result += 0.5*A0(Mhh)*CpVZVZhhhh();

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VWp_heavy(double p ) const
{
   std::complex<double> result;

   result += -4*AbsSqr(CpconjVWpHpAh())*B00(p,MAh,MHp);
   result += -4*AbsSqr(CpconjVWpHphh())*B00(p,Mhh,MHp);
   result += AbsSqr(CpconjVWpVPHp())*B0(p,0,MHp);
   result += AbsSqr(CpconjVWpVWphh())*B0(p,MVWp,Mhh);
   result += AbsSqr(CpconjVWpVZHp())*B0(p,MVZ,MHp);
   result += 0.5*A0(MAh)*CpVWpconjVWpAhAh();
   result += A0(MHp)*CpVWpconjVWpconjHpHp();
   result += 0.5*A0(Mhh)*CpVWpconjVWphhhh();

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_329;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_329 += B0(p,MFd(gI1),MAh)*Conj(CpbarFdFdAhPL(gO2,gI1))*
         CpbarFdFdAhPR(gO1,gI1)*MFd(gI1);
   }
   result += tmp_329;
   std::complex<double> tmp_330;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_330 += B0(p,MFd(gI2),Mhh)*Conj(CpbarFdhhFdPL(gO2,gI2))*
         CpbarFdhhFdPR(gO1,gI2)*MFd(gI2);
   }
   result += tmp_330;
   std::complex<double> tmp_331;
   std::complex<double> tmp_332;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_332 += B0(p,MFd(gI2),MVZ)*Conj(CpbarFdVZFdPR(gO2,gI2))*
         CpbarFdVZFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_331 += tmp_332;
   result += (-4) * tmp_331;
   std::complex<double> tmp_333;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_333 += B0(p,MFu(gI2),MHp)*Conj(CpbarFdconjHpFuPL(gO2,gI2))*
         CpbarFdconjHpFuPR(gO1,gI2)*MFu(gI2);
   }
   result += tmp_333;
   std::complex<double> tmp_334;
   std::complex<double> tmp_335;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_335 += B0(p,MFu(gI2),MVWp)*Conj(CpbarFdconjVWpFuPR(gO2,gI2))*
         CpbarFdconjVWpFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_334 += tmp_335;
   result += (-4) * tmp_334;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_336;
   std::complex<double> tmp_337;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_337 += B1(p,MFd(gI1),MAh)*Conj(CpbarFdFdAhPR(gO2,gI1))*
         CpbarFdFdAhPR(gO1,gI1);
   }
   tmp_336 += tmp_337;
   result += (-0.5) * tmp_336;
   std::complex<double> tmp_338;
   std::complex<double> tmp_339;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_339 += B1(p,MFu(gI2),MHp)*Conj(CpbarFdconjHpFuPR(gO2,gI2))*
         CpbarFdconjHpFuPR(gO1,gI2);
   }
   tmp_338 += tmp_339;
   result += (-0.5) * tmp_338;
   std::complex<double> tmp_340;
   std::complex<double> tmp_341;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_341 += B1(p,MFu(gI2),MVWp)*Conj(CpbarFdconjVWpFuPL(gO2,gI2))*
         CpbarFdconjVWpFuPL(gO1,gI2);
   }
   tmp_340 += tmp_341;
   result += (-1) * tmp_340;
   std::complex<double> tmp_342;
   std::complex<double> tmp_343;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_343 += B1(p,MFd(gI2),Mhh)*Conj(CpbarFdhhFdPR(gO2,gI2))*
         CpbarFdhhFdPR(gO1,gI2);
   }
   tmp_342 += tmp_343;
   result += (-0.5) * tmp_342;
   std::complex<double> tmp_344;
   std::complex<double> tmp_345;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_345 += B1(p,MFd(gI2),MVZ)*Conj(CpbarFdVZFdPL(gO2,gI2))*
         CpbarFdVZFdPL(gO1,gI2);
   }
   tmp_344 += tmp_345;
   result += (-1) * tmp_344;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_346;
   std::complex<double> tmp_347;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_347 += B1(p,MFd(gI1),MAh)*Conj(CpbarFdFdAhPL(gO2,gI1))*
         CpbarFdFdAhPL(gO1,gI1);
   }
   tmp_346 += tmp_347;
   result += (-0.5) * tmp_346;
   std::complex<double> tmp_348;
   std::complex<double> tmp_349;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_349 += B1(p,MFu(gI2),MHp)*Conj(CpbarFdconjHpFuPL(gO2,gI2))*
         CpbarFdconjHpFuPL(gO1,gI2);
   }
   tmp_348 += tmp_349;
   result += (-0.5) * tmp_348;
   std::complex<double> tmp_350;
   std::complex<double> tmp_351;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_351 += B1(p,MFu(gI2),MVWp)*Conj(CpbarFdconjVWpFuPR(gO2,gI2))*
         CpbarFdconjVWpFuPR(gO1,gI2);
   }
   tmp_350 += tmp_351;
   result += (-1) * tmp_350;
   std::complex<double> tmp_352;
   std::complex<double> tmp_353;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_353 += B1(p,MFd(gI2),Mhh)*Conj(CpbarFdhhFdPL(gO2,gI2))*
         CpbarFdhhFdPL(gO1,gI2);
   }
   tmp_352 += tmp_353;
   result += (-0.5) * tmp_352;
   std::complex<double> tmp_354;
   std::complex<double> tmp_355;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_355 += B1(p,MFd(gI2),MVZ)*Conj(CpbarFdVZFdPR(gO2,gI2))*
         CpbarFdVZFdPR(gO1,gI2);
   }
   tmp_354 += tmp_355;
   result += (-1) * tmp_354;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_356;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_356 += B0(p,MFe(gI1),MAh)*Conj(CpbarFeFeAhPL(gO2,gI1))*
         CpbarFeFeAhPR(gO1,gI1)*MFe(gI1);
   }
   result += tmp_356;
   std::complex<double> tmp_357;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_357 += B0(p,MFe(gI2),Mhh)*Conj(CpbarFehhFePL(gO2,gI2))*
         CpbarFehhFePR(gO1,gI2)*MFe(gI2);
   }
   result += tmp_357;
   std::complex<double> tmp_358;
   std::complex<double> tmp_359;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_359 += B0(p,MFe(gI2),MVZ)*Conj(CpbarFeVZFePR(gO2,gI2))*
         CpbarFeVZFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_358 += tmp_359;
   result += (-4) * tmp_358;
   std::complex<double> tmp_360;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_360 += B0(p,MFv(gI2),MHp)*Conj(CpbarFeconjHpFvPL(gO2,gI2))*
         CpbarFeconjHpFvPR(gO1,gI2)*MFv(gI2);
   }
   result += tmp_360;
   std::complex<double> tmp_361;
   std::complex<double> tmp_362;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_362 += B0(p,MFv(gI2),MVWp)*Conj(CpbarFeconjVWpFvPR(gO2,gI2))*
         CpbarFeconjVWpFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_361 += tmp_362;
   result += (-4) * tmp_361;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_363;
   std::complex<double> tmp_364;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_364 += B1(p,MFe(gI1),MAh)*Conj(CpbarFeFeAhPR(gO2,gI1))*
         CpbarFeFeAhPR(gO1,gI1);
   }
   tmp_363 += tmp_364;
   result += (-0.5) * tmp_363;
   std::complex<double> tmp_365;
   std::complex<double> tmp_366;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_366 += B1(p,MFv(gI2),MHp)*Conj(CpbarFeconjHpFvPR(gO2,gI2))*
         CpbarFeconjHpFvPR(gO1,gI2);
   }
   tmp_365 += tmp_366;
   result += (-0.5) * tmp_365;
   std::complex<double> tmp_367;
   std::complex<double> tmp_368;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_368 += B1(p,MFv(gI2),MVWp)*Conj(CpbarFeconjVWpFvPL(gO2,gI2))*
         CpbarFeconjVWpFvPL(gO1,gI2);
   }
   tmp_367 += tmp_368;
   result += (-1) * tmp_367;
   std::complex<double> tmp_369;
   std::complex<double> tmp_370;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_370 += B1(p,MFe(gI2),Mhh)*Conj(CpbarFehhFePR(gO2,gI2))*
         CpbarFehhFePR(gO1,gI2);
   }
   tmp_369 += tmp_370;
   result += (-0.5) * tmp_369;
   std::complex<double> tmp_371;
   std::complex<double> tmp_372;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_372 += B1(p,MFe(gI2),MVZ)*Conj(CpbarFeVZFePL(gO2,gI2))*
         CpbarFeVZFePL(gO1,gI2);
   }
   tmp_371 += tmp_372;
   result += (-1) * tmp_371;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_373;
   std::complex<double> tmp_374;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_374 += B1(p,MFe(gI1),MAh)*Conj(CpbarFeFeAhPL(gO2,gI1))*
         CpbarFeFeAhPL(gO1,gI1);
   }
   tmp_373 += tmp_374;
   result += (-0.5) * tmp_373;
   std::complex<double> tmp_375;
   std::complex<double> tmp_376;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_376 += B1(p,MFv(gI2),MHp)*Conj(CpbarFeconjHpFvPL(gO2,gI2))*
         CpbarFeconjHpFvPL(gO1,gI2);
   }
   tmp_375 += tmp_376;
   result += (-0.5) * tmp_375;
   std::complex<double> tmp_377;
   std::complex<double> tmp_378;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_378 += B1(p,MFv(gI2),MVWp)*Conj(CpbarFeconjVWpFvPR(gO2,gI2))*
         CpbarFeconjVWpFvPR(gO1,gI2);
   }
   tmp_377 += tmp_378;
   result += (-1) * tmp_377;
   std::complex<double> tmp_379;
   std::complex<double> tmp_380;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_380 += B1(p,MFe(gI2),Mhh)*Conj(CpbarFehhFePL(gO2,gI2))*
         CpbarFehhFePL(gO1,gI2);
   }
   tmp_379 += tmp_380;
   result += (-0.5) * tmp_379;
   std::complex<double> tmp_381;
   std::complex<double> tmp_382;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_382 += B1(p,MFe(gI2),MVZ)*Conj(CpbarFeVZFePR(gO2,gI2))*
         CpbarFeVZFePR(gO1,gI2);
   }
   tmp_381 += tmp_382;
   result += (-1) * tmp_381;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_383;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_383 += B0(p,MFu(gI1),MAh)*Conj(CpbarFuFuAhPL(gO2,gI1))*
         CpbarFuFuAhPR(gO1,gI1)*MFu(gI1);
   }
   result += tmp_383;
   std::complex<double> tmp_384;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_384 += B0(p,MFd(gI2),MHp)*Conj(CpbarFuHpFdPL(gO2,gI2))*
         CpbarFuHpFdPR(gO1,gI2)*MFd(gI2);
   }
   result += tmp_384;
   std::complex<double> tmp_385;
   std::complex<double> tmp_386;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_386 += B0(p,MFd(gI2),MVWp)*Conj(CpbarFuVWpFdPR(gO2,gI2))*
         CpbarFuVWpFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_385 += tmp_386;
   result += (-4) * tmp_385;
   std::complex<double> tmp_387;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_387 += B0(p,MFu(gI2),Mhh)*Conj(CpbarFuhhFuPL(gO2,gI2))*
         CpbarFuhhFuPR(gO1,gI2)*MFu(gI2);
   }
   result += tmp_387;
   std::complex<double> tmp_388;
   std::complex<double> tmp_389;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_389 += B0(p,MFu(gI2),0)*Conj(CpbarFuVPFuPR(gO2,gI2))*CpbarFuVPFuPL
         (gO1,gI2)*MFu(gI2);
   }
   tmp_388 += tmp_389;
   result += (-4) * tmp_388;
   std::complex<double> tmp_390;
   std::complex<double> tmp_391;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_391 += B0(p,MFu(gI2),MVZ)*Conj(CpbarFuVZFuPR(gO2,gI2))*
         CpbarFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_390 += tmp_391;
   result += (-4) * tmp_390;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_392;
   std::complex<double> tmp_393;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_393 += B1(p,MFu(gI1),MAh)*Conj(CpbarFuFuAhPR(gO2,gI1))*
         CpbarFuFuAhPR(gO1,gI1);
   }
   tmp_392 += tmp_393;
   result += (-0.5) * tmp_392;
   std::complex<double> tmp_394;
   std::complex<double> tmp_395;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_395 += B1(p,MFu(gI2),Mhh)*Conj(CpbarFuhhFuPR(gO2,gI2))*
         CpbarFuhhFuPR(gO1,gI2);
   }
   tmp_394 += tmp_395;
   result += (-0.5) * tmp_394;
   std::complex<double> tmp_396;
   std::complex<double> tmp_397;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_397 += B1(p,MFd(gI2),MHp)*Conj(CpbarFuHpFdPR(gO2,gI2))*
         CpbarFuHpFdPR(gO1,gI2);
   }
   tmp_396 += tmp_397;
   result += (-0.5) * tmp_396;
   std::complex<double> tmp_398;
   std::complex<double> tmp_399;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_399 += B1(p,MFu(gI2),0)*Conj(CpbarFuVPFuPL(gO2,gI2))*CpbarFuVPFuPL
         (gO1,gI2);
   }
   tmp_398 += tmp_399;
   result += (-1) * tmp_398;
   std::complex<double> tmp_400;
   std::complex<double> tmp_401;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_401 += B1(p,MFd(gI2),MVWp)*Conj(CpbarFuVWpFdPL(gO2,gI2))*
         CpbarFuVWpFdPL(gO1,gI2);
   }
   tmp_400 += tmp_401;
   result += (-1) * tmp_400;
   std::complex<double> tmp_402;
   std::complex<double> tmp_403;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_403 += B1(p,MFu(gI2),MVZ)*Conj(CpbarFuVZFuPL(gO2,gI2))*
         CpbarFuVZFuPL(gO1,gI2);
   }
   tmp_402 += tmp_403;
   result += (-1) * tmp_402;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_404;
   std::complex<double> tmp_405;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_405 += B1(p,MFu(gI1),MAh)*Conj(CpbarFuFuAhPL(gO2,gI1))*
         CpbarFuFuAhPL(gO1,gI1);
   }
   tmp_404 += tmp_405;
   result += (-0.5) * tmp_404;
   std::complex<double> tmp_406;
   std::complex<double> tmp_407;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_407 += B1(p,MFu(gI2),Mhh)*Conj(CpbarFuhhFuPL(gO2,gI2))*
         CpbarFuhhFuPL(gO1,gI2);
   }
   tmp_406 += tmp_407;
   result += (-0.5) * tmp_406;
   std::complex<double> tmp_408;
   std::complex<double> tmp_409;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_409 += B1(p,MFd(gI2),MHp)*Conj(CpbarFuHpFdPL(gO2,gI2))*
         CpbarFuHpFdPL(gO1,gI2);
   }
   tmp_408 += tmp_409;
   result += (-0.5) * tmp_408;
   std::complex<double> tmp_410;
   std::complex<double> tmp_411;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_411 += B1(p,MFu(gI2),0)*Conj(CpbarFuVPFuPR(gO2,gI2))*CpbarFuVPFuPR
         (gO1,gI2);
   }
   tmp_410 += tmp_411;
   result += (-1) * tmp_410;
   std::complex<double> tmp_412;
   std::complex<double> tmp_413;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_413 += B1(p,MFd(gI2),MVWp)*Conj(CpbarFuVWpFdPR(gO2,gI2))*
         CpbarFuVWpFdPR(gO1,gI2);
   }
   tmp_412 += tmp_413;
   result += (-1) * tmp_412;
   std::complex<double> tmp_414;
   std::complex<double> tmp_415;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_415 += B1(p,MFu(gI2),MVZ)*Conj(CpbarFuVZFuPR(gO2,gI2))*
         CpbarFuVZFuPR(gO1,gI2);
   }
   tmp_414 += tmp_415;
   result += (-1) * tmp_414;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::tadpole_hh() const
{
   std::complex<double> result;

   result += -0.5*A0(MAh)*CphhAhAh();
   result += A0(MVWp)*CphhbargWpCgWpC();
   result += A0(MVWp)*CphhbargWpgWp();
   result += A0(MVZ)*CphhbargZgZ();
   result += -(A0(MHp)*CphhconjHpHp());
   result += 4*A0(MVWp)*CphhconjVWpVWp();
   result += -0.5*A0(Mhh)*Cphhhhhh();
   result += 2*A0(MVZ)*CphhVZVZ();
   std::complex<double> tmp_416;
   std::complex<double> tmp_417;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_417 += A0(MFd(gI1))*(CphhbarFdFdPL(gI1,gI1) + CphhbarFdFdPR(gI1,
         gI1))*MFd(gI1);
   }
   tmp_416 += tmp_417;
   result += (6) * tmp_416;
   std::complex<double> tmp_418;
   std::complex<double> tmp_419;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_419 += A0(MFe(gI1))*(CphhbarFeFePL(gI1,gI1) + CphhbarFeFePR(gI1,
         gI1))*MFe(gI1);
   }
   tmp_418 += tmp_419;
   result += (2) * tmp_418;
   std::complex<double> tmp_420;
   std::complex<double> tmp_421;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_421 += A0(MFu(gI1))*(CphhbarFuFuPL(gI1,gI1) + CphhbarFuFuPR(gI1,
         gI1))*MFu(gI1);
   }
   tmp_420 += tmp_421;
   result += (6) * tmp_420;

   return result * oneOver16PiSqr;

}








void CLASSNAME::calculate_MVG_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MVG) = 0.;
}

void CLASSNAME::calculate_MHp_pole()
{
   if (problems.is_tachyon(Hp))
      return;
   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MHp) old_MHp(MHp), new_MHp(MHp);

   do {
      const double p = old_MHp;
      const double self_energy = Re(self_energy_Hp(p));
      const double mass_sqr = Sqr(MHp) - self_energy;

      if (mass_sqr < 0.)
         problems.flag_tachyon(Hp);

      PHYSICAL(MHp) = AbsSqrt(mass_sqr);

      new_MHp = PHYSICAL(MHp);
      diff = MaxRelDiff(new_MHp, old_MHp);
      old_MHp = new_MHp;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);
}

void CLASSNAME::calculate_MFv_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MFv).setConstant(0.);
}

void CLASSNAME::calculate_MAh_pole()
{
   if (problems.is_tachyon(Ah))
      return;
   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MAh) old_MAh(MAh), new_MAh(MAh);

   do {
      const double p = old_MAh;
      const double self_energy = Re(self_energy_Ah(p));
      const double mass_sqr = Sqr(MAh) - self_energy;

      if (mass_sqr < 0.)
         problems.flag_tachyon(Ah);

      PHYSICAL(MAh) = AbsSqrt(mass_sqr);

      new_MAh = PHYSICAL(MAh);
      diff = MaxRelDiff(new_MAh, old_MAh);
      old_MAh = new_MAh;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);
}

void CLASSNAME::calculate_Mhh_pole()
{
   if (problems.is_tachyon(hh))
      return;
   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(Mhh) old_Mhh(Mhh), new_Mhh(Mhh);

   do {
      const double p = old_Mhh;
      const double self_energy = Re(self_energy_hh(p));
      const double mass_sqr = Sqr(Mhh) - self_energy;

      if (mass_sqr < 0.)
         problems.flag_tachyon(hh);

      PHYSICAL(Mhh) = AbsSqrt(mass_sqr);

      new_Mhh = PHYSICAL(Mhh);
      diff = MaxRelDiff(new_Mhh, old_Mhh);
      old_Mhh = new_Mhh;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);
}

void CLASSNAME::calculate_MVP_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MVP) = 0.;
}

void CLASSNAME::calculate_MVZ_pole()
{
   if (problems.is_tachyon(VZ))
      return;
   // diagonalization with medium precision
   const double p = MVZ;
   const double self_energy = Re(self_energy_VZ(p));
   const double mass_sqr = Sqr(MVZ) - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(VZ);

   PHYSICAL(MVZ) = AbsSqrt(mass_sqr);
}

void CLASSNAME::calculate_MFd_pole()
{
   // diagonalization with medium precision
   Eigen::Matrix<double,3,3> self_energy_1;
   Eigen::Matrix<double,3,3> self_energy_PL;
   Eigen::Matrix<double,3,3> self_energy_PR;
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fd());
   for (unsigned es = 0; es < 3; ++es) {
      const double p = Abs(MFd(es));
      for (unsigned i1 = 0; i1 < 3; ++i1) {
         for (unsigned i2 = 0; i2 < 3; ++i2) {
            self_energy_1(i1,i2)  = Re(self_energy_Fd_1(p,i1,i2)
               );
            self_energy_PL(i1,i2) = Re(self_energy_Fd_PL(p,i1,i2
               ));
            self_energy_PR(i1,i2) = Re(self_energy_Fd_PR(p,i1,i2
               ));
         }
      }
      const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,3,3> M_1loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(Vd) mix_Vd;
      decltype(Ud) mix_Ud;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_1loop, eigen_values, mix_Vd, mix_Ud, eigenvalue_error);
      problems.flag_bad_mass(StandardModel_info::Fd, eigenvalue_error
         > precision * Abs(eigen_values(0)));
   #else
      fs_svd(M_1loop, eigen_values, mix_Vd, mix_Ud);
   #endif
      if (es == 0) {
         PHYSICAL(Vd) = mix_Vd;
         PHYSICAL(Ud) = mix_Ud;
      }
      PHYSICAL(MFd(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MFu_pole()
{
   // diagonalization with medium precision
   Eigen::Matrix<double,3,3> self_energy_1;
   Eigen::Matrix<double,3,3> self_energy_PL;
   Eigen::Matrix<double,3,3> self_energy_PR;
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fu());
   for (unsigned es = 0; es < 3; ++es) {
      const double p = Abs(MFu(es));
      for (unsigned i1 = 0; i1 < 3; ++i1) {
         for (unsigned i2 = 0; i2 < 3; ++i2) {
            self_energy_1(i1,i2)  = Re(self_energy_Fu_1(p,i1,i2)
               );
            self_energy_PL(i1,i2) = Re(self_energy_Fu_PL(p,i1,i2
               ));
            self_energy_PR(i1,i2) = Re(self_energy_Fu_PR(p,i1,i2
               ));
         }
      }
      const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,3,3> M_1loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(Vu) mix_Vu;
      decltype(Uu) mix_Uu;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_1loop, eigen_values, mix_Vu, mix_Uu, eigenvalue_error);
      problems.flag_bad_mass(StandardModel_info::Fu, eigenvalue_error
         > precision * Abs(eigen_values(0)));
   #else
      fs_svd(M_1loop, eigen_values, mix_Vu, mix_Uu);
   #endif
      if (es == 0) {
         PHYSICAL(Vu) = mix_Vu;
         PHYSICAL(Uu) = mix_Uu;
      }
      PHYSICAL(MFu(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MFe_pole()
{
   // diagonalization with medium precision
   Eigen::Matrix<double,3,3> self_energy_1;
   Eigen::Matrix<double,3,3> self_energy_PL;
   Eigen::Matrix<double,3,3> self_energy_PR;
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fe());
   for (unsigned es = 0; es < 3; ++es) {
      const double p = Abs(MFe(es));
      for (unsigned i1 = 0; i1 < 3; ++i1) {
         for (unsigned i2 = 0; i2 < 3; ++i2) {
            self_energy_1(i1,i2)  = Re(self_energy_Fe_1(p,i1,i2)
               );
            self_energy_PL(i1,i2) = Re(self_energy_Fe_PL(p,i1,i2
               ));
            self_energy_PR(i1,i2) = Re(self_energy_Fe_PR(p,i1,i2
               ));
         }
      }
      const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,3,3> M_1loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(Ve) mix_Ve;
      decltype(Ue) mix_Ue;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_1loop, eigen_values, mix_Ve, mix_Ue, eigenvalue_error);
      problems.flag_bad_mass(StandardModel_info::Fe, eigenvalue_error
         > precision * Abs(eigen_values(0)));
   #else
      fs_svd(M_1loop, eigen_values, mix_Ve, mix_Ue);
   #endif
      if (es == 0) {
         PHYSICAL(Ve) = mix_Ve;
         PHYSICAL(Ue) = mix_Ue;
      }
      PHYSICAL(MFe(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MVWp_pole()
{
   if (problems.is_tachyon(VWp))
      return;
   // diagonalization with medium precision
   const double p = MVWp;
   const double self_energy = Re(self_energy_VWp(p));
   const double mass_sqr = Sqr(MVWp) - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(VWp);

   PHYSICAL(MVWp) = AbsSqrt(mass_sqr);
}


double CLASSNAME::calculate_MFu_DRbar(double m_pole, int idx) const
{
   const double p = m_pole;
   const double self_energy_1  = Re(self_energy_Fu_1_heavy_rotated(p, idx
      , idx));
   const double self_energy_PL = Re(self_energy_Fu_PL_heavy_rotated(p,
      idx, idx));
   const double self_energy_PR = Re(self_energy_Fu_PR_heavy_rotated(p,
      idx, idx));

   const double currentScale = get_scale();
   const double qcd_1l = -0.008443431970194815*(5. - 3.*Log(Sqr(MFu(2))
      /Sqr(currentScale)))*Sqr(g3);
   const double qcd_2l = -0.003408916029785599*Power(g3,4) +
      0.0011495761378943394*Power(g3,4)*Log(Sqr(MFu(2))/Sqr(currentScale)) -
      0.00024060895909416413*Power(g3,4)*Sqr(Log(Power(MFu(2),2)/Sqr(
      currentScale)));

   const double m_susy_drbar = m_pole + self_energy_1 + m_pole * (
      self_energy_PL + self_energy_PR + qcd_1l + qcd_2l);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MFd_DRbar(double m_sm_msbar, int idx) const
{
   const double p = m_sm_msbar;
   const double self_energy_1  = Re(self_energy_Fd_1_heavy_rotated(p, idx
      , idx));
   const double self_energy_PL = Re(self_energy_Fd_PL_heavy_rotated(p,
      idx, idx));
   const double self_energy_PR = Re(self_energy_Fd_PR_heavy_rotated(p,
      idx, idx));
   const double m_tree = MFd(2);
   const double m_sm_drbar = m_sm_msbar * (1 - 0.00020496318737651018*
      Power(g3,4) + 0.0006860288475783287*Sqr(g1) + 0.0023747152416172916*Sqr(
      g2) - 0.008443431970194815*Sqr(g3));

   const double m_susy_drbar = m_sm_drbar / (1.0 - self_energy_1/m_tree -
      self_energy_PL - self_energy_PR);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MFe_DRbar(double m_sm_msbar, int idx) const
{
   const double p = m_sm_msbar;
   const double self_energy_1  = Re(self_energy_Fe_1_heavy_rotated(p, idx
      , idx));
   const double self_energy_PL = Re(self_energy_Fe_PL_heavy_rotated(p,
      idx, idx));
   const double self_energy_PR = Re(self_energy_Fe_PR_heavy_rotated(p,
      idx, idx));
   const double m_sm_drbar = m_sm_msbar * (1 - 0.0023747152416172916*(0.6
      *Sqr(g1) - Sqr(g2)));

   const double m_susy_drbar = m_sm_drbar + self_energy_1 + m_sm_drbar *
      (self_energy_PL + self_energy_PR);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MFv_DRbar(double, int) const
{
   return 0.0;
}

double CLASSNAME::calculate_MVP_DRbar(double)
{
   return 0.0;
}

double CLASSNAME::calculate_MVZ_DRbar(double m_pole)
{
   const double p = m_pole;
   const double self_energy = Re(self_energy_VZ(p));
   const double mass_sqr = Sqr(m_pole) + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_tachyon(VZ);
      return m_pole;
   }

   return AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVWp_DRbar(double m_pole)
{
   const double p = m_pole;
   const double self_energy = Re(self_energy_VWp(p));
   const double mass_sqr = Sqr(m_pole) + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_tachyon(VWp);
      return m_pole;
   }

   return AbsSqrt(mass_sqr);
}


double CLASSNAME::ThetaW() const
{
   return ArcTan((0.7745966692414834*g1)/g2);
}


std::ostream& operator<<(std::ostream& ostr, const StandardModel<Two_scale>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
