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

// File generated at Sat 29 Nov 2014 17:55:07

#include "StandardModel_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

#include <iostream>

namespace flexiblesusy {

#define CLASSNAME StandardModel_susy_parameters
#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

StandardModel_susy_parameters::StandardModel_susy_parameters(const StandardModel_input_parameters& input_)
   : Beta_function()
   , Lambdax(0), Yu(Eigen::Matrix<double,3,3>::Zero()), Yd(Eigen::Matrix<double
   ,3,3>::Zero()), Ye(Eigen::Matrix<double,3,3>::Zero()), g1(0), g2(0), g3(0),
   v(0)

   , input(input_)
{
   set_number_of_parameters(numberOfParameters);
}

StandardModel_susy_parameters::StandardModel_susy_parameters(
   double scale_, double loops_, double thresholds_,
   const StandardModel_input_parameters& input_
   , double Lambdax_, const Eigen::Matrix<double,3,3>& Yu_, const Eigen::Matrix
   <double,3,3>& Yd_, const Eigen::Matrix<double,3,3>& Ye_, double g1_, double
   g2_, double g3_, double v_

)
   : Beta_function()
   , Lambdax(Lambdax_), Yu(Yu_), Yd(Yd_), Ye(Ye_), g1(g1_), g2(g2_), g3(g3_), v
   (v_)

   , input(input_)
{
   set_number_of_parameters(numberOfParameters);
   set_scale(scale_);
   set_loops(loops_);
   set_thresholds(thresholds_);
}

Eigen::ArrayXd StandardModel_susy_parameters::beta() const
{
   return calc_beta().get();
}

StandardModel_susy_parameters StandardModel_susy_parameters::calc_beta() const
{
   Susy_traces susy_traces;
   calc_susy_traces(susy_traces);

   double beta_Lambdax(calc_beta_Lambdax_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_Yu(calc_beta_Yu_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_Yd(calc_beta_Yd_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_Ye(calc_beta_Ye_one_loop(TRACE_STRUCT));
   double beta_g1(calc_beta_g1_one_loop(TRACE_STRUCT));
   double beta_g2(calc_beta_g2_one_loop(TRACE_STRUCT));
   double beta_g3(calc_beta_g3_one_loop(TRACE_STRUCT));
   double beta_v(calc_beta_v_one_loop(TRACE_STRUCT));

   if (get_loops() > 1) {
      beta_Lambdax += calc_beta_Lambdax_two_loop(TRACE_STRUCT);
      beta_Yu += calc_beta_Yu_two_loop(TRACE_STRUCT);
      beta_Yd += calc_beta_Yd_two_loop(TRACE_STRUCT);
      beta_Ye += calc_beta_Ye_two_loop(TRACE_STRUCT);
      beta_g1 += calc_beta_g1_two_loop(TRACE_STRUCT);
      beta_g2 += calc_beta_g2_two_loop(TRACE_STRUCT);
      beta_g3 += calc_beta_g3_two_loop(TRACE_STRUCT);
      beta_v += calc_beta_v_two_loop(TRACE_STRUCT);

   }


   return StandardModel_susy_parameters(get_scale(), get_loops(), get_thresholds(), input,
                    beta_Lambdax, beta_Yu, beta_Yd, beta_Ye, beta_g1, beta_g2, beta_g3, beta_v);
}

void StandardModel_susy_parameters::clear()
{
   reset();
   Lambdax = 0.;
   Yu = Eigen::Matrix<double,3,3>::Zero();
   Yd = Eigen::Matrix<double,3,3>::Zero();
   Ye = Eigen::Matrix<double,3,3>::Zero();
   g1 = 0.;
   g2 = 0.;
   g3 = 0.;
   v = 0.;

}



const Eigen::ArrayXd StandardModel_susy_parameters::get() const
{
   Eigen::ArrayXd pars(numberOfParameters);

   pars(0) = Lambdax;
   pars(1) = Yu(0,0);
   pars(2) = Yu(0,1);
   pars(3) = Yu(0,2);
   pars(4) = Yu(1,0);
   pars(5) = Yu(1,1);
   pars(6) = Yu(1,2);
   pars(7) = Yu(2,0);
   pars(8) = Yu(2,1);
   pars(9) = Yu(2,2);
   pars(10) = Yd(0,0);
   pars(11) = Yd(0,1);
   pars(12) = Yd(0,2);
   pars(13) = Yd(1,0);
   pars(14) = Yd(1,1);
   pars(15) = Yd(1,2);
   pars(16) = Yd(2,0);
   pars(17) = Yd(2,1);
   pars(18) = Yd(2,2);
   pars(19) = Ye(0,0);
   pars(20) = Ye(0,1);
   pars(21) = Ye(0,2);
   pars(22) = Ye(1,0);
   pars(23) = Ye(1,1);
   pars(24) = Ye(1,2);
   pars(25) = Ye(2,0);
   pars(26) = Ye(2,1);
   pars(27) = Ye(2,2);
   pars(28) = g1;
   pars(29) = g2;
   pars(30) = g3;
   pars(31) = v;


   return pars;
}

void StandardModel_susy_parameters::print(std::ostream& ostr) const
{
   ostr << "susy parameters:\n";
   ostr << "Lambdax = " << Lambdax << '\n';
   ostr << "Yu = " << Yu << '\n';
   ostr << "Yd = " << Yd << '\n';
   ostr << "Ye = " << Ye << '\n';
   ostr << "g1 = " << g1 << '\n';
   ostr << "g2 = " << g2 << '\n';
   ostr << "g3 = " << g3 << '\n';
   ostr << "v = " << v << '\n';

}

void StandardModel_susy_parameters::set(const Eigen::ArrayXd& pars)
{
   Lambdax = pars(0);
   Yu(0,0) = pars(1);
   Yu(0,1) = pars(2);
   Yu(0,2) = pars(3);
   Yu(1,0) = pars(4);
   Yu(1,1) = pars(5);
   Yu(1,2) = pars(6);
   Yu(2,0) = pars(7);
   Yu(2,1) = pars(8);
   Yu(2,2) = pars(9);
   Yd(0,0) = pars(10);
   Yd(0,1) = pars(11);
   Yd(0,2) = pars(12);
   Yd(1,0) = pars(13);
   Yd(1,1) = pars(14);
   Yd(1,2) = pars(15);
   Yd(2,0) = pars(16);
   Yd(2,1) = pars(17);
   Yd(2,2) = pars(18);
   Ye(0,0) = pars(19);
   Ye(0,1) = pars(20);
   Ye(0,2) = pars(21);
   Ye(1,0) = pars(22);
   Ye(1,1) = pars(23);
   Ye(1,2) = pars(24);
   Ye(2,0) = pars(25);
   Ye(2,1) = pars(26);
   Ye(2,2) = pars(27);
   g1 = pars(28);
   g2 = pars(29);
   g3 = pars(30);
   v = pars(31);

}

const StandardModel_input_parameters& StandardModel_susy_parameters::get_input() const
{
   return input;
}

void StandardModel_susy_parameters::set_input_parameters(const StandardModel_input_parameters& input_)
{
   input = input_;
}

void StandardModel_susy_parameters::calc_susy_traces(Susy_traces& susy_traces) const
{
   TRACE_STRUCT.traceYdAdjYd = (Yd*Yd.adjoint()).trace();
   TRACE_STRUCT.traceYeAdjYe = (Ye*Ye.adjoint()).trace();
   TRACE_STRUCT.traceYuAdjYu = (Yu*Yu.adjoint()).trace();
   TRACE_STRUCT.traceYdAdjYdYdAdjYd = (Yd*Yd.adjoint()*Yd*Yd.adjoint()).trace()
      ;
   TRACE_STRUCT.traceYeAdjYeYeAdjYe = (Ye*Ye.adjoint()*Ye*Ye.adjoint()).trace()
      ;
   TRACE_STRUCT.traceYuAdjYuYuAdjYu = (Yu*Yu.adjoint()*Yu*Yu.adjoint()).trace()
      ;
   TRACE_STRUCT.traceYdAdjYuYuAdjYd = (Yd*Yu.adjoint()*Yu*Yd.adjoint()).trace()
      ;
   TRACE_STRUCT.traceYdAdjYdYdAdjYdYdAdjYd = (Yd*Yd.adjoint()*Yd*Yd.adjoint()*
      Yd*Yd.adjoint()).trace();
   TRACE_STRUCT.traceYdAdjYdYdAdjYuYuAdjYd = (Yd*Yd.adjoint()*Yd*Yu.adjoint()*
      Yu*Yd.adjoint()).trace();
   TRACE_STRUCT.traceYdAdjYuYuAdjYdYdAdjYd = (Yd*Yu.adjoint()*Yu*Yd.adjoint()*
      Yd*Yd.adjoint()).trace();
   TRACE_STRUCT.traceYdAdjYuYuAdjYuYuAdjYd = (Yd*Yu.adjoint()*Yu*Yu.adjoint()*
      Yu*Yd.adjoint()).trace();
   TRACE_STRUCT.traceYeAdjYeYeAdjYeYeAdjYe = (Ye*Ye.adjoint()*Ye*Ye.adjoint()*
      Ye*Ye.adjoint()).trace();
   TRACE_STRUCT.traceYuAdjYuYuAdjYuYuAdjYu = (Yu*Yu.adjoint()*Yu*Yu.adjoint()*
      Yu*Yu.adjoint()).trace();

}

std::ostream& operator<<(std::ostream& ostr, const StandardModel_susy_parameters& susy_pars)
{
   susy_pars.print(std::cout);
   return ostr;
}

} // namespace flexiblesusy
