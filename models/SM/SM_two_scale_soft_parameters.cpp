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

// File generated at Fri 22 Aug 2014 11:56:05

#include "SM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

#include <iostream>

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

SM_soft_parameters::SM_soft_parameters(const SM_input_parameters& input_)
   : SM_susy_parameters(input_)
   , mu2(0)

{
   set_number_of_parameters(numberOfParameters);
}

SM_soft_parameters::SM_soft_parameters(
   const SM_susy_parameters& susy_model
   , double mu2_

)
   : SM_susy_parameters(susy_model)
   , mu2(mu2_)

{
   set_number_of_parameters(numberOfParameters);
}

Eigen::ArrayXd SM_soft_parameters::beta() const
{
   return calc_beta().get();
}

SM_soft_parameters SM_soft_parameters::calc_beta() const
{
   Soft_traces soft_traces;
   calc_soft_traces(soft_traces);

   double beta_mu2(calc_beta_mu2_one_loop(TRACE_STRUCT));

   if (get_loops() > 1) {
      beta_mu2 += calc_beta_mu2_two_loop(TRACE_STRUCT);

   }


   const SM_susy_parameters susy_betas(SM_susy_parameters::calc_beta());

   return SM_soft_parameters(susy_betas, beta_mu2);
}

void SM_soft_parameters::clear()
{
   SM_susy_parameters::clear();

   mu2 = 0.;

}

const Eigen::ArrayXd SM_soft_parameters::get() const
{
   Eigen::ArrayXd pars(SM_susy_parameters::get());
   pars.conservativeResize(numberOfParameters);

   pars(32) = mu2;


   return pars;
}

void SM_soft_parameters::print(std::ostream& ostr) const
{
   SM_susy_parameters::print(ostr);
   ostr << "soft parameters:\n";
   ostr << "mu2 = " << mu2 << '\n';

}

void SM_soft_parameters::set(const Eigen::ArrayXd& pars)
{
   SM_susy_parameters::set(pars);

   mu2 = pars(32);

}

void SM_soft_parameters::calc_soft_traces(Soft_traces& soft_traces) const
{
   TRACE_STRUCT.traceYdAdjYd = (Yd*Yd.adjoint()).trace();
   TRACE_STRUCT.traceYeAdjYe = (Ye*Ye.adjoint()).trace();
   TRACE_STRUCT.traceYuAdjYu = (Yu*Yu.adjoint()).trace();
   TRACE_STRUCT.traceYdAdjYdYdAdjYd = (Yd*Yd.adjoint()*Yd*Yd.adjoint()).trace()
      ;
   TRACE_STRUCT.traceYdAdjYuYuAdjYd = (Yd*Yu.adjoint()*Yu*Yd.adjoint()).trace()
      ;
   TRACE_STRUCT.traceYeAdjYeYeAdjYe = (Ye*Ye.adjoint()*Ye*Ye.adjoint()).trace()
      ;
   TRACE_STRUCT.traceYuAdjYuYuAdjYu = (Yu*Yu.adjoint()*Yu*Yu.adjoint()).trace()
      ;

}

std::ostream& operator<<(std::ostream& ostr, const SM_soft_parameters& soft_pars)
{
   soft_pars.print(std::cout);
   return ostr;
}

} // namespace flexiblesusy
