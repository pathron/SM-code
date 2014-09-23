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

#include "SM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of v.
 *
 * @return one-loop beta function
 */
double SM_susy_parameters::calc_beta_v_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_v;

   beta_v = -(oneOver16PiSqr*(3*traceYdAdjYd + traceYeAdjYe + 3*
      traceYuAdjYu)*v);


   return beta_v;
}

/**
 * Calculates the two-loop beta function of v.
 *
 * @return two-loop beta function
 */
double SM_susy_parameters::calc_beta_v_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_v;

   beta_v = -0.00125*twoLoop*v*(1257*Power(g1,4) - 8575*Power(g2,4) -
      5400*traceYdAdjYdYdAdjYd + 1200*traceYdAdjYuYuAdjYd - 1800*
      traceYeAdjYeYeAdjYe - 5400*traceYuAdjYuYuAdjYu + 2060*traceYuAdjYu*Sqr(g1
      ) + 6300*traceYuAdjYu*Sqr(g2) - 90*Sqr(g1)*Sqr(g2) + 60*traceYeAdjYe*(27*
      Sqr(g1) + 35*Sqr(g2)) + 16000*traceYuAdjYu*Sqr(g3) + 20*traceYdAdjYd*(43*
      Sqr(g1) + 315*Sqr(g2) + 800*Sqr(g3)) + 1200*Sqr(Lambdax));


   return beta_v;
}

} // namespace flexiblesusy
