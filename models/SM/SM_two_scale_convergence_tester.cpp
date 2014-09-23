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

#include "SM_two_scale_convergence_tester.hpp"
#include <cmath>
#include <algorithm>
#include "wrappers.hpp"

namespace flexiblesusy {

#define OLD1(p) ol.get_##p()
#define NEW1(p) ne.get_##p()

#define OLD(p,i) ol.get_##p()(i)
#define NEW(p,i) ne.get_##p()(i)

SM_convergence_tester<Two_scale>::SM_convergence_tester(SM<Two_scale>* model, double accuracy_goal)
   : Convergence_tester_DRbar<SM<Two_scale> >(model, accuracy_goal)
{
}

SM_convergence_tester<Two_scale>::~SM_convergence_tester()
{
}

double SM_convergence_tester<Two_scale>::max_rel_diff() const
{
   const SM<Two_scale>& ol = get_last_iteration_model();
   const SM<Two_scale>& ne = get_model();

   double diff[3] = { 0 };

   diff[0] = MaxRelDiff(OLD1(MHp),NEW1(MHp));
   diff[1] = MaxRelDiff(OLD1(MAh),NEW1(MAh));
   diff[2] = MaxRelDiff(OLD1(Mhh),NEW1(Mhh));

   return *std::max_element(diff, diff + 3);

}

} // namespace flexiblesusy
