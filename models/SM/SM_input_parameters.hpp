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

#ifndef SM_INPUT_PARAMETERS_H
#define SM_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct SM_input_parameters {
   double LambdaIN;
   double LambdaxInput;
   Eigen::Matrix<double,3,3> YuInput;
   Eigen::Matrix<double,3,3> YdInput;
   Eigen::Matrix<double,3,3> YeInput;
   double vInput;

   SM_input_parameters()
      : LambdaIN(0), LambdaxInput(0), YuInput(Eigen::Matrix<double,3,3>::Zero()),
   YdInput(Eigen::Matrix<double,3,3>::Zero()), YeInput(Eigen::Matrix<double,3,3
   >::Zero()), vInput(0)

   {}
};

} // namespace flexiblesusy

#endif
