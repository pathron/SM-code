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

// File generated at Fri 28 Nov 2014 22:36:33

#include "StandardModel_input_parameters.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

std::ostream& operator<<(std::ostream& ostr, const StandardModel_input_parameters& input)
{
   ostr << "LambdaIN = " << INPUT(LambdaIN) << ", ";
   ostr << "LambdaxInput = " << INPUT(LambdaxInput) << ", ";
   ostr << "YuInput = " << INPUT(YuInput) << ", ";
   ostr << "YdInput = " << INPUT(YdInput) << ", ";
   ostr << "YeInput = " << INPUT(YeInput) << ", ";
   ostr << "vInput = " << INPUT(vInput) << ", ";

   return ostr;
}

} // namespace flexiblesusy