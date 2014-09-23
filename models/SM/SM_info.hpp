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

#ifndef SM_INFO_H
#define SM_INFO_H

#include <iosfwd>

namespace flexiblesusy {

namespace SM_info {
   enum Particles : unsigned {Hp, Fv, Ah, hh, VZ, Fd, Fu, Fe, VG, VP, VWp,
      NUMBER_OF_PARTICLES};

   enum Parameters : unsigned {Lambdax, Yu00, Yu01, Yu02, Yu10, Yu11, Yu12,
      Yu20, Yu21, Yu22, Yd00, Yd01, Yd02, Yd10, Yd11, Yd12, Yd20, Yd21, Yd22, Ye00
      , Ye01, Ye02, Ye10, Ye11, Ye12, Ye20, Ye21, Ye22, g1, g2, g3, v, mu2,
      NUMBER_OF_PARAMETERS};

   extern const unsigned particle_multiplicities[NUMBER_OF_PARTICLES];
   extern const char* particle_names[NUMBER_OF_PARTICLES];
   extern const char* particle_latex_names[NUMBER_OF_PARTICLES];
   extern const char* parameter_names[NUMBER_OF_PARAMETERS];
   extern const char* model_name;
   extern const bool is_low_energy_model;

   void print(std::ostream&);
}

} // namespace flexiblesusy

#endif
