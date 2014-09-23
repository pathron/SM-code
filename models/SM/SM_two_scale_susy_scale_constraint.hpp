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

// File generated at Fri 22 Aug 2014 11:56:06

#ifndef SM_TWO_SCALE_SUSY_SCALE_CONSTRAINT_H
#define SM_TWO_SCALE_SUSY_SCALE_CONSTRAINT_H

#include "SM_susy_scale_constraint.hpp"
#include "SM_input_parameters.hpp"
#include "two_scale_constraint.hpp"

namespace flexiblesusy {

template <class T>
class SM;

class Two_scale;

template<>
class SM_susy_scale_constraint<Two_scale> : public Constraint<Two_scale> {
public:
   SM_susy_scale_constraint();
   SM_susy_scale_constraint(const SM_input_parameters&);
   virtual ~SM_susy_scale_constraint();
   virtual void apply();
   virtual double get_scale() const;
   virtual void set_model(Two_scale_model*);

   void clear();
   double get_initial_scale_guess() const;
   void initialize();
   void set_input_parameters(const SM_input_parameters&);

private:
   double scale;
   double initial_scale_guess;
   SM<Two_scale>* model;
   SM_input_parameters inputPars;

   void update_scale();
};

} // namespace flexiblesusy

#endif
