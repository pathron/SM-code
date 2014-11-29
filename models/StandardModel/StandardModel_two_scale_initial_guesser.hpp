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

// File generated at Sat 29 Nov 2014 17:55:09

#ifndef StandardModel_TWO_SCALE_INITIAL_GUESSER_H
#define StandardModel_TWO_SCALE_INITIAL_GUESSER_H

#include "StandardModel_initial_guesser.hpp"
#include "StandardModel_input_parameters.hpp"
#include "StandardModel_two_scale_low_scale_constraint.hpp"
#include "StandardModel_two_scale_susy_scale_constraint.hpp"
#include "two_scale_initial_guesser.hpp"
#include "error.hpp"
#include "lowe.h"

#include <sstream>

namespace flexiblesusy {

template <class T>
class StandardModel;

class Two_scale;

/**
 * @class StandardModel_initial_guesser<Two_scale>
 * @brief initial guesser for the StandardModel
 */

template<>
class StandardModel_initial_guesser<Two_scale> : public Initial_guesser<Two_scale> {
public:
   StandardModel_initial_guesser(StandardModel<Two_scale>*,
                               const StandardModel_input_parameters&,
                               const QedQcd&,
                               const StandardModel_low_scale_constraint<Two_scale>&,
                               const StandardModel_susy_scale_constraint<Two_scale>&);
   virtual ~StandardModel_initial_guesser();
   virtual void guess(); ///< initial guess

   void set_running_precision(double p) { running_precision = p; }

private:
   StandardModel<Two_scale>* model; ///< pointer to model class
   StandardModel_input_parameters input_pars;
   QedQcd oneset;   ///< Standard Model low-energy data
   double mu_guess; ///< guessed DR-bar mass of up-quark
   double mc_guess; ///< guessed DR-bar mass of charm-quark
   double mt_guess; ///< guessed DR-bar mass of top-quark
   double md_guess; ///< guessed DR-bar mass of down-quark
   double ms_guess; ///< guessed DR-bar mass of strange-quark
   double mb_guess; ///< guessed DR-bar mass of bottom-quark
   double me_guess; ///< guessed DR-bar mass of electron
   double mm_guess; ///< guessed DR-bar mass of muon
   double mtau_guess; ///< guessed DR-bar mass of tau
   double running_precision; ///< Runge-Kutta RG running precision
   StandardModel_low_scale_constraint<Two_scale> low_constraint;
   StandardModel_susy_scale_constraint<Two_scale> susy_constraint;

   void guess_susy_parameters();
   void guess_soft_parameters();
   void calculate_DRbar_yukawa_couplings();
   void calculate_Yu_DRbar();
   void calculate_Yd_DRbar();
   void calculate_Ye_DRbar();
};

} // namespace flexiblesusy

#endif
