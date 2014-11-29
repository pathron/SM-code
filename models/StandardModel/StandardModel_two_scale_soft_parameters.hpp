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

// File generated at Fri 28 Nov 2014 22:36:32

#ifndef StandardModel_TWO_SCALE_soft_parameters_H
#define StandardModel_TWO_SCALE_soft_parameters_H

#include "rge.h"
#include "StandardModel_two_scale_susy_parameters.hpp"

#include <iosfwd>

namespace flexiblesusy {

#ifdef TRACE_STRUCT_TYPE
   #undef TRACE_STRUCT_TYPE
#endif
#define TRACE_STRUCT_TYPE Soft_traces

class StandardModel_soft_parameters : public StandardModel_susy_parameters {
public:
   explicit StandardModel_soft_parameters(const StandardModel_input_parameters& input_ = StandardModel_input_parameters());
   StandardModel_soft_parameters(const StandardModel_susy_parameters& , double mu2_
);
   virtual ~StandardModel_soft_parameters() {}
   virtual Eigen::ArrayXd beta() const;
   virtual const Eigen::ArrayXd get() const;
   virtual void print(std::ostream&) const;
   virtual void set(const Eigen::ArrayXd&);

   StandardModel_soft_parameters calc_beta() const;
   virtual void clear();

   void set_mu2(double mu2_) { mu2 = mu2_; }

   double get_mu2() const { return mu2; }


protected:
   double mu2;


private:
   static const int numberOfParameters = 33;

   struct Soft_traces {
      double traceYdAdjYd;
      double traceYeAdjYe;
      double traceYuAdjYu;
      double traceYdAdjYdYdAdjYd;
      double traceYdAdjYuYuAdjYd;
      double traceYeAdjYeYeAdjYe;
      double traceYuAdjYuYuAdjYu;

   };
   void calc_soft_traces(Soft_traces&) const;

   double calc_beta_mu2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mu2_two_loop(const TRACE_STRUCT_TYPE&) const;

};

std::ostream& operator<<(std::ostream&, const StandardModel_soft_parameters&);

} // namespace flexiblesusy

#endif
