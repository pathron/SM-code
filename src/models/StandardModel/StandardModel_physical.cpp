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

// File generated at Fri 28 Nov 2014 22:36:38

#include "StandardModel_physical.hpp"

#include <iostream>

namespace flexiblesusy {

StandardModel_physical::StandardModel_physical()
   :
    MVG(0), MHp(0), MFv(Eigen::Array<double,3,1>::Zero()), MAh(0), Mhh(0), MVP
       (0), MVZ(0), MFd(Eigen::Array<double,3,1>::Zero()), MFu(Eigen::Array<double
       ,3,1>::Zero()), MFe(Eigen::Array<double,3,1>::Zero()), MVWp(0)

   , Vd(Eigen::Matrix<std::complex<double>,3,3>::Zero()), Ud(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), Vu(Eigen::Matrix<std::complex<double>,3,
      3>::Zero()), Uu(Eigen::Matrix<std::complex<double>,3,3>::Zero()), Ve(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), Ue(Eigen::Matrix<
      std::complex<double>,3,3>::Zero())

{
}

void StandardModel_physical::clear()
{
   MVG = 0.;
   MHp = 0.;
   MFv = Eigen::Matrix<double,3,1>::Zero();
   MAh = 0.;
   Mhh = 0.;
   MVP = 0.;
   MVZ = 0.;
   MFd = Eigen::Matrix<double,3,1>::Zero();
   Vd = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ud = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFu = Eigen::Matrix<double,3,1>::Zero();
   Vu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Uu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFe = Eigen::Matrix<double,3,1>::Zero();
   Ve = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ue = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MVWp = 0.;

}

void StandardModel_physical::print(std::ostream& ostr) const
{
   ostr << "----------------------------------------\n"
           "pole masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MHp = " << MHp << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MAh = " << MAh << '\n';
   ostr << "Mhh = " << Mhh << '\n';
   ostr << "MVP = " << MVP << '\n';
   ostr << "MVZ = " << MVZ << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MVWp = " << MVWp << '\n';

   ostr << "----------------------------------------\n"
           "pole mass mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "Vd = " << Vd << '\n';
   ostr << "Ud = " << Ud << '\n';
   ostr << "Vu = " << Vu << '\n';
   ostr << "Uu = " << Uu << '\n';
   ostr << "Ve = " << Ve << '\n';
   ostr << "Ue = " << Ue << '\n';

}

std::ostream& operator<<(std::ostream& ostr, const StandardModel_physical& phys_pars)
{
   phys_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
