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

#ifndef StandardModel_SLHA_IO_H
#define StandardModel_SLHA_IO_H

#include "StandardModel_two_scale_model_slha.hpp"
#include "StandardModel_info.hpp"
#include "StandardModel_physical.hpp"
#include "slha_io.hpp"
#include "ew_input.hpp"

#include <Eigen/Core>
#include <string>
#include <utility>

#define Pole(p) physical.p
#define PHYSICAL(p) model.get_physical().p
#define LOCALPHYSICAL(p) physical.p
#define MODELPARAMETER(p) model.get_##p()
#define DEFINE_PARAMETER(p)                                            \
   typename std::remove_const<typename std::remove_reference<decltype(MODELPARAMETER(p))>::type>::type p;
#define DEFINE_POLE_MASS(p)                                            \
   typename std::remove_const<typename std::remove_reference<decltype(PHYSICAL(p))>::type>::type p;
#define SM(p) Electroweak_constants::p

namespace flexiblesusy {

struct StandardModel_input_parameters;
class Spectrum_generator_settings;

class StandardModel_slha_io {
public:
   StandardModel_slha_io();
   ~StandardModel_slha_io() {}

   void clear();

   void fill(QedQcd& qedqcd) const { slha_io.fill(qedqcd); }
   void fill(StandardModel_input_parameters&) const;
   template <class T> void fill(StandardModel_slha<T>&) const;
   void fill(Spectrum_generator_settings&) const;
   double get_input_scale() const;
   double get_parameter_output_scale() const;
   const SLHA_io& get_slha_io() const { return slha_io; }
   void read_from_file(const std::string&);
   void set_extpar(const StandardModel_input_parameters&);
   template <class T> void set_extra(const StandardModel_slha<T>&);
   void set_minpar(const StandardModel_input_parameters&);
   void set_sminputs(const softsusy::QedQcd&);
   template <class T> void set_spectrum(const StandardModel_slha<T>&);
   template <class T> void set_spectrum(const StandardModel<T>&);
   void set_spinfo(const Problems<StandardModel_info::NUMBER_OF_PARTICLES>&);
   void write_to_file(const std::string&);
   void write_to_stream(std::ostream& ostr = std::cout) { slha_io.write_to_stream(ostr); }

   static void fill_minpar_tuple(StandardModel_input_parameters&, int, double);
   static void fill_extpar_tuple(StandardModel_input_parameters&, int, double);
   static void fill_flexiblesusy_tuple(Spectrum_generator_settings&, int, double);

   static void convert_to_slha_convention(StandardModel_physical&);

   template <class T>
   static void fill_slhaea(SLHAea::Coll&, const StandardModel_slha<T>&, const QedQcd& qedqcd);

   template <class T>
   static SLHAea::Coll fill_slhaea(const StandardModel_slha<T>&, const QedQcd& qedqcd);

private:
   SLHA_io slha_io; ///< SLHA io class
   static unsigned const NUMBER_OF_DRBAR_BLOCKS = 6;
   static char const * const drbar_blocks[NUMBER_OF_DRBAR_BLOCKS];

   void set_mass(const StandardModel_physical&, bool);
   void set_mixing_matrices(const StandardModel_physical&, bool);
   template <class T> void set_model_parameters(const StandardModel_slha<T>&);
   double read_scale() const;
   template <class T> void fill_drbar_parameters(StandardModel_slha<T>&) const;
   template <class T> void fill_physical(StandardModel_slha<T>&) const;
};

/**
 * Reads DR-bar parameters, pole masses and mixing matrices from a
 * SLHA output file.
 */
template <class T>
void StandardModel_slha_io::fill(StandardModel_slha<T>& model) const
{
   fill_drbar_parameters(model);
   fill_physical(model);
}

/**
 * Reads DR-bar parameters from a SLHA output file.
 */
template <class T>
void StandardModel_slha_io::fill_drbar_parameters(StandardModel_slha<T>& model) const
{
   model.set_g1(slha_io.read_entry("gauge", 1) * 1.2909944487358056);
   model.set_g2(slha_io.read_entry("gauge", 2));
   model.set_g3(slha_io.read_entry("gauge", 3));
   {
      DEFINE_PARAMETER(Yu);
      slha_io.read_block("Yu", Yu);
      model.set_Yu(Yu);
   }
   {
      DEFINE_PARAMETER(Yd);
      slha_io.read_block("Yd", Yd);
      model.set_Yd(Yd);
   }
   {
      DEFINE_PARAMETER(Ye);
      slha_io.read_block("Ye", Ye);
      model.set_Ye(Ye);
   }
   model.set_mu2(slha_io.read_entry("SM", 1));
   model.set_Lambdax(slha_io.read_entry("SM", 2));
   model.set_v(slha_io.read_entry("HMIX", 3));


   model.set_scale(read_scale());
}

/**
 * Reads pole masses and mixing matrices from a SLHA output file.
 */
template <class T>
void StandardModel_slha_io::fill_physical(StandardModel_slha<T>& model) const
{
   {
      DEFINE_PARAMETER(Vu);
      slha_io.read_block("UULMIX", Vu);
      PHYSICAL(Vu) = Vu;
   }
   {
      DEFINE_PARAMETER(Vd);
      slha_io.read_block("UDLMIX", Vd);
      PHYSICAL(Vd) = Vd;
   }
   {
      DEFINE_PARAMETER(Uu);
      slha_io.read_block("UURMIX", Uu);
      PHYSICAL(Uu) = Uu;
   }
   {
      DEFINE_PARAMETER(Ud);
      slha_io.read_block("UDRMIX", Ud);
      PHYSICAL(Ud) = Ud;
   }
   {
      DEFINE_PARAMETER(Ve);
      slha_io.read_block("UELMIX", Ve);
      PHYSICAL(Ve) = Ve;
   }
   {
      DEFINE_PARAMETER(Ue);
      slha_io.read_block("UERMIX", Ue);
      PHYSICAL(Ue) = Ue;
   }

   PHYSICAL(MVG) = slha_io.read_entry("MASS", 21);
   PHYSICAL(MFv)(0) = slha_io.read_entry("MASS", 12);
   PHYSICAL(MFv)(1) = slha_io.read_entry("MASS", 14);
   PHYSICAL(MFv)(2) = slha_io.read_entry("MASS", 16);
   PHYSICAL(Mhh) = slha_io.read_entry("MASS", 25);
   PHYSICAL(MVP) = slha_io.read_entry("MASS", 22);
   PHYSICAL(MVZ) = slha_io.read_entry("MASS", 23);
   PHYSICAL(MFd)(0) = slha_io.read_entry("MASS", 1);
   PHYSICAL(MFd)(1) = slha_io.read_entry("MASS", 3);
   PHYSICAL(MFd)(2) = slha_io.read_entry("MASS", 5);
   PHYSICAL(MFu)(0) = slha_io.read_entry("MASS", 2);
   PHYSICAL(MFu)(1) = slha_io.read_entry("MASS", 4);
   PHYSICAL(MFu)(2) = slha_io.read_entry("MASS", 6);
   PHYSICAL(MFe)(0) = slha_io.read_entry("MASS", 11);
   PHYSICAL(MFe)(1) = slha_io.read_entry("MASS", 13);
   PHYSICAL(MFe)(2) = slha_io.read_entry("MASS", 15);
   PHYSICAL(MVWp) = slha_io.read_entry("MASS", 24);

}

template <class T>
void StandardModel_slha_io::fill_slhaea(SLHAea::Coll& slhaea, const StandardModel_slha<T>& model, const QedQcd& qedqcd)
{
   StandardModel_slha_io slha_io;
   const StandardModel_input_parameters& input = model.get_input();
   const Problems<StandardModel_info::NUMBER_OF_PARTICLES>& problems
      = model.get_problems();
   const bool error = problems.have_problem();

   slha_io.set_spinfo(problems);
   slha_io.set_sminputs(qedqcd);
   slha_io.set_minpar(input);
   slha_io.set_extpar(input);
   if (!error) {
      slha_io.set_spectrum(model);
      slha_io.set_extra(model);
   }

   slhaea = slha_io.get_slha_io().get_data();
}

template <class T>
SLHAea::Coll StandardModel_slha_io::fill_slhaea(const StandardModel_slha<T>& model, const QedQcd& qedqcd)
{
   SLHAea::Coll slhaea;
   StandardModel_slha_io::fill_slhaea(slhaea, model, qedqcd);

   return slhaea;
}

/**
 * Stores the model (DR-bar) parameters in the SLHA object.
 *
 * @param model model class
 */
template <class T>
void StandardModel_slha_io::set_model_parameters(const StandardModel_slha<T>& model)
{
   {
      std::ostringstream block;
      block << "Block gauge Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(g1) * 0.7745966692414834), "gY")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(g2)), "g2")
            << FORMAT_ELEMENT(3, (MODELPARAMETER(g3)), "g3")
      ;
      slha_io.set_block(block);
   }
   slha_io.set_block("Yu", MODELPARAMETER(Yu), "Yu", model.get_scale());
   slha_io.set_block("Yd", MODELPARAMETER(Yd), "Yd", model.get_scale());
   slha_io.set_block("Ye", MODELPARAMETER(Ye), "Ye", model.get_scale());
   {
      std::ostringstream block;
      block << "Block SM Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(1, (MODELPARAMETER(mu2)), "mu2")
            << FORMAT_ELEMENT(2, (MODELPARAMETER(Lambdax)), "Lambdax")
      ;
      slha_io.set_block(block);
   }
   {
      std::ostringstream block;
      block << "Block HMIX Q= " << FORMAT_SCALE(model.get_scale()) << '\n'
            << FORMAT_ELEMENT(3, (MODELPARAMETER(v)), "v")
      ;
      slha_io.set_block(block);
   }

}

/**
 * Writes extra SLHA blocks
 *
 * @param model model class
 */
template <class T>
void StandardModel_slha_io::set_extra(const StandardModel_slha<T>& model)
{
   const StandardModel_physical physical(model.get_physical_slha());


}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in BPMZ convention
 */
template <class T>
void StandardModel_slha_io::set_spectrum(const StandardModel<T>& model)
{
   const StandardModel_slha<T> model_slha(model);
   set_spectrum(model_slha);
}

/**
 * Stores the model (DR-bar) parameters, masses and mixing matrices in
 * the SLHA object.
 *
 * @param model model class in SLHA convention
 */
template <class T>
void StandardModel_slha_io::set_spectrum(const StandardModel_slha<T>& model)
{
   const StandardModel_physical physical(model.get_physical_slha());
   const bool write_sm_masses = model.do_calculate_sm_pole_masses();

   set_model_parameters(model);
   set_mass(physical, write_sm_masses);
   set_mixing_matrices(physical, write_sm_masses);
}

} // namespace flexiblesusy

#endif
