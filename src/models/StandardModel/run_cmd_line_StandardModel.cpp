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

#include "StandardModel_input_parameters.hpp"
#include "StandardModel_spectrum_generator.hpp"
#include "StandardModel_slha_io.hpp"

#include "command_line_options.hpp"
#include "lowe.h"
#include "logger.hpp"

#include <iostream>

namespace flexiblesusy {

void print_usage()
{
   std::cout <<
      "Usage: run_cmd_line_StandardModel.x [options]\n"
      "Options:\n"
      "  --LambdaIN=<value>\n"

      "  --help,-h                         print this help message"
             << std::endl;
}

void set_command_line_parameters(int argc, char* argv[],
                                 StandardModel_input_parameters& input)
{
   for (int i = 1; i < argc; ++i) {
      const std::string option(argv[i]);

      if(Command_line_options::get_parameter_value(option, "--LambdaIN=", input.LambdaIN))
         continue;

      
      if (option == "--help" || option == "-h") {
         print_usage();
         exit(EXIT_SUCCESS);
      }

      ERROR("Unrecognized command line option: " << option);
      exit(EXIT_FAILURE);
   }
}

} // namespace flexiblesusy


int main(int argc, char* argv[])
{
   using namespace flexiblesusy;
   using namespace softsusy;
   typedef Two_scale algorithm_type;

   StandardModel_input_parameters input;
   set_command_line_parameters(argc, argv, input);

   QedQcd oneset;
   oneset.toMz();

   StandardModel_spectrum_generator<algorithm_type> spectrum_generator;
   spectrum_generator.set_precision_goal(1.0e-4);
   spectrum_generator.set_max_iterations(0);         // 0 == automatic
   spectrum_generator.set_calculate_sm_masses(0);    // 0 == no
   spectrum_generator.set_parameter_output_scale(0); // 0 == susy scale
   spectrum_generator.set_pole_mass_loop_order(2);   // 2-loop
   spectrum_generator.set_ewsb_loop_order(2);        // 2-loop
   spectrum_generator.set_beta_loop_order(2);        // 2-loop
   spectrum_generator.set_threshold_corrections_loop_order(1); // 1-loop

   spectrum_generator.run(oneset, input);

   const int exit_code = spectrum_generator.get_exit_code();
   const StandardModel_slha<algorithm_type> model(spectrum_generator.get_model());

   // SLHA output
   SLHAea::Coll slhaea(StandardModel_slha_io::fill_slhaea(model, oneset));

   std::cout << slhaea;

   return exit_code;
}
