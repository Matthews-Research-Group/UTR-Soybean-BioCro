#include <vector>
#include <string>
#include <unordered_map>
#include "../framework/module.h"                  // for differential_module and update
#include "../framework/module_helper_functions.h"  // for get_ip
#include "../framework/constants.h"                // for getting physical constants
#include "thornley_nutrient_dynamics.h"  // for the thornley_nutrient_dynamics namespace
#include "thornley_biomass_calculator.h"
// #include <Rcpp.h> // for debugging

// using namespace Rcpp; // for debugging
using thornley_nutrient_dynamics::generate_multi_organ_quantity_names;
using thornley_nutrient_dynamics::get_multi_organ_ips;
using thornley_nutrient_dynamics::get_multi_organ_ops;
using thornley_nutrient_dynamics::organ;

thornley_biomass_calculator::thornley_biomass_calculator(
    std::vector<organ> const& organs,
    state_map const& input_quantities,
    state_map* output_quantities)
    :  // Define basic module properties by passing its name to its parent class
      direct_module(),

      // Store the organ names
      organs(organs),

      // Get pointers to organ input quantities
      substrate_carbon_ips(get_multi_organ_ips(input_quantities, organs, "substrate_carbon")),
      structural_carbon_ips(get_multi_organ_ips(input_quantities, organs, "structural_carbon")),

      // Get pointers to input parameters
      carbon_to_mass_factor_ips(get_multi_organ_ips(input_quantities, organs, "carbon_to_mass_factor")),

      // Get pointers to output parameters
      mass_ops(get_multi_organ_ops(output_quantities, organs, ""))
{
}

std::vector<std::string> thornley_biomass_calculator::get_inputs(
    std::vector<organ> const& organs)
{
    // Add the external substrate sources
    std::vector<std::string> inputs = get_external_substrate_quantity_names(organs);  // mol / m^2

    // List the quantity names that are guaranteed to exist for each organ
    std::vector<std::string> quantities_for_each_organ = {
        "substrate_carbon",         // mol / m^2
        "structural_carbon",        // mol / m^2
        "carbon_to_mass_factor"     // Mg/ha / (mol/m^2)
    };

    return inputs;
}

std::vector<std::string> thornley_biomass_calculator::get_outputs(std::vector<organ> const& organs)
{
    // List the quantity names that exist for each organ
    std::vector<std::string> quantities_for_each_organ = {
        ""    // Mg / ha
    };
    // Append the organ names as prefixes
    return generate_multi_organ_quantity_names(organs, quantities_for_each_organ);
}

void thornley_biomass_calculator::do_multi_organ_operation() const
{
    for (size_t i = 0; i < organs.size(); ++i) {
        update(mass_ops[i], *carbon_to_mass_factor_ips[i] * (*substrate_carbon_ips[i] + *structural_carbon_ips[i]));
    }
}
