#include <vector>
#include <string>
#include <cmath>                         // for pow
#include "../framework/module.h"                  // for direct_module and update
#include "../framework/module_helper_functions.h"  // for get_ip, get_input
#include "thornley_nutrient_dynamics.h"  // for the thornley_nutrient_dynamics namespace
#include "thornley_transport_calculator.h"

using thornley_nutrient_dynamics::generate_quantity_names_from_transport_links;
using thornley_nutrient_dynamics::get_multi_organ_ips;
using thornley_nutrient_dynamics::get_transport_link_ips;
using thornley_nutrient_dynamics::transport_link;
using thornley_nutrient_dynamics::organ_library;

thornley_transport_calculator::thornley_transport_calculator(
    std::vector<transport_link> const& organ_links,
    state_map const& input_quantities,
    state_map* output_quantities)
    :  // Define basic module properties by passing its name to its parent class
      direct_module(),

      // Store the organ links
      organ_links(organ_links),

      // Get pointers to input parameters
      substrate_carbon_ips(get_transport_link_ips(input_quantities, organ_links, "substrate_carbon")),
      structural_carbon_ips(get_transport_link_ips(input_quantities, organ_links, "structural_carbon")),
      
      // Get references to input parameters
      substrate_conductance_ips(get_ip(input_quantities, generate_pairwise_names("substrate_conductance", organ_links))),
      transportation_gamma_exponent(get_input(input_quantities, "transportation_gamma_exponent")),
      Pod_start_dvi(get_input(input_quantities, "Pod_start_dvi")),
      stop_growth_dvi(get_input(input_quantities, "stop_growth_dvi")),
      DVI(get_input(input_quantities, "DVI")), // Q: any difference from {}?

      // Get pointers to output parameters
      substrate_transport_ops(get_op(output_quantities, generate_pairwise_names("substrate_transport", organ_links)))

{   
}

std::vector<std::string> thornley_transport_calculator::get_inputs(std::vector<transport_link> const& organ_links)
{
    // List the quantity names that exist for each organ in each link
    std::vector<std::string> quantities_for_each_organ = {
        "",                  // Mg / ha
        "substrate_carbon",  // mol / m^2
        "structural_carbon"  // mol / m^2
    };

    // Append the organ names as prefixes
    std::vector<std::string> inputs = generate_quantity_names_from_transport_links(organ_links, quantities_for_each_organ);

    // Add the substrate conductances [Mg / hr / [Mg / ha]^beta]
    std::vector<std::string> substrate_conductance_names = generate_pairwise_names("substrate_conductance", organ_links);   
    inputs.insert(inputs.end(), substrate_conductance_names.begin(), substrate_conductance_names.end());

    // Add other inputs that don't depend on the organ_links
    inputs.push_back("transportation_gamma_exponent");  // dimensionless 
    inputs.push_back("Pod_start_dvi");
    inputs.push_back("stop_growth_dvi");
    inputs.push_back("DVI");
    return inputs;
}

std::vector<std::string> thornley_transport_calculator::get_outputs(std::vector<transport_link> const& organ_links)
{
    // Get the substrate transport rate names
    return generate_pairwise_names("substrate_transport", organ_links);  // mol / m^2 / hr
}

void thornley_transport_calculator::do_multi_organ_operation() const
{

    // Calculate transport rates between organs and update the relevant outputs
    for (size_t i = 0; i < organ_links.size(); ++i) {
        double const pairwise_mass = std::min(*structural_carbon_ips[i].first, *structural_carbon_ips[i].second);  ;  // mol / m^2
        double const beta_factor = pow(pairwise_mass, transportation_gamma_exponent);   // [mol/m^2]^beta
        double const substrate_gradient = (*substrate_carbon_ips[i].first / *structural_carbon_ips[i].first -
                                          *substrate_carbon_ips[i].second / *structural_carbon_ips[i].second);   // [dimensionless]
        double transport_rate = beta_factor * *substrate_conductance_ips[i] * substrate_gradient;  // mol / m^2 / hr
        
        if ((organ_links[i].second.name() == "Pod" && (DVI < Pod_start_dvi)) || (DVI>stop_growth_dvi)){
            transport_rate = 0;
        }

        update(substrate_transport_ops[i], transport_rate);
    }
}
