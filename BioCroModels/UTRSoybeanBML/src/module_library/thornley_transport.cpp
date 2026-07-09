#include <vector>
#include <string>
#include "../framework/module.h"                  // for direct_module and update
#include "../framework/module_helper_functions.h"  // for get_ip
#include "thornley_nutrient_dynamics.h"  // for the thornley_nutrient_dynamics namespace
#include "thornley_transport.h"

using thornley_nutrient_dynamics::generate_pairwise_names;
using thornley_nutrient_dynamics::generate_quantity_names_from_transport_links;
using thornley_nutrient_dynamics::get_transport_link_ips;
using thornley_nutrient_dynamics::get_transport_link_ops;
using thornley_nutrient_dynamics::organ;
using thornley_nutrient_dynamics::transport_link;

thornley_transport::thornley_transport(
    std::vector<transport_link> const& organ_links,
    state_map const& input_quantities,
    state_map* output_quantities)
    :  // Define basic module properties by passing its name to its parent class
      differential_module(),

      // Store the organ links
      organ_links(organ_links),

      // Get pointers to input parameters
      substrate_transport_ips(get_ip(input_quantities, generate_pairwise_names("substrate_transport", organ_links))),

      // Get pointers to output parameters
      substrate_carbon_ops(get_transport_link_ops(output_quantities, organ_links, "substrate_carbon"))
{
}

std::vector<std::string> thornley_transport::get_inputs(std::vector<transport_link> const& organ_links)
{
    // List the quantity names that exist for each organ in each link
    std::vector<std::string> quantities_for_each_organ = {
        "" // Mg / ha
    };

    // Append the organ names as prefixes
    std::vector<std::string> inputs = generate_quantity_names_from_transport_links(organ_links, quantities_for_each_organ);
    
    // Add the substrate transports
    std::vector<std::string> substrate_transport_names = generate_pairwise_names("substrate_transport", organ_links); // mol / m2

    // Append substrate transports
    inputs.insert(inputs.end(), substrate_transport_names.begin(), substrate_transport_names.end());

    return inputs;  
}

std::vector<std::string> thornley_transport::get_outputs(std::vector<transport_link> const& organ_links)
{
    // List the quantity names that exist for each organ in each link
    std::vector<std::string> quantities_for_each_organ = {
        "substrate_carbon"  // mol / m^2
    };

    // Append the organ names as prefixes
    std::vector<std::string> outputs = generate_quantity_names_from_transport_links(organ_links, quantities_for_each_organ);

    return outputs;
}

void thornley_transport::do_multi_organ_operation() const
{
    // Substrate carbon is transferred pairwise between organ links.
    // Note that we are relying on the additive property of differential_module::update
    for (size_t i = 0; i < organ_links.size(); ++i) {
        // A positive value for the transport rate means that substrate flows
        // from the first organ to the second
        update(substrate_carbon_ops[i].first, -*substrate_transport_ips[i]);  // mol / m^2
        update(substrate_carbon_ops[i].second, *substrate_transport_ips[i]);  // mol / m^2
    }
}
