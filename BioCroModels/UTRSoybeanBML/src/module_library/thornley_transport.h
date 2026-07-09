#ifndef THORNLEY_TRANSPORT_H
#define THORNLEY_TRANSPORT_H

#include <vector>
#include <string>
#include "../framework/module.h"
#include "../framework/state_map.h"
#include "thornley_nutrient_dynamics.h"

/**
 * @class thornley_transport
 *
 * @brief Part of Thornley model for plant growth; this module uses the outputs
 * from the transport calculator to compute derivatives. See that module
 * for more information.
 */
class thornley_transport : public differential_module
{
   public:
    thornley_transport(
        std::vector<thornley_nutrient_dynamics::transport_link> const& organ_links,
        state_map const& input_quantities,
        state_map* output_quantities);

    static std::vector<std::string> get_inputs();
    static std::vector<std::string> get_outputs();
    static std::string get_name() { return "thornley_tranport"; }

   private:
    std::vector<thornley_nutrient_dynamics::transport_link> const organ_links;

    // Pointers to input parameters
    std::vector<double const*> const substrate_transport_ips;

    // Pointers to output parameters
    std::vector<std::pair<double*, double*>> const substrate_carbon_ops;

   protected:
    void do_multi_organ_operation() const;
    static std::vector<std::string> get_inputs(std::vector<thornley_nutrient_dynamics::transport_link> const& organ_links);
    static std::vector<std::string> get_outputs(std::vector<thornley_nutrient_dynamics::transport_link> const& organ_links);
};

#endif
