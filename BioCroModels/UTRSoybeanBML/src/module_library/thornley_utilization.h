#ifndef THORNLEY_UTILIZATION_H
#define THORNLEY_UTILIZATION_H

#include <vector>
#include <string>
#include "../framework/module.h"
#include "../framework/state_map.h"
#include "thornley_nutrient_dynamics.h"

/**
 * @class thornley_utilization
 * 
 * @brief Part of Thornley model for plant growth; this module uses the outputs
 * from the utilization growth calculator to compute derivatives. See that module
 * for more information.
 */
class thornley_utilization : public differential_module
{
   public:
    thornley_utilization(
        std::vector<thornley_nutrient_dynamics::organ> const& organs,
        state_map const& input_quantities,
        state_map* output_quantities);

    static std::vector<std::string> get_inputs();
    static std::vector<std::string> get_outputs();
    static std::string get_name() { return "thornley_utilization"; }

   private:
    std::vector<thornley_nutrient_dynamics::organ> const organs;

    // Pointers to input parameters
    std::vector<double const*> const carbon_to_mass_factor_ips;
    std::vector<double const*> const substrate_carbon_source_rate_ips;
    std::vector<double const*> const utilization_rate_ips;
    std::vector<double const*> const respiration_factor_ips;
    std::vector<double const*> const structural_senescence_rate_ips;
    std::vector<double const*> const substrate_senescence_rate_ips;

    // References to input parameters
    double const& stop_growth_dvi;
    double const& DVI;

    // Pointers to output parameters
    std::vector<double*> const structural_carbon_ops;
    std::vector<double*> const substrate_carbon_ops;
    std::vector<double*> const respiration_loss_ops;
    std::vector<double*> const senescence_loss_ops;

   protected:
    void do_multi_organ_operation() const;
    static std::vector<std::string> get_inputs(std::vector<thornley_nutrient_dynamics::organ> const& organs);
    static std::vector<std::string> get_outputs(std::vector<thornley_nutrient_dynamics::organ> const& organs);
};

#endif
