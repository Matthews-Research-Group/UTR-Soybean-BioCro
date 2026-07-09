#ifndef THORNLEY_BIOMASS_CALCULATOR_H
#define THORNLEY_BIOMASS_CALCULATOR_H

#include <vector>
#include <string>
#include "../framework/module.h"
#include "../framework/state_map.h"
#include "thornley_nutrient_dynamics.h"

/**
 * @class thornley_biomass
 * 
 * @brief 
 * 
 * This module calculates and updates the organ biomass
 * by converting total C change in mol / m^2 to 
 * organ dry mass change in Mg / ha
 */
class thornley_biomass_calculator : public direct_module
{
   public:
    thornley_biomass_calculator(
        std::vector<thornley_nutrient_dynamics::organ> const& organs,
        state_map const& input_quantities,
        state_map* output_quantities);

    static std::vector<std::string> get_inputs();
    static std::vector<std::string> get_outputs();
    static std::string get_name() {return "thornley_biomass_calculator"; }

   private:
    std::vector<thornley_nutrient_dynamics::organ> const organs;

    // Pointers to input parameters
    std::vector<double const*> const substrate_carbon_ips;
    std::vector<double const*> const structural_carbon_ips;
    
    // References to input parameters
    std::vector<double const*> const carbon_to_mass_factor_ips;

    // Pointers to output parameters
    std::vector<double*> const mass_ops;

   protected:
    void do_multi_organ_operation() const;
    static std::vector<std::string> get_inputs(std::vector<thornley_nutrient_dynamics::organ> const& organs);
    static std::vector<std::string> get_outputs(std::vector<thornley_nutrient_dynamics::organ> const& organs);
};

#endif
