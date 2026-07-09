#ifndef THORNLEY_UTILIZATION_CALCULATOR_H
#define THORNLEY_UTILIZATION_CALCULATOR_H

#include <vector>
#include <string>
#include "../framework/module.h"
#include "../framework/state_map.h"
#include "thornley_nutrient_dynamics.h"

/**
 * @class thornley_utilization_calculator
 * 
 * @brief Part of the Thornley model for plant growth; although these calculations
 * could be calculated by the utilization_growth derivative module in principle, we
 * perform them here in a steady state model so they are available as outputs of a
 * BioCro simulation.
 * 
 * This model and its associated modules are primarily based on Thornley, J. H. M.
 * "A Model to Describe the Partitioning of Photosynthate during Vegetative Plant
 * Growth." Ann Bot 36, 419–430 (1972).
 * 
 * In this model, two types of carbon are distinguished (`substrate` carbon and
 * `structural` carbon) and each organ has its own separate pools of substrate and
 * structural carbon.
 * 
 * Substrate carbon is ultimately derived from photosynthesis. It is utilized for
 * growth and respiration, and it is transported between organs.
 * 
 * Structural carbon is lost during respiration, with some of the lost mass going to
 * the litter and the remainder being reclaimed as substrate carbon. Moreover, some
 * organs (called `storage organs`) can convert structural carbon to substrate carbon
 * in a process that is distinct from senescence.
 * 
 * A key assumption of this model is that substrate carbon is converted to structural
 * carbon via a catalyzed reaction following Hill kinetics. We refer to this process as
 * "substrate utilization".  In other words, the utilization rate for structural carbon
 * within a particular plant organ is determined by
 * 
 * rate = max_rate / (1 + (km / [substrate_carbon])^n),
 * 
 * where `n` is the Hill coefficient (usually 1 or 2), `[substrate_carbon]` is the
 * concentration of substrate carbon in the organ, and `km` is the value of
 * `[substrate_carbon]` that produces a rate at half of the maximum value.
 * 
 * The substrate carbon concentration is generally unavailable, but this issue can be
 * addressed by noting that 
 * 
 * [substrate_carbon] = (amount of substrate carbon in organ) / (volume of organ)
 *                    ~ (amount of substrate carbon in organ) / (mass of organ)
 *                    ~ (amount of substrate carbon in organ) / (amount of structural carbon in organ),
 * 
 * where we have assumed that the volume of the organ should be proportional to its mass,
 * and its mass should be proprtional to the amount of structural carbon it contains. In
 * his paper, Thornley refers to `(amount of substrate carbon in organ) / (mass of organ)`
 * as the "mass fraction" and uses it in place of the substrate carbon concentration.
 * Following this replacement, the reaction rate becomes
 * 
 * rate = max_rate / (1 + (km / mass_fraction)^n),
 * 
 * where `km` is now understood to be the value of `mass_fraction` that produces a rate at
 * half of the maximum value.
 * 
 * In BioCro, we are generally more interested in the total dry mass of an organ rather than
 * just the amount of carbon it contains, since dry mass is easier to measure. For this reason,
 * a user-supplied conversion factor `carbon_to_mass_factor` is used to convert
 * structural carbon to a real mass value (expressed on a ground area basis).
 */
class thornley_utilization_calculator : public direct_module
{
   public:
    thornley_utilization_calculator(
        std::vector<thornley_nutrient_dynamics::organ> const& organs,
        state_map const& input_quantities,
        state_map* output_quantities);

    static std::vector<std::string> get_inputs();
    static std::vector<std::string> get_outputs();
    static std::string get_name() { return "thornley_utilization_calculator"; }

   private:
    std::vector<thornley_nutrient_dynamics::organ> const organs;

    // Pointers to input parameters
    std::vector<double const*> const structural_carbon_ips;
    std::vector<double const*> const substrate_carbon_ips;
    std::vector<double const*> const utilization_rate_constant_ips;
    std::vector<double const*> const utilization_km_ips;
    std::vector<double const*> const senescence_fraction_max_ips;
    std::vector<double const*> const senescence_alpha_ips;
    std::vector<double const*> const senescence_beta_ips;
    std::vector<double const*> const senescence_reuse_factor_ips;
    
    // References to input parameters
    double const& Pod_start_dvi;
    double const& stop_growth_dvi;
    double const& DVI;

    // Pointers to output parameters
    std::vector<double*> const utilization_rate_ops;
    std::vector<double*> const structural_senescence_rate_ops;
    std::vector<double*> const substrate_senescence_rate_ops;

   protected:
    void do_multi_organ_operation() const;
    static std::vector<std::string> get_inputs(std::vector<thornley_nutrient_dynamics::organ> const& organs);
    static std::vector<std::string> get_outputs(std::vector<thornley_nutrient_dynamics::organ> const& organs);
};

#endif
