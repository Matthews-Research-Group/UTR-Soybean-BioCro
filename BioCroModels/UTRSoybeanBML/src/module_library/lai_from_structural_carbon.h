#ifndef UTRSoybeanBML_LAI_CALCULATOR_H
#define UTRSoybeanBML_LAI_CALCULATOR_H

#include "../framework/module.h"
#include "../framework/state_map.h"

namespace UTRSoybeanBML
{
/**
 * @class lai_from_structural_carbon
 *
 * @brief A UTR version of the original BioCro parameter calculator
 *
 * This module calculates SLA and LAI. 
 * LAI = Specific Leaf Area (SLA) *
 * Leaf structural Carbon *
 * and Carbon to mass conversion factor
 *
 */
class lai_from_structural_carbon : public direct_module
{
   public:
    lai_from_structural_carbon(state_map const& input_quantities, state_map* output_quantities)
        : direct_module{},

          // Get pointers to input quantities
          iSp_ip{get_ip(input_quantities, "iSp")},
          TTc_ip{get_ip(input_quantities, "TTc")},
          Sp_thermal_time_decay_ip{get_ip(input_quantities, "Sp_thermal_time_decay")},
          Leaf_structural_carbon_ip{get_ip(input_quantities, "Leaf_structural_carbon")},
          Leaf_carbon_to_mass_factor_ip(get_ip(input_quantities, "Leaf_carbon_to_mass_factor")),
          // Get pointers to output quantities
          Sp_op{get_op(output_quantities, "Sp")},
          lai_op{get_op(output_quantities, "lai")}   
    {
    }
    static string_vector get_inputs();
    static string_vector get_outputs();
    static std::string get_name() { return "lai_from_structural_carbon"; }

   private:
    // Pointers to input quantities
    const double* iSp_ip;
    const double* TTc_ip;
    const double* Sp_thermal_time_decay_ip;
    const double* Leaf_structural_carbon_ip;
    const double* Leaf_carbon_to_mass_factor_ip;

    // Pointers to output quantities
    double* Sp_op;
    double* lai_op;

    // Main operation
    void do_operation() const;
};

string_vector lai_from_structural_carbon::get_inputs()
{
    return {
        "iSp",
        "TTc",
        "Sp_thermal_time_decay",
        "Leaf_structural_carbon",
        "Leaf_carbon_to_mass_factor"
    };
}

string_vector lai_from_structural_carbon::get_outputs()
{
    return {
        "Sp",
        "lai"
    };
}

void lai_from_structural_carbon::do_operation() const
{
    // Collect inputs and make calculations
    double Sp = (*iSp_ip) - (*TTc_ip) * (*Sp_thermal_time_decay_ip);

    // Update the output quantity list
    update(Sp_op, Sp);
    update(lai_op, *Leaf_structural_carbon_ip * *Leaf_carbon_to_mass_factor_ip * Sp);
}

}  // namespace UTRSoybeanBML
#endif
