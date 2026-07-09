#ifndef VARYING_SLA_H
#define VARYING_SLA_H

#include "../framework/module.h"
#include "../framework/state_map.h"

namespace BioCroWater 
{
/**
 *  @class varying_SLA
 *
 *  @brief Allows SLA to vary based on development index.
 *
 *  respectively.
 */
class varying_SLA: public direct_module
{
   public:
    varying_SLA(
        state_map const& input_quantities,
        state_map* output_quantities)
        : direct_module{},

          // Get references to input quantities
          DVI{get_input(input_quantities, "DVI")},

          // Get pointers to output quantities
          Sp_op{get_op(output_quantities, "Sp")}
    {
    }
    static string_vector get_inputs();
    static string_vector get_outputs();
    static std::string get_name() { return "varying_SLA"; }

   private:
    // Pointers to input quantities
    double const& DVI;

    // Pointers to output quantities
    double* Sp_op;

    // Main operation
    void do_operation() const;
};

string_vector varying_SLA::get_inputs()
{
    return {
        "DVI"               // dimensionless
    };
}

string_vector varying_SLA::get_outputs()
{
    return {
        "Sp"  // micromol / m^2 / s
    };
}

void varying_SLA::do_operation() const
{
    double sla;
    if(DVI<0.55){
      sla = 2.58;
    }else if(DVI>1.55){
      sla = 2;
    }else{
      sla = -2.392922+12.452842*DVI-6.202293*pow(DVI,2.0);
    }

    update(Sp_op, sla);
}

}  // namespace 
#endif
