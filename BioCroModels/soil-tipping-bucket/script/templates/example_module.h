#ifndef %1$s_EXAMPLE_MODULE_H
#define %1$s_EXAMPLE_MODULE_H

#include "../framework/module.h"
#include "../framework/state_map.h"

namespace %1$s
{
/**
 * @class example_module
 *
 * @brief An example of a very simple `direct_module`.
 *
 * This module has two inputs (`A` and `B`) and two outputs (`C` and `D`) whose
 * values are calculated according to the following equations:
 *
 * `C = A + B`              [Equation (1)]
 *
 * `D = A - 0.5 * B`        [Equation (2)]
 *
 */
class example_module : public direct_module
{
   public:
    example_module(
        state_map const& input_quantities,
        state_map* output_quantities)
        : direct_module(),

          // Get references to input quantities
          A{get_input(input_quantities, "A")},
          B{get_input(input_quantities, "B")},

          // Get pointers to output quantities
          C_op{get_op(output_quantities, "C")},
          D_op{get_op(output_quantities, "D")}
    {
    }
    static string_vector get_inputs();
    static string_vector get_outputs();
    static std::string get_name() { return "example_module"; }

   private:
    // References to input quantities
    double const& A;
    double const& B;

    // Pointers to output quantities
    double* C_op;
    double* D_op;

    // Main operation
    void do_operation() const;
};

string_vector example_module::get_inputs()
{
    return {
        "A",  // dimensionless
        "B"   // dimensionless
    };
}

string_vector example_module::get_outputs()
{
    return {
        "C",  // dimensionless
        "D"   // dimensionless
    };
}

void example_module::do_operation() const
{
    update(C_op, A + B);
    update(D_op, A - 0.5 * B);
}

}  // namespace %1$s
#endif
