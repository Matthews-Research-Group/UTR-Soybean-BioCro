#include "module_library.h"
#include "../framework/module_creator.h"  // for create_mc

// Include all the header files that define the modules.
// #include "example_module.h"
#include "lai_from_structural_carbon.h"
#include "thornley_nutrient_dynamics_lsrp.h"

creator_map UTRSoybeanBML::module_library::library_entries =
{
//    {"example_module",                                        &create_mc<example_module>},
    {"lai_from_structural_carbon",                            &create_mc<lai_from_structural_carbon>},
    {"thornley_utilization_calculator_lsrp",                  &create_mc<thornley_utilization_calculator_lsrp>},
    {"thornley_utilization_lsrp",                             &create_mc<thornley_utilization_lsrp>},
    {"thornley_transport_calculator_lsrp",                    &create_mc<thornley_transport_calculator_lsrp>},
    {"thornley_transport_lsrp",                               &create_mc<thornley_transport_lsrp>},
    {"thornley_biomass_calculator_lsrp",                      &create_mc<thornley_biomass_calculator_lsrp>}
};
