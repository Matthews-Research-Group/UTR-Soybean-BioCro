#include "module_library.h"
#include "../framework/module_creator.h"  // for create_mc

// Include all the header files that define the modules.
#include "example_module.h"
#include "soil_evaporation.h"
#include "soil_surface_runoff.h"
#include "soil_type_selector.hpp"
#include "soil_water_downflow.h"
#include "soil_water_tiledrain.h"
#include "soil_water_upflow.h"
#include "multi_layer_soil_profile.hpp"
#include "multilayer_soil_profile_avg.hpp"
#include "soil_water_uptake.h"
#include "varying_SLA.h"
#include "surface_excess_water.hpp"

creator_map BioCroWater::module_library::library_entries =
{
    {"example_module", &create_mc<example_module>},
    {"soil_evaporation", &create_mc<soil_evaporation>},
    {"soil_surface_runoff", &create_mc<soil_surface_runoff>},
    {"soil_type_selector", &create_mc<soil_type_selector>},
    {"soil_water_downflow", &create_mc<soil_water_downflow>},
    {"soil_water_tiledrain", &create_mc<soil_water_tiledrain>},
    {"soil_water_upflow", &create_mc<soil_water_upflow>},
    {"multi_layer_soil_profile", &create_mc<multi_layer_soil_profile>},
    {"surface_excess_water", &create_mc<surface_excess_water>},
    {"multilayer_soil_profile_avg", &create_mc<multilayer_soil_profile_avg>},
    {"soil_water_uptake", &create_mc<soil_water_uptake>},
    {"varying_SLA", &create_mc<varying_SLA>}
};
