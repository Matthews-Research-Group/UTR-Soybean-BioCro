#include "module_library.h"
#include "../framework/module_creator.h"  // for create_mc

// Include all the header files that define the modules.
#include "example_module.h"

creator_map %1$s::module_library::library_entries =
{
    {"example_module", &create_mc<example_module>}
};
