// WARNING: This file was included in this package by the BioCro skeleton module
// library and should not be manually edited.

#include <string>
#include <exception>                       // for std::exception
#include <Rinternals.h>                    // for Rf_error
#include "framework/R_helper_functions.h"  // for r_string_vector_from_vector
#include "framework/state_map.h"           // for string_vector
#include "R_skeleton_version.h"
#include "skeleton_version.h"

using std::string;
using %1$s::skeleton_version;

extern "C" {
/**
 *  @brief Returns the skeleton module library version number as a string.
 */
SEXP R_skeleton_version()
{
    try {
        string_vector result = {skeleton_version};
        return r_string_vector_from_vector(result);
    } catch (std::exception const& e) {
        Rf_error((string("Caught exception in R_skeleton_version: ") + e.what()).c_str());
    } catch (...) {
        Rf_error("Caught unhandled exception in R_skeleton_version.");
    }
}
}
