// WARNING: This file was included in this package by the BioCro skeleton module
// library and should not be manually edited.

#ifndef R_MODULE_LIBRARY_H
#define R_MODULE_LIBRARY_H

#include <Rinternals.h>  // for SEXP

extern "C" SEXP R_module_creators(SEXP module_names);
extern "C" SEXP R_get_all_modules();
extern "C" SEXP R_get_all_quantities();

#endif
