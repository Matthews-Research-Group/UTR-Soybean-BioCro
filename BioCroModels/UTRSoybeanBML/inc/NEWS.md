# biocro/boost version 1.0.0

- This version is based on version 1.71.0 of the boost libraries, but changes
  have been made to address several issues identified by `R CMD check`:

  - Some paths have been shortened to comply with CRAN requirements for path
    lengths

  - Replaced instances of `sprintf` in boost files with `snprintf` to avoid
    `R CMD check` warnings

  - Removed usage of deprecated C++ standard library functions
    `std::unary_function` and `std::binary_function` to address `R CMD check`
    warnings; this is a temporary manual fix that should not be required after
    version 1.84 of the boost libraries is released.

- This version is used by version 3.0.2 of the BioCro R package
