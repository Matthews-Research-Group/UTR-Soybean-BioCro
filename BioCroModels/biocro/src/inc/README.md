# biocro/inc

## Overview

This repository contains external libraries that are included with BioCro. At
the moment, it consists of a subset of the boost C++ libraries. See `NEWS.md`
for more details.

## Boost

Boost does not assure backward compatibility, so changes to boost could break
BioCro. Thus, we don't want to link our code to a user supplied Boost
installation, and we include a version with BioCro instead. Boost is very large,
so we want to include only the necessary parts.

BioCro developers should see `adding_the_boost_libraries.md` for
instructions on extracting parts of the boost libraries.

### License

The boost library is licensed under version 1.0 of the boost license. A copy
of this license is included here in the file `LICENSE_1_0.txt`.
