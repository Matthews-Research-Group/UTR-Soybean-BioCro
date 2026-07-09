## %1$s
PUT AN OVERVIEW OF THE MODULE LIBRARY HERE

### Installation

#### Requirements

- The [R environment](https://cran.r-project.org/).
- On Windows, a version of [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
  appropriate for your version of R.
- On Linux, gcc and g++ version 4.9.3 or greater (consult documentation for your
  distribution for installation instructions).
- On MacOS, Xcode.

#### Installation steps

First, obtain a local copy of this repository, which can be accomplished using
either of two methods:
1. If you are new to Git, the easiest way to get a local copy is to install
   GitHub Desktop and use the "Open with GitHub Desktop option in the "Code"
   dropdown on the GitHub repository page.
2. Alternatively, clone this repository using Git on the command line in the
   usual fashion by running `git clone <repository URL>` where
   `<repository URL>` is the URL for this repository. This repository contains a
   Git submodule, so you will need to take the additional step of running
   `git submodule update --init` to obtain it.

After obtaining a local copy of the source code, install the package from the
command line or from within R using one of the following sets of commands. These
assume that the source files are in a directory named  '%1$s'.

- From the command line
```
cd path_to_unzipped_directory
R CMD INSTALL %1$s
```

- Or from within R
```
setwd('path_to_unzipped_directory')
install.packages('%1$s', repos=NULL, type='SOURCE')
```

### An example

The following code will print a list of all modules available in this library,
return information about one of them (`example_module`), and then run that
module:
```
library(BioCro)
library(%1$s)
get_all_modules('%1$s')
module_info('%1$s:example_module')
evaluate_module('%1$s:example_module', list(A = 1, B = 2))
```
For more information about using BioCro modules in R, please see the
[BioCro framework R package](https://github.com/biocro/biocro).

### Source

This package was derived from the
[BioCro skeleton module library](https://github.com/biocro/skelBML). See
`skelBML_description` for more information about the version of the skeleton
library that was used.

### License

The `%1$s` R package is licensed under the MIT license, while the BioCro C++
framework is licensed under version 3 or greater of the GNU Lesser General
Public License (LGPL). This scheme allows people to freely develop models for
any use (public or private) under the MIT license, but any changes to the
framework that assembles and solves models must make source code changes
available to all users under the LGPL. See `LICENSE.note` for more details.
