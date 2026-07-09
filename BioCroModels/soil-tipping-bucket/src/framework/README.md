## BioCro C++ framework

The **BioCro C++ framework** defines a set of C++ classes that can be used to
create sets of reusable equations, and to define and solve dynamical systems
formed from such equations. For more information, please see our publication in
_in silico Plants_
([Lochocki _et al._, 2022](https://doi.org/10.1093/insilicoplants/diac003)).

### Using the BioCro framework

The BioCro framework is typically accessed from R using the BioCro R package,
rather than directly from its C++ code. Detailed information about installing
and using the R package is available at the
[biocro/biocro](https://github.com/biocro/biocro) GitHub repository.
It is also possible to use BioCro without the R package; however, instructions
for using BioCro without the R package are not available at this time.

### Requirements

The BioCro C++ framework requires portions of version 1.71.0 of the Boost
software library. Other versions may also work, although the code has only been
tested with that version. The relevant portions of the Boost library are
included with the BioCro R package discussed
[above](#using-the-biocro-framework).

### Making contributions

Please see the contribution guidelines before submitting changes.
These may be found in Chapter One of the _Developer's Manual_ on the
[public BioCro Documentation web site](https://biocro.github.io).

### Software Documentation

See the [public BioCro Documentation web
site](https://biocro.github.io), which includes documentation for the C++
framework.

### License

The BioCro C++ framework is licensed under version 3 or greater of the GNU
Lesser General Public License (LGPL), while the BioCro R package is licensed
under the MIT license. This scheme allows people to freely develop models for
any use (public or private) under the MIT license, but any changes to the
framework that assembles and solves models must make source code changes
available to all users under the LGPL.
