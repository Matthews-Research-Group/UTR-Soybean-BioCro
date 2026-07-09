<!--
This file should document all pull requests and all user-visible changes.

When a pull request is completed, changes made should be added to a section at
the top of this file called "# Unreleased". All changes should be categorized
under "## MAJOR CHANGES", "## MINOR CHANGES", or "## BUG FIXES" following the
major.minor.patch structure of semantic versioning. When applicable, entries
should include direct links to the relevant pull requests.

Then, when a new release is made, "# Unreleased" should be replaced by a heading
with the new version number, such as "# CHANGES IN BioCro C++ Framework VERSION
2.0.0." This section will combine all of the release notes from all of the pull
requests merged in since the previous release.

Subsequent commits will then include a new "Unreleased" section in preparation
for the next release.
-->

# BioCro C++ Framework VERSION 1.1.3

- Replaced `fabs` by `std::abs` in `quadratic_root.cpp` and
  `numerical_jacobian.h` since `std::abs` is more correct for C++ code

# BioCro C++ Framework VERSION 1.1.2

- Replaced `std::abs` by `fabs` in `quadratic_root.cpp` since `std::abs` is for
  integer types
- Removed un-implemented overloaded `abs` defined in `state_map.h`

# BioCro C++ Framework VERSION 1.1.1

- Updated README with links to the new stable location for the BioCro R
  package and its documentation

# BioCro C++ Framework VERSION 1.1.0

- Two new values were added to the `physical_constants` namespace (defined in
  `constants.h`: `dr_stomata` and `dr_boundary`.
- Several new functions for finding the roots of quadratic equations were added
  (see `quadratic_root.h` and `quadratic_root.cpp`).
- Both overloaded versions of `dynamical_system::update_drivers` were updated to
  prevent out-of-bounds errors when accessing vector elements, an important bug
  fix.

# BioCro C++ Framework VERSION 1.0.0

- This is the first version of the standalone BioCro C++ framework. Previously,
  this code had been included as part of the BioCro R package source code.
