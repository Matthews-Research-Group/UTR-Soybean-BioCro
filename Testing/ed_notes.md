## soybean_debug

The `soybean_debug` model definition in `debugging.R` was created from the model
used in `SoyFACE_simulation_testing.R`. To do this, the source code for
`run_biocro` was temporarily modified to include the following lines:

```r
write.csv(as.data.frame(initial_values),            file = 'initial_values.csv',            row.names = FALSE)
write.csv(as.data.frame(parameters),                file = 'parameters.csv',                row.names = FALSE)
write.csv(as.data.frame(direct_module_names),       file = 'direct_module_names.csv',       row.names = FALSE)
write.csv(as.data.frame(differential_module_names), file = 'differential_module_names.csv', row.names = FALSE)
write.csv(as.data.frame(ode_solver),                file = 'ode_solver.csv',                row.names = FALSE)
write.csv(drivers,                                  file = 'drivers.csv',                   row.names = FALSE)
```

Then, the testing script was run. All the inputs to `run_biocro` were saved to
text files so they could be accessed and included in `debugging.R`.

## verbose output

The simulations for 2002, 2004, and 2006 do not report any error messages when
`verbose` is set to `TRUE` in `run_biocro`. The simulation for 2005 results in
the following message:

```
The ODE solver reports the following:
boost::numeric::odeint::integrate_const encountered an error and has returned a partial result:
Max number of iterations exceeded (200).
The dynamical system reports the following:
93300 derivatives were calculated
```

So it looks like the solver gets stuck for this year. This happens before any
NA values are produced.

## checking for NaN

We can figure out whether NaNs first appear as the result of a module
calculation or whether they are passed from the ODE solver. This can be done by
adding a few lines to `dynamical_system::calculate_derivative`, which is found
in `src/framework/dynamical_system.h`. This function should be replaced by the
following, where the lines checking for NaN are new:

```c++
template <typename vector_type, typename time_type>
void dynamical_system::calculate_derivative(const vector_type& x, vector_type& dxdt, const time_type& t)
{
    // Check for NaN in x
    if (!std::all_of(x.begin(), x.end(), [](double v) { return !std::isnan(v); })) {
        throw std::range_error("Thrown in calculate_derivative: NaN found in x (input).");
    }

    ++ncalls;
    update_all_quantities(x, t);
    run_differential_modules(dxdt);

    // Check for NaN in dxdt
    if (!std::all_of(dxdt.begin(), dxdt.end(), [](double v) { return !std::isnan(v); })) {
        throw std::range_error("Thrown in calculate_derivative: NaN found in dxdt (output).");
    }
}
```

Note: it is also necessary to include two libraries:

```c++
#include <algorithm>         // for std::all_of
#include <cmath>             // for std::isnan
```

Errors thrown by this function will be printed to R if `verbose` is `TRUE` when
calling `run_biocro`. When running 2002, 2004, 2005, and 2006 (in that order),
the results are:

```
System startup information:

[pass] No quantities were defined multiple times in the inputs

[pass] All module inputs were properly defined

[pass] All differential module outputs were included in the initial values

[pass] There are no cyclic dependencies among the direct modules.

ODE solver description:
Name: rkck54
Output step size: 1.000000
Relative error tolerance: 0.000100
Absolute error tolerance: 0.000100
Maximum attempts to find a new step size: 200

The ODE solver reports the following:
boost::numeric::odeint::integrate_const encountered an error and has returned a partial result:
Thrown in calculate_derivative: NaN found in dxdt (output).
The dynamical system reports the following:
84407 derivatives were calculated


System startup information:

[pass] No quantities were defined multiple times in the inputs

[pass] All module inputs were properly defined

[pass] All differential module outputs were included in the initial values

[pass] There are no cyclic dependencies among the direct modules.

ODE solver description:
Name: rkck54
Output step size: 1.000000
Relative error tolerance: 0.000100
Absolute error tolerance: 0.000100
Maximum attempts to find a new step size: 200

The ODE solver reports the following:
boost::numeric::odeint::integrate_const encountered an error and has returned a partial result:
Thrown in calculate_derivative: NaN found in dxdt (output).
The dynamical system reports the following:
86549 derivatives were calculated


System startup information:

[pass] No quantities were defined multiple times in the inputs

[pass] All module inputs were properly defined

[pass] All differential module outputs were included in the initial values

[pass] There are no cyclic dependencies among the direct modules.

ODE solver description:
Name: rkck54
Output step size: 1.000000
Relative error tolerance: 0.000100
Absolute error tolerance: 0.000100
Maximum attempts to find a new step size: 200

The ODE solver reports the following:
boost::numeric::odeint::integrate_const encountered an error and has returned a partial result:
Max number of iterations exceeded (200).
The dynamical system reports the following:
93300 derivatives were calculated


System startup information:

[pass] No quantities were defined multiple times in the inputs

[pass] All module inputs were properly defined

[pass] All differential module outputs were included in the initial values

[pass] There are no cyclic dependencies among the direct modules.

ODE solver description:
Name: rkck54
Output step size: 1.000000
Relative error tolerance: 0.000100
Absolute error tolerance: 0.000100
Maximum attempts to find a new step size: 200

The ODE solver reports the following:
boost::numeric::odeint::integrate_const encountered an error and has returned a partial result:
Thrown in calculate_derivative: NaN found in dxdt (output).
The dynamical system reports the following:
84407 derivatives were calculated
```

So it looks like at least one module is returning NaN for a derivative before
the ODE solver starts setting differential quantities to NaN.

## checking for NaN by module outputs

We can figure out which quantity first has NaN for its derivative. This can be
done by adding a few lines to `dynamical_system::run_differential_modules`,
which is found in `src/framework/dynamical_system.h`. This function should be
replaced by the following, where the for loop for running each module and
checking for NaN is new:

```c++
template <class vector_type>
void dynamical_system::run_differential_modules(vector_type& dxdt)
{
    // Reset the derivatives of the differential quantities
    for (auto& x : differential_quantity_derivatives) {
        x.second = 0.0;
    }

    // Run the modules
    for (auto const& m : differential_modules) {
        m->run();

        // Check for NaN
        for (auto& x : differential_quantity_derivatives) {
            if (std::isnan(x.second)) {
                throw std::range_error(x.first);
            }
        }
    }

    // Store the output in the derivative vector
    for (size_t i = 0; i < dxdt.size(); i++) {
        dxdt[i] = *(differential_quantity_ptr_pairs[i].second) * (*timestep_ptr);
    }
}
```

Note: it is also necessary to include one library:

```c++
#include <cmath>             // for std::isnan
```

Errors thrown by this function will be printed to R if `verbose` is `TRUE` when
calling `run_biocro`. When running 2002, 2004, 2005, and 2006 (in that order),
the results are:

```
System startup information:

[pass] No quantities were defined multiple times in the inputs

[pass] All module inputs were properly defined

[pass] All differential module outputs were included in the initial values

[pass] There are no cyclic dependencies among the direct modules.

ODE solver description:
Name: rkck54
Output step size: 1.000000
Relative error tolerance: 0.000100
Absolute error tolerance: 0.000100
Maximum attempts to find a new step size: 200

The ODE solver reports the following:
boost::numeric::odeint::integrate_const encountered an error and has returned a partial result:
Leaf_substrate_carbon
The dynamical system reports the following:
84407 derivatives were calculated


System startup information:

[pass] No quantities were defined multiple times in the inputs

[pass] All module inputs were properly defined

[pass] All differential module outputs were included in the initial values

[pass] There are no cyclic dependencies among the direct modules.

ODE solver description:
Name: rkck54
Output step size: 1.000000
Relative error tolerance: 0.000100
Absolute error tolerance: 0.000100
Maximum attempts to find a new step size: 200

The ODE solver reports the following:
boost::numeric::odeint::integrate_const encountered an error and has returned a partial result:
Leaf_substrate_carbon
The dynamical system reports the following:
86549 derivatives were calculated


System startup information:

[pass] No quantities were defined multiple times in the inputs

[pass] All module inputs were properly defined

[pass] All differential module outputs were included in the initial values

[pass] There are no cyclic dependencies among the direct modules.

ODE solver description:
Name: rkck54
Output step size: 1.000000
Relative error tolerance: 0.000100
Absolute error tolerance: 0.000100
Maximum attempts to find a new step size: 200

The ODE solver reports the following:
boost::numeric::odeint::integrate_const encountered an error and has returned a partial result:
Max number of iterations exceeded (200).
The dynamical system reports the following:
93300 derivatives were calculated


System startup information:

[pass] No quantities were defined multiple times in the inputs

[pass] All module inputs were properly defined

[pass] All differential module outputs were included in the initial values

[pass] There are no cyclic dependencies among the direct modules.

ODE solver description:
Name: rkck54
Output step size: 1.000000
Relative error tolerance: 0.000100
Absolute error tolerance: 0.000100
Maximum attempts to find a new step size: 200

The ODE solver reports the following:
boost::numeric::odeint::integrate_const encountered an error and has returned a partial result:
Leaf_substrate_carbon
The dynamical system reports the following:
84407 derivatives were calculated
```

So it looks like NaN is returned for `Leaf_substrate_carbon` (possibly among
others; this is just the first one encountered). We can get a list of modules
that have this quantity as an output:

```r
all_quantities <- get_all_quantities('UTRSoybeanBML')
all_quantities[all_quantities$quantity_name == 'Leaf_substrate_carbon' & all_quantities$quantity_type == 'output', ]
```

Looks like there are two: `UTRSoybeanBML:thornley_transport_lsrp` and
`UTRSoybeanBML:thornley_utilization_lsrp`. Based on the order of the modules in
the model definition, I think the utilization module gets run first.

## checking UTRSoybeanBML:thornley_utilization_lsrp

Looks like the derivative for `Leaf_substrate_carbon` becomes NaN because
`substrate_carbon_source_rate` is NaN. The value of this quantity is taken from
`canopy_assimilation_rate`.

This was determined by using a modified version of
`src/module_library/thornley_utilization.cpp` in the `UTRSoybeanBML` repository.
A copy of the modified file is included in this directory, which has been
renamed to `DEBUG_thornley_utilization.cpp`.

## checking BioCro:c3_leaf_photosynthesis

Looks like there are several sets of inputs that cause an NaN assimilation rate.
The following are three examples, which are the first times there is an NaN
assimilation rate during 2002, 2004, and 2006, respectively:

```
c3_leaf_photosynthesis:
  absorbed_ppfd = 1.786528e+02
  temp = 2.135132e+01
  rh = 3.850781e-01
  vmax1 = 1.112000e+02
  jmax = 1.950000e+02
  tpu_rate_max = 1.300000e+01
  Rd = 1.280000e+00
  b0 = 8.000000e-03
  b1 = 1.060000e+01
  Gs_min = 1.000000e-03
  Catm = 3.720000e+02
  atmospheric_pressure = 1.013250e+05
  O2 = 2.100000e+02
  theta = 7.600000e-01
  StomataWS = 7.553973e-01
  electrons_per_carboxylation = 4.500000e+00
  electrons_per_oxygenation = 5.250000e+00
  average_absorbed_shortwave = -3.114077e-03
  windspeed = 5.457125e+00
  height = -9.394983e-05
  specific_heat_of_air = 1.010000e+03
  minimum_gbw = 8.000000e-02
  windspeed_height = 5.000000e+00
  beta_PSII = 5.000000e-01

c3_leaf_photosynthesis:
  absorbed_ppfd = 4.484838e+02
  temp = 2.820970e+01
  rh = 2.914652e-01
  vmax1 = 1.112000e+02
  jmax = 1.950000e+02
  tpu_rate_max = 1.300000e+01
  Rd = 1.280000e+00
  b0 = 8.000000e-03
  b1 = 1.060000e+01
  Gs_min = 1.000000e-03
  Catm = 3.720000e+02
  atmospheric_pressure = 1.013250e+05
  O2 = 2.100000e+02
  theta = 7.600000e-01
  StomataWS = 7.500526e-01
  electrons_per_carboxylation = 4.500000e+00
  electrons_per_oxygenation = 5.250000e+00
  average_absorbed_shortwave = -8.330642e-02
  windspeed = 6.027481e+00
  height = -1.000932e-03
  specific_heat_of_air = 1.010000e+03
  minimum_gbw = 8.000000e-02
  windspeed_height = 5.000000e+00
  beta_PSII = 5.000000e-01

c3_leaf_photosynthesis:
  absorbed_ppfd = 1.786528e+02
  temp = 2.135132e+01
  rh = 3.850781e-01
  vmax1 = 1.112000e+02
  jmax = 1.950000e+02
  tpu_rate_max = 1.300000e+01
  Rd = 1.280000e+00
  b0 = 8.000000e-03
  b1 = 1.060000e+01
  Gs_min = 1.000000e-03
  Catm = 3.720000e+02
  atmospheric_pressure = 1.013250e+05
  O2 = 2.100000e+02
  theta = 7.600000e-01
  StomataWS = 7.553973e-01
  electrons_per_carboxylation = 4.500000e+00
  electrons_per_oxygenation = 5.250000e+00
  average_absorbed_shortwave = -3.114077e-03
  windspeed = 5.457125e+00
  height = -9.394983e-05
  specific_heat_of_air = 1.010000e+03
  minimum_gbw = 8.000000e-02
  windspeed_height = 5.000000e+00
  beta_PSII = 5.000000e-01
```

These sets of inputs all have a negative value for `average_absorbed_shortwave`.
This value comes from `sunML`.

This was determined by using a modified version of
`src/module_library/c3_leaf_photosynthesis.cpp` in the `BioCro` repository.
A copy of the modified file is included in this directory, which has been
renamed to `DEBUG_c3_leaf_photosynthesis.cpp`.
