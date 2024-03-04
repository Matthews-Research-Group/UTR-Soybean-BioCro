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