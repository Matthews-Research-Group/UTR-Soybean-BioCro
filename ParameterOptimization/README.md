# Parameter optimization

To reproduce the UTR parameters, source the `multiyear_utr_params_optimization.R` script. 

This script calls the `multiyear_UTR_optim()` function which calculates the weighted predicted RMSEs described in Piao et al.

Certain parameters, i.e. Stem, Root, Pod respiration rates, and Leaf, Stem, Root senescence reuse factors, are assumed to be the same to reduce the number of parameters to optimize.

The `soybean_parameter_expansion.R` script has the `optim_params_conversion` function, which expands the reduced number of parameters to optimize back to the original size.

Note: This optimization can take several hours to days to run. 

The optimization results in `Optimization_output_2026_03_02.txt` was completed on University of Illinois Urbana-Champaign Campus Cluster.

Service provider: National Center for Supercomputing Applications (NCSA) 
Service: Illinois Computes
Resource: Illinois Campus Cluster Program (ICCP)
Hardware: Dual AMD EPYC 7713 CPU, 512GB RAM, 25G node

Required R packages (the version included in Models folder):
- BioCro v3.3.1 
- BioCroWater
- UTRSoybeanBML
- DEoptim (tested on version 2.2-8)
