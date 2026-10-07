# UTR-Soybean-BioCro

Code and data to reproduce the results in:

> Piao, Ximin, Edward B. Lochocki, Justin M. McGrath, and Megan L. Matthews. “Integrating Carbon Utilization and Transport Processes into a Crop Growth Model Enables the Prediction of Emergent Soybean Carbon Allocation Behavior.” bioRxiv, Aug 28, 2026, 2026.08.27.747615. https://doi.org/10.64898/2026.08.27.747615.

This work applies a utilization–transport–resistance (UTR) carbon allocation framework to [Soybean-BioCro](https://doi.org/10.1093/insilicoplants/diab032). 
The UTR processes are included the **UTRSoybeanBML** BioCro Module Library (BML). 
Soil water is simulated with the **BioCroWater** BML. 
The model is calibrated against 2002 and 2005 SoyFACE biomass data, and the scripts in the "ReproduceFigures" folder regenerate every figure, RMSE statistic, and sensitivity analysis in the manuscript.

## Repository structure

```
UTR-Soybean-BioCro/
├── BioCroModels/                   # R packages used by the model
│   ├── biocro/                     #   BioCro v3.3.1 (core framework)
│   └── soil-tipping-bucket/        #   BioCroWater v1.0.0 (soil water BML)
├── Data/
│   ├── SoyFACE_data/               # 2002–2006 Pioneer93B15 biomass at ambient and elevated CO2
│   ├── EnergyFarm_data/            # 2021–2024 LD11 biomass (Energy Farm, Urbana, IL)
│   ├── 2022_Carb_data/             # 2022 LD11 leaf and stem carbohydrate measurements
│   ├── Weather_data/               # Hourly weather for SoyFACE and Energy Farm growing seasons
│   └── Soybean-BioCro_Parameters/  # Model setup scripts and UTR parameters
├── ParameterOptimization/          # DEoptim calibration of the UTR parameters
└── ReproduceFigures/               # Scripts that produce all manuscript figures
```

## Installation

### Requirements

- [R](https://cran.r-project.org/)
  [Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows, Xcode on macOS, or gcc/g++ ≥ 4.9.3 on Linux.
- [RStudio](https://posit.co/download/rstudio-desktop/) is recommended.

### 1. Install the BioCro packages

From the repository root:

```sh
R CMD INSTALL BioCroModels/biocro
R CMD INSTALL BioCroModels/soil-tipping-bucket
```

You also need [UTRSoybeanBML](https://github.com/Matthews-Research-Group/UTRSoybeanBML). Install it the same way, with `R CMD INSTALL path/to/UTRSoybeanBML`.

### 2. Install CRAN dependencies

```r
install.packages(c(
  "DEoptim",                                     # parameter optimization
  "ggplot2", "grid", "gridExtra", "patchwork",   # plotting
  "dplyr", "tidyr", "purrr", "reshape2", "scales",
  "ggtext", "rstudioapi"
))
```

## Usage

### Reproduce the figures

Source `ReproduceFigures/generate_all_figures.R`. Figures are written to `ReproduceFigures/GeneratedFigures/`.

### Re-run the parameter optimization (optional)

Source `ParameterOptimization/multiyear_utr_params_optimization.R`.

> **Note:** The optimization can take several hours to days. The results used in the manuscript are already included in `ParameterOptimization/Optmization_output_2026-08-04.txt`. They were produced on the Illinois Campus Cluster (NCSA), and the figure scripts read them directly, so you do not need to re-run the optimization to reproduce the figures.

## Citation

If you use this code or data, please cite 
> Piao, Ximin, Edward B. Lochocki, Justin M. McGrath, and Megan L. Matthews. “Integrating Carbon Utilization and Transport Processes into a Crop Growth Model Enables the Prediction of Emergent Soybean Carbon Allocation Behavior.” bioRxiv, January 1, 2026, 2026.08.27.747615. https://doi.org/10.64898/2026.08.27.747615.

