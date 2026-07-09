# Reproduce figures

Scripts to reproduce all figures and model-fit statistics (MSE) from
Piao et al. 2026, *"Integrating carbon utilization and transport processes
into a crop growth model enables the prediction of emergent soybean carbon
allocation behavior."*

## Usage

Open `generate_all_figures.R` in RStudio (it uses `rstudioapi` to set the
working directory to its own location) and source the whole file. All
figures are written to `GeneratedFigures/`, which is created automatically
if it doesn't already exist. MSE and sensitivity results are printed to the
console as the script runs.

Near the top of the script, three flags control which optional (slower)
figures get produced:

- `RUN_SENSITIVITY` — pod-mass sensitivity to +/-10% input perturbations
- `PLOT_SOURCE_SINK_MANIPULATION` — shade / pod-removal / defoliation / hail treatment figures
- `PLOT_LAI` — simulated-vs-observed LAI figures (SoyFACE and LD11)

## Layout

- `generate_all_figures.R` — main driver script; run this to reproduce everything.
- `PlotScripts/` — plotting scripts/functions sourced by `generate_all_figures.R`
  (one file per figure group: biomass, partitioning, allocation, LAI, and the
  shade/pod-removal/defoliation/hail source-sink manipulation treatments).
- `CalculateScripts/` — `calculate_mse.R` and `calculate_sensitivity.R`, helper
  functions for model-fit (MSE) and sensitivity calculations, sourced directly
  from `generate_all_figures.R`.
- `GeneratedFigures/` — output folder; all `.png` figures land here. Not
  tracked in git — regenerated on every run.

Scripts under `PlotScripts/` and `CalculateScripts/` read data using paths
relative to `ReproduceFigures/` (e.g. `../Data/...`), not their own location
— this works because `source()` doesn't change R's working directory, so it
stays at `ReproduceFigures/` for the whole run.

## Required packages

- BioCro, UTRSoybeanBML, BioCroWater (soybean model and modules)
- ggplot2, grid, gridExtra, patchwork (plotting and layout)
- dplyr, tidyr, purrr, reshape2, scales (data wrangling and formatting)
- ggtext (rich text labels in the defoliation figure)
- rstudioapi (used once, to set the working directory)

## Data dependencies

Expects the sibling `../Data/` (weather, SoyFACE, LD11/EnergyFarm, and 2022
carbohydrate data) and `../ParameterOptimization/` directories from this
repository to be present.
