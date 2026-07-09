# Reproduce figures
Scripts to reproduce all figures, statistics (MSE), and sensitivity analysis from Piao et al. 2026.

## Usage
Run `generate_all_figures.R`. 
All figures are written to `GeneratedFigures/`, which is created automatically if it doesn't already exist. 
MSE and sensitivity results are printed to the console as the script runs.

Three flags control which optional (slower) figures and analysis get produced:
- `RUN_SENSITIVITY` — final pod mass sensitivity to +/-10% input perturbations
- `PLOT_SOURCE_SINK_MANIPULATION` — shade / pod-removal / defoliation / hail figures
- `PLOT_LAI` — simulated-vs-observed LAI figures (SoyFACE and LD11)

## File Structure
- `generate_all_figures.R` — main script; run this to reproduce everything.
- `PlotScripts/` — plotting scripts/functions sourced by `generate_all_figures.R`
- `CalculateScripts/` — helper functions for mean squared errors (MSEs) and sensitivity calculations.
- `GeneratedFigures/` — figures output folder; all `.png` figures saved here. Regenerated on every run.

## Required packages
- BioCro, UTRSoybeanBML, BioCroWater (BioCro module libraries)
- ggplot2, grid, gridExtra, patchwork (plotting and layout)
- dplyr, tidyr, purrr, reshape2, scales (data wrangling and formatting)
- ggtext (rich text labels)
- rstudioapi (used to set the working directory)