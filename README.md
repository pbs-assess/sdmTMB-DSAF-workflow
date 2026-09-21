# sdmTMB DSAF workflow

An sdmTMB modelling workflow for DFO's Distribution Shifts Assessment Framework
(DSAF). It fits spatiotemporal species distribution models to survey data,
projects them forward under climate scenarios, and calculates a standard set of
distribution-shift indicators.

**Status:** IN ACTIVE DEVELOPMENT - expect frequent, possibly code-breaking refactors. The pipeline is currently being built and tested on
placeholder data: Pacific spiny dogfish (WCVI synoptic survey) and a
placeholder downscaled bottom-temperature projection (SSP1-2.6, ensemble
median).

## Workflow

| Step | File | Output (in `data-generated/`) |
|------|------|-------------------------------|
| 0 | `R/00-settings.R` | Directory paths and plot theme, sourced by every stage |
| 1 | `R/01-data-prep.qmd` | Survey data with covariates, a prediction grid covering historical + projection years, and cell-by-region polygons |
| 2 | `R/02-fit-model.qmd` | Fitted sdmTMB models (`model-outputs/fits/*.rds`) |
| 3 | `R/03-get-shift-indicators.qmd` | Simulation draws from `project()` (`model-outputs/sims_*.qs2`) and indicator summaries; renders to PDF |

**Step 1: data prep.** Aligns survey sets and gridded covariates (bottom
temperature, depth) onto a common grid and CRS. Also splits WCVI into two
example regions (north/south at Y = 5420 km). Cells that straddle the boundary
are split, and each part gets its own area.

**Step 2: model fitting.** Fits spatiotemporal models with random-walk or
AR(1) structure (needed for forecasting). Diagnostics are still a placeholder.

**Step 3: shift indicators.** Draws `project()` simulations onto the prediction
grid and calculates these indicators, both by region and for the whole domain:

- abundance trend (as mean log density)
- proportional abundance
- centre of gravity (E/W, N/S)
- effective area occupied
- range edges (leading/trailing)
- mean density-weighted depth

Each indicator gets per-year estimates with intervals, plus a linear trend for
the historical and projection periods.

The calculations work directly on the matrix of simulation draws
(cells × sims).

## Functions

- `R/functions.R`: everything the pipeline sources -- the matrix-based
  indicator functions (`make_sim_draws()`, `abundance_by_group()`,
  `summarise_estimates()`, `summarise_trend()`, `get_area_occupied_sims()`,
  `get_range_edge_sims()`, ...) plus plotting and spatial helpers
  (`plot_field()`, `get_region_borders()`, ...)
- `R/functions-legacy.R`: the original long-format indicator functions, used
  only by `R/verify-shift-indicators-matrix.R`

## Checks

- `R/tests.R`: checks that split-cell area and abundance are counted correctly
- `R/verify-shift-indicators-matrix.R`: checks the matrix pipeline against the
  long-format pipeline
- `R/verify-shift-indicators-tinyvast.R`: runs the indicators on an independent
  (tinyVAST) density matrix

## Data

`data-raw/` and `data-generated/` are git-ignored and not included in the repo.
The helper scripts `R/demo-00-*.R`, `R/demo-01-make-prediction-grid.R` and
`R/make-*-polygons.R` document how some inputs (GLORYS, the survey grid,
region polygons) were obtained or built.

## Dependencies

sdmTMB (development version and project2 branch), dplyr, tidyr, purrr, ggplot2, sf, terra, qs2,
gfplot, data.table, matrixStats.
