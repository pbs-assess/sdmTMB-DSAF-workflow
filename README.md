# sdmTMB DSAF workflow

An sdmTMB modelling workflow for DFO's Distribution Shifts Assessment Framework
(DSAF). It fits spatiotemporal species distribution models to survey data,
projects them forward under climate scenarios, and calculates a standard set of
distribution-shift indicators.

**Status:** IN ACTIVE DEVELOPMENT - expect frequent, possibly code-breaking
refactors. The pipeline is currently being built and tested on example data:
Pacific spiny dogfish (WCVI synoptic survey, from `surveyjoin`) with GLORYS or a
placeholder downscaled bottom-temperature projection (SSP1-2.6, ensemble
median).

## Setup

1. Install the development version of sdmTMB: `pak::pak("sdmTMB/sdmTMB")`.
2. The `dshift` package (shift-indicator and covariate functions) lives in this
   repo (`dshift/`). You don't need to install it: `R/functions.R` loads it from
   source with `pkgload::load_all()` every time, so the pipeline always uses
   the current version. Just make sure `pkgload` is installed (it comes with
   `devtools`). If you want `library(dshift)` outside the pipeline, run
   `pak::local_install("dshift")`.
3. Install the other packages listed under [Dependencies](#dependencies).

## Workflow

| Step | File | Output |
|------|------|--------|
| 0 | `R/00-settings.R` | Directory paths and plot theme, sourced by every stage |
| 1 | `R/01-data-prep.R` (spin script) | Survey sets with covariates, and raw covariates for each grid cell over the historical + projection years (`grid-covariates.rds`) |
| 2 | `R/02-fit-model.qmd` | Fitted sdmTMB models (`model-outputs/fits/*.rds`); historical prediction maps |
| 3 | `R/03-get-shift-indicators.qmd` | Simulation draws from `project()` (`model-outputs/sims_*.qs2`) and indicator summaries; renders to PDF |

**Step 1: data prep.** Attaches raw covariates (bottom temperature, depth) to
the survey sets and to the grid cells with `attach_covariate()`, and puts both
on a common grid and CRS. Nothing here depends on the model: the prediction grid
itself (year range, covariate scaling, factor reference levels) is built in
steps 2 and 3 by joining `grid-covariates.rds` to the survey grid by `cell`. All inputs are set at the top of the document, so this
is the step to adapt for your own data (see [Using your own data](#using-your-own-data)).
It is a spin script: run it line by line as plain R, or render the report with
`quarto render R/01-data-prep.R`.
With a monthly covariate file, each set gets the layer for its tow date, and the
grid gets the same month in each year (`predict_month`). The dogfish demo
splits WCVI into two example regions (north/south at Y = 5420 km). Cells that
straddle the boundary are split, and each part gets its own area.

**Step 2: model fitting.** Selects the surveys to model (`surveys_use`; step 1
prepares every set in the sets file), centres and scales the covariates, fits
spatiotemporal models with random-walk or AR(1) structure (needed for
forecasting) using `poly(x, 2, raw = TRUE)` on the scaled covariates, and
compares them by AIC. New data (prediction grids) are scaled with the means and
SDs of `fit$data`. Diagnostics are still a placeholder.

**Step 3: shift indicators.** Builds the prediction grid for the forecast
horizon, draws `project()` simulations onto it, and calculates these
indicators, both by region and for the whole domain:

- abundance trend (as mean log density)
- proportional abundance
- centre of gravity (E/W, N/S)
- effective area occupied
- range edges (leading/trailing)
- mean density-weighted depth

Each indicator gets per-year estimates with intervals, plus a linear trend for
the historical and projection periods. The calculations work directly on the
matrix of simulation draws (cells × sims).

Scripts outside the 01 → 03 chain that read the old `prediction-grid.rds` or
`wcvi-cell-polygons.rds` (`prep-slides.R`, `jensens-inequality.R`,
`dogfish-retro.R`, `verify-shift-indicators-*.R`, `presentations/`) need
updating: build the grid as in 03 instead.

The prediction grid can be historical-only (no projection years). A region
shapefile (e.g. NAFO divisions) is optional and only used to outline regions on
maps; it needs a `region` column matching `region` in the grid.

## Using your own data

`01-data-prep.R` needs four inputs. The demo scripts show how each was made
for dogfish:

| Input | Dogfish demo | Notes |
|-------|--------------|-------|
| Survey sets | `R/download-surveyjoin-data.R` | One row per set, with tow date, catch, effort, lon/lat, depth. Map your columns to the names used downstream in the first chunk. |
| Survey grid | `R/demo-00-get-survey-grid.R` | sf polygons, one row per `cell` x `region`, with `X`/`Y` (km, shared by split cells) and `area` (km²). Footprint clipped to the surveyed GEBCO depth range, so run `demo-00-get-GEBCO.R` first. Covariates are attached in step 1. |
| Bottom temperature | `R/demo-00-get-GLORYS.R` | Downloads monthly GLORYS bottom temperature for your bbox. Needs a free [Copernicus Marine](https://marine.copernicus.eu) account; setup steps are at the top of the script. |
| Depth | `R/demo-00-get-GEBCO.R` | Downloads GEBCO bathymetry for your bbox (no login needed) and saves it as a single layer of positive depth in m. Used for both the sets (fit) and the grid (prediction). |

Covariate files must follow the standard in `?dshift::attach_covariate`: a
NetCDF (or anything `terra::rast()` reads) with a defined CRS, one variable, and
time stored in `terra::time()` (yearly or monthly layers; static covariates like
depth have one layer and no time). ~~Aggregate finer cadences (daily, 5-day) to
monthly first.~~

## Functions

- `dshift/` (package): the matrix-based indicator functions (`make_sim_draws()`,
  `abundance_by_group()`, `summarise_estimates()`, `summarise_trend()`,
  `get_area_occupied_sims()`, `get_range_edge_sims()`, ...) and
  `attach_covariate()`
- `R/functions.R`: loads `dshift` and holds workflow helpers:
  `get_region_borders()`, `check_data_sources()`

## Checks

- `R/tests.R`: checks that split-cell area and abundance are counted correctly
- `R/verify-shift-indicators-tinyvast.R`: runs the indicators on an independent
  (tinyVAST) density matrix
- `R/verify-shift-indicators-matrix.R`: checked the matrix pipeline against the
  original long-format pipeline; needs `R/functions-legacy.R`, which is no
  longer in the repo

## Data

`data-raw/` and `data-generated/` are git-ignored, except
`data-generated/dogfish-syn-wcvi-sets.rds` (the demo survey sets). The helper
scripts `R/demo-00-*.R`, `R/download-*.R`, `R/format-*.R` and
`R/make-*-polygons.R` document how the other inputs (GLORYS, the survey grid,
depth, region polygons) were obtained or built.

### Data access

`R/download-hpfx-data.R` reads its (non-public) download URL from the `HPFX_URL`
environment variable. To set it once:

1. Run `file.edit("~/.Renviron")` and add `HPFX_URL=<the URL>` (ask the DSAF
   team for the URL), ending with a new line.
2. Save, then run `readRenviron("~/.Renviron")` (or restart R).

## Dependencies

- **Pipeline:** sdmTMB (development version), dshift (this repo, loaded via
  pkgload), dplyr, tidyr, purrr, tibble, rlang, ggplot2, scales, sf, terra,
  qs2, data.table, matrixStats, cli, here, pkgload
- **GitHub-only:** gfplot (`pbs-assess/gfplot`; plot theme, `coord_sf_auto()`,
  `XY_to_sf()`), pacea (`pbs-assess/pacea`; BC coastline on the dogfish maps)
- **Reports and plots:** knitr, patchwork, broom, kableExtra, viridis, ggthemes,
  DHARMa
- **Data download (optional):** CopernicusMarine, stars, CFtime, keyring
  (GLORYS); ncdf4 (GEBCO); surveyjoin (GitHub: `DFO-NOAA-Pacific/surveyjoin`)
