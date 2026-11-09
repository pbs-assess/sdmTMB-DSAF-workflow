# dshift

Distribution-shift indicators (abundance, log density, range edges, area
occupied) computed from posterior simulation draws of sdmTMB projections.

Draft package under active development; the API may change.

## Installation

dshift lives in a subfolder of the `sdmTMB-DSAF-workflow` repo, so pak needs
the package name given explicitly (otherwise it assumes the package is named
after the repo):

```r
pak::pkg_install("pbs-assess/sdmTMB-DSAF-workflow/dshift")
```

Or with remotes:

```r
remotes::install_github("pbs-assess/sdmTMB-DSAF-workflow", subdir = "dshift")
```

If another package depends on dshift, use the same explicit form in its
`DESCRIPTION`:

```
Imports:
    dshift
Remotes:
    dshift=github::pbs-assess/sdmTMB-DSAF-workflow/dshift
```
