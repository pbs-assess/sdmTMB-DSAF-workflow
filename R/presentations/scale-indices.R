#' Scale indices independently within each level of a grouping variable
#'
#' Splits `dat` by `.groups` and scales each split to its own reference
#' (max upr or geomean of est) -- e.g. to remove an intercept offset
#' between model variants before comparing shape.
#'
#' @param min_year,max_year If supplied, the scaling reference is computed
#'   from years in `[min_year, max_year]` only -- e.g. the historical
#'   period, or a single year if `min_year == max_year` -- rather than
#'   including projected years.
scale_indices <- function(dat, .groups,
  scale_by = c("upr", "geomean"),
  space = c("natural", "log"),
  min_year = NULL,
  max_year = NULL
) {

  scale_by <- match.arg(scale_by, c("upr", "geomean"))
  space <- match.arg(space, c("natural", "log"))
  combine <- if (space == "log") `-` else `/`
  avg_fun <- if (space == "log") {
    \(x, ...) mean(x, ...)
  } else {
    \(x, ...) exp(mean(log(x), ...))
  }

  dat_split <- if (is.null(.groups)) {
    list(dat)
  } else {
    dat |> group_by(across(all_of(.groups))) |> group_split()
  }

  dat_split |>
    purrr::map_dfr(\(x) {
      # reference is computed from years in [min_year, max_year] only, so
      # projected years don't bias the scaling baseline; combine() still
      # applies to the full series
      ref <- rep(TRUE, nrow(x))
      if (!is.null(min_year)) ref <- ref & x$year >= min_year
      if (!is.null(max_year)) ref <- ref & x$year <= max_year

      if (sum(!is.na(x$est[ref]))) {
        reference <- switch(scale_by,
          upr = max(x$upr[ref], na.rm = TRUE),
          geomean = avg_fun(x$est[ref], na.rm = TRUE)
        )
        mutate(x,
          est = combine(est, reference),
          lwr = combine(lwr, reference),
          upr = combine(upr, reference)
        )
      } else {
        mutate(x, est = NA_real_, lwr = NA_real_, upr = NA_real_)
      }
    })
}
