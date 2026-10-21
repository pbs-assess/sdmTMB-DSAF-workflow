library(dplyr)
library(purrr)
library(stringr)

# To keep the HPFX link private, save it in your .Renviron:
# 1. Open ~/.Renviron with file.edit('~/.Renviron').
# 2. Add the URL `HPFX_URL=https://...`.
# 3. Reload with readRenviron('~/.Renviron') (or restart R).
base_url <- Sys.getenv("HPFX_URL")
if (base_url == "") {
  cli::cli_abort(c(
    "{.envvar HPFX_URL} isn't set.",
    "i" = "Open {.file ~/.Renviron} with {.run file.edit('~/.Renviron')}.",
    " " = "Add the line {.code HPFX_URL=https://...}.",
    "i" = "Reload with {.run readRenviron('~/.Renviron')} (or restart R)."
  ))
}
dest_dir <- here::here("data-raw", "hpfx")

# Recursively list every file under an Apache directory index
list_hpfx_files <- function(url, path = "") {
  hrefs <- readLines(paste0(url, path), warn = FALSE) |>
    str_extract_all('(?<=href=")[^"]+(?=")') |>
    unlist()
  # Drop column-sort links ("?C=N;O=D") and the parent directory ("/dfo/...")
  hrefs <- hrefs[!str_detect(hrefs, "^[?/]")]

  is_dir <- str_ends(hrefs, "/")
  # str_c(), not paste0(): paste0("a/", character(0)) returns "a/", which recurses forever
  files <- tibble(path = str_c(path, hrefs[!is_dir]))
  subdirs <- map(str_c(path, hrefs[is_dir]), \(p) list_hpfx_files(url, p))
  bind_rows(files, subdirs)
}

download_hpfx_files <- function(files, url = base_url, dest = dest_dir) {
  dest_files <- file.path(dest, files$path)
  fs::dir_create(unique(dirname(dest_files)))
  cli::cli_inform("Downloading {nrow(files)} file{?s} to {.file {dest}}")
  # Queued requests count against curl's 10 s connect timeout, so only queue as
  # many as curl will run at once (6 per host); multiplex = FALSE lets all 6 run
  # in parallel instead of waiting on one connection
  batch <- ceiling(seq_along(dest_files) / 6)
  map2(
    split(paste0(url, files$path), batch), split(dest_files, batch),
    \(u, d) curl::multi_download(u, d, resume = TRUE, progress = TRUE, multiplex = FALSE)
  ) |>
    list_rbind()
}

# --- List files ---------------------------------------------------------------
hpfx_files <- list_hpfx_files(base_url) |>
  mutate(folder = dirname(path), file = basename(path))

count(hpfx_files, folder)

# --- Choose what to download -------------------------------------------------
# Use this to download all the data:
# folders <- c(
#   "NAODS-CanOE_Hindcast_bottom_T_O2",
#   "data_anomalies_from_hist-clim/HadGEM3-ES_SSP585"
# )
# file_pattern <- "\\.nc$"

# Use this to pick and choose based on folder and filename
folders <- "NAODS-CanOE_Hindcast_bottom_T_O2"
file_pattern <- "\\.nc$"
# file_pattern <- "_bottomT\\.nc$"
# file_pattern <- "_bottomO2\\.nc$"

to_download <- hpfx_files |>
  filter(
    map(folders, \(f) str_starts(paste0(folder, "/"), fixed(paste0(f, "/")))) |>
      reduce(`|`),
    str_detect(file, file_pattern)
  )

cli::cli_inform("{nrow(to_download)} file{?s} selected from {n_distinct(to_download$folder)} folder{?s}")

results <- download_hpfx_files(to_download)
filter(results, !success | status_code >= 400)
