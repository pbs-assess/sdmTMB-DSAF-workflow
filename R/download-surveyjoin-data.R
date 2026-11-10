# Download survey sets (catch + effort, with set dates) from surveyjoin
#
# surveyjoin combines DFO Pacific and NOAA (AFSC, NWFSC) trawl surveys in a
# common format, zeros included. The example uses Pacific Spiny Dogfish; swap in
# your own survey sets here, keeping one row per set with its tow date.
#
# Outputs:
# - data-generated/<species>-survey-sets.rds: all surveys
# - data-generated/<species>-syn-wcvi-sets.rds: SYN WCVI only, committed for the demos

source(here::here("R", "00-settings.R"))

# surveyjoin is GitHub-only:
# pak::pkg_install("DFO-NOAA-Pacific/surveyjoin")

species <- "dogfish"
species_common_name <- "pacific spiny dogfish" # surveyjoin uses lowercase common names

# One-time download to a local cache, then load into the local SQL database
surveyjoin::cache_data()
surveyjoin::load_sql_data()

sets <- surveyjoin::get_data(species_common_name)

dplyr::count(sets, survey_name, effort_units)

saveRDS(sets, file.path(data_gen_dir, paste0(species, "-survey-sets.rds")))

# SYN WCVI subset used by the demos. This one file is committed (git add -f,
# since data-generated/ is ignored) so the demos run without surveyjoin.
sets |>
  dplyr::filter(survey_name == "SYN WCVI") |>
  saveRDS(file.path(data_gen_dir, paste0(species, "-syn-wcvi-sets.rds")))
