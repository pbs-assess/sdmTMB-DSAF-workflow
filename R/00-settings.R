# Define directories in once place to make it easier to change downstream pointers
data_raw_dir <- here::here("data-raw")

data_gen_dir <- here::here("data-generated")
dir.create(data_gen_dir, showWarnings = FALSE, recursive = TRUE)

model_input_dir <- file.path(data_gen_dir, "data-for-sdmTMB")
dir.create(model_input_dir, showWarnings = FALSE, recursive = TRUE)

# Note - in the end this isn't something we will cache for all species; but it
# is handy to save them while we are testing.
model_out_dir <- file.path(data_gen_dir, "model-outputs")
dir.create(model_out_dir, showWarnings = FALSE, recursive = TRUE)

fit_dir <- file.path(model_out_dir, "fits")
dir.create(fit_dir, showWarnings = FALSE)


ggplot2::theme_set(gfplot::theme_pbs())

blues <- bayesplot::color_scheme_get("blue")

# @Check: this is the wording to use?
pvalue_labels <- dplyr::tribble(
  ~p.value, ~p.significant,
  0.0005,   "Very strong",
  0.005,    "Strong",
  0.03,     "Moderate",
  0.5,      "Not significant"
)
