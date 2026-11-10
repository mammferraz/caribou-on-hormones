# ============================================================
# 04_climate_road_analyses.R
#
# PART 4 — COMBINED CLIMATE × WINTER ROAD EFFECTS
# Purpose:
#   • Additive and tensor-product interaction models for
#     snow depth + road proximity on cortisol and T3
#   • Year as fixed effect; survey as random intercept
# ============================================================

source("00_data_processing.R")
source("00b_road_layers.R")

outdir <- "caribou_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# I/O — re-load pre-computed lags and road distances
# (run 02 and 03 first to generate these files, OR re-derive
#  them here by sourcing the relevant sections)
# ============================================================

infile_climate <- "climatedaily_all_surveys.csv"
infile_roads   <- "MNR_Road_Segments.geojson"
infile_lookup  <- "winter_roads_all2.csv"
infile_open    <- "Roads_open.csv"

# ============================================================
# NEAREST CLIMATE STATION
# NOTE: climate_group is assigned PER SAMPLE in 00_data_processing.R,
# based on the nearest station to each sample's own GPS coordinates —
# NOT by field campaign (survey_field) and NOT by boundary-assigned
# range (survey). Per-campaign/per-range lookup was checked and found
# too coarse (e.g. ~38% of Missisa-range samples are geographically
# closer to James Bay's station than to Missisa's own nominal
# station). See 00_data_processing.R for the nearest-station logic.
# ============================================================

# ============================================================
# CLIMATE LAGS
# FT3/T3 models below use the 5-day short-lag window (unprefixed
# columns, e.g. snow_on_ground). FCM/cortisol models use a 2-day
# short-lag window (columns suffixed _fcm2, e.g. snow_on_ground_fcm2),
# per the FCM/FT3 lag-window split adopted in 02_climate_analyses.R
# (see that script's header for the rationale: Ashley et al. 2011
# ACTH-challenge data for FCM vs. no caribou-specific FT3 study).
# ============================================================

clim_all <- read_csv(infile_climate, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(date = as.Date(local_date), local_year = as.integer(local_year),
         doy = yday(date), month = month(date))

clim_winter <- clim_all %>% filter(month %in% c(1, 2, 3))
clim_split  <- split(clim_winter, clim_winter$climate_group)

sample_dates <- dat %>%
  filter(!is.na(date), !is.na(climate_group)) %>%
  distinct(pcid, survey, survey_field, year, date, climate_group)

calc_lag_mean <- function(target_date, clim_df, lag_start, lag_end) {
  clim_df %>%
    filter(date <  target_date - days(lag_start - 1),
           date >= target_date - days(lag_end)) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
}

climate_lags <- sample_dates %>%
  rowwise() %>%
  mutate(lag = list(calc_lag_mean(date, clim_split[[climate_group]], 1, 5))) %>%
  ungroup() %>%
  unnest(lag)

# FCM (cortisol) 2-day short-lag window — suffixed _fcm2 to stay
# distinct from the 5-day (unprefixed) columns above, which FT3 uses.
climate_lags_fcm2 <- sample_dates %>%
  rowwise() %>%
  mutate(lag = list(calc_lag_mean(date, clim_split[[climate_group]], 1, 2))) %>%
  ungroup() %>%
  unnest(lag) %>%
  rename_with(~ paste0(.x, "_fcm2"), -c(pcid, survey, survey_field, year, date, climate_group))

dat_clim <- dat %>%
  left_join(
    climate_lags %>% select(-survey, -survey_field, -climate_group),
    by = c("pcid", "year", "date")
  ) %>%
  left_join(
    climate_lags_fcm2 %>% select(-survey, -survey_field, -climate_group),
    by = c("pcid", "year", "date")
  )

# ============================================================
# ROAD DISTANCES
# ============================================================

# Winter roads now come from the provincial layer via 00b_road_layers.R
# (the old MNR_Road_Segments.geojson is a 541-segment spatial clip of the
# full 968; 7 of the missing segments sit inside the sampling window).
roads_sf <- road_layer("winter") %>% clean_names()

road_lookup <- read_csv(infile_lookup, show_col_types = FALSE) %>% clean_names()
roads_open  <- read_csv(infile_open,   show_col_types = FALSE) %>% clean_names()
roads_named <- roads_sf %>% left_join(road_lookup, by = "ogf_id")

dat_sf <- dat_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

dat_sf$dist_road_km <- nearest_road_km(dat_sf, roads_named)

# NOTE: (pcid, date) is NOT a unique row key — 18 individuals were
# sampled twice on the same day (repeat pellets, e.g. one unassayed),
# so joining on (pcid, date) alone with relationship = "many-to-many"
# silently produced a cross-product for those 18 individuals, inflating
# n beyond the true sample size (previously surfaced as models
# reporting n = 1,400/1,393 despite only 1,380 total samples). Adding
# tz_sample_number makes the join key unique per row, since
# (pcid, date, tz_sample_number) has exactly one match in the data.
dat_roads <- dat_sf %>%
  st_drop_geometry() %>%
  select(pcid, date, tz_sample_number, dist_road_km) %>%
  right_join(dat, by = c("pcid", "date", "tz_sample_number"), relationship = "one-to-one")

# ============================================================
# COMBINED DATASET
# ============================================================

dat_combined <- dat_clim %>%
  left_join(dat_roads %>% select(pcid, date, tz_sample_number, dist_road_km),
            by = c("pcid", "date", "tz_sample_number"), relationship = "one-to-one") %>%
  filter(is.finite(dist_road_km)) %>%
  mutate(survey = factor(survey), year_f = factor(year))
# NOTE: climate finiteness (snow_on_ground for FT3, snow_on_ground_fcm2
# for FCM) is now checked per-hormone below, since the two hormones use
# different short-lag windows and therefore different missingness.

cat("dat_combined row count after fixing the join key (should not exceed 1,380):",
    nrow(dat_combined), "\n")

# ============================================================
# 1) CORTISOL — ADDITIVE AND INTERACTION MODELS
# ============================================================

dat_comb_cort <- dat_combined %>%
  filter(is.finite(cortisol), !is.na(julian), is.finite(snow_on_ground_fcm2))

# Additive
# NOTE: road distance is log-transformed here (k = 5), matching the fix
# applied in 03_winter_road_analyses.R — raw dist_road_km with k = 8 was
# overfitting a right-skewed distance measure, the same issue already
# handled with log() for mine distance (07) and wolf suitability (08).
# NOTE: FCM uses the 2-day short-lag snow variable (snow_on_ground_fcm2);
# FT3 below keeps the original 5-day snow_on_ground (see header note).
m_cort_add <- gam(
  cortisol ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(snow_on_ground_fcm2, k = 4, bs = "tp") +
    s(log(dist_road_km),  k = 5, bs = "tp") +
    s(julian,             k = 8, bs = "tp"),
  data   = dat_comb_cort,
  family = Gamma(link = "log"),
  method = "REML"
)

# Tensor-product interaction (snow × log road distance)
m_cort_int <- gam(
  cortisol ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    te(snow_on_ground_fcm2, log(dist_road_km), k = c(6, 6), bs = "tp") +
    s(julian, k = 8, bs = "tp"),
  data   = dat_comb_cort,
  family = Gamma(link = "log"),
  method = "REML"
)

cat("Cortisol combined model (AIC, BIC):\n")
aic_cort_comb <- AIC(m_cort_add, m_cort_int)
bic_cort_comb <- BIC(m_cort_add, m_cort_int)
print(aic_cort_comb); print(bic_cort_comb)
capture.output(print(aic_cort_comb), print(bic_cort_comb),
               file = file.path(outdir, "GAM_AIC_04_cortisol_combined.txt"))

summary(m_cort_add)
summary(m_cort_int)
capture.output(summary(m_cort_add), file = file.path(outdir, "GAM_summary_04_cortisol_combined_add.txt"))
capture.output(summary(m_cort_int), file = file.path(outdir, "GAM_summary_04_cortisol_combined_int.txt"))

png(file.path(outdir, "Fig_GAM_Cortisol_Combined_Climate_Road.png"), 2400, 1600, res = 300)
plot(m_cort_add, pages = 1)
dev.off()

cort_comb_year_contrasts <- emmeans(m_cort_add, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(cort_comb_year_contrasts), file = file.path(outdir, "GAM_contrasts_04_cortisol_combined_year.txt"))

# ============================================================
# 2) T3 — ADDITIVE AND INTERACTION MODELS
# ============================================================

dat_comb_t3 <- dat_combined %>%
  filter(is.finite(t3), t3 > 0, !is.na(julian), is.finite(snow_on_ground))
# T3/FT3 keeps the original 5-day short-lag snow_on_ground (unchanged).

m_t3_add <- gam(
  t3 ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(snow_on_ground,     k = 4, bs = "tp") +
    s(log(dist_road_km),  k = 5, bs = "tp") +
    s(julian,             k = 8, bs = "tp"),
  data   = dat_comb_t3,
  family = Gamma(link = "log"),
  method = "REML"
)

m_t3_int <- gam(
  t3 ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    te(snow_on_ground, log(dist_road_km), k = c(6, 6), bs = "tp") +
    s(julian, k = 8, bs = "tp"),
  data   = dat_comb_t3,
  family = Gamma(link = "log"),
  method = "REML"
)

cat("T3 combined model (AIC, BIC):\n")
aic_t3_comb <- AIC(m_t3_add, m_t3_int)
bic_t3_comb <- BIC(m_t3_add, m_t3_int)
print(aic_t3_comb); print(bic_t3_comb)
capture.output(print(aic_t3_comb), print(bic_t3_comb),
               file = file.path(outdir, "GAM_AIC_04_t3_combined.txt"))

summary(m_t3_add)
summary(m_t3_int)
capture.output(summary(m_t3_add), file = file.path(outdir, "GAM_summary_04_t3_combined_add.txt"))
capture.output(summary(m_t3_int), file = file.path(outdir, "GAM_summary_04_t3_combined_int.txt"))

png(file.path(outdir, "Fig_GAM_T3_Combined_Climate_Road.png"), 2400, 1600, res = 300)
plot(m_t3_add, pages = 1)
dev.off()

t3_comb_year_contrasts <- emmeans(m_t3_add, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(t3_comb_year_contrasts), file = file.path(outdir, "GAM_contrasts_04_t3_combined_year.txt"))

# ============================================================
# 3) SPATIAL AUTOCORRELATION DIAGNOSTICS (variogram on residuals)
# Purpose: same diagnostic applied to the road-only models in
#   03_winter_road_analyses.R, now extended to the combined
#   climate+road models. Residuals are averaged to one value per unique
#   sampling location (rounded projected coordinates, using the same
#   roads_named CRS as above) before the variogram is computed, since
#   several samples share identical GPS coordinates (repeat visits).
# ============================================================

if (!requireNamespace("gstat", quietly = TRUE)) install.packages("gstat", repos = "https://cloud.r-project.org")
library(gstat)
# Spatial-smooth coordinates are held in EPSG:4326 (degrees) — what every
# s(x, y) term in this project has always used. Stated explicitly rather
# than inherited from a road layer's CRS, so swapping road layers can never
# silently re-project the smooths. Road DISTANCES are separate, in metres.
SPATIAL_SMOOTH_CRS <- 4326


sf_cort_comb <- dat_comb_cort %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(SPATIAL_SMOOTH_CRS)
sf_t3_comb <- dat_comb_t3 %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(SPATIAL_SMOOTH_CRS)

check_spatial_autocorr_comb <- function(resid_vec, data_sf, label) {
  stopifnot(nrow(data_sf) == length(resid_vec))
  data_sf$resid <- resid_vec

  coords <- st_coordinates(data_sf)
  data_sf$x_round <- round(coords[, 1])
  data_sf$y_round <- round(coords[, 2])

  site_resid <- data_sf %>%
    st_drop_geometry() %>%
    group_by(x_round, y_round) %>%
    summarise(resid = mean(resid), n_samples = n(), .groups = "drop")

  cat(sprintf("\n%s: %d rows collapsed to %d unique sampling locations for the variogram\n",
              label, nrow(data_sf), nrow(site_resid)))

  site_sf <- st_as_sf(site_resid, coords = c("x_round", "y_round"), crs = st_crs(data_sf))

  vgram <- gstat::variogram(resid ~ 1, data = site_sf)

  png(file.path(outdir, paste0("Variogram_", label, "_CombinedModel.png")), 1600, 1200, res = 200)
  print(plot(vgram, main = paste0(label, " — variogram of combined climate+road model deviance residuals"),
             xlab = "Distance (m)", ylab = "Semivariance"))
  dev.off()

  vgram_fit <- tryCatch(gstat::fit.variogram(vgram, gstat::vgm("Sph")), error = function(e) NULL)
  if (!is.null(vgram_fit)) {
    cat(sprintf("%s variogram model (nugget / partial sill / range):\n", label))
    print(vgram_fit)
    capture.output(print(vgram_fit),
                   file = file.path(outdir, paste0("Variogram_fit_", label, "_CombinedModel.txt")))
  } else {
    message(label, ": variogram model fit did not converge — inspect the empirical variogram plot directly.")
  }
  list(vgram = vgram, vgram_fit = vgram_fit)
}

vres_cort_comb <- check_spatial_autocorr_comb(residuals(m_cort_add, type = "deviance"), sf_cort_comb, "FCM")
vres_t3_comb   <- check_spatial_autocorr_comb(residuals(m_t3_add,   type = "deviance"), sf_t3_comb,   "FT3")

# ---- Refit additive AND interaction models with s(x, y) added ----

dat_comb_cort_xy <- dat_comb_cort %>%
  bind_cols(st_coordinates(sf_cort_comb) %>% as.data.frame() %>% rename(x_proj = X, y_proj = Y))
dat_comb_t3_xy <- dat_comb_t3 %>%
  bind_cols(st_coordinates(sf_t3_comb) %>% as.data.frame() %>% rename(x_proj = X, y_proj = Y))

m_cort_add_xy <- gam(
  cortisol ~
    s(survey, bs = "re") + year_f + preg_group +
    s(snow_on_ground_fcm2, k = 4, bs = "tp") +
    s(log(dist_road_km), k = 5, bs = "tp") +
    s(julian,            k = 8, bs = "tp") +
    s(x_proj, y_proj, bs = "tp", k = 30),
  data = dat_comb_cort_xy, family = Gamma(link = "log"), method = "REML"
)
m_cort_int_xy <- gam(
  cortisol ~
    s(survey, bs = "re") + year_f + preg_group +
    te(snow_on_ground_fcm2, log(dist_road_km), k = c(6, 6), bs = "tp") +
    s(julian, k = 8, bs = "tp") +
    s(x_proj, y_proj, bs = "tp", k = 30),
  data = dat_comb_cort_xy, family = Gamma(link = "log"), method = "REML"
)
m_t3_add_xy <- gam(
  t3 ~
    s(survey, bs = "re") + year_f + preg_group +
    s(snow_on_ground,    k = 4, bs = "tp") +
    s(log(dist_road_km), k = 5, bs = "tp") +
    s(julian,            k = 8, bs = "tp") +
    s(x_proj, y_proj, bs = "tp", k = 30),
  data = dat_comb_t3_xy, family = Gamma(link = "log"), method = "REML"
)
m_t3_int_xy <- gam(
  t3 ~
    s(survey, bs = "re") + year_f + preg_group +
    te(snow_on_ground, log(dist_road_km), k = c(6, 6), bs = "tp") +
    s(julian, k = 8, bs = "tp") +
    s(x_proj, y_proj, bs = "tp", k = 30),
  data = dat_comb_t3_xy, family = Gamma(link = "log"), method = "REML"
)

cat("\nCombined models with s(x, y) added — AIC/BIC comparison against the originally reported (no-spatial-term) models:\n")
aic_cort_comb_xy <- AIC(m_cort_add, m_cort_add_xy, m_cort_int, m_cort_int_xy)
bic_cort_comb_xy <- BIC(m_cort_add, m_cort_add_xy, m_cort_int, m_cort_int_xy)
aic_t3_comb_xy   <- AIC(m_t3_add,   m_t3_add_xy,   m_t3_int,   m_t3_int_xy)
bic_t3_comb_xy   <- BIC(m_t3_add,   m_t3_add_xy,   m_t3_int,   m_t3_int_xy)
print(aic_cort_comb_xy); print(bic_cort_comb_xy)
print(aic_t3_comb_xy); print(bic_t3_comb_xy)
capture.output(print(aic_cort_comb_xy), print(bic_cort_comb_xy), print(aic_t3_comb_xy), print(bic_t3_comb_xy),
               file = file.path(outdir, "GAM_AIC_04_spatial_xy_comparison.txt"))

vres_cort_comb_after <- check_spatial_autocorr_comb(residuals(m_cort_add_xy, type = "deviance"), sf_cort_comb, "FCM_afterXY")
vres_t3_comb_after   <- check_spatial_autocorr_comb(residuals(m_t3_add_xy,   type = "deviance"), sf_t3_comb,   "FT3_afterXY")

message("Part 3 (spatial autocorrelation diagnostics) complete — outputs in: ", outdir)

# ============================================================
# 4) PRIMARY MODEL DECISION — spatially-adjusted model reported
# Project decision: for landscape covariates (road, mine, wolf, land
# cover, population), the s(x, y)-augmented model is reported as the
# PRIMARY model, regardless of BIC (see 03_winter_road_analyses.R
# Part 2E for full rationale). This applies here to the combined
# climate+road models: m_cort_add_xy / m_cort_int_xy / m_t3_add_xy /
# m_t3_int_xy are the primary reported models. The outputs below
# OVERWRITE the original (non-spatial) GAM_summary/Fig/contrasts files
# from Parts 1-2 so every file on disk reflects the model actually
# reported in the manuscript.
# ============================================================

summary(m_cort_add_xy)
summary(m_cort_int_xy)
capture.output(summary(m_cort_add_xy), file = file.path(outdir, "GAM_summary_04_cortisol_combined_add.txt"))
capture.output(summary(m_cort_int_xy), file = file.path(outdir, "GAM_summary_04_cortisol_combined_int.txt"))

png(file.path(outdir, "Fig_GAM_Cortisol_Combined_Climate_Road.png"), 2400, 1600, res = 300)
plot(m_cort_add_xy, pages = 1)
dev.off()

cort_comb_year_contrasts <- emmeans(m_cort_add_xy, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(cort_comb_year_contrasts), file = file.path(outdir, "GAM_contrasts_04_cortisol_combined_year.txt"))

summary(m_t3_add_xy)
summary(m_t3_int_xy)
capture.output(summary(m_t3_add_xy), file = file.path(outdir, "GAM_summary_04_t3_combined_add.txt"))
capture.output(summary(m_t3_int_xy), file = file.path(outdir, "GAM_summary_04_t3_combined_int.txt"))

png(file.path(outdir, "Fig_GAM_T3_Combined_Climate_Road.png"), 2400, 1600, res = 300)
plot(m_t3_add_xy, pages = 1)
dev.off()

t3_comb_year_contrasts <- emmeans(m_t3_add_xy, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(t3_comb_year_contrasts), file = file.path(outdir, "GAM_contrasts_04_t3_combined_year.txt"))

# Re-run the AIC/BIC model-selection table (additive vs. interaction)
# using the spatially-adjusted models, since that comparison is also
# reported in the manuscript.
cat("Cortisol combined model, spatially-adjusted (AIC, BIC):\n")
aic_cort_comb <- AIC(m_cort_add_xy, m_cort_int_xy)
bic_cort_comb <- BIC(m_cort_add_xy, m_cort_int_xy)
print(aic_cort_comb); print(bic_cort_comb)
capture.output(print(aic_cort_comb), print(bic_cort_comb),
               file = file.path(outdir, "GAM_AIC_04_cortisol_combined.txt"))

cat("T3 combined model, spatially-adjusted (AIC, BIC):\n")
aic_t3_comb <- AIC(m_t3_add_xy, m_t3_int_xy)
bic_t3_comb <- BIC(m_t3_add_xy, m_t3_int_xy)
print(aic_t3_comb); print(bic_t3_comb)
capture.output(print(aic_t3_comb), print(bic_t3_comb),
               file = file.path(outdir, "GAM_AIC_04_t3_combined.txt"))

message("Part 4 (primary model = spatially-adjusted) complete — GAM_summary_04_*_combined_add.txt / GAM_summary_04_*_combined_int.txt now reflect the s(x,y)-augmented models.")

message("Part 4 complete — outputs in: ", outdir)

