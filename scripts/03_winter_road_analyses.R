# ============================================================
# 03_winter_road_analyses.R
#
# PART 3 — SPATIOTEMPORAL WINTER ROAD EXPOSURE
# Purpose:
#   • Compute distance from each sample to nearest winter road
#   • Characterise road openness by year (days_open)
#   • Test effects on cortisol and T3
#   • Primary axis: year; survey as random effect and colour
# ============================================================

source("00_data_processing.R")
source("00b_road_layers.R")

outdir <- "caribou_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# I/O
# ============================================================

# WINTER ROAD SOURCE — CORRECTED 2026-07-26
#   Winter roads were previously read from MNR_Road_Segments.geojson,
#   which contains 541 Winter segments. The provincial file
#   MNR_Road_Segments2.geojson contains 968, and the 541 are a strict
#   SUBSET of them — the smaller file is a spatial clip, not a full
#   download. Seven of the 427 absent segments fall inside the sampling
#   window, so distance-to-nearest-winter-road was biased upward for an
#   unknown number of samples. Winter roads now come from the provincial
#   file via road_layer("winter"), matching the all-season layer.
infile_lookup <- "winter_roads_all2.csv"
infile_open   <- "Roads_open.csv"

# ============================================================
# LOAD SPATIAL ROADS
# ============================================================

roads_sf <- road_layer("winter") %>% clean_names()

road_lookup <- read_csv(infile_lookup, show_col_types = FALSE) %>% clean_names()
roads_open  <- read_csv(infile_open,   show_col_types = FALSE) %>% clean_names()
# roads_open: columns road, year, days_open (pre-computed)

# The name lookup only covers the originally-downloaded segments, so the
# newly included ones join to NA. That is fine: road names are used for
# reporting only, never for the distance calculation.
roads_named <- roads_sf %>% left_join(road_lookup, by = "ogf_id")
cat("Winter road segments:", nrow(roads_named),
    "| with a name in the lookup:", sum(!is.na(roads_named$road)), "\n")

# ============================================================
# DISTANCE FROM EACH SAMPLE TO NEAREST WINTER ROAD
# Memory-safe nearest-feature search in metres — see 00b_road_layers.R.
# ============================================================

dat_sf <- dat_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

dat_sf$dist_road_km <- nearest_road_km(dat_sf, roads_named)

# NOTE: (pcid, date) is NOT a unique row key — 18 individuals were
# sampled twice on the same day (repeat pellets), so a many-to-many
# join on (pcid, date) alone silently cross-products those rows,
# inflating n. Adding tz_sample_number makes the key unique per row
# (same fix as 04_climate_road_analyses.R and 10_supplementary_figures.R).
dat_roads <- dat_sf %>%
  st_drop_geometry() %>%
  select(pcid, date, tz_sample_number, dist_road_km) %>%
  right_join(dat, by = c("pcid", "date", "tz_sample_number"), relationship = "one-to-one") %>%
  left_join(
    roads_open %>%
      group_by(year) %>%
      summarise(road_days_open = mean(days_open, na.rm = TRUE), .groups = "drop"),
    by = "year"
  ) %>%
  mutate(
    road_days_open = replace_na(road_days_open, 0),
    survey         = factor(survey),
    year_f         = factor(year)
  )

# ============================================================
# 1) DESCRIPTIVE: ROAD DISTANCE BY YEAR × SURVEY
# ============================================================

p_road_dist <- ggplot(
  dat_roads %>% filter(is.finite(dist_road_km)),
  aes(x = factor(year), y = dist_road_km, colour = survey)
) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
              alpha = 0.5, size = 1.8) +
  geom_boxplot(aes(group = interaction(year, survey)),
               position = position_dodge(width = 0.6),
               outlier.shape = NA, fill = NA, linewidth = 0.7, width = 0.4) +
  scale_colour_manual(values = survey_colours, na.value = "grey50") +
  scale_y_log10() +
  labs(x = "Year", y = "Distance to nearest winter road (km, log scale)",
       colour = "Survey/Range",
       title = "Winter road proximity by year and survey") +
  theme_pub

ggsave(file.path(outdir, "Fig_RoadDistance_by_Year_Survey.png"),
       p_road_dist, width = 8, height = 5, dpi = 300)

write_csv(
  dat_roads %>%
    filter(is.finite(dist_road_km)) %>%
    group_by(year, survey) %>%
    summarise(
      n             = n(),
      dist_mean_km  = mean(dist_road_km),
      dist_median_km= median(dist_road_km),
      dist_sd_km    = sd(dist_road_km),
      .groups       = "drop"
    ),
  file.path(outdir, "Summary_RoadDistance_by_Year_Survey.csv")
)

# ============================================================
# 2) GAMs — ROAD PROXIMITY EFFECTS ON CORTISOL AND T3
# ============================================================

# --- 2A) Cortisol ---

dat_cort_road <- dat_roads %>%
  filter(is.finite(dist_road_km), is.finite(cortisol), !is.na(julian))

m_cort_road <- gam(
  cortisol ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(log(dist_road_km), k = 5, bs = "tp") +
    s(julian, k = 5, bs = "tp"),
  data   = dat_cort_road,
  family = Gamma(link = "log"),
  method = "REML"
)

m_cort_noroad <- gam(
  cortisol ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(julian, k = 8, bs = "tp"),
  data   = dat_cort_road,
  family = Gamma(link = "log"),
  method = "REML"
)

cat("Cortisol — road model selection (AIC, BIC):\n")
aic_cort_road <- AIC(m_cort_noroad, m_cort_road)
bic_cort_road <- BIC(m_cort_noroad, m_cort_road)
print(aic_cort_road); print(bic_cort_road)
capture.output(print(aic_cort_road), print(bic_cort_road),
               file = file.path(outdir, "GAM_AIC_03_cortisol_road.txt"))

summary(m_cort_road)
capture.output(summary(m_cort_road), file = file.path(outdir, "GAM_summary_03_cortisol_road.txt"))

png(file.path(outdir, "Fig_GAM_Cortisol_RoadDistance.png"), 2000, 1400, res = 300)
plot(m_cort_road, pages = 1)
dev.off()

cort_road_year_contrasts <- emmeans(m_cort_road, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(cort_road_year_contrasts), file = file.path(outdir, "GAM_contrasts_03_cortisol_road_year.txt"))

# --- 2B) T3 ---

dat_t3_road <- dat_roads %>%
  filter(is.finite(dist_road_km), is.finite(t3), t3 > 0, !is.na(julian))

m_t3_road <- gam(
  t3 ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(log(dist_road_km), k = 5, bs = "tp") +
    s(julian,       k = 5, bs = "tp"),
  data   = dat_t3_road,
  family = Gamma(link = "log"),
  method = "REML"
)

m_t3_noroad <- gam(
  t3 ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(julian, k = 8, bs = "tp"),
  data   = dat_t3_road,
  family = Gamma(link = "log"),
  method = "REML"
)

cat("T3 — road model selection (AIC, BIC):\n")
aic_t3_road <- AIC(m_t3_noroad, m_t3_road)
bic_t3_road <- BIC(m_t3_noroad, m_t3_road)
print(aic_t3_road); print(bic_t3_road)
capture.output(print(aic_t3_road), print(bic_t3_road),
               file = file.path(outdir, "GAM_AIC_03_t3_road.txt"))

summary(m_t3_road)
capture.output(summary(m_t3_road), file = file.path(outdir, "GAM_summary_03_t3_road.txt"))

png(file.path(outdir, "Fig_GAM_T3_RoadDistance.png"), 2000, 1400, res = 300)
plot(m_t3_road, pages = 1)
dev.off()

t3_road_year_contrasts <- emmeans(m_t3_road, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(t3_road_year_contrasts), file = file.path(outdir, "GAM_contrasts_03_t3_road_year.txt"))

# --- 2C) Road openness effect ---
# road_days_open is a year-level aggregate (one value per year), so it is
# perfectly collinear with year_f. We therefore drop year_f and use
# road_days_open as a *linear* parametric predictor. This tests whether the
# inter-annual variation in road openness predicts cortisol, controlling for
# sex/preg group, survey range (RE), and day-of-year.

dat_open_cort <- dat_roads %>%
  filter(is.finite(cortisol), is.finite(road_days_open), !is.na(julian))

m_cort_open <- gam(
  cortisol ~
    s(survey, bs = "re") +
    road_days_open +          # linear — only 4 unique values (one per year)
    preg_group +
    s(julian, k = 8, bs = "tp"),
  data   = dat_open_cort,
  family = Gamma(link = "log"),
  method = "REML"
)

summary(m_cort_open)
capture.output(summary(m_cort_open), file = file.path(outdir, "GAM_summary_03_cortisol_roadopen.txt"))

# Also test T3 ~ road openness
dat_open_t3 <- dat_roads %>%
  filter(is.finite(t3), t3 > 0, is.finite(road_days_open), !is.na(julian))

m_t3_open <- gam(
  t3 ~
    s(survey, bs = "re") +
    road_days_open +
    preg_group +
    s(julian, k = 8, bs = "tp"),
  data   = dat_open_t3,
  family = Gamma(link = "log"),
  method = "REML"
)

summary(m_t3_open)
capture.output(summary(m_t3_open), file = file.path(outdir, "GAM_summary_03_t3_roadopen.txt"))

# ============================================================
# 2D) SPATIAL AUTOCORRELATION DIAGNOSTICS (variogram on residuals)
# Purpose: check whether residual spatial autocorrelation could be
#   inflating the effective sample size (and therefore the significance)
#   of the road-distance smooths — flagged as a concern especially for
#   FT3~road, where the road-distance effect is large (ΔAIC = 186.7).
#
# Some samples share IDENTICAL GPS coordinates (repeat visits to the
# same collection site). A raw per-sample variogram would treat those
# as separate zero-distance observations and read as pure nugget noise
# rather than genuine small-scale variation, so residuals are averaged
# to one value per unique sampling location before the variogram is
# computed. If the variogram shows real spatial structure (sill clearly
# above the nugget, with a non-trivial range), the road models are
# refit with an added s(x, y) smooth on the projected coordinates, and
# the road-distance effect is re-checked with that term included.
# ============================================================

if (!requireNamespace("gstat", quietly = TRUE)) install.packages("gstat", repos = "https://cloud.r-project.org")
library(gstat)

# Build sf points matching each model's row order exactly.
#
# CORRECTION: the previous comment here claimed these coordinates were in
# "a projected, metres-based CRS". They were not — roads_named was read as
# EPSG:4326, so x_proj/y_proj (and therefore the variogram distance axis)
# have always been in DEGREES, not metres. The variogram plots are still
# valid for judging whether spatial structure exists, but their x-axis
# label should read degrees; treat the fitted "range" values as relative,
# not as metres.
#
# The CRS is now stated explicitly rather than inherited from
# st_crs(roads_named), so that changing the road layer can never silently
# re-project the spatial smooths and alter every model in the project.
# Road DISTANCES are computed separately, in metres (00b_road_layers.R).
SPATIAL_SMOOTH_CRS <- 4326

sf_cort_road <- dat_cort_road %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(SPATIAL_SMOOTH_CRS)
sf_t3_road <- dat_t3_road %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(SPATIAL_SMOOTH_CRS)

# Reusable helper: attach a residual vector to an sf points object,
# collapse to one averaged residual per unique (rounded) location,
# compute + plot + save the empirical variogram, and fit a spherical
# variogram model to summarise nugget vs. sill vs. range.
check_spatial_autocorr <- function(resid_vec, data_sf, label) {
  stopifnot(nrow(data_sf) == length(resid_vec))
  data_sf$resid <- resid_vec

  coords <- st_coordinates(data_sf)
  data_sf$x_round <- round(coords[, 1])   # nearest metre — collapses
  data_sf$y_round <- round(coords[, 2])   # floating-point-only duplicates

  site_resid <- data_sf %>%
    st_drop_geometry() %>%
    group_by(x_round, y_round) %>%
    summarise(resid = mean(resid), n_samples = n(), .groups = "drop")

  cat(sprintf("\n%s: %d rows collapsed to %d unique sampling locations for the variogram\n",
              label, nrow(data_sf), nrow(site_resid)))

  site_sf <- st_as_sf(site_resid, coords = c("x_round", "y_round"), crs = st_crs(data_sf))

  vgram <- gstat::variogram(resid ~ 1, data = site_sf)

  png(file.path(outdir, paste0("Variogram_", label, "_RoadModel.png")), 1600, 1200, res = 200)
  print(plot(vgram, main = paste0(label, " — variogram of road-model deviance residuals"),
             xlab = "Distance (m)", ylab = "Semivariance"))
  dev.off()

  vgram_fit <- tryCatch(gstat::fit.variogram(vgram, gstat::vgm("Sph")), error = function(e) NULL)
  if (!is.null(vgram_fit)) {
    cat(sprintf("%s variogram model (nugget / partial sill / range):\n", label))
    print(vgram_fit)
    capture.output(print(vgram_fit),
                   file = file.path(outdir, paste0("Variogram_fit_", label, "_RoadModel.txt")))
  } else {
    message(label, ": variogram model fit did not converge — inspect the empirical variogram plot directly.")
  }
  list(vgram = vgram, vgram_fit = vgram_fit)
}

vres_cort <- check_spatial_autocorr(residuals(m_cort_road, type = "deviance"), sf_cort_road, "FCM")
vres_t3   <- check_spatial_autocorr(residuals(m_t3_road,   type = "deviance"), sf_t3_road,   "FT3")

# ---- Refit with s(x, y) added, and re-check whether dist_road_km survives ----

dat_cort_road_xy <- dat_cort_road %>%
  bind_cols(st_coordinates(sf_cort_road) %>% as.data.frame() %>% rename(x_proj = X, y_proj = Y))
dat_t3_road_xy <- dat_t3_road %>%
  bind_cols(st_coordinates(sf_t3_road) %>% as.data.frame() %>% rename(x_proj = X, y_proj = Y))

m_cort_road_xy <- gam(
  cortisol ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(log(dist_road_km), k = 5, bs = "tp") +
    s(julian,       k = 5, bs = "tp") +
    s(x_proj, y_proj, bs = "tp", k = 30),
  data   = dat_cort_road_xy,
  family = Gamma(link = "log"),
  method = "REML"
)

plot(m_cort_road_xy, pages =1, scheme=2)

m_t3_road_xy <- gam(
  t3 ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(log(dist_road_km), k = 5, bs = "tp") +
    s(julian,       k = 5, bs = "tp") +
    s(x_proj, y_proj, bs = "tp", k = 30),
  data   = dat_t3_road_xy,
  family = Gamma(link = "log"),
  method = "REML"
)

plot(m_t3_road_xy, pages =1, scheme=2)

cat("\nRoad models with s(x, y) added — AIC/BIC comparison against the originally reported (no-spatial-term) models:\n")
aic_cort_xy <- AIC(m_cort_road, m_cort_road_xy)
aic_t3_xy   <- AIC(m_t3_road,   m_t3_road_xy)
bic_cort_xy <- BIC(m_cort_road, m_cort_road_xy)
bic_t3_xy   <- BIC(m_t3_road,   m_t3_road_xy)
print(aic_cort_xy); print(bic_cort_xy)
print(aic_t3_xy); print(bic_t3_xy)
capture.output(print(aic_cort_xy), print(bic_cort_xy), print(aic_t3_xy), print(bic_t3_xy),
               file = file.path(outdir, "GAM_AIC_03_spatial_xy_comparison.txt"))

cat("\nDoes dist_road_km remain significant once s(x, y) is included?\n")
print(summary(m_cort_road_xy)$s.table)
print(summary(m_t3_road_xy)$s.table)
capture.output(summary(m_cort_road_xy), file = file.path(outdir, "GAM_summary_03_cortisol_road_xy.txt"))
capture.output(summary(m_t3_road_xy),   file = file.path(outdir, "GAM_summary_03_t3_road_xy.txt"))

# Re-check the variogram on the xy-augmented models' residuals — if
# s(x, y) has done its job, this should look flatter than the original.
vres_cort_after <- check_spatial_autocorr(residuals(m_cort_road_xy, type = "deviance"), sf_cort_road, "FCM_afterXY")
vres_t3_after   <- check_spatial_autocorr(residuals(m_t3_road_xy,   type = "deviance"), sf_t3_road,   "FT3_afterXY")

message("Part 3D (spatial autocorrelation diagnostics) complete — outputs in: ", outdir)

# ============================================================
# 2E) PRIMARY MODEL DECISION — spatially-adjusted model reported
# Project decision: for landscape covariates (road, mine, wolf, land
# cover, population — see 05/06/07/08), the s(x, y)-augmented model is
# reported as the PRIMARY model, regardless of BIC. Rationale: AIC
# favours adding s(x, y) in every one of these covariates; BIC opposes
# it for FCM but agrees for FT3. Given the known risk that repeat
# visits to a small number of GPS locations can inflate the effective
# sample size behind these smooths, we treat the spatial term as a
# necessary correction rather than an optional complexity penalty, and
# report the s(x, y)-augmented model throughout — for both hormones,
# not only where AIC and BIC agree. This overrides the general
# AIC/BIC-discordance policy stated in Methods for THESE FIVE
# COVARIATES specifically; all other reported models (climate, road
# openness, etc.) are unaffected.
#
# m_cort_road_xy / m_t3_road_xy (built above) are therefore the primary
# reported road models. The outputs below OVERWRITE the original
# (non-spatial) GAM_summary/Fig/contrasts files from Part 2A/2B so that
# every file on disk reflects the model actually reported in the
# manuscript.
# ============================================================

summary(m_cort_road_xy)
capture.output(summary(m_cort_road_xy), file = file.path(outdir, "GAM_summary_03_cortisol_road.txt"))

png(file.path(outdir, "Fig_GAM_Cortisol_RoadDistance.png"), 2000, 1400, res = 300)
plot(m_cort_road_xy, pages = 1)
dev.off()

cort_road_year_contrasts <- emmeans(m_cort_road_xy, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(cort_road_year_contrasts), file = file.path(outdir, "GAM_contrasts_03_cortisol_road_year.txt"))

summary(m_t3_road_xy)
capture.output(summary(m_t3_road_xy), file = file.path(outdir, "GAM_summary_03_t3_road.txt"))

png(file.path(outdir, "Fig_GAM_T3_RoadDistance.png"), 2000, 1400, res = 300)
plot(m_t3_road_xy, pages = 1)
dev.off()

t3_road_year_contrasts <- emmeans(m_t3_road_xy, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(t3_road_year_contrasts), file = file.path(outdir, "GAM_contrasts_03_t3_road_year.txt"))

message("Part 2E (primary model = spatially-adjusted) complete — GAM_summary_03_cortisol_road.txt / GAM_summary_03_t3_road.txt now reflect the s(x,y)-augmented models.")

message("Part 3 complete — outputs in: ", outdir)

