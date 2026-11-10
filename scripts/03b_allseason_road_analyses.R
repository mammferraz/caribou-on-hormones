# ============================================================
# 03b_allseason_road_analyses.R
#
# PART 3B — ALL-SEASON ROAD PROXIMITY
# Purpose:
#   • Compute distance from each sample to nearest all-season road
#   • Test effects on cortisol and T3
#   • Primary axis: year; survey as random effect and colour
#
# Data: MNR_Road_Segments2.geojson
#   Downloaded 2026-07-12 from:
#   https://geohub.lio.gov.on.ca/datasets/lio::mnr-road-segments
#   Full provincial dataset (not pre-filtered).
#
# ALL-SEASON ROAD DEFINITION — CORRECTED 2026-07-26
#   Previously this script defined all-season roads as
#   NATIONAL_ROAD_CLASS == "Local / Strata" ONLY (8,648 of 517,836
#   segments, 1.7%), excluding "Resource / Recreation" on the stated
#   assumption that those were "largely seasonal/closed forest management
#   roads". That assumption was wrong: of the 508,220 Resource/Recreation
#   segments, 343,960 have STATUS = "Open" and 312,616 have
#   PASSABLE_IND = "Yes". Local/Strata is a narrow ownership class, not a
#   road network, and produced scattered fragments rather than the dense
#   network visible in the region.
#
#   All-season roads are now defined as ALL NON-WINTER segments
#   (Resource/Recreation + Local/Strata = 516,868), loaded via the shared
#   helper in 00b_road_layers.R. Winter roads are handled in 03.
#
#   NOTE: every all-season road result predating 2026-07-26 (distances,
#   GAM summaries, SupFig 7a-b, SupFig 11k-l, SupFig 12/13, and the
#   full-model term in 11/12) was computed on the old 8,648-segment
#   definition and must be regenerated.
# ============================================================

library(tidyverse)
library(sf)
library(mgcv)
library(emmeans)
library(janitor)

source("00_data_processing.R")
source("00b_road_layers.R")

outdir <- "caribou_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# LOAD SPATIAL ROADS (all non-winter segments)
# Cached to caribou_outputs/roads_allseason.gpkg on first run.
# ============================================================

roads_all_sf <- road_layer("allseason")
cat("All-season road features:", nrow(roads_all_sf), "\n")

# ============================================================
# DISTANCE FROM EACH SAMPLE TO NEAREST ALL-SEASON ROAD
# Uses nearest_road_km(): spatially indexed, O(n) memory. The old
# st_distance(points, roads) built a full n x m matrix, which with
# ~517k road segments is ~5.7 GB and exceeds R's vector memory limit.
# ============================================================

dat_sf <- dat_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

cat("Computing distances to nearest all-season road...\n")
dat_sf$dist_allroad_km <- nearest_road_km(dat_sf, roads_all_sf)

dat_allroads <- dat_sf %>%
  st_drop_geometry() %>%
  select(pcid, date, tz_sample_number, dist_allroad_km) %>%
  right_join(dat, by = c("pcid", "date", "tz_sample_number"), relationship = "one-to-one") %>%
  mutate(
    survey = factor(survey),
    year_f = factor(year)
  )

cat("Road distance summary (km):\n")
print(summary(dat_allroads$dist_allroad_km))

write_csv(
  dat_allroads %>%
    filter(is.finite(dist_allroad_km)) %>%
    select(pcid, date, tz_sample_number, survey, year, dist_allroad_km),
  file.path(outdir, "AllSeasonRoad_distances.csv")
)

# ============================================================
# 1) DESCRIPTIVE: ROAD DISTANCE BY YEAR × SURVEY
# ============================================================

p_allroad_dist <- ggplot(
  dat_allroads %>% filter(is.finite(dist_allroad_km)),
  aes(x = factor(year), y = dist_allroad_km, colour = survey)
) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
              alpha = 0.5, size = 1.8) +
  geom_boxplot(aes(group = interaction(year, survey)),
               position = position_dodge(width = 0.6),
               outlier.shape = NA, fill = NA, linewidth = 0.7, width = 0.4) +
  scale_colour_manual(values = survey_colours, na.value = "grey50") +
  scale_y_log10() +
  labs(x = "Year", y = "Distance to nearest all-season road (km, log scale)",
       colour = "Survey/Range",
       title = "All-season road proximity by year and survey") +
  theme_pub

ggsave(file.path(outdir, "Fig_AllSeasonRoad_Distance_by_Year_Survey.png"),
       p_allroad_dist, width = 8, height = 5, dpi = 300)

write_csv(
  dat_allroads %>%
    filter(is.finite(dist_allroad_km)) %>%
    group_by(year, survey) %>%
    summarise(
      n              = n(),
      dist_mean_km   = mean(dist_allroad_km),
      dist_median_km = median(dist_allroad_km),
      dist_sd_km     = sd(dist_allroad_km),
      .groups        = "drop"
    ),
  file.path(outdir, "Summary_AllSeasonRoad_Distance_by_Year_Survey.csv")
)

# ============================================================
# 2) GAMs — ALL-SEASON ROAD PROXIMITY EFFECTS
# ============================================================

# --- 2A) Cortisol ---

dat_cort_allroad <- dat_allroads %>%
  filter(is.finite(dist_allroad_km), is.finite(cortisol), !is.na(julian))

m_cort_allroad <- gam(
  cortisol ~
    s(survey,              bs = "re") +
    year_f +
    preg_group +
    s(log(dist_allroad_km), k = 5, bs = "tp") +
    s(julian,               k = 5, bs = "tp"),
  data   = dat_cort_allroad,
  family = Gamma(link = "log"),
  method = "REML"
)

m_cort_noallroad <- gam(
  cortisol ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(julian, k = 8, bs = "tp"),
  data   = dat_cort_allroad,
  family = Gamma(link = "log"),
  method = "REML"
)

cat("Cortisol — all-season road model selection (AIC, BIC):\n")
aic_cort_allroad <- AIC(m_cort_noallroad, m_cort_allroad)
bic_cort_allroad <- BIC(m_cort_noallroad, m_cort_allroad)
print(aic_cort_allroad); print(bic_cort_allroad)
capture.output(print(aic_cort_allroad), print(bic_cort_allroad),
               file = file.path(outdir, "GAM_AIC_03b_cortisol_allroad.txt"))

summary(m_cort_allroad)
capture.output(summary(m_cort_allroad),
               file = file.path(outdir, "GAM_summary_03b_cortisol_allroad.txt"))

png(file.path(outdir, "Fig_GAM_Cortisol_AllSeasonRoad.png"), 2000, 1400, res = 300)
plot(m_cort_allroad, pages = 1)
dev.off()

cort_allroad_year_contrasts <- emmeans(m_cort_allroad, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(cort_allroad_year_contrasts),
               file = file.path(outdir, "GAM_contrasts_03b_cortisol_allroad_year.txt"))

# --- 2B) T3 ---

dat_t3_allroad <- dat_allroads %>%
  filter(is.finite(dist_allroad_km), is.finite(t3), t3 > 0, !is.na(julian))

m_t3_allroad <- gam(
  t3 ~
    s(survey,               bs = "re") +
    year_f +
    preg_group +
    s(log(dist_allroad_km), k = 5, bs = "tp") +
    s(julian,               k = 5, bs = "tp"),
  data   = dat_t3_allroad,
  family = Gamma(link = "log"),
  method = "REML"
)

m_t3_noallroad <- gam(
  t3 ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(julian, k = 8, bs = "tp"),
  data   = dat_t3_allroad,
  family = Gamma(link = "log"),
  method = "REML"
)

cat("T3 — all-season road model selection (AIC, BIC):\n")
aic_t3_allroad <- AIC(m_t3_noallroad, m_t3_allroad)
bic_t3_allroad <- BIC(m_t3_noallroad, m_t3_allroad)
print(aic_t3_allroad); print(bic_t3_allroad)
capture.output(print(aic_t3_allroad), print(bic_t3_allroad),
               file = file.path(outdir, "GAM_AIC_03b_t3_allroad.txt"))

summary(m_t3_allroad)
capture.output(summary(m_t3_allroad),
               file = file.path(outdir, "GAM_summary_03b_t3_allroad.txt"))

png(file.path(outdir, "Fig_GAM_T3_AllSeasonRoad.png"), 2000, 1400, res = 300)
plot(m_t3_allroad, pages = 1)
dev.off()

t3_allroad_year_contrasts <- emmeans(m_t3_allroad, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(t3_allroad_year_contrasts),
               file = file.path(outdir, "GAM_contrasts_03b_t3_allroad_year.txt"))

# ============================================================
# 2C) SPATIAL AUTOCORRELATION CHECK
#   Same approach as 03: refit with s(x, y) and compare
# ============================================================

if (!requireNamespace("gstat", quietly = TRUE)) install.packages("gstat", repos = "https://cloud.r-project.org")
library(gstat)

# NOTE: the coordinates behind s(x_proj, y_proj) are held in EPSG:4326
# (degrees), which is what every spatial smooth in this project has always
# used. This is stated explicitly rather than inherited from
# st_crs(roads_all_sf) so that changing the road layer can never silently
# re-project the spatial smooths and alter every model in the project.
# Road DISTANCES are computed separately, in metres (see 00b_road_layers.R).
SPATIAL_SMOOTH_CRS <- 4326

sf_cort_allroad <- dat_cort_allroad %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(SPATIAL_SMOOTH_CRS)

sf_t3_allroad <- dat_t3_allroad %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(SPATIAL_SMOOTH_CRS)

check_spatial_autocorr <- function(resid_vec, data_sf, label) {
  stopifnot(nrow(data_sf) == length(resid_vec))
  data_sf$resid <- resid_vec
  coords <- st_coordinates(data_sf)
  data_sf$x_round <- round(coords[, 1])
  data_sf$y_round <- round(coords[, 2])
  site_resid <- data_sf %>%
    st_drop_geometry() %>%
    group_by(x_round, y_round) %>%
    summarise(resid = mean(resid), n_samples = n(), .groups = "drop")
  cat(sprintf("\n%s: %d rows → %d unique locations\n",
              label, nrow(data_sf), nrow(site_resid)))
  site_sf <- st_as_sf(site_resid, coords = c("x_round", "y_round"), crs = st_crs(data_sf))
  vgram <- gstat::variogram(resid ~ 1, data = site_sf)
  png(file.path(outdir, paste0("Variogram_", label, "_AllRoadModel.png")), 1600, 1200, res = 200)
  print(plot(vgram, main = paste0(label, " — variogram of all-season road model residuals"),
             xlab = "Distance (m)", ylab = "Semivariance"))
  dev.off()
  vgram_fit <- tryCatch(gstat::fit.variogram(vgram, gstat::vgm("Sph")), error = function(e) NULL)
  if (!is.null(vgram_fit)) {
    print(vgram_fit)
    capture.output(print(vgram_fit),
                   file = file.path(outdir, paste0("Variogram_fit_", label, "_AllRoadModel.txt")))
  }
  list(vgram = vgram, vgram_fit = vgram_fit)
}

vres_cort <- check_spatial_autocorr(residuals(m_cort_allroad, type = "deviance"), sf_cort_allroad, "FCM_allroad")
vres_t3   <- check_spatial_autocorr(residuals(m_t3_allroad,   type = "deviance"), sf_t3_allroad,   "FT3_allroad")

# Refit with s(x, y)
dat_cort_allroad_xy <- dat_cort_allroad %>%
  bind_cols(st_coordinates(sf_cort_allroad) %>% as.data.frame() %>% rename(x_proj = X, y_proj = Y))
dat_t3_allroad_xy <- dat_t3_allroad %>%
  bind_cols(st_coordinates(sf_t3_allroad) %>% as.data.frame() %>% rename(x_proj = X, y_proj = Y))

m_cort_allroad_xy <- gam(
  cortisol ~
    s(survey,               bs = "re") +
    year_f +
    preg_group +
    s(log(dist_allroad_km), k = 5, bs = "tp") +
    s(julian,               k = 5, bs = "tp") +
    s(x_proj, y_proj,       bs = "tp", k = 30),
  data   = dat_cort_allroad_xy,
  family = Gamma(link = "log"),
  method = "REML"
)

m_t3_allroad_xy <- gam(
  t3 ~
    s(survey,               bs = "re") +
    year_f +
    preg_group +
    s(log(dist_allroad_km), k = 5, bs = "tp") +
    s(julian,               k = 5, bs = "tp") +
    s(x_proj, y_proj,       bs = "tp", k = 30),
  data   = dat_t3_allroad_xy,
  family = Gamma(link = "log"),
  method = "REML"
)

cat("\nAll-season road models with s(x, y) — AIC/BIC vs. no spatial term:\n")
print(AIC(m_cort_allroad, m_cort_allroad_xy))
print(BIC(m_cort_allroad, m_cort_allroad_xy))
print(AIC(m_t3_allroad,   m_t3_allroad_xy))
print(BIC(m_t3_allroad,   m_t3_allroad_xy))

capture.output(
  print(AIC(m_cort_allroad, m_cort_allroad_xy)),
  print(BIC(m_cort_allroad, m_cort_allroad_xy)),
  print(AIC(m_t3_allroad,   m_t3_allroad_xy)),
  print(BIC(m_t3_allroad,   m_t3_allroad_xy)),
  file = file.path(outdir, "GAM_AIC_03b_spatial_xy_comparison.txt")
)

cat("\nDoes dist_allroad_km remain significant once s(x, y) is included?\n")
print(summary(m_cort_allroad_xy)$s.table)
print(summary(m_t3_allroad_xy)$s.table)
capture.output(summary(m_cort_allroad_xy), file = file.path(outdir, "GAM_summary_03b_cortisol_allroad_xy.txt"))
capture.output(summary(m_t3_allroad_xy),   file = file.path(outdir, "GAM_summary_03b_t3_allroad_xy.txt"))

# Post-spatial variograms
check_spatial_autocorr(residuals(m_cort_allroad_xy, type = "deviance"), sf_cort_allroad, "FCM_allroad_afterXY")
check_spatial_autocorr(residuals(m_t3_allroad_xy,   type = "deviance"), sf_t3_allroad,   "FT3_allroad_afterXY")

message("Part 3B complete — outputs in: ", outdir)
