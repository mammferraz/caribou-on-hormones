# ============================================================
# 07_mines_analyses.R
#
# PART 7 — MINING PROXIMITY
# Purpose:
#   • Map mine proximity raster
#   • Extract distance-to-nearest-mine at sample locations
#   • Test effects on cortisol and T3 by year and survey
# ============================================================

source("00_data_processing.R")

outdir <- "caribou_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# I/O
# ============================================================

infile_mines <- "mines.tif"

# ============================================================
# LOAD RASTER & EXTRACT AT SAMPLE LOCATIONS
# ============================================================

mines_r <- terra::rast(infile_mines)

pts_sf_mn <- dat_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(mines_r))

# mines.tif encodes mine presence/density (0 = no mine, 1–8 = mine activity).
# Compute distance from each cell to the nearest mine cell, then extract.
# Crop to a generous buffer around sample locations first for speed.
buf_mn     <- st_buffer(pts_sf_mn, dist = 500000)   # 500 km buffer
mines_crop <- terra::crop(mines_r, vect(buf_mn))

# Binary: mine cells = 1, all others = NA
mine_binary <- mines_crop
mine_binary[mine_binary == 0] <- NA

# Distance (metres) from every cell to nearest mine cell
cat("Computing distance to nearest mine raster (may take ~1 min)...\n")
dist_rast <- terra::distance(mine_binary)

# Extract and convert to km
mines_vals <- terra::extract(dist_rast, vect(pts_sf_mn))[, 2] / 1000

dat_mines <- dat_map %>%
  mutate(dist_mine_km = mines_vals) %>%
  filter(is.finite(dist_mine_km), !is.na(preg_group)) %>%
  mutate(survey = factor(survey), year_f = factor(year),
         # +0.1 km offset guards against log(0) if any sample sits exactly on
         # a mine cell (dist_mine_km = 0); also matches the transform used in
         # 09_manuscript_figures.R so the reported model and the figure agree.
         log_dist_mine = log(dist_mine_km + 0.1))

cat("Distance to nearest mine (km):\n")
print(summary(dat_mines$dist_mine_km))

# ============================================================
# 1) MAP: MINE DISTANCE + CARIBOU LOCATIONS
# ============================================================

pts_sf_proj <- dat_mines %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(mines_r))

buf           <- st_buffer(pts_sf_proj, dist = 50000)
dist_crop     <- terra::crop(dist_rast / 1000, vect(buf))  # km
mines_df      <- as.data.frame(dist_crop, xy = TRUE, na.rm = FALSE)
names(mines_df) <- c("x", "y", "dist_mine_km")

# Mine cluster centroids — group adjacent non-zero cells into patches,
# then take the centroid of each patch as the "mine location" symbol
mine_binary_map <- terra::crop(mine_binary, vect(buf))
mine_patches    <- terra::patches(mine_binary_map, directions = 8)
mine_polys      <- terra::as.polygons(mine_patches)   # one polygon per patch
mine_centers    <- terra::centroids(mine_polys) |> st_as_sf()

p_mines_map <- ggplot() +
  geom_raster(data = mines_df, aes(x, y, fill = dist_mine_km)) +
  geom_sf(data = pts_sf_proj,
          aes(colour = survey, shape = preg_group),
          size = 1.8, alpha = 0.95) +
  geom_sf(data = mine_centers,
          shape = 23, size = 3, fill = "yellow", colour = "black", stroke = 0.7) +
  scale_fill_viridis_c(name = "Distance to\nnearest mine (km)", option = "magma",
                        direction = -1, na.value = "white") +
  scale_colour_manual(values = survey_colours, na.value = "grey50") +
  scale_shape_manual(values = preg_shapes) +
  coord_sf(expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.35) +
  labs(title  = "Distance to nearest mine and caribou sampling locations",
       x = "Easting", y = "Northing",
       colour = "Survey/Range", shape = "Group") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(outdir, "Fig_MineDistance_Map.png"),
       p_mines_map, width = 8, height = 6, dpi = 300)


# ============================================================
# 2) MINE DISTANCE BY YEAR × SURVEY
# ============================================================

p_mine_year <- ggplot(
  dat_mines,
  aes(x = factor(year), y = dist_mine_km, colour = survey)
) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
              alpha = 0.5, size = 1.8) +
  geom_boxplot(aes(group = interaction(year, survey)),
               position = position_dodge(0.6),
               outlier.shape = NA, fill = NA, linewidth = 0.7, width = 0.4) +
  scale_colour_manual(values = survey_colours, na.value = "grey50") +
  labs(x = "Year", y = "Distance to nearest mine (km)",
       colour = "Survey/Range",
       title  = "Mining proximity at caribou locations by year") +
  theme_pub

ggsave(file.path(outdir, "Fig_MineDistance_by_Year_Survey.png"),
       p_mine_year, width = 8, height = 5, dpi = 300)


# ============================================================
# 3) GAMs — MINE PROXIMITY EFFECTS
# ============================================================

# --- 3A) Cortisol ---

dat_mines_cort <- dat_mines %>%
  filter(is.finite(cortisol), !is.na(julian))

m_cort_mine <- gam(
  cortisol ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(log_dist_mine, k = 6, bs = "tp") +
    s(julian,         k = 8, bs = "tp"),
  data   = dat_mines_cort,
  family = Gamma(link = "log"),
  method = "REML"
)

summary(m_cort_mine)
capture.output(summary(m_cort_mine), file = file.path(outdir, "GAM_summary_07_cortisol_mine.txt"))
gam.check(m_cort_mine)

png(file.path(outdir, "Fig_GAM_Cortisol_MineProximity.png"), 2000, 1400, res = 300)
plot(m_cort_mine, pages = 1)
dev.off()

cort_mine_year_contrasts <- emmeans(m_cort_mine, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(cort_mine_year_contrasts), file = file.path(outdir, "GAM_contrasts_07_cortisol_mine_year.txt"))

# --- 3B) T3 ---

dat_mines_t3 <- dat_mines %>%
  filter(is.finite(t3), t3 > 0, !is.na(julian))

m_t3_mine <- gam(
  t3 ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(log_dist_mine, k = 6, bs = "tp") +
    s(julian,         k = 8, bs = "tp"),
  data   = dat_mines_t3,
  family = Gamma(link = "log"),
  method = "REML"
)

summary(m_t3_mine)
capture.output(summary(m_t3_mine), file = file.path(outdir, "GAM_summary_07_t3_mine.txt"))

png(file.path(outdir, "Fig_GAM_T3_MineProximity.png"), 2000, 1400, res = 300)
plot(m_t3_mine, pages = 1)
dev.off()

t3_mine_year_contrasts <- emmeans(m_t3_mine, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(t3_mine_year_contrasts), file = file.path(outdir, "GAM_contrasts_07_t3_mine_year.txt"))

# ============================================================
# 4) SPATIAL AUTOCORRELATION DIAGNOSTICS (variogram on residuals)
# Purpose: same diagnostic applied to the road-distance models in
#   03_winter_road_analyses.R, now extended to mine proximity — check
#   whether residual spatial autocorrelation could be inflating the
#   apparent significance of the log_dist_mine smooth.
#
# Some samples share IDENTICAL GPS coordinates (repeat visits to the
# same collection site). A raw per-sample variogram would treat those
# as separate zero-distance observations and read as pure nugget noise
# rather than genuine small-scale variation, so residuals are averaged
# to one value per unique sampling location before the variogram is
# computed. If the variogram shows real spatial structure (sill clearly
# above the nugget, with a non-trivial range), the mine models are
# refit with an added s(x, y) smooth on the projected coordinates, and
# the log_dist_mine effect is re-checked with that term included.
# ============================================================

if (!requireNamespace("gstat", quietly = TRUE)) install.packages("gstat", repos = "https://cloud.r-project.org")
library(gstat)

sf_cort_mine <- dat_mines_cort %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(mines_r))
sf_t3_mine <- dat_mines_t3 %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(mines_r))

check_spatial_autocorr_mine <- function(resid_vec, data_sf, label) {
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

  png(file.path(outdir, paste0("Variogram_", label, "_MineModel.png")), 1600, 1200, res = 200)
  print(plot(vgram, main = paste0(label, " — variogram of mine-model deviance residuals"),
             xlab = "Distance (m)", ylab = "Semivariance"))
  dev.off()

  vgram_fit <- tryCatch(gstat::fit.variogram(vgram, gstat::vgm("Sph")), error = function(e) NULL)
  if (!is.null(vgram_fit)) {
    cat(sprintf("%s variogram model (nugget / partial sill / range):\n", label))
    print(vgram_fit)
    capture.output(print(vgram_fit),
                   file = file.path(outdir, paste0("Variogram_fit_", label, "_MineModel.txt")))
  } else {
    message(label, ": variogram model fit did not converge — inspect the empirical variogram plot directly.")
  }
  list(vgram = vgram, vgram_fit = vgram_fit)
}

vres_cort_mine <- check_spatial_autocorr_mine(residuals(m_cort_mine, type = "deviance"), sf_cort_mine, "FCM")
vres_t3_mine   <- check_spatial_autocorr_mine(residuals(m_t3_mine,   type = "deviance"), sf_t3_mine,   "FT3")

# ---- Refit with s(x, y) added, and re-check whether log_dist_mine survives ----

dat_cort_mine_xy <- dat_mines_cort %>%
  bind_cols(st_coordinates(sf_cort_mine) %>% as.data.frame() %>% rename(x_proj = X, y_proj = Y))
dat_t3_mine_xy <- dat_mines_t3 %>%
  bind_cols(st_coordinates(sf_t3_mine) %>% as.data.frame() %>% rename(x_proj = X, y_proj = Y))

m_cort_mine_xy <- gam(
  cortisol ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(log_dist_mine, k = 6, bs = "tp") +
    s(julian,        k = 8, bs = "tp") +
    s(x_proj, y_proj, bs = "tp", k = 30),
  data   = dat_cort_mine_xy,
  family = Gamma(link = "log"),
  method = "REML"
)

plot(m_cort_mine_xy, pages = 1, scheme = 2)

m_t3_mine_xy <- gam(
  t3 ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(log_dist_mine, k = 6, bs = "tp") +
    s(julian,        k = 8, bs = "tp") +
    s(x_proj, y_proj, bs = "tp", k = 30),
  data   = dat_t3_mine_xy,
  family = Gamma(link = "log"),
  method = "REML"
)

plot(m_t3_mine_xy, pages = 1, scheme = 2)

cat("\nMine models with s(x, y) added — AIC/BIC comparison against the originally reported (no-spatial-term) models:\n")
aic_cort_mine_xy <- AIC(m_cort_mine, m_cort_mine_xy)
aic_t3_mine_xy   <- AIC(m_t3_mine,   m_t3_mine_xy)
bic_cort_mine_xy <- BIC(m_cort_mine, m_cort_mine_xy)
bic_t3_mine_xy   <- BIC(m_t3_mine,   m_t3_mine_xy)
print(aic_cort_mine_xy); print(bic_cort_mine_xy)
print(aic_t3_mine_xy); print(bic_t3_mine_xy)
capture.output(print(aic_cort_mine_xy), print(bic_cort_mine_xy), print(aic_t3_mine_xy), print(bic_t3_mine_xy),
               file = file.path(outdir, "GAM_AIC_07_spatial_xy_comparison.txt"))

cat("\nDoes log_dist_mine remain significant once s(x, y) is included?\n")
print(summary(m_cort_mine_xy)$s.table)
print(summary(m_t3_mine_xy)$s.table)
capture.output(summary(m_cort_mine_xy), file = file.path(outdir, "GAM_summary_07_cortisol_mine_xy.txt"))
capture.output(summary(m_t3_mine_xy),   file = file.path(outdir, "GAM_summary_07_t3_mine_xy.txt"))

# Re-check the variogram on the xy-augmented models' residuals — if
# s(x, y) has done its job, this should look flatter than the original.
vres_cort_mine_after <- check_spatial_autocorr_mine(residuals(m_cort_mine_xy, type = "deviance"), sf_cort_mine, "FCM_afterXY")
vres_t3_mine_after   <- check_spatial_autocorr_mine(residuals(m_t3_mine_xy,   type = "deviance"), sf_t3_mine,   "FT3_afterXY")

message("Part 7 spatial autocorrelation diagnostics complete — outputs in: ", outdir)

# ============================================================
# 5) PRIMARY MODEL DECISION — spatially-adjusted model reported
# Project decision: for landscape covariates (road, mine, wolf, land
# cover, population), the s(x, y)-augmented model is reported as the
# PRIMARY model, regardless of BIC (see 03_winter_road_analyses.R
# Part 2E for full rationale). m_cort_mine_xy / m_t3_mine_xy are
# therefore the primary reported mine models. The outputs below
# OVERWRITE the original (non-spatial) GAM_summary/Fig/contrasts files
# from Part 3 so every file on disk reflects the model actually
# reported in the manuscript.
# ============================================================

summary(m_cort_mine_xy)
capture.output(summary(m_cort_mine_xy), file = file.path(outdir, "GAM_summary_07_cortisol_mine.txt"))

png(file.path(outdir, "Fig_GAM_Cortisol_MineProximity.png"), 2000, 1400, res = 300)
plot(m_cort_mine_xy, pages = 1)
dev.off()

cort_mine_year_contrasts <- emmeans(m_cort_mine_xy, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(cort_mine_year_contrasts), file = file.path(outdir, "GAM_contrasts_07_cortisol_mine_year.txt"))

summary(m_t3_mine_xy)
capture.output(summary(m_t3_mine_xy), file = file.path(outdir, "GAM_summary_07_t3_mine.txt"))

png(file.path(outdir, "Fig_GAM_T3_MineProximity.png"), 2000, 1400, res = 300)
plot(m_t3_mine_xy, pages = 1)
dev.off()

t3_mine_year_contrasts <- emmeans(m_t3_mine_xy, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(t3_mine_year_contrasts), file = file.path(outdir, "GAM_contrasts_07_t3_mine_year.txt"))

message("Part 5 (primary model = spatially-adjusted) complete — GAM_summary_07_cortisol_mine.txt / GAM_summary_07_t3_mine.txt now reflect the s(x,y)-augmented models.")

message("Part 7 complete — outputs in: ", outdir)
