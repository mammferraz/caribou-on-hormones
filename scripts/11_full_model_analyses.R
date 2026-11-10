# ============================================================
# 11_full_model_analyses.R
#
# PART 11 — FULL COMBINED MODEL (ALL LANDSCAPE + CLIMATE COVARIATES)
# Purpose:
#   • Requested check: fit ONE model per hormone containing all of
#     short-lag snow-on-ground, winter road distance, wolf habitat
#     suitability, mine distance, distance to nearest community, AND
#     land cover composition (PCA, 4 axes) simultaneously (plus the usual
#     survey/year/pregnancy/Julian day terms), fit BOTH without and
#     with a spatial s(x, y) term, and report the comparison.
#   • This is an additional robustness check requested after the
#     single-covariate spatial-autocorrelation work in 03/05/06/07/08 —
#     it is NOT a replacement for those individual analyses, which
#     remain the primary reported models. This script asks: do the
#     individual-covariate conclusions hold up when every landscape
#     covariate competes for the same variance simultaneously, and does
#     the spatial term still matter once all of them are in the model
#     together?
#   • Land cover is summarised by the same PCA used in 05/09/10
#     (LC_PC1-4, fit on the land-cover-complete sample before
#     intersecting with the other five covariates below). Including it
#     further shrinks the complete-case sample size — check the
#     reported n against 05/06/07/08's individual n's before
#     interpreting effect sizes.
# ============================================================

source("00_data_processing.R")
source("00b_road_layers.R")

outdir <- "caribou_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# I/O
# ============================================================

infile_climate <- "climatedaily_all_surveys.csv"
infile_roads        <- "MNR_Road_Segments.geojson"
infile_roads_all    <- "MNR_Road_Segments2.geojson"   # all classes; see 00b_road_layers.R
infile_lookup  <- "winter_roads_all2.csv"
infile_wolf      <- "wolfmean_t.tif"
infile_mines     <- "mines.tif"
infile_pop       <- "population_density.tif"
infile_wolfstack <- "wolfstack.tif"

# ============================================================
# 1) SHORT-LAG CLIMATE (snow-on-ground)
# Same lag-window logic as 02/04/09/10 (see those scripts for full
# rationale on per-sample nearest-station assignment, computed once in
# 00_data_processing.R as climate_group).
# FCM (cortisol) full model below uses a 2-day short-lag window
# (snow_on_ground_fcm2); FT3 full model keeps the 5-day window
# (snow_on_ground) — see 02_climate_analyses.R header for rationale.
# ============================================================

clim_all <- read_csv(infile_climate, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(date = as.Date(local_date), local_year = as.integer(local_year),
         month = month(date))

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

climate_lags_fcm2 <- sample_dates %>%
  rowwise() %>%
  mutate(lag = list(calc_lag_mean(date, clim_split[[climate_group]], 1, 2))) %>%
  ungroup() %>%
  unnest(lag)

# climate_lags is one row per unique (pcid, year, date) — joining without
# tz_sample_number is safe here because the RIGHT-hand table has no
# duplicate keys, so this left_join cannot multiply rows (unlike the
# road-distance join fixed elsewhere, where the duplicate-prone table
# was on the right).
clim_for_join <- climate_lags %>%
  select(pcid, year, date, snow_on_ground)

# FCM's 2-day window, suffixed _fcm2 to stay distinct when joined
# alongside the 5-day columns above.
clim_for_join_fcm2 <- climate_lags_fcm2 %>%
  select(pcid, year, date, snow_on_ground) %>%
  rename(snow_on_ground_fcm2 = snow_on_ground)

# ============================================================
# 2) WINTER ROAD DISTANCE (direct row-order extraction on dat_map,
# same approach already used for mine/wolf/population in 06/07/08 —
# avoids the join entirely, so the (pcid, date) duplicate-key issue
# fixed elsewhere for 03/04/10 never arises here)
# ============================================================

# Winter roads from the provincial layer via 00b_road_layers.R. The old
# MNR_Road_Segments.geojson is a 541-segment spatial clip of the full 968;
# 7 of the absent segments fall inside the sampling window.
roads_sf <- road_layer("winter") %>% clean_names()
road_lookup <- read_csv(infile_lookup, show_col_types = FALSE) %>% clean_names()
roads_named <- roads_sf %>% left_join(road_lookup, by = "ogf_id")
# Spatial-smooth coordinates are held in EPSG:4326 (degrees) — what every
# s(x, y) term in this project has always used. Stated explicitly rather
# than inherited from a road layer's CRS, so swapping road layers can never
# silently re-project the smooths. Road DISTANCES are separate, in metres.
SPATIAL_SMOOTH_CRS <- 4326


pts_sf_road <- dat_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
road_vals <- nearest_road_km(pts_sf_road, roads_named)

# ============================================================
# 2b) ALL-SEASON ROAD DISTANCE (same approach as winter road above,
#   using all non-winter segments of MNR_Road_Segments2.geojson)
#   NOTE: file is ~1.3 GB — loading takes ~30–60 s
# ============================================================

# ALL-SEASON ROADS — CORRECTED 2026-07-26. Previously "Local / Strata"
# only (8,648 of 517,836 segments, an ownership class, not a network).
# Now ALL NON-WINTER segments (516,868) via 00b_road_layers.R.
roads_all_sf <- road_layer("allseason") %>% clean_names()

pts_sf_allroad <- dat_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
allroad_vals <- nearest_road_km(pts_sf_allroad, roads_all_sf)

# ============================================================
# 3) WOLF HABITAT SUITABILITY (same as 08)
# ============================================================

wolf_r <- terra::rast(infile_wolf)
pts_sf_wolf <- dat_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(wolf_r))
wolf_vals <- terra::extract(wolf_r, vect(pts_sf_wolf))[, 2]

# ============================================================
# 4) MINE DISTANCE (same as 07)
# ============================================================

mines_r <- terra::rast(infile_mines)
pts_sf_mn <- dat_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(mines_r))
buf_mn     <- st_buffer(pts_sf_mn, dist = 500000)
mines_crop <- terra::crop(mines_r, vect(buf_mn))
mine_binary <- mines_crop
mine_binary[mine_binary == 0] <- NA
message("Computing distance to nearest mine raster (may take ~1 min)...")
dist_rast_mine <- terra::distance(mine_binary)
mine_vals <- terra::extract(dist_rast_mine, vect(pts_sf_mn))[, 2] / 1000

# ============================================================
# 5) DISTANCE TO NEAREST COMMUNITY (same as 06)
# Communities are derived from the population density raster (no
# independent settlements dataset was available for this study region —
# see 06_population_analyses.R for full rationale): contiguous
# populated cells are clustered into patches, small (< ~1 km^2)
# clusters are dropped as likely raster noise, and each remaining
# cluster's centroid is treated as a candidate community location.
# ============================================================

pop_r <- terra::rast(infile_pop)
pts_sf_pop <- dat_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(pop_r))

buf_pop_11  <- st_buffer(pts_sf_pop, dist = 500000)
pop_crop_11 <- terra::crop(pop_r, vect(buf_pop_11))

cell_area_km2_11     <- prod(res(pop_crop_11)) / 1e6
min_cluster_cells_11 <- max(1, round(1 / cell_area_km2_11))

pop_binary_11 <- pop_crop_11
pop_binary_11[pop_binary_11 <= 0] <- NA
pop_binary_11[!is.na(pop_binary_11)] <- 1

pop_patches_11 <- terra::patches(pop_binary_11, directions = 8)
patch_sizes_11 <- terra::freq(pop_patches_11)
small_patches_11 <- patch_sizes_11$value[patch_sizes_11$count < min_cluster_cells_11]
if (length(small_patches_11) > 0) pop_patches_11[pop_patches_11 %in% small_patches_11] <- NA

pop_polys_11 <- terra::as.polygons(pop_patches_11)
community_centers_11 <- terra::centroids(pop_polys_11) %>% st_as_sf()
message("Part 11: ", nrow(community_centers_11), " derived community locations.")

dist_mat_pop_11 <- st_distance(pts_sf_pop, community_centers_11)
pop_vals <- as.numeric(apply(dist_mat_pop_11, 1, min)) / 1000   # dist_community_km, kept as "pop_vals" for minimal downstream renaming

# ============================================================
# 5b) LAND COVER (PCA, same as 05/09/10)
# PCA is fit on the land-cover-complete sample (all dat_map rows with a
# full set of active land cover classes), matching how 05/09/10 build
# LC_PC1-4 — NOT restricted to the other five covariates' complete
# cases, so the PCA loadings are stable and comparable to the primary
# land-cover models. LC_PC1-4 scores are then carried forward and
# intersected with everything else in the final assembly step below.
# ============================================================

wolfstack <- terra::rast(infile_wolfstack)

pts_sf_lc <- dat_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(wolfstack))

lc_vals_raw <- terra::extract(wolfstack, vect(pts_sf_lc))

dat_lc_full <- dat_map %>%
  bind_cols(lc_vals_raw %>% select(-ID))

lc_class_cols_all <- c(
  "Water", "Wetlands", "OpenPeat", "TreedPeat",
  "Coniferous", "Sparse", "Deciduous", "NonTreed",
  "Anthropogenic", "Disturbance",
  "OpenBog", "TreedBog", "OpenFen", "TreedFen",
  "ThicketSwamp", "ConiferousSwamp", "DeciduousSwamp", "Marsh"
)
lc_bands_present <- intersect(lc_class_cols_all, names(wolfstack))

buf_lc <- st_buffer(pts_sf_lc, dist = 50000)
lc_stack_crop <- terra::crop(wolfstack[[lc_bands_present]], vect(buf_lc))
band_means   <- terra::global(lc_stack_crop, "mean", na.rm = TRUE)[, 1]
active_bands <- lc_bands_present[band_means > 0.005]
message("Land cover classes retained for PCA: ", paste(active_bands, collapse = ", "))

lc_complete <- dat_lc_full %>%
  filter(!apply(select(dat_lc_full, all_of(active_bands)), 1, anyNA))

pca_lc_full <- prcomp(lc_complete[, active_bands], scale. = TRUE)
cat("Land cover PCA variance explained (first 4 PCs):\n")
print(round(summary(pca_lc_full)$importance[2, 1:4] * 100, 1))

lc_scores <- lc_complete %>%
  select(pcid, date, tz_sample_number) %>%
  bind_cols(as_tibble(pca_lc_full$x[, 1:4])) %>%
  rename(LC_PC1 = PC1, LC_PC2 = PC2, LC_PC3 = PC3, LC_PC4 = PC4)

# ============================================================
# 6) PROJECTED COORDINATES FOR THE SPATIAL TERM
# Uses roads_named's CRS as the single reference projection for
# s(x_proj, y_proj) — the same projected, metres-based CRS already used
# for the spatial diagnostics in 03/04/10. This is independent of
# whichever CRS each raster happened to use for its own extraction
# above; it only matters that x_proj/y_proj are self-consistent across
# samples.
# ============================================================

pts_sf_xy <- dat_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(SPATIAL_SMOOTH_CRS)
xy_coords <- st_coordinates(pts_sf_xy)

# ============================================================
# 7) ASSEMBLE THE FULL COMBINED DATASET
# Road/wolf/mine/population attached by row order (dat_map row order is
# preserved throughout every extraction above); climate joined by
# (pcid, year, date) as noted above; land cover PCA scores joined by
# (pcid, date, tz_sample_number) — safe because lc_scores has at most
# one row per unique sample (no duplicate keys on the right-hand side).
# ============================================================

dat_full <- dat_map %>%
  mutate(
    dist_road_km      = road_vals,
    dist_allroad_km   = allroad_vals,
    wolf_suit         = wolf_vals,
    dist_mine_km      = mine_vals,
    dist_community_km = pop_vals,
    x_proj            = xy_coords[, 1],
    y_proj            = xy_coords[, 2]
  ) %>%
  left_join(clim_for_join, by = c("pcid", "year", "date")) %>%
  left_join(clim_for_join_fcm2, by = c("pcid", "year", "date")) %>%
  left_join(lc_scores, by = c("pcid", "date", "tz_sample_number"), relationship = "one-to-one") %>%
  filter(
    is.finite(snow_on_ground),
    is.finite(snow_on_ground_fcm2),
    is.finite(dist_road_km),
    is.finite(dist_allroad_km),
    is.finite(wolf_suit),
    is.finite(dist_mine_km),
    is.finite(dist_community_km),
    is.finite(LC_PC1), is.finite(LC_PC2), is.finite(LC_PC3), is.finite(LC_PC4),
    !is.na(preg_group)
  ) %>%
  mutate(
    survey               = factor(survey),
    year_f               = factor(year),
    log_dist_road_km     = log(dist_road_km),
    log_dist_allroad_km  = log(dist_allroad_km),
    log_wolf_suit        = log(wolf_suit + 0.001),
    log_dist_mine        = log(dist_mine_km + 0.1),
    log_dist_community   = log(dist_community_km + 0.1)
  )

cat("Full combined dataset — complete-case n (all 6 covariates present):",
    nrow(dat_full), "\n")

# ============================================================
# 8) CORTISOL — FULL MODEL, WITHOUT AND WITH SPATIAL TERM
# ============================================================

dat_full_cort <- dat_full %>% filter(is.finite(cortisol), !is.na(julian))
cat("Cortisol full-model complete-case n:", nrow(dat_full_cort), "\n")

m_cort_full <- gam(
  cortisol ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(snow_on_ground_fcm2,   k = 4, bs = "tp") +
    s(log_dist_road_km,      k = 5, bs = "tp") +
    s(log_dist_allroad_km,   k = 5, bs = "tp") +
    s(log_wolf_suit,         k = 6, bs = "tp") +
    s(log_dist_mine,         k = 6, bs = "tp") +
    s(log_dist_community,    k = 6, bs = "tp") +
    s(LC_PC1,                k = 3, bs = "tp") +
    s(LC_PC2,                k = 3, bs = "tp") +
    s(LC_PC3,                k = 3, bs = "tp") +
    s(LC_PC4,                k = 3, bs = "tp") +
    s(julian,                k = 8, bs = "tp"),
  data   = dat_full_cort,
  family = Gamma(link = "log"),
  method = "REML"
)

cat("\n=== Cortisol — FULL model, WITHOUT spatial term ===\n")
summary(m_cort_full)
capture.output(summary(m_cort_full), file = file.path(outdir, "GAM_summary_11_cortisol_full_nospatial.txt"))

png(file.path(outdir, "Fig_GAM_Cortisol_FullModel_NoSpatial.png"), 2600, 1800, res = 300)
plot(m_cort_full, pages = 1)
dev.off()

# ---- Spatial autocorrelation check on the full model's residuals ----

sf_cort_full <- dat_full_cort %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(SPATIAL_SMOOTH_CRS)

check_spatial_autocorr_full <- function(resid_vec, data_sf, label) {
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

  png(file.path(outdir, paste0("Variogram_", label, "_FullModel.png")), 1600, 1200, res = 200)
  print(plot(vgram, main = paste0(label, " — variogram of full-model deviance residuals"),
             xlab = "Distance (m)", ylab = "Semivariance"))
  dev.off()

  vgram_fit <- tryCatch(gstat::fit.variogram(vgram, gstat::vgm("Sph")), error = function(e) NULL)
  if (!is.null(vgram_fit)) {
    cat(sprintf("%s variogram model (nugget / partial sill / range):\n", label))
    print(vgram_fit)
    capture.output(print(vgram_fit),
                   file = file.path(outdir, paste0("Variogram_fit_", label, "_FullModel.txt")))
  } else {
    message(label, ": variogram model fit did not converge — inspect the empirical variogram plot directly.")
  }
  list(vgram = vgram, vgram_fit = vgram_fit)
}

if (!requireNamespace("gstat", quietly = TRUE)) install.packages("gstat", repos = "https://cloud.r-project.org")
library(gstat)

vres_cort_full <- check_spatial_autocorr_full(residuals(m_cort_full, type = "deviance"), sf_cort_full, "FCM")

m_cort_full_xy <- gam(
  cortisol ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(snow_on_ground_fcm2,   k = 4, bs = "tp") +
    s(log_dist_road_km,      k = 5, bs = "tp") +
    s(log_dist_allroad_km,   k = 5, bs = "tp") +
    s(log_wolf_suit,         k = 6, bs = "tp") +
    s(log_dist_mine,         k = 6, bs = "tp") +
    s(log_dist_community,    k = 6, bs = "tp") +
    s(LC_PC1,                k = 3, bs = "tp") +
    s(LC_PC2,                k = 3, bs = "tp") +
    s(LC_PC3,                k = 3, bs = "tp") +
    s(LC_PC4,                k = 3, bs = "tp") +
    s(julian,                k = 8, bs = "tp") +
    s(x_proj, y_proj,        bs = "tp", k = 30),
  data   = dat_full_cort,
  family = Gamma(link = "log"),
  method = "REML"
)

cat("\n=== Cortisol — FULL model, WITH spatial term s(x,y) ===\n")
summary(m_cort_full_xy)
capture.output(summary(m_cort_full_xy), file = file.path(outdir, "GAM_summary_11_cortisol_full_spatial.txt"))

png(file.path(outdir, "Fig_GAM_Cortisol_FullModel_Spatial.png"), 2600, 1800, res = 300)
plot(m_cort_full_xy, pages = 1)
dev.off()

vres_cort_full_after <- check_spatial_autocorr_full(residuals(m_cort_full_xy, type = "deviance"), sf_cort_full, "FCM_afterXY")

cat("\nCortisol full model — AIC/BIC, without vs. with spatial term:\n")
aic_cort_full <- AIC(m_cort_full, m_cort_full_xy)
bic_cort_full <- BIC(m_cort_full, m_cort_full_xy)
print(aic_cort_full); print(bic_cort_full)
capture.output(print(aic_cort_full), print(bic_cort_full),
               file = file.path(outdir, "GAM_AIC_11_cortisol_full_spatial_comparison.txt"))

cort_full_year_contrasts <- emmeans(m_cort_full_xy, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(cort_full_year_contrasts), file = file.path(outdir, "GAM_contrasts_11_cortisol_full_year.txt"))

# ============================================================
# 9) T3 — FULL MODEL, WITHOUT AND WITH SPATIAL TERM
# ============================================================

dat_full_t3 <- dat_full %>% filter(is.finite(t3), t3 > 0, !is.na(julian))
cat("T3 full-model complete-case n:", nrow(dat_full_t3), "\n")

m_t3_full <- gam(
  t3 ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(snow_on_ground,        k = 4, bs = "tp") +
    s(log_dist_road_km,      k = 5, bs = "tp") +
    s(log_dist_allroad_km,   k = 5, bs = "tp") +
    s(log_wolf_suit,         k = 6, bs = "tp") +
    s(log_dist_mine,         k = 6, bs = "tp") +
    s(log_dist_community,    k = 6, bs = "tp") +
    s(LC_PC1,                k = 3, bs = "tp") +
    s(LC_PC2,                k = 3, bs = "tp") +
    s(LC_PC3,                k = 3, bs = "tp") +
    s(LC_PC4,                k = 3, bs = "tp") +
    s(julian,                k = 8, bs = "tp"),
  data   = dat_full_t3,
  family = Gamma(link = "log"),
  method = "REML"
)

cat("\n=== T3 — FULL model, WITHOUT spatial term ===\n")
summary(m_t3_full)
capture.output(summary(m_t3_full), file = file.path(outdir, "GAM_summary_11_t3_full_nospatial.txt"))

png(file.path(outdir, "Fig_GAM_T3_FullModel_NoSpatial.png"), 2600, 1800, res = 300)
plot(m_t3_full, pages = 1)
dev.off()

sf_t3_full <- dat_full_t3 %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(SPATIAL_SMOOTH_CRS)

vres_t3_full <- check_spatial_autocorr_full(residuals(m_t3_full, type = "deviance"), sf_t3_full, "FT3")

m_t3_full_xy <- gam(
  t3 ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(snow_on_ground,        k = 4, bs = "tp") +
    s(log_dist_road_km,      k = 5, bs = "tp") +
    s(log_dist_allroad_km,   k = 5, bs = "tp") +
    s(log_wolf_suit,         k = 6, bs = "tp") +
    s(log_dist_mine,         k = 6, bs = "tp") +
    s(log_dist_community,    k = 6, bs = "tp") +
    s(LC_PC1,                k = 3, bs = "tp") +
    s(LC_PC2,                k = 3, bs = "tp") +
    s(LC_PC3,                k = 3, bs = "tp") +
    s(LC_PC4,                k = 3, bs = "tp") +
    s(julian,                k = 8, bs = "tp") +
    s(x_proj, y_proj,        bs = "tp", k = 30),
  data   = dat_full_t3,
  family = Gamma(link = "log"),
  method = "REML"
)

cat("\n=== T3 — FULL model, WITH spatial term s(x,y) ===\n")
summary(m_t3_full_xy)
capture.output(summary(m_t3_full_xy), file = file.path(outdir, "GAM_summary_11_t3_full_spatial.txt"))

png(file.path(outdir, "Fig_GAM_T3_FullModel_Spatial.png"), 2600, 1800, res = 300)
plot(m_t3_full_xy, pages = 1)
dev.off()

vres_t3_full_after <- check_spatial_autocorr_full(residuals(m_t3_full_xy, type = "deviance"), sf_t3_full, "FT3_afterXY")

cat("\nT3 full model — AIC/BIC, without vs. with spatial term:\n")
aic_t3_full <- AIC(m_t3_full, m_t3_full_xy)
bic_t3_full <- BIC(m_t3_full, m_t3_full_xy)
print(aic_t3_full); print(bic_t3_full)
capture.output(print(aic_t3_full), print(bic_t3_full),
               file = file.path(outdir, "GAM_AIC_11_t3_full_spatial_comparison.txt"))

t3_full_year_contrasts <- emmeans(m_t3_full_xy, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(t3_full_year_contrasts), file = file.path(outdir, "GAM_contrasts_11_t3_full_year.txt"))

# ============================================================
# 10) SIDE-BY-SIDE SIGNIFICANCE SUMMARY (without vs. with spatial term)
# ============================================================

cat("\n=== SUMMARY: does each covariate remain significant, without vs. with s(x,y)? ===\n")
cat("\n--- Cortisol (FCM) ---\n")
cat("Without spatial term:\n"); print(summary(m_cort_full)$s.table)
cat("With spatial term:\n");    print(summary(m_cort_full_xy)$s.table)
cat("\n--- T3 (FT3) ---\n")
cat("Without spatial term:\n"); print(summary(m_t3_full)$s.table)
cat("With spatial term:\n");    print(summary(m_t3_full_xy)$s.table)

# ============================================================
# 11) SAVE FITTED MODELS + DATA FOR REUSE
# 12_full_model_curve_comparison.R loads these .rds files rather than
# refitting the full model from scratch. This matters because wolf,
# mine, road, and land cover are all spatially correlated with one
# another — with covariates this collinear, mgcv's REML optimiser can
# land on slightly different (comparably valid) splits of wiggliness
# between correlated smooths depending on tiny floating-point/row-order
# differences between two independently-rebuilt datasets. Saving and
# reloading the actual fitted objects guarantees 12's comparison
# figure plots the exact same model shown here, not a re-fit that
# merely shares the same formula.
# ============================================================

saveRDS(m_cort_full,    file.path(outdir, "model_11_m_cort_full.rds"))
saveRDS(m_cort_full_xy, file.path(outdir, "model_11_m_cort_full_xy.rds"))
saveRDS(m_t3_full,      file.path(outdir, "model_11_m_t3_full.rds"))
saveRDS(m_t3_full_xy,   file.path(outdir, "model_11_m_t3_full_xy.rds"))
saveRDS(dat_full_cort,  file.path(outdir, "model_11_dat_full_cort.rds"))
saveRDS(dat_full_t3,    file.path(outdir, "model_11_dat_full_t3.rds"))

message("\nPart 11 (full combined model, with/without spatial term) complete — outputs in: ", outdir)
message("Fitted models saved to ", outdir, " as model_11_*.rds for reuse by 12_full_model_curve_comparison.R")
