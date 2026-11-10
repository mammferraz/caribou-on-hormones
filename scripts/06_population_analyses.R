# ============================================================
# 06_population_analyses.R
#
# PART 6 — HUMAN SETTLEMENT PROXIMITY (DISTANCE TO NEAREST COMMUNITY)
# Purpose:
#   • Derive candidate "community" locations from the population
#     density raster (no independent settlements dataset was available
#     for this study region — see note below), then compute distance
#     from each caribou sample to the nearest one
#   • Test effects on cortisol and T3 by year and survey
#
# CHANGE FROM EARLIER VERSION: this script previously modelled raw
# population density directly (log1p(pop_density)). That covariate was
# almost entirely zero (median = 0, 3rd quartile = 0 people/km² across
# samples) since these ranges are remote wilderness, which is a poor
# basis for a smooth term. Distance to the nearest human settlement is
# a more standard, more interpretable exposure metric, and is
# structurally identical to the winter-road and mine-distance
# covariates already used elsewhere in this project (03/07).
#
# DERIVING "COMMUNITIES" FROM THE POPULATION RASTER: an actual
# settlements point/polygon layer for northern Ontario was not
# available in this project (the only settlements file on hand,
# Greenbelt_towns_and_villages.geojson, covers the southern Ontario
# Greenbelt — hundreds of km from these ranges — and is not used here).
# Instead, contiguous cells with detectable population density are
# clustered into patches (same terra::patches() technique already used
# for mine clusters in 07_mines_analyses.R), small clusters likely
# representing raster/interpolation noise rather than a real settlement
# are dropped, and each remaining cluster's centroid is treated as a
# candidate community location. The cluster-size distribution and a
# diagnostic map are saved so this choice can be visually sanity-checked
# and the threshold/minimum-size parameters adjusted if needed.
# ============================================================

source("00_data_processing.R")

outdir <- "caribou_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# I/O
# ============================================================

infile_pop <- "population_density.tif"

# ============================================================
# 1) DERIVE COMMUNITY LOCATIONS FROM THE POPULATION RASTER
# ============================================================

pop_r <- terra::rast(infile_pop)

cat("Population density CRS:\n");   print(crs(pop_r))
cat("Population density resolution:\n"); print(res(pop_r))

pts_sf <- dat_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(pop_r))

buf_pop   <- st_buffer(pts_sf, dist = 500000)   # 500 km buffer, matching 07's mine-distance convention
pop_crop  <- terra::crop(pop_r, vect(buf_pop))

# Threshold: any cell with detectable population density (> 0) is
# flagged as "settled". Contiguous settled cells are clustered into
# patches; patches smaller than ~1 km^2 are dropped as likely noise
# before taking cluster centroids as candidate community locations.
# These two parameters are the main tunable choices in this derivation
# — adjust and re-inspect Fig_DerivedCommunities_Diagnostic.png and the
# printed cluster-size summary below if the result looks wrong (e.g.
# far too many tiny clusters, or everything merging into one blob).
pop_threshold     <- 0
min_cluster_km2   <- 1
cell_area_km2     <- prod(res(pop_crop)) / 1e6
min_cluster_cells <- max(1, round(min_cluster_km2 / cell_area_km2))
cat(sprintf("Cell area: %.4f km^2 | minimum cluster size: %d cells (~%.1f km^2)\n",
            cell_area_km2, min_cluster_cells, min_cluster_cells * cell_area_km2))

pop_binary <- pop_crop
pop_binary[pop_binary <= pop_threshold] <- NA
pop_binary[!is.na(pop_binary)] <- 1

pop_patches <- terra::patches(pop_binary, directions = 8)

patch_sizes <- terra::freq(pop_patches)
cat("Population-cluster size distribution (cells), before size filtering:\n")
print(summary(patch_sizes$count))

small_patches <- patch_sizes$value[patch_sizes$count < min_cluster_cells]
if (length(small_patches) > 0) {
  cat(sprintf("Dropping %d of %d clusters smaller than %d cells (likely noise)\n",
              length(small_patches), nrow(patch_sizes), min_cluster_cells))
  pop_patches[pop_patches %in% small_patches] <- NA
}

pop_polys <- terra::as.polygons(pop_patches)
community_centers <- terra::centroids(pop_polys) %>% st_as_sf()
cat("Number of derived community locations after filtering small clusters:", nrow(community_centers), "\n")

# Diagnostic map — inspect this before trusting the derived distances.
png(file.path(outdir, "Fig_DerivedCommunities_Diagnostic.png"), 2400, 1800, res = 200)
plot(pop_crop, main = "Population density with derived 'community' centroids (red x)")
plot(vect(community_centers), add = TRUE, col = "red", pch = 4, cex = 1.2, lwd = 2)
dev.off()

# ============================================================
# 2) DISTANCE FROM EACH SAMPLE TO NEAREST DERIVED COMMUNITY
# ============================================================

dist_mat_pop      <- st_distance(pts_sf, community_centers)
dist_community_km <- as.numeric(apply(dist_mat_pop, 1, min)) / 1000

dat_pop <- dat_map %>%
  mutate(dist_community_km = dist_community_km) %>%
  filter(is.finite(dist_community_km), !is.na(preg_group)) %>%
  mutate(survey = factor(survey), year_f = factor(year),
         # +0.1 km offset guards against log(0) for any sample that falls
         # inside a derived community cluster; matches the offset
         # convention used for mine distance (07).
         log_dist_community = log(dist_community_km + 0.1))

cat("Distance to nearest derived community (km):\n")
print(summary(dat_pop$dist_community_km))

# ============================================================
# 3) MAP: DISTANCE TO NEAREST COMMUNITY + CARIBOU LOCATIONS
# ============================================================

pts_sf_proj <- dat_pop %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(pop_r))

# Map is cropped tightly to the sampling area (50 km buffer around the
# caribou points, same convention as the mine-distance map in 07), NOT
# the full 500 km buffer used for the distance calculation itself.
# Plotting all 557 derived community centroids from that much larger
# buffer previously blew out the map's automatic extent to cover most
# of the raster's bounding box instead of just the study area.
buf           <- st_buffer(pts_sf_proj, dist = 50000)
pop_crop_map  <- terra::crop(pop_crop, vect(buf))
dist_rast_pop <- terra::distance(pop_crop_map, vect(community_centers))
dist_crop     <- dist_rast_pop / 1000  # km
pop_df        <- as.data.frame(dist_crop, xy = TRUE, na.rm = FALSE)
names(pop_df) <- c("x", "y", "dist_community_km")

# Only show community centroids that actually fall within the mapped
# area — matches how 07_mines_analyses.R derives its mine-cluster
# centroids from the cropped (not full-buffer) raster.
community_centers_map <- community_centers[st_intersects(community_centers, buf, sparse = FALSE), ]

buf_bbox <- st_bbox(buf)

p_pop_map <- ggplot() +
  geom_raster(data = pop_df, aes(x, y, fill = dist_community_km)) +
  geom_sf(data = pts_sf_proj,
          aes(colour = survey, shape = preg_group),
          size = 1.8, alpha = 0.95) +
  geom_sf(data = community_centers_map, shape = 4, size = 2, colour = "red", stroke = 1) +
  scale_fill_viridis_c(name = "Distance to\nnearest community (km)", option = "magma",
                        direction = -1, na.value = "white") +
  scale_colour_manual(values = survey_colours, na.value = "grey50") +
  scale_shape_manual(values = preg_shapes) +
  coord_sf(xlim = c(buf_bbox["xmin"], buf_bbox["xmax"]),
           ylim = c(buf_bbox["ymin"], buf_bbox["ymax"]), expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.35) +
  labs(title  = "Distance to nearest community and caribou sampling locations",
       x = "Easting", y = "Northing",
       colour = "Survey/Range", shape = "Group") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"), legend.position = "right")

ggsave(file.path(outdir, "Fig_CommunityDistance_Map.png"),
       p_pop_map, width = 8, height = 6, dpi = 300)


# ============================================================
# 4) DISTANCE TO NEAREST COMMUNITY BY YEAR × SURVEY
# ============================================================

p_popden_year <- ggplot(
  dat_pop,
  aes(x = factor(year), y = dist_community_km, colour = survey)
) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
              alpha = 0.5, size = 1.8) +
  geom_boxplot(aes(group = interaction(year, survey)),
               position = position_dodge(0.6),
               outlier.shape = NA, fill = NA, linewidth = 0.7, width = 0.4) +
  scale_colour_manual(values = survey_colours, na.value = "grey50") +
  scale_y_log10() +
  labs(x = "Year", y = "Distance to nearest community (km, log scale)",
       colour = "Survey/Range",
       title  = "Distance to nearest community at caribou locations by year") +
  theme_pub

ggsave(file.path(outdir, "Fig_CommunityDistance_by_Year_Survey.png"),
       p_popden_year, width = 8, height = 5, dpi = 300)


# ============================================================
# 5) GAMs — DISTANCE-TO-COMMUNITY EFFECTS
# ============================================================

# --- 5A) Cortisol ---

dat_pop_cort <- dat_pop %>%
  filter(is.finite(cortisol), !is.na(julian))

m_cort_pop <- gam(
  cortisol ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(log_dist_community, k = 6, bs = "tp") +
    s(julian,              k = 8, bs = "tp"),
  data   = dat_pop_cort,
  family = Gamma(link = "log"),
  method = "REML"
)

summary(m_cort_pop)
capture.output(summary(m_cort_pop), file = file.path(outdir, "GAM_summary_06_cortisol_pop.txt"))
gam.check(m_cort_pop)

p_cort_pop <- draw(m_cort_pop, select = "s(log_dist_community)") +
  ggtitle("Cortisol — partial effect of distance to nearest community")

ggsave(file.path(outdir, "Fig_Cortisol_CommunityDistance_PartialEffect.png"),
       p_cort_pop, width = 6, height = 5, dpi = 300)

png(file.path(outdir, "Fig_GAM_Cortisol_CommunityDistance.png"), 2000, 1400, res = 300)
plot(m_cort_pop, pages = 1)
dev.off()

cort_pop_year_contrasts <- emmeans(m_cort_pop, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(cort_pop_year_contrasts), file = file.path(outdir, "GAM_contrasts_06_cortisol_pop_year.txt"))

# --- 5B) T3 ---

dat_pop_t3 <- dat_pop %>%
  filter(is.finite(t3), t3 > 0, !is.na(julian))

m_t3_pop <- gam(
  t3 ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(log_dist_community, k = 6, bs = "tp") +
    s(julian,              k = 8, bs = "tp"),
  data   = dat_pop_t3,
  family = Gamma(link = "log"),
  method = "REML"
)

summary(m_t3_pop)
capture.output(summary(m_t3_pop), file = file.path(outdir, "GAM_summary_06_t3_pop.txt"))

p_t3_pop <- draw(m_t3_pop, select = "s(log_dist_community)") +
  ggtitle("T3 — partial effect of distance to nearest community")

ggsave(file.path(outdir, "Fig_T3_CommunityDistance_PartialEffect.png"),
       p_t3_pop, width = 6, height = 5, dpi = 300)

png(file.path(outdir, "Fig_GAM_T3_CommunityDistance.png"), 2000, 1400, res = 300)
plot(m_t3_pop, pages = 1)
dev.off()

t3_pop_year_contrasts <- emmeans(m_t3_pop, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(t3_pop_year_contrasts), file = file.path(outdir, "GAM_contrasts_06_t3_pop_year.txt"))

# ============================================================
# 6) SPATIAL AUTOCORRELATION DIAGNOSTICS (variogram on residuals)
# Purpose: same diagnostic applied to the road-distance models in
#   03_winter_road_analyses.R, now extended to distance-to-community —
#   check whether residual spatial autocorrelation could be inflating
#   the apparent significance of the log_dist_community smooth.
#
# Some samples share IDENTICAL GPS coordinates (repeat visits to the
# same collection site). A raw per-sample variogram would treat those
# as separate zero-distance observations and read as pure nugget noise
# rather than genuine small-scale variation, so residuals are averaged
# to one value per unique sampling location before the variogram is
# computed. If the variogram shows real spatial structure (sill clearly
# above the nugget, with a non-trivial range), the community-distance
# models are refit with an added s(x, y) smooth on the projected
# coordinates, and the log_dist_community effect is re-checked with
# that term included.
# ============================================================

if (!requireNamespace("gstat", quietly = TRUE)) install.packages("gstat", repos = "https://cloud.r-project.org")
library(gstat)

sf_cort_pop <- dat_pop_cort %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(pop_r))
sf_t3_pop <- dat_pop_t3 %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(pop_r))

check_spatial_autocorr_pop <- function(resid_vec, data_sf, label) {
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

  png(file.path(outdir, paste0("Variogram_", label, "_CommunityModel.png")), 1600, 1200, res = 200)
  print(plot(vgram, main = paste0(label, " — variogram of community-distance-model deviance residuals"),
             xlab = "Distance (m)", ylab = "Semivariance"))
  dev.off()

  vgram_fit <- tryCatch(gstat::fit.variogram(vgram, gstat::vgm("Sph")), error = function(e) NULL)
  if (!is.null(vgram_fit)) {
    cat(sprintf("%s variogram model (nugget / partial sill / range):\n", label))
    print(vgram_fit)
    capture.output(print(vgram_fit),
                   file = file.path(outdir, paste0("Variogram_fit_", label, "_CommunityModel.txt")))
  } else {
    message(label, ": variogram model fit did not converge — inspect the empirical variogram plot directly.")
  }
  list(vgram = vgram, vgram_fit = vgram_fit)
}

vres_cort_pop <- check_spatial_autocorr_pop(residuals(m_cort_pop, type = "deviance"), sf_cort_pop, "FCM")
vres_t3_pop   <- check_spatial_autocorr_pop(residuals(m_t3_pop,   type = "deviance"), sf_t3_pop,   "FT3")

# ---- Refit with s(x, y) added, and re-check whether log_dist_community survives ----

dat_cort_pop_xy <- dat_pop_cort %>%
  bind_cols(st_coordinates(sf_cort_pop) %>% as.data.frame() %>% rename(x_proj = X, y_proj = Y))
dat_t3_pop_xy <- dat_pop_t3 %>%
  bind_cols(st_coordinates(sf_t3_pop) %>% as.data.frame() %>% rename(x_proj = X, y_proj = Y))

m_cort_pop_xy <- gam(
  cortisol ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(log_dist_community, k = 6, bs = "tp") +
    s(julian,              k = 8, bs = "tp") +
    s(x_proj, y_proj, bs = "tp", k = 30),
  data   = dat_cort_pop_xy,
  family = Gamma(link = "log"),
  method = "REML"
)

plot(m_cort_pop_xy, pages = 1, scheme = 2)

m_t3_pop_xy <- gam(
  t3 ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(log_dist_community, k = 6, bs = "tp") +
    s(julian,              k = 8, bs = "tp") +
    s(x_proj, y_proj, bs = "tp", k = 30),
  data   = dat_t3_pop_xy,
  family = Gamma(link = "log"),
  method = "REML"
)

plot(m_t3_pop_xy, pages = 1, scheme = 2)

cat("\nCommunity-distance models with s(x, y) added — AIC/BIC comparison against the originally reported (no-spatial-term) models:\n")
aic_cort_pop_xy <- AIC(m_cort_pop, m_cort_pop_xy)
aic_t3_pop_xy   <- AIC(m_t3_pop,   m_t3_pop_xy)
bic_cort_pop_xy <- BIC(m_cort_pop, m_cort_pop_xy)
bic_t3_pop_xy   <- BIC(m_t3_pop,   m_t3_pop_xy)
print(aic_cort_pop_xy); print(bic_cort_pop_xy)
print(aic_t3_pop_xy); print(bic_t3_pop_xy)
capture.output(print(aic_cort_pop_xy), print(bic_cort_pop_xy), print(aic_t3_pop_xy), print(bic_t3_pop_xy),
               file = file.path(outdir, "GAM_AIC_06_spatial_xy_comparison.txt"))

cat("\nDoes log_dist_community remain significant once s(x, y) is included?\n")
print(summary(m_cort_pop_xy)$s.table)
print(summary(m_t3_pop_xy)$s.table)
capture.output(summary(m_cort_pop_xy), file = file.path(outdir, "GAM_summary_06_cortisol_pop_xy.txt"))
capture.output(summary(m_t3_pop_xy),   file = file.path(outdir, "GAM_summary_06_t3_pop_xy.txt"))

# Re-check the variogram on the xy-augmented models' residuals — if
# s(x, y) has done its job, this should look flatter than the original.
vres_cort_pop_after <- check_spatial_autocorr_pop(residuals(m_cort_pop_xy, type = "deviance"), sf_cort_pop, "FCM_afterXY")
vres_t3_pop_after   <- check_spatial_autocorr_pop(residuals(m_t3_pop_xy,   type = "deviance"), sf_t3_pop,   "FT3_afterXY")

message("Part 6 spatial autocorrelation diagnostics complete — outputs in: ", outdir)

# ============================================================
# 7) PRIMARY MODEL DECISION — spatially-adjusted model reported
# Project decision: for landscape covariates (road, mine, wolf, land
# cover, community distance), the s(x, y)-augmented model is reported
# as the PRIMARY model, regardless of BIC (see 03_winter_road_analyses.R
# Part 2E for full rationale). m_cort_pop_xy / m_t3_pop_xy are
# therefore the primary reported community-distance models. The outputs
# below OVERWRITE the original (non-spatial) GAM_summary/Fig/contrasts
# files from Part 5 so every file on disk reflects the model actually
# reported in the manuscript.
# ============================================================

summary(m_cort_pop_xy)
capture.output(summary(m_cort_pop_xy), file = file.path(outdir, "GAM_summary_06_cortisol_pop.txt"))

p_cort_pop <- draw(m_cort_pop_xy, select = "s(log_dist_community)") +
  ggtitle("Cortisol — partial effect of distance to nearest community")
ggsave(file.path(outdir, "Fig_Cortisol_CommunityDistance_PartialEffect.png"),
       p_cort_pop, width = 6, height = 5, dpi = 300)

png(file.path(outdir, "Fig_GAM_Cortisol_CommunityDistance.png"), 2000, 1400, res = 300)
plot(m_cort_pop_xy, pages = 1)
dev.off()

cort_pop_year_contrasts <- emmeans(m_cort_pop_xy, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(cort_pop_year_contrasts), file = file.path(outdir, "GAM_contrasts_06_cortisol_pop_year.txt"))

summary(m_t3_pop_xy)
capture.output(summary(m_t3_pop_xy), file = file.path(outdir, "GAM_summary_06_t3_pop.txt"))

p_t3_pop <- draw(m_t3_pop_xy, select = "s(log_dist_community)") +
  ggtitle("T3 — partial effect of distance to nearest community")
ggsave(file.path(outdir, "Fig_T3_CommunityDistance_PartialEffect.png"),
       p_t3_pop, width = 6, height = 5, dpi = 300)

png(file.path(outdir, "Fig_GAM_T3_CommunityDistance.png"), 2000, 1400, res = 300)
plot(m_t3_pop_xy, pages = 1)
dev.off()

t3_pop_year_contrasts <- emmeans(m_t3_pop_xy, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(t3_pop_year_contrasts), file = file.path(outdir, "GAM_contrasts_06_t3_pop_year.txt"))

message("Part 7 (primary model = spatially-adjusted) complete — GAM_summary_06_cortisol_pop.txt / GAM_summary_06_t3_pop.txt now reflect the s(x,y)-augmented models.")

message("Part 6 complete — outputs in: ", outdir)
