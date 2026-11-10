# ============================================================
# 12_full_model_curve_comparison.R
#
# PART 12 — COMPARING PARTIAL-EFFECT CURVE SHAPES ACROSS MODELS
# Purpose:
#   • Rather than compare p-values, directly compare the SHAPE of each
#     covariate's partial-effect curve across three model
#     specifications:
#       (a) "Individual"     — that covariate's own primary spatially-
#                               adjusted model (matching 03/05/06/07/08),
#                               fit on its own complete-case sample.
#       (b) "Full, no space" — the seven-covariate combined model from
#                               11_full_model_analyses.R, WITHOUT s(x,y).
#       (c) "Full + space"   — the same combined model WITH s(x,y).
#   • This produces a supplementary figure so we can visually judge
#     whether a covariate's estimated relationship changes shape (not
#     just significance) once it competes with the other covariates,
#     and once spatial location is accounted for.
#   • Snow-on-ground does not have its own "individual" model elsewhere
#     in this project (02's climate model includes minimum temperature
#     too, with a cyclic Julian-day spline, a different spec). For
#     comparability with the other covariates here, a fresh
#     snow-only + s(x,y) model is fit as its "individual" reference —
#     this is NOT the same as 02's primary reported climate model, and
#     should be labelled as such in the figure/legend.
#   • Road is now split into two panels: winter roads (as reported
#     elsewhere in the manuscript) and all-year roads (Local/Strata),
#     matching the two road terms fit in 11_full_model_analyses.R.
#   • FCM's snow term uses the 2-day short-lag window; FT3's uses the
#     5-day window — see 02_climate_analyses.R header for rationale.
# ============================================================

source("00_data_processing.R")
source("00b_road_layers.R")

library(gratia)
library(patchwork)

outdir    <- "caribou_outputs"
outdir_ms <- file.path(outdir, "manuscript")
dir.create(outdir_ms, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("gstat", quietly = TRUE)) install.packages("gstat", repos = "https://cloud.r-project.org")

# ============================================================
# I/O
# ============================================================

infile_climate   <- "climatedaily_all_surveys.csv"
infile_roads     <- "MNR_Road_Segments.geojson"
infile_roads_all <- "MNR_Road_Segments2.geojson"   # all classes; see 00b_road_layers.R
infile_lookup    <- "winter_roads_all2.csv"
infile_wolf      <- "wolfmean_t.tif"
infile_mines     <- "mines.tif"
infile_pop       <- "population_density.tif"
infile_wolfstack <- "wolfstack.tif"

# ============================================================
# 1) EXTRACT ALL COVARIATES (identical to 11_full_model_analyses.R —
# see that script for the full rationale on each join/extraction step)
# ============================================================

# --- Climate ---
# FCM (cortisol) uses a 2-day short-lag window (snow_on_ground_fcm2);
# FT3 keeps the 5-day window (snow_on_ground) — see
# 02_climate_analyses.R header for rationale.
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

clim_for_join <- climate_lags %>% select(pcid, year, date, snow_on_ground)
clim_for_join_fcm2 <- climate_lags_fcm2 %>%
  select(pcid, year, date, snow_on_ground) %>%
  rename(snow_on_ground_fcm2 = snow_on_ground)

# --- Winter road distance ---
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

# --- All-season road distance (all non-winter segments) ---
# NOTE: file is ~1.3 GB — loading takes ~30-60 s. Matches
# 11_full_model_analyses.R's "2b) ALL-SEASON ROAD DISTANCE" block exactly.
# ALL-SEASON ROADS — CORRECTED 2026-07-26. Previously "Local / Strata"
# only (8,648 of 517,836 segments, an ownership class, not a network).
# Now ALL NON-WINTER segments (516,868) via 00b_road_layers.R.
roads_all_sf <- road_layer("allseason") %>% clean_names()

pts_sf_allroad <- dat_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
allroad_vals <- nearest_road_km(pts_sf_allroad, roads_all_sf)

# --- Wolf habitat suitability ---
wolf_r <- terra::rast(infile_wolf)
pts_sf_wolf <- dat_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(wolf_r))
wolf_vals <- terra::extract(wolf_r, vect(pts_sf_wolf))[, 2]

# --- Mine distance ---
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

# --- Distance to nearest community (derived from population raster —
# no independent settlements dataset was available for this region;
# see 06_population_analyses.R for full rationale) ---
pop_r <- terra::rast(infile_pop)
pts_sf_pop <- dat_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(pop_r))

buf_pop_12  <- st_buffer(pts_sf_pop, dist = 500000)
pop_crop_12 <- terra::crop(pop_r, vect(buf_pop_12))

cell_area_km2_12     <- prod(res(pop_crop_12)) / 1e6
min_cluster_cells_12 <- max(1, round(1 / cell_area_km2_12))

pop_binary_12 <- pop_crop_12
pop_binary_12[pop_binary_12 <= 0] <- NA
pop_binary_12[!is.na(pop_binary_12)] <- 1

pop_patches_12 <- terra::patches(pop_binary_12, directions = 8)
patch_sizes_12 <- terra::freq(pop_patches_12)
small_patches_12 <- patch_sizes_12$value[patch_sizes_12$count < min_cluster_cells_12]
if (length(small_patches_12) > 0) pop_patches_12[pop_patches_12 %in% small_patches_12] <- NA

pop_polys_12 <- terra::as.polygons(pop_patches_12)
community_centers_12 <- terra::centroids(pop_polys_12) %>% st_as_sf()

dist_mat_pop_12 <- st_distance(pts_sf_pop, community_centers_12)
pop_vals <- as.numeric(apply(dist_mat_pop_12, 1, min)) / 1000   # dist_community_km

# --- Land cover (PCA) ---
wolfstack <- terra::rast(infile_wolfstack)
pts_sf_lc <- dat_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(wolfstack))
lc_vals_raw <- terra::extract(wolfstack, vect(pts_sf_lc))
dat_lc_full <- dat_map %>% bind_cols(lc_vals_raw %>% select(-ID))

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

lc_complete <- dat_lc_full %>%
  filter(!apply(select(dat_lc_full, all_of(active_bands)), 1, anyNA))
pca_lc_full <- prcomp(lc_complete[, active_bands], scale. = TRUE)
lc_scores <- lc_complete %>%
  select(pcid, date, tz_sample_number) %>%
  bind_cols(as_tibble(pca_lc_full$x[, 1:4])) %>%
  rename(LC_PC1 = PC1, LC_PC2 = PC2, LC_PC3 = PC3, LC_PC4 = PC4)

# --- Projected coordinates (single reference CRS for every s(x,y) term) ---
pts_sf_xy <- dat_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(SPATIAL_SMOOTH_CRS)
xy_coords <- st_coordinates(pts_sf_xy)

dat_all <- dat_map %>%
  mutate(
    dist_road_km      = road_vals,
    dist_allroad_km   = allroad_vals,
    wolf_suit         = wolf_vals,
    dist_mine_km      = mine_vals,
    dist_community_km = pop_vals,
    x_proj            = xy_coords[, 1],
    y_proj            = xy_coords[, 2],
    survey            = factor(survey),
    year_f            = factor(year)
  ) %>%
  left_join(clim_for_join, by = c("pcid", "year", "date")) %>%
  left_join(clim_for_join_fcm2, by = c("pcid", "year", "date")) %>%
  left_join(lc_scores, by = c("pcid", "date", "tz_sample_number"), relationship = "one-to-one") %>%
  mutate(
    log_dist_road_km    = log(dist_road_km),
    log_dist_allroad_km = log(dist_allroad_km),
    log_wolf_suit       = log(wolf_suit + 0.001),
    log_dist_mine       = log(dist_mine_km + 0.1),
    log_dist_community  = log(dist_community_km + 0.1)
  )

# ============================================================
# 2) INDIVIDUAL MODELS — each covariate's own primary spatially-
# adjusted model, fit on its own complete-case sample (matching
# 03/05/06/07/08's k values exactly; snow is a fresh snow-only + s(x,y)
# comparator, see header note above)
# ============================================================

fit_individual <- function(response, covariate_term, covariate_k, extra_filter,
                            julian_k = 8) {
  d <- dat_all %>%
    filter(is.finite(.data[[response]]), !is.na(julian), !is.na(preg_group))
  if (response == "t3") d <- d %>% filter(t3 > 0)
  d <- extra_filter(d)
  form <- as.formula(sprintf(
    "%s ~ s(survey, bs = 're') + year_f + preg_group + s(%s, k = %d, bs = 'tp') + s(julian, k = %d, bs = 'tp') + s(x_proj, y_proj, bs = 'tp', k = 30)",
    response, covariate_term, covariate_k, julian_k
  ))
  gam(form, data = d, family = Gamma(link = "log"), method = "REML")
}

# NOTE: snow's individual model uses hormone-specific short-lag windows —
# FCM: 2-day (snow_on_ground_fcm2); FT3: 5-day (snow_on_ground) — see
# 02_climate_analyses.R header for rationale.
filt_snow_fcm2 <- function(d) d %>% filter(is.finite(snow_on_ground_fcm2))
filt_snow      <- function(d) d %>% filter(is.finite(snow_on_ground))
filt_road      <- function(d) d %>% filter(is.finite(dist_road_km))
filt_allroad   <- function(d) d %>% filter(is.finite(dist_allroad_km))
filt_wolf      <- function(d) d %>% filter(is.finite(wolf_suit))
filt_mine      <- function(d) d %>% filter(is.finite(dist_mine_km))
filt_pop       <- function(d) d %>% filter(is.finite(dist_community_km))
filt_lc        <- function(d) d %>% filter(is.finite(LC_PC1), is.finite(LC_PC2), is.finite(LC_PC3), is.finite(LC_PC4))

m_cort_snow_ind <- fit_individual("cortisol", "snow_on_ground_fcm2", 4, filt_snow_fcm2, julian_k = 8)
m_t3_snow_ind   <- fit_individual("t3",       "snow_on_ground",      4, filt_snow,      julian_k = 8)
m_cort_road_ind <- fit_individual("cortisol", "log_dist_road_km", 5, filt_road, julian_k = 5)
m_t3_road_ind   <- fit_individual("t3",       "log_dist_road_km", 5, filt_road, julian_k = 5)
m_cort_allroad_ind <- fit_individual("cortisol", "log_dist_allroad_km", 5, filt_allroad, julian_k = 5)
m_t3_allroad_ind   <- fit_individual("t3",       "log_dist_allroad_km", 5, filt_allroad, julian_k = 5)
m_cort_wolf_ind <- fit_individual("cortisol", "log_wolf_suit",    6, filt_wolf, julian_k = 8)
m_t3_wolf_ind   <- fit_individual("t3",       "log_wolf_suit",    6, filt_wolf, julian_k = 8)
m_cort_mine_ind <- fit_individual("cortisol", "log_dist_mine",    6, filt_mine, julian_k = 8)
m_t3_mine_ind   <- fit_individual("t3",       "log_dist_mine",    6, filt_mine, julian_k = 8)

m_cort_pop_ind <- fit_individual("cortisol", "log_dist_community", 6, filt_pop, julian_k = 8)
m_t3_pop_ind   <- fit_individual("t3",       "log_dist_community", 6, filt_pop, julian_k = 8)

# Land cover's individual model uses 4 PCs, matching 05, fit directly.
dat_lc_cort_ind <- dat_all %>% filter(is.finite(cortisol), !is.na(julian), !is.na(preg_group)) %>% filt_lc()
dat_lc_t3_ind   <- dat_all %>% filter(is.finite(t3), t3 > 0, !is.na(julian), !is.na(preg_group)) %>% filt_lc()
m_cort_lc_ind <- gam(cortisol ~ s(survey, bs = "re") + year_f + preg_group +
                       s(LC_PC1, k = 3, bs = "tp") + s(LC_PC2, k = 3, bs = "tp") +
                       s(LC_PC3, k = 3, bs = "tp") + s(LC_PC4, k = 3, bs = "tp") +
                       s(julian, k = 8, bs = "tp") + s(x_proj, y_proj, bs = "tp", k = 30),
                     data = dat_lc_cort_ind, family = Gamma(link = "log"), method = "REML")
m_t3_lc_ind   <- gam(t3 ~ s(survey, bs = "re") + year_f + preg_group +
                       s(LC_PC1, k = 3, bs = "tp") + s(LC_PC2, k = 3, bs = "tp") +
                       s(LC_PC3, k = 3, bs = "tp") + s(LC_PC4, k = 3, bs = "tp") +
                       s(julian, k = 8, bs = "tp") + s(x_proj, y_proj, bs = "tp", k = 30),
                     data = dat_lc_t3_ind, family = Gamma(link = "log"), method = "REML")

message("Individual models fit — n's: snow(FCM,2d)=", nrow(filt_snow_fcm2(dat_all %>% filter(is.finite(cortisol)))),
        ", snow(FT3,5d)=", nrow(filt_snow(dat_all %>% filter(is.finite(t3), t3 > 0))),
        ", winter_road=", nrow(filt_road(dat_all %>% filter(is.finite(cortisol)))),
        ", all_road=", nrow(filt_allroad(dat_all %>% filter(is.finite(cortisol)))),
        ", wolf=", nrow(filt_wolf(dat_all %>% filter(is.finite(cortisol)))),
        ", mine=", nrow(filt_mine(dat_all %>% filter(is.finite(cortisol)))),
        ", pop=",  nrow(filt_pop(dat_all %>% filter(is.finite(cortisol)))),
        ", lc=",   nrow(dat_lc_cort_ind))

# ============================================================
# 3) FULL COMBINED MODEL — without and with s(x,y)
# (identical specification to 11_full_model_analyses.R)
#
# IMPORTANT: this section REUSES the exact fitted model objects saved
# by 11_full_model_analyses.R (as caribou_outputs/model_11_*.rds)
# rather than refitting them from a freshly rebuilt dataset. Wolf,
# mine, road, and land cover are all spatially correlated with one
# another; with covariates this collinear, mgcv's REML optimiser can
# land on slightly different (comparably valid) splits of wiggliness
# between correlated smooths depending on tiny floating-point/row-order
# differences between two independently-rebuilt versions of "the same"
# dataset. Loading 11's actual saved objects guarantees this script's
# comparison figure shows the SAME full-model curves you already
# inspected in 11's own diagnostic plots, not a re-fit that merely
# shares the same formula. If 11 has not been run yet (no .rds files
# found), this falls back to refitting here, with a warning — run 11
# first if you want a guaranteed match.
# ============================================================

rds_paths <- list(
  m_cort_full    = file.path(outdir, "model_11_m_cort_full.rds"),
  m_cort_full_xy = file.path(outdir, "model_11_m_cort_full_xy.rds"),
  m_t3_full      = file.path(outdir, "model_11_m_t3_full.rds"),
  m_t3_full_xy   = file.path(outdir, "model_11_m_t3_full_xy.rds"),
  dat_full_cort  = file.path(outdir, "model_11_dat_full_cort.rds"),
  dat_full_t3    = file.path(outdir, "model_11_dat_full_t3.rds")
)

if (all(file.exists(unlist(rds_paths)))) {

  message("Loading full-model objects saved by 11_full_model_analyses.R — guarantees an exact match with that script's diagnostic plots.")
  m_cort_full    <- readRDS(rds_paths$m_cort_full)
  m_cort_full_xy <- readRDS(rds_paths$m_cort_full_xy)
  m_t3_full      <- readRDS(rds_paths$m_t3_full)
  m_t3_full_xy   <- readRDS(rds_paths$m_t3_full_xy)
  dat_full_cort  <- readRDS(rds_paths$dat_full_cort)
  dat_full_t3    <- readRDS(rds_paths$dat_full_t3)

} else {

  warning("11_full_model_analyses.R's saved model files were not found in ", outdir,
          " — refitting the full model here instead. Run 11 first (it now saves ",
          "model_11_*.rds) if you want this script's full-model curves to exactly ",
          "match 11's own diagnostic plots; a fresh refit CAN differ slightly from ",
          "11's fit due to collinearity among wolf/mine/road/land cover (see note above).")

  dat_full <- dat_all %>%
    filter(
      is.finite(snow_on_ground), is.finite(snow_on_ground_fcm2),
      is.finite(dist_road_km), is.finite(dist_allroad_km), is.finite(wolf_suit),
      is.finite(dist_mine_km), is.finite(dist_community_km),
      is.finite(LC_PC1), is.finite(LC_PC2), is.finite(LC_PC3), is.finite(LC_PC4),
      !is.na(preg_group)
    )

  # snow_col differs by hormone: FCM uses the 2-day window
  # (snow_on_ground_fcm2), FT3 keeps the 5-day window (snow_on_ground).
  full_formula <- function(response, spatial, snow_col) {
    base <- sprintf(
      "%s ~ s(survey, bs = 're') + year_f + preg_group + s(%s, k = 4, bs = 'tp') + s(log_dist_road_km, k = 5, bs = 'tp') + s(log_dist_allroad_km, k = 5, bs = 'tp') + s(log_wolf_suit, k = 6, bs = 'tp') + s(log_dist_mine, k = 6, bs = 'tp') + s(log_dist_community, k = 6, bs = 'tp') + s(LC_PC1, k = 3, bs = 'tp') + s(LC_PC2, k = 3, bs = 'tp') + s(LC_PC3, k = 3, bs = 'tp') + s(LC_PC4, k = 3, bs = 'tp') + s(julian, k = 8, bs = 'tp')",
      response, snow_col
    )
    if (spatial) base <- paste(base, "+ s(x_proj, y_proj, bs = 'tp', k = 30)")
    as.formula(base)
  }

  dat_full_cort <- dat_full %>% filter(is.finite(cortisol), !is.na(julian))
  dat_full_t3   <- dat_full %>% filter(is.finite(t3), t3 > 0, !is.na(julian))

  m_cort_full    <- gam(full_formula("cortisol", FALSE, "snow_on_ground_fcm2"), data = dat_full_cort, family = Gamma(link = "log"), method = "REML")
  m_cort_full_xy <- gam(full_formula("cortisol", TRUE,  "snow_on_ground_fcm2"), data = dat_full_cort, family = Gamma(link = "log"), method = "REML")
  m_t3_full      <- gam(full_formula("t3",       FALSE, "snow_on_ground"),      data = dat_full_t3,   family = Gamma(link = "log"), method = "REML")
  m_t3_full_xy   <- gam(full_formula("t3",       TRUE,  "snow_on_ground"),      data = dat_full_t3,   family = Gamma(link = "log"), method = "REML")
}

cat("Full-model complete-case n — cortisol:", nrow(dat_full_cort), " | T3:", nrow(dat_full_t3), "\n")

# ============================================================
# 4) EXTRACT AND COMPARE PARTIAL-EFFECT CURVES
# ============================================================

extract_curve <- function(model, smooth_name, xvar, model_label) {
  sm <- smooth_estimates(model, smooth = smooth_name) |> add_confint()
  sm %>%
    transmute(x = .data[[xvar]], estimate = .estimate,
              lower = .lower_ci, upper = .upper_ci,
              model = model_label)
}

build_comparison <- function(xvar_cort, xvar_t3, ind_cort, ind_t3,
                              term_ind_cort, term_ind_t3, term_full_cort, term_full_t3,
                              response_model_cort, response_model_t3) {
  # response_model_cort/t3 are the full (no-space / +space) models.
  # xvar/term are split per hormone since snow uses different columns for
  # FCM (2-day, snow_on_ground_fcm2) vs FT3 (5-day, snow_on_ground); for
  # every other covariate, pass the same value for both _cort and _t3.
  d_cort <- bind_rows(
    extract_curve(ind_cort,           term_ind_cort,  xvar_cort, "Individual model"),
    extract_curve(response_model_cort[[1]], term_full_cort, xvar_cort, "Full model, no space"),
    extract_curve(response_model_cort[[2]], term_full_cort, xvar_cort, "Full model, + space")
  )
  d_t3 <- bind_rows(
    extract_curve(ind_t3,             term_ind_t3,  xvar_t3, "Individual model"),
    extract_curve(response_model_t3[[1]],   term_full_t3, xvar_t3, "Full model, no space"),
    extract_curve(response_model_t3[[2]],   term_full_t3, xvar_t3, "Full model, + space")
  )
  list(cort = d_cort, t3 = d_t3)
}

model_colours <- c("Individual model" = "#000000",
                    "Full model, no space" = "#0072B2",
                    "Full model, + space"  = "#D55E00")

make_panel <- function(curve_df, xlab, ylab, title) {
  ggplot(curve_df, aes(x = x, y = estimate, colour = model, fill = model)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = model_colours, name = "Model") +
    scale_fill_manual(values = model_colours, guide = "none") +
    labs(x = xlab, y = ylab, title = title) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(size = 9, face = "bold"),
          legend.position = "bottom")
}

# Snow: FCM uses the 2-day short-lag column/term, FT3 uses the 5-day one.
curves_snow <- build_comparison("snow_on_ground_fcm2", "snow_on_ground", m_cort_snow_ind, m_t3_snow_ind,
                                 "s(snow_on_ground_fcm2)", "s(snow_on_ground)",
                                 "s(snow_on_ground_fcm2)", "s(snow_on_ground)",
                                 list(m_cort_full, m_cort_full_xy), list(m_t3_full, m_t3_full_xy))
curves_winter_road <- build_comparison("log_dist_road_km", "log_dist_road_km", m_cort_road_ind, m_t3_road_ind,
                                 "s(log_dist_road_km)", "s(log_dist_road_km)",
                                 "s(log_dist_road_km)", "s(log_dist_road_km)",
                                 list(m_cort_full, m_cort_full_xy), list(m_t3_full, m_t3_full_xy))
curves_all_road <- build_comparison("log_dist_allroad_km", "log_dist_allroad_km", m_cort_allroad_ind, m_t3_allroad_ind,
                                 "s(log_dist_allroad_km)", "s(log_dist_allroad_km)",
                                 "s(log_dist_allroad_km)", "s(log_dist_allroad_km)",
                                 list(m_cort_full, m_cort_full_xy), list(m_t3_full, m_t3_full_xy))
curves_wolf <- build_comparison("log_wolf_suit", "log_wolf_suit", m_cort_wolf_ind, m_t3_wolf_ind,
                                 "s(log_wolf_suit)", "s(log_wolf_suit)",
                                 "s(log_wolf_suit)", "s(log_wolf_suit)",
                                 list(m_cort_full, m_cort_full_xy), list(m_t3_full, m_t3_full_xy))
curves_mine <- build_comparison("log_dist_mine", "log_dist_mine", m_cort_mine_ind, m_t3_mine_ind,
                                 "s(log_dist_mine)", "s(log_dist_mine)",
                                 "s(log_dist_mine)", "s(log_dist_mine)",
                                 list(m_cort_full, m_cort_full_xy), list(m_t3_full, m_t3_full_xy))
curves_pop  <- build_comparison("log_dist_community", "log_dist_community", m_cort_pop_ind, m_t3_pop_ind,
                                 "s(log_dist_community)", "s(log_dist_community)",
                                 "s(log_dist_community)", "s(log_dist_community)",
                                 list(m_cort_full, m_cort_full_xy), list(m_t3_full, m_t3_full_xy))
curves_lc1  <- build_comparison("LC_PC1", "LC_PC1", m_cort_lc_ind, m_t3_lc_ind,
                                 "s(LC_PC1)", "s(LC_PC1)", "s(LC_PC1)", "s(LC_PC1)",
                                 list(m_cort_full, m_cort_full_xy), list(m_t3_full, m_t3_full_xy))
curves_lc2  <- build_comparison("LC_PC2", "LC_PC2", m_cort_lc_ind, m_t3_lc_ind,
                                 "s(LC_PC2)", "s(LC_PC2)", "s(LC_PC2)", "s(LC_PC2)",
                                 list(m_cort_full, m_cort_full_xy), list(m_t3_full, m_t3_full_xy))
curves_lc3  <- build_comparison("LC_PC3", "LC_PC3", m_cort_lc_ind, m_t3_lc_ind,
                                 "s(LC_PC3)", "s(LC_PC3)", "s(LC_PC3)", "s(LC_PC3)",
                                 list(m_cort_full, m_cort_full_xy), list(m_t3_full, m_t3_full_xy))
curves_lc4  <- build_comparison("LC_PC4", "LC_PC4", m_cort_lc_ind, m_t3_lc_ind,
                                 "s(LC_PC4)", "s(LC_PC4)", "s(LC_PC4)", "s(LC_PC4)",
                                 list(m_cort_full, m_cort_full_xy), list(m_t3_full, m_t3_full_xy))

all_curve_sets <- list(
  snow = curves_snow, winter_road = curves_winter_road, all_road = curves_all_road,
  wolf = curves_wolf, mine = curves_mine,
  population = curves_pop, LC_PC1 = curves_lc1, LC_PC2 = curves_lc2,
  LC_PC3 = curves_lc3, LC_PC4 = curves_lc4
)

# Snow's xlab differs by hormone (FCM: 2-day window; FT3: 5-day window);
# everything else is shared between hormones.
xlabs_cort <- c(snow = "2-day snow-on-ground (cm)", winter_road = "log(winter road distance, km)",
           all_road = "log(all-year road distance, km)",
           wolf = "log(wolf suitability)", mine = "log(mine distance, km)",
           population = "log(community distance, km)",
           LC_PC1 = "Land cover PC1", LC_PC2 = "Land cover PC2",
           LC_PC3 = "Land cover PC3", LC_PC4 = "Land cover PC4")
xlabs_t3 <- xlabs_cort
xlabs_t3["snow"] <- "5-day snow-on-ground (cm)"

panel_titles <- c(snow = "snow", winter_road = "winter road", all_road = "all-year road",
                   wolf = "wolf", mine = "mine", population = "population",
                   LC_PC1 = "LC_PC1", LC_PC2 = "LC_PC2", LC_PC3 = "LC_PC3", LC_PC4 = "LC_PC4")

panels_cort <- lapply(names(all_curve_sets), function(nm) {
  make_panel(all_curve_sets[[nm]]$cort, xlabs_cort[[nm]], "Partial effect on log(FCM)", panel_titles[[nm]])
})
panels_t3 <- lapply(names(all_curve_sets), function(nm) {
  make_panel(all_curve_sets[[nm]]$t3, xlabs_t3[[nm]], "Partial effect on log(FT3)", panel_titles[[nm]])
})

figS_cort <- wrap_plots(panels_cort, ncol = 3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Supplementary Figure 12. FCM partial-effect curve shape: individual vs. full-model specifications",
    caption = "Each panel compares the same covariate's partial-effect curve from (i) its own individual spatially-adjusted model (03/05/06/07/08 spec; a fresh snow-only+space model for snow, see Methods note), (ii) the seven-covariate full model without a spatial term, and (iii) the full model with s(x, y). Winter and all-year roads are shown as separate panels. Shaded bands are 95% CI.",
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.caption = element_text(size = 8, colour = "grey40"))
  ) &
  theme(legend.position = "bottom")

figS_t3 <- wrap_plots(panels_t3, ncol = 3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Supplementary Figure 13. FT3 partial-effect curve shape: individual vs. full-model specifications",
    caption = "Same comparison as Supplementary Figure 12a, for FT3.",
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.caption = element_text(size = 8, colour = "grey40"))
  ) &
  theme(legend.position = "bottom")

# Height increased from 14 to 18: 10 panels at ncol=3 now wrap to 4 rows
# (was 9 panels / 3 rows), after adding the all-year road panel.
ggsave(file.path(outdir_ms, "SupplementaryFigure12_FCM_CurveComparison.png"),
       figS_cort, width = 13, height = 18, dpi = 300)
ggsave(file.path(outdir_ms, "SupplementaryFigure13_FT3_CurveComparison.png"),
       figS_t3, width = 13, height = 18, dpi = 300)

# Also save the underlying curve data, in case a different plotting
# style is wanted later.
curve_data_long <- bind_rows(
  lapply(names(all_curve_sets), function(nm) {
    bind_rows(
      all_curve_sets[[nm]]$cort %>% mutate(covariate = nm, hormone = "FCM"),
      all_curve_sets[[nm]]$t3   %>% mutate(covariate = nm, hormone = "FT3")
    )
  })
)
write_csv(curve_data_long, file.path(outdir, "Curve_comparison_data_12.csv"))

message("\nPart 12 (curve-shape comparison figures) complete — outputs in: ", outdir_ms)
