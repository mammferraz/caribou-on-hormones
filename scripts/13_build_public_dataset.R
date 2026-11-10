# ============================================================
# 13_build_public_dataset.R
#
# BUILD THE DE-IDENTIFIED, ANALYSIS-READY DATASET FOR PUBLIC RELEASE
#
# WHY
# ---
# Boreal caribou are listed as Threatened under Canada's Species at Risk
# Act, and precise faecal collection coordinates could be used to locate
# animals. Latitude and longitude are therefore WITHHELD from the public
# release.
#
# Simply deleting the coordinate columns would make the pipeline
# irreproducible, because every landscape covariate is derived FROM those
# coordinates. This script instead computes every spatially derived
# covariate at the true coordinates and then writes out only the derived
# values, so all reported models can be refit exactly without the
# coordinates ever leaving this machine.
#
# WHAT IS AND IS NOT REPRODUCIBLE FROM THE PUBLIC FILE
#   Reproducible : every model in scripts 01-08 and 11-12 that uses the
#                  derived covariates (road / mine / settlement distances,
#                  wolf suitability, land cover PCs, lagged climate),
#                  plus all year, range, reproductive-group and Julian-day
#                  terms.
#   NOT reproducible : the two-dimensional spatial smooth s(x, y), which
#                  needs the coordinates themselves. Fitted mgcv objects are
#                  deliberately NOT published either: gam() stores its model
#                  frame inside the object, so any model containing
#                  s(x_proj, y_proj) carries the full coordinate table and
#                  m$model$x_proj would recover every collection location.
#                  The rendered partial-effect surfaces and the .txt model
#                  summaries are published instead. This limitation is stated
#                  in the README and in the Data Availability statement.
#
# OUTPUT
#   caribou_outputs/public/caribou_hormones_analysis_ready.csv
#   caribou_outputs/public/DATA_DICTIONARY.md
#
# Run AFTER 00-08 so the cached road layers and rasters are available.
# ============================================================

source("00_data_processing.R")
source("00b_road_layers.R")

suppressPackageStartupMessages({
  library(terra)
})

outdir_pub <- file.path("caribou_outputs", "public")
dir.create(outdir_pub, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1) SPATIAL COVARIATES (computed at the TRUE coordinates)
# ============================================================

pts <- dat_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

message("Computing road distances...")
dist_road_km    <- nearest_road_km(pts, road_layer("winter"))
dist_allroad_km <- nearest_road_km(pts, road_layer("allseason"))

message("Extracting wolf habitat suitability...")
wolf_r    <- terra::rast("wolfmean_t.tif")
wolf_suit <- terra::extract(wolf_r, vect(st_transform(pts, crs(wolf_r))))[, 2]

message("Computing mine distances...")
mines_r    <- terra::rast("mines.tif")
pts_mn     <- st_transform(pts, crs(mines_r))
mines_crop <- terra::crop(mines_r, vect(st_buffer(pts_mn, 500000)))
mine_bin   <- mines_crop; mine_bin[mine_bin == 0] <- NA
dist_mine_km <- terra::extract(terra::distance(mine_bin), vect(pts_mn))[, 2] / 1000

message("Computing settlement distances...")
pop_r    <- terra::rast("population_density.tif")
pts_pop  <- st_transform(pts, crs(pop_r))
pop_crop <- terra::crop(pop_r, vect(st_buffer(pts_pop, 500000)))
pop_bin  <- pop_crop
pop_bin[pop_bin <= 0] <- NA
pop_bin[!is.na(pop_bin)] <- 1
patches  <- terra::patches(pop_bin, directions = 8)
sizes    <- terra::freq(patches)
min_cells <- max(1, round(1 / (prod(res(pop_crop)) / 1e6)))
small    <- sizes$value[sizes$count < min_cells]
if (length(small) > 0) patches[patches %in% small] <- NA
comm     <- terra::centroids(terra::as.polygons(patches)) %>% st_as_sf()
dist_community_km <- as.numeric(apply(st_distance(pts_pop, comm), 1, min)) / 1000

message("Extracting land cover and computing PCA axes...")
wolfstack <- terra::rast("wolfstack.tif")
pts_lc    <- st_transform(pts, crs(wolfstack))
lc_vals   <- terra::extract(wolfstack, vect(pts_lc)) %>% select(-ID)
lc_all <- c("Water","Wetlands","OpenPeat","TreedPeat","Coniferous","Sparse",
            "Deciduous","NonTreed","Anthropogenic","Disturbance","OpenBog",
            "TreedBog","OpenFen","TreedFen","ThicketSwamp","ConiferousSwamp",
            "DeciduousSwamp","Marsh")
bands_present <- intersect(lc_all, names(wolfstack))
crop_lc <- terra::crop(wolfstack[[bands_present]], vect(st_buffer(pts_lc, 50000)))
active  <- bands_present[terra::global(crop_lc, "mean", na.rm = TRUE)[, 1] > 0.005]

lc_ok  <- !apply(lc_vals[, active, drop = FALSE], 1, anyNA)
pca_lc <- prcomp(lc_vals[lc_ok, active], scale. = TRUE)
LC <- matrix(NA_real_, nrow = nrow(lc_vals), ncol = 4,
             dimnames = list(NULL, paste0("LC_PC", 1:4)))
LC[lc_ok, ] <- pca_lc$x[, 1:4]

# ============================================================
# 2) CLIMATE LAGS (2-day for FCM, 5-day for FT3)
# ============================================================

message("Computing climate lags...")
clim <- read_csv("climatedaily_all_surveys.csv", show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(date = as.Date(local_date), month = month(date)) %>%
  filter(month %in% c(1, 2, 3))
clim_split <- split(clim, clim$climate_group)

lag_mean <- function(d, cdf, a, b) {
  cdf %>% filter(date < d - days(a - 1), date >= d - days(b)) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
}
lags <- dat_map %>%
  select(pcid, date, tz_sample_number, climate_group) %>%
  rowwise() %>%
  mutate(l5 = list(lag_mean(date, clim_split[[climate_group]], 1, 5)),
         l2 = list(lag_mean(date, clim_split[[climate_group]], 1, 2))) %>%
  ungroup() %>%
  unnest(l5, names_sep = "_") %>%
  unnest(l2, names_sep = "_") %>%
  transmute(pcid, date, tz_sample_number,
            min_temp_lag5_c   = l5_min_temperature,
            snow_lag5_cm      = l5_snow_on_ground,
            min_temp_lag2_c   = l2_min_temperature,
            snow_lag2_cm      = l2_snow_on_ground)

# ============================================================
# 3) ASSEMBLE — NO COORDINATES
# ============================================================

public <- dat_map %>%
  mutate(
    dist_winterroad_km = dist_road_km,
    dist_allroad_km    = dist_allroad_km,
    wolf_suitability   = wolf_suit,
    dist_mine_km       = dist_mine_km,
    dist_settlement_km = dist_community_km
  ) %>%
  bind_cols(as_tibble(LC)) %>%
  left_join(lags, by = c("pcid", "date", "tz_sample_number")) %>%
  transmute(
    year, julian,
    range              = as.character(survey),
    reproductive_group = as.character(preg_group),
    fcm_ng_g           = cortisol,
    ft3_ng_g           = t3,
    fpm_ng_g           = progesterone,
    dist_winterroad_km, dist_allroad_km, dist_mine_km,
    dist_settlement_km, wolf_suitability,
    LC_PC1, LC_PC2, LC_PC3, LC_PC4,
    min_temp_lag2_c, snow_lag2_cm, min_temp_lag5_c, snow_lag5_cm
  )

# ---- LOCATION PROTECTION: coarsen derived covariates and shuffle rows ----
# Distances to four known features at full precision would let anyone
# triangulate each sample's location, so every derived covariate is
# coarsened before release. Rows are shuffled (fixed seed) so that
# consecutive sample IDs do not mark samples from the same site.
public <- public %>%
  mutate(
    across(c(dist_winterroad_km, dist_allroad_km, dist_mine_km,
             dist_settlement_km), ~ round(.x, 0)),
    wolf_suitability = round(wolf_suitability, 3),
    across(c(LC_PC1, LC_PC2, LC_PC3, LC_PC4), ~ round(.x, 1))
  )
set.seed(2026)
public <- public[sample(nrow(public)), ]
public$sample_id <- paste0("S", sprintf("%04d", seq_len(nrow(public))))
public <- public %>% select(sample_id, everything())

# ---- HARD GUARD: nothing coordinate-like may escape ----
banned <- c("latitude", "longitude", "lat", "lon", "long", "x_proj", "y_proj",
            "coords", "geometry", "easting", "northing", "pcid",
            "tz_sample_number", "site_number", "archive_label", "date")
leaked <- names(public)[tolower(names(public)) %in% banned]
if (length(leaked) > 0) {
  stop("Refusing to write: identifying columns present -> ",
       paste(leaked, collapse = ", "))
}
if (any(grepl("lat|lon|coord|geom", tolower(names(public))))) {
  stop("Refusing to write: a column name still looks coordinate-like.")
}
message("Guard passed — no coordinate or individual-identifier columns present.")

write_csv(public, file.path(outdir_pub, "caribou_hormones_analysis_ready.csv"))
message("Wrote ", nrow(public), " rows x ", ncol(public), " columns to ",
        file.path(outdir_pub, "caribou_hormones_analysis_ready.csv"))

# ============================================================
# 4) DATA DICTIONARY
# ============================================================

dict <- c(
"# Data dictionary — caribou_hormones_analysis_ready.csv", "",
"One row per faecal sample with GPS coordinates (n = 1,375 of 1,380 samples; the five",
"samples without coordinates have no derived covariates). Collection coordinates are",
"deliberately omitted; see README.md and the manuscript's Data Availability statement.", "",
"**Coarsened to protect animal locations:** distances are rounded to the nearest km,",
"wolf suitability to 3 decimals and land cover PC scores to 1 decimal; rows are in",
"random order. Refits from this file therefore differ slightly from published values.", "",
"| Column | Units | Description |",
"|---|---|---|",
"| sample_id | — | Anonymous sample identifier. Does not encode location or individual. |",
"| year | — | Sampling year (2021-2024). |",
"| julian | day | Day of year of collection. |",
"| range | — | Caribou range, assigned by spatial overlap with published boundaries. |",
"| reproductive_group | — | Male, Female_nonpregnant or Female_pregnant. |",
"| fcm_ng_g | ng/g | Faecal cortisol metabolites. |",
"| ft3_ng_g | ng/g | Faecal triiodothyronine metabolites. |",
"| fpm_ng_g | ng/g | Faecal progesterone metabolites (used for pregnancy assignment). |",
"| dist_winterroad_km | km | Distance to nearest winter road (rounded to 1 km). |",
"| dist_allroad_km | km | Distance to nearest all-season road (all non-winter classes; rounded to 1 km). |",
"| dist_mine_km | km | Distance to nearest mine (rounded to 1 km). |",
"| dist_settlement_km | km | Distance to nearest derived settlement centroid (rounded to 1 km). |",
"| wolf_suitability | 0-1 | Predicted wolf habitat suitability (3 decimals). |",
"| LC_PC1-LC_PC4 | — | Land cover principal component scores (14 classes, 1-km window; 1 decimal). |",
"| min_temp_lag2_c | °C | Mean daily minimum temperature, 1-2 d before collection (FCM models). |",
"| snow_lag2_cm | cm | Mean snow on ground, 1-2 d before collection (FCM models). |",
"| min_temp_lag5_c | °C | Mean daily minimum temperature, 1-5 d before collection (FT3 models). |",
"| snow_lag5_cm | cm | Mean snow on ground, 1-5 d before collection (FT3 models). |",
"",
"## Not included", "",
"Latitude, longitude, projected coordinates, individual genetic identifiers and",
"collection dates are withheld. The two-dimensional spatial smooth s(x, y) used",
"in the published models therefore cannot be refit from this file. Fitted mgcv",
"objects are NOT published either: gam() stores its model frame inside the object,",
"so any model containing s(x_proj, y_proj) would carry the full coordinate table.",
"The .txt model summaries and the rendered partial-effect surfaces are published",
"instead, in outputs/.")
writeLines(dict, file.path(outdir_pub, "DATA_DICTIONARY.md"))
message("Wrote DATA_DICTIONARY.md")
