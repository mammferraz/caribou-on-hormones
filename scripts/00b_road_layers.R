# ============================================================
# 00b_road_layers.R
#
# SHARED ROAD LAYERS + MEMORY-SAFE DISTANCE HELPER
#
# WHY THIS FILE EXISTS
# --------------------
# Road loading and distance-to-road code was previously duplicated across
# 03, 03b, 04, 09, 10, 11 and 12. Those copies drifted apart, and two
# errors propagated silently through all of them:
#
#   1. WRONG ALL-SEASON ROAD DEFINITION.
#      "All-season road" was defined as NATIONAL_ROAD_CLASS ==
#      "Local / Strata" only. That is 8,648 of 517,836 segments (1.7%) —
#      a narrow ownership class, not a road network. The actual
#      non-winter network in this dataset is "Resource / Recreation"
#      (508,220 segments), which had been excluded on the untested
#      assumption that it was "largely seasonal/closed". It is not:
#      343,960 Resource/Recreation segments have STATUS = "Open" and
#      312,616 have PASSABLE_IND = "Yes".
#      All-season roads are now ALL NON-WINTER segments (516,868).
#
#   2. INCOMPLETE WINTER ROAD LAYER.
#      Winter roads were read from MNR_Road_Segments.geojson, which holds
#      541 Winter segments. The provincial file MNR_Road_Segments2.geojson
#      holds 968, and the 541 are a strict subset of them. Seven of the
#      427 missing segments fall inside the sampling window, so winter-road
#      distances were biased. BOTH layers now come from the provincial file.
#
# MEMORY SAFETY
# -------------
# The old code called st_distance(points, roads), which materialises an
# n x m matrix. With 1,375 points and 516,868 roads that is ~5.7 GB and
# will hit R's vector memory limit. nearest_road_km() instead uses
# st_nearest_feature() (spatially indexed) and then computes only the
# n pairwise distances — O(n) memory.
#
# CRS NOTE
# --------
# Distances are computed in EPSG:3161 (NAD83 / Ontario MNR Lambert,
# metres), which is the appropriate projected CRS for this region.
# This file does NOT touch the x_proj / y_proj coordinates used by the
# s(x, y) smooths in the analysis scripts — those are computed
# independently in each script and are deliberately left unchanged so
# that swapping the road layers does not silently alter every spatial
# smooth in the project.
#
# USAGE
#   source("00b_road_layers.R")
#   winter <- road_layer("winter")
#   allrd  <- road_layer("allseason")
#   d_km   <- nearest_road_km(pts_sf, winter)
#
# The first call parses the 1.3 GB GeoJSON and writes a GeoPackage cache
# to caribou_outputs/. Later calls read the cache (seconds, not minutes).
# Use road_layer(..., refresh = TRUE) to rebuild.
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
})

ROADS_SOURCE   <- "MNR_Road_Segments2.geojson"   # provincial, all classes
ROADS_CACHE_DIR <- "caribou_outputs"
ROADS_DIST_CRS <- 3161                            # NAD83 / Ontario MNR Lambert (m)

# Buffer (degrees) around the sample point cloud used when caching roads.
# ~2 deg latitude is roughly 220 km, far beyond any plausible
# nearest-road distance in a network this dense, so clipping here cannot
# bias the computed distances. It exists purely to keep the cached layers
# a workable size.
ROADS_CLIP_BUFFER_DEG <- 2

.roads_bbox <- function() {
  if (!exists("dat_map", envir = .GlobalEnv)) {
    stop("dat_map not found — source('00_data_processing.R') before 00b_road_layers.R")
  }
  dm <- get("dat_map", envir = .GlobalEnv)
  b  <- sf::st_bbox(sf::st_as_sf(dm, coords = c("longitude", "latitude"), crs = 4326))
  k  <- ROADS_CLIP_BUFFER_DEG
  sf::st_as_sfc(sf::st_bbox(c(xmin = unname(b["xmin"]) - k, xmax = unname(b["xmax"]) + k,
                              ymin = unname(b["ymin"]) - k, ymax = unname(b["ymax"]) + k),
                            crs = 4326))
}

#' Load the winter or all-season road layer (cached).
#'
#' @param type "winter" (NATIONAL_ROAD_CLASS == "Winter") or "allseason"
#'   (every non-Winter class: Resource / Recreation + Local / Strata).
#' @param refresh rebuild the cache from the source GeoJSON.
road_layer <- function(type = c("winter", "allseason"), refresh = FALSE) {
  type  <- match.arg(type)
  dir.create(ROADS_CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
  cache <- file.path(ROADS_CACHE_DIR, paste0("roads_", type, ".gpkg"))

  if (!refresh && file.exists(cache)) {
    out <- sf::st_read(cache, quiet = TRUE)
    message(sprintf("road_layer('%s'): %s segments (from cache)", type,
                    format(nrow(out), big.mark = ",")))
    return(out)
  }

  if (!file.exists(ROADS_SOURCE)) {
    stop("Road source not found: ", ROADS_SOURCE)
  }
  message("road_layer('", type, "'): parsing ", ROADS_SOURCE,
          " — this takes a minute the first time...")
  rds <- sf::st_read(ROADS_SOURCE, quiet = TRUE)

  rc_col <- names(rds)[grepl("national_road_class", names(rds), ignore.case = TRUE)][1]
  if (is.na(rc_col)) stop("NATIONAL_ROAD_CLASS column not found in ", ROADS_SOURCE)

  is_winter <- rds[[rc_col]] == "Winter"
  rds <- if (type == "winter") rds[is_winter, ] else rds[!is_winter, ]

  rds <- sf::st_transform(rds, 4326)
  rds <- rds[!sf::st_is_empty(sf::st_geometry(rds)), ]

  # Clip to the sampling window + buffer (see ROADS_CLIP_BUFFER_DEG note).
  n_before <- nrow(rds)
  clipped <- tryCatch({
    old_s2 <- sf::sf_use_s2(); suppressMessages(sf::sf_use_s2(FALSE))
    keep <- lengths(suppressWarnings(sf::st_intersects(rds, .roads_bbox()))) > 0
    suppressMessages(sf::sf_use_s2(old_s2))
    rds[keep, ]
  }, error = function(e) {
    message("  ⚠ clip failed (", e$message, ") — caching the full layer instead.")
    NULL
  })
  if (!is.null(clipped) && nrow(clipped) > 0) rds <- clipped

  # Keep only the few attributes anything downstream actually uses, so the
  # cache stays small and reads fast.
  keep_cols <- intersect(c("OGF_ID", "ogf_id", rc_col, "STATUS", "status",
                           "MAINTENANCE_LEVEL", "maintenance_level",
                           "PASSABLE_IND", "passable_ind", "ROAD_NAME", "road_name"),
                         names(rds))
  rds <- rds[, keep_cols]

  sf::st_write(rds, cache, delete_dsn = TRUE, quiet = TRUE)
  message(sprintf("road_layer('%s'): %s segments kept of %s (clipped to study window); cached to %s",
                  type, format(nrow(rds), big.mark = ","),
                  format(n_before, big.mark = ","), cache))
  rds
}

#' Distance (km) from each point to the nearest road segment.
#'
#' Memory-safe: uses a spatial index to find the nearest segment per point,
#' then computes only n pairwise distances. Never builds an n x m matrix.
#'
#' @param pts_sf sf POINT object (any CRS).
#' @param roads_sf sf LINESTRING object from road_layer().
#' @return numeric vector of distances in km, same length/order as pts_sf.
nearest_road_km <- function(pts_sf, roads_sf) {
  stopifnot(inherits(pts_sf, "sf"), inherits(roads_sf, "sf"))
  if (nrow(roads_sf) == 0) stop("nearest_road_km(): road layer is empty.")

  p <- sf::st_transform(pts_sf,  ROADS_DIST_CRS)
  r <- sf::st_transform(roads_sf, ROADS_DIST_CRS)

  idx <- sf::st_nearest_feature(p, r)
  d   <- sf::st_distance(p, r[idx, ], by_element = TRUE)
  as.numeric(d) / 1000
}

message("00b_road_layers.R loaded — road_layer() and nearest_road_km() available.")
