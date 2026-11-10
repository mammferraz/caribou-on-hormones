# ============================================================
# 00_data_processing.R
#
# PURPOSE: Shared data loading, cleaning, and preparation for
#   all Caribou Ontario analyses. Source this file at the top
#   of every analysis script:
#
#     source("00_data_processing.R")
#
# After sourcing, the following objects are available:
#   dat          — full cleaned dataset (all animals)
#   dat_map      — subset with valid GPS coordinates
#   caribou_ranges — sf object: boreal caribou range polygons
#   outdir       — output directory (created if absent)
#
# Primary temporal variable: year (integer)
# Survey/range: assigned from spatial boundary overlap, used
#   as random effect in GAMs and colour in plots.
# ============================================================

setwd("/Users/marciaferraz/Desktop/Toronto Zoo/Collaboration/Caribou_ON/Final analyses")

# ============================================================
# PACKAGES
# ============================================================

library(tidyverse)
library(lubridate)
library(janitor)
library(sf)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(scales)
library(viridis)
library(units)
library(mgcv)
library(emmeans)
library(gratia)
library(patchwork)
library(multcompView)

# ============================================================
# INPUT FILES
# ============================================================

infile        <- "All_data_updated.csv"
infile_ranges <- "Caribou_range_boundary.shp"   # boreal caribou range polygons (NAD83)

# Output directory — each analysis script may override this
outdir <- "caribou_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# LOAD & CLEAN BIOLOGICAL DATA
# ============================================================

dat_raw <- suppressWarnings(read_csv(infile, show_col_types = FALSE)) %>%
  janitor::clean_names() %>%
  select(-any_of("...24"))      # drop trailing empty column if present

dat <- dat_raw %>%
  mutate(
    sex                  = toupper(trimws(sex)),
    survey               = trimws(survey),
    survey_field         = trimws(survey),  # original field campaign label — used for climate station matching
    pcid                 = as.character(pcid),
    date                 = suppressWarnings(as.Date(date)),
    year                 = as.integer(year),
    month                = as.integer(month),
    day                  = as.integer(day),
    cortisol             = as.numeric(cortisol),
    progesterone         = as.numeric(progesterone),
    t3                   = as.numeric(t3),
    latitude             = as.numeric(latitude),
    longitude            = as.numeric(longitude),
    pregnancy_status_num = suppressWarnings(as.numeric(pregnancy_status))
  ) %>%
  # Fill julian from Year/Month/Day where Date is NA
  mutate(
    julian = if_else(
      is.na(julian) & !is.na(year) & !is.na(month) & !is.na(day),
      as.numeric(yday(make_date(year, month, day))),
      as.numeric(julian)
    )
  )

# ============================================================
# PREGNANCY GROUP
# ============================================================

dat <- dat %>%
  mutate(
    preg_group = case_when(
      sex == "M"                            ~ "Male",
      sex == "F" & pregnancy_status_num == 1 ~ "Female_pregnant",
      sex == "F"                            ~ "Female_nonpregnant",
      TRUE                                  ~ NA_character_
    ),
    preg_group = factor(
      preg_group,
      levels = c("Male", "Female_nonpregnant", "Female_pregnant")
    )
  )

# Spatial subset — animals with GPS coordinates
dat_map <- dat %>% filter(!is.na(latitude), !is.na(longitude))

# ============================================================
# SPATIAL SURVEY REASSIGNMENT
# Each animal's survey label is replaced by the name of the
# boreal caribou range polygon its coordinates fall within.
# This ensures survey groupings reflect boundary geography,
# not field-campaign labels. Animals outside all polygons
# retain their original survey label.
# ============================================================

caribou_ranges <- st_read(infile_ranges, quiet = TRUE) %>%
  st_transform(4326)

# Identify the range name column
range_name_col <- intersect(
  c("RANGE_NAME", "Range_Name", "range_name", "RANGENAME",
    "NAME", "name", "HERD_NAME", "herd_name"),
  names(caribou_ranges)
)
if (length(range_name_col) == 0) {
  stop("No range name column found in ", infile_ranges,
       ". Available: ", paste(names(caribou_ranges), collapse = ", "))
}
range_name_col <- range_name_col[1]
message("Range name column: ", range_name_col)

# Spatial join: each animal point → range polygon it falls within
dat_sf <- dat_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

joined <- st_join(
  dat_sf,
  caribou_ranges %>% select(range_name_assigned = all_of(range_name_col)),
  join = st_within,
  left = TRUE
)

survey_lookup <- joined %>%
  st_drop_geometry() %>%
  select(pcid, latitude, longitude, range_name_assigned) %>%
  # Remove spaces so "James Bay" matches "JamesBay" etc.
  mutate(range_name_assigned = gsub(" ", "", range_name_assigned)) %>%
  # Deduplicate: a point may fall within overlapping polygons; keep first match
  distinct(pcid, latitude, longitude, .keep_all = TRUE)

# Report animals outside all polygons
outside <- survey_lookup %>% filter(is.na(range_name_assigned))
if (nrow(outside) > 0) {
  message(nrow(outside), " animal(s) outside all range polygons — retaining original survey label.")
}

# Report any label reassignments
reassigned <- dat_map %>%
  left_join(survey_lookup, by = c("pcid", "latitude", "longitude")) %>%
  filter(!is.na(range_name_assigned),
         gsub(" ", "", survey) != range_name_assigned) %>%
  count(survey, range_name_assigned, name = "n")
if (nrow(reassigned) > 0) {
  message("Survey reassignments (original → boundary-assigned):")
  print(reassigned)
}

# Apply reassignment to full dataset
dat <- dat %>%
  left_join(survey_lookup, by = c("pcid", "latitude", "longitude")) %>%
  mutate(
    survey = factor(
      gsub(" ", "",
        if_else(!is.na(range_name_assigned),
                range_name_assigned,
                gsub(" ", "", survey))
      )
    )
  ) %>%
  select(-range_name_assigned)

# Rebuild dat_map with corrected labels
dat_map <- dat %>% filter(!is.na(latitude), !is.na(longitude))

message("Sample counts by survey and year:")
print(dat %>% count(survey, year) %>% arrange(year, survey))

# ============================================================
# NEAREST CLIMATE STATION (per-sample, by GPS distance)
#
# Climate is matched to whichever station is geographically nearest to
# each individual sample's own coordinates — NOT by field campaign
# (survey_field) and NOT by boundary-assigned range (survey). Per-range
# or per-campaign lookup was checked and found too coarse: e.g. Missisa
# is large enough that ~38% of its samples are geographically closer to
# James Bay's station (Ogoki Post) than to its own nominal station
# (Lansdowne House), and all of Pagwachuan's samples are closer to Ogoki
# Post than to Armstrong, the station its field campaign nominally used.
# This produces one climate_group value per sample, computed directly,
# rather than a per-campaign/per-range table.
#
# climatedaily_all_surveys.csv station coordinates (from the file itself):
#   Ogoki Post            -85.90029639, 51.65779361
#   Lansdowne House (AUT) -87.93611111, 52.19611111
#   Pickle Lake (AUT)     -90.21777972, 51.44916694
#   Armstrong (AUT)       -88.90527778, 50.29361111
# Geraldton A (proxy for Missisa 2021 snow-on-ground only; see
# 02_climate_analyses.R) is NOT a candidate station here — it only
# supplies snow-on-ground for Lansdowne House's own climate_group.
#
# 5 samples lack GPS coordinates; these fall back to their field
# campaign's nominal station.
# ============================================================

climate_stations <- tribble(
  ~climate_group,              ~stn_lon,      ~stn_lat,
  "JamesBay",                  -85.90029639,  51.65779361,
  "Missisa",                   -87.93611111,  52.19611111,
  "Ozhiski",                   -90.21777972,  51.44916694,
  "NipigonKinlochPagwachuan",  -88.90527778,  50.29361111
)
# NOTE: "NipigonKinlochPagwachuan" has no hyphens deliberately —
# multcompLetters() (used in 02_climate_analyses.R for Tukey-group
# lettering) splits contrast names on "-", so a group label containing
# hyphens breaks it (e.g. "JamesBay-Nipigon-Kinloch-Pagwachuan" parses
# as too many groups). Matches the renamed value now used in
# climatedaily_all_surveys.csv's climate_group column.

haversine_km <- function(lon1, lat1, lon2, lat2) {
  R <- 6371
  p1 <- lat1 * pi / 180; p2 <- lat2 * pi / 180
  dphi    <- (lat2 - lat1) * pi / 180
  dlambda <- (lon2 - lon1) * pi / 180
  a <- sin(dphi / 2)^2 + cos(p1) * cos(p2) * sin(dlambda / 2)^2
  2 * R * asin(sqrt(a))
}

# Station-year availability (Jan-Mar winter window), from
# climatedaily_all_surveys.csv as actually assembled. Ogoki Post
# (JamesBay's station) only began reporting in Nov 2022, so it has NO
# winter data for 2021 or 2022 — a genuine station-history gap, not a
# download error. All other stations have complete 2021-2024 winter
# coverage. nearest_climate_station() is therefore year-aware: for a
# given sample-year, it picks the nearest station THAT HAS DATA that
# year, rather than the nearest station unconditionally (which would
# silently drop 224 samples — 169 in 2021, 55 in 2022 — that are
# geographically nearest to Ogoki Post but have no year-matched data).
station_year_available <- tribble(
  ~climate_group,              ~year_available,
  "JamesBay",                  2023,
  "JamesBay",                  2024,
  "Missisa",                   2021,
  "Missisa",                   2022,
  "Missisa",                   2023,
  "Missisa",                   2024,
  "Ozhiski",                   2021,
  "Ozhiski",                   2022,
  "Ozhiski",                   2023,
  "Ozhiski",                   2024,
  "NipigonKinlochPagwachuan",  2021,
  "NipigonKinlochPagwachuan",  2022,
  "NipigonKinlochPagwachuan",  2023,
  "NipigonKinlochPagwachuan",  2024
)

nearest_climate_station <- function(lon, lat, yr) {
  if (is.na(lon) || is.na(lat) || is.na(yr)) return(NA_character_)
  available <- station_year_available$climate_group[station_year_available$year_available == yr]
  candidates <- climate_stations[climate_stations$climate_group %in% available, ]
  if (nrow(candidates) == 0) return(NA_character_)
  d <- vapply(
    seq_len(nrow(candidates)),
    function(i) haversine_km(lon, lat, candidates$stn_lon[i], candidates$stn_lat[i]),
    numeric(1)
  )
  candidates$climate_group[which.min(d)]
}

# Fallback for the handful of samples with no GPS coordinates: use their
# original field campaign's nominal station.
fallback_climate_map <- tribble(
  ~survey_field, ~climate_group_fallback,
  "Ozhiski",     "Ozhiski",
  "Missisa",     "Missisa",
  "Nipigon",     "NipigonKinlochPagwachuan",
  "Pagwachuan",  "NipigonKinlochPagwachuan",
  "JamesBay",    "JamesBay"
)

dat <- dat %>%
  rowwise() %>%
  mutate(climate_group = nearest_climate_station(longitude, latitude, year)) %>%
  ungroup() %>%
  left_join(fallback_climate_map, by = "survey_field") %>%
  mutate(climate_group = if_else(is.na(climate_group), climate_group_fallback, climate_group)) %>%
  select(-climate_group_fallback)

dat_map <- dat %>% filter(!is.na(latitude), !is.na(longitude))

message("Samples per nearest climate station (climate_group):")
print(dat %>% count(climate_group))
message("Cross-tab: assigned range (survey) x nearest climate station:")
print(dat %>% count(survey, climate_group) %>% arrange(survey, climate_group))

# ============================================================
# SHARED PLOT THEME & AESTHETICS
# ============================================================

theme_pub <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold"),
    legend.position  = "right"
  )

# Colour palette — Okabe-Ito colorblind-safe palette
# doi:10.1038/nmeth.1618 / https://jfly.uni-koeln.de/color/
survey_colours <- c(
  "Missisa"    = "#0072B2",   # blue
  "Ozhiski"   = "#E69F00",   # orange
  "JamesBay"  = "#009E73",   # bluish green
  "Nipigon"   = "#CC79A7",   # mauve/pink
  "Pagwachuan"= "#56B4E9",   # sky blue
  "Kinloch"   = "#D55E00"    # vermillion
)

preg_shapes <- c(
  "Male"               = 15,
  "Female_nonpregnant" = 17,
  "Female_pregnant"    = 16
)

preg_colours <- c(
  "Male"               = "#999999",   # grey   (Okabe-Ito)
  "Female_nonpregnant" = "#0072B2",   # blue
  "Female_pregnant"    = "#D55E00"    # vermillion
)

