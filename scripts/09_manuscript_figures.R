# ============================================================
# 09_manuscript_figures.R
#
# PURPOSE: All publication-ready main manuscript figures
#   Figure 1 — Sample structure (map + hormone profiles + baselines)
#   Figure 2 — Short-lag climate effects on FCM and FT3
#   Figure 3 — Land cover effects (2×4 named-class models; biplot → FigS6)
#   Figure 4 — Distance to nearest community + mining proximity
#   Figure 5 — Wolf habitat suitability
#
# Run from:  Final analyses/   (setwd handled by 00_data_processing.R)
# Output:    caribou_outputs/manuscript/Figure1.png … Figure5.png
# ============================================================

source("00_data_processing.R")
source("00b_road_layers.R")

library(gratia)
library(patchwork)
library(ggrepel)

outdir_ms <- file.path("caribou_outputs", "manuscript")
dir.create(outdir_ms, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# SHARED THEME & HELPERS
# ============================================================

theme_ms <- theme_classic(base_size = 10) +
  theme(
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold", size = 9),
    legend.key.size   = unit(0.4, "cm"),
    legend.text       = element_text(size = 8),
    legend.title      = element_text(size = 9, face = "bold"),
    axis.title        = element_text(size = 9),
    axis.text         = element_text(size = 8),
    plot.tag          = element_text(face = "bold", size = 11),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank()
  )

preg_labels <- c(
  "Male"               = "Male",
  "Female_nonpregnant" = "Non-pregnant female",
  "Female_pregnant"    = "Pregnant female"
)

year_colours <- c(
  "2021" = "#2166ac",
  "2022" = "#4dac26",
  "2023" = "#d6604d",
  "2024" = "#762a83"
)

# Helper: extract smooth estimate from a GAM as ggplot-ready data frame
smooth_df <- function(model, smooth_name, newcol = NULL) {
  sm <- smooth_estimates(model, smooth = smooth_name) |> add_confint()
  if (!is.null(newcol)) names(sm)[names(sm) == smooth_name] <- newcol
  sm
}

# Helper: make a partial-effect panel
partial_panel <- function(sm_data, xvar, xlab, ylab = "Partial effect",
                          raw_data = NULL, raw_xvar = NULL) {
  p <- ggplot(sm_data, aes(x = .data[[xvar]])) +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = "grey80", alpha = 0.6) +
    geom_line(aes(y = .estimate), colour = "black", linewidth = 0.8) +
    labs(x = xlab, y = ylab) +
    theme_ms
  if (!is.null(raw_data) && !is.null(raw_xvar)) {
    p <- p + geom_rug(data = raw_data, aes(x = .data[[raw_xvar]]),
                      inherit.aes = FALSE,
                      sides = "b", alpha = 0.40, length = unit(0.03, "npc"))
  }
  p
}

# Helper: raster crop + data frame (for map panels)
rast_df <- function(rast_obj, pts_sf, buffer_m = 50000, varname = "value") {
  buf    <- st_buffer(pts_sf, dist = buffer_m)
  cr     <- terra::crop(rast_obj, vect(buf))
  df     <- as.data.frame(cr, xy = TRUE, na.rm = FALSE)
  names(df) <- c("x", "y", varname)
  df
}

dat_gam <- dat %>%
  mutate(
    survey     = factor(survey),
    year_f     = factor(year),
    preg_group = factor(preg_group,
                        levels = c("Male", "Female_nonpregnant", "Female_pregnant")),
    preg_label = factor(
      recode(as.character(preg_group), !!!preg_labels),
      levels = unname(preg_labels)
    )
  )

# ============================================================
# FIGURE 1 — SAMPLE STRUCTURE
# Panels: a) map  b) FCM by year × range  c) T3 by year × range
# ============================================================

message("\n── Figure 1: Sample structure ──")

# ---- Panel a: Map ------------------------------------------

land      <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf") |> st_transform(4326)
provinces <- rnaturalearth::ne_states(country = "canada", returnclass = "sf") |> st_transform(4326)
lakes     <- rnaturalearth::ne_download(scale = "large", type = "lakes",
                         category = "physical", returnclass = "sf") |> st_transform(4326)

# Define bb FIRST — allseason_roads clipping depends on it
pts_map <- st_as_sf(dat_map, coords = c("longitude", "latitude"), crs = 4326)
bb <- st_bbox(pts_map)
bb["xmin"] <- bb["xmin"] - 1.2; bb["xmax"] <- bb["xmax"] + 1.2
bb["ymin"] <- bb["ymin"] - 0.8; bb["ymax"] <- bb["ymax"] + 0.8

# ---- Road layers -------------------------------------------
# Both networks come from the shared, cached loader in 00b_road_layers.R,
# so the map shows exactly the road networks used to compute the
# distance-to-road covariates. Previously this block loaded roads
# independently and used a different (wrong) all-season definition than
# the analysis scripts, so the map and the models disagreed.
#
# Note the layers are already clipped to the sampling window + 2 deg by
# road_layer(), and coord_sf() clips again for display, so no geometry
# cutting is needed here.
winter_roads    <- tryCatch(road_layer("winter"),
                            error = function(e) { message("⚠ winter roads: ", e$message); NULL })
allseason_roads <- tryCatch(road_layer("allseason"),
                            error = function(e) { message("⚠ all-season roads: ", e$message); NULL })


# Fail loudly rather than quietly producing a map with no roads on it.
if (is.null(allseason_roads) || nrow(allseason_roads) == 0) {
  warning("Figure 1a will have NO all-season roads — allseason_roads is empty. ",
          "Check the messages above before using this figure.", call. = FALSE)
}
if (is.null(winter_roads) || nrow(winter_roads) == 0) {
  warning("Figure 1a will have NO winter roads — winter_roads is empty.", call. = FALSE)
}

dat_map_lbl <- dat_map %>%
  mutate(preg_label = factor(recode(as.character(preg_group), !!!preg_labels),
                              levels = unname(preg_labels)))

# ---- Road layer with its own colour legend ------------------
# Both road types are combined into ONE sf object with a road_type column so
# they can be mapped to a colour aesthetic and get a real legend. Previously
# the roads were drawn as two fixed-colour layers with an invisible
# annotate() placeholder and the note "legend handled by caption" — which
# meant the reader had no way to tell the two networks apart on the figure
# itself. ggnewscale lets the roads have their own colour scale independent
# of the range colours used by the sample points; if it isn't installed we
# fall back to a linetype legend, which uses an aesthetic the points don't.
road_colours <- c("All-season road" = "#1A9850",  # strong green
                  "Winter road"     = "#CC5500")  # orange

# ggnewscale is needed to give the roads a colour legend that is independent
# of the range colours used by the sample points. It is a small, pure-R
# package, so install it automatically if missing rather than silently
# producing a figure with no road legend.
if (!requireNamespace("ggnewscale", quietly = TRUE)) {
  message("Installing 'ggnewscale' (needed for the road legend on Figure 1a)...")
  try(install.packages("ggnewscale", repos = "https://cloud.r-project.org"),
      silent = TRUE)
}
have_ggnewscale <- requireNamespace("ggnewscale", quietly = TRUE)

roads_comb <- NULL
if (!is.null(winter_roads) && nrow(winter_roads) > 0) {
  roads_comb <- st_sf(road_type = "Winter road",
                      geometry = st_geometry(winter_roads))
}
if (!is.null(allseason_roads) && nrow(allseason_roads) > 0) {
  ar <- st_sf(road_type = "All-season road",
              geometry = st_geometry(allseason_roads))
  # All-season roads appended LAST so they draw on top: the Local/Strata
  # network inside this window is sparse and fragmented, and was previously
  # being lost under the denser winter network.
  roads_comb <- if (is.null(roads_comb)) ar else rbind(roads_comb, ar)
}

# LINE WEIGHTS — the two networks differ in density by ~400x (798 winter
# segments vs ~328,000 all-season), so they cannot share a line width.
# Drawn at a common linewidth = 0.9 the all-season network rendered as a
# solid green mass across the southern half of the map, hiding the
# province boundaries, lakes and sample points beneath it. All-season
# roads are therefore drawn hairline and semi-transparent so that
# network DENSITY reads as colour intensity (which is the informative
# signal — the network thins sharply northwards into the ranges), while
# the sparse winter network stays legible as individual routes.
# Legend keys are drawn at a uniform, readable width regardless.
# LINE WEIGHTS — the two networks differ in density by ~400x (798 winter
# segments vs ~328,000 all-season), so they cannot share a line width.
# Drawn at a common width the all-season network rendered as a solid green
# mass across the southern half of the map. All-season roads are therefore
# drawn hairline and semi-transparent so that network DENSITY reads as
# colour intensity (the informative signal — the network thins sharply
# northwards into the ranges), while the sparse winter network stays
# legible as individual routes.
#
# LEGEND BUG (fixed): these were previously two separate geom_sf layers,
# each containing only one road_type. ggplot then built the legend across
# both layers and applied the LAST layer's appearance to every key, so both
# keys rendered in the all-season green and the winter/orange key was wrong.
# Everything is now drawn from ONE geom_sf with colour, linewidth and alpha
# all mapped to road_type, so each key is derived from its own data and the
# colours are guaranteed to match the map.
road_widths <- c("All-season road" = 0.12, "Winter road" = 0.85)
road_alphas <- c("All-season road" = 0.45, "Winter road" = 0.95)

if (!is.null(roads_comb) && have_ggnewscale) {
  # All-season first so the dense southern network sits UNDER the winter
  # routes; factor level order still controls the legend order.
  roads_comb$road_type <- factor(roads_comb$road_type,
                                 levels = c("Winter road", "All-season road"))
  roads_comb <- roads_comb[order(roads_comb$road_type == "Winter road"), ]

  road_layers <- list(
    geom_sf(data = roads_comb,
            aes(colour = road_type, linewidth = road_type, alpha = road_type)),
    scale_colour_manual(values = road_colours, name = "Roads",
                        guide = guide_legend(order = 1,
                                             override.aes = list(linewidth = 1.4,
                                                                 alpha = 1))),
    scale_linewidth_manual(values = road_widths, guide = "none"),
    scale_alpha_manual(values = road_alphas, guide = "none"),
    ggnewscale::new_scale_colour()
  )
} else {
  # FALLBACK (ggnewscale unavailable): fixed-colour layers, no legend key.
  # Never map roads to `colour` here — the points' scale_colour_manual()
  # would override it and blank them out.
  road_layers <- c(
    if (!is.null(allseason_roads) && nrow(allseason_roads) > 0)
      list(geom_sf(data = allseason_roads, colour = road_colours[["All-season road"]],
                   linewidth = road_widths[["All-season road"]],
                   alpha     = road_alphas[["All-season road"]])) else list(),
    if (!is.null(winter_roads) && nrow(winter_roads) > 0)
      list(geom_sf(data = winter_roads, colour = road_colours[["Winter road"]],
                   linewidth = road_widths[["Winter road"]],
                   alpha     = road_alphas[["Winter road"]])) else list()
  )
  if (length(road_layers) > 0) {
    message("NOTE: 'ggnewscale' unavailable — Figure 1a roads are drawn but ",
            "have no legend key. Identify them in the caption: ",
            "orange = winter roads, green = all-season roads.")
  }
}

p1_a <- ggplot() +
  geom_sf(data = land,      fill = "grey96", colour = NA) +
  geom_sf(data = lakes,     fill = "#cce5ff", colour = NA) +
  geom_sf(data = provinces, fill = NA, colour = "grey55", linewidth = 0.25) +
  geom_sf(data = caribou_ranges, fill = "grey80", colour = "grey40",
          alpha = 0.30, linewidth = 0.35) +
  road_layers +
  geom_point(
    data  = dat_map_lbl,
    aes(longitude, latitude, colour = survey, shape = preg_label),
    alpha = 0.85, size = 1.6
  ) +
  scale_colour_manual(values = survey_colours, name = "Range",
                      guide = guide_legend(order = 2, ncol = 3, byrow = TRUE)) +
  scale_shape_manual(
    values = c("Male" = 15, "Non-pregnant female" = 17, "Pregnant female" = 16),
    name = "Group", guide = guide_legend(order = 3, nrow = 1)
  ) +
  annotation_scale(location = "bl", width_hint = 0.20, text_cex = 0.6) +
  coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
           ylim = c(bb["ymin"], bb["ymax"]), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_ms +
  theme(panel.background = element_rect(fill = "#e8f4f8", colour = NA),
        legend.position  = "right", legend.box = "vertical")

# ---- Panels b & c: hormone by year × range -----------------

make_panel_bc <- function(data, yvar, ylab, ybreaks, show_legend = FALSE) {
  d <- data |> filter(is.finite(.data[[yvar]]), !is.na(preg_label))
  p <- ggplot(d, aes(x = year_f, y = .data[[yvar]], colour = survey)) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.65),
                alpha = 0.35, size = 0.65) +
    geom_boxplot(aes(group = interaction(year_f, survey)),
                 position = position_dodge(width = 0.65),
                 outlier.shape = NA, fill = NA, linewidth = 0.45, width = 0.40) +
    stat_summary(aes(group = interaction(year_f, survey)),
                 fun = mean, geom = "point", shape = 23, size = 1.6,
                 fill = "white", colour = "black",
                 position = position_dodge(width = 0.65)) +
    facet_wrap(~ preg_label, ncol = 1) +
    # Must match panel a's Range guide exactly, otherwise patchwork's
    # guides = "collect" treats them as two different guides and draws both.
    scale_colour_manual(values = survey_colours, name = "Range",
                        guide = guide_legend(order = 2, ncol = 3, byrow = TRUE)) +
    scale_y_log10(breaks = ybreaks,
                  labels = scales::label_number(drop0trailing = TRUE)) +
    labs(x = "Year", y = ylab) +
    theme_ms +
    theme(legend.position = if (show_legend) "right" else "none")
  p
}

p1_b <- make_panel_bc(dat_gam, "cortisol", "FCM (ng/g, log scale)",
                       c(30, 100, 300, 1000), show_legend = TRUE)
p1_c <- make_panel_bc(dat_gam |> filter(t3 > 0), "t3", "FT3 (ng/g, log scale)",
                       c(10, 30, 100, 300), show_legend = FALSE)

# ---- Assemble Figure 1 -------------------------------------
# Layout:  [       a: map (full width)      ]
#          [  b: FCM boxes  |  c: T3 boxes  ]

# LEGEND LAYOUT — the three collected guides (Roads, Range, Group) were
# previously laid out side by side (legend.box = "horizontal"). Their
# combined width exceeded the figure width, so the outermost keys were
# clipped at both ends: "Winter road" lost its "Wi" on the left and
# "Pagwachuan" lost its final "n" on the right. Stacking the guide boxes
# vertically removes the horizontal overflow entirely, and each individual
# guide is then wrapped onto a fixed number of columns (set on each scale
# individually, NOT via a blanket guides() on the patchwork: panel a carries
# two separate colour scales because of ggnewscale, and a blanket
# guides(colour = ...) would silently retarget the Roads guide instead of
# the Range one).
fig1 <- (p1_a) /
  (p1_b | p1_c) +
  plot_layout(heights = c(5, 4), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag             = element_text(face = "bold", size = 12),
    plot.margin          = margin(4, 5, 4, 5, "pt"),
    legend.position      = "bottom",
    legend.justification = "center",
    legend.box           = "vertical",     # stack Roads / Range / Group
    legend.box.just      = "center",
    legend.margin        = margin(2, 4, 2, 4, "pt"),
    legend.spacing.y     = unit(2, "pt")
  )

ggsave(file.path(outdir_ms, "Figure1.png"), fig1,
       width = 10, height = 12.6, dpi = 300)
message("Figure 1 saved.")


# ============================================================
# FIGURE 2 — SHORT-LAG CLIMATE EFFECTS
# 2 rows × 3 cols
# Row 1 (FCM): julian day | min temp (2-day lag) | snow (2-day lag)
# Row 2 (FT3): julian day | min temp (5-day lag) | snow (5-day lag)
# (Range RE before/after climate → FigS2 in supplementary)
# FCM uses a 2-day short-lag window (Ashley et al. 2011 ACTH-challenge
# data in caribou/reindeer); FT3 keeps the 5-day window (no
# caribou-specific FT3 excretion-lag study exists). See
# 02_climate_analyses.R header for full rationale.
# ============================================================

message("\n── Figure 2: Climate effects ──")

# ---- Load and merge climate lags ---------------------------

infile_climate  <- "climatedaily_all_surveys.csv"

# NOTE: climate_group is assigned PER SAMPLE in 00_data_processing.R,
# based on the nearest station to each sample's own GPS coordinates —
# NOT by field campaign (survey_field) and NOT by boundary-assigned
# range (survey). See 00_data_processing.R for the nearest-station
# logic; see 02_climate_analyses.R for the full rationale.

clim_all <- read_csv(infile_climate, show_col_types = FALSE) |>
  janitor::clean_names() |>
  mutate(date = as.Date(local_date), local_year = as.integer(local_year),
         month = month(date))

clim_winter <- clim_all |> filter(month %in% c(1, 2, 3))

sample_survey_dates <- dat |>
  filter(!is.na(date), !is.na(climate_group)) |>
  distinct(pcid, survey, survey_field, year, date, climate_group)

clim_winter_split <- split(clim_winter, clim_winter$climate_group)

calc_lag_mean <- function(target_date, clim_df, lag_start, lag_end) {
  clim_df |>
    filter(date <  target_date - lubridate::days(lag_start - 1),
           date >= target_date - lubridate::days(lag_end)) |>
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
}

climate_lags <- sample_survey_dates |>
  rowwise() |>
  mutate(short  = list(calc_lag_mean(date, clim_winter_split[[climate_group]], 1, 5)),
         short2 = list(calc_lag_mean(date, clim_winter_split[[climate_group]], 1, 2))) |>
  ungroup() |>
  unnest(short, names_sep = "_") |>
  unnest(short2, names_sep = "_")

# Rename the nested columns to have clean names matching 02.
# short_*  (1-5 d) -> *_lag5, used for FT3 (unchanged)
# short2_* (1-2 d) -> *_lag2, used for FCM
clim_lag_clean <- climate_lags |>
  select(pcid, year, date,
         short_min_temperature, short_snow_on_ground,
         short2_min_temperature, short2_snow_on_ground) |>
  rename(min_lag5  = short_min_temperature,
         snow_lag5 = short_snow_on_ground,
         min_lag2  = short2_min_temperature,
         snow_lag2 = short2_snow_on_ground)

dat_clim <- dat |>
  left_join(clim_lag_clean, by = c("pcid", "year", "date")) |>
  mutate(survey = factor(survey), year_f = factor(year))

# ---- Climate GAMs ------------------------------------------

dat_cort_clim <- dat_clim |>
  filter(is.finite(cortisol), is.finite(min_lag2),
         is.finite(snow_lag2), !is.na(julian)) |>
  mutate(preg_group = factor(preg_group,
                              levels = c("Male","Female_nonpregnant","Female_pregnant")))

dat_t3_clim <- dat_clim |>
  filter(is.finite(t3), t3 > 0, is.finite(min_lag5),
         is.finite(snow_lag5), !is.na(julian)) |>
  mutate(preg_group = factor(preg_group,
                              levels = c("Male","Female_nonpregnant","Female_pregnant")))

m_cort_clim <- gam(
  cortisol ~ s(survey, bs = "re") + year_f + preg_group +
    s(min_lag2,  k = 6, bs = "tp") +
    s(snow_lag2, k = 4, bs = "tp") +
    s(julian, bs = "cc", k = 12),
  data = dat_cort_clim, family = Gamma(link = "log"), method = "REML"
)

# PRIMARY MODEL for FT3 (matches 02_climate_analyses.R): combined
# (snow + temperature). UPDATED — after reducing snow's k from 8 to 4,
# the AIC comparison in 02_climate_analyses.R flips: the combined model
# now has the clearly lower/better AIC (12317.83 vs 12368.42, ΔAIC ≈
# 50.6 favouring combined; previously snow-only won by only ΔAIC ≈ 2.1
# at k = 8), and snow alone is no longer significant on its own
# (p = 0.555). Panels d (julian), e (temperature), and f (snow) below
# are therefore all drawn from this combined model, matching FCM's
# convention (all three of its panels come from the single m_cort_clim
# model). The snow-only model is fit ONLY for reference/comparison — it
# is NOT the reported model and its curves are NOT used in any panel.
m_t3_clim <- gam(
  t3 ~ s(survey, bs = "re") + year_f + preg_group +
    s(min_lag5,  k = 4, bs = "tp") +
    s(snow_lag5, k = 4, bs = "tp") +
    s(julian, bs = "cc", k = 7),
  data = dat_t3_clim, family = Gamma(link = "log"), method = "REML"
)

m_t3_snowonly <- gam(
  t3 ~ s(survey, bs = "re") + year_f + preg_group +
    s(snow_lag5, k = 4, bs = "tp") +
    s(julian, bs = "cc", k = 7),
  data = dat_t3_clim, family = Gamma(link = "log"), method = "REML"
)

message("Climate GAM summaries:")
print(summary(m_cort_clim)$s.table)
message("FT3 primary (combined snow + temperature):")
print(summary(m_t3_clim)$s.table)
message("FT3 snow-only (comparison only, not reported):")
print(summary(m_t3_snowonly)$s.table)

# ---- Partial effects ----------------------------------------

sm_cort_jd   <- smooth_df(m_cort_clim, "s(julian)")
sm_cort_temp <- smooth_df(m_cort_clim, "s(min_lag2)")
sm_cort_snow <- smooth_df(m_cort_clim, "s(snow_lag2)")
sm_t3_jd     <- smooth_df(m_t3_clim, "s(julian)")
sm_t3_temp   <- smooth_df(m_t3_clim, "s(min_lag5)")
sm_t3_snow   <- smooth_df(m_t3_clim, "s(snow_lag5)")

p2_a <- partial_panel(sm_cort_jd,   "julian",    "Julian day",
                       "Partial effect on log(FCM)", dat_cort_clim, "julian")
p2_b <- partial_panel(sm_cort_temp, "min_lag2",  "2-day min. temperature (°C)",
                       "Partial effect on log(FCM)", dat_cort_clim, "min_lag2")
p2_c <- partial_panel(sm_cort_snow, "snow_lag2", "2-day snow-on-ground (cm)",
                       "Partial effect on log(FCM)", dat_cort_clim, "snow_lag2")
p2_d <- partial_panel(sm_t3_jd,    "julian",    "Julian day",
                       "Partial effect on log(FT3)", dat_t3_clim,  "julian")
p2_e <- partial_panel(sm_t3_temp,   "min_lag5",  "5-day min. temperature (°C)",
                       "Partial effect on log(FT3)", dat_t3_clim,  "min_lag5")
p2_f <- partial_panel(sm_t3_snow,   "snow_lag5", "5-day snow-on-ground (cm)",
                       "Partial effect on log(FT3)", dat_t3_clim,  "snow_lag5")

# ---- Assemble Figure 2 — 2 rows × 3 cols -------------------

fig2 <- (p2_a | p2_b | p2_c) /
        (p2_d | p2_e | p2_f) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 12),
        plot.margin = margin(3, 5, 3, 5, "pt"))

ggsave(file.path(outdir_ms, "Figure2.png"), fig2,
       width = 12, height = 7, dpi = 300)
message("Figure 2 saved.")


# ============================================================
# FIGURE 3 — LAND COVER EFFECTS (named-class interpretive models)
# 2 rows × 4 cols
# Row 1 (FCM): OpenBog | TreedPeat | ConiferousSwamp | Deciduous
# Row 2 (FT3): OpenBog | TreedPeat | ConiferousSwamp | Deciduous
# (PCA biplot → FigS6 in supplementary)
# ============================================================

message("\n── Figure 3: Land cover ──")

# ---- Extract land cover proportions from wolfstack ----------

wolfstack <- terra::rast("wolfstack.tif")

dat_sf_wgs <- dat_map |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(crs(wolfstack))

lc_vals <- terra::extract(wolfstack, vect(dat_sf_wgs))

dat_lc <- dat_map |>
  bind_cols(lc_vals |> select(-ID)) |>
  filter(!is.na(preg_group)) |>
  mutate(survey = factor(survey), year_f = factor(year),
         preg_group = factor(preg_group,
                              levels = c("Male","Female_nonpregnant","Female_pregnant")))

# Active bands: mean > 0.005 (matching 05_land_cover_analyses.R)
lc_class_cols_all <- c(
  "Water", "Wetlands", "OpenPeat", "TreedPeat",
  "Coniferous", "Sparse", "Deciduous", "NonTreed",
  "Anthropogenic", "Disturbance",
  "OpenBog", "TreedBog", "OpenFen", "TreedFen",
  "ThicketSwamp", "ConiferousSwamp", "DeciduousSwamp", "Marsh"
)
lc_bands_present <- intersect(lc_class_cols_all, names(wolfstack))

buf_lc <- st_buffer(dat_sf_wgs, dist = 50000)
lc_stack_crop <- terra::crop(wolfstack[[lc_bands_present]], vect(buf_lc))
band_means    <- terra::global(lc_stack_crop, "mean", na.rm = TRUE)[, 1]
active_bands  <- lc_bands_present[band_means > 0.005]
message("Active land cover bands: ", paste(active_bands, collapse = ", "))

# PCA on active bands (from sample-location extractions)
lc_mat <- dat_lc |> select(all_of(active_bands)) |> na.omit()
lc_complete <- dat_lc |> filter(!apply(select(dat_lc, all_of(active_bands)), 1, anyNA))
pca <- prcomp(lc_complete[, active_bands], scale. = TRUE)

message("PCA variance explained:")
print(round(summary(pca)$importance[2, 1:6] * 100, 1))

pc_scores <- as.data.frame(pca$x[, 1:4])
names(pc_scores) <- c("PC1", "PC2", "PC3", "PC4")
dat_lc_pca <- bind_cols(lc_complete, pc_scores)

# Note: PCA biplot goes to FigS6 (supplementary). Figure 3 is a clean 2×4 named-class grid.

# ---- Interpretive named-class GAMs -------------------------

named_interp <- c("OpenBog", "TreedPeat", "ConiferousSwamp", "Deciduous")

dat_lc_cort <- dat_lc_pca |> filter(is.finite(cortisol))
dat_lc_t3   <- dat_lc_pca |> filter(is.finite(t3), t3 > 0)

# sqrt-transform OpenBog and Deciduous (right-skewed proportions) so knots
# are placed at informative quantiles and the x axis is more evenly spread
dat_lc_cort <- dat_lc_cort |>
  mutate(sqrt_OpenBog   = sqrt(OpenBog),
         sqrt_Deciduous = sqrt(Deciduous))
dat_lc_t3 <- dat_lc_t3 |>
  mutate(sqrt_OpenBog   = sqrt(OpenBog),
         sqrt_Deciduous = sqrt(Deciduous))

# PROJECT DECISION: for landscape covariates (road, mine, wolf, land
# cover, population), the s(x, y)-augmented model is reported as the
# PRIMARY model, regardless of BIC (see 03_winter_road_analyses.R
# Part 2E for full rationale). Attach projected coordinates so the
# figure-generating models here match the primary models fitted in
# 05_land_cover_analyses.R.
dat_lc_cort <- dat_lc_cort |>
  bind_cols(st_coordinates(st_transform(
    st_as_sf(dat_lc_cort, coords = c("longitude", "latitude"), crs = 4326),
    crs(wolfstack)
  )) |> as.data.frame() |> rename(x_proj = X, y_proj = Y))
dat_lc_t3 <- dat_lc_t3 |>
  bind_cols(st_coordinates(st_transform(
    st_as_sf(dat_lc_t3, coords = c("longitude", "latitude"), crs = 4326),
    crs(wolfstack)
  )) |> as.data.frame() |> rename(x_proj = X, y_proj = Y))

m_cort_interp <- gam(
  cortisol ~ s(survey, bs = "re") + year_f + preg_group +
    s(sqrt_OpenBog,    k = 6, bs = "tp") +
    s(TreedPeat,       k = 6, bs = "tp") +
    s(ConiferousSwamp, k = 6, bs = "tp") +
    s(sqrt_Deciduous,  k = 6, bs = "tp") +
    s(julian, k = 8, bs = "tp") +
    s(x_proj, y_proj, bs = "tp", k = 30),
  data = dat_lc_cort, family = Gamma(link = "log"), method = "REML"
)

m_t3_interp <- gam(
  t3 ~ s(survey, bs = "re") + year_f + preg_group +
    s(sqrt_OpenBog,    k = 6, bs = "tp") +
    s(TreedPeat,       k = 6, bs = "tp") +
    s(ConiferousSwamp, k = 6, bs = "tp") +
    s(sqrt_Deciduous,  k = 6, bs = "tp") +
    s(julian, k = 8, bs = "tp") +
    s(x_proj, y_proj, bs = "tp", k = 30),
  data = dat_lc_t3, family = Gamma(link = "log"), method = "REML"
)

# Helper: build sqrt-scale x-axis breaks labelled in original proportion units
sqrt_prop_axis <- function(breaks_orig = c(0, 0.05, 0.1, 0.2, 0.3, 0.5)) {
  scale_x_continuous(
    breaks = sqrt(breaks_orig),
    labels = as.character(breaks_orig),
    expand = expansion(mult = 0.02)
  )
}

message("Land cover GAM summaries:")
print(summary(m_cort_interp)$s.table)
print(summary(m_t3_interp)$s.table)

# ---- Panels b-e (FCM) and f-i (FT3) -----------------------

make_lc_panel <- function(model, raw_data, varname,
                           xlab, ylab = "Partial effect") {
  sm <- smooth_df(model, paste0("s(", varname, ")"))
  partial_panel(sm, varname, xlab, ylab, raw_data, varname) +
    theme(plot.margin = margin(3, 6, 3, 6, "pt"))
}

cort_panels <- list(
  make_lc_panel(m_cort_interp, dat_lc_cort, "sqrt_OpenBog",
                "Open bog proportion", "Partial effect on log(FCM)") + sqrt_prop_axis(),
  make_lc_panel(m_cort_interp, dat_lc_cort, "TreedPeat",
                "Treed peat proportion", ""),
  make_lc_panel(m_cort_interp, dat_lc_cort, "ConiferousSwamp",
                "Coniferous swamp proportion", ""),
  make_lc_panel(m_cort_interp, dat_lc_cort, "sqrt_Deciduous",
                "Deciduous proportion", "") + sqrt_prop_axis()
)

t3_panels <- list(
  make_lc_panel(m_t3_interp, dat_lc_t3, "sqrt_OpenBog",
                "Open bog proportion", "Partial effect on log(FT3)") + sqrt_prop_axis(),
  make_lc_panel(m_t3_interp, dat_lc_t3, "TreedPeat",
                "Treed peat proportion", ""),
  make_lc_panel(m_t3_interp, dat_lc_t3, "ConiferousSwamp",
                "Coniferous swamp proportion", ""),
  make_lc_panel(m_t3_interp, dat_lc_t3, "sqrt_Deciduous",
                "Deciduous proportion", "") + sqrt_prop_axis()
)

row_cort <- wrap_plots(cort_panels, nrow = 1)
row_t3   <- wrap_plots(t3_panels,   nrow = 1)

fig3 <- (row_cort) /
        (row_t3) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 12),
        plot.margin = margin(3, 5, 3, 5, "pt"))

ggsave(file.path(outdir_ms, "Figure3.png"), fig3,
       width = 12, height = 7, dpi = 300)
message("Figure 3 saved.")


# ============================================================
# FIGURE 4 — COMMUNITY PROXIMITY & MINING PROXIMITY
# Panels: a) community distance map  b) distance ~ FCM  c) distance ~ FT3
#         d) mine distance map  e) log(mine dist) ~ FCM  f) log(mine dist) ~ FT3
# ============================================================

message("\n── Figure 4: Community distance + Mining ──")

# ---- Distance to nearest community ---------------------------
# Communities are derived from the population density raster (no
# independent settlements dataset was available for this region — see
# 06_population_analyses.R for full rationale): contiguous cells with
# detectable population density are clustered into patches, small
# clusters (< ~1 km^2, likely raster noise) are dropped, and each
# remaining cluster's centroid is treated as a candidate community
# location. Distance to the nearest one replaces raw population density
# as the human-settlement exposure metric, since population density was
# almost entirely zero across samples.

pop_r  <- terra::rast("population_density.tif")
pts_sf_pop <- dat_map |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(crs(pop_r))

buf_pop_fig4  <- st_buffer(pts_sf_pop, dist = 500000)
pop_crop_fig4 <- terra::crop(pop_r, vect(buf_pop_fig4))

pop_threshold_fig4     <- 0
min_cluster_km2_fig4   <- 1
cell_area_km2_fig4     <- prod(res(pop_crop_fig4)) / 1e6
min_cluster_cells_fig4 <- max(1, round(min_cluster_km2_fig4 / cell_area_km2_fig4))

pop_binary_fig4 <- pop_crop_fig4
pop_binary_fig4[pop_binary_fig4 <= pop_threshold_fig4] <- NA
pop_binary_fig4[!is.na(pop_binary_fig4)] <- 1

pop_patches_fig4 <- terra::patches(pop_binary_fig4, directions = 8)
patch_sizes_fig4 <- terra::freq(pop_patches_fig4)
small_patches_fig4 <- patch_sizes_fig4$value[patch_sizes_fig4$count < min_cluster_cells_fig4]
if (length(small_patches_fig4) > 0) pop_patches_fig4[pop_patches_fig4 %in% small_patches_fig4] <- NA

pop_polys_fig4 <- terra::as.polygons(pop_patches_fig4)
community_centers_fig4 <- terra::centroids(pop_polys_fig4) |> st_as_sf()
message("Figure 4: ", nrow(community_centers_fig4), " derived community locations.")

dist_mat_pop_fig4 <- st_distance(pts_sf_pop, community_centers_fig4)
dist_community_km_fig4 <- as.numeric(apply(dist_mat_pop_fig4, 1, min)) / 1000

dat_pop <- dat_map |>
  mutate(dist_community_km = dist_community_km_fig4) |>
  filter(is.finite(dist_community_km), !is.na(preg_group)) |>
  mutate(survey = factor(survey), year_f = factor(year),
         log_dist_community = log(dist_community_km + 0.1),
         preg_group = factor(preg_group,
                              levels = c("Male","Female_nonpregnant","Female_pregnant")))

# Map is cropped tightly to the sampling area (60 km buffer around the
# caribou points), NOT the full 500 km buffer used for the distance
# calculation itself. Plotting all derived community centroids from
# that much larger buffer previously blew out the map's automatic
# extent to cover most of the raster's bounding box instead of just the
# study area (same fix applied in 06_population_analyses.R).
buf_pop_fig4_map <- st_buffer(pts_sf_pop, dist = 60000)
pop_crop_fig4_map <- terra::crop(pop_crop_fig4, vect(buf_pop_fig4_map))
dist_rast_pop_fig4 <- terra::distance(pop_crop_fig4_map, vect(community_centers_fig4))
pop_df <- as.data.frame(dist_rast_pop_fig4 / 1000, xy = TRUE, na.rm = FALSE)
names(pop_df) <- c("x", "y", "dist_community_km")

community_centers_fig4_map <- community_centers_fig4[st_intersects(community_centers_fig4, buf_pop_fig4_map, sparse = FALSE), ]
buf_pop_fig4_bbox <- st_bbox(buf_pop_fig4_map)

pts_pop_plot <- dat_pop |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(crs(pop_r))

p4_a <- ggplot() +
  geom_raster(data = pop_df |> filter(!is.na(dist_community_km)),
              aes(x, y, fill = pmin(dist_community_km, 300))) +
  geom_sf(data = pts_pop_plot, aes(colour = survey), size = 0.9, alpha = 0.8) +
  geom_sf(data = community_centers_fig4_map, shape = 4, size = 1.3, colour = "red", stroke = 0.8) +
  scale_fill_gradient(low = "#d7301f", high = "white",
                      name = "Dist. to\ncommunity (km)", na.value = "grey90",
                      breaks = c(0, 100, 200, 300),
                      labels = c("0","100","200","≥300")) +
  scale_colour_manual(values = survey_colours, name = "Range") +
  annotation_scale(location = "bl", width_hint = 0.20, text_cex = 0.6) +
  labs(x = NULL, y = NULL) +
  coord_sf(xlim = c(buf_pop_fig4_bbox["xmin"], buf_pop_fig4_bbox["xmax"]),
           ylim = c(buf_pop_fig4_bbox["ymin"], buf_pop_fig4_bbox["ymax"]), expand = FALSE) +
  theme_ms +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_rect(fill = "#e8f4f8", colour = NA))

# Community-distance GAMs
dat_pop_cort <- dat_pop |> filter(is.finite(cortisol))
dat_pop_t3   <- dat_pop |> filter(is.finite(t3), t3 > 0)

# PROJECT DECISION: for landscape covariates, the s(x, y)-augmented
# model is reported as the PRIMARY model, regardless of BIC (see
# 03_winter_road_analyses.R Part 2E). Attach projected coordinates so
# this figure model matches the primary model in 06_population_analyses.R.
dat_pop_cort <- dat_pop_cort |>
  bind_cols(st_coordinates(st_transform(
    st_as_sf(dat_pop_cort, coords = c("longitude", "latitude"), crs = 4326),
    crs(pop_r)
  )) |> as.data.frame() |> rename(x_proj = X, y_proj = Y))
dat_pop_t3 <- dat_pop_t3 |>
  bind_cols(st_coordinates(st_transform(
    st_as_sf(dat_pop_t3, coords = c("longitude", "latitude"), crs = 4326),
    crs(pop_r)
  )) |> as.data.frame() |> rename(x_proj = X, y_proj = Y))

m_pop_cort <- gam(cortisol ~ s(survey, bs = "re") + year_f + preg_group +
                    s(log_dist_community, k = 6, bs = "tp") + s(julian, k = 8, bs = "tp") +
                    s(x_proj, y_proj, bs = "tp", k = 30),
                  data = dat_pop_cort, family = Gamma(link = "log"), method = "REML")
m_pop_t3   <- gam(t3 ~ s(survey, bs = "re") + year_f + preg_group +
                    s(log_dist_community, k = 6, bs = "tp") + s(julian, k = 8, bs = "tp") +
                    s(x_proj, y_proj, bs = "tp", k = 30),
                  data = dat_pop_t3, family = Gamma(link = "log"), method = "REML")

sm_pop_cort <- smooth_df(m_pop_cort, "s(log_dist_community)")
sm_pop_t3   <- smooth_df(m_pop_t3,   "s(log_dist_community)")

p4_b <- partial_panel(sm_pop_cort, "log_dist_community", "log(distance to nearest community, km)",
                       "Partial effect on log(FCM)", dat_pop_cort, "log_dist_community")
p4_c <- partial_panel(sm_pop_t3,   "log_dist_community", "log(distance to nearest community, km)",
                       "Partial effect on log(FT3)", dat_pop_t3,  "log_dist_community")

# ---- Mine distance -----------------------------------------

mines_r <- terra::rast("mines.tif")
pts_sf_mn <- dat_map |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(crs(mines_r))

buf_mn     <- st_buffer(pts_sf_mn, dist = 500000)
mines_crop <- terra::crop(mines_r, vect(buf_mn))
mine_binary <- mines_crop; mine_binary[mine_binary == 0] <- NA

message("Computing mine distance raster (~1 min)...")
dist_rast <- terra::distance(mine_binary)

mines_vals <- terra::extract(dist_rast, vect(pts_sf_mn))[, 2] / 1000

dat_mines <- dat_map |>
  mutate(dist_mine_km = mines_vals) |>
  filter(is.finite(dist_mine_km), !is.na(preg_group)) |>
  mutate(survey = factor(survey), year_f = factor(year),
         log_dist_mine = log(dist_mine_km + 0.1),
         preg_group = factor(preg_group,
                              levels = c("Male","Female_nonpregnant","Female_pregnant")))

# Mine distance map
mine_df <- rast_df(dist_rast, pts_sf_mn, buffer_m = 100000, varname = "dist_mine_m")
mine_df$dist_mine_km <- mine_df$dist_mine_m / 1000

pts_mine_sf <- dat_mines |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(crs(mines_r))

p4_d <- ggplot() +
  geom_raster(data = mine_df |> filter(!is.na(dist_mine_km)),
              aes(x, y, fill = pmin(dist_mine_km, 300))) +
  geom_sf(data = pts_mine_sf, aes(colour = survey), size = 0.9, alpha = 0.8) +
  scale_fill_gradient(low = "#d7301f", high = "white",
                      name = "Dist. to mine (km)",
                      na.value = "grey90",
                      breaks = c(0, 100, 200, 300),
                      labels = c("0","100","200","≥300")) +
  scale_colour_manual(values = survey_colours, name = "Range") +
  annotation_scale(location = "bl", width_hint = 0.20, text_cex = 0.6) +
  labs(x = NULL, y = NULL) +
  coord_sf(expand = FALSE) +
  theme_ms +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_rect(fill = "#e8f4f8", colour = NA))

dat_mines_cort <- dat_mines |> filter(is.finite(cortisol))
dat_mines_t3   <- dat_mines |> filter(is.finite(t3), t3 > 0)

# PROJECT DECISION: for landscape covariates, the s(x, y)-augmented
# model is reported as the PRIMARY model, regardless of BIC (see
# 03_winter_road_analyses.R Part 2E). Attach projected coordinates so
# this figure model matches the primary model in 07_mines_analyses.R.
dat_mines_cort <- dat_mines_cort |>
  bind_cols(st_coordinates(st_transform(
    st_as_sf(dat_mines_cort, coords = c("longitude", "latitude"), crs = 4326),
    crs(mines_r)
  )) |> as.data.frame() |> rename(x_proj = X, y_proj = Y))
dat_mines_t3 <- dat_mines_t3 |>
  bind_cols(st_coordinates(st_transform(
    st_as_sf(dat_mines_t3, coords = c("longitude", "latitude"), crs = 4326),
    crs(mines_r)
  )) |> as.data.frame() |> rename(x_proj = X, y_proj = Y))

m_mine_cort <- gam(cortisol ~ s(survey, bs = "re") + year_f + preg_group +
                     s(log_dist_mine, k = 6, bs = "tp") + s(julian, k = 8, bs = "tp") +
                     s(x_proj, y_proj, bs = "tp", k = 30),
                   data = dat_mines_cort, family = Gamma(link = "log"), method = "REML")
m_mine_t3   <- gam(t3 ~ s(survey, bs = "re") + year_f + preg_group +
                     s(log_dist_mine, k = 6, bs = "tp") + s(julian, k = 8, bs = "tp") +
                     s(x_proj, y_proj, bs = "tp", k = 30),
                   data = dat_mines_t3, family = Gamma(link = "log"), method = "REML")

sm_mine_cort <- smooth_df(m_mine_cort, "s(log_dist_mine)")
sm_mine_t3   <- smooth_df(m_mine_t3,   "s(log_dist_mine)")

p4_e <- partial_panel(sm_mine_cort, "log_dist_mine", "log(distance to nearest mine, km)",
                       "Partial effect on log(FCM)", dat_mines_cort, "log_dist_mine")
p4_f <- partial_panel(sm_mine_t3,   "log_dist_mine", "log(distance to nearest mine, km)",
                       "Partial effect on log(FT3)", dat_mines_t3,  "log_dist_mine")

# ---- Assemble Figure 4 -------------------------------------

fig4 <- (p4_a | p4_b | p4_c) /
        (p4_d | p4_e | p4_f) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 12),
        plot.margin = margin(3, 5, 3, 5, "pt"))

ggsave(file.path(outdir_ms, "Figure4.png"), fig4,
       width = 12, height = 8, dpi = 300)
message("Figure 4 saved.")


# ============================================================
# FIGURE 5 — WOLF HABITAT SUITABILITY
# Panels: a) wolf map  b) wolf ~ FCM  c) wolf ~ FT3
# ============================================================

message("\n── Figure 5: Wolf suitability ──")

wolf_r <- terra::rast("wolfmean_t.tif")

pts_sf_wf <- dat_map |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(crs(wolf_r))

wolf_vals <- terra::extract(wolf_r, vect(pts_sf_wf))[, 2]

dat_wolf <- dat_map |>
  mutate(wolf_suit = wolf_vals) |>
  filter(is.finite(wolf_suit), !is.na(preg_group)) |>
  mutate(survey = factor(survey), year_f = factor(year),
         preg_group = factor(preg_group,
                              levels = c("Male","Female_nonpregnant","Female_pregnant")),
         log_wolf_suit = log(wolf_suit + 0.001))   # log-transform: suitability is right-skewed

wolf_df <- rast_df(wolf_r, pts_sf_wf, buffer_m = 60000, varname = "wolf_suit")
pts_wolf_sf <- dat_wolf |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(crs(wolf_r))

p5_a <- ggplot() +
  geom_raster(data = wolf_df |> filter(!is.na(wolf_suit)),
              aes(x, y, fill = wolf_suit)) +
  geom_sf(data = pts_wolf_sf,
          aes(colour = survey), size = 0.9, alpha = 0.8) +
  scale_fill_gradient(low = "white", high = "#1a6930",
                      name = "Wolf suitability", na.value = "grey90") +
  scale_colour_manual(values = survey_colours, name = "Range") +
  annotation_scale(location = "bl", width_hint = 0.20, text_cex = 0.6) +
  labs(x = NULL, y = NULL) +
  coord_sf(expand = FALSE) +
  theme_ms +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.background = element_rect(fill = "#e8f4f8", colour = NA))

dat_wolf_cort <- dat_wolf |> filter(is.finite(cortisol))
dat_wolf_t3   <- dat_wolf |> filter(is.finite(t3), t3 > 0)

# PROJECT DECISION: for landscape covariates, the s(x, y)-augmented
# model is reported as the PRIMARY model, regardless of BIC (see
# 03_winter_road_analyses.R Part 2E). Attach projected coordinates so
# this figure model matches the primary model in 08_predation_analyses.R.
dat_wolf_cort <- dat_wolf_cort |>
  bind_cols(st_coordinates(st_transform(
    st_as_sf(dat_wolf_cort, coords = c("longitude", "latitude"), crs = 4326),
    crs(wolf_r)
  )) |> as.data.frame() |> rename(x_proj = X, y_proj = Y))
dat_wolf_t3 <- dat_wolf_t3 |>
  bind_cols(st_coordinates(st_transform(
    st_as_sf(dat_wolf_t3, coords = c("longitude", "latitude"), crs = 4326),
    crs(wolf_r)
  )) |> as.data.frame() |> rename(x_proj = X, y_proj = Y))

m_wolf_cort <- gam(cortisol ~ s(survey, bs = "re") + year_f + preg_group +
                     s(log_wolf_suit, k = 6, bs = "tp") + s(julian, k = 8, bs = "tp") +
                     s(x_proj, y_proj, bs = "tp", k = 30),
                   data = dat_wolf_cort, family = Gamma(link = "log"), method = "REML")
m_wolf_t3   <- gam(t3 ~ s(survey, bs = "re") + year_f + preg_group +
                     s(log_wolf_suit, k = 6, bs = "tp") + s(julian, k = 8, bs = "tp") +
                     s(x_proj, y_proj, bs = "tp", k = 30),
                   data = dat_wolf_t3, family = Gamma(link = "log"), method = "REML")

sm_wolf_cort <- smooth_df(m_wolf_cort, "s(log_wolf_suit)")
sm_wolf_t3   <- smooth_df(m_wolf_t3,   "s(log_wolf_suit)")

# X-axis: log scale internally, labelled in original suitability units
wolf_log_axis <- scale_x_continuous(
  breaks = log(c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 1) + 0.001),
  labels = c("0", "0.01", "0.05", "0.1", "0.25", "0.5", "1"),
  name   = "Wolf habitat suitability"
)

p5_b <- partial_panel(sm_wolf_cort, "log_wolf_suit", "Wolf habitat suitability (log scale)",
                       "Partial effect on log(FCM)", dat_wolf_cort, "log_wolf_suit") +
        wolf_log_axis
p5_c <- partial_panel(sm_wolf_t3,   "log_wolf_suit", "Wolf habitat suitability (log scale)",
                       "Partial effect on log(FT3)", dat_wolf_t3,  "log_wolf_suit") +
        wolf_log_axis

fig5 <- (p5_a) /
        (p5_b | p5_c) +
  plot_layout(heights = c(1.5, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 12),
        plot.margin = margin(3, 5, 3, 5, "pt"))

ggsave(file.path(outdir_ms, "Figure5.png"), fig5,
       width = 9, height = 9, dpi = 300)
message("Figure 5 saved.")

message("\n✓ All main manuscript figures saved to: ", outdir_ms)
