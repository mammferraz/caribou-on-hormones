# ============================================================
# 08_predation_analyses.R
#
# PART 8 — PREDATION RISK (WOLF HABITAT SUITABILITY)
# Purpose:
#   • Map predicted wolf habitat suitability
#   • Extract values at caribou sample locations
#   • Test effects on cortisol and T3 by year and survey
# ============================================================

source("00_data_processing.R")

outdir <- "caribou_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# I/O
# ============================================================

infile_wolf <- "wolfmean_t.tif"

# ============================================================
# LOAD RASTER & EXTRACT AT SAMPLE LOCATIONS
# ============================================================

wolf_r <- terra::rast(infile_wolf)

pts_sf_wf <- dat_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(wolf_r))

wolf_vals <- terra::extract(wolf_r, vect(pts_sf_wf))[, 2]

dat_wolf <- dat_map %>%
  mutate(wolf_suit = wolf_vals) %>%
  filter(is.finite(wolf_suit), !is.na(preg_group)) %>%
  mutate(survey = factor(survey), year_f = factor(year),
         log_wolf_suit = log(wolf_suit + 0.001))   # log-transform: suitability is right-skewed

cat("Wolf habitat suitability summary:\n")
print(summary(dat_wolf$wolf_suit))

# ============================================================
# 1) MAP: WOLF SUITABILITY + CARIBOU LOCATIONS
# ============================================================

pts_sf_proj <- dat_wolf %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(wolf_r))

buf       <- st_buffer(pts_sf_proj, dist = 50000)
wolf_crop <- terra::crop(wolf_r, vect(buf))
wolf_df   <- as.data.frame(wolf_crop, xy = TRUE, na.rm = FALSE)
names(wolf_df) <- c("x", "y", "wolf_suit")

p_wolf_map <- ggplot() +
  geom_raster(data = wolf_df, aes(x, y, fill = wolf_suit)) +
  geom_sf(data = pts_sf_proj,
          aes(colour = survey, shape = preg_group),
          size = 1.8, alpha = 0.95) +
  scale_fill_gradient(
    name     = "Wolf habitat\nsuitability",
    low      = "white",
    high     = "#0072B2",   # Okabe-Ito blue
    limits   = c(0, NA),
    na.value = "white"
  ) +
  scale_colour_manual(values = survey_colours, na.value = "grey50") +
  scale_shape_manual(values = preg_shapes) +
  coord_sf(expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.35) +
  labs(title  = "Wolf habitat suitability and caribou sampling locations",
       x = "Easting", y = "Northing",
       colour = "Survey/Range", shape = "Group") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(outdir, "Fig_WolfSuitability_Map.png"),
       p_wolf_map, width = 8, height = 6, dpi = 300)


# ============================================================
# 2) WOLF SUITABILITY BY YEAR × SURVEY
# ============================================================

p_wolf_year <- ggplot(
  dat_wolf,
  aes(x = factor(year), y = wolf_suit, colour = survey)
) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
              alpha = 0.5, size = 1.8) +
  geom_boxplot(aes(group = interaction(year, survey)),
               position = position_dodge(0.6),
               outlier.shape = NA, fill = NA, linewidth = 0.7, width = 0.4) +
  scale_colour_manual(values = survey_colours, na.value = "grey50") +
  labs(x = "Year", y = "Wolf habitat suitability",
       colour = "Survey/Range",
       title  = "Wolf habitat suitability at caribou locations by year") +
  theme_pub

ggsave(file.path(outdir, "Fig_WolfSuitability_by_Year_Survey.png"),
       p_wolf_year, width = 8, height = 5, dpi = 300)

write_csv(
  dat_wolf %>%
    group_by(year, survey) %>%
    summarise(n = n(), wolf_mean = mean(wolf_suit), wolf_sd = sd(wolf_suit),
              .groups = "drop"),
  file.path(outdir, "Summary_WolfSuitability_by_Year_Survey.csv")
)


# ============================================================
# 3) GAMs — WOLF SUITABILITY EFFECTS
# ============================================================

# --- 3A) Cortisol ---

dat_wolf_cort <- dat_wolf %>%
  filter(is.finite(cortisol), !is.na(julian))

m_cort_wolf <- gam(
  cortisol ~
    s(survey,         bs = "re") +
    year_f +
    preg_group +
    s(log_wolf_suit,  k = 6, bs = "tp") +   # log-transformed; k reduced from 8
    s(julian,         k = 8, bs = "tp"),    # k = 8 to match 09_manuscript_figures.R
  data   = dat_wolf_cort,
  family = Gamma(link = "log"),
  method = "REML"
)

summary(m_cort_wolf)
capture.output(summary(m_cort_wolf), file = file.path(outdir, "GAM_summary_08_cortisol_wolf.txt"))
gam.check(m_cort_wolf)

p_cort_wolf <- draw(m_cort_wolf, select = "s(log_wolf_suit)") +
  ggtitle("Cortisol — partial effect of wolf habitat suitability (log scale)")

ggsave(file.path(outdir, "Fig_Cortisol_WolfSuitability_PartialEffect.png"),
       p_cort_wolf, width = 6, height = 5, dpi = 300)

png(file.path(outdir, "Fig_GAM_Cortisol_WolfSuitability.png"), 2000, 1400, res = 300)
plot(m_cort_wolf, pages = 1)
dev.off()

cort_wolf_year_contrasts <- emmeans(m_cort_wolf, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(cort_wolf_year_contrasts), file = file.path(outdir, "GAM_contrasts_08_cortisol_wolf_year.txt"))

# --- 3B) T3 ---

dat_wolf_t3 <- dat_wolf %>%
  filter(is.finite(t3), t3 > 0, !is.na(julian))

m_t3_wolf <- gam(
  t3 ~
    s(survey,         bs = "re") +
    year_f +
    preg_group +
    s(log_wolf_suit,  k = 6, bs = "tp") +   # log-transformed; k reduced from 8
    s(julian,         k = 8, bs = "tp"),    # k = 8 to match 09_manuscript_figures.R
  data   = dat_wolf_t3,
  family = Gamma(link = "log"),
  method = "REML"
)

summary(m_t3_wolf)
capture.output(summary(m_t3_wolf), file = file.path(outdir, "GAM_summary_08_t3_wolf.txt"))

p_t3_wolf <- draw(m_t3_wolf, select = "s(log_wolf_suit)") +
  ggtitle("T3 — partial effect of wolf habitat suitability (log scale)")

ggsave(file.path(outdir, "Fig_T3_WolfSuitability_PartialEffect.png"),
       p_t3_wolf, width = 6, height = 5, dpi = 300)

png(file.path(outdir, "Fig_GAM_T3_WolfSuitability.png"), 2000, 1400, res = 300)
plot(m_t3_wolf, pages = 1)
dev.off()

t3_wolf_year_contrasts <- emmeans(m_t3_wolf, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(t3_wolf_year_contrasts), file = file.path(outdir, "GAM_contrasts_08_t3_wolf_year.txt"))

# ============================================================
# 4) SPATIAL AUTOCORRELATION DIAGNOSTICS (variogram on residuals)
# Purpose: same diagnostic applied to the road-distance models in
#   03_winter_road_analyses.R, now extended to wolf habitat suitability —
#   check whether residual spatial autocorrelation could be inflating
#   the apparent significance of the log_wolf_suit smooth.
#
# Some samples share IDENTICAL GPS coordinates (repeat visits to the
# same collection site). A raw per-sample variogram would treat those
# as separate zero-distance observations and read as pure nugget noise
# rather than genuine small-scale variation, so residuals are averaged
# to one value per unique sampling location before the variogram is
# computed. If the variogram shows real spatial structure (sill clearly
# above the nugget, with a non-trivial range), the wolf models are
# refit with an added s(x, y) smooth on the projected coordinates, and
# the log_wolf_suit effect is re-checked with that term included.
# ============================================================

if (!requireNamespace("gstat", quietly = TRUE)) install.packages("gstat", repos = "https://cloud.r-project.org")
library(gstat)

sf_cort_wolf <- dat_wolf_cort %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(wolf_r))
sf_t3_wolf <- dat_wolf_t3 %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(wolf_r))

check_spatial_autocorr_wolf <- function(resid_vec, data_sf, label) {
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

  png(file.path(outdir, paste0("Variogram_", label, "_WolfModel.png")), 1600, 1200, res = 200)
  print(plot(vgram, main = paste0(label, " — variogram of wolf-model deviance residuals"),
             xlab = "Distance (m)", ylab = "Semivariance"))
  dev.off()

  vgram_fit <- tryCatch(gstat::fit.variogram(vgram, gstat::vgm("Sph")), error = function(e) NULL)
  if (!is.null(vgram_fit)) {
    cat(sprintf("%s variogram model (nugget / partial sill / range):\n", label))
    print(vgram_fit)
    capture.output(print(vgram_fit),
                   file = file.path(outdir, paste0("Variogram_fit_", label, "_WolfModel.txt")))
  } else {
    message(label, ": variogram model fit did not converge — inspect the empirical variogram plot directly.")
  }
  list(vgram = vgram, vgram_fit = vgram_fit)
}

vres_cort_wolf <- check_spatial_autocorr_wolf(residuals(m_cort_wolf, type = "deviance"), sf_cort_wolf, "FCM")
vres_t3_wolf   <- check_spatial_autocorr_wolf(residuals(m_t3_wolf,   type = "deviance"), sf_t3_wolf,   "FT3")

# ---- Refit with s(x, y) added, and re-check whether log_wolf_suit survives ----

dat_cort_wolf_xy <- dat_wolf_cort %>%
  bind_cols(st_coordinates(sf_cort_wolf) %>% as.data.frame() %>% rename(x_proj = X, y_proj = Y))
dat_t3_wolf_xy <- dat_wolf_t3 %>%
  bind_cols(st_coordinates(sf_t3_wolf) %>% as.data.frame() %>% rename(x_proj = X, y_proj = Y))

m_cort_wolf_xy <- gam(
  cortisol ~
    s(survey,        bs = "re") +
    year_f +
    preg_group +
    s(log_wolf_suit, k = 6, bs = "tp") +
    s(julian,        k = 8, bs = "tp") +
    s(x_proj, y_proj, bs = "tp", k = 30),
  data   = dat_cort_wolf_xy,
  family = Gamma(link = "log"),
  method = "REML"
)

plot(m_cort_wolf_xy, pages = 1, scheme = 2)

m_t3_wolf_xy <- gam(
  t3 ~
    s(survey,        bs = "re") +
    year_f +
    preg_group +
    s(log_wolf_suit, k = 6, bs = "tp") +
    s(julian,        k = 8, bs = "tp") +
    s(x_proj, y_proj, bs = "tp", k = 30),
  data   = dat_t3_wolf_xy,
  family = Gamma(link = "log"),
  method = "REML"
)

plot(m_t3_wolf_xy, pages = 1, scheme = 2)

cat("\nWolf models with s(x, y) added — AIC/BIC comparison against the originally reported (no-spatial-term) models:\n")
aic_cort_wolf_xy <- AIC(m_cort_wolf, m_cort_wolf_xy)
aic_t3_wolf_xy   <- AIC(m_t3_wolf,   m_t3_wolf_xy)
bic_cort_wolf_xy <- BIC(m_cort_wolf, m_cort_wolf_xy)
bic_t3_wolf_xy   <- BIC(m_t3_wolf,   m_t3_wolf_xy)
print(aic_cort_wolf_xy); print(bic_cort_wolf_xy)
print(aic_t3_wolf_xy); print(bic_t3_wolf_xy)
capture.output(print(aic_cort_wolf_xy), print(bic_cort_wolf_xy), print(aic_t3_wolf_xy), print(bic_t3_wolf_xy),
               file = file.path(outdir, "GAM_AIC_08_spatial_xy_comparison.txt"))

cat("\nDoes log_wolf_suit remain significant once s(x, y) is included?\n")
print(summary(m_cort_wolf_xy)$s.table)
print(summary(m_t3_wolf_xy)$s.table)
capture.output(summary(m_cort_wolf_xy), file = file.path(outdir, "GAM_summary_08_cortisol_wolf_xy.txt"))
capture.output(summary(m_t3_wolf_xy),   file = file.path(outdir, "GAM_summary_08_t3_wolf_xy.txt"))

# Re-check the variogram on the xy-augmented models' residuals — if
# s(x, y) has done its job, this should look flatter than the original.
vres_cort_wolf_after <- check_spatial_autocorr_wolf(residuals(m_cort_wolf_xy, type = "deviance"), sf_cort_wolf, "FCM_afterXY")
vres_t3_wolf_after   <- check_spatial_autocorr_wolf(residuals(m_t3_wolf_xy,   type = "deviance"), sf_t3_wolf,   "FT3_afterXY")

message("Part 8 spatial autocorrelation diagnostics complete — outputs in: ", outdir)

# ============================================================
# 5) PRIMARY MODEL DECISION — spatially-adjusted model reported
# Project decision: for landscape covariates (road, mine, wolf, land
# cover, population), the s(x, y)-augmented model is reported as the
# PRIMARY model, regardless of BIC (see 03_winter_road_analyses.R
# Part 2E for full rationale). m_cort_wolf_xy / m_t3_wolf_xy are
# therefore the primary reported wolf models. The outputs below
# OVERWRITE the original (non-spatial) GAM_summary/Fig/contrasts files
# from Part 3 so every file on disk reflects the model actually
# reported in the manuscript.
# ============================================================

summary(m_cort_wolf_xy)
capture.output(summary(m_cort_wolf_xy), file = file.path(outdir, "GAM_summary_08_cortisol_wolf.txt"))

p_cort_wolf <- draw(m_cort_wolf_xy, select = "s(log_wolf_suit)") +
  ggtitle("Cortisol — partial effect of wolf habitat suitability (log scale)")
ggsave(file.path(outdir, "Fig_Cortisol_WolfSuitability_PartialEffect.png"),
       p_cort_wolf, width = 6, height = 5, dpi = 300)

png(file.path(outdir, "Fig_GAM_Cortisol_WolfSuitability.png"), 2000, 1400, res = 300)
plot(m_cort_wolf_xy, pages = 1)
dev.off()

cort_wolf_year_contrasts <- emmeans(m_cort_wolf_xy, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(cort_wolf_year_contrasts), file = file.path(outdir, "GAM_contrasts_08_cortisol_wolf_year.txt"))

summary(m_t3_wolf_xy)
capture.output(summary(m_t3_wolf_xy), file = file.path(outdir, "GAM_summary_08_t3_wolf.txt"))

p_t3_wolf <- draw(m_t3_wolf_xy, select = "s(log_wolf_suit)") +
  ggtitle("T3 — partial effect of wolf habitat suitability (log scale)")
ggsave(file.path(outdir, "Fig_T3_WolfSuitability_PartialEffect.png"),
       p_t3_wolf, width = 6, height = 5, dpi = 300)

png(file.path(outdir, "Fig_GAM_T3_WolfSuitability.png"), 2000, 1400, res = 300)
plot(m_t3_wolf_xy, pages = 1)
dev.off()

t3_wolf_year_contrasts <- emmeans(m_t3_wolf_xy, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(t3_wolf_year_contrasts), file = file.path(outdir, "GAM_contrasts_08_t3_wolf_year.txt"))

message("Part 5 (primary model = spatially-adjusted) complete — GAM_summary_08_cortisol_wolf.txt / GAM_summary_08_t3_wolf.txt now reflect the s(x,y)-augmented models.")

message("Part 8 complete — outputs in: ", outdir)
