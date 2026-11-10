# ============================================================
# 05_land_cover_analyses.R
#
# PART 5 — LAND COVER CONTEXT (WOLFSTACK PCA-BASED GAMs)
# Purpose:
#   • Map Wolfstack 1-km moving-window land cover proportions
#   • Extract at caribou sample locations
#   • Reduce to PCA axes
#   • Test whether land cover gradients explain cortisol/T3
#     controlling for pregnancy group, year, julian day, survey
# ============================================================

source("00_data_processing.R")

outdir <- "caribou_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# I/O
# ============================================================

infile_wolfstack <- "wolfstack.tif"

# ============================================================
# LOAD RASTER & EXTRACT AT SAMPLE LOCATIONS
# ============================================================

wolfstack <- terra::rast(infile_wolfstack)

dat_sf_wgs <- dat_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(wolfstack))

lc_vals <- terra::extract(wolfstack, vect(dat_sf_wgs))
lc_cols <- names(lc_vals)[-1]  # all extracted bands (drop ID)

dat_lc <- dat_map %>%
  bind_cols(lc_vals %>% select(-ID)) %>%
  filter(is.finite(cortisol) | is.finite(t3), !is.na(preg_group)) %>%
  mutate(survey = factor(survey), year_f = factor(year))

# ============================================================
# 1) LAND COVER MAP — all classes in faceted multi-panel
# ============================================================

# All land cover class bands we expect from wolfstack
lc_class_cols_all <- c(
  "Water", "Wetlands", "OpenPeat", "TreedPeat",
  "Coniferous", "Sparse", "Deciduous", "NonTreed",
  "Anthropogenic", "Disturbance",
  "OpenBog", "TreedBog", "OpenFen", "TreedFen",
  "ThicketSwamp", "ConiferousSwamp", "DeciduousSwamp", "Marsh"
)
lc_bands_present <- intersect(lc_class_cols_all, names(wolfstack))

# Project sample points and define bounding box
pts_sf_lc <- dat_map %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(wolfstack))

buf_lc <- st_buffer(pts_sf_lc, dist = 50000)

# Crop all bands and convert to a single long-format data frame
lc_stack_crop <- terra::crop(wolfstack[[lc_bands_present]], vect(buf_lc))

# Keep only bands with at least some non-trivial signal (mean > 0.005)
band_means <- global(lc_stack_crop, "mean", na.rm = TRUE)[, 1]
active_bands <- lc_bands_present[band_means > 0.005]
lc_stack_crop <- lc_stack_crop[[active_bands]]
message("Land cover classes retained: ", paste(active_bands, collapse = ", "))

# Convert raster to long data frame — one row per pixel per class
lc_df_list <- lapply(active_bands, function(b) {
  df <- as.data.frame(lc_stack_crop[[b]], xy = TRUE, na.rm = FALSE)
  names(df) <- c("x", "y", "proportion")
  df$class   <- b
  df
})
lc_df_long <- do.call(rbind, lc_df_list)
lc_df_long$class <- factor(lc_df_long$class, levels = active_bands)

# Sample points data frame for overlaying
pts_df <- as.data.frame(st_coordinates(pts_sf_lc))
names(pts_df) <- c("x", "y")
pts_df$survey     <- dat_map$survey[!is.na(dat_map$latitude) & !is.na(dat_map$longitude)]
pts_df$preg_group <- dat_map$preg_group[!is.na(dat_map$latitude) & !is.na(dat_map$longitude)]

# Build faceted figure — one panel per land cover class
n_cols  <- 3
n_bands <- length(active_bands)
fig_h   <- ceiling(n_bands / n_cols) * 3.5 + 1   # dynamic height

p_lc_map <- ggplot() +
  geom_raster(data = lc_df_long,
              aes(x, y, fill = proportion)) +
  geom_point(data = pts_df,
             aes(x, y, colour = survey, shape = preg_group),
             size = 0.7, alpha = 0.8) +
  facet_wrap(~ class, ncol = n_cols) +
  scale_fill_gradientn(
    colours  = c("#f7f7f7", "#a8ddb5", "#41ab5d", "#238b45", "#005a32"),
    na.value = "white",
    name     = "Proportion",
    limits   = c(0, 1),
    breaks   = c(0, 0.25, 0.5, 0.75, 1),
    guide    = guide_colourbar(barwidth = 0.8, barheight = 5)
  ) +
  scale_colour_manual(values = survey_colours, na.value = "grey50",
                      name = "Survey") +
  scale_shape_manual(values = preg_shapes, name = "Group") +
  coord_equal(expand = FALSE) +
  labs(
    title = "Land cover class proportions (1-km moving window) and caribou sampling locations",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title     = element_text(face = "bold", size = 10),
    axis.text      = element_blank(),
    axis.ticks     = element_blank(),
    panel.spacing  = unit(0.3, "lines"),
    strip.text     = element_text(face = "bold", size = 8),
    legend.position = "right"
  )

ggsave(file.path(outdir, "Fig_LandCover_Map.png"),
       p_lc_map, width = 11, height = fig_h, dpi = 300, limitsize = FALSE)

message("Land cover map saved (", n_bands, " classes).")

# ============================================================
# 2) PCA ON LAND COVER CLASS PROPORTIONS ONLY
# ============================================================

# Land cover class columns only (exclude climate, topography,
# vegetation indices, species codes, wildlife density, distances)
lc_class_cols <- c(
  "Water", "Wetlands", "OpenPeat", "TreedPeat",
  "Coniferous", "Sparse", "Deciduous", "NonTreed",
  "Anthropogenic", "Disturbance",
  "OpenBog", "TreedBog", "OpenFen", "TreedFen",
  "ThicketSwamp", "ConiferousSwamp", "DeciduousSwamp", "Marsh"
)
# Keep only those present in the extracted data
lc_class_cols <- intersect(lc_class_cols, names(dat_lc))

lc_matrix <- dat_lc %>%
  select(all_of(lc_class_cols)) %>%
  select(where(~ var(., na.rm = TRUE) > 0)) %>%
  replace(is.na(.), 0)

pca_lc <- prcomp(lc_matrix, scale. = TRUE, center = TRUE)
cat("PCA variance explained (first 5 PCs) — land cover classes only:\n")
print(summary(pca_lc)$importance[, 1:5])

# Top loadings per PC for interpretation
rot <- as.data.frame(pca_lc$rotation[, 1:4])
cat("\nTop 5 loadings per PC:\n")
for (pc in paste0("PC", 1:4)) {
  vals <- setNames(rot[[pc]], rownames(rot))
  top_idx <- order(abs(vals), decreasing = TRUE)[1:5]
  cat("\n---", pc, "---\n")
  print(round(vals[top_idx], 3))
}

dat_lc <- dat_lc %>%
  bind_cols(
    as_tibble(pca_lc$x[, 1:4]) %>%
      rename(LC_PC1 = PC1, LC_PC2 = PC2, LC_PC3 = PC3, LC_PC4 = PC4)
  )

# PCA biplot (first two axes)
# Scale arrows so the longest reaches ~75% of the score range
pca_df  <- as.data.frame(pca_lc$x[, 1:2]) %>%
  bind_cols(dat_lc %>% select(survey, year))
load_df <- as.data.frame(pca_lc$rotation[, 1:2]) %>%
  rownames_to_column("land_class")

arrow_scale <- 0.75 * max(abs(pca_lc$x[, 1:2])) /
               max(sqrt(load_df$PC1^2 + load_df$PC2^2))

# Keep only top 6 loadings by magnitude to avoid label clutter
load_df_top <- load_df %>%
  mutate(mag = sqrt(PC1^2 + PC2^2)) %>%
  slice_max(mag, n = 6)

p_biplot <- ggplot(pca_df, aes(PC1, PC2, colour = survey)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_segment(
    data = load_df_top,
    aes(x = 0, y = 0,
        xend = PC1 * arrow_scale,
        yend = PC2 * arrow_scale),
    arrow = arrow(length = unit(0.25, "cm")),
    colour = "grey20", linewidth = 0.6, inherit.aes = FALSE
  ) +
  geom_label(
    data = load_df_top,
    aes(PC1 * arrow_scale * 1.1,
        PC2 * arrow_scale * 1.1,
        label = land_class),
    size = 2.5, colour = "grey10", fill = alpha("white", 0.7),
    label.size = NA, inherit.aes = FALSE
  ) +
  scale_colour_manual(values = survey_colours, na.value = "grey50") +
  labs(x = paste0("PC1 (", round(100 * summary(pca_lc)$importance[2, 1], 1), "%)"),
       y = paste0("PC2 (", round(100 * summary(pca_lc)$importance[2, 2], 1), "%)"),
       colour = "Survey/Range",
       title  = "Land cover PCA — sample locations") +
  theme_pub

ggsave(file.path(outdir, "Fig_LandCover_PCA_Biplot.png"),
       p_biplot, width = 8, height = 6, dpi = 300)


# ============================================================
# 3) GAMs — CORTISOL ~ LAND COVER
#
# Step 1 — model selection (AIC):
#   (a) Null
#   (b) Full PCA: LC_PC1–PC4  ← primary reported model
#
# Step 2 — interpretive named-class model:
#   One top-loading class per PC, chosen to minimise collinearity:
#     PC1 (+): OpenBog     (loading 0.343; open peatland gradient)
#     PC2 (−): TreedPeat   (loading −0.461; treed peatland)
#     PC3 (−): ConiferousSwamp (loading −0.468; conifer swamp)
#     PC4 (+): Deciduous    (loading 0.402; deciduous/upland; Marsh filtered out as near-zero)
#   Used for partial-effect figures and ecological interpretation only.
# ============================================================

# sqrt-transform OpenBog and Deciduous (right-skewed proportions) for better knot placement
dat_lc <- dat_lc %>%
  mutate(sqrt_OpenBog   = sqrt(OpenBog),
         sqrt_Deciduous = sqrt(Deciduous))

dat_lc_cort <- dat_lc %>% filter(is.finite(cortisol), !is.na(julian))

# (a) Null
m_cort_nolc <- gam(
  cortisol ~
    s(survey, bs = "re") +
    year_f + preg_group +
    s(julian, k = 8, bs = "tp"),
  data   = dat_lc_cort,
  family = Gamma(link = "log"),
  method = "REML"
)

# (b) Full PCA — primary model
m_cort_pca <- gam(
  cortisol ~
    s(survey, bs = "re") +
    year_f + preg_group +
    s(LC_PC1, k = 3, bs = "tp") +
    s(LC_PC2, k = 3, bs = "tp") +
    s(LC_PC3, k = 3, bs = "tp") +
    s(LC_PC4, k = 3, bs = "tp") +
    s(julian, k = 8, bs = "tp"),
  data   = dat_lc_cort,
  family = Gamma(link = "log"),
  method = "REML"
)

cat("\n=== Cortisol: null vs PCA model (AIC) ===\n")
aic_cort <- AIC(m_cort_nolc, m_cort_pca)
aic_cort$deltaAIC <- aic_cort$AIC - min(aic_cort$AIC)
print(round(aic_cort, 2))
capture.output(print(round(aic_cort, 2)), file = file.path(outdir, "GAM_AIC_05_cortisol_lc.txt"))
cat("\n--- PCA model summary ---\n"); summary(m_cort_pca)
capture.output(summary(m_cort_pca), file = file.path(outdir, "GAM_summary_05_cortisol_lc_pca.txt"))

# GAM plot for primary model
png(file.path(outdir, "Fig_GAM_Cortisol_LandCover.png"), 2400, 1200, res = 300)
plot(m_cort_pca, pages = 1, scheme = 1,
     main = "Cortisol ~ Land cover PCs (partial effects)")
dev.off()

cort_lc_year_contrasts <- emmeans(m_cort_pca, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(cort_lc_year_contrasts), file = file.path(outdir, "GAM_contrasts_05_cortisol_lc_year.txt"))

# (c) Interpretive named-class model — one representative per PC
# PC1: sqrt(OpenBog) (+), PC2: TreedPeat (−), PC3: ConiferousSwamp (−), PC4: sqrt(Deciduous) (+)
# OpenBog and Deciduous are sqrt-transformed (right-skewed proportions)
m_cort_interp <- gam(
  cortisol ~
    s(survey, bs = "re") +
    year_f + preg_group +
    s(sqrt_OpenBog,    k = 6, bs = "tp") +
    s(TreedPeat,       k = 6, bs = "tp") +
    s(ConiferousSwamp, k = 6, bs = "tp") +
    s(sqrt_Deciduous,  k = 6, bs = "tp") +
    s(julian,          k = 8, bs = "tp"),
  data   = dat_lc_cort,
  family = Gamma(link = "log"),
  method = "REML"
)
cat("\n--- Interpretive named-class model (cortisol) ---\n")
summary(m_cort_interp)
capture.output(summary(m_cort_interp), file = file.path(outdir, "GAM_summary_05_cortisol_lc_interp.txt"))

# Partial effect plots for each named class (gratia for cleaner output)
# Note: x-axis is on sqrt scale; original proportion labels added for publication in 09_manuscript_figures.R
p_cort_interp <- draw(m_cort_interp,
                      select = c("s(sqrt_OpenBog)", "s(TreedPeat)",
                                 "s(ConiferousSwamp)", "s(sqrt_Deciduous)"),
                      residuals = FALSE) &
  theme_bw(base_size = 10) &
  labs(y = "Partial effect on log(cortisol)")

ggsave(file.path(outdir, "Fig_GAM_Cortisol_LandCover_Interp.png"),
       p_cort_interp, width = 10, height = 8, dpi = 300)

# ============================================================
# 4) GAMs — T3 ~ LAND COVER
# ============================================================

dat_lc_t3 <- dat_lc %>% filter(is.finite(t3), t3 > 0, !is.na(julian))

# (a) Null
m_t3_nolc <- gam(
  t3 ~
    s(survey, bs = "re") +
    year_f + preg_group +
    s(julian, k = 8, bs = "tp"),
  data   = dat_lc_t3,
  family = Gamma(link = "log"),
  method = "REML"
)

# (b) Full PCA — primary model
m_t3_pca <- gam(
  t3 ~
    s(survey, bs = "re") +
    year_f + preg_group +
    s(LC_PC1, k = 3, bs = "tp") +
    s(LC_PC2, k = 3, bs = "tp") +
    s(LC_PC3, k = 3, bs = "tp") +
    s(LC_PC4, k = 3, bs = "tp") +
    s(julian, k = 8, bs = "tp"),
  data   = dat_lc_t3,
  family = Gamma(link = "log"),
  method = "REML"
)

cat("\n=== T3: null vs PCA model (AIC) ===\n")
aic_t3 <- AIC(m_t3_nolc, m_t3_pca)
aic_t3$deltaAIC <- aic_t3$AIC - min(aic_t3$AIC)
print(round(aic_t3, 2))
capture.output(print(round(aic_t3, 2)), file = file.path(outdir, "GAM_AIC_05_t3_lc.txt"))
cat("\n--- PCA model summary ---\n"); summary(m_t3_pca)
capture.output(summary(m_t3_pca), file = file.path(outdir, "GAM_summary_05_t3_lc_pca.txt"))

# GAM plot for primary model
png(file.path(outdir, "Fig_GAM_T3_LandCover.png"), 2400, 1200, res = 300)
plot(m_t3_pca, pages = 1, scheme = 1,
     main = "T3 ~ Land cover PCs (partial effects)")
dev.off()

t3_lc_year_contrasts <- emmeans(m_t3_pca, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(t3_lc_year_contrasts), file = file.path(outdir, "GAM_contrasts_05_t3_lc_year.txt"))

# (c) Interpretive named-class model — same representatives as cortisol (with sqrt transforms)
m_t3_interp <- gam(
  t3 ~
    s(survey, bs = "re") +
    year_f + preg_group +
    s(sqrt_OpenBog,    k = 6, bs = "tp") +
    s(TreedPeat,       k = 6, bs = "tp") +
    s(ConiferousSwamp, k = 6, bs = "tp") +
    s(sqrt_Deciduous,  k = 6, bs = "tp") +
    s(julian,          k = 8, bs = "tp"),
  data   = dat_lc_t3,
  family = Gamma(link = "log"),
  method = "REML"
)
cat("\n--- Interpretive named-class model (T3) ---\n")
summary(m_t3_interp)
capture.output(summary(m_t3_interp), file = file.path(outdir, "GAM_summary_05_t3_lc_interp.txt"))

p_t3_interp <- draw(m_t3_interp,
                    select = c("s(sqrt_OpenBog)", "s(TreedPeat)",
                               "s(ConiferousSwamp)", "s(sqrt_Deciduous)"),
                    residuals = FALSE) &
  theme_bw(base_size = 10) &
  labs(y = "Partial effect on log(T3)")

ggsave(file.path(outdir, "Fig_GAM_T3_LandCover_Interp.png"),
       p_t3_interp, width = 10, height = 8, dpi = 300)

# ============================================================
# 5) SPATIAL AUTOCORRELATION DIAGNOSTICS (variogram on residuals)
# Purpose: same diagnostic applied to the road-distance models in
#   03_winter_road_analyses.R (and now also to mine distance, wolf
#   suitability, and population density in 07/08/06), extended here to
#   the primary reported land-cover PCA models (m_cort_pca / m_t3_pca)
#   — flagged in the external QC report as a concern especially for
#   FT3~land-cover, where the effect is large (ΔAIC = 486.7).
#
# Some samples share IDENTICAL GPS coordinates (repeat visits to the
# same collection site). A raw per-sample variogram would treat those
# as separate zero-distance observations and read as pure nugget noise
# rather than genuine small-scale variation, so residuals are averaged
# to one value per unique sampling location before the variogram is
# computed. If the variogram shows real spatial structure (sill clearly
# above the nugget, with a non-trivial range), the land-cover PCA
# models are refit with an added s(x, y) smooth on the projected
# coordinates, and the LC_PC1-4 effects are re-checked with that term
# included.
# ============================================================

if (!requireNamespace("gstat", quietly = TRUE)) install.packages("gstat", repos = "https://cloud.r-project.org")
library(gstat)

sf_cort_lc <- dat_lc_cort %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(wolfstack))
sf_t3_lc <- dat_lc_t3 %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs(wolfstack))

check_spatial_autocorr_lc <- function(resid_vec, data_sf, label) {
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

  png(file.path(outdir, paste0("Variogram_", label, "_LandCoverModel.png")), 1600, 1200, res = 200)
  print(plot(vgram, main = paste0(label, " — variogram of land-cover-model deviance residuals"),
             xlab = "Distance (m)", ylab = "Semivariance"))
  dev.off()

  vgram_fit <- tryCatch(gstat::fit.variogram(vgram, gstat::vgm("Sph")), error = function(e) NULL)
  if (!is.null(vgram_fit)) {
    cat(sprintf("%s variogram model (nugget / partial sill / range):\n", label))
    print(vgram_fit)
    capture.output(print(vgram_fit),
                   file = file.path(outdir, paste0("Variogram_fit_", label, "_LandCoverModel.txt")))
  } else {
    message(label, ": variogram model fit did not converge — inspect the empirical variogram plot directly.")
  }
  list(vgram = vgram, vgram_fit = vgram_fit)
}

vres_cort_lc <- check_spatial_autocorr_lc(residuals(m_cort_pca, type = "deviance"), sf_cort_lc, "FCM")
vres_t3_lc   <- check_spatial_autocorr_lc(residuals(m_t3_pca,   type = "deviance"), sf_t3_lc,   "FT3")

# ---- Refit with s(x, y) added, and re-check whether LC_PC1-4 survive ----

dat_cort_lc_xy <- dat_lc_cort %>%
  bind_cols(st_coordinates(sf_cort_lc) %>% as.data.frame() %>% rename(x_proj = X, y_proj = Y))
dat_t3_lc_xy <- dat_lc_t3 %>%
  bind_cols(st_coordinates(sf_t3_lc) %>% as.data.frame() %>% rename(x_proj = X, y_proj = Y))

m_cort_pca_xy <- gam(
  cortisol ~
    s(survey, bs = "re") +
    year_f + preg_group +
    s(LC_PC1, k = 3, bs = "tp") +
    s(LC_PC2, k = 3, bs = "tp") +
    s(LC_PC3, k = 3, bs = "tp") +
    s(LC_PC4, k = 3, bs = "tp") +
    s(julian, k = 8, bs = "tp") +
    s(x_proj, y_proj, bs = "tp", k = 30),
  data   = dat_cort_lc_xy,
  family = Gamma(link = "log"),
  method = "REML"
)

plot(m_cort_pca_xy, pages = 1, scheme = 2)

m_t3_pca_xy <- gam(
  t3 ~
    s(survey, bs = "re") +
    year_f + preg_group +
    s(LC_PC1, k = 3, bs = "tp") +
    s(LC_PC2, k = 3, bs = "tp") +
    s(LC_PC3, k = 3, bs = "tp") +
    s(LC_PC4, k = 3, bs = "tp") +
    s(julian, k = 8, bs = "tp") +
    s(x_proj, y_proj, bs = "tp", k = 30),
  data   = dat_t3_lc_xy,
  family = Gamma(link = "log"),
  method = "REML"
)

plot(m_t3_pca_xy, pages = 1, scheme = 2)

cat("\nLand-cover PCA models with s(x, y) added — AIC/BIC comparison against the originally reported (no-spatial-term) models:\n")
aic_cort_lc_xy <- AIC(m_cort_pca, m_cort_pca_xy)
aic_t3_lc_xy   <- AIC(m_t3_pca,   m_t3_pca_xy)
bic_cort_lc_xy <- BIC(m_cort_pca, m_cort_pca_xy)
bic_t3_lc_xy   <- BIC(m_t3_pca,   m_t3_pca_xy)
print(aic_cort_lc_xy); print(bic_cort_lc_xy)
print(aic_t3_lc_xy); print(bic_t3_lc_xy)
capture.output(print(aic_cort_lc_xy), print(bic_cort_lc_xy), print(aic_t3_lc_xy), print(bic_t3_lc_xy),
               file = file.path(outdir, "GAM_AIC_05_spatial_xy_comparison.txt"))

cat("\nDo LC_PC1-4 remain significant once s(x, y) is included?\n")
print(summary(m_cort_pca_xy)$s.table)
print(summary(m_t3_pca_xy)$s.table)
capture.output(summary(m_cort_pca_xy), file = file.path(outdir, "GAM_summary_05_cortisol_lc_pca_xy.txt"))
capture.output(summary(m_t3_pca_xy),   file = file.path(outdir, "GAM_summary_05_t3_lc_pca_xy.txt"))

# Re-check the variogram on the xy-augmented models' residuals — if
# s(x, y) has done its job, this should look flatter than the original.
vres_cort_lc_after <- check_spatial_autocorr_lc(residuals(m_cort_pca_xy, type = "deviance"), sf_cort_lc, "FCM_afterXY")
vres_t3_lc_after   <- check_spatial_autocorr_lc(residuals(m_t3_pca_xy,   type = "deviance"), sf_t3_lc,   "FT3_afterXY")

message("Part 5 spatial autocorrelation diagnostics complete — outputs in: ", outdir)

# ============================================================
# 6) PRIMARY MODEL DECISION — spatially-adjusted model reported
# Project decision: for landscape covariates (road, mine, wolf, land
# cover, population), the s(x, y)-augmented model is reported as the
# PRIMARY model, regardless of BIC (see 03_winter_road_analyses.R
# Part 2E for full rationale). This applies to BOTH land-cover models
# used in this project: the full PCA model (m_cort_pca / m_t3_pca —
# used for AIC-based model selection) AND the interpretive named-class
# model (m_cort_interp / m_t3_interp — used for the partial-effect
# figures reported in 09_manuscript_figures.R Figure 3). The
# interpretive model did not yet have an xy-augmented version, so it is
# built here using the same projected coordinates already attached to
# dat_cort_lc_xy / dat_t3_lc_xy above.
#
# The outputs below OVERWRITE the original (non-spatial) GAM_summary/
# Fig/contrasts files from Parts 3-4 so every file on disk — including
# the ones 09_manuscript_figures.R draws from — reflects the model
# actually reported in the manuscript.
# ============================================================

# --- 6A) PCA model (AIC-selection reporting) ---

summary(m_cort_pca_xy)
capture.output(summary(m_cort_pca_xy), file = file.path(outdir, "GAM_summary_05_cortisol_lc_pca.txt"))

png(file.path(outdir, "Fig_GAM_Cortisol_LandCover.png"), 2400, 1200, res = 300)
plot(m_cort_pca_xy, pages = 1, scheme = 1,
     main = "Cortisol ~ Land cover PCs + s(x,y) (partial effects)")
dev.off()

cort_lc_year_contrasts <- emmeans(m_cort_pca_xy, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(cort_lc_year_contrasts), file = file.path(outdir, "GAM_contrasts_05_cortisol_lc_year.txt"))

summary(m_t3_pca_xy)
capture.output(summary(m_t3_pca_xy), file = file.path(outdir, "GAM_summary_05_t3_lc_pca.txt"))

png(file.path(outdir, "Fig_GAM_T3_LandCover.png"), 2400, 1200, res = 300)
plot(m_t3_pca_xy, pages = 1, scheme = 1,
     main = "T3 ~ Land cover PCs + s(x,y) (partial effects)")
dev.off()

t3_lc_year_contrasts <- emmeans(m_t3_pca_xy, ~ year_f) %>% pairs(adjust = "tukey")
capture.output(print(t3_lc_year_contrasts), file = file.path(outdir, "GAM_contrasts_05_t3_lc_year.txt"))

# --- 6B) Interpretive named-class model (figure reporting) ---
# Built fresh here with s(x, y) added, using the same representative
# classes as the original interpretive model (sqrt_OpenBog, TreedPeat,
# ConiferousSwamp, sqrt_Deciduous).

m_cort_interp_xy <- gam(
  cortisol ~
    s(survey, bs = "re") +
    year_f + preg_group +
    s(sqrt_OpenBog,    k = 6, bs = "tp") +
    s(TreedPeat,       k = 6, bs = "tp") +
    s(ConiferousSwamp, k = 6, bs = "tp") +
    s(sqrt_Deciduous,  k = 6, bs = "tp") +
    s(julian,          k = 8, bs = "tp") +
    s(x_proj, y_proj,  bs = "tp", k = 30),
  data   = dat_cort_lc_xy,
  family = Gamma(link = "log"),
  method = "REML"
)
cat("\n--- Interpretive named-class model + s(x,y) (cortisol) ---\n")
summary(m_cort_interp_xy)
capture.output(summary(m_cort_interp_xy), file = file.path(outdir, "GAM_summary_05_cortisol_lc_interp.txt"))

p_cort_interp <- draw(m_cort_interp_xy,
                      select = c("s(sqrt_OpenBog)", "s(TreedPeat)",
                                 "s(ConiferousSwamp)", "s(sqrt_Deciduous)"),
                      residuals = FALSE) &
  theme_bw(base_size = 10) &
  labs(y = "Partial effect on log(cortisol)")

ggsave(file.path(outdir, "Fig_GAM_Cortisol_LandCover_Interp.png"),
       p_cort_interp, width = 10, height = 8, dpi = 300)

m_t3_interp_xy <- gam(
  t3 ~
    s(survey, bs = "re") +
    year_f + preg_group +
    s(sqrt_OpenBog,    k = 6, bs = "tp") +
    s(TreedPeat,       k = 6, bs = "tp") +
    s(ConiferousSwamp, k = 6, bs = "tp") +
    s(sqrt_Deciduous,  k = 6, bs = "tp") +
    s(julian,          k = 8, bs = "tp") +
    s(x_proj, y_proj,  bs = "tp", k = 30),
  data   = dat_t3_lc_xy,
  family = Gamma(link = "log"),
  method = "REML"
)
cat("\n--- Interpretive named-class model + s(x,y) (T3) ---\n")
summary(m_t3_interp_xy)
capture.output(summary(m_t3_interp_xy), file = file.path(outdir, "GAM_summary_05_t3_lc_interp.txt"))

p_t3_interp <- draw(m_t3_interp_xy,
                    select = c("s(sqrt_OpenBog)", "s(TreedPeat)",
                               "s(ConiferousSwamp)", "s(sqrt_Deciduous)"),
                    residuals = FALSE) &
  theme_bw(base_size = 10) &
  labs(y = "Partial effect on log(T3)")

ggsave(file.path(outdir, "Fig_GAM_T3_LandCover_Interp.png"),
       p_t3_interp, width = 10, height = 8, dpi = 300)

message("Part 6 (primary model = spatially-adjusted) complete — GAM_summary_05_*_lc_pca.txt and GAM_summary_05_*_lc_interp.txt now reflect the s(x,y)-augmented models.")

message("Part 5 complete — outputs in: ", outdir)

