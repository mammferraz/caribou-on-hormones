# ============================================================
# 01_general_analyses.R
#
# PART 1 — GENERAL DESCRIPTIVE ANALYSES
# Purpose:
#   • Dataset QC and sample structure summaries
#   • Hormone distribution plots (by year, coloured by survey)
#   • Sampling maps
#   • GAMs for hormone variation by pregnancy group and year
#   • Model-based baseline hormone predictions
#   • Baseline exceedance summaries
# ============================================================

source("00_data_processing.R")

outdir <- "caribou_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1) SAMPLE STRUCTURE SUMMARIES
# ============================================================

# 1A. Samples per year × survey × sex
counts_year_survey <- dat %>%
  count(year, survey, sex, name = "n_samples") %>%
  arrange(year, survey, sex)
write_csv(counts_year_survey, file.path(outdir, "Summary_counts_by_year_survey_sex.csv"))

# 1B. Unique animals per year × survey
unique_animals <- dat %>%
  distinct(year, survey, sex, pcid) %>%
  count(year, survey, sex, name = "n_unique_animals")
write_csv(unique_animals, file.path(outdir, "Summary_unique_animals_by_year_survey.csv"))

# 1C. Repeated sampling
repeats <- dat %>%
  filter(!is.na(pcid)) %>%
  count(pcid, year, survey, name = "n_samples_this_animal") %>%
  group_by(year, survey) %>%
  summarise(
    n_unique_animals = n(),
    n_repeated       = sum(n_samples_this_animal > 1),
    max_samples      = max(n_samples_this_animal),
    .groups          = "drop"
  )
write_csv(repeats, file.path(outdir, "Summary_repeated_animals_by_year_survey.csv"))

# 1D. Pregnant females
pregnant_summary <- dat %>%
  filter(sex == "F") %>%
  group_by(year, survey) %>%
  summarise(
    n_female_samples   = n(),
    n_preg_samples     = sum(pregnancy_status_num == 1, na.rm = TRUE),
    n_unique_females   = n_distinct(pcid),
    n_unique_pregnant  = n_distinct(pcid[pregnancy_status_num == 1]),
    .groups            = "drop"
  )
write_csv(pregnant_summary, file.path(outdir, "Summary_pregnant_females_by_year_survey.csv"))

# 1E. Hormone summaries by year × survey × pregnancy group
hormone_summary <- dat %>%
  group_by(year, survey, preg_group) %>%
  summarise(
    n                   = n(),
    cortisol_mean       = mean(cortisol,     na.rm = TRUE),
    cortisol_sd         = sd(cortisol,       na.rm = TRUE),
    cortisol_median     = median(cortisol,   na.rm = TRUE),
    progesterone_mean   = mean(progesterone, na.rm = TRUE),
    progesterone_sd     = sd(progesterone,   na.rm = TRUE),
    t3_mean             = mean(t3,           na.rm = TRUE),
    t3_sd               = sd(t3,             na.rm = TRUE),
    .groups             = "drop"
  )
write_csv(hormone_summary, file.path(outdir, "Summary_hormones_by_year_survey_group.csv"))

# ============================================================
# 2) HORMONE DISTRIBUTION PLOTS
# ============================================================

# 2A) Overall distributions by pregnancy group (all data pooled)
plot_hormone_overall <- function(data, yvar, ylab, ybreaks, fname) {
  p <- ggplot(
    data %>% filter(is.finite(.data[[yvar]])),
    aes(x = preg_group, y = .data[[yvar]],
        colour = preg_group, shape = preg_group)
  ) +
    geom_violin(fill = "grey85", alpha = 0.6, trim = FALSE, colour = NA) +
    geom_jitter(width = 0.18, alpha = 0.5, size = 2) +
    geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.8, fill = NA,
                 colour = "black") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                 fill = "red", colour = "red") +
    scale_shape_manual(values  = preg_shapes) +
    scale_colour_manual(values = preg_colours) +
    scale_y_log10(breaks = ybreaks) +
    labs(x = NULL, y = ylab) +
    theme_pub + theme(legend.position = "none")
  ggsave(file.path(outdir, fname), p, width = 5, height = 6, dpi = 300)
  invisible(p)
}

plot_hormone_overall(dat, "cortisol",    "Cortisol (log scale)",
                     c(30, 100, 300, 1000),       "Fig_Cortisol_by_Group.png")
plot_hormone_overall(dat, "t3",          "T3 (log scale)",
                     c(10, 30, 100, 300),          "Fig_T3_by_Group.png")
plot_hormone_overall(dat, "progesterone","Progesterone (log scale)",
                     c(100, 500, 1000, 5000, 10000), "Fig_Progesterone_by_Group.png")


# 2B) Distributions by YEAR (x-axis), coloured by survey
#     Faceted by pregnancy group
plot_hormone_by_year <- function(data, yvar, ylab, ybreaks, fname) {
  d <- data %>%
    filter(is.finite(.data[[yvar]]), !is.na(preg_group)) %>%
    mutate(year_f = factor(year))

  p <- ggplot(d, aes(x = year_f, y = .data[[yvar]], colour = survey)) +
    geom_jitter(
      position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.7),
      alpha = 0.6, size = 1.8
    ) +
    geom_boxplot(
      aes(group = interaction(year_f, survey)),
      position = position_dodge(width = 0.7),
      outlier.shape = NA, fill = NA, linewidth = 0.7, width = 0.5
    ) +
    stat_summary(
      aes(group = interaction(year_f, survey)),
      fun = mean, geom = "point", shape = 23, size = 3,
      fill = "white", colour = "black",
      position = position_dodge(width = 0.7)
    ) +
    facet_wrap(~ preg_group, ncol = 1) +
    scale_colour_manual(values = survey_colours, na.value = "grey50") +
    scale_y_log10(breaks = ybreaks) +
    labs(x = "Year", y = ylab, colour = "Survey/Range") +
    theme_pub +
    theme(axis.text.x = element_text(size = 10))
  ggsave(file.path(outdir, fname), p,
         width = max(6, 1.5 * n_distinct(d$year_f)), height = 10, dpi = 300)
  invisible(p)
}

plot_hormone_by_year(dat, "cortisol",    "Cortisol (log scale)",
                     c(30, 100, 300, 1000),        "Fig_Cortisol_by_Year_Survey.png")
plot_hormone_by_year(dat, "t3",          "T3 (log scale)",
                     c(10, 30, 100, 300),           "Fig_T3_by_Year_Survey.png")
plot_hormone_by_year(dat, "progesterone","Progesterone (log scale)",
                     c(100, 500, 1000, 5000, 10000),"Fig_Progesterone_by_Year_Survey.png")


# ============================================================
# 3) SAMPLING MAP
# ============================================================

land      <- ne_countries(scale = "large", returnclass = "sf") %>% st_transform(4326)
provinces <- ne_states(country = "canada", returnclass = "sf") %>% st_transform(4326)
lakes     <- ne_download(scale = "large", type = "lakes",
                         category = "physical", returnclass = "sf") %>% st_transform(4326)
rivers    <- ne_download(scale = "large", type = "rivers_lake_centerlines",
                         category = "physical", returnclass = "sf") %>% st_transform(4326)

pts <- st_as_sf(dat_map, coords = c("longitude", "latitude"),
                crs = 4326, remove = FALSE)
bb <- st_bbox(pts)
bb[c("xmin","xmax")] <- bb[c("xmin","xmax")] + c(-1, 1)
bb[c("ymin","ymax")] <- bb[c("ymin","ymax")] + c(-1, 1)

p_map <- ggplot() +
  geom_sf(data = land,      fill = "grey95",    colour = NA) +
  geom_sf(data = lakes,     fill = "lightblue", colour = NA, alpha = 0.7) +
  geom_sf(data = rivers,    colour = "lightblue3", linewidth = 0.3) +
  geom_sf(data = provinces, fill = NA, colour = "grey40", linewidth = 0.3) +
  geom_sf(data = caribou_ranges, fill = "forestgreen", colour = "darkgreen",
          alpha = 0.15, linewidth = 0.5) +
  geom_sf_text(
    data          = caribou_ranges,
    aes(label     = .data[[range_name_col]]),
    size          = 2.5, colour = "darkgreen",
    fontface      = "italic", check_overlap = TRUE
  ) +
  geom_point(
    data  = st_drop_geometry(pts),
    aes(longitude, latitude, colour = survey, shape = preg_group),
    alpha = 0.7, size = 2.5
  ) +
  scale_shape_manual(values = preg_shapes) +
  scale_colour_manual(values = survey_colours, na.value = "grey50") +
  coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
           ylim = c(bb["ymin"], bb["ymax"]), expand = FALSE) +
  labs(
    title  = "Caribou sampling locations and boreal caribou range boundaries",
    x = "Longitude", y = "Latitude",
    colour = "Survey/Range", shape = "Group"
  ) +
  theme_minimal() + theme(panel.grid = element_blank())

ggsave(file.path(outdir, "Map_Sampling_Locations.png"),
       p_map, width = 9, height = 8, dpi = 300)


# ============================================================
# 4) GAMs — HORMONE ~ PREGNANCY GROUP + YEAR + s(JULIAN)
#    Year is a fixed numeric effect capturing interannual trends.
#    Survey is a random effect capturing range-level baseline.
# ============================================================

dat_gam <- dat %>%
  mutate(
    survey     = factor(survey),
    year_num   = as.numeric(year),   # continuous year for trend smooth
    year_f     = factor(year),       # factor year for pairwise contrasts
    preg_group = factor(preg_group,
                        levels = c("Male", "Female_nonpregnant", "Female_pregnant"))
  )

# ------------------------------------------------------------
# 4A) Cortisol
# ------------------------------------------------------------

dat_cort <- dat_gam %>%
  filter(is.finite(cortisol), !is.na(survey), !is.na(preg_group), !is.na(julian))

m_cortisol <- gam(
  cortisol ~
    s(survey,   bs = "re") +   # range-level random intercept
    year_f +                    # interannual variation
    preg_group +                # sex/pregnancy main effect
    s(julian,   k = 20),       # seasonal smooth (k raised from 8: edf near limit)
  data   = dat_cort,
  family = Gamma(link = "log"),
  method = "REML"
)

summary(m_cortisol)
capture.output(summary(m_cortisol), file = file.path(outdir, "GAM_summary_01_cortisol_preggroup.txt"))
gam.check(m_cortisol)

# Pairwise contrasts: pregnancy group and year
emm_cort_group <- emmeans(m_cortisol, ~ preg_group)
emm_cort_year  <- emmeans(m_cortisol, ~ year_f)
cort_group_contrasts <- pairs(emm_cort_group, adjust = "tukey")
cort_year_contrasts  <- pairs(emm_cort_year,  adjust = "tukey")
capture.output(print(cort_group_contrasts), file = file.path(outdir, "GAM_contrasts_01_cortisol_preggroup.txt"))
capture.output(print(cort_year_contrasts),  file = file.path(outdir, "GAM_contrasts_01_cortisol_year.txt"))

write_csv(as.data.frame(emm_cort_group),
          file.path(outdir, "EMM_Cortisol_by_PregGroup.csv"))
write_csv(as.data.frame(emm_cort_year),
          file.path(outdir, "EMM_Cortisol_by_Year.csv"))

# ------------------------------------------------------------
# 4B) T3
# ------------------------------------------------------------

dat_t3 <- dat_gam %>%
  filter(is.finite(t3), t3 > 0, !is.na(survey), !is.na(preg_group), !is.na(julian))

m_t3 <- gam(
  t3 ~
    s(survey,   bs = "re") +
    year_f +
    preg_group +
    s(julian,   k = 20),
  data   = dat_t3,
  family = Gamma(link = "log"),
  method = "REML"
)

# Note: a preg_group × julian interaction model (bs="sz") was tested and
# rejected — additive model beat it by >25 AIC units once k was sufficient.

summary(m_t3)
capture.output(summary(m_t3), file = file.path(outdir, "GAM_summary_01_t3_preggroup.txt"))
emm_t3_group <- emmeans(m_t3, ~ preg_group)
emm_t3_year  <- emmeans(m_t3, ~ year_f)
t3_group_contrasts <- pairs(emm_t3_group, adjust = "tukey")
t3_year_contrasts  <- pairs(emm_t3_year,  adjust = "tukey")
capture.output(print(t3_group_contrasts), file = file.path(outdir, "GAM_contrasts_01_t3_preggroup.txt"))
capture.output(print(t3_year_contrasts),  file = file.path(outdir, "GAM_contrasts_01_t3_year.txt"))

write_csv(as.data.frame(emm_t3_group),
          file.path(outdir, "EMM_T3_by_PregGroup.csv"))
write_csv(as.data.frame(emm_t3_year),
          file.path(outdir, "EMM_T3_by_Year.csv"))

# ------------------------------------------------------------
# 4C) Progesterone (females only)
# ------------------------------------------------------------

dat_prog <- dat_gam %>%
  filter(sex == "F",
         preg_group %in% c("Female_nonpregnant", "Female_pregnant"),
         is.finite(progesterone), progesterone > 0,
         !is.na(survey), !is.na(julian)) %>%
  mutate(preg_group = droplevels(preg_group))

m_prog <- gam(
  progesterone ~
    s(survey,   bs = "re") +
    year_f +
    preg_group +
    s(julian,   k = 20),
  data   = dat_prog,
  family = Gamma(link = "log"),
  method = "REML"
)

summary(m_prog)
capture.output(summary(m_prog), file = file.path(outdir, "GAM_summary_01_prog_preggroup.txt"))
emm_prog_group <- emmeans(m_prog, ~ preg_group)
prog_group_contrasts <- pairs(emm_prog_group, adjust = "tukey")
capture.output(print(prog_group_contrasts), file = file.path(outdir, "GAM_contrasts_01_prog_preggroup.txt"))

write_csv(as.data.frame(emm_prog_group),
          file.path(outdir, "EMM_Progesterone_by_PregGroup.csv"))


# ============================================================
# 5) MODEL-BASED BASELINE PREDICTIONS (marginalised over
#    survey and julian day)
# ============================================================

ref_year <- levels(dat_gam$year_f)[1]

pred_df <- expand.grid(
  preg_group = levels(dat_gam$preg_group),
  julian     = median(dat_gam$julian, na.rm = TRUE),
  survey     = factor(dat_gam$survey[1], levels = levels(factor(dat_gam$survey))),
  year_f     = factor(ref_year, levels = levels(dat_gam$year_f))
)

# Cortisol predictions
pred_cort <- predict(m_cortisol, newdata = pred_df,
                     exclude = c("s(survey)", "s(julian)"),
                     type = "link", se.fit = TRUE)
# T3 predictions
pred_t3 <- predict(m_t3, newdata = pred_df,
                   exclude = c("s(survey)", "s(julian)"),
                   type = "link", se.fit = TRUE)

hormone_predictions <- pred_df %>%
  mutate(
    cortisol_hat = exp(as.numeric(pred_cort$fit)),
    cortisol_lo  = exp(as.numeric(pred_cort$fit) - 1.96 * as.numeric(pred_cort$se.fit)),
    cortisol_hi  = exp(as.numeric(pred_cort$fit) + 1.96 * as.numeric(pred_cort$se.fit)),
    t3_hat       = exp(as.numeric(pred_t3$fit)),
    t3_lo        = exp(as.numeric(pred_t3$fit) - 1.96 * as.numeric(pred_t3$se.fit)),
    t3_hi        = exp(as.numeric(pred_t3$fit) + 1.96 * as.numeric(pred_t3$se.fit))
  ) %>%
  select(preg_group, starts_with("cortisol"), starts_with("t3"))

write_csv(hormone_predictions,
          file.path(outdir, "Baseline_hormone_predictions.csv"))

# Plot predictions
plot_df <- hormone_predictions %>%
  pivot_longer(-preg_group,
               names_to  = c("hormone", ".value"),
               names_pattern = "(cortisol|t3)_(.*)") %>%
  mutate(hormone = recode(hormone, cortisol = "Cortisol", t3 = "T3"))

p_baseline <- ggplot(plot_df, aes(x = preg_group, y = hat)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.9) +
  facet_wrap(~ hormone, scales = "fixed") +
  scale_y_continuous(limits = c(0, 200)) +
  labs(x = NULL, y = "Predicted baseline concentration",
       title = "Model-based baseline hormones by pregnancy group",
       subtitle = "GAM predictions marginalised over survey and Julian day") +
  theme_classic(base_size = 13) +
  theme(strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(outdir, "Fig_Baseline_Hormone_Predictions.png"),
       p_baseline, width = 8, height = 4.5, dpi = 300)


# ============================================================
# 6) BASELINE EXCEEDANCE
# ============================================================

dat_baseline <- dat %>%
  filter(is.finite(cortisol) | is.finite(t3),
         !is.na(preg_group), !is.na(survey), !is.na(julian)) %>%
  left_join(hormone_predictions, by = "preg_group") %>%
  mutate(
    delta_cortisol  = if_else(is.finite(cortisol), cortisol - cortisol_hat, NA_real_),
    delta_t3        = if_else(is.finite(t3),        t3 - t3_hat,            NA_real_),
    cort_above_base = if_else(is.finite(cortisol),  cortisol > cortisol_hat, NA),
    cort_above_ci   = if_else(is.finite(cortisol),  cortisol > cortisol_hi,  NA),
    t3_above_base   = if_else(is.finite(t3),        t3 > t3_hat,             NA),
    t3_above_ci     = if_else(is.finite(t3),        t3 > t3_hi,              NA)
  )

# Summary by year × survey
summary_exceedance <- dat_baseline %>%
  group_by(year, survey, preg_group) %>%
  summarise(
    n                   = n(),
    n_cort              = sum(is.finite(cortisol)),
    n_t3                = sum(is.finite(t3)),
    pct_cort_above_base = mean(cort_above_base, na.rm = TRUE) * 100,
    pct_cort_above_ci   = mean(cort_above_ci,   na.rm = TRUE) * 100,
    pct_t3_above_base   = mean(t3_above_base,   na.rm = TRUE) * 100,
    pct_t3_above_ci     = mean(t3_above_ci,     na.rm = TRUE) * 100,
    .groups             = "drop"
  )
write_csv(summary_exceedance,
          file.path(outdir, "Summary_Baseline_Exceedance_by_Year_Survey.csv"))

# Plot: observed cortisol vs baseline, x = julian day, colour = survey, facet = year × preg_group
p_cort_baseline <- ggplot(
  dat_baseline %>% filter(is.finite(cortisol)),
  aes(julian, cortisol, colour = survey)
) +
  geom_rect(aes(xmin = -Inf, xmax = Inf,
                ymin = cortisol_lo, ymax = cortisol_hi),
            fill = "grey80", alpha = 0.3, inherit.aes = FALSE) +
  geom_hline(aes(yintercept = cortisol_hat), linewidth = 1) +
  geom_point(alpha = 0.6, size = 2) +
  facet_grid(preg_group ~ year) +
  scale_y_log10(breaks = c(10, 100, 200, 800)) +
  scale_colour_manual(values = survey_colours, na.value = "grey50") +
  labs(title = "Observed cortisol relative to model-based baseline",
       x = "Julian day", y = "Cortisol (log scale)", colour = "Survey/Range") +
  theme_classic(base_size = 12)

ggsave(file.path(outdir, "Fig_Cortisol_vs_Baseline_by_Year.png"),
       p_cort_baseline, width = 3 * n_distinct(dat_baseline$year), height = 8, dpi = 300)

p_t3_baseline <- ggplot(
  dat_baseline %>% filter(is.finite(t3)),
  aes(julian, t3, colour = survey)
) +
  geom_rect(aes(xmin = -Inf, xmax = Inf,
                ymin = t3_lo, ymax = t3_hi),
            fill = "grey80", alpha = 0.3, inherit.aes = FALSE) +
  geom_hline(aes(yintercept = t3_hat), linewidth = 1) +
  geom_point(alpha = 0.6, size = 2) +
  facet_grid(preg_group ~ year) +
  scale_y_log10(breaks = c(10, 100, 200, 1300)) +
  scale_colour_manual(values = survey_colours, na.value = "grey50") +
  labs(title = "Observed T3 relative to model-based baseline",
       x = "Julian day", y = "T3 (log scale)", colour = "Survey/Range") +
  theme_classic(base_size = 12)

ggsave(file.path(outdir, "Fig_T3_vs_Baseline_by_Year.png"),
       p_t3_baseline, width = 3 * n_distinct(dat_baseline$year), height = 8, dpi = 300)

message("Part 1 complete — outputs in: ", outdir)
