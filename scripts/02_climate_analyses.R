# ============================================================
# 02_climate_analyses.R
#
# PART 2 — CLIMATE INFLUENCE
# Purpose:
#   • Characterise inter-annual climate structure (Jan–Mar)
#   • Link short-lag winter climate to cortisol (FCM) and T3 (FT3),
#     plus medium-lag (6–15 d) and long-lag (16–30 d) windows
#
# SHORT-LAG WINDOW, SPLIT BY HORMONE (updated 2026-07-19):
#   FCM (cortisol) now uses a 2-day short-lag window (1–2 d), based on
#   Ashley et al. 2011 (Gen Comp Endocrinol 172:382-391) — an ACTH-
#   challenge validation study in captive caribou/reindeer showing
#   fecal GC peaks ~8 h post-stimulus and returns to baseline by ~48 h,
#   i.e. a much shorter true excretion/integration window than the
#   previously-used 5-day default.
#   FT3 keeps the original 5-day short-lag window (1–5 d): no
#   caribou/Rangifer-specific FT3 excretion-lag study exists in the
#   literature, and the best available proxy (Pasciu et al., roe deer)
#   suggests FT3 integrates temperature over a longer (~7-day) window,
#   arguing against shortening it to match FCM.
#   See `short2_*` (FCM, 1–2 d) vs `short_*` (FT3, 1–5 d) columns below.
#
# IMPORTANT: climate station (climate_group) is assigned PER SAMPLE in
#   00_data_processing.R, based on the nearest station to each sample's
#   own GPS coordinates — NOT by field campaign (survey_field) and NOT
#   by boundary-assigned range (survey). Rigorous per-sample nearest-
#   station assignment was checked against both simpler alternatives
#   and found necessary: ~38% of Missisa-range samples are
#   geographically closer to James Bay's station (Ogoki Post) than to
#   Missisa's own nominal station (Lansdowne House), and all
#   Pagwachuan samples are closer to Ogoki Post than to Armstrong, the
#   station its field campaign nominally used. See
#   00_data_processing.R for station coordinates and the
#   nearest-station logic (climate_stations, nearest_climate_station()).
#
# DATA FIX (see climatedaily_all_surveys.csv): Lansdowne House (AUT)
#   has no snow-on-ground records for Jan-Mar 2021. 90 daily rows from
#   Geraldton A (~136 km SE, same boreal ecotype; climate ID 6042717)
#   were appended to the CSV under climate_group = "Missisa" for
#   Jan-Mar 2021 only, one row per date, to supply snow-on-ground for
#   that period ONLY. All other Geraldton A fields (temperature,
#   precipitation, humidity, wind) are deliberately blanked in the CSV
#   so calc_lag_mean()'s na.rm = TRUE averaging cannot blend them with
#   Lansdowne House's values.
# ============================================================

source("00_data_processing.R")

outdir <- "caribou_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# LOAD CLIMATE DATA
# ============================================================

infile_climate <- "climatedaily_all_surveys.csv"

clim_all <- read_csv(infile_climate, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(
    date       = as.Date(local_date),
    local_year = as.integer(local_year),
    doy        = yday(date),
    month      = month(date)
  )

clim_winter <- clim_all %>% filter(month %in% c(1, 2, 3))

climate_vars <- c(
  "mean_temperature", "min_temperature", "max_temperature",
  "total_precipitation", "snow_on_ground",
  "min_rel_humidity", "max_rel_humidity", "speed_max_gust"
)

# ============================================================
# SAMPLE DATES x NEAREST CLIMATE STATION (per-sample, from 00)
# ============================================================

sample_survey_dates <- dat %>%
  filter(!is.na(date), !is.na(climate_group)) %>%
  distinct(pcid, survey, survey_field, year, date, climate_group)

# ============================================================
# 1) CLIMATE OVERLAY PLOT (FEB–MAR) BY NEAREST CLIMATE STATION
#    x-axis: calendar date (Feb–Mar)
#    All years shown as thin context lines; sampled years bold
#    Collection dates shown as dots
# ============================================================

clim_lines <- clim_all %>%
  filter(month %in% c(2, 3)) %>%
  select(climate_group, date, local_year,
         mean_temperature, min_temperature, snow_on_ground,
         total_precipitation) %>%
  pivot_longer(
    cols     = c(mean_temperature, min_temperature, snow_on_ground, total_precipitation),
    names_to = "variable", values_to = "value"
  ) %>%
  mutate(
    season_date = as.Date(sprintf("2000-%02d-%02d", month(date), day(date))),
    year_f      = factor(local_year),
    variable    = recode(variable,
      mean_temperature    = "Mean Temp (°C)",
      min_temperature     = "Min Temp (°C)",
      snow_on_ground      = "Snow on Ground (cm)",
      total_precipitation = "Total Precipitation (mm)"
    )
  )

# Sample collection dates (Feb–Mar only) for dot overlay
sample_points <- sample_survey_dates %>%
  filter(month(date) %in% c(2, 3)) %>%
  rename(local_year = year) %>%
  left_join(
    clim_all %>%
      filter(month %in% c(2, 3)) %>%
      select(climate_group, date, local_year,
             mean_temperature, min_temperature, snow_on_ground,
             total_precipitation),
    by = c("climate_group", "date", "local_year")
  ) %>%
  pivot_longer(
    cols     = c(mean_temperature, min_temperature, snow_on_ground, total_precipitation),
    names_to = "variable", values_to = "value"
  ) %>%
  mutate(
    season_date = as.Date(sprintf("2000-%02d-%02d", month(date), day(date))),
    year_f      = factor(local_year),
    variable    = recode(variable,
      mean_temperature    = "Mean Temp (°C)",
      min_temperature     = "Min Temp (°C)",
      snow_on_ground      = "Snow on Ground (cm)",
      total_precipitation = "Total Precipitation (mm)"
    )
  )

# Years actually sampled at each nearest-climate-station group
# (for bolding those lines)
station_years_set <- sample_survey_dates %>% distinct(climate_group, year)

p_clim <- ggplot(clim_lines,
                 aes(season_date, value, group = year_f, colour = year_f)) +
  geom_line(linewidth = 0.35, alpha = 0.3) +
  geom_line(
    data = clim_lines %>%
      semi_join(station_years_set, by = c("climate_group", "local_year" = "year")),
    linewidth = 1.1, alpha = 0.9
  ) +
  geom_point(data = sample_points, size = 2.2, alpha = 0.95) +
  facet_grid(climate_group ~ variable, scales = "free_y") +
  scale_x_date(date_labels = "%b %d", date_breaks = "2 weeks") +
  scale_colour_viridis_d(option = "turbo", name = "Year") +
  labs(x = "Date (Feb–Mar)", y = NULL,
       title = "Winter climate (Feb–Mar) by nearest climate station",
       subtitle = "Bold = years actually sampled at that station; thin = context years; dots = collection dates") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text.x     = element_text(face = "bold", size = 9),
    strip.text.y     = element_text(face = "bold", size = 9, angle = 0),
    axis.text.x      = element_text(angle = 35, hjust = 1, size = 7),
    legend.position  = "bottom"
  )

ggsave(file.path(outdir, "Climate_FebMar_ByCampaign_WithSampleDots.png"),
       p_clim,
       width  = 14,
       height = 2.5 * n_distinct(clim_lines$climate_group),
       dpi    = 300)


# ============================================================
# 2) INTER-ANNUAL CLIMATE COMPARISONS (GAMs + emmeans)
#    Grouping factor: nearest climate station (climate_group);
#    year_f included where a station covers multiple years.
# ============================================================

# Unique nearest-station x year combinations actually sampled
station_years <- sample_survey_dates %>%
  distinct(climate_group, year)

clim_gam <- clim_all %>%
  filter(month %in% c(1, 2, 3)) %>%
  mutate(doy = ifelse(doy == 366, 365, doy)) %>%
  inner_join(station_years,
             by = c("climate_group", "local_year" = "year"),
             relationship = "many-to-many") %>%
  mutate(year_f = factor(local_year))

# GAM per climate variable; primary factor = climate_group (nearest station)
fit_climate_gam <- function(var, data) {
  df <- data %>% filter(is.finite(.data[[var]]))
  if (n_distinct(df$climate_group) < 2) return(NULL)

  has_multi_year <- df %>%
    distinct(climate_group, year_f) %>%
    count(climate_group) %>%
    filter(n > 1) %>% nrow() > 0

  fml <- if (has_multi_year) {
    as.formula(paste0(var, " ~ climate_group + year_f + s(doy, bs = 'cc', k = 10)"))
  } else {
    as.formula(paste0(var, " ~ climate_group + s(doy, bs = 'cc', k = 10)"))
  }
  m   <- gam(fml, data = df, method = "REML")
  emm <- emmeans(m, ~ climate_group)
  list(variable = var, model = m, emmeans = emm,
       contrasts = pairs(emm, adjust = "tukey"))
}

# Inter-annual (year-based) comparison, pooling across whichever
# stations were relevant in each year — primary comparison for
# reporting climate differences by year, per project convention
# (organise by year, not station/range). NOTE: because most stations
# are only relevant in a subset of years, year and station are only
# partially separable here. Report accordingly.
fit_climate_gam_by_year <- function(var, data) {
  df <- data %>% filter(is.finite(.data[[var]]))
  if (n_distinct(df$year_f) < 2) return(NULL)
  fml <- as.formula(paste0(var, " ~ year_f + s(doy, bs = 'cc', k = 10)"))
  m   <- gam(fml, data = df, method = "REML")
  emm <- emmeans(m, ~ year_f)
  list(variable = var, model = m, emmeans = emm,
       contrasts = pairs(emm, adjust = "tukey"))
}

climate_results <- map(climate_vars, fit_climate_gam, data = clim_gam)
names(climate_results) <- climate_vars

climate_results_byyear <- map(climate_vars, fit_climate_gam_by_year, data = clim_gam)
names(climate_results_byyear) <- climate_vars

climate_tukey_byyear <- map_dfr(
  climate_results_byyear,
  function(r) if (!is.null(r)) as.data.frame(r$contrasts) else NULL,
  .id = "climate_variable"
)
write_csv(climate_tukey_byyear,
          file.path(outdir, "Climate_GAM_Tukey_Comparisons_ByYear_JanMar.csv"))

climate_tukey <- map_dfr(
  climate_results,
  function(r) if (!is.null(r)) as.data.frame(r$contrasts) else NULL,
  .id = "climate_variable"
)
write_csv(climate_tukey,
          file.path(outdir, "Climate_GAM_Tukey_Comparisons_ByCampaign_JanMar.csv"))

# Boxplot: variables with significant climate_group (station) differences
sig_vars <- climate_tukey %>% filter(p.value < 0.05) %>%
  distinct(climate_variable) %>% pull()

if (length(sig_vars) > 0) {
  clim_box <- clim_gam %>%
    select(date, climate_group, year_f, all_of(sig_vars)) %>%
    pivot_longer(all_of(sig_vars), names_to = "variable", values_to = "value") %>%
    mutate(variable_label = recode(variable,
             mean_temperature    = "Mean temperature (°C)",
             min_temperature     = "Minimum temperature (°C)",
             max_temperature     = "Maximum temperature (°C)",
             snow_on_ground      = "Snow on ground (cm)",
             min_rel_humidity    = "Minimum relative humidity (%)",
             max_rel_humidity    = "Maximum relative humidity (%)",
             total_precipitation = "Total precipitation (mm)",
             speed_max_gust      = "Max gust speed (km/h)"))

  cld_df <- map_dfr(sig_vars, function(v) {
    r <- climate_results[[v]]
    if (is.null(r)) return(NULL)
    cmp    <- as.data.frame(r$contrasts)
    p_vals <- setNames(cmp$p.value, gsub(" ", "", cmp$contrast))
    ltrs   <- multcompLetters(p_vals)$Letters
    tibble(climate_group = names(ltrs), cld = ltrs, variable = v)
  })

  cld_pos <- clim_box %>%
    group_by(variable, climate_group) %>%
    summarise(y_pos = max(value, na.rm = TRUE), .groups = "drop") %>%
    left_join(cld_df, by = c("variable", "climate_group")) %>%
    left_join(clim_box %>% distinct(variable, variable_label), by = "variable")

  p_climate_box <- ggplot(clim_box,
                           aes(x = year_f, y = value, fill = climate_group)) +
    geom_boxplot(alpha = 0.7, width = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.2, size = 0.6, colour = "grey40") +
    stat_summary(fun = mean, geom = "point", shape = 23,
                 size = 3, fill = "white", colour = "black", stroke = 0.8) +
    geom_text(
      data = cld_pos,
      aes(x = climate_group, y = y_pos, label = cld),
      inherit.aes = FALSE, vjust = -0.5, size = 3.2, fontface = "bold"
    ) +
    facet_grid(variable_label ~ climate_group, scales = "free_y") +
    labs(x = "Year", y = NULL, fill = "Nearest station",
         title = "Winter climate by year and nearest climate station (Jan–Mar)",
         subtitle = "Significant differences only; diamond = mean; letters = Tukey groups") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.text        = element_text(face = "bold"),
          axis.text.x       = element_text(angle = 35, hjust = 1, size = 9))

  ggsave(file.path(outdir, "Climate_Boxplots_by_Year_Campaign.png"),
         p_climate_box,
         width  = max(7, 1.5 * n_distinct(clim_gam$climate_group)),
         height = 3 + 2.5 * length(sig_vars),
         dpi    = 300)
}


# ============================================================
# 3) CLIMATE LAGS — MERGE WITH BIOLOGICAL DATA
#    Lags computed per sample (pcid) x date, using each sample's own
#    nearest climate station (climate_group). Joined back to dat via
#    pcid + year + date (NOT survey_field/survey, since climate_group
#    is now per-sample, not per-campaign or per-range).
# ============================================================

clim_winter_split <- split(clim_winter, clim_winter$climate_group)

calc_lag_mean <- function(target_date, clim_df, lag_start, lag_end) {
  clim_df %>%
    filter(date <  target_date - days(lag_start - 1),
           date >= target_date - days(lag_end)) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
}

climate_lags <- sample_survey_dates %>%
  rowwise() %>%
  mutate(
    short  = list(calc_lag_mean(date, clim_winter_split[[climate_group]], 1,  5)),
    short2 = list(calc_lag_mean(date, clim_winter_split[[climate_group]], 1,  2)),
    medium = list(calc_lag_mean(date, clim_winter_split[[climate_group]], 6,  15)),
    long   = list(calc_lag_mean(date, clim_winter_split[[climate_group]], 16, 30))
  ) %>%
  ungroup() %>%
  pivot_longer(c(short, short2, medium, long), names_to = "period", values_to = "climate") %>%
  unnest(climate) %>%
  pivot_wider(names_from = period,
              values_from = -c(pcid, survey, survey_field, year, date, climate_group, period),
              names_glue  = "{period}_{.value}")
# short_*  = 1-5 day window (used for FT3 / T3 models, unchanged)
# short2_* = 1-2 day window (used for FCM / cortisol models)

write_csv(climate_lags,
          file.path(outdir, "Climate_Lagged_Summaries_byCampaignDate.csv"))

# Join to biological data via pcid + date (each sample's own nearest station)
dat_clim <- dat %>%
  left_join(
    climate_lags %>% select(-survey, -survey_field, -climate_group),
    by = c("pcid", "year", "date")
  )

cat("Rows with climate lag data (FT3, 5-day short lag):", sum(!is.na(dat_clim$short_mean_temperature)), "\n")
cat("Rows without climate lag data (FT3, 5-day short lag):", sum(is.na(dat_clim$short_mean_temperature)), "\n")
cat("Rows with climate lag data (FCM, 2-day short lag):", sum(!is.na(dat_clim$short2_mean_temperature)), "\n")
cat("Rows without climate lag data (FCM, 2-day short lag):", sum(is.na(dat_clim$short2_mean_temperature)), "\n")


# ============================================================
# 4) GAMs — CLIMATE EFFECTS ON CORTISOL AND T3
#    survey (range-assigned) as random intercept; year_f as fixed
#    effect; climate_group used only for computing climate lags.
# ============================================================

dat_gam_clim <- dat_clim %>%
  mutate(survey = factor(survey), year_f = factor(year))

# --- 4A) Cortisol (FCM) — 2-day short-lag window (short2_*) ---

dat_cort_clim <- dat_gam_clim %>%
  filter(is.finite(cortisol),
         is.finite(short2_min_temperature),
         is.finite(short2_snow_on_ground),
         !is.na(julian))

m_cort_clim <- gam(
  cortisol ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(short2_min_temperature, k = 6) +
    s(short2_snow_on_ground,  k = 4) +
    s(julian, bs = "cc", k = 8),
  data   = dat_cort_clim,
  family = Gamma(link = "log"),
  method = "REML"
)

summary(m_cort_clim)
capture.output(summary(m_cort_clim), file = file.path(outdir, "GAM_summary_02_cortisol_climate.txt"))
gam.check(m_cort_clim)

# Model selection: temperature only vs snow only vs both
m_cort_snow <- update(m_cort_clim, . ~ . - s(short2_min_temperature, k = 6))
m_cort_temp <- update(m_cort_clim, . ~ . - s(short2_snow_on_ground,  k = 4))

cat("Cortisol climate model AICs:\n")
aic_cort_clim <- AIC(m_cort_clim, m_cort_snow, m_cort_temp)
print(aic_cort_clim)
capture.output(print(aic_cort_clim), file = file.path(outdir, "GAM_AIC_02_cortisol_climate.txt"))

png(file.path(outdir, "Fig_GAM_Cortisol_ShortLag_Climate.png"), 2000, 1400, res = 300)
plot(m_cort_clim, pages = 1)
dev.off()

# Year pairwise contrasts (survey is RE, not suitable for emmeans)
cort_clim_year_contrasts  <- emmeans(m_cort_clim, ~ year_f)    %>% pairs(adjust = "tukey")
cort_clim_group_contrasts <- emmeans(m_cort_clim, ~ preg_group) %>% pairs(adjust = "tukey")
capture.output(print(cort_clim_year_contrasts),  file = file.path(outdir, "GAM_contrasts_02_cortisol_year.txt"))
capture.output(print(cort_clim_group_contrasts), file = file.path(outdir, "GAM_contrasts_02_cortisol_preggroup.txt"))

# --- 4B) T3 (FT3) — unchanged 5-day short-lag window (short_*) ---
# Kept at 5 days: no caribou/Rangifer-specific FT3 excretion-lag study
# exists; best available proxy (Pasciu et al., roe deer) suggests a
# longer (~7-day) FT3 integration window, so FT3's window is NOT
# shortened to match FCM's new 2-day window (see header note above).

dat_t3_clim <- dat_gam_clim %>%
  filter(is.finite(t3), t3 > 0,
         is.finite(short_snow_on_ground),
         is.finite(short_min_temperature),
         !is.na(julian))

# PRIMARY MODEL for FT3: combined (snow + temperature).
# UPDATED after reducing short_snow_on_ground's k from 8 to 4 (matching
# the harmonized k used throughout the project — see 11/12): the AIC
# comparison FLIPS. With k = 4, the combined model now has the clearly
# lower (better) AIC (see GAM_AIC_02_t3_climate.txt): m_t3_clim
# (combined) AIC = 12317.83 vs m_t3_snowonly AIC = 12368.42 — a decisive
# ΔAIC of ~50.6 favouring the combined model (previously, at k = 8,
# snow-only had won by a thin ΔAIC of ~2.1). Snow on its own is also no
# longer significant in the snow-only model (p = 0.555). The combined
# model is therefore now the primary/reported FT3 climate model,
# matching FCM's convention (m_cort_clim, also combined). The snow-only
# model is still fit below purely for the AIC comparison and as a
# historical reference — it is NOT the reported model.
m_t3_clim <- gam(
  t3 ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(short_min_temperature, k = 4) +
    s(short_snow_on_ground,  k = 4) +
    s(julian, bs = "cc", k = 8),
  data   = dat_t3_clim,
  family = Gamma(link = "log"),
  method = "REML"
)

summary(m_t3_clim)
capture.output(summary(m_t3_clim), file = file.path(outdir, "GAM_summary_02_t3_climate.txt"))
gam.check(m_t3_clim)

# Snow-only model — fit ONLY for the AIC comparison below, NOT used as
# the reported/primary FT3 climate model (see note above).
m_t3_snowonly <- gam(
  t3 ~
    s(survey, bs = "re") +
    year_f +
    preg_group +
    s(short_snow_on_ground, k = 4) +
    s(julian, bs = "cc", k = 8),
  data   = dat_t3_clim,
  family = Gamma(link = "log"),
  method = "REML"
)
capture.output(summary(m_t3_snowonly), file = file.path(outdir, "GAM_summary_02_t3_climate_snowonly_for_comparison.txt"))

cat("T3 climate model AICs (snow-only vs snow+temp) — combined is now primary:\n")
aic_t3_clim <- AIC(m_t3_snowonly, m_t3_clim)
print(aic_t3_clim)
capture.output(print(aic_t3_clim), file = file.path(outdir, "GAM_AIC_02_t3_climate.txt"))

png(file.path(outdir, "Fig_GAM_T3_ShortLag_Climate.png"), 2000, 1400, res = 300)
plot(m_t3_clim, pages = 1)
dev.off()

t3_clim_year_contrasts  <- emmeans(m_t3_clim, ~ year_f)    %>% pairs(adjust = "tukey")
t3_clim_group_contrasts <- emmeans(m_t3_clim, ~ preg_group) %>% pairs(adjust = "tukey")
capture.output(print(t3_clim_year_contrasts),  file = file.path(outdir, "GAM_contrasts_02_t3_year.txt"))
capture.output(print(t3_clim_group_contrasts), file = file.path(outdir, "GAM_contrasts_02_t3_preggroup.txt"))

message("Part 2 complete — outputs in: ", outdir)
