# Dual-hormone monitoring in boreal woodland caribou, northern Ontario

Analysis code and de-identified data for:

> de Almeida Monteiro Melo Ferraz M, Mastromonaco GF, McFarlane S.
> *Winter climate and anthropogenic disturbance jointly structure glucocorticoid
> and thyroid hormone variation in boreal woodland caribou (*Rangifer tarandus
> caribou*).

Faecal cortisol (FCM) and faecal triiodothyronine (FT3) metabolites from
1,380 samples (1,375 with GPS coordinates, and therefore derived covariates, in
the public table) collected across six caribou ranges over four winters
(2021–2024), modelled against winter climate, road proximity, land cover,
mining, settlement proximity and predation risk using generalised additive
models.

---

## ⚠ Collection coordinates are deliberately withheld

Boreal caribou are listed as **Threatened** under Canada's *Species at Risk
Act*, and precise faecal collection coordinates could be used to locate
animals. **Latitude, longitude, projected coordinates, individual genetic
identifiers and exact collection dates are not published in the dataset.**
`scripts/13_build_public_dataset.R` computes every spatially derived
covariate at the true coordinates on a local machine and writes out **only
the derived values**. **Derived covariates are also deliberately coarsened.** Distances to four
known features at full precision would let anyone triangulate each sample's
location, so distances are rounded to the nearest km, wolf suitability to 3
decimals and land cover PC scores to 1 decimal, and rows are in random order.
The published table is therefore analysis-ready:

| Reproducible from this repository | Not reproducible |
|---|---|
| The non-spatial components of the models in scripts `01`–`08`, `11`, `12`: road, mine and settlement distances; wolf habitat suitability; land cover principal components; lagged climate; year, range, reproductive-group and Julian-day terms. Because covariates are coarsened, refits approximate the published statistics rather than reproduce them exactly. | The two-dimensional spatial smooth `s(x, y)`, which requires the coordinates themselves and is part of every primary landscape model. Its fitted surfaces are shown in `outputs/manuscript figures/SupplementaryFigure11.png` and its statistics in `outputs/GAM_summary_*.txt`. |

The distance-computation scripts (`00b`, `03`, `03b`, `13`) are included for
transparency, but cannot be re-run from this repository alone because they
need the coordinates and the large source rasters. They are the code that
*produced* the published covariates.

Coordinates may be made available to researchers for legitimate scientific
purposes on reasonable request to the corresponding author, subject to a
data-sharing agreement.

---

## Repository layout

```
data/
  caribou_hormones_analysis_ready.csv   # 1,375 rows, no coordinates, coarsened covariates
  DATA_DICTIONARY.md                    # every column defined
scripts/
  00_data_processing.R … 13_build_public_dataset.R
outputs/
  GAM_summary_*.txt                     # model summaries as reported
  GAM_AIC_*.txt                         # model-selection tables
  manuscript figures/                   # Figures 1-5, Supplementary 5-13
```

Supplementary Figures 1–4 are the immunoassay validation figures
(parallelism, accuracy, precision). They are not produced by these scripts
and are distributed with the article rather than here.

## Reproducing the reported models

```r
install.packages(c("tidyverse","mgcv","gratia","emmeans","sf","terra",
                   "patchwork","janitor","gstat","ggnewscale"))

dat <- readr::read_csv("data/caribou_hormones_analysis_ready.csv")
dat$range <- factor(dat$range)   # required for the random-effect smooth

# e.g. the winter-road FT3 model, minus its spatial term
mgcv::gam(
  ft3_ng_g ~ s(range, bs = "re") + factor(year) + reproductive_group +
             s(log(dist_winterroad_km), k = 5) + s(julian, k = 5),
  data = dat, family = Gamma(link = "log"), method = "REML")
```

Reported models additionally include `s(x_proj, y_proj, k = 30)`. Because
that term needs coordinates, and covariates are coarsened, refits from this
repository will differ from the published statistics; the published values
are in `outputs/`.

**Why no fitted model objects.** Serialised `mgcv` objects are deliberately
*not* distributed. `gam()` stores its model frame inside the fitted object,
so any model containing `s(x_proj, y_proj)` carries the full coordinate
table — `m$model$x_proj` would recover every collection location in two
lines of R. The `.txt` summaries and rendered partial-effect surfaces are
published instead.

## Road network definitions

Road layers come from the MNR Road Segments provincial dataset (Ontario
GeoHub), 517,836 segments in three classes:

| `NATIONAL_ROAD_CLASS` | Segments | Used as |
|---|---|---|
| Resource / Recreation | 508,220 | all-season |
| Local / Strata | 8,648 | all-season |
| Winter | 968 | winter |

**All-season roads = every non-winter segment (516,868).** An earlier
version of this analysis used `Local / Strata` alone (1.7% of the layer, an
ownership class rather than a road network); that definition was wrong and
has been corrected. Winter roads come from the full provincial file, not a
pre-filtered extract, which was previously missing 427 of the 968 segments.

## Software

R 4.5.2; `mgcv` 1.9-4 for GAMs, `sf` and `terra` for spatial work,
`gratia` for partial effects, `emmeans` for contrasts, `gstat` for
variograms.

## Licence

Code: MIT. Data: CC BY 4.0. Please cite the paper above.

## Contact

Samantha McFarlane — samantha.mcfarlane@ec.gc.ca
