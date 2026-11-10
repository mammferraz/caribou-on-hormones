# Data dictionary — caribou_hormones_analysis_ready.csv

One row per faecal sample with GPS coordinates (n = 1,375 of 1,380 samples; the five
samples without coordinates have no derived covariates). Collection coordinates are
deliberately omitted; see README.md and the manuscript's Data Availability statement.

**Coarsened to protect animal locations:** distances are rounded to the nearest km,
wolf suitability to 3 decimals and land cover PC scores to 1 decimal; rows are in
random order. Refits from this file therefore differ slightly from published values.

| Column | Units | Description |
|---|---|---|
| sample_id | — | Anonymous sample identifier. Does not encode location or individual. |
| year | — | Sampling year (2021-2024). |
| julian | day | Day of year of collection. |
| range | — | Caribou range, assigned by spatial overlap with published boundaries. |
| reproductive_group | — | Male, Female_nonpregnant or Female_pregnant. |
| fcm_ng_g | ng/g | Faecal cortisol metabolites. |
| ft3_ng_g | ng/g | Faecal triiodothyronine metabolites. |
| fpm_ng_g | ng/g | Faecal progesterone metabolites (used for pregnancy assignment). |
| dist_winterroad_km | km | Distance to nearest winter road (rounded to 1 km). |
| dist_allroad_km | km | Distance to nearest all-season road (all non-winter classes; rounded to 1 km). |
| dist_mine_km | km | Distance to nearest mine (rounded to 1 km). |
| dist_settlement_km | km | Distance to nearest derived settlement centroid (rounded to 1 km). |
| wolf_suitability | 0-1 | Predicted wolf habitat suitability (3 decimals). |
| LC_PC1-LC_PC4 | — | Land cover principal component scores (14 classes, 1-km window; 1 decimal). |
| min_temp_lag2_c | °C | Mean daily minimum temperature, 1-2 d before collection (FCM models). |
| snow_lag2_cm | cm | Mean snow on ground, 1-2 d before collection (FCM models). |
| min_temp_lag5_c | °C | Mean daily minimum temperature, 1-5 d before collection (FT3 models). |
| snow_lag5_cm | cm | Mean snow on ground, 1-5 d before collection (FT3 models). |

## Not included

Latitude, longitude, projected coordinates, individual genetic identifiers and
collection dates are withheld. The two-dimensional spatial smooth s(x, y) used
in the published models therefore cannot be refit from this file. Fitted mgcv
objects are NOT published either: gam() stores its model frame inside the object,
so any model containing s(x_proj, y_proj) would carry the full coordinate table.
The .txt model summaries and the rendered partial-effect surfaces are published
instead, in outputs/.
