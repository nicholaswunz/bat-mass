# Bat seasonal fat use
[![license](https://img.shields.io/badge/license-MIT%20+%20file%20LICENSE-lightgrey.svg)](https://choosealicense.com/)
![Open Source
Love](https://badges.frapsoft.com/os/v2/open-source.svg?v=103)

This repository contains code and data needed to reproduce the article:

**Wu N. C., Villada-Cadavid, T., Welbergen J. A., & Turbill C.** (In submission) Seasonal fattening among bat populations globally: storing energy for survival in a changing world.


**Raw data**
- `raw_dat.csv` - Raw data used for the analysis.

***R*** **code**
- `batmass_analysis.R` - Data cleaning, analysis and figure production.

**External files**
- `wc2.1_2.5m_bio_1.tif` - [WorldClim bio1](https://www.worldclim.org/data/bioclim.html) - Mean annual surface temperature.
- `wc2.1_2.5m_bio_4.tif` - [WorldClim bio4](https://www.worldclim.org/data/bioclim.html) - Temperature seasonality.
- `wc2.1_2.5m_bio_15.tif` - [WorldClim bio15](https://www.worldclim.org/data/bioclim.html) - Rainfall seasonality.
- `NIRv_sd_rast.tif` - Gross primary production seasonality as NIRv.

## Abstract
Seasonality is a fundamental challenge for life on Earth and energy storage prior to colder and drier periods by fattening is a common strategy of animals for survival. Fattening should reflect a trade-off between an expected seasonal energy deficit and the costs of increased body mass, which are particularly important to flying endotherms. We applied a phylogenetic, multilevel meta-analysis to published data on body mass change (Δ*M*<sub>b</sub>), a proxy of fat storage, among bat populations over low productivity periods (winter in temperate regions and dry periods in the tropics) to examine its relationship with global variation in yearly average and seasonality of local climates. We found that Δ*M*<sub>b</sub> increased with decreasing mean annual surface temperature (MAST) but Δ*M*<sub>b</sub> also increased at higher MAST with higher seasonality of rainfall. Seasonal use of body energy reserves by bats is predicted to be widespread in warm, seasonal climates at low latitudes but is poorly studied compared to cold temperate regions. We also found that, in colder climates only, overwinter mass loss was less for females than males, as predicted by the 'thrifty females' hypothesis, and furthermore, that Δ*M*<sub>b</sub> has increased with year of study in warm climates, possibly linked to effects of global climate change on their energetics. Our quantitative synthesis describes the intrinsic and environmental factors determining the extent of seasonal reliance on stored body energy reserves among bat populations, and highlights the global importance of seasonal fattening for survival in this diverse and widespread mammal group.

**Keywords:** Chiroptera, hibernation, life-history, reproduction, torpor, energetics, migration

## Meta-data
Column descriptors for `raw_dat.csv` are as follows:

- pub_ID: Unique identifier of the study.
- pub_year: Year of the publication
- first_author: Last name of the first author.
- source: Where the values were taken (text, figures, tables, raw data). 
- suborder: Suborder of the species.
- family: Family of the species.
- species: Full species name as stated in the study.
- species_OTL: Full species name as stated in the [Open Tree of Life](https://tree.opentreeoflife.org/opentree/argus/opentree15.1@ott93302).
- roost_type: Recorded wintering roost of the study site (cave, tree, broad).
- pd: Whether Pd was detected (yes, blank).
- country: Location of the study site at country level.
- hemisphere: Location of the study site relative to the equator (northern, southern)
- lat: Latitude of the study site (decimal degrees).
- long: Longitude of the study site. (decimal degrees).
- clim_ID: Unique identifier of the study site climate following the [Köppen climate classification](https://en.wikipedia.org/wiki/K%C3%B6ppen_climate_classification).
- study_year: The year/s when the study was conducted.
- mean_year_study: The average year when the study was conducted (for studies with multiple years).
- sex: Sex of the species (Female, Male, Both).
- age: Age of the species (Adult, Juvenile).
- date_pre: The approximate date when the pre-winter body mass was taken.
- date_max: The approximate date when the maximum body mass was taken.
- date_post: The approximate date when the post-winter body mass was taken.
- month_num_pre: The month when the pre-winter body mass was taken (numeric).
- month_num_max: The month when the maximum body mass was taken (numeric).
- month_num_post: The month when the post-winter body mass was taken (numeric).
- days_dif: The number of days between measurements.
- month_pre: The month when the pre-winter body mass was taken (Gregorian).
- month_max: The month when the maximum body mass was taken (Gregorian).
- month_post: The month when the post-winter body mass was taken (Gregorian).
- value_pre_mean: Average body mass value pre-winter (grams).
- value_pre_SD: Standard deviation of the body mass value pre-winter.
- pre_n: Sample size when the pre-winter body mass was taken (numeric).
- value_max_mean: Average body mass value at the start of winter (grams).
- value_max_SD:  Standard deviation of the body mass value at the start of winter.
- max_n: Sample size when the maximum body mass was taken (numeric).
- value_post_mean: Average body mass value post-winter (grams).
- value_post_SD: Standard deviation of the body mass value post-winter.
- post_n: Sample size when the post-winter body mass was taken (numeric).
- pre_post: Categorised name whether change in body mass was pre-winter or post-winter (pre, post).
- mass_g: Average body mass value post-winter (grams).
- abs_diff: The absolute difference between value_max_mean and either value_pre_mean or value_post_mean.
- abs_vi: The absolute variance between value_max_mean and either value_pre_mean or value_post_mean.
- rel_diff: The relative difference between value_max_mean and either value_pre_mean or value_post_mean.
- lnRR: The log response ratio between value_max_mean and either value_pre_mean or value_post_mean.
- lnRR_vi: The log response ratio variance between value_max_mean and either value_pre_mean or value_post_mean. 
- notes: Extra notes regarding the study.
- title: Title of the study.
- ref: Link to the study.

## License
This repository is provided by the authors under the MIT License ([MIT](http://opensource.org/licenses/MIT)).
