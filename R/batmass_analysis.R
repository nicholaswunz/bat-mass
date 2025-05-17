## SET UP ## -------------------------------------------------------------------
# Load library
#devtools::install_github("daniel1noble/orchaRd", force = TRUE)
pacman::p_load(tidyverse, ggplot2, cowplot, MuMIn, 
               rotl, ape, metafor, orchaRd, performance, 
               rsm, geodata, terra, sp, raster)

# Functions
mytheme <- function() {
  theme_bw() +
    theme(panel.border          = element_rect(fill = NA, colour = "black"), # set border around plot.
          panel.grid.major      = element_blank(), # remove major grid lines
          panel.grid.minor      = element_blank(), # remove minor grid lines
          axis.line             = element_blank(), # remove axis lines
          axis.ticks            = element_line(colour = "black"),
          axis.text             = element_text(size = 10, colour = "black"), # axis text size
          axis.title            = element_text(size = 10), # axis title size
          axis.title.y          = element_text(vjust = 3), # increase distance from the y-axis
          axis.title.x          = element_text(vjust = -1), # increase distance from the x-axis
          panel.background      = element_rect(fill = NA),
          plot.background       = element_rect(fill = NA, color = NA), # remove background colour
          plot.margin           = unit(c(0.2, 0.2, 0.2, 0.2), units = , "cm"),
          legend.background     = element_rect(fill = NA, color = NA), # get rid of legend bg
          legend.box.background = element_rect(fill = NA, color = NA), # get rid of legend panel bg
          strip.text.x          = element_text(size = 10, color = "black", face = "bold"), # for facet plots
          strip.background      = element_rect(fill = NA, color = NA)
    )
} # set up plot theme

## RAW DATA ## -----------------------------------------------------------------
# Load raw data
raw_dat <- read.csv("https://github.com/nicholaswunz/bat-mass/raw/refs/heads/main/data/raw_data.csv") %>%
  dplyr::select(-c(notes, title, ref)) %>%
  dplyr::mutate(lat_abs     = abs(lat),
                lnMass      = log(mass_g),
                pre_post    = factor(pre_post, levels = c("pre", "post")),
                sex         = factor(sex, levels = c("Female", "Male", "Both")),
                roost_type  = factor(roost_type, levels = c("cave", "tree", "broad"))
                )

## CLIMATE DATA ## -------------------------------------------------------------
# Download WorldCLim .tif files () to local folder

# Import all files as a list from your local folder
rastlist <- list.files(path = "/LOCAL_FOLDER", 
                       pattern = '.tif$', all.files = TRUE, full.names = TRUE)

# Stack all layers
allrasters <- raster::stack(rastlist)
bioclim <- allrasters
bioclim <- bioclim[[c(1, 7, 14)]] # extract relevent variables
names(bioclim) <- c("MAST", "rain_season_cv", "T_season")

# Extract coordinates from raw_dat
xy <- sp::SpatialPoints(cbind(raw_dat$lon, raw_dat$lat), 
                        proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs"))

# Extract BioClim values from coordinates
values <- data.frame(terra::extract(bioclim, xy)) %>%
  dplyr::mutate(T_season_sd    = T_season / 100) # Temperature seasonality (standard deviation *100)

# Load GPP data
NIRv_sd_rast     <- raster::raster("https://github.com/nicholaswunz/bat-mass/raw/refs/heads/main/files/NIRv_sd_rast.tif")
NIRv_sd_rast_0.5 <- terra::aggregate(NIRv_sd_rast, 10) # resolution to 0.5
GPP_values       <- data.frame(raster::extract(NIRv_sd_rast_0.5, xy)) %>% rename(GPP_sd = "raster..extract.NIRv_sd_rast_0.5..xy.")

# Combine climate data with raw_dat
clim_dat <- cbind(raw_dat, values, GPP_values)

## DATA SUMMARY ## -------------------------------------------------------------
clim_dat %>%
  group_by(family) %>%
  summarise(n = length(unique(species)))

## PHYLOGENY RECONSTRUCTION ## -------------------------------------------------
sp_all <- sort(unique(as.character(raw_dat$species_OTL)))  # generate list of species (as character format)
taxa   <- rotl::tnrs_match_names(names = sp_all)  # match taxonomic names to the OTL

# check if species list match OT identifier 
#taxa[taxa$approximate_match == TRUE,] # none so far

# retrieving phylogenetic relationships among taxa in the form of a trimmed
# sub-tree
tree <- rotl::tol_induced_subtree(ott_ids = rotl::ott_id(taxa), label_format = "name")

# Change OTL species name to match raw_dat
tree$tip.label[tree$tip.label == "mrcaott307132ott864246"] <- "Myotis_lucifugus"
tree$tip.label[tree$tip.label == "mrcaott62482ott62486"] <- "Artibeus_planirostris"

# Compute branch lengths
tree <- ape::compute.brlen(tree, method = "Grafen", power = 1)
tree <- ape::multi2di(tree, random = TRUE)  # use a randomization approach to deal with polytomies

# Create the phylogenetic correlation matrix from our phylogenetic tree.
Acov <- ape::vcv.phylo(tree, corr = FALSE)

## 3 INITIAL MODELS ## ------------------------------------------------------------------
# Check for multicollinearity. VIF >5 is a problem
car::vif(lm(lnRR ~ MAST + T_season_sd + rain_season_cv + GPP_sd + mean_year_study + pre_post + lnMass + sex + roost_type + days_dif, data = clim_dat))
# no values above 5

MAST_mod <- metafor::rma.mv(lnRR ~ MAST * rain_season_cv + MAST * mean_year_study + lnMass + pre_post * sex + roost_type + days_dif,
                            V = lnRR_vi,
                            random = list(~1 | species, 
                                          ~1 | species_OTL,
                                          ~1 | pub_ID), 
                            R = list(species_OTL = Acov),
                            test = "t", dfs = "contain",
                            data = clim_dat %>%
                              dplyr::mutate(species_OTL = stringr::str_replace_all(species_OTL, " ", "_"))
)

season_mod <- metafor::rma.mv(lnRR ~ T_season_sd * rain_season_cv + MAST * mean_year_study + lnMass + pre_post * sex + roost_type + days_dif,
                              V = lnRR_vi,
                              random = list(~1 | species, 
                                            ~1 | species_OTL,
                                            ~1 | pub_ID), 
                              R = list(species_OTL = Acov),
                              test = "t", dfs = "contain",
                              data = clim_dat %>%
                                dplyr::mutate(species_OTL = stringr::str_replace_all(species_OTL, " ", "_"))
)

GPP_mod <- metafor::rma.mv(lnRR ~ MAST * GPP_sd + MAST * mean_year_study + lnMass + pre_post * sex + roost_type + days_dif,
                           V = lnRR_vi,
                           random = list(~1 | species, 
                                         ~1 | species_OTL,
                                         ~1 | pub_ID), 
                           R = list(species_OTL = Acov),
                           test = "t", dfs = "contain",
                           data = clim_dat %>%
                             dplyr::mutate(species_OTL = stringr::str_replace_all(species_OTL, " ", "_"))
)

performance::compare_performance(MAST_mod, season_mod, GPP_mod, rank = TRUE) # Compare models

## MODEL SELECTION ## ----------------------------------------------------------
eval(metafor:::.MuMIn)
options(na.action = "na.fail")  # required for dredge to run

mumin_dat <- clim_dat %>% dplyr::mutate(species_OTL = stringr::str_replace_all(species_OTL, " ", "_")) 
mumin_dat <- mumin_dat[!apply(mumin_dat[, c("lnRR_vi", "MAST", "T_season_sd", "rain_season_cv", 
                                            "mean_year_study", "lnMass","pre_post", "sex",
                                            "roost_type", "days_dif")],
                              1, anyNA), ]  # remove rows with any NA value on one of these columns

# Run multi-level model
mumin_model <- metafor::rma.mv(lnRR, lnRR_vi, mods = ~ MAST * rain_season_cv + MAST * mean_year_study + lnMass + pre_post * sex + roost_type + days_dif,
                               random = list(~1 | species, 
                                             ~1 | species_OTL,
                                             ~1 | pub_ID), 
                               R = list(species_OTL = Acov),
                               test = "t", dfs = "contain",
                               data = mumin_dat)

candidate_models <- MuMIn::dredge(mumin_model, trace = 2) # generate all possible combinations of moderators 

options(na.action = "na.omit")  # set back to default
subset(candidate_models)  # display all models within 2 values of AICc
subset(candidate_models, delta <= 6)  # display all models within 2 values of AICc
MuMIn::sw(candidate_models) # display weight of each moderators

## FINAL MODEL ## --------------------------------------------------------------
final_mod <- metafor::rma.mv(lnRR, lnRR_vi, mods = ~ MAST * rain_season_cv + MAST * mean_year_study + pre_post + sex,
                             random = list(~1 | species, 
                                           ~1 | species_OTL,
                                           ~1 | pub_ID), 
                             R = list(species_OTL = Acov),
                             test = "t", dfs = "contain",
                             data = clim_dat %>%
                               dplyr::mutate(species_OTL = stringr::str_replace_all(species_OTL, " ", "_"))
)

summary(final_mod)
final_mod$sigma2[2] / sum(final_mod$sigma2) # Phylogenetic signal
orchaRd::r2_ml(final_mod)

## MAST-SEX INT MODEL ## -----------------------------------------------------------
sex_pre_mod <- metafor::rma.mv(lnRR, lnRR_vi, mods = ~ MAST * sex,
                                random = list(~1 | species, 
                                              ~1 | species_OTL,
                                              ~1 | pub_ID), 
                                R = list(species_OTL = Acov),
                                test = "t", dfs = "contain",
                                data = clim_dat %>%
                                  dplyr::filter(sex != "Both" & pre_post == "pre") %>%
                                  dplyr::mutate(species_OTL = stringr::str_replace_all(species_OTL, " ", "_"))
)

sex_post_mod <- metafor::rma.mv(lnRR, lnRR_vi, mods = ~ MAST * sex,
                             random = list(~1 | species, 
                                           ~1 | species_OTL,
                                           ~1 | pub_ID), 
                             R = list(species_OTL = Acov),
                             test = "t", dfs = "contain",
                             data = clim_dat %>%
                             dplyr::filter(sex != "Both" & pre_post == "post") %>%
                             dplyr::mutate(species_OTL = stringr::str_replace_all(species_OTL, " ", "_"))
)

summary(sex_pre_mod)
summary(sex_post_mod)

## FIGURE PRODUCTION ## --------------------------------------------------------
# Figures were cleaned up for publication using Adobe Illustrator.

# Fig 1
clim_mod <- lm(rel_diff ~  MAST * rain_season_cv, data = clim_dat)
summary(clim_mod)

par(mfrow = c(1,2))
contour(clim_mod, rain_season_cv ~ MAST, 
        image = T, 
        img.col = colorspace::sequential_hcl(12, palette = "YlOrRd", rev = T),
        labcex = 1,
        xlabs = c("Mean annual surface temperature (°C)", "Rainfall seasonality (CV)"))
persp(clim_mod, rain_season_cv ~ MAST, 
      zlab = "Mb (%)",
      col = colorspace::sequential_hcl(12, palette = "YlOrRd", rev = T), 
      contours = "colors",
      xlabs = c("Mean annual surface temperature (°C)", "Rainfall seasonality (CV)"))
par(mfrow = c(1,1))

# Fig 2
# Sex
sex_data <- orchaRd::mod_results(final_mod, mod = "sex", group = "species", weights = "prop")$mod_table %>%
  mutate(estimate =  (exp(estimate) - 1) * 100,
         lowerCL = (exp(lowerCL) - 1) * 100,
         upperCL = (exp(upperCL) - 1) * 100,
  )

sex_plot <- ggplot(sex_data, aes(x = name, y = estimate)) +
  ggforce::geom_sina(data = clim_dat, aes(x = sex, y = rel_diff), 
                     colour = "black", alpha = 0.1) +
  geom_linerange(aes(ymin = lowerCL, 
                     ymax = upperCL), color = "black") + 
  geom_point(size = 3, shape = 21, colour = "black", fill = "black") +
  labs(y = expression(Delta*italic("M")["b"]*" (%)"), x = NULL) +
  mytheme() 


# Year
year_MAST_data <- orchaRd::mod_results(final_mod, mod = "mean_year_study", group = "species", weights = "prop",
                                       by = "MAST", 
                                       at = list("MAST" = quantile(clim_dat$MAST, probs = c(0.25, 0.5, 0.75))))$mod_table %>%
  dplyr::mutate(estimate =  (exp(estimate) - 1) * 100,
                lowerCL = (exp(lowerCL) - 1) * 100,
                upperCL = (exp(upperCL) - 1) * 100, 
                condition  = case_when(
                  condition <9 ~ "9 °C",
                  condition >9 & condition <19  ~ "13 °C",
                  condition >19 ~ "19 °C"),
                condition  = factor(condition, levels = c("9 °C", "13 °C", "19 °C"))
  )

year_plot <- ggplot(year_MAST_data, aes(x = moderator)) +
  geom_point(data = clim_dat, aes(x = mean_year_study, y = rel_diff), 
             colour = "black", alpha = 0.1) +
  geom_line(aes(y = estimate, linetype = condition, colour = condition)) +
  labs(y = expression(Delta*italic("M")["b"]*" (%)"), x = "Study year", fill = "MAST", colour = "MAST") +
  scale_color_manual(values = c("#4a6fe3", "black", "#d33f6a")) +
  scale_fill_manual(values = c("#4a6fe3", "black", "#d33f6a")) +
  mytheme()

cowplot::plot_grid(sex_plot, year_plot + theme(legend.position = "bottom"), 
                   align = "h", rel_widths = c(0.7, 1), 
                   ncol = 2, labels = c("a", "b"))

# Fig 3
# MAST * sex
sex_post_data <- orchaRd::mod_results(sex_post_mod, mod = "MAST", group = "species", weights = "prop", 
                                 by = "sex")$mod_table %>%
  dplyr::mutate(estimate =  (exp(estimate) - 1) * 100,
                lowerCL = (exp(lowerCL) - 1) * 100,
                upperCL = (exp(upperCL) - 1) * 100
                )

ggplot(sex_post_data, aes(x = moderator)) +
  geom_point(data = clim_dat, aes(x = MAST, y = rel_diff), 
             colour = "grey") +
  geom_ribbon(aes(ymin = lowerCL, ymax = upperCL, linetype = condition), alpha = 0.1) +
  geom_line(aes(y = estimate, linetype = condition)) +
  labs(y = expression(Delta*italic("M")["b"]*" (%)"), x = "MAST  (°C)") +
  mytheme()

# Fig 4
# Bat distribution
world  <- map_data("world")
bat_sr <- readRDS(file.path(data_path, "bat_all_sr.rds"))
bioclim_ag <- raster::projectRaster(bioclim, crs = crs(bat_sr), res = 0.25)

bat_ag           <- raster::resample(bat_sr, bioclim_ag)
bioclim_ag <- raster::stack(bioclim_ag) # convert SpatRaster to RasterStack
bioclim_bat_rast <- raster::stack(bioclim_ag, bat_ag) # merge raster files together

# Categorise mass loss
bioclim_bat_df  <- raster::as.data.frame(bioclim_bat_rast, xy = T) %>%
  dplyr::rename(species_n = layer) %>%
  dplyr::mutate(intrcpt_est    = 0.4924,
                intrcpt_se     = 0.0443,
                MAST_est       = -0.0181,
                MAST_se        = 0.0015,
                rain_est       = -0.004203,
                rain_se        = 0.000501,
                int_est        = 0.00019,
                int_se         = 0.000024,
                lnRR           = intrcpt_est + (MAST * MAST_est) + (rain_season_cv * rain_est) + (MAST * rain_season_cv * int_est),
                delta_Mb       = (exp(lnRR) - 1) * 100,
                lnRR_low       = intrcpt_est + (MAST * MAST_est) + (rain_season_cv * rain_est) + (MAST * rain_season_cv * (int_est - int_se)),
                lnRR_hig       = intrcpt_est + (MAST * MAST_est) + (rain_season_cv * rain_est) + (MAST * rain_season_cv * (int_est + int_se)),
                delta_Mb_low   = (exp(lnRR_low) - 1) * 100,
                delta_Mb_hig   = (exp(lnRR_hig) - 1) * 100,
                se             = abs(delta_Mb_hig - delta_Mb_low)
  ) %>%
  na.omit(species_n) %>%
  # Recode
  dplyr::mutate(category = case_when(
    delta_Mb <= 0 ~ '<1',                       
    delta_Mb >0 & delta_Mb < 10 ~ '1-10', 
    delta_Mb >= 10 & delta_Mb < 20 ~ '10-20',      
    delta_Mb >= 20 & delta_Mb < 30 ~ '20-30',
    delta_Mb >= 30 & delta_Mb < 40 ~ '30-40', 
    delta_Mb >= 40 & delta_Mb < 50 ~ '40-50',  
    delta_Mb >=50  ~ '>50'                    
  )) %>% 
  # Convert to ordered factor
  dplyr::mutate(category = factor(category,
                                  levels = c('<1', '1-10', '10-20', '20-30',
                                             '30-40', '40-50', '>50'),
                                  ordered = TRUE))

# Fig 4a
map_plot <- bioclim_bat_df %>%
  ggplot() +
  geom_map(data = world %>% dplyr::filter(region != "Antarctica"), map = world, aes(long, lat, map_id = region),
           fill = "lightgrey", colour = NA) +
  geom_raster(aes(y = y, x = x, fill = category)) +
  geom_point(data = clim_dat, aes(x = lon, y = lat), colour = "black", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  colorspace::scale_fill_discrete_sequential(palette = "YlOrRd", rev = T) +
  colorspace::scale_colour_continuous_sequential(palette = "YlOrRd", rev = T) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, fill = expression(Delta*italic("M")["b"]*" (%)")) +
  theme_void() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 10),
        legend.title = element_text(size = 8), 
        legend.position = "bottom") +
  coord_fixed(ratio = 1)

# Fig 4b
se_plot <- bioclim_bat_df %>%
  ggplot() +
  geom_map(data = world %>% dplyr::filter(region != "Antarctica"), map = world, aes(long, lat, map_id = region),
           fill = "lightgrey", colour = NA) +
  geom_raster(aes(y = y, x = x, fill = se)) +
  colorspace::scale_fill_continuous_sequential(palette = "Blues 3", rev = T) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, fill = "se") +
  theme_void() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 10),
        legend.title = element_text(size = 8), 
        legend.position = "bottom") +
  coord_fixed(ratio = 1)


bioclim_df <- raster::as.data.frame(bioclim, xy = T) %>%
  drop_na() %>%
  dplyr::filter(y >-60)

# Fig 4c
temp_map <- ggplot() +
  geom_map(data = world %>% dplyr::filter(region != "Antarctica"), map = world, aes(long, lat, map_id = region),
           fill = NA, colour = "lightgrey") +
  geom_raster(data = bioclim_df, aes(y = y, x = x, fill = MAST)) +
  #colorspace::scale_fill_continuous_diverging(palette = "Blue Red 3", mid = 20, rev = F) +
  colorspace::scale_fill_continuous_sequential(palette = "Purple-Blue", rev = F) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "bottom") +
  coord_fixed(ratio = 1)

# Fig 4d
rain_map <- ggplot() +
  geom_map(data = world %>% dplyr::filter(region != "Antarctica"), map = world, aes(long, lat, map_id = region),
           fill = NA, colour = "lightgrey") +
  geom_raster(data = bioclim_df, aes(y = y, x = x, fill = rain_season_cv)) +
  colorspace::scale_fill_continuous_sequential(palette = "Purple-Blue", rev = T) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "bottom") +
  coord_fixed(ratio = 1)

# Fig 4e
cat_plot <- data.frame(bioclim_bat_df %>% 
                         dplyr::group_by(category) %>% 
                         dplyr::summarise(delta_Mb   = length(delta_Mb[!is.na(delta_Mb)]))) %>%
  dplyr::mutate(category = as.factor(category),
                freq = delta_Mb / sum(delta_Mb) * 100) %>%
  ggplot(aes(x = category, y = freq, fill = category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(freq, 1)), vjust = -1, size = 3) +
  colorspace::scale_fill_discrete_sequential(palette = "YlOrRd", rev = T) +
  ylim(0, 50) +
  labs(y = "% occupancy", x = expression(Delta*italic("M")["b"]*" (%)")) +
  mytheme()

## SUPPLEMENTARY FIGURES ## ----------------------------------------------------
# Fig S4
MAST_tseason_plot <- bioclim_bat_df %>%
  mutate(MAST = round(MAST, digits = 1),
         T_season = round(T_season / 100, digits = 1),) %>%
  group_by(MAST, T_season) %>%
  summarise(species_n = median(species_n)) %>%
  data.frame() %>%
  ggplot() +
  geom_hex(aes(y = T_season, x = MAST), bins = 150, fill = "lightgrey") +
  geom_point(data = clim_dat, aes(x = MAST, y = T_season / 100), shape = 21, fill = "white", size = 2) +
  ylab("Temperature seasonality (σ)") + xlab("Mean annual surface temperature (°C)") +
  mytheme()

cor(clim_dat$MAST, clim_dat$T_season_sd)

MAST_rainseason_plot <- bioclim_bat_df %>%
  mutate(MAST = round(MAST, digits = 1),
         rain_season_cv = round(rain_season_cv, digits = 2)) %>%
  group_by(MAST, rain_season_cv) %>%
  summarise(species_n = median(species_n)) %>%
  data.frame() %>%
  ggplot() +
  geom_hex(aes(y = rain_season_cv, x = MAST), bins = 150, fill = "lightgrey") +
  geom_point(data = clim_dat, aes(x = MAST, y = rain_season_cv), shape = 21, fill = "white", size = 2) +
  ylab("Rainfall seasonality (CV)") + xlab("Mean annual surface temperature (°C)") +
  mytheme()

t_rainseason_plot <- bioclim_bat_df %>%
  mutate(T_season = round(T_season / 100, digits = 1),
         rain_season_cv = round(rain_season_cv, digits = 2)) %>%
  group_by(T_season, rain_season_cv) %>%
  summarise(species_n = median(species_n)) %>%
  data.frame() %>%
  ggplot() +
  geom_hex(aes(y = rain_season_cv, x = T_season), bins = 150, fill = "lightgrey") +
  geom_point(data = clim_dat, aes(x = T_season / 100, y = rain_season_cv), shape = 21, fill = "white", size = 2) +
  ylab("Rainfall seasonality (CV)") + xlab("Temperature seasonality (σ)") +
  mytheme()

cowplot::plot_grid(MAST_tseason_plot, MAST_rainseason_plot, t_rainseason_plot, 
                   labels = c("a", "b", "c"),  align = 'v', axis = 'l')

## Reviewer suggestions ## -----------------------------------------------------
# Run final model with Northern American bats removed (potential influence of WNS). 
final_mod_NA_removed <- metafor::rma.mv(lnRR, lnRR_vi, mods = ~ MAST * rain_season_cv + MAST * mean_year_study + pre_post + sex,
                             random = list(~1 | species, 
                                           ~1 | species_OTL,
                                           ~1 | pub_ID), 
                             R = list(species_OTL = Acov),
                             test = "t", dfs = "contain",
                             data = clim_dat %>%
                               dplyr::mutate(species_OTL = stringr::str_replace_all(species_OTL, " ", "_")) %>%
                               dplyr::filter(!country %in% c("USA", "Canada")) # Remove NA
)

summary(final_mod_NA_removed)
