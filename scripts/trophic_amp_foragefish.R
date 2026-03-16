################################################################################
## Script:  trophic_amp_foragefish.R
## Project: NES-LTER Trophic Amplification - Forage Fish
##          Cross-Site LTER Pelagic Synthesis Working Group 
## Data:    NOAA NEFSC Bottom Trawl Survey (Fall & Spring)
##          https://www.fisheries.noaa.gov/inport/item/22557
##          Data Downloaded 09-MAR-2026
## Author:  Alexandra C. Cabanelas Bermudez
## Created: August 2025  |  Updated: March 2026
##
## Purpose: Prepares forage fish biomass from NEFSC bottom trawl survey for 
##          trophic amplification analysis.        
##          Applies vessel/gear correction factors and produces log-transformed
##          annual means and 5-year running means by region and season for
##          for 1998–2023 (to match chl-a/sat records)
##       1) Apply vessel/gear corrections (DCF, GCF, VCF, Bigelow rhoW)
##       2) log10(x + (min/2)) for each station; x = summed forage fish CPUE
##       3) Average across stations per region/season/year
##       4) 5-year running mean (k=5, index-aware)
##       5) Compute SD of time series
##     Troph amp code starts line 447
##
## Inputs (data/raw/):
##   - 22560_NEFSCFallFisheriesIndependentBottomTrawlData/
##       22560_UNION_FSCS_SVCAT.csv                 #catch data
##       22560_UNION_FSCS_SVSTA.csv                 #metadata
##       Fall_SVDBS_SupportTables/SVDBS_SVSPECIES_LIST.csv
##       Fall_SVDBS_SupportTables/SVDBS_SVMSTRATA.csv
##   - 22561_NEFSCSpringFisheriesIndependentBottomTrawlData/
##       22561_UNION_FSCS_SVCAT.csv
##       22561_UNION_FSCS_SVSTA.csv
##   - NEFSC_conversion_factors.csv   (OceanAdapt/GitHub)
##     https://github.com/pinskylab/OceanAdapt/blob/master/data_raw/NEFSC_conversion_factors.csv
##   - Miller2010_Bigelow_calibration_factors.csv   (Miller et al. 2010)
##     https://repository.library.noaa.gov/view/noaa/3726
##   - EcomonStrata_v4.shp  (from EcoMon)
##   - EcomonStrata_v4b.shp (from EcoMon)
##
## Outputs (data/output/):
##   - trophamp_ff_1998_2023_full = all station rows with sd_ff and run_mean_ff joined on
##   - trophamp_ff_1998_2023_sum  = one row per region/season/year (annual means only)
##   - trophamp_ff_1998_2023_sd   = sd_ff per region/season
################################################################################

## ------------------------------------------ ##
#            Packages
## ------------------------------------------ ##
library(tidyverse)  
library(runner)     # running mean
library(janitor)
library(data.table)
library(sf)

## ------------------------------------------ ##
#            Data  
## ------------------------------------------ ##
## --- Fall ---
fall_catch <- read_csv(here::here("data", "raw",
                       "22560_NEFSCFallFisheriesIndependentBottomTrawlData",
                       "22560_UNION_FSCS_SVCAT.csv"),
              # read as character to avoid floating point precision loss (id ~ 2e17)
                       col_types = cols(ID = col_character())) %>%
  clean_names()

fall_meta <- read_csv(here::here("data", "raw",
                      "22560_NEFSCFallFisheriesIndependentBottomTrawlData",
                      "22560_UNION_FSCS_SVSTA.csv"),
                      col_types = cols(ID = col_character())) %>%
  clean_names()

## --- Spring ---
spring_catch <- read_csv(here::here("data", "raw",
                         "22561_NEFSCSpringFisheriesIndependentBottomTrawlData",
                         "22561_UNION_FSCS_SVCAT.csv"),
                         col_types = cols(ID = col_character())) %>%
  clean_names()

spring_meta <- read_csv(here::here("data", "raw",
                        "22561_NEFSCSpringFisheriesIndependentBottomTrawlData",
                        "22561_UNION_FSCS_SVSTA.csv"),
                        col_types = cols(ID = col_character())) %>%
  clean_names() %>%
  # fix typo
  rename(end_est_towdate = being_est_towdate)

## --- spp codes --- 
spp_codes <- read_csv(here::here("data", "raw",
                      "22560_NEFSCFallFisheriesIndependentBottomTrawlData",
                      "Fall_SVDBS_SupportTables",
                      "SVDBS_SVSPECIES_LIST.csv")) %>%
  clean_names()

## --- region info --- 
svmstrata <- read_csv(here::here("data", "raw",
                      "22560_NEFSCFallFisheriesIndependentBottomTrawlData",
                      "Fall_SVDBS_SupportTables",
                      "SVDBS_SVMSTRATA.csv")) %>%
  clean_names() %>%
  select(stratum, stratum_name) %>%
  distinct()

## --- Conversion Factors ---
# Source: NEFSC_conversion_factors.csv (OceanAdapt GitHub)
#https://github.com/pinskylab/OceanAdapt/raw/master/data_raw/NEFSC_conversion_factors.csv
nefsc_cf_full <- read_csv(here::here("data", "raw",
                          "NEFSC_conversion_factors.csv")) %>%
  clean_names() 

## --- Bigelow / Pisces (rhoW) ---
# Source: Miller et al. 2010, NEFSC Ref Doc 10-05, Tables 56-58
#https://repository.library.noaa.gov/view/noaa/3726
bigelow_cf <- read_csv(here::here("data", "raw",
                       "Miller2010_Bigelow_calibration_factors.csv")) %>%
  clean_names() 

## ------------------------------------------ ##
#            Species codes
## ------------------------------------------ ##
## --- Spp Codes ---
#     SVDBS_SVSPECIES_LIST.csv
# Atlantic herring    == Clupea harengus      == 032
# Atlantic mackerel   == Scomber scombrus     == 121 
# Atlantic butterfish == Peprilus spp         == 131
#                        Ammodytes spp      
# Northern sand lance == Ammodytes dubius     == 181
# American sand lance == Ammodytes americanus == 734

# svspp codes for target forage fish species
forage_spp <- spp_codes %>%
  filter(svspp %in% c("032", "121", "131", "181", "734")) %>%
  mutate(target_taxa = case_when(
    svspp == "032"              ~ "Atlantic herring",
    svspp == "121"              ~ "Atlantic mackerel",
    svspp == "131"              ~ "Butterfish",
    svspp %in% c("181", "734")  ~ "Sand lance"
  )) 

forage_codes <- forage_spp$svspp

## ------------------------------------------ ##
#            Conversion Factors
## ------------------------------------------ ##
## --- DCF / GCF / VCF ---
# DCF = Door Conversion Factor      applied: year < 1985
# GCF = Gear Conversion Factor      applied: spring, year 1973-1981
# VCF = Vessel Conversion Factor    applied: vessel == "DE" (Delaware II)
nefsc_cf <- nefsc_cf_full %>%
  mutate(svspp = as.character(svspp)) %>%
  # only keep species relevant to forage fish
  filter(svspp %in% as.character(as.integer(forage_codes))) %>%
  select(svspp, dcf_wt, gcf_wt, vcf_wt) %>%
  mutate(svspp = str_pad(svspp, width = 3, pad = "0"))  # 32 -> 032

## --- Bigelow / Pisces (rhoW) ---
# Bigelow (HB) and Pisces (PC) vessel corrections post-2008
# Season-specific where available; combined estimate where not
# Atlantic mackerel and N. sand lance only appear in Table 56 (combined seasons)
#there weren't enough season-specific observations to estimate separate spring/fall factors
# Ammodytes americanus (svspp 734) has no correction factor in this document
# using 181 value as proxy

# Tables 56, 57, 58
# Columns 5 & 6
# ρ_W = ratio of the total survey-wide biomass estimate between vessels
bigelow_cf <- bigelow_cf %>%
  #mutate(svspp = as.character(svspp)) %>%
  select(svspp, season, rho_w) %>%
  mutate(season = str_to_title(season))

## ------------------------------------------ ##
#            Tidy & Join
## ------------------------------------------ ##
# join catch to haul metadata, filter to forage spp only
# sum across sexes (catchsex) for same species at same haul
# one row per haul (id) x species (svspp), with station-level metadata
# catch summed across sexes (don't need sex-specific catch)
# need to add metadata from (SVSTA) to catch df with location, depth, vessel, etc

prep_survey <- function(catch, meta, season_label) {
  catch %>%
    # keep only target forage fish species
    filter(svspp %in% forage_codes) %>%
    # join haul-level metadata 
    # id = Concatenation of Cruise6, Stratum, Tow and Station values
    left_join(
      meta %>%
        select(id, svvessel, est_year, est_month, est_day,
               begin_est_towdate, end_est_towdate,
               decdeg_beglat, decdeg_beglon, avgdepth,
               surftemp, surfsalin, bottemp, botsalin),
      by = "id"
    ) %>%
    rename(
      year   = est_year,
      month  = est_month,
      day    = est_day,
      lat    = decdeg_beglat,
      lon    = decdeg_beglon,
      depth  = avgdepth,
      wtcpue = expcatchwt  # expanded catch weight (kg per tow, standardized)
    ) %>%
    mutate(
      season_survey = season_label,   # survey season (Fall/Spring)
      wtcpue        = as.numeric(wtcpue),
      svvessel      = as.character(svvessel)
    ) %>%
    # sum across sexes for same species at same haul
    # this collapses male + female + unknown into one row per haul x svspp
    group_by(id, svspp, year, month, day, begin_est_towdate, end_est_towdate,
             lat, lon, depth, stratum,
             svvessel, season_survey, surftemp, surfsalin, bottemp, botsalin) %>%
    summarise(wtcpue = sum(wtcpue, na.rm = TRUE), .groups = "drop")
}

fall_df   <- prep_survey(fall_catch,   fall_meta,   "Fall")
spring_df <- prep_survey(spring_catch, spring_meta, "Spring")

fall_df %>% count(svspp)
fall_df %>% count(id, svspp) %>% filter(n > 1) # no duplicate haul x svspp rows

## ------------------------------------------ ##
#            Vessel & Gear Corrections
## ------------------------------------------ ##
## Goal: standardize wtcpue across time so that catch from different vessels
# and gear configurations are comparable to a common baseline
# Albatross IV, the primary survey vessel 1963-2008

## Step 1: DCF - Door conversion (pre-1985, all vessels)
#   The survey switched otter trawl door types in 1985.
#   For affected species, catches before 1985 are multiplied by DCF
#   to make them comparable to post-1985 catches.
#   Source: NEFSC_conversion_factors.csv, column dcf_wt

## Step 2: GCF - Gear conversion (spring only, 1973-1981)
#   A gear change affected spring surveys during this period.
#   Catches are divided by GCF to standardize.
#   Only applied to spring survey, years 1973-1981.
#   Source: NEFSC_conversion_factors.csv, column gcf_wt

## Step 3: VCF - Vessel conversion (Delaware II only)
#   The Delaware II ("DE") had different catchability than Albatross IV.
#   Catches are multiplied by VCF to standardize.
#   Source: NEFSC_conversion_factors.csv, column vcf_wt

## Step 4: rhoW - Bigelow/Pisces correction (HB & PC vessels)
#   In 2008 the survey switched from Albatross IV to Henry B. Bigelow (HB).
#   Pisces (PC) uses the same corrections as Bigelow.
#   The Bigelow catches MORE of most species due to net design differences,
#   so we DIVIDE by rhoW to bring Bigelow catches back down to
#   Albatross IV equivalent levels.
#   rhoW is the biomass-based calibration factor (ρ_W) from Tables 56-58
#   of Miller et al. 2010 (NEFSC Ref Doc 10-05).
#   Season-specific values used where available (herring, butterfish);
#   combined season estimate used where not (mackerel, N. sand lance).
#   svspp 734 (A. americanus) has no published rhoW -- no correction applied.

# corrections are species-specific
# Values == 1.0 = no correction needed

apply_corrections <- function(df, nefsc_cf, bigelow_cf) {
  # join DCF/GCF/VCF factors onto data by species
  df <- df %>%
    left_join(nefsc_cf, by = "svspp")
  
  # Step 1: DCF = multiply pre-1985 catches by door conversion factor
  # corrects for net door change
  df <- df %>%
    mutate(wtcpue = case_when(
      year < 1985 & dcf_wt != 1 ~ wtcpue * dcf_wt,
      TRUE                      ~ wtcpue
    ))
  
  # Step 2: GCF (spring only, 1973-1981) = divide catches by gear conversion factor
  # corrects for gear change in spring surveys only
  df <- df %>%
    mutate(wtcpue = case_when(
      season_survey == "Spring" & year > 1972 & year < 1982 & gcf_wt != 1 ~ wtcpue / gcf_wt,
      TRUE                                                                ~ wtcpue
    ))
  
  # Step 3: VCF - multiply Delaware II catches by vessel conversion factor
  # corrects for catchability difference vs Albatross IV
  df <- df %>%
    mutate(wtcpue = case_when(
      svvessel == "DE" & vcf_wt != 1 ~ wtcpue * vcf_wt,
      TRUE                           ~ wtcpue
    ))
  
  df <- df %>% select(-dcf_wt, -gcf_wt, -vcf_wt)
  
  # Step 4: Bigelow/Pisces rhoW
  # join season-specific rhoW where available, "both" as fallback
  # svspp 734 (A. americanus) has no Bigelow CF - no correction applied
  rho_specific <- bigelow_cf %>% 
    filter(season != "Both")
  rho_both     <- bigelow_cf %>% 
    filter(season == "Both") %>% select(svspp, rho_w)
  
  df <- df %>%
    left_join(rho_specific %>% rename(survey_season = season),
              by = c("svspp", "season_survey" = "survey_season")) %>%
    left_join(rho_both %>% rename(rho_w_both = rho_w), by = "svspp") %>%
    mutate(
      rho_final = coalesce(rho_w, rho_w_both),
      wtcpue = case_when(
        svvessel %in% c("HB", "PC") & !is.na(rho_final) ~ wtcpue / rho_final,
        TRUE                                             ~ wtcpue
      )
    ) %>%
    select(-rho_w, -rho_w_both, -rho_final)
  
  return(df)
}

fall_df   <- apply_corrections(fall_df,   nefsc_cf, bigelow_cf)
spring_df <- apply_corrections(spring_df, nefsc_cf, bigelow_cf)

fall_df %>% 
  filter(svvessel %in% c("HB", "PC")) %>% 
  summarise(min_year = min(year, na.rm = TRUE), 
            max_year = max(year, na.rm = TRUE))
# should be post-2008 only

## ------------------------------------------ ##
#            Assign Regions
## ------------------------------------------ ##
## --- read in EcoMon shapefiles to get region names ---
ecomap1 <- st_read(here::here("data", "raw", "EcomonStrata_v4.shp"))
ecomap2 <- st_read(here::here("data", "raw", "EcomonStrata_v4b.shp"))

ecomap <- rbind(ecomap1, ecomap2)
st_crs(ecomap) <- 4326 #set crs
ecomap_valid <- st_make_valid(ecomap)

# function to assign region to any df with lat/lon columns
assign_region <- function(df, ecomap_valid) {
  df_sf <- df %>%
    filter(!is.na(lat), !is.na(lon)) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326)
  
  df_sf %>%
    bind_cols(
      ecomap_valid %>%
        st_drop_geometry() %>%
        slice(st_nearest_feature(df_sf, ecomap_valid))
    ) %>%
    mutate(
      lon = st_coordinates(.)[,1],
      lat = st_coordinates(.)[,2]
    ) %>%
    st_drop_geometry()
}

fall_df_region   <- assign_region(fall_df,   ecomap_valid) %>%
  clean_names() %>%
  # add region info from svmstrata
  left_join(svmstrata, by = "stratum")
spring_df_region <- assign_region(spring_df, ecomap_valid) %>%
  clean_names() %>%
  # add region info from svmstrata
  left_join(svmstrata, by = "stratum")

# check NAs
sum(is.na(fall_df_region$region))
sum(is.na(spring_df_region$region))

# check region assignment
table(fall_df_region$region)
table(spring_df_region$region)

fall_df_region %>%
  count(region, stratum_name) %>%
  arrange(region, stratum_name) %>%
  print(n = 50)

## --- plot --- 
ggplot(fall_df_region, aes(x = lon, 
                           y = lat, color = region)) +
  geom_point(size = 2) +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "none")

## --- filter by regions of interest --- 
fall_df_region <- fall_df_region %>%
  mutate(
    region = case_when(
      region == "SNE" ~ "MAB",
      TRUE ~ region
    )
  ) %>%
  filter(region %in% c("MAB", "GB", "GOM")) %>%
  select(-c(name, numof_poly, numof_sta, area))

spring_df_region <- spring_df_region %>%
  mutate(
    region = case_when(
      region == "SNE" ~ "MAB",
      TRUE ~ region
    )
  ) %>%
  filter(region %in% c("MAB", "GB", "GOM")) %>%
  select(-c(name, numof_poly,numof_sta, area))

table(fall_df_region$region, useNA = "ifany")
table(spring_df_region$region, useNA = "ifany")

## ------------------------------------------ ##
#            Combine & Assign Calendar Season
## ------------------------------------------ ##
ff <- bind_rows(fall_df_region, spring_df_region) %>%
  left_join(forage_spp %>% select(svspp, target_taxa), by = "svspp") %>%
  mutate(season = case_when(
    month %in% 3:5   ~ "Spring",
    month %in% 6:8   ~ "Summer",
    month %in% 9:11  ~ "Fall",
    TRUE             ~ "Winter"
  ))

## ------------------------------------------ ##
#            Total Forage Fish Abundance
## ------------------------------------------ ##
# sum cpue across species per haul

## --- Sum ALL 4 species into single forage fish index per haul ---  
unique(ff$target_taxa)
ff <- ff %>%
  group_by(id, year, month, day, begin_est_towdate, end_est_towdate,
           season, season_survey, lat, lon, depth, stratum, stratum_name, 
           region, svvessel, surftemp, surfsalin, bottemp, botsalin) %>%
  summarise(wtcpue = sum(wtcpue, na.rm = TRUE), .groups = "drop")

## --- plot raw ---
ggplot(ff, aes(x = year, y = wtcpue)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~season, scales = "free") +
  theme_bw() +
  labs(title = "Raw forage fish CPUE (corrected)")

ggplot(ff, aes(x = year, y = wtcpue)) +
  geom_point(alpha = 0.3) +
  facet_grid(season ~ region, scales = "free") +
  theme_bw() +
  labs(title = "Raw forage fish CPUE (corrected)")

ggplot() +
  geom_sf(data = ecomap_valid, aes(fill = Region), 
          alpha = 0.3, color = "gray40") +
  borders("world", colour = "gray50", fill = "gray90") +
  geom_point(data = ff, aes(x = lon, y = lat, color = region), 
             alpha = 0.1, size = 2) +
  coord_sf(xlim = c(-80, -60), ylim = c(30, 50)) +
  theme_bw() 

## ------------------------------------------ ##
#    Split into Two Time Series
## ------------------------------------------ ##

## --- Time series 1 = 1977-1987 ---
ts1_ff <- ff %>% 
  filter(between(year, 1977, 1987))

## --- Time series 2 = 1998-2023 ---
# matched to chla time series
ts2_ff <- ff %>% 
  filter(between(year, 1998, 2023))

## ------------------------------------------ ##
#    1) Log Transformation
## ------------------------------------------ ##
# log10(x + min_nonzero/2) per region/season/species group
# find the minimum non-zero value for each region and season

log_transform_ff <- function(df) {
  df %>%
    group_by(region, season) %>% 
    mutate(
      min_nonzero = min(wtcpue[wtcpue > 0], na.rm = TRUE),
      log10_ff    = log10(wtcpue + min_nonzero / 2)
    ) %>%
    ungroup()
}

ts1_ff <- log_transform_ff(ts1_ff) # TS1: 1977–1987 
ts2_ff <- log_transform_ff(ts2_ff) # TS2: 1998–2023

## --- plot ---
ggplot(ts2_ff, aes(x = year, y = log10_ff)) +
  geom_point() +
  facet_grid(season ~ region, scales = "free") +
  theme_bw()

## ------------------------------------------ ##
#    2) Annual Means by Region, Season, Species
## ------------------------------------------ ##
# average across stations for a cruise or year/season

# two versions, needed downstream:
# ts_ff     - station-level, log_mean_ff broadcast as repeated column
# ts_ff_sum - one row per region/season/year/species, used for running mean
## ------------------------------------------ ##

## --- Full station-level data with annual mean column attached ---
annual_means_ff <- function(df) {
  df %>%
    group_by(region, season, year) %>% 
    mutate(
      log_mean_ff   = mean(log10_ff,  na.rm = TRUE),
      mean_sfc_temp = mean(surftemp,  na.rm = TRUE),
      mean_sfc_salt = mean(surfsalin, na.rm = TRUE),
      mean_btm_temp = mean(bottemp,   na.rm = TRUE),
      mean_btm_salt = mean(botsalin,  na.rm = TRUE)
    ) %>%
    ungroup()
}

ts1_ff <- annual_means_ff(ts1_ff) # TS1: 1977–1987 
ts2_ff <- annual_means_ff(ts2_ff) # TS2: 1998–2023

## --- Summarize (one row per region/season/year) ---
# used for running mean
## --- Time series 1 = 1977-1987 ---
ts1_ff_sum <- ts1_ff %>%
  group_by(region, season, year) %>% 
  summarize(
    log_mean_ff   = mean(log_mean_ff,   na.rm = TRUE),
    mean_sfc_temp = mean(mean_sfc_temp, na.rm = TRUE),
    mean_sfc_salt = mean(mean_sfc_salt, na.rm = TRUE),
    mean_btm_temp = mean(mean_btm_temp, na.rm = TRUE),
    mean_btm_salt = mean(mean_btm_salt, na.rm = TRUE),
    .groups = "drop"
  )

## --- Time series 2 = 1998-2023 ---
ts2_ff_sum <- ts2_ff %>%
  group_by(region, season, year) %>% 
  summarize(
    log_mean_ff   = mean(log_mean_ff,   na.rm = TRUE),
    mean_sfc_temp = mean(mean_sfc_temp, na.rm = TRUE),
    mean_sfc_salt = mean(mean_sfc_salt, na.rm = TRUE),
    mean_btm_temp = mean(mean_btm_temp, na.rm = TRUE),
    mean_btm_salt = mean(mean_btm_salt, na.rm = TRUE),
    .groups = "drop"
  )

## --- plot ---
ggplot(ts2_ff_sum, aes(x = year, y = log_mean_ff)) +
  geom_point() +
  facet_grid(season~region) 

## ------------------------------------------ ##
#    3) 5-Year Running Mean
## ------------------------------------------ ##
# take a running mean with timespan ~ longest lived taxon
# smoothing over life span of organisms 
# k = 5 
# run_mean_zp uses all available years up to k (right-aligned, na_rm = TRUE)

## --- Time series 1 = 1977-1987 ---
ts1_rm_ff <- ts1_ff_sum %>%
  arrange(region, season, year) %>% 
  group_by(region, season) %>% 
  mutate(run_mean_ff = runner::mean_run(log_mean_ff, k = 5,
                                        idx = year, na_rm = TRUE)) %>%
  ungroup()

## --- Time series 2 = 1998-2023 ---
ts2_rm_ff <- ts2_ff_sum %>%
  arrange(region, season, year) %>% 
  group_by(region, season) %>% 
  mutate(run_mean_ff = runner::mean_run(log_mean_ff, k = 5,
                                        idx = year, na_rm  = TRUE)) %>%
  ungroup()

## --- plot ---
ggplot(ts2_rm_ff, aes(x = year, y = run_mean_ff)) +
  geom_point() +
  facet_grid(season~region) 

## ------------------------------------------ ##
#    4) Standard Deviation
## ------------------------------------------ ##
# SD of 5-yr log annual means across each time series
# = interannual variability metric for trophic amplification test
# does variability increase at higher trophic levels?

## --- Time series 1 = 1977-1987 ---
ts1_ff <- ts1_ff %>%
  group_by(region, season) %>%
  mutate(sd_ff = sd(log_mean_ff, na.rm = TRUE)) %>%
  ungroup()

## --- Time series 2 = 1998-2023 ---
ts2_ff <- ts2_ff %>%
  group_by(region, season) %>%
  mutate(sd_ff = sd(log_mean_ff, na.rm = TRUE)) %>%
  ungroup()

# final SD table
ff_sd <- ts2_ff %>%
  group_by(region, season) %>%
  summarize(sd_ff = sd(log_mean_ff, na.rm = TRUE), .groups = "drop")

## --- plot ---
ggplot(ff_sd, aes(x = season, y = sd_ff, color = region)) +
  geom_point(size = 3) +
  geom_line(aes(group = region)) +
  theme_bw() +
  labs(title = "SD FF 1998-2023")

## ------------------------------------------ ##
#    5) Final Output Tables
## ------------------------------------------ ##
# FULL: all station rows + run_mean_ff joined on
# SUM:  one row per region/season/year/species

## --- Time series 1 --
ts1_ff_full <- ts1_ff %>%
  left_join(ts1_rm_ff %>% select(region, season, year, run_mean_ff),
            by = c("region", "season", "year"))
# ts1_ff_full <- ts1_ff %>%
#   left_join(ts1_rm_ff %>% select(region, season, year, run_mean_ff),
#             by = c("region", "season", "year")) %>%
#   left_join(ts1_sd, by = c("region", "season"))

ts1_ff_sum <- ts1_ff_full %>%
  distinct(region, season, year, .keep_all = TRUE) %>%
  select(region, season, year, min_nonzero, log10_ff, log_mean_ff,
         mean_sfc_temp, mean_sfc_salt, mean_btm_temp, mean_btm_salt,
         sd_ff, run_mean_ff)

## --- Time series 2 ---
ts2_ff_full <- ts2_ff %>%
  left_join(ts2_rm_ff %>% select(region, season, year, run_mean_ff),
            by = c("region", "season", "year"))
# ts2_ff_full <- ts2_ff %>%
#   left_join(ts2_rm_ff %>% select(region, season, year, run_mean_ff),
#             by = c("region", "season", "year")) %>%
#   left_join(ts2_sd, by = c("region", "season"))

ts2_ff_sum <- ts2_ff_full %>%
  distinct(region, season, year, .keep_all = TRUE) %>%
  select(region, season, year, min_nonzero, log10_ff, log_mean_ff,
         mean_sfc_temp, mean_sfc_salt, mean_btm_temp, mean_btm_salt,
         sd_ff, run_mean_ff)

# write.csv(ts1_ff_full, "data/output/trophamp_ff_1977_1987_full.csv",   row.names = FALSE)
# write.csv(ts1_ff_sum,  "data/output/trophamp_ff_1977_1987_sum.csv",    row.names = FALSE)
# write.csv(ts2_ff_full, "data/output/trophamp_ff_1998_2023_full.csv",row.names = FALSE)
# write.csv(ts2_ff_sum,  "data/output/trophamp_ff_1998_2023_sum.csv", row.names = FALSE)
# write.csv(ff_sd, "data/output/trophamp_ff_1998_2023_sd.csv", row.names = FALSE)

## ------------------------------------------ ##
#    Plots
## ------------------------------------------ ##

## --- Annual log means ---
ggplot(ts2_ff_sum, aes(x = year, y = log_mean_ff)) +
  geom_point() +
  geom_line(data = ts2_rm_ff, aes(x = year, y = run_mean_ff),
            color = "red", linewidth = 1.2) +
  facet_grid(season ~ region) +
  scale_color_brewer(palette = "Set2") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12,
                                   color = "black"),
        panel.grid = element_blank()) +
  labs(title = "Forage Fish (1998-2023)",
       y = "log10 mean CPUE", x = NULL)

## --- SD ---
# ggplot(ts2_ff_sum, aes(x = region, y = sd_ff, color = season)) +
#   geom_point(size = 4, position = position_dodge(width = 0.3)) +
#   scale_color_brewer(palette = "Set2") +
#   theme_bw() +
#   labs(title = "Forage Fish Interannual Variability (1998-present)",
#        x = "Region", y = "SD of log10 annual mean CPUE", color = "Season")

## --- anomaly --- 
ts2_anom <- ts2_ff_sum %>%
  group_by(region, season) %>%
  mutate(ff_anom = log_mean_ff - mean(log_mean_ff, na.rm = TRUE)) %>%
  ungroup()

ggplot(ts2_anom, aes(year, ff_anom)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_col() +
  facet_grid(season ~ region) +
  theme_bw() +
  labs(y = "Anomaly in log10 annual mean FF", x = NULL)

ggplot(ts2_anom, aes(x = year, y = ff_anom, color = region)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~season) +
  theme_bw() +
  labs(
    y = "Anomaly in log10 annual mean FF",
    x = NULL,
    color = "Region"
  )

## --- distribution --- 
ggplot(ts2_ff_sum, aes(region, log_mean_ff, fill = season)) +
  geom_boxplot() +
  theme_bw() +
  labs(y = "log10 annual mean FF", x = NULL)

## --- 
ts2_check <- ts2_ff_sum %>%
  arrange(region, season, year) %>%
  group_by(region, season) %>%
  mutate(
    n_in_window = runner::runner(
      x = log_mean_ff,
      k = 5,
      idx = year,
      f = \(x) sum(!is.na(x))
    )
  ) %>%
  ungroup()

ggplot(ts2_check, aes(year, n_in_window)) +
  geom_point() +
  geom_line() +
  facet_grid(season ~ region) +
  theme_bw() +
  labs(y = "Observations in 5-year window", x = NULL)


## --- look at tow duration --- 
# pre-2009: intended duration = 0.5 hr; post-2009: 0.333 hr
# TOWDUR column is in minutes in the raw metadata

fall_meta <- fall_meta %>%
  mutate(haul_dur = as.numeric(towdur) / 60)

fall_meta_filtered <- fall_meta %>%
  filter(
    (est_year <  2009 & between(haul_dur,   0.5   - 0.083, 0.5   + 0.083)) |
      (est_year >= 2009 & between(haul_dur, 0.333 - 0.083, 0.333 + 0.083))
  )

# check how many hauls removed
cat("Removed", nrow(fall_meta) - nrow(fall_meta_filtered), 
    "fall hauls (", 
    round((nrow(fall_meta) - nrow(fall_meta_filtered))/nrow(fall_meta)*100, 1),
    "%)\n")

## --- Spring ---
spring_meta <- spring_meta %>%
  mutate(haul_dur = as.numeric(towdur) / 60)

spring_meta_filtered <- spring_meta %>%
  filter(
    (est_year <  2009 & between(haul_dur, 0.417, 0.583)) |
      (est_year >= 2009 & between(haul_dur, 0.250, 0.416))
  )

cat("Removed", nrow(spring_meta) - nrow(spring_meta_filtered),
    "spring hauls (",
    round((nrow(spring_meta) - nrow(spring_meta_filtered)) / nrow(spring_meta) * 100, 1),
    "%)\n")

## --- Plot function
plot_haul_dur <- function(meta, survey_label) {
  
  hline_pre  <- data.frame(period = factor("pre-2009 (target: 0.5 hr)",
                                           levels = c("pre-2009 (target: 0.5 hr)",
                                                      "post-2009 (target: 0.333 hr)")),
                           yint = c(0.417, 0.5, 0.583))
  
  hline_post <- data.frame(period = factor("post-2009 (target: 0.333 hr)",
                                           levels = c("pre-2009 (target: 0.5 hr)",
                                                      "post-2009 (target: 0.333 hr)")),
                           yint = c(0.250, 0.333, 0.416))
  meta %>%
    filter(!is.na(haul_dur)) %>%
    mutate(
      est_year = as.numeric(est_year),
      period   = factor(
        if_else(est_year < 2009,
                "pre-2009 (target: 0.5 hr)",
                "post-2009 (target: 0.333 hr)"),
        levels = c("pre-2009 (target: 0.5 hr)",
                   "post-2009 (target: 0.333 hr)")
      ),
      status = case_when(
        est_year <  2009 & between(haul_dur, 0.417, 0.583) ~ "Within window",
        est_year >= 2009 & between(haul_dur, 0.250, 0.416) ~ "Within window",
        TRUE                                               ~ "Outside window"
      )
    ) %>%
    filter(!is.na(period)) %>%
    ggplot(aes(x = est_year, y = haul_dur, color = status)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_hline(data = hline_pre,  aes(yintercept = yint),
               linetype = c("dashed", "solid", "dashed"), linewidth = 0.6) +
    geom_hline(data = hline_post, aes(yintercept = yint),
               linetype = c("dashed", "solid", "dashed"), linewidth = 0.6) +
    scale_color_manual(values = c("Within window"  = "steelblue",
                                  "Outside window" = "firebrick")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    scale_x_continuous(breaks = seq(1963, 2023, by = 4)) +
    facet_wrap(~period, scales = "free_x") +
    theme_bw() +
    labs(title    = paste(survey_label, "haul durations"),
         subtitle = "Dashed lines = ±5 min window; solid = target duration",
         x = NULL, y = "Haul duration (hrs)", color = NULL)
}

plot_haul_dur(fall_meta,   "Fall")
plot_haul_dur(spring_meta, "Spring")

plot_haul_dur <- function(meta, survey_label) {
  
  x_min <- 1963
  x_cut <- 2009
  x_max <- 2025
  
  meta %>%
    filter(!is.na(haul_dur)) %>%
    mutate(
      est_year = as.numeric(est_year),
      status = case_when(
        est_year <  2009 & between(haul_dur, 0.417, 0.583) ~ "Within window",
        est_year >= 2009 & between(haul_dur, 0.250, 0.416) ~ "Within window",
        TRUE                                               ~ "Outside window"
      )
    ) %>%
    ggplot(aes(x = est_year, y = haul_dur, color = status)) +
    geom_point(alpha = 0.5, size = 1.5) +
    # vertical cutoff line
    geom_vline(xintercept = x_cut, linetype = "dotted", 
               linewidth = 0.8, color = "gray30") +
    # pre-2009 reference lines (solid = target, dashed = window)
    geom_segment(aes(x = x_min, xend = x_cut, y = 0.500, yend = 0.500),
                 linetype = "solid",  linewidth = 0.6, color = "gray30") +
    geom_segment(aes(x = x_min, xend = x_cut, y = 0.417, yend = 0.417),
                 linetype = "dashed", linewidth = 0.6, color = "gray30") +
    geom_segment(aes(x = x_min, xend = x_cut, y = 0.583, yend = 0.583),
                 linetype = "dashed", linewidth = 0.6, color = "gray30") +
    # post-2009 reference lines
    geom_segment(aes(x = x_cut, xend = x_max, y = 0.333, yend = 0.333),
                 linetype = "solid",  linewidth = 0.6, color = "gray30") +
    geom_segment(aes(x = x_cut, xend = x_max, y = 0.250, yend = 0.250),
                 linetype = "dashed", linewidth = 0.6, color = "gray30") +
    geom_segment(aes(x = x_cut, xend = x_max, y = 0.416, yend = 0.416),
                 linetype = "dashed", linewidth = 0.6, color = "gray30") +
    scale_color_manual(values = c("Within window"  = "steelblue",
                                  "Outside window" = "firebrick")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    scale_x_continuous(breaks = seq(1963, 2025, by = 4)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title    = paste(survey_label, "haul durations"),
         subtitle = "Vertical line = 2009 vessel transition | Dashed = ±5 min window",
         x = NULL, y = "Haul duration (hrs)", color = NULL)
}

plot_haul_dur(fall_meta,   "Fall")
plot_haul_dur(spring_meta, "Spring")

## --- Count table ---
dur_table <- function(meta, survey_label) {
  meta %>%
    #filter(as.numeric(est_year) >= 1998) %>%
    mutate(
      est_year = as.numeric(est_year),
      status   = case_when(
        is.na(haul_dur)                                     ~ "NA",
        est_year <  2009 & between(haul_dur, 0.417, 0.583)  ~ "Within window",
        est_year >= 2009 & between(haul_dur, 0.250, 0.416)  ~ "Within window",
        TRUE                                                ~ "Outside window"
      )
    ) %>%
    count(est_year, status) %>%
    pivot_wider(names_from = status, values_from = n, values_fill = 0) %>%
    mutate(pct_drop  = round(Drop / (Keep + Drop + NA) * 100, 1),
           survey    = survey_label)
}

bind_rows(
  dur_table(fall_meta,   "Fall"),
  dur_table(spring_meta, "Spring")
) %>%
  arrange(survey, est_year) %>%
  print(n = 60)

## --- Outlier ---
outlier_inspect <- function(meta, survey_label) {
  meta %>%
    #filter(as.numeric(est_year) >= 1998) %>%
    mutate(est_year = as.numeric(est_year)) %>%
    filter(haul_dur > 1 | haul_dur < 0.05) %>%
    select(id, est_year, svvessel, stratum, haul_dur,
           decdeg_beglat, decdeg_beglon) %>%
    arrange(desc(haul_dur)) %>%
    mutate(survey = survey_label)
}

bind_rows(
  outlier_inspect(fall_meta,   "Fall"),
  outlier_inspect(spring_meta, "Spring")
) %>%
  print(n = 50)