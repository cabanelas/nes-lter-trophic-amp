################################################################################
## Script:  
## Project: NES-LTER Trophic Amplification - Forage Fish
##          Cross-Site LTER Pelagic Synthesis Working Group 
## Data:    NEFSC Bottom Trawl Survey (Fall & Spring)
##          https://www.fisheries.noaa.gov/inport/item/22557
##          Data Downloaded 09-MAR-2026
## Author:  Alexandra C. Cabanelas Bermudez
## Created: August 2025  |  Updated: March 2026
##
## Purpose: Prepares forage fish biomass (CPUE) from NEFSC bottom trawl 
##          survey for trophic amplification analysis.        
##          Applies vessel/gear correction factors and produces log-transformed
##          annual means and 5-year running means by region and season.              
##          TS1: 1978–1987                           
##          TS2: 1998–present (matched to chl-a) 
##       1) Apply vessel/gear corrections (DCF, GCF, VCF, Bigelow rhoW)
##       2) log10(x + (min/2)) for each station; x = summed forage fish CPUE
##       3) Average across stations per region/season/year
##       4) 5-year running mean (k=5, index-aware)
##       5) Compute SD of time series
##
## Inputs (data/raw/):
##   - 22560_NEFSCFallFisheriesIndependentBottomTrawlData/
##       22560_UNION_FSCS_SVCAT.csv
##       22560_UNION_FSCS_SVSTA.csv
##       Fall_SVDBS_SupportTables/SVDBS_SVSPECIES_LIST.csv
##   - 22561_NEFSCSpringFisheriesIndependentBottomTrawlData/
##       22561_UNION_FSCS_SVCAT.csv
##       22561_UNION_FSCS_SVSTA.csv
##   - correction_factors/NEFSC_conversion_factors.csv   (OceanAdapt/GitHub)
##  https://github.com/pinskylab/OceanAdapt/blob/master/data_raw/NEFSC_conversion_factors.csv
##   - correction_factors/bigelow_corrections.csv        (Miller et al. 2010)
##   - neus_strata.csv
##  https://github.com/pinskylab/OceanAdapt/blob/master/data_raw/neus_strata.csv
##
## Outputs (data/processed/):
##   - 
################################################################################

## ------------------------------------------ ##
#            Packages
## ------------------------------------------ ##
library(tidyverse)  # v2.0.0
library(runner)     # running mean; v0.4.3
library(janitor)
library(data.table)
library(sf)

## ------------------------------------------ ##
#            Data  
## ------------------------------------------ ##
## --- Fall ---
fall_catch <- read_csv(here::here("raw",
                                  "22560_NEFSCFallFisheriesIndependentBottomTrawlData",
                                  "22560_UNION_FSCS_SVCAT.csv"),
              # read as character to avoid floating point precision loss (id ~ 2e17)
                       col_types = cols(ID = col_character())) %>%
  clean_names()

fall_meta <- read_csv(here::here("raw",
                                 "22560_NEFSCFallFisheriesIndependentBottomTrawlData",
                                 "22560_UNION_FSCS_SVSTA.csv"),
                      col_types = cols(ID = col_character())) %>%
  clean_names()

## --- Spring ---
spring_catch <- read_csv(here::here("raw",
                                    "22561_NEFSCSpringFisheriesIndependentBottomTrawlData",
                                    "22561_UNION_FSCS_SVCAT.csv"),
                         col_types = cols(ID = col_character())) %>%
  clean_names()

spring_meta <- read_csv(here::here("raw",
                                   "22561_NEFSCSpringFisheriesIndependentBottomTrawlData",
                                   "22561_UNION_FSCS_SVSTA.csv"),
                        col_types = cols(ID = col_character())) %>%
  clean_names()

## --- spp codes --- 
spp_codes <- read_csv(here::here("raw",
                                 "22560_NEFSCFallFisheriesIndependentBottomTrawlData",
                                 "Fall_SVDBS_SupportTables",
                                 "SVDBS_SVSPECIES_LIST.csv")) %>%
  clean_names()

## --- Conversion Factors ---
# Source: NEFSC_conversion_factors.csv (OceanAdapt GitHub)
nefsc_cf_full <- read_csv(here::here("raw",
                                     "NEFSC_conversion_factors.csv")) %>%
  clean_names() 

## --- Bigelow / Pisces (rhoW) ---
# Source: Miller et al. 2010, NEFSC Ref Doc 10-05, Tables 56-58
bigelow_cf <- read_csv(here::here("raw",
                                  "Miller2010_Bigelow_calibration_factors.csv")) %>%
  clean_names() 

## --- NEUS Strata ---
nes_strata <- read_csv(here::here("raw",
                                  "neus_strata.csv"))

## ------------------------------------------ ##
#            Species codes
## ------------------------------------------ ##
## --- Spp Codes ---
# based on SVDBS_SVSPECIES_LIST.csv
# Atlantic herring    == Clupea harengus      == 032
# Atlantic mackerel   == Scomber scombrus     == 121 
# Atlantic butterfish == Peprilus spp         == 131
#                        Ammodytes spp      
# American sand lance == Ammodytes americanus == 734
# Northern sand lance == Ammodytes dubius     == 181

# svspp codes for target forage fish species
forage_spp <- spp_codes %>%
  filter(svspp %in% c("032", "121", "131", "181", "734")) %>%
  mutate(target_taxa = case_when(
    svspp == "032"              ~ "Atlantic herring",
    svspp == "121"              ~ "Atlantic mackerel",
    svspp == "131"              ~ "Butterfish",
    svspp %in% c("181", "734")  ~ "Sand lance"
  )) #### CHECK IF OKAY TO GROUP SAND LANCE

forage_codes <- forage_spp$svspp

## ------------------------------------------ ##
#            Conversion Factors
## ------------------------------------------ ##
## --- DCF / GCF / VCF ---
# DCF = Door Conversion Factor          applied: year < 1985
# GCF = General Conversion Factor       applied: spring, year 1973-1981
# VCF = Vessel Conversion Factor        applied: vessel == "DE" (Delaware II)
nefsc_cf <- nefsc_cf_full %>%
  mutate(svspp = as.character(svspp)) %>%
  # only keep species relevant to forage fish
  filter(svspp %in% as.character(as.integer(forage_codes))) %>%
  select(svspp, dcf_wt, gcf_wt, vcf_wt) %>%
  mutate(svspp = str_pad(svspp, width = 3, pad = "0"))  # "32" -> "032"

## --- Bigelow / Pisces (rhoW) ---
# Bigelow (HB) and Pisces (PC) vessel corrections post-2008
# Applied to: svvessel %in% c("HB", "PC")
# Season-specific where available; combined estimate where not
# Atlantic mackerel and N. sand lance only appear in Table 56 (combined seasons)
# meaning there weren't enough season-specific observations to estimate separate spring/fall factors
# Ammodytes dubius has no correction factor in this document
# svspp 734 not in Miller 2010 - using 181 value as proxy

# Tables 56, 57, 58
# Columns 5&6 
# ρ_W = biomass-based, total survey estimate; ratio of the total survey-wide biomass estimate between vessels
# SE (ρ_W)
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
# attached and catch summed across sexes.

# We don't need sex-specific catch, so we sum across catchsex within each
# haul x species combination.

# need to add metadata from (SVSTA) to catch df with location, depth, vessel, etc

prep_survey <- function(catch, meta, season_label) {
  catch %>%
    # keep only target forage fish species
    filter(svspp %in% forage_codes) %>%
    # join haul-level metadata (location, vessel, env conditions)
    # id = Concatenation of Cruise6, Stratum, Tow and Station values
    left_join(
      meta %>%
        select(id, svvessel, est_year, est_month,
               decdeg_beglat, decdeg_beglon, avgdepth,
               surftemp, surfsalin, bottemp, botsalin),
      by = "id"
    ) %>%
    rename(
      year   = est_year,
      month  = est_month,
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
    # this collapses e.g. male + female + unknown into one row per haul x svspp
    group_by(id, svspp, year, month, lat, lon, depth, stratum,
             svvessel, season_survey, surftemp, surfsalin, bottemp, botsalin) %>%
    summarise(wtcpue = sum(wtcpue, na.rm = TRUE), .groups = "drop")
}

fall_df   <- prep_survey(fall_catch,   fall_meta,   "Fall")
spring_df <- prep_survey(spring_catch, spring_meta, "Spring")

# verify: should have 4-5 species per survey, no duplicate haul x svspp rows
fall_df %>% count(svspp)
fall_df %>% count(id, svspp) %>% filter(n > 1)

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
#   svspp 734 (A. dubius) has no published rhoW -- no correction applied.

# corrections are species-specific.
# Values == 1.0 in the correction factor files = no correction needed.

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
# svspp 734 (A. dubius) has no Bigelow CF - no correction applied
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
ecomap1 <- st_read(here::here("raw", "EcomonStrata_v4.shp"))
ecomap2 <- st_read(here::here("raw", "EcomonStrata_v4b.shp"))

ecomap <- rbind(ecomap1, ecomap2)
st_crs(ecomap) <- 4326
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
  clean_names()
spring_df_region <- assign_region(spring_df, ecomap_valid) %>%
  clean_names()

# check NAs
sum(is.na(fall_df_region$region))
sum(is.na(spring_df_region$region))

# check region assignment
table(fall_df_region$region)
table(spring_df_region$region)

# plot
ggplot(fall_df_region, aes(x = lon, 
                           y = lat, color = region)) +
  geom_point(size = 2) +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "none")
##

### OR
# svmstrata <- read_csv(here::here("raw",
#                                  "22560_NEFSCFallFisheriesIndependentBottomTrawlData",
#                                  "Fall_SVDBS_SupportTables",
#                                  "SVDBS_SVMSTRATA.csv")) %>%
#   clean_names()
# 
# ggplot(svmstrata, aes(x = -midlon/100, y = midlat/100, color = stratum_name)) +
#   geom_point(size = 2) +
#   coord_fixed() +
#   theme_bw() +
#   theme(legend.position = "none")
# 
# svm_pts <- svmstrata %>%
#   mutate(
#     lon = -midlon / 100,
#     lat =  midlat / 100
#   ) %>%
#   filter(!is.na(lon), !is.na(lat)) %>%
#   st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
# 
# svm_region_lookup <- svm_pts %>%
#   bind_cols(
#     ecomap_valid %>%
#       st_drop_geometry() %>%
#       slice(st_nearest_feature(svm_pts, ecomap_valid))
#   ) %>%
#   st_drop_geometry() %>%
#   select(stratum, stratum_name, Region) %>%
#   distinct()
# 
# svmstrata_region <- svm_pts %>%
#   bind_cols(
#     ecomap_valid %>%
#       st_drop_geometry() %>%
#       slice(st_nearest_feature(svm_pts, ecomap_valid))
#   ) %>%
#   st_drop_geometry()
# 
# ggplot(svmstrata_region, aes(x = lon, y = lat, color = Region)) +
#   geom_point(size = 2) +
#   coord_fixed() +
#   theme_bw() +
#   theme(legend.position = "none")
# 
# fall_df1 <- fall_df %>%
#   left_join(svm_region_lookup %>% select(stratum, Region), by = "stratum")
# 
# ggplot(fall_df1, aes(x = lon, y = lat, color = Region)) +
#   geom_point(size = 2) +
#   coord_fixed() +
#   theme_bw() +
#   theme(legend.position = "none")








# ff <- ff %>%
#   left_join(strata_lookup %>% select(stratum, region),
#             by = "stratum") %>%
#   filter(region %in% c("MAB", "GB", "GOM"))

fall_df_region <- fall_df_region %>%
  mutate(
    region = case_when(
      region == "SNE" ~ "MAB",
      TRUE ~ region
    )
  ) %>%
  filter(region %in% c("MAB", "GB", "GOM"))

spring_df_region <- spring_df_region %>%
  mutate(
    region = case_when(
      region == "SNE" ~ "MAB",
      TRUE ~ region
    )
  ) %>%
  filter(region %in% c("MAB", "GB", "GOM"))

table(fall_df_region$region, useNA = "ifany")
table(spring_df_region$region, useNA = "ifany")