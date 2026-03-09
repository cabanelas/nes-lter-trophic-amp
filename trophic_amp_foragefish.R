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
#        https://github.com/pinskylab/OceanAdapt/blob/master/data_raw/NEFSC_conversion_factors.csv
##   - correction_factors/bigelow_corrections.csv        (Miller et al. 2010)
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

## ------------------------------------------ ##
#            Data  
## ------------------------------------------ ##
## --- Fall ---
fall_catch <- read_csv(here::here("raw",
                                  "22560_NEFSCFallFisheriesIndependentBottomTrawlData",
                                  "22560_UNION_FSCS_SVCAT.csv")) %>%
  clean_names()

fall_meta <- read_csv(here::here("raw",
                                 "22560_NEFSCFallFisheriesIndependentBottomTrawlData",
                                 "22560_UNION_FSCS_SVSTA.csv")) %>%
  clean_names()

## --- Spring ---
spring_catch <- read_csv(here::here("raw",
                                    "22561_NEFSCSpringFisheriesIndependentBottomTrawlData",
                                    "22561_UNION_FSCS_SVCAT.csv")) %>%
  clean_names()

spring_meta <- read_csv(here::here("raw",
                                   "22561_NEFSCSpringFisheriesIndependentBottomTrawlData",
                                   "22561_UNION_FSCS_SVSTA.csv")) %>%
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
  select(svspp, season, rho_w)

## ------------------------------------------ ##
#            Tidy & Join
## ------------------------------------------ ##
# join catch to haul metadata, filter to forage spp only
# sum across sexes (catchsex) for same species at same haul

prep_survey <- function(catch, meta, season_label) {
  catch %>%
    filter(svspp %in% forage_codes) %>%
    #mutate(svspp = as.character(svspp)) %>%
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
      wtcpue = expcatchwt
    ) %>%
    mutate(
      season_survey = season_label,   # survey season (Fall/Spring)
      wtcpue        = as.numeric(wtcpue),
      svvessel      = as.character(svvessel)
    ) %>%
    # sum across sexes for same species at same haul
    group_by(id, svspp, year, month, lat, lon, depth, stratum,
             svvessel, season_survey, surftemp, surfsalin, bottemp, botsalin) %>%
    summarise(wtcpue = sum(wtcpue, na.rm = TRUE), .groups = "drop")
}

fall_df   <- prep_survey(fall_catch,   fall_meta,   "Fall")
spring_df <- prep_survey(spring_catch, spring_meta, "Spring")
