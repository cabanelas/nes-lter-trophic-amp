################################################################################
## Script:  
## Project: NES-LTER Trophic Amplification - Zooplankton
##          Cross-Site LTER Pelagic Synthesis Working Group 
## Data:    EcoMon Plankton Survey 
## Author:  Alexandra C. Cabanelas Bermudez
## Created: November 2023  |  Updated: March 2026
##
## Purpose: Prepares zooplankton displacement volume data from the NES for 
##          trophic amplification analysis.        
##          Produces log-transformed annual means and 5-year running means 
##          by region and season for two time periods:              
##          TS1: 1978–1987                           
##          TS2: 1998–present (matched to chl-a)     
##
## Inputs (data/raw/):
##   - EcoMon_plankton_v3_10.csv  (187513.4.4.tar.gz)
## https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:0187513
## Data available for cruises through 2023 - As of 2024_11_06
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

## ------------------------------------------ ##
#            Data  
## ------------------------------------------ ##
#EcoMon_Plankton_Data_v3_10_wStrataMeta.csv
zp_full <- read_csv(here::here("raw",
                               "EcoMon_plankton_v3_10.csv")) %>%
  clean_names()

## ------------------------------------------ ##
#            Tidy Data  
## ------------------------------------------ ##
## --- Add season column --- 
zp <- zp_full %>%
  mutate(season = case_when(
    between(month_num, 3, 5)  ~ "Spring",
    between(month_num, 6, 8)  ~ "Summer",
    between(month_num, 9, 11) ~ "Fall",
    TRUE                      ~ "Winter"
  ))

## --- Filter regions of interest ---
# exclude CC and NS [Region 0]
# rename SNE to MAB to align with phytoplankton time series
# Region 1 = MAB
# Region 2 = SNE (will become MAB)
# Region 3 = GB
# Region 4 = GOM
zp <- zp %>%
  filter(region %in% c("SNE", "MAB", "GOM", "GB")) %>%
  # rename SNE as MAB to align w phyto ts
  mutate(region = case_when(
    region == "SNE" ~ "MAB",
    TRUE            ~ region
  ))

# some observations have zooplankton abundance but volume == NA
# may not have measured volume on all tows
zp_full %>%
  summarize(across(ctyp_10m2:pnepau_10m2, ~ sum(is.na(.))))

## --- Dropping NA volumes ---
zp <- zp %>%
  select(
    cruise_name, station, zoo_gear, ich_gear, date, day, month, month_num, year,
    time, depth, season, region, area, type, lon, lat,
    sfc_temp, sfc_salt, btm_temp, btm_salt,
    volume_1m2   # Zooplankton Displacement Volume (mL per 1m2 surface area)
  ) %>%
  drop_na(volume_1m2)

## ------------------------------------------ ##
#    Diagnostic Plot 
## ------------------------------------------ ##
ggplot(zp, aes(x = year, y = volume_1m2)) +
  geom_point(alpha = 0.3) +
  facet_grid(season ~ region, scales = "free") +
  theme_bw() +
  labs(title = "Raw ZP displacement volume")

## ------------------------------------------ ##
#    Split into Two Time Series
## ------------------------------------------ ##

# TS1: 1978–1987 
ts1_zp <- zp %>%
  filter(between(year, 1978, 1987)) #it starts 1977

# TS2: 1998–present (matched to chla time series)
ts2_zp <- zp %>%
  filter(year >= 1998)

## ------------------------------------------ ##
#    1) Log Transformation
## ------------------------------------------ ##
# find the minimum non-zero value for each region and season

# computed per region/season, so values are not strictly
# comparable across groups on the raw log scale
# may want to switch to a single global offset before final analysis
# if cross-region/season SD comparisons are made

log_transform <- function(df) {
  df %>%
    group_by(region, season) %>%
    mutate(
      min_nonzero = min(volume_1m2[volume_1m2 > 0], na.rm = TRUE),
      log10_zp    = log10(volume_1m2 + min_nonzero / 2)
    ) %>%
    ungroup()
}

ts1_zp <- log_transform(ts1_zp) # TS1: 1978–1987 
ts2_zp <- log_transform(ts2_zp) # TS2: 1998–present

ggplot(ts2_zp, aes(x = year, y = log10_zp)) +
  geom_point() +
  facet_grid(season~region, scales = "free")

## ------------------------------------------ ##
#    2) Annual Means by Region & Season
## ------------------------------------------ ##
# average across stations for a cruise or year/season

# two versions, needed downstream:
# ts_zp (mutate)      - keeps all station-level rows, adds log_mean_zp as a 
#                       repeated column. Needed later to re-attach
#                       station metadata (lat, lon, temp, salt).
# ts_zp_a (summarize) - collapses to one row per region/season/year.
#                       Used for running mean (step 3).
# log_mean_zp values are identical between the two - only shape differs

# --- Full station-level data with annual mean column attached
# station-level rows retained 
# Time series 1 = 1978-1987
#ts1_zp <- ts1_zp %>%
#  group_by(region, season, year) %>%
#  mutate(log_mean_zp = mean(log10_zp, na.rm = TRUE)) %>%
#  ungroup()

# Time series 2 = 1998-now
#ts2_zp <- ts2_zp %>%
#  group_by(region, season, year) %>%
#  mutate(log_mean_zp = mean(log10_zp, na.rm = TRUE)) %>%
#  ungroup()

## --- Summarize (one row per region/season/year) ---
# used for running mean
# Time series 1 = 1978-1987
ts1_zp_a <- ts1_zp %>%
  group_by(region, season, year) %>%
  summarize(log_mean_zp = mean(log10_zp, na.rm = TRUE), .groups = "drop")

# Time series 2 = 1998-now
ts2_zp_a <- ts2_zp %>%
  group_by(region, season, year) %>%
  summarize(log_mean_zp = mean(log10_zp, na.rm = TRUE), .groups = "drop")

# these two plots should look identical
#ggplot(ts2_zp, aes(x = year, y = log_mean_zp)) +
#  geom_point() +
#  facet_grid(season~region) 

ggplot(ts2_zp_a, aes(x = year, y = log_mean_zp)) +
  geom_point() +
  facet_grid(season~region) 

##NEED TO CHECK
## --- Summarize environmental conditions (one row per region/season/year) ---
ts1_env <- ts1_zp %>%
  group_by(region, season, year) %>%
  summarize(
    mean_sfc_temp = mean(sfc_temp, na.rm = TRUE),
    mean_sfc_salt = mean(sfc_salt, na.rm = TRUE),
    mean_btm_temp = mean(btm_temp, na.rm = TRUE),
    mean_btm_salt = mean(btm_salt, na.rm = TRUE),
    .groups = "drop"
  )

ts2_env <- ts2_zp %>%
  group_by(region, season, year) %>%
  summarize(
    mean_sfc_temp = mean(sfc_temp, na.rm = TRUE),
    mean_sfc_salt = mean(sfc_salt, na.rm = TRUE),
    mean_btm_temp = mean(btm_temp, na.rm = TRUE),
    mean_btm_salt = mean(btm_salt, na.rm = TRUE),
    .groups = "drop"
  )

## ------------------------------------------ ##
#    3) 5-Year Running Mean
## ------------------------------------------ ##
# take a running mean with timespan ~ longest lived taxon
# smoothing over life span of organisms 
# k = 5 
# run_mean_zp uses all available years up to k (right-aligned, na_rm = TRUE)

# Time series 1 = 1978-1987
ts1_rm <- ts1_zp_a %>%
  arrange(region, season, year) %>%
  group_by(region, season) %>%
  mutate(run_mean_zp = runner::mean_run(log_mean_zp, na_rm = TRUE, k = 5)) %>%
  ungroup()

# Time series 2 = 1998-now
ts2_rm <- ts2_zp_a %>%
  arrange(region, season, year) %>%
  group_by(region, season) %>%
  mutate(run_mean_zp = runner::mean_run(log_mean_zp, na_rm = TRUE, k = 5)) %>%
  ungroup()

## ------------------------------------------ ##
#    4) Standard Deviation
## ------------------------------------------ ##

# SD of annual log-means across each time series 
# = interannual variability metric for trophic amplification test
# does variability increase at higher trophic levels?
#
# TODO: confirm whether SD should be:
#   (a) pooled across all regions and seasons (current approach keeps separate)
#   (b) standardized (e.g., CV) for cross-trophic-level comparison

# Time series 1 = 1978-1987
ts1_zp <- ts1_zp %>%
  group_by(region, season) %>%
  mutate(sd_zp = sd(log_mean_zp, na.rm = TRUE)) %>%
  ungroup()

# Time series 2 = 1998-now
ts2_zp <- ts2_zp %>%
  group_by(region, season) %>%
  mutate(sd_zp = sd(log_mean_zp, na.rm = TRUE)) %>%
  ungroup()

## ------------------------------------------ ##
#    5) Final Output Tables           ----
## ------------------------------------------ ##

# add station-level metadata (temp, salt, lat, lon) to summarized output
# metadata comes from distinct() = one row per year/season/region 
# so represents a single station, not a spatial average. CHECK

assemble_final <- function(rm_df, full_df) {
  distinct_meta <- full_df %>%
    distinct(region, season, year, .keep_all = TRUE) %>%
    select(-log_mean_zp)

  rm_df %>%
    left_join(distinct_meta, by = c("region", "season", "year")) %>%
    relocate(region, season, year, log_mean_zp, run_mean_zp, sd_zp,
             log10_zp, volume_1m2, min_nonzero,
             lat, lon, sfc_temp, sfc_salt, btm_temp, btm_salt)
}

ts1_final <- assemble_final(ts1_rm, ts1_zp)
ts2_final <- assemble_final(ts2_rm, ts2_zp)

# both log_mean_zp (annual mean) and run_mean_zp (5-yr running
# mean) are retained in the output 

# write.csv(ts1_final, "output/trophamp_zp_runmean_NES_1978_1987.csv", row.names = FALSE)
# write.csv(ts2_final, "output/trophamp_zp_runmean_NES_1998_present.csv", row.names = FALSE)


## CHECK IF USING ENV AVERAGE
## --- Summarize environmental conditions (one row per region/season/year) ---
assemble_final <- function(rm_df, env_df) {
  rm_df %>%
    left_join(env_df, by = c("region", "season", "year")) %>%
    relocate(region, season, year, log_mean_zp, run_mean_zp, sd_zp,
             mean_sfc_temp, mean_sfc_salt, mean_btm_temp, mean_btm_salt)
}

ts1_final <- assemble_final(ts1_rm, ts1_env)
ts2_final <- assemble_final(ts2_rm, ts2_env)

## ------------------------------------------ ##
#    Plots                                     ----
## ------------------------------------------ ##

## --- Annual log means TS2 ---
ggplot(ts2_zp_a, aes(x = year, y = log_mean_zp)) +
  geom_point() +
  geom_line(data = ts2_rm, aes(x = year, y = run_mean_zp),
            color = "red", linewidth = 1.2) +
  facet_grid(season ~ region) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12,
                                   color = "black"),
        panel.grid = element_blank()) +
  labs(title = "Zooplankton Displacement Volume — TS2 (1998–present)",
       y = "log10 ZP (mL/m2)", x = NULL)

# Annual log means — TS1
ggplot(ts1_zp_a, aes(x = year, y = log_mean_zp)) +
  geom_point() +
  geom_line(data = ts1_rm, aes(x = year, y = run_mean_zp),
            color = "red", linewidth = 1.2) +
  facet_grid(season ~ region) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12,
                                   color = "black"),
        panel.grid = element_blank()) +
  labs(title = "Zooplankton Displacement Volume — TS1 (1978–1987)",
       y = "log10 ZP (mL/m2)", x = NULL)
