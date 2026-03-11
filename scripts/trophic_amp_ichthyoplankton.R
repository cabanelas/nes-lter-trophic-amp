################################################################################
## Script:  trophic_amp_ichthyoplankton.R
## Project: NES-LTER Trophic Amplification - ichthyoplankton
##          Cross-Site LTER Pelagic Synthesis Working Group 
## Data:    EcoMon Plankton Survey 
## Author:  Alexandra C. Cabanelas Bermudez
## Created: November 2023  |  Updated: March 2026
##
## Purpose: Prepares ichthyoplankton data from the NES for trophic amplification
##          analysis. Produces log-transformed annual means and 5-year running 
##          means by region and season for two time periods:              
##          TS1: 1978–1987                           
##          TS2: 1998–present (matched to chl-a)     
##       1) log10(x+ (min/2))    for each station ; x = sum trophic level 
##       2) average across stations for a cruise or year/season
##       3) take a running mean with timespan ~ longest lived taxon 
##       4) compute st.dev. of time-series 
##
## Inputs (data/raw/):
##   - EcoMon_plankton_v3_10.csv  (187513.4.4.tar.gz)
## https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:0187513
## Data available for cruises through 2023 - As of 2024_11_06
##
## Outputs (data/output/):
##   - 
################################################################################

## ------------------------------------------ ##
#            Packages
## ------------------------------------------ ##
library(tidyverse)  
library(runner)     # running mean
library(janitor)

## ------------------------------------------ ##
#            Data  
## ------------------------------------------ ##
#EcoMon_Plankton_Data_v3_10_wStrataMeta.csv
ichthyo_full <- read_csv(here::here("data", "raw",
                                    "EcoMon_plankton_v3_10.csv")) %>%
  clean_names()

## ------------------------------------------ ##
#            Tidy Data  
## ------------------------------------------ ##

## --- Filter regions of interest ---
# exclude CC and NS [Region 0]
# rename SNE to MAB to align with phytoplankton time series
# Region 1 = MAB
# Region 2 = SNE (will become MAB)
# Region 3 = GB
# Region 4 = GOM
ichthyo <- ichthyo_full %>%
  filter(region %in% c("SNE", "MAB", "GOM", "GB")) %>%
  # rename SNE as MAB to align w phyto ts
  mutate(region = case_when(
    region == "SNE" ~ "MAB",
    TRUE            ~ region
  ))

## --- Larval Forage Fish ---

# when the ich_gear cell is empty it seems like they did not sample for larval f
ichthyo_full %>%
  summarise(across(bretyr_10m2:lopame_10m2, ~ sum(is.na(.))))
sum(is.na(ichthyo_full$ich_gear))
# 3500 rows that dont have icht data

# atlantic herring == Clupea harengus
# atlantic mackerel == scomber scombrus
# atlantic butterfish == peprilus spp
# sand lance == ammodytes spp

# alewife == not in data = river h = alosa
# blueback herring == not in data = river h = alosa
ichthyo <- ichthyo %>%
  select(
    cruise_name, station, zoo_gear, ich_gear, date, day, month, month_num, year,
    time, depth, season, region, area, type, lon, lat, 
    sfc_temp, sfc_salt, btm_temp, btm_salt, 
    cluhar_10m2, # atlantic herring == Clupea harengus
    scosco_10m2, # atlantic mackerel == scomber scombrus
    pepspp_10m2, # atlantic butterfish == peprilus spp
    ammspp_10m2 # sand lance == ammodytes spp
    #cluhar_100m3, scosco_100m3, pepspp_100m3, ammspp_100m3
  ) %>%
  filter(!is.na(ich_gear))
# nofish_abnd; bretyr_abnd::lopame_abnd

ichthyo %>%
  count(cruise_name, station, region, date) %>%
  filter(n > 1)
# duplicate entries for:
#cruise == HB2303   station ==  14
#cruise == PC2104   station ==  122 
#cruise == PC2106   station ==  88
# fish values identical; zoop values differ
ichthyo <- ichthyo %>%
  arrange(cruise_name, station, region, date,
          desc(zoo_gear == "6B3Z")) %>%
  distinct(cruise_name, station, region, date, .keep_all = TRUE)

rm(ichthyo_full)
## ------------------------------------------ ##
#            Total ichthyo Abundance -----
## ------------------------------------------ ##
# abundance 10m2
# --- sum abundances of diff spp. to get ichthyo total ---
ichthyo <- ichthyo %>%
  group_by(cruise_name, station, region, date) %>% #grouping by sample
  mutate(ichthyosum = sum(cluhar_10m2, 
                          scosco_10m2, 
                          pepspp_10m2, 
                          ammspp_10m2, na.rm = T)) %>%
  ungroup()

## --- plot --- 
ggplot(ichthyo, aes(x = year, y = ichthyosum)) +
  geom_point() +
  facet_wrap(~region, scales = "free")

ggplot(ichthyo, aes(x = year, y = ichthyosum)) +
  geom_point() +
  facet_grid(season~region, scales = "free")

## ------------------------------------------ ##
#    Split into Two Time Series
## ------------------------------------------ ##

## --- Time series 1 = 1978-1987 ---
ts1_ichthyo <- ichthyo %>%
  filter(between(year, 1978, 1987)) #starts 1977

## --- Time series 2 = 1998-now ---
# matched to chla time series
ts2_ichthyo <- ichthyo %>%
  filter(year >= 1998)

## ------------------------------------------ ##
#    1) Log Transformation
## ------------------------------------------ ##
# find the minimum non-zero value for each region and season

log_transform <- function(df) {
  df %>%
    group_by(region, season) %>%
    mutate(
      min_nonzero      = min(ichthyosum[ichthyosum > 0], na.rm = TRUE),
      log10_ichthyo    = log10(ichthyosum + min_nonzero / 2)
    ) %>%
    ungroup()
}

ts1_ichthyo <- log_transform(ts1_ichthyo) # TS1: 1978–1987 
ts2_ichthyo <- log_transform(ts2_ichthyo) # TS2: 1998–present

## --- plot ---
ggplot(ts2_ichthyo, aes(x = year, y = log10_ichthyo)) +
  geom_point() +
  facet_grid(season~region, scales = "free")

## ------------------------------------------ ##
#    2) Annual Means by Region & Season
## ------------------------------------------ ##
# average across stations for a cruise or year/season

# two versions, needed downstream:
# ts_ichthyo (mutate)        - keeps all station-level rows, log_mean_ichthyo as 
#                              repeated column. 
# ts_ichthyo_sum (summarize) - collapses to one row per region/season/year.
#                              Used for running mean (step 3).
# log_mean_ichthyo values are identical between the two - only shape differs
## ------------------------------------------ ##

## --- Time series 1 = 1978-1987 ---
## --- Full station-level data with annual mean column attached ---
# station-level rows retained 
ts1_ichthyo <- ts1_ichthyo %>%
  group_by(region, season, year) %>%
  mutate(
    log_mean_ichthyo   = mean(log10_ichthyo, na.rm = TRUE),
    mean_sfc_temp      = mean(sfc_temp, na.rm = TRUE),
    mean_sfc_salt      = mean(sfc_salt, na.rm = TRUE),
    mean_btm_temp      = mean(btm_temp, na.rm = TRUE),
    mean_btm_salt      = mean(btm_salt, na.rm = TRUE)
  ) %>%
  ungroup()

## --- Summarize (one row per region/season/year) ---
# used for running mean
ts1_ichthyo_sum <- ts1_ichthyo %>%
  group_by(region, season, year) %>%
  summarize(
    log_mean_ichthyo = mean(log10_ichthyo, na.rm = TRUE),
    mean_sfc_temp    = mean(sfc_temp, na.rm = TRUE),
    mean_sfc_salt    = mean(sfc_salt, na.rm = TRUE),
    mean_btm_temp    = mean(btm_temp, na.rm = TRUE),
    mean_btm_salt    = mean(btm_salt, na.rm = TRUE),
    .groups = "drop" 
  )

## --- Time series 2 = 1998-now ---
## --- Full station-level data with annual mean column attached ---
# station-level rows retained 
ts2_ichthyo <- ts2_ichthyo %>%
  group_by(region, season, year) %>%
  mutate(
    log_mean_ichthyo   = mean(log10_ichthyo, na.rm = TRUE),
    mean_sfc_temp      = mean(sfc_temp, na.rm = TRUE),
    mean_sfc_salt      = mean(sfc_salt, na.rm = TRUE),
    mean_btm_temp      = mean(btm_temp, na.rm = TRUE),
    mean_btm_salt      = mean(btm_salt, na.rm = TRUE)
  ) %>%
  ungroup()

## --- Summarize (one row per region/season/year) ---
# used for running mean
ts2_ichthyo_sum <- ts2_ichthyo %>%
  group_by(region, season, year) %>%
  summarize(
    log_mean_ichthyo   = mean(log10_ichthyo, na.rm = TRUE),
    mean_sfc_temp      = mean(sfc_temp, na.rm = TRUE),
    mean_sfc_salt      = mean(sfc_salt, na.rm = TRUE),
    mean_btm_temp      = mean(btm_temp, na.rm = TRUE),
    mean_btm_salt      = mean(btm_salt, na.rm = TRUE),
    .groups = "drop"
  )

## --- plot ---
ggplot(ts2_ichthyo_sum, aes(x = year, y = log_mean_ichthyo)) +
  geom_point() +
  facet_grid(season~region) 

## ------------------------------------------ ##
#    3) 5-Year Running Mean
## ------------------------------------------ ##
# take a running mean with timespan ~ longest lived taxon
# smoothing over life span of organisms 
# k = 5 
# run_mean_zp uses all available years up to k (right-aligned, na_rm = TRUE)

## --- Time series 1 = 1978-1987 ---
ts1_rm <- ts1_ichthyo_sum %>%
  arrange(region, season, year) %>%
  group_by(region, season) %>%
  mutate(run_mean_ichthyo = runner::mean_run(log_mean_ichthyo, k = 5, 
                                             idx = year, na_rm = TRUE)) %>%
  ungroup()

## --- Time series 2 = 1998-now ---
ts2_rm <- ts2_ichthyo_sum %>%
  arrange(region, season, year) %>%
  group_by(region, season) %>%
  mutate(run_mean_ichthyo = runner::mean_run(log_mean_ichthyo, k = 5, 
                                             idx = year, na_rm = TRUE)) %>%
  ungroup()

## --- plot ---
ggplot(ts2_rm, aes(x = year, y = run_mean_ichthyo)) +
  geom_point() +
  facet_grid(season~region) 

## ------------------------------------------ ##
#    4) Standard Deviation
## ------------------------------------------ ##
# SD of annual log-means across each time series 
# = interannual variability metric for trophic amplification test
# does variability increase at higher trophic levels?

## --- Time series 1 = 1978-1987 ---
ts1_ichthyo <- ts1_ichthyo %>%
  group_by(region, season) %>%
  mutate(sd_ichthyo = sd(log_mean_ichthyo, na.rm = TRUE)) %>%
  ungroup()

## --- Time series 2 = 1998-now ---
ts2_ichthyo <- ts2_ichthyo %>%
  group_by(region, season) %>%
  mutate(sd_ichthyo = sd(log_mean_ichthyo, na.rm = TRUE)) %>%
  ungroup()

## --- plot ---
ggplot(ts2_ichthyo, aes(x = year, y = sd_ichthyo)) +
  geom_point() +
  facet_grid(season~region) 

## ------------------------------------------ ##
#    Final Output Tables        
## ------------------------------------------ ##
# ts_ichthyo = station-level, has all cruise info + sd_ichthyo
# ts_rm      = one row per region/season/year, has run_mean_ichthyo
# 
# FULL:  all station rows, run_mean_ichthyo joined on
# SUM:   one row per region/season/year (distinct), run_mean_ichthyo joined on

## --- Time series 1 ---
ts1_ichthyo_full <- ts1_ichthyo %>%
  left_join(ts1_rm %>% select(region, season, year, run_mean_ichthyo),
            by = c("region", "season", "year"))

ts1_ichthyo_sum <- ts1_ichthyo_full %>%
  distinct(region, season, year, .keep_all = TRUE) %>%
  select(region, season, year, min_nonzero, log10_ichthyo, log_mean_ichthyo, 
         mean_sfc_temp, mean_sfc_salt, mean_btm_temp, mean_btm_salt, 
         sd_ichthyo, run_mean_ichthyo
  )

## --- Time series 2 ---
ts2_ichthyo_full <- ts2_ichthyo %>%
  left_join(ts2_rm %>% select(region, season, year, run_mean_ichthyo),
            by = c("region", "season", "year"))

ts2_ichthyo_sum <- ts2_ichthyo_full %>%
  distinct(region, season, year, .keep_all = TRUE) %>%
  select(region, season, year, min_nonzero, log10_ichthyo, log_mean_ichthyo, 
         mean_sfc_temp, mean_sfc_salt, mean_btm_temp, mean_btm_salt, 
         sd_ichthyo, run_mean_ichthyo
  )

# write.csv(ts1_ichthyo_full, "output/trophamp_ichthyo_1978_1987_full.csv", row.names = FALSE)
# write.csv(ts1_ichthyo_sum,  "output/trophamp_ichthyo_1978_1987_sum.csv",  row.names = FALSE)
# write.csv(ts2_ichthyo_full, "output/trophamp_ichthyo_1998_2023_full.csv", row.names = FALSE)
# write.csv(ts2_ichthyo_sum,  "output/trophamp_ichthyo_1998_2023_sum.csv",  row.names = FALSE)

## --- SD ---
ts2_sd <- ts2_ichthyo_sum %>%
  distinct(region, season, sd_ichthyo)
#write.csv(ts2_sd, "output/trophamp_ichthyo_1998_2023_sd.csv", row.names = FALSE)

## ------------------------------------------ ##
#    Plots                                 
## ------------------------------------------ ##

## --- Annual log means TS2 ---
ggplot(ts2_ichthyo_sum, aes(x = year, y = log_mean_ichthyo)) +
  geom_point() +
  geom_line(data = ts2_rm, aes(x = year, y = run_mean_ichthyo),
            color = "red", linewidth = 1.2) +
  facet_grid(season ~ region) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12,
                                   color = "black"),
        panel.grid = element_blank()) +
  labs(title = "Ichthyoplankton — TS2 (1998–2023)",
       y = "log10", x = NULL)

## --- SD TS2 ---
ggplot(ts2_ichthyo_sum, aes(x = region, y = sd_ichthyo, color = season)) +
  geom_point(size = 4) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "Ichthyoplankton",
    subtitle = "Time series 2: 1998–2023",
    x = "Region",
    y = "SD of log10 annual mean",
    color = "Season"
  ) +
  theme_bw()

## --- anomaly --- 
ts2_anom <- ts2_ichthyo_sum %>%
  group_by(region, season) %>%
  mutate(ichthyo_anom = log_mean_ichthyo - mean(log_mean_ichthyo, na.rm = TRUE)) %>%
  ungroup()

ggplot(ts2_anom, aes(year, ichthyo_anom)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_col() +
  facet_grid(season ~ region) +
  theme_bw() +
  labs(y = "Anomaly in log10 annual mean ichthyo", x = NULL)

ggplot(ts2_anom, aes(x = year, y = ichthyo_anom, color = region)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~season) +
  theme_bw() +
  labs(
    y = "Anomaly in log10 annual mean ichthyo",
    x = NULL,
    color = "Region"
  )

## --- distribution --- 
ggplot(ts2_ichthyo_sum, aes(region, log_mean_ichthyo, fill = season)) +
  geom_boxplot() +
  theme_bw() +
  labs(y = "log10 annual mean ichthyo", x = NULL)

## --- 
ts2_check <- ts2_ichthyo_sum %>%
  arrange(region, season, year) %>%
  group_by(region, season) %>%
  mutate(
    n_in_window = runner::runner(
      x = log_mean_ichthyo,
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
