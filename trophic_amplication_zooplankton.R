

# trophic amp zooplankton
## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##

library(tidyverse) #v2.0.0
library(runner) #for running mean calc; v0.4.3
#library(zoo) has other runmean functions
library(janitor)

# ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##
zp_full <- read_csv(here::here("raw",
                               "EcoMon_Plankton_Data_v3_10_wStrataMeta.csv")) %>%
  clean_names()

## ------------------------------------------ ##
#            Tidy Data -----
## ------------------------------------------ ##

# add season column  
zp <- zp_full %>%
  mutate(season = case_when(between(month_num, 3, 5) ~ "Spring",
                            between(month_num, 6, 8) ~ "Summer",
                            between(month_num, 9, 11) ~ "Fall",
                            TRUE ~ "Winter"))

# select regions (excluding CC and NS) [Region 0]
# assign region names
zp <- zp %>%
  filter(region %in% c("SNE","MAB","GOM","GB")) %>%
  # rename SNE as MAB to align w phyto ts
  mutate(region = case_when(
         region == "SNE" ~ "MAB",
         TRUE ~ region
  ))
# Region 1 = MAB
# Region 2 = SNE
# Region 3 = GB
# Region 4 = GOM

## ------------------------------------------ ##
#      Zooplankton DISPLACEMENT VOLUME -----
## ------------------------------------------ ##
#############################################################
sum(is.na(zp$volume_1m2))
sum(is.na(zp$volume_100m3)) 
#3700 

zp_full %>%
  summarize(across(ctyp_10m2:pnepau_10m2, ~ sum(is.na(.))))
# for some reason there are observations that have zp abundance but the 
# volume is still NA - so maybe they didnt measure vol sometimes


# subset columns for zooplankton
zp <- zp %>%
  select(
    cruise_name, station, zoo_gear, ich_gear, date, day, month, month_num, year,
    time, depth, season, region, area, type, lon, lat, sfc_temp, sfc_salt, 
    btm_temp, btm_salt, volume_1m2) %>% 
  drop_na(volume_1m2)
# volume_1m2 = Zooplankton Displacement Volume (ml) per 1m2 of surface area
# does this not include the fish displacement...?? 


## ------------------------------------------ ##
#           1) Data Transformation -----
## ------------------------------------------ ##
# find the minimum non-zero value for each region and season

zp <- zp %>%
  group_by(region, season) %>%
  mutate(min_nonzero = min(volume_1m2[volume_1m2 > 0],
                           na.rm = TRUE)) %>%
  ungroup()

# Log transformation 
zp$log10_zp <- log10(zp$volume_1m2 + zp$min_nonzero/2)

ggplot(zp, aes(x = year, y = log10_zp)) +
  geom_point() +
  facet_grid(season~region, scales = "free")

## 1978-1988
## split into two 
ts1_zp <- zp %>%
  filter(between(year, 1978, 1987))

ts2_zp <- zp %>%
  filter(between(year, 1998, max(year)))


## ------------------------------------------ ##
#           2) Average -----
## ------------------------------------------ ##
# average across stations for a cruise or year/season

## ----
# Entire time series
## ----
zp <- zp %>%
  group_by(region, year, season) %>% 
  mutate(log_mean_zp = mean(log10_zp, na.rm = T)) %>%
  ungroup()

zp_a <- zp %>%
  group_by(region, year, season) %>%
  summarize(log_mean_zp = mean(log10_zp, na.rm = T), 
            .groups = "drop")

## ----
# Time series 1 = 1978-1987
## ----
ts1_zp <- ts1_zp %>%
  group_by(region, year, season) %>% 
  mutate(log_mean_zp = mean(log10_zp, na.rm = T)) %>%
  ungroup()

## ----
# Time series 2 = 1998-now
## ----
ts2_zp <- ts2_zp %>%
  group_by(region, year, season) %>% 
  mutate(log_mean_zp = mean(log10_zp, na.rm = T)) %>%
  ungroup()

ts2_zp_a <- ts2_zp %>%
  group_by(region, year, season) %>%
  summarize(log_mean_zp = mean(log10_zp, na.rm = T),
            .groups = "drop")

##
ggplot(ts1_zp, aes(x = year, y = log_mean_zp)) +
  geom_point() +
  facet_grid(season~region) 


## ------------------------------------------ ##
#    3) Running mean = 5 years  -----
## ------------------------------------------ ##
# take a running mean with timespan ~ longest lived taxon
# smoothing over life span of organisms 

## ----
# Entire time series 
## ----
rm_zp_tsFULL <- zp %>%
  distinct(year, region, season, .keep_all = T) %>%
  arrange(region, season, year) %>% 
  group_by(region, season) %>% 
  mutate(run_mean_zp = runner::mean_run(log_mean_zp, na_rm = T, k=5))
#mutate(running_mean = rollmean(log_mean_ZP, k = 5, fill = NA, align = "right"))

# same as above but less columns 
rm_zp_tsFULL_a <- zp_a %>%
  #distinct(year, region_name, season, .keep_all = T) %>%
  arrange(region, season, year) %>% 
  group_by(region, season) %>% 
  #summarize(runMean_FF = runner::mean_run(log_mean_FF, na_rm = TRUE, k=5))
  mutate(run_mean_zp = runner::mean_run(log_mean_zp, na_rm = T, k = 5))

## ----
# Time series 1 = 1978-1987 
## ----
rm_zp_ts1 <- ts1_zp %>%
  distinct(year, region, season, .keep_all = T) %>%
  arrange(region, season, year) %>% 
  group_by(region, season) %>%
  mutate(run_mean_zp = runner::mean_run(log_mean_zp, na_rm = T, k=5))

## ----
# Time series 2 = 1998-now
## ----
rm_zp_ts2 <- ts2_zp %>%
  distinct(year, region, season, .keep_all = T) %>%
  arrange(region, season, year) %>%
  group_by(region, season) %>%
  mutate(run_mean_zp = runner::mean_run(log_mean_zp, na_rm = T, k=5))

rm_zp_ts2_a <- ts2_zp_a %>%
  arrange(region, season, year) %>% 
  group_by(region, season) %>% 
  mutate(run_mean_zp = runner::mean_run(log_mean_zp, na_rm = T, k = 5))


## ------------------------------------------ ##
#    4) Compute SD   -----
## ------------------------------------------ ##
# compute st.dev. of time-series 

## ----
# Entire time series
## ----
zp <- zp %>%
  group_by(region, season) %>%
  mutate(sd_zp = sd(log_mean_zp, na.rm = T))
#write.csv(zp, "output/trophamp_zp_NES.csv")

#SD only
sd_zp <- zp %>% 
  select(region, season, sd_zp) %>%
  distinct(region, season, .keep_all = T) %>%
  arrange(desc(region))
#write.csv(sd_zp, "output/SD_trophamp_zp_NES.csv")

distinct_tsFull_zp <- zp %>%
  distinct(region, season, year, .keep_all = T)

rm_zp_tsFull_FINAL <- rm_zp_tsFULL_a %>%
  left_join(distinct_tsFull_zp %>% 
              select(-log_mean_zp), 
            by = c("region", "season", "year")) %>%
  select(all_of(setdiff(names(.), names(rm_zp_tsFULL_a))), 
         everything()) %>%
  #reorder col
  relocate(region, .after = region) %>%
  relocate(year, .after = day) %>%
  relocate(season, .after = depth) %>%
  relocate(log_mean_zp, .after = log10_zp)
#write.csv(rm_zp_tsFull_FINAL, "output/trophamp_zp_runmean_NES.csv")

## ----
# Time series 1 = 1978-1987
## ----
ts1_zp <- ts1_zp %>%
  #distinct(region, season, .keep_all = TRUE) %>%
  group_by(region, season) %>%
  mutate(SD_zp = sd(log_mean_zp, na.rm = T)) 

## ----
# Time series 2 = 1998-now
## ----
ts2_zp <- ts2_zp %>%
  group_by(region, season) %>%
  mutate(sd_zp = sd(log_mean_zp, na.rm = T)) 
#write.csv(ts2_zp, "output/trophamp_zp_NES_1998_2021.csv")

distinct_ts2_zp <- ts2_zp %>%
  distinct(region, season, year, .keep_all = T)

rm_zp_ts2_FINAL <- rm_zp_ts2_a %>%
  left_join(distinct_ts2_zp %>% 
              select(-log_mean_zp), 
            by = c("region", "season", "year")) %>%
  select(all_of(setdiff(names(.), names(rm_zp_ts2_a))), 
         everything()) %>%
  #reorder col
  relocate(region, .after = region) %>%
  relocate(year, .after = day) %>%
  relocate(season, .after = depth) %>%
  relocate(log_mean_zp, .after = log10_zp)
#write.csv(rm_zp_ts2_FINAL, "output/trophamp_zp_runmean_NES_1998_2021.csv")

## ------------------------------------------ ##
#    Plots   -----
## ------------------------------------------ ##
ggplot() +
  geom_point(data = ZP, aes(x = year, y = log10_ZP)) +
  geom_point(data = ZP, aes(x = year, y = log_mean_ZP),
             color = "red", shape = 1, size = 2) +
  geom_line(data = rm_zp_tsFULL, aes(x = year, y = runMean_ZP), 
            color = "red", 
            size = 2) +
  ggtitle("ZP") + 
  facet_grid(season~region_name, scales = "free") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12,
                                   color = "black"),
        panel.grid = element_blank())

ggplot() +
  #geom_point(data = sub_ZP_v2, aes(x = year, y = log10Sum)) +
  #geom_point(data = ZP, aes(x = year, y = aver),
  #           color = "red", shape = 1, size = 2) +
  geom_line(data = ZP, aes(x = year, y = log_mean_ZP),
            color = "black") +
  geom_line(data = rm_zp_tsFULL, aes(x = year, y = runMean_ZP), 
            color = "red", 
            size = 2) +
  ggtitle("ZP") + 
  facet_grid(season~region_name, scales = "free") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12,
                                   color = "black"),
        panel.grid = element_blank())
