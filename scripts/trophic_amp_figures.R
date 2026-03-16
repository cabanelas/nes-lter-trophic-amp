################################################################################
## Script:  trophic_amp_figures.R
## Project: NES-LTER Trophic Amplification
##          Cross-Site LTER Pelagic Synthesis Working Group
## Author:  Alexandra C. Cabanelas Bermudez
## Created: March 2026
##
## Purpose: Plotting.
##
## Inputs (data/output/):
##   Zooplankton (from trophic_amp_zooplankton.R):
##     - trophamp_zp_1998_2023_full.csv
##   Ichthyoplankton (from trophic_amp_ichthyoplankton.R):
##     - trophamp_ichthyo_1998_2023_full.csv
##   Forage Fish (from trophic_amp_foragefish.R):
##     - trophamp_ff_1998_2023_full.csv
##
##   Required for Figure 1 (station map):
##     - EcomonStrata_v4.shp  (data/raw/)
##     - EcomonStrata_v4b.shp (data/raw/)
##
##   Required for Figure 2 (raw vs corrected CPUE):
##     NEFSC BTS catch + metadata CSVs (data/raw/) and the
##     prep_survey() helper.
##
##   Required for Figure 7, 8, 9 
##     - EcoMon_plankton_v3_10.csv (data/raw/)
##
## Outputs (figures/):
##   - fig01_station_map.png
##   - fig02_corrections_combined.png
##   - fig02_corrections_by_season.png
##   - fig03_ts_zp.png
##   - fig03_ts_ichthyo.png
##   - fig03_ts_ff.png
##   - fig04_anom_zp.png
##   - fig04_anom_ichthyo.png
##   - fig04_anom_ff.png
##   - fig05_nwindow_zp.png
##   - fig05_nwindow_ichthyo.png
##   - fig05_nwindow_ff.png
##   - fig06a_zp_regional_ts.png
##   - fig06b_zp_regional_dist.png
##   - fig07a_zp_mab_sne_ts.png
##   - fig07b_zp_regional_boxplot.png
##   - fig08a_ip_mab_sne_ts.png
##   - fig08b_ip_regional_boxplot.png
##   - fig09a_ip_bytaxa_allregions_ts.png
##   - fig09b_ip_bytaxa_allregions_boxplot.png
################################################################################

SAVE_FIGURES <- TRUE   # set FALSE to only print to screen
FIG_DIR      <- "figures"

## ------------------------------------------ ##
#   Packages
## ------------------------------------------ ##
library(tidyverse)
library(patchwork)
library(gt)
library(sf)
library(runner)
library(janitor)
library(rnaturalearth)
library(rnaturalearthdata)

if (SAVE_FIGURES) dir.create(FIG_DIR, showWarnings = FALSE)

## ------------------------------------------ ##
#   Helper: save wrapper
## ------------------------------------------ ##
save_fig <- function(plot, filename, width = 10, height = 7, dpi = 300) {
  if (SAVE_FIGURES) {
    ggsave(
      filename = file.path(FIG_DIR, filename),
      plot     = plot,
      width    = width,
      height   = height,
      dpi      = dpi
    )
    message("Saved: ", file.path(FIG_DIR, filename))
  }
  invisible(plot)
}

## ------------------------------------------ ##
#   Data
## ------------------------------------------ ##

## --- Zooplankton ---
zp_full <- read_csv(here::here("data", "output",
                               "trophamp_zp_1998_2023_full.csv"))

## --- Ichthyoplankton ---
ichthyo_full <- read_csv(here::here("data", "output",
                                    "trophamp_ichthyo_1998_2023_full.csv"))
## --- Forage Fish ---
ff_full <- read_csv(here::here("data", "output",
                               "trophamp_ff_1998_2023_full.csv"))

## --- Collapse to one row per region/season/year for time series plots ---
zp <- zp_full %>%
  distinct(region, season, year, .keep_all = TRUE) %>%
  select(region, season, year, log_mean_zp, run_mean_zp, sd_zp)

ichthyo <- ichthyo_full %>%
  distinct(region, season, year, .keep_all = TRUE) %>%
  select(region, season, year, log_mean_ichthyo, run_mean_ichthyo, sd_ichthyo)

ff <- ff_full %>%
  distinct(region, season, year, .keep_all = TRUE) %>%
  select(region, season, year, log_mean_ff, run_mean_ff, sd_ff)

## --- Combined running mean for Fig 3 ---
rm <- zp %>%
  left_join(
    ichthyo %>% select(region, season, year, run_mean_ichthyo),
    by = c("region", "season", "year")
  )

## ------------------------------------------ ##
#   FIGURE 1: Station Coverage Map
## ------------------------------------------ ##

## --- Load EcoMon shapefile for basemap ---
ecomap1 <- st_read(here::here("data", "raw", "EcomonStrata_v4.shp"),
                   quiet = TRUE)
ecomap2 <- st_read(here::here("data", "raw", "EcomonStrata_v4b.shp"),
                   quiet = TRUE)
ecomap         <- rbind(ecomap1, ecomap2)
st_crs(ecomap) <- 4326
ecomap_valid   <- st_make_valid(ecomap)

## --- Build combined station layer ---
map_zp      <- zp_full      %>%
  select(lon, lat, region, season) %>% mutate(source = "Zooplankton")
map_ichthyo <- ichthyo_full %>%
  select(lon, lat, region, season) %>% mutate(source = "Ichthyoplankton")
map_ff      <- ff_full      %>%
  select(lon, lat, region, season) %>% mutate(source = "Forage fish")

map_all <- bind_rows(map_zp, map_ichthyo, map_ff) %>%
  mutate(source = factor(source,
                         levels = c("Zooplankton", "Ichthyoplankton", 
                                    "Forage fish")))

world <- ne_countries(scale = "medium", returnclass = "sf")

fig01 <- ggplot() +
  geom_sf(data = ecomap_valid, fill = "gray92", 
          color = "gray60", linewidth = 0.3) +
  geom_sf(data = world, fill = "gray80", color = "gray40") +
  geom_point(data = map_all,
             aes(x = lon, y = lat, color = region),
             alpha = 0.3, size = 0.8) +
  coord_sf(xlim = c(-77, -61), ylim = c(35, 46)) +
  facet_wrap(~source, ncol = 3) +
  scale_color_brewer(palette = "Set1", name = "Region") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  labs(title = "Station Coverage (1998–present)", x = NULL, y = NULL) +
  theme_bw() +
  theme(strip.text      = element_text(face = "bold"),
        axis.text       = element_text(size = 7),
        legend.position = "bottom")

print(fig01)
# save_fig(fig01, "fig01_station_map.png", width = 12, height = 5)
rm(ecomap, ecomap_valid, ecomap1, ecomap2)

## ------------------------------------------ ##
#   FIGURE 2: Raw vs. Corrected Forage Fish CPUE 
## ------------------------------------------ ##

## --- prep_survey: join catch to haul metadata, sum across sexes ---
prep_survey <- function(catch, meta, season_label) {
  catch %>%
    filter(svspp %in% forage_codes) %>%
    left_join(
      meta %>% select(id, svvessel, est_year, est_month,
                      decdeg_beglat, decdeg_beglon, avgdepth,
                      surftemp, surfsalin, bottemp, botsalin),
      by = "id"
    ) %>%
    rename(year         = est_year,
           month        = est_month,
           lat          = decdeg_beglat,
           lon          = decdeg_beglon,
           depth        = avgdepth,
           wtcpue       = expcatchwt) %>%
    mutate(season_survey = season_label,
           wtcpue        = as.numeric(wtcpue),
           svvessel      = as.character(svvessel)) %>%
    group_by(id, svspp, year, month, lat, lon, depth, stratum,
             svvessel, season_survey, surftemp, surfsalin, 
             bottemp, botsalin) %>%
    summarise(wtcpue = sum(wtcpue, na.rm = TRUE), .groups = "drop")
}
  
## --- Load raw files ---
fall_catch <- read_csv(here::here("data", "raw",
  "22560_NEFSCFallFisheriesIndependentBottomTrawlData",
                                  "22560_UNION_FSCS_SVCAT.csv"),
                       col_types = cols(ID = col_character())) %>% 
  clean_names()
  
fall_meta <- read_csv(here::here("data", "raw",
  "22560_NEFSCFallFisheriesIndependentBottomTrawlData",
                                 "22560_UNION_FSCS_SVSTA.csv"),
                      col_types = cols(ID = col_character())) %>% 
  clean_names()
  
spring_catch <- read_csv(here::here("data", "raw",
  "22561_NEFSCSpringFisheriesIndependentBottomTrawlData",
                                    "22561_UNION_FSCS_SVCAT.csv"),
                         col_types = cols(ID = col_character())) %>% 
  clean_names()
  
spring_meta <- read_csv(here::here("data", "raw",
  "22561_NEFSCSpringFisheriesIndependentBottomTrawlData",
                                   "22561_UNION_FSCS_SVSTA.csv"),
                        col_types = cols(ID = col_character())) %>% 
  clean_names()
  
spp_codes <- read_csv(here::here("data", "raw",
  "22560_NEFSCFallFisheriesIndependentBottomTrawlData",
                                 "Fall_SVDBS_SupportTables", 
                                 "SVDBS_SVSPECIES_LIST.csv")) %>%
  clean_names()
  
forage_spp <- spp_codes %>%
  filter(svspp %in% c("032", "121", "131", "181", "734")) %>%
  mutate(target_taxa = case_when(
      svspp == "032"             ~ "Atlantic herring",
      svspp == "121"             ~ "Atlantic mackerel",
      svspp == "131"             ~ "Butterfish",
      svspp %in% c("181", "734") ~ "Sand lance"
    ))
forage_codes <- forage_spp$svspp

## --- Build raw (uncorrected) CPUE ---
raw_ff <- bind_rows(
  prep_survey(fall_catch,   fall_meta,   "Fall"),
  prep_survey(spring_catch, spring_meta, "Spring")
) %>%
  filter(between(year, 1998, 2023)) %>%
  group_by(year) %>%
  summarise(wtcpue = mean(wtcpue, na.rm = TRUE), .groups = "drop") %>%
  mutate(version = "Raw")
  
raw_ff_s <- bind_rows(
  prep_survey(fall_catch,   fall_meta,   "Fall"),
  prep_survey(spring_catch, spring_meta, "Spring")
) %>%
  filter(between(year, 1998, 2023)) %>%
  group_by(season_survey, year) %>%
  summarise(wtcpue = mean(wtcpue, na.rm = TRUE), .groups = "drop") %>%
  mutate(version = "Raw")
  
## --- Corrected CPUE comes straight from ff_full ---
corrected_ff <- ff_full %>%
  group_by(year) %>%
  summarise(wtcpue = mean(wtcpue, na.rm = TRUE), .groups = "drop") %>%
  mutate(version = "Corrected")
  
corrected_ff_s <- ff_full %>%
  group_by(season_survey, year) %>%
  summarise(wtcpue = mean(wtcpue, na.rm = TRUE), .groups = "drop") %>%
  mutate(version = "Corrected")

## --- Plot: combined seasons ---
fig02a <- bind_rows(raw_ff, corrected_ff) %>%
  filter(year >= 1998) %>%
  ggplot(aes(x = year, y = wtcpue, color = version, linetype = version)) +
  geom_vline(xintercept = 2008.5, linetype = "dashed",
             color = "gray40", linewidth = 0.6) +
  annotate("text", x = 2009.2, y = Inf, label = "Bigelow",
           hjust = 0, vjust = 1.3, size = 2.8, color = "gray30") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  scale_color_manual(values = c("Raw"       = "tomato", 
                                "Corrected" = "steelblue")) +
  scale_linetype_manual(values = c("Raw"       = "dashed", 
                                   "Corrected" = "solid")) +
  labs(title    = "Effect of Vessel & Gear Corrections on Forage Fish CPUE",
       subtitle = "1998–2023 | Dashed line = Bigelow transition (2009)",
       x = NULL, y = "Mean CPUE (kg tow\u207B\u00B9)", 
       color = NULL, linetype = NULL) +
  scale_x_continuous(breaks = seq(1998, 2023, by = 5)) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold"), 
        legend.position = c(0.15, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.7), 
                                         color = NA),
        axis.text.x      = element_text(size = 10)
        )

## --- Plot: by survey season ---
fig02b <- bind_rows(raw_ff_s, corrected_ff_s) %>%
  filter(year >= 1998) %>%
  ggplot(aes(x = year, y = wtcpue, color = version, linetype = version)) +
  geom_vline(xintercept = 2008.5, linetype = "dashed",
             color = "gray40", linewidth = 0.5) +
  annotate("text", x = 2009.3, y = Inf, label = "Bigelow",
           hjust = 0, vjust = 1.4, size = 2.5, color = "gray30") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1) +
  scale_color_manual(values = c("Raw"       = "tomato", 
                                "Corrected" = "steelblue")) +
  scale_linetype_manual(values = c("Raw"       = "dashed", 
                                   "Corrected" = "solid")) +
  facet_wrap(~season_survey, scales = "free_y") +
  labs(title    = "Effect of Vessel & Gear Corrections on Forage Fish CPUE",
       subtitle = "1998–2023 | Dashed line = Bigelow transition (2009)",
       x = NULL, y = "Mean CPUE (kg tow\u207B\u00B9)",
       color = NULL, linetype = NULL) +
  theme_bw() +
  theme(strip.text      = element_text(face = "bold"),
        legend.position = "bottom",
        axis.text.x     = element_text(size = 10))

print(fig02a)
print(fig02b)
# save_fig(fig02a, "fig02_corrections_combined.png",  width = 10, height = 7)
# save_fig(fig02b, "fig02_corrections_by_season.png", width = 12, height = 7)

fig02a + 
  bind_rows(raw_ff_s, corrected_ff_s) %>%
  filter(year >= 1998) +
  facet_wrap(~season_survey) +
  theme(legend.position = "bottom")
  
rm(fall_catch, fall_meta, spring_catch, spring_meta,
   spp_codes, forage_spp, forage_codes,
   raw_ff, raw_ff_s, corrected_ff, corrected_ff_s, prep_survey)

## ------------------------------------------ ##
#   FIGURE 3: Annual Log Means + 5-yr Running Mean 
## ------------------------------------------ ##

plot_ts <- function(sum_df, rm_df, y_col, rm_col, title_label, y_label) {
  ggplot(sum_df, aes(x = year, y = .data[[y_col]])) +
    geom_point(size = 1.5, alpha = 0.7, color = "gray30") +
    geom_line(data = rm_df, aes(x = year, y = .data[[rm_col]]),
              color = "tomato", linewidth = 1.1) +
    facet_grid(season ~ region) +
    labs(title    = title_label,
         subtitle = "Points = annual log mean | Red line = 5-year running mean",
         y = y_label, x = NULL) +
    theme_bw() +
    theme(axis.text.x  = element_text(size = 9),
          strip.text   = element_text(face = "bold"),
          panel.grid   = element_blank())
}

fig03_zp <- plot_ts(
  zp, rm,
  "log_mean_zp",   "run_mean_zp",
  "Zooplankton Displacement Volume (1998–2023)",
  "log\u2081\u2080 ZP (mL m\u207B\u00B2)"
)

fig03_ichthyo <- plot_ts(
  ichthyo, rm,
  "log_mean_ichthyo", "run_mean_ichthyo",
  "Ichthyoplankton (1998–2023)",
  "log\u2081\u2080 larval abundance (ind. 10 m\u207B\u00B2)"
)

fig03_ff <- plot_ts(
  ff, ff,
  "log_mean_ff",   "run_mean_ff",
  "Forage Fish CPUE (1998–2023)",
  "log\u2081\u2080 CPUE (kg tow\u207B\u00B9)"
)

print(fig03_zp)
print(fig03_ichthyo)
print(fig03_ff)
# save_fig(fig03_zp,      "fig03_ts_zp.png")
# save_fig(fig03_ichthyo, "fig03_ts_ichthyo.png")
# save_fig(fig03_ff,      "fig03_ts_ff.png")

## ------------------------------------------ ##
#   FIGURE 4: Anomaly Bar Plots
## ------------------------------------------ ##

plot_anom <- function(sum_df, y_col, title_label, y_label) {
  df <- sum_df %>%
    group_by(region, season) %>%
    mutate(anom = .data[[y_col]] - mean(.data[[y_col]], na.rm = TRUE)) %>%
    ungroup()
  
  ggplot(df, aes(x = year, y = anom, fill = anom > 0)) +
    geom_col(show.legend = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    scale_fill_manual(values = c("TRUE" = "tomato", "FALSE" = "steelblue")) +
    facet_grid(season ~ region) +
    labs(title    = title_label,
         subtitle = "Anomaly relative to full time series mean",
         y = expression(log[10]~"annual mean abundance (per 10 m"^2*")"), 
         x = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          strip.text  = element_text(face = "bold"),
          panel.grid  = element_blank())
}

fig04_zp <- plot_anom(
  zp, "log_mean_zp",
  "Zooplankton Anomaly",
  "Anomaly in log\u2081\u2080 annual mean ZP"
)

fig04_ichthyo <- plot_anom(
  ichthyo, "log_mean_ichthyo",
  "Ichthyoplankton Anomaly",
  "Anomaly in log\u2081\u2080 annual mean ichthyo"
)

fig04_ff <- plot_anom(
  ff, "log_mean_ff",
  "Forage Fish Anomaly",
  "Anomaly in log\u2081\u2080 annual mean CPUE"
)

print(fig04_zp)
print(fig04_ichthyo)
print(fig04_ff)
# save_fig(fig04_zp,      "fig04_anom_zp.png")
# save_fig(fig04_ichthyo, "fig04_anom_ichthyo.png")
# save_fig(fig04_ff,      "fig04_anom_ff.png")

## ------------------------------------------ ##
#   FIGURE 5: N Observations in 5-Year Running Mean Window
## ------------------------------------------ ##

make_nwindow <- function(sum_df, y_col) {
  sum_df %>%
    arrange(region, season, year) %>%
    group_by(region, season) %>%
    mutate(
      n_in_window = runner::runner(
        x   = .data[[y_col]],
        k   = 5,
        idx = year,
        f   = \(x) sum(!is.na(x))
      )
    ) %>%
    ungroup()
}

plot_nwindow <- function(check_df, title_label) {
  ggplot(check_df, aes(x = year, y = n_in_window)) +
    geom_line(linewidth = 0.6, color = "steelblue") +
    geom_point(size = 1.5, color = "steelblue") +
    facet_grid(season ~ region) +
    scale_y_continuous(breaks = 1:5) +
    labs(title    = title_label,
         y = "N years in window", x = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          strip.text  = element_text(face = "bold"),
          panel.grid  = element_blank())
}

check_zp      <- make_nwindow(zp,      "log_mean_zp")
check_ichthyo <- make_nwindow(ichthyo, "log_mean_ichthyo")
check_ff      <- make_nwindow(ff,      "log_mean_ff")

fig05_zp      <- plot_nwindow(check_zp,     "Zooplankton — N years in 5-yr window")
fig05_ichthyo <- plot_nwindow(check_ichthyo,"Ichthyoplankton — N years in 5-yr window")
fig05_ff      <- plot_nwindow(check_ff,     "Forage Fish — N years in 5-yr window")

print(fig05_zp)
print(fig05_ichthyo)
print(fig05_ff)
# save_fig(fig05_zp,      "fig05_nwindow_zp.png")
# save_fig(fig05_ichthyo, "fig05_nwindow_ichthyo.png")
# save_fig(fig05_ff,      "fig05_nwindow_ff.png")

## ------------------------------------------ ##
#   FIGURE 6: Zooplankton Regional Time Series & Distribution
## ------------------------------------------ ##

## --- Time series by region ---
fig06a <- zp_full %>%
  group_by(region, season, year) %>%
  summarize(log_mean_zp = mean(log10_zp, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = year, y = log_mean_zp, color = region)) +
  geom_line(linewidth = 1) +
  geom_point() +
  facet_wrap(~season) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  labs(title = "Zooplankton Displacement Volume (1998–present)",
       y = "log\u2081\u2080 ZP (mL m\u207B\u00B2)", x = NULL, color = "Region")

## --- Distribution by region ---
fig06b <- zp %>%
  ggplot(aes(x = region, y = log_mean_zp, fill = region)) +
  geom_boxplot() +
  facet_wrap(~season) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(title = "Zooplankton Distribution by Region",
       y = "log\u2081\u2080 annual mean ZP", x = NULL) +
  theme(legend.position = "none")

print(fig06a)
print(fig06b)
# save_fig(fig06a, "fig06a_zp_regional_ts.png")
# save_fig(fig06b, "fig06b_zp_regional_dist.png")

## ------------------------------------------ ##
#   FIGURE 7: MAB vs SNE Zooplankton Comparison
## ------------------------------------------ ##

zp_raw <- read_csv(here::here("data", "raw", "EcoMon_plankton_v3_10.csv")) %>%
  clean_names() %>%
  select(region, season, year, volume_1m2) %>%
  filter(region %in% c("MAB", "SNE", "GOM", "GB"), year >= 1998) %>%
  drop_na(volume_1m2) %>%
  group_by(region, season) %>%
  mutate(
    min_nonzero = min(volume_1m2[volume_1m2 > 0], na.rm = TRUE),
    log10_zp    = log10(volume_1m2 + min_nonzero / 2)
  ) %>%
  group_by(region, season, year) %>%
  summarize(log_mean_zp = mean(log10_zp, na.rm = TRUE), .groups = "drop")

## --- SD by region (MAB and SNE separate) ---
zp_sd_compare <- zp_raw %>%
  arrange(region, season, year) %>%
  group_by(region, season) %>%
  summarize(sd_zp = sd(log_mean_zp, na.rm = TRUE), .groups = "drop")

print(zp_sd_compare)

## --- Correlation MAB vs SNE ---
zp_cor <- zp_raw %>%
  filter(region %in% c("MAB", "SNE")) %>%
  pivot_wider(names_from = region, values_from = log_mean_zp) %>%
  group_by(season) %>%
  summarize(r = cor(MAB, SNE, use = "complete.obs"))

print(zp_cor)

## --- Plot: time series MAB vs SNE ---
fig07a <- zp_raw %>%
  filter(region %in% c("MAB", "SNE")) %>%
  ggplot(aes(x = year, y = log_mean_zp, color = region)) +
  geom_line(linewidth = 1) +
  geom_point() +
  facet_wrap(~season) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  labs(title    = "Zooplankton: MAB vs SNE (1998–present)",
       subtitle = paste("Pearson r by season:",
                        zp_cor %>%
                          mutate(label = paste0(season, ": r = ", 
                                                round(r, 2))) %>%
                          pull(label) %>% paste(collapse = " | ")),
       y = expression(log[10]~"annual mean abundance (per 10 m"^2*")"),
       x = NULL, color = "Region")

## --- Plot: boxplot all 4 regions ---
fig07b <- zp_raw %>%
  ggplot(aes(x = region, y = log_mean_zp, fill = region)) +
  geom_boxplot() +
  facet_wrap(~season) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(title = "Zooplankton Distribution by Region (1998–present)",
       y = expression(log[10]~"annual mean abundance (per 10 m"^2*")"),
       x = NULL) +
  theme(legend.position = "none")

print(fig07a)
print(fig07b)
# save_fig(fig07a, "fig07a_zp_mab_sne_ts.png")
# save_fig(fig07b, "fig07b_zp_regional_boxplot.png")

rm(zp_raw, zp_cor, zp_sd_compare)

## ------------------------------------------ ##
#   FIGURE 8: MAB vs SNE Ichthyo Comparison
## ------------------------------------------ ##

ip_raw <- read_csv(here::here("data", "raw", "EcoMon_plankton_v3_10.csv")) %>%
  clean_names() %>%
  select(region, season, year,
         cluhar_10m2, scosco_10m2, pepspp_10m2, ammspp_10m2) %>%
  filter(region %in% c("MAB", "SNE", "GOM", "GB"), year >= 1998) %>%
  mutate(ichthyosum = cluhar_10m2 + scosco_10m2 + pepspp_10m2 + ammspp_10m2) %>%
  drop_na(ichthyosum) %>%
  group_by(region, season) %>%
  mutate(
    min_nonzero   = min(ichthyosum[ichthyosum > 0], na.rm = TRUE),
    log10_ichthyo = log10(ichthyosum + min_nonzero / 2)
  ) %>%
  group_by(region, season, year) %>%
  summarize(log_mean_ichthyo = mean(log10_ichthyo, na.rm = TRUE),
            .groups = "drop")

## --- SD by region (MAB and SNE separate) ---
ip_sd_compare <- ip_raw %>%
  arrange(region, season, year) %>%
  group_by(region, season) %>%
  summarize(sd_ichthyo = sd(log_mean_ichthyo, na.rm = TRUE), .groups = "drop")
print(ip_sd_compare)

## --- Correlation MAB vs SNE ---
ip_cor <- ip_raw %>%
  filter(region %in% c("MAB", "SNE")) %>%
  pivot_wider(names_from = region, values_from = log_mean_ichthyo) %>%
  group_by(season) %>%
  summarize(r = cor(MAB, SNE, use = "complete.obs"))
print(ip_cor)

## --- Plot: time series MAB vs SNE ---
fig08a <- ip_raw %>%
  filter(region %in% c("MAB", "SNE")) %>%
  ggplot(aes(x = year, y = log_mean_ichthyo, color = region)) +
  geom_line(linewidth = 1) +
  geom_point() +
  facet_wrap(~season) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  labs(title    = "Ichthyoplankton: MAB vs SNE (1998–2023)",
       subtitle = paste("Pearson r by season:",
                        ip_cor %>%
                          mutate(label = paste0(season, ": r = ",
                                                round(r, 2))) %>%
                          pull(label) %>% paste(collapse = " | ")),
       y = expression(log[10]~"annual mean abundance (per 10 m"^2*")"), 
       x = NULL, color = "Region")

## --- Plot: boxplot all 4 regions ---
fig08b <- ip_raw %>%
  ggplot(aes(x = region, y = log_mean_ichthyo, fill = region)) +
  geom_boxplot() +
  facet_wrap(~season) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(title = "Ichthyoplankton Distribution by Region (1998–2023)",
       y = expression(log[10]~"annual mean abundance (per 10 m"^2*")"), 
       x = NULL) +
  theme(legend.position = "none")

print(fig08a)
print(fig08b)

# save_fig(fig08a, "fig08a_ip_mab_sne_ts.png")
# save_fig(fig08b, "fig08b_ip_regional_boxplot.png")

rm(ip_raw, ip_cor, ip_sd_compare)

## ------------------------------------------ ##
#   FIGURE 9: ICHTHYO by taxa
## ------------------------------------------ ##

ip_taxa <- read_csv(here::here("data", "raw", "EcoMon_plankton_v3_10.csv")) %>%
  clean_names() %>%
  select(region, season, year,
         cluhar_10m2, scosco_10m2, pepspp_10m2, ammspp_10m2) %>%
  filter(region %in% c("MAB", "SNE", "GOM", "GB"), year >= 1998) %>%
  pivot_longer(cols = c(cluhar_10m2, scosco_10m2, pepspp_10m2, ammspp_10m2),
               names_to = "taxon", values_to = "abundance") %>%
  drop_na(abundance) %>%
  group_by(region, season, taxon) %>%
  mutate(
    min_nonzero = min(abundance[abundance > 0], na.rm = TRUE),
    log10_abund = log10(abundance + min_nonzero / 2)
  ) %>%
  group_by(region, season, year, taxon) %>%
  summarize(log_mean = mean(log10_abund, na.rm = TRUE), .groups = "drop") %>%
  mutate(taxon = recode(taxon,
                        "cluhar_10m2" = "Atlantic herring",
                        "scosco_10m2" = "Atlantic mackerel",
                        "pepspp_10m2" = "Butterfish",
                        "ammspp_10m2" = "Sand lance"
  ))

## --- taxa time series MAB vs SNE ---
fig09a <- ip_taxa %>%
  filter(region %in% c("MAB", "SNE")) %>%
  ggplot(aes(x = year, y = log_mean, color = region)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 0.9) +
  facet_grid(taxon ~ season) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  labs(title = "Ichthyoplankton by Taxa: MAB vs SNE (1998–2023)",
       y = expression(log[10]~"annual mean abundance (per 10 m"^2*")"),
       x = NULL, color = "Region")

## --- taxa boxplot all 4 regions ---
fig09b <- ip_taxa %>%
  ggplot(aes(x = region, y = log_mean, fill = region)) +
  geom_boxplot() +
  facet_grid(taxon ~ season) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(title = "Ichthyoplankton Taxa Distribution by Region (1998–2023)",
       y = expression(log[10]~"annual mean abundance (per 10 m"^2*")"),
       x = NULL) +
  theme(legend.position = "none")

print(fig09a)
print(fig09b)

# save_fig(fig09a, "fig09a_ip_bytaxa_allregions_ts.png")
# save_fig(fig09b, "fig09b_ip_bytaxa_allregions_boxplot.png")

rm(ip_raw, ip_cor, ip_sd_compare, ip_taxa)

## ------------------------------------------ ##
#   TABLE 1: SD Summary — all three trophic levels side by side
## ------------------------------------------ ##

sd_summary <- zp %>%
  distinct(region, season, sd_zp) %>%
  rename(SD_zooplankton = sd_zp) %>%
  left_join(
    ichthyo %>% distinct(region, season, sd_ichthyo) %>%
      rename(SD_ichthyoplankton = sd_ichthyo),
    by = c("region", "season")
  ) %>%
  left_join(
    ff %>% distinct(region, season, sd_ff) %>%
      rename(SD_foragefish = sd_ff),
    by = c("region", "season")
  ) %>%
  arrange(region, season)

print(sd_summary)

sd_summary %>%
  mutate(region = factor(region,
                         levels = c("GOM", "GB", "MAB"),
                         labels = c("Gulf of Maine", "Georges Bank", 
                                    "Mid-Atlantic Bight"))) %>%
  arrange(region, season) %>%
  gt(rowname_col = "season", groupname_col = "region") %>%
  tab_header(
    title    = "Interannual Variability by Trophic Level",
    subtitle = "SD of log\u2081\u2080 annual mean: 1998\u20132023"
  ) %>%
  cols_label(
    SD_zooplankton     = "Zooplankton",
    SD_ichthyoplankton = "Ichthyoplankton",
    SD_foragefish      = "Forage Fish"
  ) %>%
  fmt_number(columns = starts_with("SD_"), decimals = 3) %>%
  data_color(columns = starts_with("SD_"), palette = "Blues") %>%
  tab_style(
    style    = cell_text(align = "center"),
    locations = cells_row_groups()
  ) %>%
  tab_footnote(
    "SD computed across log\u2081\u2080 annual means within each region \u00D7 season group"
  ) %>%
  tab_source_note(
    "ZP & Ichthyo: EcoMon | Forage fish: NEFSC BTS"
  )


## ------------------------------------------ ##
#   FIGURE 10: SD across Trophic Levels
## ------------------------------------------ ##

sd_long <- sd_summary %>%
  pivot_longer(
    cols      = c(SD_zooplankton, SD_ichthyoplankton, SD_foragefish),
    names_to  = "trophic_level",
    values_to = "sd"
  ) %>%
  mutate(
    trophic_level = factor(trophic_level,
                           levels = c("SD_zooplankton",
                                      "SD_ichthyoplankton",
                                      "SD_foragefish"),
                           labels = c("Zooplankton",
                                      "Ichthyoplankton",
                                      "Forage Fish")),
    region = factor(region,
                    levels = c("GOM", "GB", "MAB"),
                    labels = c("Gulf of Maine", "Georges Bank",
                               "Mid-Atlantic Bight")),
    season = factor(season,
                    levels = c("Spring", "Summer", "Fall", "Winter"))
  )

fig10 <- ggplot(sd_long,
                aes(x = trophic_level, y = sd,
                    color = season, shape = season,
                    group = season)) +
  geom_point(size = 4, alpha = 0.9) +
  #geom_line(linewidth = 0.8, alpha = 0.7) +
  facet_wrap(~region, ncol = 3) +
  scale_color_brewer(palette = "Set2") +
  scale_shape_manual(values = c("Spring" = 16, "Summer" = 17,
                                "Fall"   = 15, "Winter" = 18)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    title    = "SD Across Trophic Levels (1998–2023)",
    x        = NULL,     y     = bquote("SD of" ~ log[10] ~ "annual mean"),
    color    = "Season", shape = "Season"
  ) +
  theme_bw() +
  theme(
    strip.text   = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 20, hjust = 1, color = "black"),
    panel.grid.x = element_blank(),
    legend.position = "bottom"
  )

print(fig10)
# save_fig(fig10, "fig10_sd_trophic_levels.png", width = 11, height = 5)
