library(tidyverse)
library(patchwork)
library(gt)
library(sf)

## ============================================================ ##
##  FIGURE 1: Station Coverage Map (TS2)
## ============================================================ ##

# combine all three trophic levels with a source label
map_zp      <- ts2_zp      %>% select(lon, lat, region, season) %>% mutate(source = "Zooplankton")
map_ichthyo <- ts2_ichthyo %>% select(lon, lat, region, season) %>% mutate(source = "Ichthyoplankton")
map_ff      <- ts2_ff      %>% select(lon, lat, region, season) %>% mutate(source = "Forage fish")

map_all <- bind_rows(map_zp, map_ichthyo, map_ff) %>%
  mutate(source = factor(source, levels = c("Zooplankton", "Ichthyoplankton", "Forage fish")))

fig_map <- ggplot() +
  geom_sf(data = ecomap_valid, fill = "gray92", color = "gray60", linewidth = 0.3) +
  borders("world", colour = "gray40", fill = "gray80") +
  geom_point(data = map_all, aes(x = lon, y = lat, color = region),
             alpha = 0.3, size = 0.8) +
  coord_sf(xlim = c(-77, -61), ylim = c(35, 46)) +
  facet_wrap(~source, ncol = 3) +
  scale_color_brewer(palette = "Set1", name = "Region") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  labs(title = "Station Coverage — TS2 (1998–present)",
       x = NULL, y = NULL) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold"),
        axis.text  = element_text(size = 7),
        legend.position = "bottom")

print(fig_map)


## ============================================================ ##
##  FIGURE 2: Raw vs. Corrected CPUE (Forage Fish)
##  Shows impact of vessel/gear corrections around 2009
## ============================================================ ##

# you'll need to re-run prep_survey before corrections for the "raw" version
# if you don't have it saved, here's how to reconstruct it:
fall_raw   <- prep_survey(fall_catch,   fall_meta,   "Fall")
spring_raw <- prep_survey(spring_catch, spring_meta, "Spring")

raw_ff <- bind_rows(fall_raw, spring_raw) %>%
  filter(svspp %in% forage_codes) %>%
  left_join(forage_spp %>% select(svspp, target_taxa), by = "svspp") %>%
  group_by(target_taxa, year) %>%
  summarise(wtcpue = mean(wtcpue, na.rm = TRUE), .groups = "drop") %>%
  mutate(version = "Raw")

corrected_ff <- bind_rows(fall_df, spring_df) %>%
  filter(svspp %in% forage_codes) %>%
  left_join(forage_spp %>% select(svspp, target_taxa), by = "svspp") %>%
  group_by(target_taxa, year) %>%
  summarise(wtcpue = mean(wtcpue, na.rm = TRUE), .groups = "drop") %>%
  mutate(version = "Corrected")

compare_ff <- bind_rows(raw_ff, corrected_ff) %>%
  filter(year >= 1998)

fig_corrections <- ggplot(compare_ff, aes(x = year, y = wtcpue,
                                          color = version, linetype = version)) +
  geom_vline(xintercept = 2008.5, linetype = "dashed",
             color = "gray40", linewidth = 0.6) +
  annotate("text", x = 2009.2, y = Inf, label = "Bigelow\nintroduced",
           hjust = 0, vjust = 1.3, size = 2.8, color = "gray30") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2) +
  scale_color_manual(values = c("Raw" = "tomato", "Corrected" = "steelblue")) +
  scale_linetype_manual(values = c("Raw" = "dashed", "Corrected" = "solid")) +
  facet_wrap(~target_taxa, scales = "free_y", ncol = 2) +
  labs(title    = "Effect of Vessel & Gear Corrections on Forage Fish CPUE",
       subtitle = "TS2 (1998–present) | Dashed line = Bigelow transition (2009)",
       y = "Mean CPUE (kg tow\u207B\u00B9)",
       color = NULL, linetype = NULL) +
  theme_bw() +
  theme(strip.text      = element_text(face = "bold"),
        legend.position = "bottom")

print(fig_corrections)

raw_ff <- bind_rows(fall_raw, spring_raw) %>%
  filter(svspp %in% forage_codes) %>%
  left_join(forage_spp %>% select(svspp, target_taxa), by = "svspp") %>%
  group_by(target_taxa, season_survey, year) %>%                    # added season_survey
  summarise(wtcpue = mean(wtcpue, na.rm = TRUE), .groups = "drop") %>%
  mutate(version = "Raw")

corrected_ff <- bind_rows(fall_df, spring_df) %>%
  filter(svspp %in% forage_codes) %>%
  left_join(forage_spp %>% select(svspp, target_taxa), by = "svspp") %>%
  group_by(target_taxa, season_survey, year) %>%                    # added season_survey
  summarise(wtcpue = mean(wtcpue, na.rm = TRUE), .groups = "drop") %>%
  mutate(version = "Corrected")

compare_ff <- bind_rows(raw_ff, corrected_ff) %>%
  filter(year >= 1998)

fig_corrections_season <- ggplot(compare_ff, aes(x = year, y = wtcpue,
                                                 color = version, linetype = version)) +
  geom_vline(xintercept = 2008.5, linetype = "dashed",
             color = "gray40", linewidth = 0.5) +
  annotate("text", x = 2009.3, y = Inf, label = "Bigelow",
           hjust = 0, vjust = 1.4, size = 2.5, color = "gray30") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1) +
  scale_color_manual(values = c("Raw" = "tomato", "Corrected" = "steelblue")) +
  scale_linetype_manual(values = c("Raw" = "dashed", "Corrected" = "solid")) +
  facet_grid(season_survey ~ target_taxa, scales = "free_y") +    # season as rows, species as cols
  labs(title    = "Effect of Vessel & Gear Corrections on Forage Fish CPUE",
       subtitle = "TS2 (1998–present) | Dashed line = Bigelow transition (2009)",
       x = "Year",
       y = "Mean CPUE (kg tow\u207B\u00B9)",
       color = NULL, linetype = NULL) +
  theme_bw() +
  theme(strip.text      = element_text(face = "bold"),
        legend.position = "bottom",
        axis.text.x     = element_text(angle = 45, hjust = 1, size = 7))

print(fig_corrections_season)
## ============================================================ ##
##  FIGURE 3: Annual Log Means + 5-yr Running Mean (TS2)
##  Clean version of what's already in your scripts
## ============================================================ ##

plot_ts <- function(sum_df, rm_df, y_col, rm_col, title_label, y_label) {
  ggplot(sum_df, aes(x = year, y = .data[[y_col]])) +
    geom_point(size = 1.5, alpha = 0.7, color = "gray30") +
    geom_line(data = rm_df, aes(x = year, y = .data[[rm_col]]),
              color = "tomato", linewidth = 1.2) +
    facet_grid(season ~ region) +
    labs(title    = title_label,
         subtitle = "Points = annual log mean | Red line = 5-year running mean",
         y = y_label, x = NULL) +
    theme_bw() +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
          strip.text   = element_text(face = "bold"),
          panel.grid   = element_blank())
}

fig_ts_zp <- plot_ts(ts2_zp_sum, ts2_rm,
                     "log_mean_zp",      "run_mean_zp",
                     "Zooplankton Displacement Volume — TS2 (1998–present)",
                     "log\u2081\u2080 ZP (mL m\u207B\u00B2)")

fig_ts_ichthyo <- plot_ts(ts2_ichthyo_sum, ts2_rm,
                          "log_mean_ichthyo", "run_mean_ichthyo",
                          "Ichthyoplankton — TS2 (1998–present)",
                          "log\u2081\u2080 larval abundance (ind. 10 m\u207B\u00B2)")

fig_ts_ff <- plot_ts(ts2_ff_sum, ts2_rm_ff,
                     "log_mean_ff",      "run_mean_ff",
                     "Forage Fish CPUE — TS2 (1998–present)",
                     "log\u2081\u2080 CPUE (kg tow\u207B\u00B9)")

print(fig_ts_zp)
print(fig_ts_ichthyo)
print(fig_ts_ff)


## ============================================================ ##
##  FIGURE 4: Anomaly Bar Plots (TS2)
## ============================================================ ##

plot_anom <- function(sum_df, y_col, anom_col, title_label, y_label) {
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
         y = y_label, x = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          strip.text  = element_text(face = "bold"),
          panel.grid  = element_blank())
}

fig_anom_zp <- plot_anom(ts2_zp_sum, "log_mean_zp", "anom",
                         "Zooplankton Anomaly — TS2",
                         "Anomaly in log\u2081\u2080 annual mean ZP")

fig_anom_ichthyo <- plot_anom(ts2_ichthyo_sum, "log_mean_ichthyo", "anom",
                              "Ichthyoplankton Anomaly — TS2",
                              "Anomaly in log\u2081\u2080 annual mean ichthyo")

fig_anom_ff <- plot_anom(ts2_ff_sum, "log_mean_ff", "anom",
                         "Forage Fish Anomaly — TS2",
                         "Anomaly in log\u2081\u2080 annual mean CPUE")

print(fig_anom_zp)
print(fig_anom_ichthyo)
print(fig_anom_ff)


## ============================================================ ##
##  FIGURE 5: N Observations in 5-Year Running Mean Window
## ============================================================ ##

plot_nwindow <- function(check_df, title_label) {
  ggplot(check_df, aes(x = year, y = n_in_window)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_hline(yintercept = 5, linetype = "dashed", color = "tomato", linewidth = 0.7) +
    facet_grid(season ~ region) +
    labs(title    = title_label,
         subtitle = "Dashed line = full 5-year window",
         y = "N years in window", x = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          strip.text  = element_text(face = "bold"),
          panel.grid  = element_blank())
}

# rebuild ts2_check for each trophic level if needed
make_nwindow <- function(sum_df, y_col) {
  sum_df %>%
    arrange(region, season, year) %>%
    group_by(region, season) %>%
    mutate(n_in_window = runner::runner(
      x   = .data[[y_col]],
      k   = 5,
      idx = year,
      f   = \(x) sum(!is.na(x))
    )) %>%
    ungroup()
}

ts2_check_zp      <- make_nwindow(ts2_zp_sum,      "log_mean_zp")
ts2_check_ichthyo <- make_nwindow(ts2_ichthyo_sum,  "log_mean_ichthyo")
ts2_check_ff      <- make_nwindow(ts2_ff_sum,       "log_mean_ff")

fig_nwin_zp      <- plot_nwindow(ts2_check_zp,      "Zooplankton — N years in 5-yr window")
fig_nwin_ichthyo <- plot_nwindow(ts2_check_ichthyo,  "Ichthyoplankton — N years in 5-yr window")
fig_nwin_ff      <- plot_nwindow(ts2_check_ff,       "Forage Fish — N years in 5-yr window")

print(fig_nwin_zp)
print(fig_nwin_ichthyo)
print(fig_nwin_ff)


## ============================================================ ##
##  TABLE: SD Summary (all three trophic levels side by side)
## ============================================================ ##

sd_zp_tbl <- ts2_zp %>%
  distinct(region, season, sd_zp) %>%
  rename(SD_zooplankton = sd_zp)

sd_ichthyo_tbl <- ts2_ichthyo %>%
  distinct(region, season, sd_ichthyo) %>%
  rename(SD_ichthyoplankton = sd_ichthyo)

sd_ff_tbl <- ts2_ff %>%
  distinct(region, season, sd_ff) %>%
  rename(SD_foragefish = sd_ff)

sd_summary <- sd_zp_tbl %>%
  left_join(sd_ichthyo_tbl, by = c("region", "season")) %>%
  left_join(sd_ff_tbl,      by = c("region", "season")) %>%
  arrange(region, season)

# print to console
print(sd_summary)

# formatted gt table
sd_summary %>%
  gt(rowname_col = "season", groupname_col = "region") %>%
  tab_header(
    title    = "Interannual Variability by Trophic Level",
    subtitle = "SD of log\u2081\u2080 annual mean CPUE — TS2 (1998\u2013present)"
  ) %>%
  cols_label(
    SD_zooplankton     = "Zooplankton",
    SD_ichthyoplankton = "Ichthyoplankton",
    SD_foragefish      = "Forage Fish"
  ) %>%
  fmt_number(columns = starts_with("SD_"), decimals = 3) %>%
  data_color(
    columns  = starts_with("SD_"),
    palette  = "Blues"
  ) %>%
  tab_footnote("SD computed across annual log\u2081\u2080 means within each region \u00D7 season group") %>%
  tab_source_note("TS2: 1998\u2013present | ZP & Ichthyo: EcoMon | Forage fish: NEFSC BTS")




ff_plot <- ff %>%
  filter(svspp %in% forage_codes) %>%
  left_join(forage_spp %>% select(svspp, target_taxa), 
            by = c("svspp","target_taxa")) %>%
  group_by(target_taxa, year, season, region) %>%
  summarise(wtcpue = mean(wtcpue, na.rm = TRUE), .groups = "drop") %>%
  group_by(target_taxa, season, region) %>%          # z-score within each taxa/season/region combo
  mutate(
    anomaly = (wtcpue - mean(wtcpue, na.rm = TRUE)) / sd(wtcpue, na.rm = TRUE)
  ) %>%
  ungroup()

plot_species <- function(sp_name) {
  df <- ff_plot %>% filter(target_taxa == sp_name)
  
  ggplot(df, aes(x = year, y = anomaly)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = pmin(anomaly, 0), ymax = 0), fill = "steelblue", alpha = 0.4) +
    geom_ribbon(aes(ymin = 0, ymax = pmax(anomaly, 0)), fill = "tomato", alpha = 0.4) +
    geom_line(color = "gray30", linewidth = 0.5) +
    facet_grid(season ~ region) +
    labs(
      title = paste("CPUE Anomaly:", sp_name),
      subtitle = "Standardized (z-score) relative to full time series mean",
      x = "Year",
      y = "Anomaly (SD)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold")
    )
}

for (sp in unique(ff_plot$target_taxa)) {
  print(plot_species(sp))
}






ff_plot <- ff %>%
  filter(svspp %in% forage_codes) %>%
  left_join(forage_spp %>% select(svspp, target_taxa), 
            by = c("svspp","target_taxa")) %>%
  group_by(target_taxa, year, season, region) %>%
  summarise(wtcpue = mean(wtcpue, na.rm = TRUE), .groups = "drop") %>%
  group_by(target_taxa, season, region) %>%
  mutate(
    log_wtcpue = log1p(wtcpue),
    anomaly_zscore = (wtcpue - mean(wtcpue, na.rm = TRUE)) / sd(wtcpue, na.rm = TRUE),
    anomaly_log    = (log_wtcpue - mean(log_wtcpue, na.rm = TRUE)) / sd(log_wtcpue, na.rm = TRUE)
  ) %>%
  ungroup()

plot_species <- function(sp_name, anomaly_type = c("zscore", "log")) {
  anomaly_type <- match.arg(anomaly_type)
  
  df <- ff_plot %>%
    filter(target_taxa == sp_name) %>%
    mutate(anomaly = if(anomaly_type == "zscore") anomaly_zscore else anomaly_log)  # base R if/else
  
  subtitle <- if(anomaly_type == "zscore") {
    "Standardized z-score relative to full time series mean"
  } else {
    "Log(x+1) standardized z-score relative to full time series mean"
  }
  
  ggplot(df, aes(x = year, y = anomaly)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = pmin(anomaly, 0), ymax = 0), fill = "steelblue", alpha = 0.4) +
    geom_ribbon(aes(ymin = 0, ymax = pmax(anomaly, 0)), fill = "tomato", alpha = 0.4) +
    geom_line(color = "gray30", linewidth = 0.5) +
    facet_grid(season ~ region) +
    labs(
      title    = paste("CPUE Anomaly:", sp_name),
      subtitle = subtitle,
      x = "Year",
      y = "Anomaly (SD)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text  = element_text(face = "bold")
    )
}

# compare both side by side for one species
library(patchwork)
plot_species("Atlantic herring", "zscore") + plot_species("Atlantic herring", "log")

# or loop through all species, printing both versions
for (sp in unique(ff_plot$target_taxa)) {
  print(plot_species(sp, "zscore") + plot_species(sp, "log") +
          plot_annotation(title = sp))
}

for (sp in unique(ff_plot$target_taxa)) {
  
  p <- plot_species(sp, "zscore") +
    plot_species(sp, "log") +
    plot_annotation(title = sp)
  
  print(p)
  
  ggsave(
    filename = paste0("figures/ff_", sp, ".png"),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
  )
}




sandlance_plot <- ff %>%
  filter(svspp %in% c("181", "734")) %>%
  mutate(species_name = case_when(
    svspp == "181" ~ "Northern sand lance (A. dubius)",
    svspp == "734" ~ "American sand lance (A. americanus)"
  )) %>%
  group_by(species_name, year, season, region) %>%
  summarise(wtcpue = mean(wtcpue, na.rm = TRUE), .groups = "drop") %>%
  group_by(species_name, season, region) %>%
  mutate(
    log_wtcpue     = log1p(wtcpue),
    anomaly_zscore = (wtcpue - mean(wtcpue, na.rm = TRUE)) / sd(wtcpue, na.rm = TRUE),
    anomaly_log    = (log_wtcpue - mean(log_wtcpue, na.rm = TRUE)) / sd(log_wtcpue, na.rm = TRUE)
  ) %>%
  ungroup()






a_am <- ff %>%
  filter(svspp == "734") %>%
  mutate(species_name = "American sand lance (A. americanus)") %>%
  group_by(year, season, region) %>%
  summarise(
    wtcpue = mean(wtcpue, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(region, season, year)

ggplot(a_am, aes(x = year, y = wtcpue)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_line(linewidth = 0.4, alpha = 0.7) +
  facet_grid(season ~ region, scales = "free_y") +
  labs(
    title = "American sand lance (A. americanus)",
    subtitle = "Raw annual mean WTCPUE",
    x = "Year",
    y = "WTCPUE"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )
