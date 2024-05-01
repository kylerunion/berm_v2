#figures
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(ggpubr)
library(RColorBrewer)
library(ggh4x)
library(ggpattern)
library(terra)
library(ggpubr)
library(forcats)
library(patchwork)
library(ggpmisc)
library(rstatix)
library(moments)
library(Hmisc)
library(car)
library(gridExtra)
library(ggtext)
library(ggthemes)
library(lubridate)
library(rstatix)
library(viridis)
library(ggsci)



#figure 2
veg<-read.csv("output/core_vegplot_combined.csv", header=T)
veg$date <- as.Date(veg$date)
veg <- veg[veg$date > as.Date("2021-01-01"),]
veg <- veg[veg$site != "fluxa",]
veg <- veg %>%
  drop_na(bgm2.core.stemscaled) %>%
  mutate(mo = month(date), yr = year(date)) %>%
  dplyr::select(mo, yr, site, bgm2.core.stemscaled, elevation) %>%
  mutate(elevation_group = ifelse(elevation > summary(veg$elevation)[5], "H", NA)) %>%
  mutate(elevation_group = ifelse(elevation <= summary(veg$elevation)[5], "M", elevation_group)) %>%
  mutate(elevation_group = ifelse(elevation <= summary(veg$elevation)[2], "L", elevation_group)) %>%
  mutate(site_group = ifelse(site %in% c("folly"), "b", "a")) %>%
  mutate(month_group = ifelse(mo == 2, 1, NA)) %>%
  mutate(month_group = ifelse(mo %in% c(5,6), 2, month_group)) %>%
  mutate(month_group = ifelse(mo == 8, 3, month_group)) %>%
  mutate(month_group = ifelse(mo == 11, 4, month_group)) %>%
  group_by(site) %>%
  mutate(norm_bgb = bgm2.core.stemscaled/max(bgm2.core.stemscaled)) %>%
  ungroup()
table(veg$elevation_group)

library(ggthemes)
plot_theme <- function(){ 
  font <- "Georgia"   #assign font family up front
  theme_hc() %+replace%    #replace elements we want to change
    theme(
      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      axis.ticks = element_blank(),          #strip axis ticks
      #text elements
      plot.title = element_text(family = font, size = 10, face = 'bold', hjust = 0.5, vjust = 2),
      axis.title = element_text(family = font, size = 10),               
      axis.text = element_text(family = font, size = 7), 
      legend.title = element_text(family = font, size = 10),
      legend.key.width = unit(15,"pt"),
      legend.text = element_text(family = font, size = 6),
      legend.box="vertical",
      panel.background = element_rect(fill = "white", color = "white")
    )
}

bgb_seasons <- ggplot(data = veg, aes(x = month_group, y = bgm2.core.stemscaled)) +
  geom_jitter(aes(col = factor(yr), shape = factor(site)), height = 0, width = 0.2) +
  geom_smooth(color = "black") +
  scale_x_continuous(breaks = c(1,2,3,4), labels = c("Feb", "May-June", "Aug", "Nov"), limits = c(1,4)) +
  coord_cartesian(ylim = c(0,2000)) +
  #scale_y_continuous(limits = c(0,2000)) +
  scale_shape_manual(values=1:nlevels(as.factor(veg$site)), labels = c(deanc = "Dean Creek", dupli = "Duplin River", fluxb = "Flux B", folly = "Folly River", huntc = "Hunt Camp", north = "North Sapelo", skida = "Skidaway", ugami = "UGAMI")) +
  scale_color_viridis_d(breaks = c(2021,2022,2023)) +
  labs(color = "Year", shape = "Site", x = "Sampling Event", y = expression(paste("Belowground Biomass (g m"^-2~")"))) +
  guides(shape = guide_legend(nrow = 9), color = guide_legend(nrow = 3)) +
  plot_theme() + 
  theme(legend.position = "right")
bgb_seasons
ggsave("results/bgb_seasons.jpg", bgb_seasons, width = 4, height = 4, dpi = 600)

#figure 3
bg_all_tr_test_results <- read.csv("results/bg_all_tr_test_results_v2.0.csv")
bg_all_tr_test_results$version <- "2.0"
bg_tr_test_results <- bg_all_tr_test_results %>%
  filter(run == 1) %>%
  dplyr::select(-run)
flatt = pal_flatui(palette = "flattastic")(9)
#p-value matrix among groups
df <- bg_tr_test_results
df2 <- veg
df2$year <- year(df2$date)
df2 <- df2 %>%
  filter(year > 2015)

bgb_site_density <- df2 %>%
  group_by(site) %>%
  ggplot(aes(x = bgm2.core.stemscaled, y = reorder(site,bgm2.core.stemscaled), fill = site)) +
  ggridges::stat_density_ridges() +  
  scale_fill_manual(name = "Site", values = flatt) +
  labs(x = expression("Belowground Biomass (g " ~ m^-~2 ~ ")"), y = "Observation Density", title = "Site") +
  theme(
    text = element_text(size = 14, family = "Georgia"),
    plot.title = element_text(size = 14, face = 'bold', hjust = 0.5, vjust = 2),
    legend.position = "none",
    legend.spacing.y = unit(0, 'in'),
    legend.key = element_rect(fill = "white"),
    legend.key.size = unit(0.15, 'in'),
    legend.text = element_text(size = 8),
    panel.background = element_rect(fill = "white", color = "white"),
    axis.text.x = element_text(angle = 45, vjust = 0.5)) 
bgb_site_density

bgb_year_density <- df2 %>%
  filter(year > 2015) %>%
  mutate(year = factor(year)) %>%
  group_by(year) %>%
  ggplot(aes(x = bgm2.core.stemscaled, y = reorder(year,bgm2.core.stemscaled), fill = year)) +
  ggridges::stat_density_ridges() +  
  scale_fill_manual(name = "Year", values = years_col) +
  labs(x = expression("Belowground Biomass (g " ~ m^-~2~")"), y = "Observation Density", title = "Year") +
  theme(
    text = element_text(size = 14, family = "Georgia"),
    plot.title = element_text(size = 14, face = 'bold', hjust = 0.5, vjust = 2),
    legend.position = "none",
    legend.spacing.y = unit(0, 'in'),
    legend.key = element_rect(fill = "white"),
    legend.key.size = unit(0.15, 'in'),
    legend.text = element_text(size = 8),
    panel.background = element_rect(fill = "white", color = "white"),
    axis.text.x = element_text(angle = 45, vjust = 0.5)) 
bgb_year_density

library(tidyverse)

p_pal <- c(" " = "#F6EDBD", "*" = "#EDBB8A", "**" = "#DE8A5A", "***" = "#CA562C")
signif_size = 3.5

bgb_site_t <- pairwise.wilcox.test(df2$bgm2.core.stemscaled, df2$site, p.adjust.method = "bonf", paired = F)
bgb_site_t_p <- as.data.frame(bgb_site_t$p.value)

dat2 <-
  bgb_site_t_p %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column("Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>%
  mutate(signif = ifelse(value < 0.05, "*", " ")) %>%
  mutate(signif = ifelse(value < 0.01, "**", signif)) %>%
  mutate(signif = ifelse(value < 0.001, "***", signif))

bgb_site_t_plot <- ggplot(data = dat2, aes(Var1, Var2, fill = signif)) +
  geom_tile() +
  geom_text(aes(label = signif), size = signif_size) +
  scale_fill_manual(values = p_pal, na.value = "white", na.translate = F) +
  labs(x = "Site", y = "Site", fill = "p-value") +
  theme(
    text = element_text(size = 14, family = "Georgia"),
    legend.position = "right",
    legend.spacing.y = unit(5, 'pt'),
    legend.key.size = unit(0.15, 'in'),
    panel.background = element_rect(fill = "white", color = "white"),
    axis.text.x = element_text(angle = 45, vjust = 0.5)) 
bgb_site_t_plot

bgb_year_t <- pairwise.wilcox.test(df2$bgm2.core.stemscaled, df2$year, p.adjust.method = "bonf", paired = F)
bgb_year_t_p <- as.data.frame(bgb_year_t$p.value)

dat2 <-
  bgb_year_t_p %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column("Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>%
  mutate(signif = ifelse(value < 0.05, "*", " ")) %>%
  mutate(signif = ifelse(value < 0.01, "**", signif)) %>%
  mutate(signif = ifelse(value < 0.001, "***", signif))

p_pal <- c(" " = "#F6EDBD", "*" = "#EDBB8A", "**" = "#DE8A5A", "***" = "#CA562C")

bgb_year_t_plot <- ggplot(data = dat2, aes(Var1, Var2, fill = signif)) +
  geom_tile() +
  geom_text(aes(label = signif), size = signif_size) +
  scale_fill_manual(values = p_pal, na.value = "white", na.translate = F) +
  labs(x = "Year", y = "Year", fill = "p-value") +
  theme(
    text = element_text(size = 14, family = "Georgia"),
    legend.position = "none",
    legend.spacing.y = unit(5, 'pt'),
    legend.key = element_rect(fill = "white"),
    legend.key.size = unit(0.15, 'in'),
    legend.text = element_text(size = 12),
    panel.background = element_rect(fill = "white", color = "white"),
    axis.text.x = element_text(angle = 45, vjust = 0.5)) 
bgb_year_t_plot

site_year_diffs <-  (bgb_site_density + bgb_year_density) / (bgb_site_t_plot + bgb_year_t_plot + plot_layout(guides = "collect")) +
  plot_annotation(tag_levels = "a")
site_year_diffs
ggsave("results/site_year_diffs_obs.jpg", site_year_diffs, width = 10, height = 8)

#figure 4

v <- 2
data<-read_csv(paste0("/home/kyle/Documents/Documents/UT/Sapelo/berm_output/processed_bg", v, ".0.csv"))

#all potential predictors
#all_feats <- read_csv("output/xgb2.0_features.csv")
#potential_feats <- all_feats$potential
#potential_feats <- potential_feats[!grepl("diff|roll|lagg|lpred|start|grow", potential_feats)]
#potential_feats
#potential_feats_full <- c("Elevation", "Greenup Day of Year", "Predicted Aboveground Biomass", "Predicted Foliar Chlorophyll", "Predicted Leaf Area Index", "Predicted Foliar N", "")

length(unique(data$elevation))

biological <- c("pred|ndvi|photo")
bio_list <- colnames(data[grepl(biological, colnames(data))])
bio_list <- bio_list[!grepl("diff|roll|lagg|lpred|start|grow", bio_list)]
climate <- c("lst|doy|prcp|vp|srad|dayl|green|par|growingday|tmin|tmax")
clim_list <- colnames(data[grepl(climate, colnames(data))])
clim_list <- clim_list[!grepl("diff|roll|lagg|lpred|start|grow|sum|greenup|.avg|greenupmo|dayl|root|tmax|tmin|par", clim_list)]
clim_list <- clim_list[!clim_list == 'doy']
hydro <- c("lochi|loclo|flood|lowater|hiwater")
hyd_list <- colnames(data[grepl(hydro, colnames(data))])
hyd_list <- hyd_list[!grepl("diff|roll|lagg|lpred|start|grow|sum", hyd_list)]
hyd_list <- hyd_list[!hyd_list == 'hiwater']
hyd_list <- hyd_list[!hyd_list == 'lowater']
physical <- c("elevation")
phys_list <- colnames(data[grepl(physical, colnames(data))])
phys_list <- phys_list[!grepl("loc", phys_list)]
bg_list <- c("bgm2.core.stemscaled")
names_list <- c(bio_list, clim_list, hyd_list, phys_list)
full_names <- c("Predicted\nAGB", "Predicted\n% Foliar N", "Predicted\n% Foliar CHL", "Predicted\nLAI", "NDVI\n", "Predicted Tot.\nFoliar N", 
                "LST\n", "Greenup\nDoY", "Avg.\nPrecip.", "Avg. Solar\nRadiation", "Avg. Vapor\nPressure",
                "FF", "II", "DI",
                "Elevation\n")
name_match <- data.frame(names_list, full_names)
bio_list2 <- name_match$full_names[match(bio_list, name_match$names_list)]
clim_list2 <- name_match$full_names[match(clim_list, name_match$names_list)]
hyd_list2 <- name_match$full_names[match(hyd_list, name_match$names_list)]
phys_list2 <- name_match$full_names[match(phys_list, name_match$names_list)]

#data <- data %>%
#  rename_at(as.vector(na.omit(name_match$names_list[match(names(data), name_match$names_list)])),
#           ~as.vector(na.omit(name_match$full_names[match(names(data), name_match$names_list)])))

#label_y <- expression(atop(x = "", y = atop(x = paste("BGB"), y = paste("(g " ~ m^-2 ~ ")"))))
label_y <- "BGB"

point_size <- 0.01
point_color <- "grey30"
smooth_width <- 2
title_size <- 1.1

plot_theme <- theme(
  #axis.text.x = element_blank(),
  #axis.ticks = element_blank(),
  axis.title.x = element_text(family = "Georgia"),
  #axis.title.y= element_blank(),
  axis.title.y = element_text(family = "Georgia", angle = 90, margin = margin(0,-15,0,0)),
  axis.line = element_line(colour = "black", size = 0.25),
  legend.position = "none",
  legend.title = element_blank(),
  plot.title = element_blank(),
  plot.margin = unit(c(1,10,1,1), "pt"),
  panel.background = element_rect(fill = "white", color = "white"))


plot_list <- list()
for (i in bio_list) {
  df <- data[,c(i,"bgm2.core.stemscaled")]
  df <- df[complete.cases(df),]
  loess_data <- stats::loess(paste("bgm2.core.stemscaled ~", i), data = df, span = 0.5)
  plot_list[[i]] <-  ggplot(data = df, aes_string(x = i, y = "bgm2.core.stemscaled")) +
    geom_point(size = point_size, color = point_color) +
    geom_smooth(se = F, color = "#209058", size = smooth_width) + 
    scale_x_continuous(breaks = c(signif(min(df[1]),2), signif(max(df[1]),2)), limits = c(signif(min(df[1]),2), signif(max(df[1]),2))) +
    scale_y_continuous(breaks = c(0,4000)) +
    labs(x = name_match$full_names[match(i, name_match$names_list)], y = label_y) +
    plot_theme
  #ggsave(plot = plot_list[[i]], file = paste0("results/bg_feat_trends/", i, ".png"))
}
plot_list <- plot_list[c("predag.allom.l8", "predchl.l8", "predn.l8", "photo", "predlai.l8", "ndvi")]
b <- grid.arrange(grobs = plot_list, nrow = 6, top = grid::textGrob(expression(underline("Biological")), gp = grid::gpar(fontfamily = "Georgia", cex = title_size)))
b <- as_ggplot(b)

plot_list <- list()
for (i in clim_list) {
  df <- data[,c(i,"bgm2.core.stemscaled")]
  df <- df[complete.cases(df),]
  loess_data <- stats::loess(paste("bgm2.core.stemscaled ~", i), data = df, span = 0.5)
  plot_list[[i]] <-  ggplot(data = df, aes_string(x = i, y = "bgm2.core.stemscaled")) +
    geom_point(size = point_size, color = point_color) +
    geom_smooth(se = F, color = "#D84840", size = smooth_width) + 
    scale_x_continuous(breaks = c(signif(min(df[1]),2), signif(max(df[1]),2)), limits = c(signif(min(df[1]),2), signif(max(df[1]),2))) +
    scale_y_continuous(breaks = c(0,4000)) +
    labs(x = name_match$full_names[match(i, name_match$names_list)], y = label_y) +
    plot_theme + 
    theme(axis.title.y = element_blank(), axis.text.y = element_blank())
  #ggsave(plot = plot_list[[i]], file = paste0("results/bg_feat_trends/", i, ".png"))
}
plot_list <- plot_list[c("lst", "greendoy", "prcp_mean", "srad_mean", "vp_mean")]
c <- grid.arrange(grobs = plot_list, nrow = 6, top = grid::textGrob(expression(underline("Climatic")), gp = grid::gpar(fontfamily = "Georgia", cex = title_size)))
c <- as_ggplot(c)

plot_list <- list()
for (i in hyd_list) {
  df <- data[,c(i,"bgm2.core.stemscaled")]
  df <- df[complete.cases(df),]
  loess_data <- stats::loess(paste("bgm2.core.stemscaled ~", i), data = df, span = 0.5)
  plot_list[[i]] <-  ggplot(data = df, aes_string(x = i, y = "bgm2.core.stemscaled")) +
    geom_point(size = point_size, color = point_color) +
    geom_smooth(se = F, color = "#4080C0", size = smooth_width) + 
    scale_x_continuous(breaks = c(signif(min(df[1]),2), signif(max(df[1]),2)), limits = c(signif(min(df[1]),2), signif(max(df[1]),2))) +
    scale_y_continuous(breaks = c(0,4000)) +
    labs(x = name_match$full_names[match(i, name_match$names_list)], y = label_y) +
    plot_theme + 
    theme(axis.title.y = element_blank(), axis.text.y = element_blank())
  #ggsave(plot = plot_list[[i]], file = paste0("results/bg_feat_trends/", i, ".png"))
}
h <- grid.arrange(grobs = plot_list, nrow = 6, top = grid::textGrob(expression(underline("Hydrologic")), gp = grid::gpar(fontfamily = "Georgia", cex = title_size)))
h <- as_ggplot(h)

plot_list <- list()
for (i in phys_list) {
  df <- data[,c(i,"bgm2.core.stemscaled")]
  df <- df[complete.cases(df),]
  loess_data <- stats::loess(paste("bgm2.core.stemscaled ~", i), data = df, span = 0.5)
  plot_list[[i]] <-  ggplot(data = df, aes_string(x = i, y = "bgm2.core.stemscaled")) +
    geom_point(size = point_size, color = point_color) +
    geom_smooth(se = F, color = "#E0B860", size = smooth_width) + 
    scale_x_continuous(breaks = c(signif(min(df[1]),2), signif(max(df[1]),2)), limits = c(signif(min(df[1]),2), signif(max(df[1]),2))) +
    scale_y_continuous(breaks = c(0,4000)) +
    labs(x = name_match$full_names[match(i, name_match$names_list)], y = label_y) +
    plot_theme + 
    theme(axis.title.y = element_blank(), axis.text.y = element_blank())
  #ggsave(plot = plot_list[[i]], file = paste0("results/bg_feat_trends/", i, ".png"))
}
p <- grid.arrange(grobs = plot_list, nrow = 6, top = grid::textGrob(expression(underline("Elevation")), gp = grid::gpar(fontfamily = "Georgia", cex = title_size)))
p <- as_ggplot(p)

fig4 <- b + c + p + h + 
  plot_layout(ncol = 4)
fig4

ggsave(file = "results/bg_candidate_predictors.jpg", fig4, width = 8, height = 12, dpi = 900)

#figure 5
feats <- read.csv("output/xgb2.0_all_import_features.csv")
feats2 <- feats %>%
  group_by(features) %>%
  pivot_longer(cols = c("train.1", "train.2", "train.3", "train.4", "train.5"), names_to = "train", values_to = "mean_importance")
feats3 <- feats2 %>%
  group_by(features) %>%
  summarise(importance_mean = mean(mean_importance), importance_sd = sd(mean_importance), importance_min = min(mean_importance), importance_max = max(mean_importance))
feats3 <- feats3[order(feats3$importance_mean, decreasing = T),]
feats_list <- feats3$features
physical <- c("elevation")
biological <- c("lai|ndvi|ag|chl|photo|growN|growing|_N_")
climate <- c("lst|doy")
hydrologic <- c("lochi|loclo|flood|lowater|hiwater")
#other <- c("doy")
feats3$category <- NA
feats3$category[grepl(physical, feats3$features)] <- "Elevation"
feats3$category[grepl(biological, feats3$features)] <- "Biological"
feats3$category[grepl(climate, feats3$features)] <- "Climatic"
feats3$category[grepl(hydrologic, feats3$features)] <- "Hydrologic"
#feats3$category[grepl(other, feats3$features)] <- "other"

feats4 <- feats3 %>%
  group_by(category) %>%
  summarise(cat_sum = sum(importance_mean))


col2hex <- function(x, alpha = FALSE) {
  args <- as.data.frame(t(col2rgb(x, alpha = alpha)))
  args <- c(args, list(names = x, maxColorValue = 255))
  do.call(rgb, args)
}
col2hex("tan1")
col2hex("steelblue1")
col2hex("palegreen")
col2hex("indianred1")

feats3$features_full <- c("FF", "Rolling prev. 5 month mean of II", "Elevation", "Rolling prev. 5 month mean of DI",
                          "Rolling prev. 4 month mean of II", "Green-up DoY", "End of prev. growing season LAI", "II",
                          "Rolling prev. 5 month mean of NDVI", "End of prev. growing season % foliar N", "End of prev. growing season AGB",
                          "Rolling prev. 3 month mean of NDVI", "Rolling prev. 5 month mean of % foliar N", "Change in LAI from end of prev. 5 month",
                          "Change in CHL from end of prev. 3 month", "Rolling prev. 3 month mean of II", "Change in AGB from prev. growing season", 
                          "Change in % foliar N from prev. 3 month", "Change in AGB from end of prev. 5 month", "End of prev. growing season total foliar N",
                          "Rolling prev. 5 month mean of LST", "DI", "Rolling prev. 5 month mean of CHL", "Change in AGB from end of prev. 2 month",
                          "Prev. month CHL", "Rolling prev. 3 month mean of CHL", "Change in NDVI from end of prev. 2 month")
write_csv(feats3, "Rmarkdown/xgb_mean_import_features.csv")



feat_importance_bar <- ggplot(data = feats3, aes(x = fct_reorder(features_full, importance_mean), y = importance_mean, fill = category)) +
  geom_col() +
  coord_flip() +
  geom_errorbar(aes(ymin = importance_min, ymax = importance_max), size = 0.35, color = "gray20", width = 0.25) +
  labs(x = "Feature", y = "Feature Importance", fill = "Category") +
  scale_fill_manual(values = c("#209058", "#D84840","#E0B860", "#4080C0")) +
  theme(
    text = element_text(size = 10, family = "Georgia"),
    axis.text = element_text(size = 8, color = "black"),
    legend.position = "none", #c(0.67,0.7),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = "white"))
feat_importance_bar
ggsave("results/bg_feat_importance_bar.png", feat_importance_bar, dpi = 300, height = 12, width = 16, units = "in")

feat_importance_bar_sum <- ggplot(data = feats4, aes(x = fct_reorder(category, cat_sum), y = cat_sum, fill = category)) +
  geom_col() +
  coord_flip() +
  labs(y = "Summed\nFeature Importance", title = "Category", x = "") +
  scale_y_continuous(breaks =c(0,0.3,0.6)) +
  scale_fill_manual(values = c("#209058", "#D84840", "#E0B860", "#4080C0")) +
  theme(
    axis.text = element_text(size = 8, color = "black"),
    legend.position = "none",
    text = element_text(size = 8, family = "Georgia"),
    panel.background = element_rect(fill = "white", color = "white"))
feat_importance_bar_sum
ggsave("results/bg_feat_importance_bar_sum.png", feat_importance_bar_sum, dpi = 300, height = 12, width = 16, units = "in")

feat_importance_combined <- feat_importance_bar + inset_element(feat_importance_bar_sum, left = 0.2, bottom = 0.1, right = 1.0, top = 0.6)
feat_importance_combined
ggsave("results/bg_feat_importance_combined.jpg", feat_importance_combined, dpi = 300, height = 6, width = 6, units = "in")

#figure 6

bg_model_output <- read_csv("output/bg_model_output.csv")
bg_model_1 <- bg_model_output %>%
  dplyr::select(bgm2.core.stemscaled, train.1, pred.1, site, group, pix) %>%
  rename(train = train.1, pred = pred.1) %>%
  mutate(run = 1) %>%
  na.omit
bg_model_2 <- bg_model_output %>%
  dplyr::select(bgm2.core.stemscaled, train.2, pred.2, site, group, pix) %>%
  rename(train = train.2, pred = pred.2) %>%
  mutate(run = 2) %>%
  na.omit
bg_model_3 <- bg_model_output %>%
  dplyr::select(bgm2.core.stemscaled, train.3, pred.3, site, group, pix) %>%
  rename(train = train.3, pred = pred.3) %>%
  mutate(run = 3) %>%
  na.omit
bg_model_4 <- bg_model_output %>%
  dplyr::select(bgm2.core.stemscaled, train.4, pred.4, site, group, pix) %>%
  rename(train = train.4, pred = pred.4) %>%
  mutate(run = 4) %>%
  na.omit
bg_model_5 <- bg_model_output %>%
  dplyr::select(bgm2.core.stemscaled, train.5, pred.5, site, group, pix) %>%
  rename(train = train.5, pred = pred.5) %>%
  mutate(run = 5) %>%
  na.omit
bg_model <- rbind(bg_model_1, bg_model_2, bg_model_3, bg_model_4, bg_model_5)

bgb_pred_theme <- function(){ 
  font <- "Georgia"   #assign font family up front
  theme_hc() %+replace%    #replace elements we want to change
    theme(
      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      axis.ticks = element_blank(),          #strip axis ticks
      #text elements
      plot.title = element_text(family = font, size = 10, face = 'bold', hjust = 0.5, vjust = 2),
      axis.title = element_text(family = font, size = 10),               
      axis.text = element_text(family = font, size = 7), 
      axis.title.x = element_blank(), 
      panel.background = element_rect(fill = "white", color = "white")
    )
}
library(ggsci)
flatt = pal_flatui(palette = "flattastic")(9)

bg_model_preds <- ggplot(bg_model, aes(x = pred, y = bgm2.core.stemscaled, shape = factor(train), color = site)) +
  geom_point() +
  labs(x = expression("Predicted Belowground Biomass (g " ~ m^-2 ~ ")"), y = expression("Observed Belowground Biomass (g " ~ m^-2 ~ ")")) +
  scale_shape_manual(name = "Split", values = c(1, 16), labels = c("Test", "Train")) +
  scale_color_manual(name = "Site", values = flatt) +
  theme(
    legend.position = "right",
    legend.spacing.y = unit(0, 'in'),
    legend.key = element_rect(fill = "white"),
    legend.key.size = unit(0.15, 'in'),
    text = element_text(size = 8, family = "Georgia"),
    panel.background = element_rect(fill = "white", color = "white"))
bg_model_preds
ggsave(file = "results/bg_model_preds.jpg", bg_model_preds, width = 3, height = 3, units = 'in')

bg_model_test_avg <- bg_model %>%
  group_by(group, pix, site) %>%
  filter(train == "0") %>%
  summarise(resp = bgm2.core.stemscaled, mean_pred = mean(pred))

### fad test
bg_model_test_avg_fad <- bg_model_test_avg %>%
  filter(pix == "fad") %>%
  mutate(diff = resp - mean_pred)
mae_fad <- (sum((bg_model_test_avg_fad$mean_pred-bg_model_test_avg_fad$resp), na.rm=T))/length(bg_model_test_avg_fad$resp[!is.na(bg_model_test_avg_fad$resp)])
#read veg from smoothed. stem density
fad_veg <- veg[veg$plot == "s1" | veg$plot == "s2" | veg$plot == "s6",]
fad_veg1 <- fad_veg[fad_veg$date > as.Date("2016-05-01") & fad_veg$date < as.Date("2016-06-01"),]
fad_veg2 <- fad_veg[fad_veg$date > as.Date("2017-08-01") & fad_veg$date < as.Date("2017-12-01"),]
fad_veg3 <- fad_veg[fad_veg$date > as.Date("2022-02-01") & fad_veg$date < as.Date("2022-04-01"),]
fad_veg_d <- rbind(fad_veg1, fad_veg2, fad_veg3)
fad_core <- mean(fad_veg_d$core_live_stem_count/fad_veg_d$core_area_cm2*1000, na.rm = T)
fad_plot <- mean(fad_veg_d$stemsm2, na.rm = T)
fad_ratio_mean <- fad_core/fad_plot
fad_ratio <- mean(((fad_veg_d$core_live_stem_count/fad_veg_d$core_area_cm2*1000)/fad_veg_d$stemsm2), na.rm = T)
entire_core <- mean(veg$core_live_stem_count/veg$core_area_cm2*1000, na.rm = T)
entire_plot <- mean(veg$stemsm2, na.rm = T)
entire_ratio_mean <- entire_core/entire_plot
entire_ratio <- mean((veg$core_live_stem_count/veg$core_area_cm2*1000/veg$stemsm2), na.rm = T)
entire_ratio_sd <- sd((veg$core_live_stem_count/veg$core_area_cm2*1000/veg$stemsm2), na.rm = T)

entire_ratio/fad_ratio


bg_model_preds_train_avg <- ggplot(bg_model_test_avg, aes(x = mean_pred, y = resp), color = "grey20") +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 1981, ymax = Inf), fill = "grey90") +
  geom_abline(intercept = 0, slope = 1, lty = 1, col = "grey") +
  geom_point(size = 0.25) +
  coord_equal(ratio = 1, xlim = c(0,4000), ylim = c(0,4000)) +
  labs(x = expression("Predicted Belowground Biomass (g " ~ m^-2 ~ ")"), y = expression("Observed Belowground Biomass (g " ~ m^-2 ~ ")")) +
  scale_color_manual(name = "Site", values = flatt) +
  theme(
    legend.position = "none",
    legend.spacing.y = unit(0, 'in'),
    legend.key = element_rect(fill = "white"),
    legend.key.size = unit(0.15, 'in'),
    legend.text = element_text(size = 10, family = "Georgia"),
    text = element_text(size = 8, family = "Georgia"),
    panel.background = element_rect(fill = "white", color = "white"),
    aspect.ratio = 1)
bg_model_preds_train_avg
ggsave(file = "results/bg_model_preds_train_avg.jpg", bg_model_preds_train_avg, width = 4, height = 3, units = 'in')

#(figure S6)
residuals <- ggplot(data = bg_model_test_avg, aes(x = resp, y = resp-mean_pred)) +
  geom_rect(aes(xmin = 1981, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "grey90") +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 0, lty = 1) +
  geom_segment(aes(x = 1981, y = -0, xend = 1981, yend = 343), col = "red", lty = 2) +
  geom_segment(aes(x = -Inf, y = 343, xend = 1981, yend = 343), col = "red", lty = 2) +
  labs(x = expression("Observed Belowground Biomass (g " ~ m^-2 ~ ")"), y = expression("Residual (g " ~ m^-2 ~ ")")) +
  theme(
    legend.position = "none",
    legend.spacing.y = unit(0, 'in'),
    legend.key = element_rect(fill = "white"),
    legend.key.size = unit(0.15, 'in'),
    legend.text = element_text(size = 10, family = "Georgia"),
    text = element_text(size = 8, family = "Georgia"),
    panel.background = element_rect(fill = "white", color = "white"),
    aspect.ratio = 1)
residuals
ggsave(file = "results/bg_model_residuals.jpg", residuals, width = 4, height = 3, units = 'in')


residual_model <- lm((bg_model_test_avg$resp - bg_model_test_avg$mean_pred) ~ bg_model_test_avg$resp)
summary(residual_model)

table <- read_csv("Rmarkdown/metrics_table.csv")
rmse_comp <- ggplot(data = table[table$Variable == "Belowground Biomass",], aes(x = Model, y = RMSE, fill = Model)) +
  geom_col() +
  scale_fill_manual(values = c("BERMv1.0" = "#A997AB", "BERMv2.0" = "#7EC488")) +
  labs(y = expression("RMSE (g " ~ m^-2 ~ ")")) +
  theme(
    legend.position = "bottom",
    legend.spacing.y = unit(0, 'in'),
    legend.key = element_rect(fill = "white"),
    legend.key.size = unit(0.15, 'in'),
    legend.text = element_text(size = 10, family = "Georgia"),
    text = element_text(size = 8, family = "Georgia"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    panel.background = element_rect(fill = "white", color = "white"),
    aspect.ratio = 2/1)
rmse_comp

nrmse_comp <- ggplot(data = table[table$Variable == "Belowground Biomass",], aes(x = Model, y = nRMSE, fill = Model)) +
  geom_col() +
  scale_fill_manual(values = c("BERMv1.0" = "#A997AB", "BERMv2.0" = "#7EC488")) +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "nRMSE") +
  theme(
    legend.position = "none",
    legend.key = element_rect(fill = "white"),
    legend.key.size = unit(0.15, 'in'),
    text = element_text(size = 8, family = "Georgia"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    panel.background = element_rect(fill = "white", color = "white"),
    aspect.ratio = 2/1)
nrmse_comp

library(patchwork)
berm_perf <- bg_model_preds_train_avg + ((rmse_comp + nrmse_comp)/guide_area() + plot_layout(guides = "collect", heights = c(1000,1))) + 
  plot_layout(widths = c(2,2)) +
  plot_annotation(tag_levels = "a")
ggsave(file = "results/berm_perf.jpg", berm_perf, width = 6, height = 3.2, dpi = 600)

#figure 7
#berm2.0 paper 2022 flux marsh
doboy_rivers <- vect("/home/kyle/Documents/Documents/UT/Sapelo/Mapping/GA_streams/Rivers_Streams_Doboy.shp")
doboy_rivers <- project(doboy_rivers, "epsg:32617")
ga_rivers <- vect("/home/kyle/Documents/Documents/UT/Sapelo/Mapping/GA_streams/Rivers_Streams_GA.shp")
ga_rivers <- project(ga_rivers, "epsg:32617")
ga_rivers2 <- terra::crop(ga_rivers, ext(464000, 484000, 3470000, 3490000))
flux_marsh <- vect("/home/kyle/Documents/Documents/UT/Sapelo/Mapping/Flux Tower/FluxOutline.shp")
flux_marsh <- project(flux_marsh, "epsg:32617")
ga_coast <- vect("/home/kyle/Documents/Documents/UT/Sapelo/Mapping/Coastline/gshhg_coastline_clipped2.shp")
ga_coast <- project(ga_coast, "epsg:32617")
cut_no <- length(list.files(paste0("/media/kyle/Seagate Expansion Drive/data/results/model_xgb_belowground_biomass/", "0306020408"), pattern = ".csv", full.names = T))
out<- data.frame()

for (c in 1:(cut_no)){
  bgb <- fread(paste0("/media/kyle/Seagate Expansion Drive/data/results/model_xgb_belowground_biomass/", "0306020408", "/modeled_xgb_belowground_biomass_", "0306020408", "_cutno", c, ".csv"), select = c("pix", "date","predbg.xgb", "predag.allom.l8", "predn.l8", "predchl.l8", "predlai.l8", "elevation", "flood_time", "local_hiwater", "local_lowater"))
  bgb <- bgb[bgb$date == as.Date("2022-11-15"),]
  bgb <- bgb %>%
    drop_na(pix)
  bgb$xcoord <- substr(bgb$pix, start = 1, stop = 6)
  bgb$ycoord <- substr(bgb$pix, start = 8, stop = 15)
  bgb2 <- as_spatraster(bgb, xycols = 12:13, crs = 32617)
  bgb3 <- terra::mask(bgb2, flux_marsh)
  bgb4 <- as.data.frame(bgb3)
  out <- rbind(out, bgb4)
}
write_csv(out, "results/berm2.0_flux_maps/data_202211.csv")

df <- read_csv("results/berm2.0_flux_maps/data_202202.csv")
df$xcoord <- substr(df$pix, start = 1, stop = 6)
df$ycoord <- substr(df$pix, start = 8, stop = 15)
df$date <- as.Date(df$date)
bgb02 <- df[,c(1,2,3)]
bgb02$xcoord <- substr(bgb02$pix, start = 1, stop = 6)
bgb02$ycoord <- substr(bgb02$pix, start = 8, stop = 15)
bgb02 <- bgb02[,c(3,4,5)]
bgb02 <- as_spatraster(bgb02, xycols = 2:3, crs = 32617)

df <- read_csv("results/berm2.0_flux_maps/data_202205.csv")
df$xcoord <- substr(df$pix, start = 1, stop = 6)
df$ycoord <- substr(df$pix, start = 8, stop = 15)
df$date <- as.Date(df$date)
bgb05 <- df[,c(1,2,3)]
bgb05$xcoord <- substr(bgb05$pix, start = 1, stop = 6)
bgb05$ycoord <- substr(bgb05$pix, start = 8, stop = 15)
bgb05 <- bgb05[,c(3,4,5)]
bgb05 <- as_spatraster(bgb05, xycols = 2:3, crs = 32617)

df <- read_csv("results/berm2.0_flux_maps/data_202208.csv")
df$xcoord <- substr(df$pix, start = 1, stop = 6)
df$ycoord <- substr(df$pix, start = 8, stop = 15)
df$date <- as.Date(df$date)
bgb08 <- df[,c(1,2,3)]
bgb08$xcoord <- substr(bgb08$pix, start = 1, stop = 6)
bgb08$ycoord <- substr(bgb08$pix, start = 8, stop = 15)
bgb08 <- bgb08[,c(3,4,5)]
bgb08 <- as_spatraster(bgb08, xycols = 2:3, crs = 32617)

df <- read_csv("results/berm2.0_flux_maps/data_202211.csv")
df$xcoord <- substr(df$pix, start = 1, stop = 6)
df$ycoord <- substr(df$pix, start = 8, stop = 15)
df$date <- as.Date(df$date)
bgb11 <- df[,c(1,2,3)]
bgb11$xcoord <- substr(bgb11$pix, start = 1, stop = 6)
bgb11$ycoord <- substr(bgb11$pix, start = 8, stop = 15)
bgb11 <- bgb11[,c(3,4,5)]
bgb11 <- as_spatraster(bgb11, xycols = 2:3, crs = 32617)

df <- read_csv("results/berm2.0_flux_maps/data_202205.csv")
df$xcoord <- substr(df$pix, start = 1, stop = 6)
df$ycoord <- substr(df$pix, start = 8, stop = 15)
df$date <- as.Date(df$date)

agb <- df[,c(1,2,4)]
agb$xcoord <- substr(agb$pix, start = 1, stop = 6)
agb$ycoord <- substr(agb$pix, start = 8, stop = 15)
agb <- agb[,c(3,4,5)]
agb2 <- as_spatraster(agb, xycols = 2:3, crs = 32617)

nit <- df[,c(1,2,5)]
nit$xcoord <- substr(nit$pix, start = 1, stop = 6)
nit$ycoord <- substr(nit$pix, start = 8, stop = 15)
nit <- nit[,c(3,4,5)]
nit2 <- as_spatraster(nit, xycols = 2:3, crs = 32617)

chl <- df[,c(1,2,6)]
chl$xcoord <- substr(chl$pix, start = 1, stop = 6)
chl$ycoord <- substr(chl$pix, start = 8, stop = 15)
chl <- chl[,c(3,4,5)]
chl2 <- as_spatraster(chl, xycols = 2:3, crs = 32617)

lai <- df[,c(1,2,7)]
lai$xcoord <- substr(lai$pix, start = 1, stop = 6)
lai$ycoord <- substr(lai$pix, start = 8, stop = 15)
lai <- lai[,c(3,4,5)]
lai2 <- as_spatraster(lai, xycols = 2:3, crs = 32617)

elv <- df[,c(1,2,8)]
elv$xcoord <- substr(elv$pix, start = 1, stop = 6)
elv$ycoord <- substr(elv$pix, start = 8, stop = 15)
elv <- elv[,c(3,4,5)]
elv2 <- as_spatraster(elv, xycols = 2:3, crs = 32617)

flo <- df[,c(1,2,9)]
flo$xcoord <- substr(flo$pix, start = 1, stop = 6)
flo$ycoord <- substr(flo$pix, start = 8, stop = 15)
flo <- flo[,c(3,4,5)]
flo2 <- as_spatraster(flo, xycols = 2:3, crs = 32617)

hig <- df[,c(1,2,10)]
hig$xcoord <- substr(hig$pix, start = 1, stop = 6)
hig$ycoord <- substr(hig$pix, start = 8, stop = 15)
hig <- hig[,c(3,4,5)]
hig2 <- as_spatraster(hig, xycols = 2:3, crs = 32617)

low <- df[,c(1,2,11)]
low$xcoord <- substr(low$pix, start = 1, stop = 6)
low$ycoord <- substr(low$pix, start = 8, stop = 15)
low <- low[,c(3,4,5)]
low2 <- as_spatraster(low, xycols = 2:3, crs = 32617)

library(ggthemes)
mapmap_theme <- function(){ 
  font <- "Georgia"   #assign font family up front
  theme_hc() %+replace%    #replace elements we want to change
    theme(
      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      #text elements
      plot.title = element_blank(), 
      plot.subtitle = element_blank(), 
      axis.text.x = element_text(family = font, size = 10, vjust =  0),
      axis.text.y = element_text(family = font, size = 10),
      axis.title = element_blank(), 
      axis.ticks.length=unit(-0.25, "cm"), 
      panel.background = element_rect(fill = NA, color = "black"),
      text = element_text(family = font),
      plot.margin = margin(0,0,0,0),
      panel.border = element_rect(color = "black", fill = NA, size = 0)
    )
} 
map_theme <- function(){ 
  font <- "Georgia"   #assign font family up front
  theme_hc() %+replace%    #replace elements we want to change
    theme(
      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      #text elements
      plot.title = element_text(family = font, size = 16, hjust = 0.5, vjust = -1),
      plot.subtitle = element_text(family = font, size = 12, hjust = 0.5, vjust = -1),
      plot.caption = element_text(family = font, size = 16, hjust = 0.5, vjust = 1),
      axis.text = element_blank(), 
      axis.title = element_blank(), 
      axis.ticks = element_blank(), 
      panel.background = element_rect(fill = "white", color = "white"),
      legend.position = "bottom",
      legend.title = element_text(size = 12, family = font),
      legend.text = element_text(size = 8, family = font, angle = 0),
      text = element_text(family = font),
      plot.margin = margin(0,0,0,0)    )
} 



bgb_leg <- expression(BGB~g~m^{-~2})
bgb_02_map <- ggplot(data = bgb02) +
  #geom_sf(data = doboy_rivers, fill = "lightblue", lwd = 0) +
  #geom_sf(data = doboy_sound, col = "red") +
  geom_raster(aes(x = x, y = y, fill = predbg.xgb)) + 
  coord_sf() + 
  scale_x_continuous(breaks = c(473000, 474000), limits = c(472605, 474075)) +
  scale_y_continuous(breaks = c(3478000, 3480000), limits = c(3478065, 3480525)) +
  scale_fill_paletteer_c("grDevices::terrain.colors", na.value = "transparent", direction = -1, limits = c(0,2000)) +
  labs(title = "", fill = bgb_leg, x = "", subtitle = "2022-02-15", y = "") +
  map_theme()
#bgb_02_map

bgb_05_map <- ggplot(data = bgb05) +
  #geom_sf(data = doboy_rivers, fill = "lightblue", lwd = 0) +
  #geom_sf(data = doboy_sound, col = "red") +
  geom_raster(aes(x = x, y = y, fill = predbg.xgb)) + 
  coord_sf() + 
  scale_x_continuous(breaks = c(473000, 474000), limits = c(472605, 474075)) +
  scale_y_continuous(breaks = c(3478000, 3480000), limits = c(3478065, 3480525)) +
  scale_fill_paletteer_c("grDevices::terrain.colors", na.value = "transparent", direction = -1, limits = c(0,2000)) +
  labs(title = "", fill = bgb_leg, x = "", subtitle = "2022-05-15", y = "") +
  map_theme()
#bgb_05_map

bgb_08_map <- ggplot(data = bgb08) +
  #geom_sf(data = doboy_rivers, fill = "lightblue", lwd = 0) +
  #geom_sf(data = doboy_sound, col = "red") +
  geom_raster(aes(x = x, y = y, fill = predbg.xgb)) + 
  coord_sf() + 
  scale_x_continuous(breaks = c(473000, 474000), limits = c(472605, 474075)) +
  scale_y_continuous(breaks = c(3478000, 3480000), limits = c(3478065, 3480525)) +
  scale_fill_paletteer_c("grDevices::terrain.colors", na.value = "transparent", direction = -1, limits = c(0,2000)) +
  labs(title = "", fill = bgb_leg, x = "", subtitle = "2022-08-15", y = "") +
  map_theme()
#bgb_08_map

bgb_11_map <- ggplot(data = bgb11) +
  #geom_sf(data = doboy_rivers, fill = "lightblue", lwd = 0) +
  #geom_sf(data = doboy_sound, col = "red") +
  geom_raster(aes(x = x, y = y, fill = predbg.xgb)) + 
  coord_sf() + 
  scale_x_continuous(breaks = c(473000, 474000), limits = c(472605, 474075)) +
  scale_y_continuous(breaks = c(3478000, 3480000), limits = c(3478065, 3480525)) +
  scale_fill_paletteer_c("grDevices::terrain.colors", na.value = "transparent", direction = -1, limits = c(0,2000)) +
  labs(title = "", fill = bgb_leg, x = "", subtitle = "2022-11-15", y = "") +
  map_theme()
#bgb_11_map


library(plyr)

agb_leg <- expression(g~m^{-~2})
leg_vals <- c(round_any(min(agb$predag.allom.l8),10, f = floor), round_any(max(agb$predag.allom.l8),10,f = ceiling))
agb_map <- ggplot(data = agb2) +
  #geom_sf(data = doboy_rivers, fill = "lightblue", lwd = 0) +
  #geom_sf(data = doboy_sound, col = "red") +
  geom_raster(aes(x = x, y = y, fill = predag.allom.l8)) + 
  coord_sf() + 
  scale_x_continuous(breaks = c(473000, 474000), limits = c(472605, 474075)) +
  scale_y_continuous(breaks = c(3478000, 3480000), limits = c(3478065, 3480525)) +
  scale_fill_paletteer_c("grDevices::terrain.colors", na.value = "transparent", direction = -1, breaks = leg_vals, limits = leg_vals, labels = leg_vals) +
  labs(caption = "AGB", fill = agb_leg, x = "", subtitle = "2022-05-15", y = "") +
  map_theme() +
  theme(legend.position = "bottom")
#agb_map

leg_vals <- c(round_any(min(lai$predlai.l8),0.1, f = floor), round_any(max(lai$predlai.l8),0.1,f = ceiling))
lai_map <- ggplot(data = lai2) +
  #geom_sf(data = doboy_rivers, fill = "lightblue", lwd = 0) +
  #geom_sf(data = doboy_sound, col = "red") +
  geom_raster(aes(x = x, y = y, fill = predlai.l8)) + 
  coord_sf() + 
  scale_x_continuous(breaks = c(473000, 474000), limits = c(472605, 474075)) +
  scale_y_continuous(breaks = c(3478000, 3480000), limits = c(3478065, 3480525)) +
  scale_fill_paletteer_c("grDevices::terrain.colors", na.value = "transparent", direction = -1, breaks = leg_vals, limits = leg_vals, labels = leg_vals) +
  labs(caption = "LAI", fill = "", x = "", subtitle = "2022-05-15", y = "") +
  map_theme()
#lai_map

chl_leg <- expression(g~g^{-~1})
leg_vals <- c((round_any(min(chl$predchl.l8),0.1, f = floor)+0.1), round_any(max(chl$predchl.l8),0.1,f = ceiling))
chl_map <- ggplot(data = chl2) +
  #geom_sf(data = doboy_rivers, fill = "lightblue", lwd = 0) +
  #geom_sf(data = doboy_sound, col = "red") +
  geom_raster(aes(x = x, y = y, fill = predchl.l8)) + 
  coord_sf() + 
  scale_x_continuous(breaks = c(473000, 474000), limits = c(472605, 474075)) +
  scale_y_continuous(breaks = c(3478000, 3480000), limits = c(3478065, 3480525)) +
  scale_fill_paletteer_c("grDevices::YlGn", na.value = "transparent", direction = -1, breaks = leg_vals, limits = leg_vals, labels = leg_vals) +
  labs(caption = "CHL", fill = chl_leg, x = "", subtitle = "2022-05-15", y = "") +
  map_theme()
#chl_map

leg_vals <- c((round_any(min(nit$predn.l8),0.1, f = floor)+0.1), (round_any(max(nit$predn.l8),0.1,f = ceiling))-0.2)
nit_map <- ggplot(data = nit2) +
  #geom_sf(data = doboy_rivers, fill = "lightblue", lwd = 0) +
  #geom_sf(data = doboy_sound, col = "red") +
  geom_raster(aes(x = x, y = y, fill = predn.l8)) + 
  coord_sf() + 
  scale_x_continuous(breaks = c(473000, 474000), limits = c(472605, 474075)) +
  scale_y_continuous(breaks = c(3478000, 3480000), limits = c(3478065, 3480525)) +
  scale_fill_paletteer_c("grDevices::Sunset", na.value = "transparent", direction = -1, breaks = leg_vals, limits = leg_vals, labels = leg_vals) +
  labs(caption = "Foliar N", fill = "%", x = "", subtitle = "2022-05-15", y = "") +
  map_theme()
#nit_map

leg_vals <- c(round_any(min(flo$flood_time),0.1, f = floor), round_any(max(flo$flood_time),0.1,f = ceiling))
flo_map <- ggplot(data = flo2) +
  #geom_sf(data = doboy_rivers, fill = "lightblue", lwd = 0) +
  #geom_sf(data = doboy_sound, col = "red") +
  geom_raster(aes(x = x, y = y, fill = flood_time)) + 
  coord_sf() + 
  scale_x_continuous(breaks = c(473000, 474000), limits = c(472605, 474075)) +
  scale_y_continuous(breaks = c(3478000, 3480000), limits = c(3478065, 3480525)) +
  scale_fill_paletteer_c("viridis::cividis", na.value = "transparent", direction = -1, breaks = leg_vals, limits = leg_vals, labels = leg_vals) +
  labs(caption = "FF", fill = "%", x = "", subtitle = "2022-05-15", y = "") +
  map_theme()
#flo_map

leg_vals <- c(round_any(min(elv$elevation),0.1, f = floor), round_any(max(elv$elevation),0.1,f = ceiling))
elv_map <- ggplot(data = elv2) +
  #geom_sf(data = doboy_rivers, fill = "lightblue", lwd = 0) +
  #geom_sf(data = doboy_sound, col = "red") +
  geom_raster(aes(x = x, y = y, fill = elevation)) + 
  coord_sf() + 
  scale_x_continuous(breaks = c(473000, 474000), limits = c(472605, 474075)) +
  scale_y_continuous(breaks = c(3478000, 3480000), limits = c(3478065, 3480525)) +
  scale_fill_paletteer_c("viridis::cividis", na.value = "transparent", direction = 1, breaks = leg_vals, limits = leg_vals, labels = leg_vals) +
  labs(caption = "Elevation", fill = "m", x = "", subtitle = "2022-05-15", y = "") +
  map_theme()
#elv_map

leg_vals <- c(round_any(min(hig$local_hiwater),0.1, f = floor), (round_any(max(hig$local_hiwater),0.1,f = ceiling))-0.2)
hig_map <- ggplot(data = hig2) +
  #geom_sf(data = doboy_rivers, fill = "lightblue", lwd = 0) +
  #geom_sf(data = doboy_sound, col = "red") +
  geom_raster(aes(x = x, y = y, fill = local_hiwater)) + 
  coord_sf() + 
  scale_x_continuous(breaks = c(473000, 474000), limits = c(472605, 474075)) +
  scale_y_continuous(breaks = c(3478000, 3480000), limits = c(3478065, 3480525)) +
  scale_fill_paletteer_c("grDevices::Geyser", na.value = "transparent", direction = -1, breaks = leg_vals, limits = leg_vals, labels = leg_vals) +
  labs(caption = "II", fill = "", x = "", subtitle = "2022-05-15", y = "") +
  map_theme()
#hig_map

leg_vals <- c((round_any(min(low$local_lowater),0.1, f = floor)+0.3), round_any(max(low$local_lowater),0.1,f = ceiling))
low_map <- ggplot(data = low2) +
  #geom_sf(data = doboy_rivers, fill = "lightblue", lwd = 0) +
  #geom_sf(data = doboy_sound, col = "red") +
  geom_raster(aes(x = x, y = y, fill = local_lowater)) + 
  coord_sf() + 
  scale_x_continuous(breaks = c(473000, 474000), limits = c(472605, 474075)) +
  scale_y_continuous(breaks = c(3478000, 3480000), limits = c(3478065, 3480525)) +
  scale_fill_paletteer_c("grDevices::Geyser", na.value = "transparent", direction = 1, breaks = leg_vals, limits = leg_vals, labels = leg_vals) +
  labs(caption = "DI", fill = "", x = "", subtitle = "2022-05-15", y = "") +
  map_theme()
#low_map

library(ggspatial)
ga_coast2 <- st_as_sf(ga_coast)
class(ga_coast2)
ext <- ext(470000, 483200, 3465000, 3495000)
ga_rivers3 <- crop(ga_rivers2,ext)

flux_point <- data.frame(x = 473362.2, y = 3479003.0)
flux_map <- ggplot(data = ga_coast2) +
  geom_rect(aes(xmin = 470000, xmax = 483200, ymin = 3465000, ymax = 3495000),fill = "lightblue") +
  geom_sf(fill = "white", color = "black", lwd = 0) +
  geom_sf(data = ga_rivers3, fill = "lightblue", lwd = 0) +
  annotation_scale(location = "br", pad_x = unit(0.2, "in"), pad_y = unit(0.58, "in"), width_hint = 0.4, height = unit(0.2, "in"), text_cex = 1.7) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(1.8, "in"),  pad_y = unit(1.0, "in"), height = unit(1, "in"), width = unit(1, "in"),
                         style = north_arrow_fancy_orienteering) +
  #geom_point(data = flux_point, aes(x = x, y = y), color = "gold3", size = 20, shape="\u2605") +
  geom_sf(data = flux_marsh, fill = "gold3", lwd = 0) +
  geom_rect(aes(xmin = 470000, xmax = 483200, ymin = 3465000, ymax = 3495000),fill = NA, color = "black") +
  coord_sf(expand = F, datum=st_crs(32617), label_graticule = "SE") + 
  scale_x_continuous(limits = c(470000, 483200), breaks = c(475000, 480000)) +
  scale_y_continuous(limits = c(3465000, 3495000), position = "right", breaks = c(3470000, 3475000, 3480000, 3485000, 3490000)) +
  labs(x = "", y = "") +
  mapmap_theme() +
  theme(plot.margin = margin(t=0))
#flux_map



library(patchwork)
bgb_maps <- bgb_02_map + bgb_05_map + bgb_08_map + bgb_11_map + plot_layout(ncol = 4, nrow = 1, guides = "collect") + 
  plot_annotation(caption = "BGB") & 
  theme(legend.position = 'bottom', legend.key.width = unit(100, 'pt'), plot.caption = element_text(family = "Georgia", size = 16, hjust = 0.5, vjust = 0)) 
#bgb_maps
ag_maps <- agb_map + lai_map + chl_map + nit_map + plot_layout(ncol = 4, nrow = 1) & theme(legend.position = 'bottom', legend.key.width = unit(8, 'pt'))
#ag_maps
other_maps <- elv_map + flo_map + hig_map + low_map + plot_layout(ncol = 4, nrow = 1) & theme(legend.position = 'bottom', legend.key.width = unit(8, 'pt'))
#other_maps

berm_demonstration_maps <- bgb_maps  /
  ((ag_maps /
      other_maps) |
     flux_map)  + 
  plot_layout(heights = c(1,1.8))
ggsave("results/berm_demonstration_maps_2022.jpg", berm_demonstration_maps, width = 10, height = 14, dpi = 600)

#figure 8
bg_pred_error_type <- read_csv("results/bg_table_by_site_year_v2.0.csv")
letters_only <- function(x) !grepl("[^A-Za-z]", x)
bg_pred_error_type$type <- ifelse(letters_only(bg_pred_error_type$value) == T, "Site", "Year")
bg_pred_error_type <- bg_pred_error_type[bg_pred_error_type$n_test >=10,]
write_csv(bg_pred_error_type, "Rmarkdown/bg_pred_error_type.csv")
#label_y <- expression(atop(x = "", y = atop(x = paste("Predicted"), y = paste("Belowground Biomass (g " ~ m^-2 ~ ")"))))
site_col = "#207890"
year_col = "#E85040"
site.year_col = "#F8C060"
t.test(bg_pred_error_type$nrmse[bg_pred_error_type$type == "Site"], bg_pred_error_type$nrmse[bg_pred_error_type$type == "Year"])
bg_pred_type_sd <- bg_pred_error_type %>%
  group_by(type) %>%
  summarise(sd_nrmse = sd(nrmse))
label_y2 <- expression(atop(x = "", y = atop(x = paste("Predicted Belowground"), y = paste("Biomass Error (nRMSE)"))))

bg_pred_error_site_year <- ggplot(bg_pred_error_type, aes(x = type, y = nrmse, fill = type)) +
  geom_boxplot(position = position_dodge(1), lwd = 0.25, outlier.size = 0.05) +
  stat_compare_means( aes(label = ..p.signif..), method = "wilcox", label.x = 1.4, label.y = 0.37) +
  scale_x_discrete(labels = c("Site", "Year")) +
  scale_y_continuous(labels = scales::percent, limits = c(0,0.4)) +
  scale_fill_manual(values = c(site_col, year_col)) +
  labs(title = "", y = label_y2) +
  bgb_obs_var_theme2() 
bg_pred_sd_site_year_colplot <- ggplot(bg_pred_type_sd, aes(x = type, y = sd_nrmse, fill = type)) +
  geom_col(position = position_dodge(1), lwd = 0.25) +
  scale_x_discrete(labels = c("Site", "Year")) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c(site_col, year_col)) +
  labs(title = "", y = "") +
  bgb_obs_var_theme() 
bg_pred_cc_col <- ggplot(data = CC_data) +
  geom_col(aes(x = type, y = coef/100, fill = type))  +
  labs(y = "Percent of Total\nExplained Variance", x = "", title = "") +
  scale_y_continuous(labels = scales::percent, breaks = c(0,0.25,0.5,0.75, 1.0), limits = c(0,1)) +
  scale_fill_manual(values = c(site_col,  year_col, site.year_col), labels = c("Site", "Year", "Site:Year")) +
  scale_x_discrete(labels = c("Site", "Year", "Site:Year"))  +
  bgb_obs_var_theme() 

### putting it all together
var_obs_pred <-  bg_mean_site_year  + bg_cc_col +
  bg_pred_error_site_year  + bg_pred_cc_col +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a") & theme(plot.margin = unit(c(10,5,0,5), "pt"), plot.tag.position = c(0.1,1))
var_obs_pred
ggsave("results/var_obs_pred.jpg", var_obs_pred, dpi = 300, height = 6, width = 6, units = "in")


