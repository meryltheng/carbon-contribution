library(terra)
library(tidyterra)
library(geodata)
library(ggplot2)
library(ggpubr)

# define directories
input_dir = 'Input/Climate/CMIP6'
output_dir = 'Output/Climate/CMIP6'
plot_dir = "Plots/Climate"

# define models and ssp_scenarios to loop over
models = c("ACCESS-CM2", "CNRM-CM6-1-HR", "MPI-ESM1-2-LR")
ssp_scenarios = c("historical", "ssp126", "ssp245", "ssp370")

# crop to land area
aus <- gadm(country="AUS", level=0, path=tempdir(), resolution=2)
#aus <- crop(aus, ext(suitability[[1]][[1]]))

"----------- Plot climate input layers ---------------"
# -- tmax ("Max. temperature of wettest quarter")
temp_max <- lapply(models, function(m){
  lapply(ssp_scenarios, function(s){
    r <- rast(paste0(file.path(input_dir, m, s),"/temp_max.tif"))
    names(r) <- s
    r = mask(r, aus)
    r
  })
})
temp_max_r <- list(rast(temp_max[[1]]), rast(temp_max[[2]]), rast(temp_max[[3]]))

tmax_plots <- lapply(temp_max_r, function(x){
  names(x) <- c("Present \n(2014)","SSP1-2.6 \n(2100)", "SSP2-4.5 \n(2100)", "SSP3-7.0 \n(2100)")
  plot <- ggplot() +
    geom_spatraster(data = x) +
    facet_wrap(~ lyr, ncol = 4) +
    scale_fill_whitebox_c(
      palette = "bl_yl_rd",
      labels = scales::label_number(suffix = "º"),
      n.breaks = 9,
      limits = c(0,50),
      guide = guide_legend(reverse = TRUE)
    ) +
    labs(fill = "tmax") +
    geom_spatvector(data = aus, fill = NA) +
    theme_void() 
  return(plot)
})

tmax_plot <- ggarrange(tmax_plots[[1]] + theme(strip.text.x = element_text(size=13)) ,
                            tmax_plots[[2]] + theme(strip.text.x = element_blank()),
                            tmax_plots[[3]] + theme(strip.text.x = element_blank()),
                            ncol = 1, nrow = 3, align = "h", common.legend = TRUE, legend = "right") 

tmax_plot <- annotate_figure(tmax_plot,
                                  left = text_grob("    MPI-ESM1-2-LR                    CNRM-CM6-1-HR             ACCESS-CM2       ", rot = 90, vjust = 1, size = 12, face = "bold"))

ggsave(filename=paste0(plot_dir, "/tmax.png"), plot=tmax_plot, width=8.5, height=6.5) 

# -- tmin ("Min. temperature of wettest quarter")
temp_min <- lapply(models, function(m){
  lapply(ssp_scenarios, function(s){
    r <- rast(paste0(file.path(input_dir, m, s),"/temp_min.tif"))
    names(r) <- s
    r = mask(r, aus)
    r
  })
})
temp_min_r <- list(rast(temp_min[[1]]), rast(temp_min[[2]]), rast(temp_min[[3]]))

tmin_plots <- lapply(temp_min_r, function(x){
  names(x) <- c("Present \n(2014)","SSP1-2.6 \n(2100)", "SSP2-4.5 \n(2100)", "SSP3-7.0 \n(2100)")
  plot <- ggplot() +
    geom_spatraster(data = x) +
    facet_wrap(~ lyr, ncol = 4) +
    scale_fill_whitebox_c(
      palette = "bl_yl_rd",
      labels = scales::label_number(suffix = "º"),
      n.breaks = 9,
      limits = c(0,50),
      guide = guide_legend(reverse = TRUE)
    ) +
    labs(fill = "tmin") +
    geom_spatvector(data = aus, fill = NA) +
    theme_void() 
  return(plot)
})

tmin_plot <- ggarrange(tmin_plots[[1]] + theme(strip.text.x = element_text(size=13)) ,
                       tmin_plots[[2]] + theme(strip.text.x = element_blank()),
                       tmin_plots[[3]] + theme(strip.text.x = element_blank()),
                       ncol = 1, nrow = 3, align = "h", common.legend = TRUE, legend = "right") 

tmin_plot <- annotate_figure(tmin_plot,
                             left = text_grob("    MPI-ESM1-2-LR                    CNRM-CM6-1-HR             ACCESS-CM2       ", rot = 90, vjust = 1, size = 12, face = "bold"))

ggsave(filename=paste0(plot_dir, "/tmin.png"), plot=tmin_plot, width=8.5, height=6.5) 

# -- srad ("Solar radiation")
solar_rad <- lapply(models, function(m){
  lapply(ssp_scenarios, function(s){
    r <- rast(paste0(file.path(input_dir, m, s),"/solar_rad.tif"))
    names(r) <- s
    r = mask(r, aus)
    r
  })
})
solar_rad_r <- list(rast(solar_rad[[1]]), rast(solar_rad[[2]]), rast(solar_rad[[3]]))

srad_plots <- lapply(solar_rad_r, function(x){
  names(x) <- c("Present \n(2014)","SSP1-2.6 \n(2100)", "SSP2-4.5 \n(2100)", "SSP3-7.0 \n(2100)")
  plot <- ggplot() +
    geom_spatraster(data = x) +
    facet_wrap(~ lyr, ncol = 4) +
    scale_fill_grass_c(
      palette = "inferno",
      n.breaks = 9,
      guide = guide_legend(reverse = TRUE)
    ) +
    labs(fill = bquote(atop("Solar radiation",
                            "(W" ~ m ^ 2 ~ h ^ -1 ~")"))) +
    geom_spatvector(data = aus, fill = NA) +
    theme_void() 
  return(plot)
})

srad_plot <- ggarrange(srad_plots[[1]] + theme(strip.text.x = element_text(size=13)) ,
                       srad_plots[[2]] + theme(strip.text.x = element_blank()),
                       srad_plots[[3]] + theme(strip.text.x = element_blank()),
                       ncol = 1, nrow = 3, align = "h", common.legend = TRUE, legend = "right") 

srad_plot <- annotate_figure(srad_plot,
                             left = text_grob("    MPI-ESM1-2-LR                    CNRM-CM6-1-HR             ACCESS-CM2       ", rot = 90, vjust = 1, size = 12, face = "bold"))

ggsave(filename=paste0(plot_dir, "/srad.png"), plot=srad_plot, width=8.5, height=6.5) 

# -- hours_of_high_RH
hours_of_high_RH <- lapply(models, function(m){
  lapply(ssp_scenarios, function(s){
    r <- rast(paste0(file.path(input_dir, m, s),"/hours_of_high_RH.tif"))
    names(r) <- s
    r = mask(r, aus)
    r
  })
})
hours_of_high_RH_r <- list(rast(hours_of_high_RH[[1]]), rast(hours_of_high_RH[[2]]), rast(hours_of_high_RH[[3]]))

RH_plots <- lapply(hours_of_high_RH_r, function(x){
  names(x) <- c("Present \n(2014)","SSP1-2.6 \n(2100)", "SSP2-4.5 \n(2100)", "SSP3-7.0 \n(2100)")
  plot <- ggplot() +
    geom_spatraster(data = x) +
    facet_wrap(~ lyr, ncol = 4) +
    scale_fill_whitebox_c(
      palette = "deep",
      n.breaks = 12,
      limits = c(0,24),
      guide = guide_legend(reverse = TRUE)
    ) +
    labs(fill = "Hours of \nhigh RH") +
    geom_spatvector(data = aus, fill = NA) +
    theme_void() 
  return(plot)
})

RH_plot <- ggarrange(RH_plots[[1]] + theme(strip.text.x = element_text(size=13)) ,
                       RH_plots[[2]] + theme(strip.text.x = element_blank()),
                       RH_plots[[3]] + theme(strip.text.x = element_blank()),
                       ncol = 1, nrow = 3, align = "h", common.legend = TRUE, legend = "right") 

RH_plot <- annotate_figure(RH_plot,
                             left = text_grob("    MPI-ESM1-2-LR                    CNRM-CM6-1-HR             ACCESS-CM2       ", rot = 90, vjust = 1, size = 12, face = "bold"))

ggsave(filename=paste0(plot_dir, "/hours_of_high_RH.png"), plot=RH_plot, width=8.5, height=6.5) 

"----------- Plot infection risk maps ---------------"
# -- inf. risk at tmin
infection_risk_tmin <- lapply(models, function(m){
  lapply(ssp_scenarios, function(s){
    r <- rast(paste0(file.path(output_dir, m, s),"/infection_risk_tmin.tif"))
    names(r) <- s
    r = mask(r, aus)
    r
  })
})
infection_risk_tmin_r <- list(rast(infection_risk_tmin[[1]]), rast(infection_risk_tmin[[2]]), rast(infection_risk_tmin[[3]]))

risk_tmin_plots <- lapply(infection_risk_tmin_r, function(x){
  names(x) <- c("Present \n(2014)","SSP1-2.6 \n(2100)", "SSP2-4.5 \n(2100)", "SSP3-7.0 \n(2100)")
  plot <- ggplot() +
    geom_spatraster(data = x) +
    facet_wrap(~ lyr, ncol = 4) +
    #scale_fill_manual(values = c("grey90", "black"), na.value = "transparent", na.translate = F) +
    scale_fill_whitebox_c(
      palette = "viridi",
      n.breaks = 10,
      limits = c(0,1),
      guide = guide_legend(reverse = TRUE)
    ) +
    labs(fill = "Infection \nrisk at tmin") +
    geom_spatvector(data = aus, fill = NA) +
    theme_void() 
  return(plot)
})

risk_tmin_plot <- ggarrange(risk_tmin_plots[[1]] + theme(strip.text.x = element_text(size=13)) ,
                            risk_tmin_plots[[2]] + theme(strip.text.x = element_blank()),
                            risk_tmin_plots[[3]] + theme(strip.text.x = element_blank()),
                            ncol = 1, nrow = 3, align = "h", common.legend = TRUE, legend = "right") 

risk_tmin_plot <- annotate_figure(risk_tmin_plot,
                                  left = text_grob("    MPI-ESM1-2-LR                    CNRM-CM6-1-HR             ACCESS-CM2       ", rot = 90, vjust = 1, size = 12, face = "bold"))

ggsave(filename=paste0(plot_dir, "/risk_tmin.png"), plot=risk_tmin_plot, width=8.5, height=6.5) 

# -- inf. risk at tmax
infection_risk_tmax <- lapply(models, function(m){
  lapply(ssp_scenarios, function(s){
    r <- rast(paste0(file.path(output_dir, m, s),"/infection_risk_tmax.tif"))
    names(r) <- s
    r = mask(r, aus)
    r
  })
})

infection_risk_tmax_r <- list(rast(infection_risk_tmax[[1]]), rast(infection_risk_tmax[[2]]), rast(infection_risk_tmax[[3]]))

risk_tmax_plots <- lapply(infection_risk_tmax_r, function(x){
  names(x) <- c("Present \n(2014)","SSP1-2.6 \n(2100)", "SSP2-4.5 \n(2100)", "SSP3-7.0 \n(2100)")
  plot <- ggplot() +
    geom_spatraster(data = x) +
    facet_wrap(~ lyr, ncol = 4) +
    #scale_fill_manual(values = c("grey90", "black"), na.value = "transparent", na.translate = F) +
    scale_fill_whitebox_c(
      palette = "viridi",
      n.breaks = 10,
      limits = c(0,1),
      guide = guide_legend(reverse = TRUE)
    ) +
    labs(fill = "Infection \nrisk at tmax") +
    geom_spatvector(data = aus, fill = NA) +
    theme_void() 
  return(plot)
})

risk_tmax_plot <- ggarrange(risk_tmax_plots[[1]] + theme(strip.text.x = element_text(size=13)) ,
                            risk_tmax_plots[[2]] + theme(strip.text.x = element_blank()),
                            risk_tmax_plots[[3]] + theme(strip.text.x = element_blank()),
                            ncol = 1, nrow = 3, align = "h", common.legend = TRUE, legend = "right") 

risk_tmax_plot <- annotate_figure(risk_tmax_plot,
                                  left = text_grob("    MPI-ESM1-2-LR                    CNRM-CM6-1-HR             ACCESS-CM2       ", rot = 90, vjust = 1, size = 12, face = "bold"))

ggsave(filename=paste0(plot_dir, "/risk_tmax.png"), plot=risk_tmax_plot, width=8.5, height=6.5) 

"----------- Plot climate suitability maps ---------------"
# climate suitability
threshold_val = 0.5 # 0.2 or 0.5
suitability <- lapply(models, function(m){
  lapply(ssp_scenarios, function(s){
    r <- rast(paste0(file.path(output_dir, m, s),"/suitability_", threshold_val, ".tif"))
    names(r) <- s
    r = mask(r, aus)
    r
  })
})
suitability_r <- list(rast(suitability[[1]]), rast(suitability[[2]]), rast(suitability[[3]]))

aus <- crop(aus, ext(suitability[[1]][[1]]))

suitability_plots <- lapply(suitability_r, function(x){
  names(x) <- c("Present \n (2014)","SSP1-2.6 \n (2100)", "SSP2-4.5 \n (2100)", "SSP3-7.0 \n (2100)")
  plot <- ggplot() +
      geom_spatraster(data = as.factor(x)) +
      facet_wrap(~ lyr, ncol = 4) +
      scale_fill_manual(values = c("grey90", "black"), na.value = "transparent", na.translate = F) +
      geom_spatvector(data = aus, fill = NA) +
      labs(fill = "Suitability") +
      theme_void()
    return(plot)
  })

suitability_plot <- ggarrange(suitability_plots[[1]] + theme(strip.text.x = element_text(size=13)) ,
                              suitability_plots[[2]] + theme(strip.text.x = element_blank()),
                              suitability_plots[[3]] + theme(strip.text.x = element_blank()),
                             ncol = 1, nrow = 3, align = "h", common.legend = TRUE, legend = "right") #

suitability_plot <- annotate_figure(suitability_plot,
                                left = text_grob("    MPI-ESM1-2-LR                    CNRM-CM6-1-HR             ACCESS-CM2       ", rot = 90, vjust = 1, size = 12, face = "bold"))


#ggsave(filename=paste0(plot_dir, "/suitability_", threshold_val,".pdf"), plot=suitability_plot, width=8.5, height=6.5) 
ggsave(filename=paste0(plot_dir, "/suitability_", threshold_val,".png"), plot=suitability_plot, width=8.5, height=6.5) 

"----------- Apply clim. suit. mask to damage map ---------------"
nvis_layer <- terra::rast('Data/Spatial/NVIS_MVG_RAS_AA.tif')
load(file = "Output/Damage_plus_tables.RData")

# Australia boundary
aus <- gadm(country="AUS", level=0, path=tempdir(), resolution=2)
aus <- project(aus, crs(nvis_layer)) # reproject climate layer to nvis layer
aus <- crop(aus, ext(nvis_layer)) # crop climate layer to nvis layer to match extents

# pick a model and scenario
m = "ACCESS-CM2"
s = "historical"
threshold = 0.5

# get impact % and carbon lost values (same for all scenarios)
Impact_perc <- strsplit(output_tables[[m]][[s]][[paste(threshold)]]$Impact_perc, " ") 
Impact_perc <- as.numeric(sapply(Impact_perc, "[[" , 1))

C_lost_ha <- strsplit(output_tables[[m]][[s]][[paste(threshold)]]$C_lost_ha, " ")
C_lost_ha <- as.numeric(sapply(C_lost_ha, "[[" , 1))
C_lost_ha[is.na(C_lost_ha)] <- 0

# create impact and carbon lost layer based on MVG classification (before mask)
impact_mat <- matrix(c(c(1:32,99), c(Impact_perc,0)), ncol = 2, byrow = F) # is-to-become
impact_layer <- classify(nvis_layer, impact_mat)
closs_mat <- matrix(c(c(1:32,99), c(C_lost_ha,0)), ncol = 2, byrow = F) # is-to-become
closs_layer <- classify(nvis_layer, closs_mat)

# load climate mask layer
output_dir = 'Output/Climate/CMIP6' # climate mask directory
clim_mask <- rast(paste0(file.path(output_dir, m, s),"/suitability_", threshold, ".tif"))
#clim_mask = mask(clim_mask, aus)

for (s in c("historical", "ssp245")){
  
  # read climate mask
  clim_mask <- rast(paste0(file.path(output_dir, m, s),"/suitability_", threshold, ".tif"))
  
  # reproject and resample climate mask layer to nvis layer
  clim_mask <- project(clim_mask, crs(nvis_layer), method = "near") # reproject climate layer to nvis layer
  clim_mask <- crop(clim_mask, nvis_layer) # crop climate layer to nvis layer to match extents
  empty_layer <- rast(nrows=nrow(nvis_layer), ncols=ncol(nvis_layer), # create empty raster with the same extent and resolution as nvis layer
                      xmin=terra::ext(nvis_layer)[1], xmax=terra::ext(nvis_layer)[2],
                      ymin=terra::ext(nvis_layer)[3], ymax=terra::ext(nvis_layer)[4])
  clim_mask_resampled <- terra::resample(clim_mask, empty_layer, method="near") # Resample climate suitability layer to resolution as nvis layer
  
  # mask impact and carbon lost layer with climate mask
  impact_masked = impact_layer
  impact_masked[clim_mask_resampled==0] = NA # set all cells that are not suitable to NA
  plot(impact_masked)
  
  closs_masked = closs_layer
  closs_masked[clim_mask_resampled==0] = NA # set all cells that are not suitable to NA
  plot(closs_masked)
  
  # create pattern mask to show unsuitable areas
  clim_mask_crop <- crop(clim_mask, aus, mask = TRUE, snap = "near")
  clim_mask_crop[clim_mask_crop==0] = NA # delineate suitable areas
  clim_mask_poly = as.polygons(clim_mask_crop)
  pts <- spatSample(aus, size = 350, method = "regular") # grid pattern
  pts <- mask(pts, clim_mask_poly, inverse = TRUE)
  
  map_layers[[s]] <- list(closs_masked,
                          impact_masked,
                          clim_mask_poly,
                          pts)
rm(clim_mask, clim_mask_resampled, clim_mask_crop, impact_masked, closs_masked, 
   pts, clim_mask_poly, empty_layer)
}

save(map_layers, file = "Output/Climate/damage_map_data.RData")

# PLOT
damage_plots <- list()
# load("Output/Climate/damage_map_data.RData")

for (i in 1:2){
  # - carbon
  closs_plot <- ggplot() +
    geom_spatraster(data = map_layers[[i]][[1]]) +
    scale_fill_whitebox_c(
      palette = "muted",
      #labels = scales::label_number(suffix = " t/ha"),
      n.breaks = 5,
      guide = guide_legend(reverse = TRUE)
    ) +
    geom_spatvector(data = aus, fill = NA) +
    geom_spatvector(data = map_layers[[i]][[3]], fill = NA) +
    geom_sf(data = map_layers[[i]][[4]], color = "grey50", shape = 4, size = 1) +
    theme_void() +
    labs(fill = "Carbon \nsequestration\npotential lost \n(t/ha)", title = NULL)
  
  # - damages
  impact_plot <- ggplot() +
    geom_spatraster(data = map_layers[[i]][[2]]) +
    scale_fill_whitebox_c(
      palette = "bl_yl_rd",
      #labels = scales::label_number(suffix = " %"),
      n.breaks = 5,
      guide = guide_legend(reverse = TRUE)
    ) +
    geom_spatvector(data = aus, fill = NA) +
    geom_spatvector(data = map_layers[[i]][[3]], fill = NA) +
    geom_sf(data = map_layers[[i]][[4]], color = "grey50", shape = 4, size = 1) +
    theme_void() +
    labs(fill = "Asset value \nreduction (%)", title = NULL)
  
  damage_plots[[i]] <- list(closs_plot,
                            impact_plot)
}

closs_maps <- ggarrange(damage_plots[[1]][[1]] + ggtitle("Present (2014)") + theme(plot.title = element_text(hjust = 0.5)), 
                        damage_plots[[2]][[1]] + ggtitle("Future (2100)") + theme(plot.title = element_text(hjust = 0.5)), 
                        ncol = 2, nrow = 1, 
                        common.legend = T, legend = "right",
                        labels = c("(a)", "(b)")) 

impact_maps <- ggarrange(damage_plots[[1]][[2]] + ggtitle("Present (2014)") + theme(plot.title = element_text(hjust = 0.5)), 
                         damage_plots[[2]][[2]] + ggtitle("Future (2100)") + theme(plot.title = element_text(hjust = 0.5)), 
                         ncol = 2, nrow = 1, 
                         common.legend = T, legend = "right",
                         labels = c("(c)", "(d)"))

damages_map <- ggarrange(closs_maps, impact_maps, 
                         ncol = 1, nrow = 2, align = "h")

ggsave(file.path(plot_dir, paste0("damage_map.png")), damages_map, width = 8, height = 6.67)


"----------- Overall damages ---------------"
# load sample results if processed
damage_total <- read.csv("Output/Damage_plus_plot.csv")

# conversions
damage_total$Damage <- damage_total$Damage/1000000 # in $m
damage_total$threshold_val <- as.factor(damage_total$threshold_val)

# plot
threshold_cols <- c("grey50", "grey80") #c("tomato", "gold")
damages_plot <- ggviolin(damage_total, x = "ssp_scenario", y = "Damage", fill = "threshold_val", facet.by = 'model',
                         palette = threshold_cols, size = 0.3, add = "median_iqr", width = 0.75, add.params = list(size=0.1)) + 
  ggtitle('') + xlab('Scenario') + ylab('Damages ($M)') + 
  #scale_x_discrete(labels = c("historical" = "Present \n (2014)", "ssp126" = "SSP1-2.6 \n (2100)", "ssp245" = "SSP2-4.5 \n (2100)", "ssp370" = "SSP3-7.0 \n (2100)")) +
  scale_x_discrete(labels = c("historical" = "Present\n(2014)", "ssp126" = "Future\n(2100)\nSSP1-2.6", "ssp245" = "Future\n(2100)\nSSP2-4.5", "ssp370" = "Future\n(2100)\nSSP3-7.0")) +
  scale_fill_manual(name = "Infection risk threshold", labels = c("0.2" = "Low (>0.2)", "0.5" = "Moderate (>0.5)"),
                    values = threshold_cols) +
  ylim(0, 1200) +
  theme(strip.background = element_blank(), 
        strip.text.x = element_text(size=11),
        axis.text.x = element_text(size=10),
        axis.title = element_text(size=12, face = "bold"))

ggsave(paste0(plot_dir, "/overall_damages.png"), damages_plot, width = 10, height = 3.5, dpi = 300)

# ggerrorplot(damage_total, x = "ssp_scenario", y = "Damage", color = "threshold_val", facet.by = 'model',
#          palette = threshold_cols, desc_stat = "median_iqr", 
#          size = 0.3, width = 0.5, add.params = list(size=0.1),
#          position = position_dodge(0.3)) + 
#   ggtitle('') + xlab('SSP scenario') + ylab('Damages ($m)') + 
#   theme(strip.background = element_blank(), 
#         strip.text.x = element_text(size=11))