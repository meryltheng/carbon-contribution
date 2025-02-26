library(terra)
library(geodata)

# define directories
input_dir = 'Input/Climate/CMIP6'
plot_dir = "Plots/Climate"

# define models and ssp_scenarios to loop over
models = c("ACCESS-CM2", "CNRM-CM6-1-HR", "MPI-ESM1-2-LR")
ssp_scenarios = c("historical", "ssp126", "ssp245", "ssp370")

# crop to land area
aus <- gadm(country="AUS", level=0, path=tempdir(), resolution=2)

"----------- Plot climate input layers ---------------"
# -- tmax
temp_max <- lapply(models, function(m){
  lapply(ssp_scenarios, function(s){
    r <- rast(paste0(file.path(input_dir, m, s),"/temp_max.tif"))
    names(r) <- s
    r = mask(r, aus)
    r
  })
})
temp_max_r <- list(rast(temp_max[[1]]), rast(temp_max[[2]]), rast(temp_max[[3]]))
temp_max_r <-  c(temp_max_r[[1]], temp_max_r[[2]], temp_max_r[[3]])
names(temp_max_r) <- apply(expand.grid(ssp_scenarios, models), 1, paste, collapse=".")
aus <- crop(aus, temp_max_r[[1]])

tmax_plot <- ggplot() +
  geom_spatraster(data = temp_max_r) +
  facet_wrap(~ lyr, ncol = 4) +
  scale_fill_whitebox_c(
    palette = "bl_yl_rd",
    labels = scales::label_number(suffix = "º"),
    n.breaks = 20,
    guide = guide_legend(reverse = TRUE)
  ) +
  geom_spatvector(data = aus, fill = NA) +
  theme_void() +
  labs(fill = "", title = "Max. temperature of wettest quarter")

ggsave(filename=paste0(plot_dir, "/tmax.pdf"), plot=tmax_plot, width=10, height=6.5) 

# -- tmin
temp_min <- lapply(models, function(m){
  lapply(ssp_scenarios, function(s){
    r <- rast(paste0(file.path(input_dir, m, s),"/temp_min.tif"))
    names(r) <- s
    r = mask(r, aus)
    r
  })
})
temp_min_r <- list(rast(temp_min[[1]]), rast(temp_min[[2]]), rast(temp_min[[3]]))
temp_min_r <-  c(temp_min_r[[1]], temp_min_r[[2]], temp_min_r[[3]])
names(temp_min_r) <- apply(expand.grid(ssp_scenarios, models), 1, paste, collapse=".")
#aus <- crop(aus, temp_min_r[[1]])

tmin_plot <- ggplot() +
  geom_spatraster(data = temp_min_r) +
  facet_wrap(~ lyr, ncol = 4) +
  scale_fill_whitebox_c(
    palette = "bl_yl_rd",
    labels = scales::label_number(suffix = "º"),
    n.breaks = 20,
    guide = guide_legend(reverse = TRUE)
  ) +
  geom_spatvector(data = aus, fill = NA) +
  theme_void() +
  labs(fill = "", title = "Min. temperature of wettest quarter")

ggsave(filename=paste0(plot_dir, "/tmin.pdf"), plot=tmax_plot, width=10, height=6.5) 

# -- srad
solar_rad <- lapply(models, function(m){
  lapply(ssp_scenarios, function(s){
    r <- rast(paste0(file.path(input_dir, m, s),"/solar_rad.tif"))
    names(r) <- s
    r = mask(r, aus)
    r
  })
})
solar_rad_r <- list(rast(solar_rad[[1]]), rast(solar_rad[[2]]), rast(solar_rad[[3]]))
solar_rad_r <-  c(solar_rad_r[[1]], solar_rad_r[[2]], solar_rad_r[[3]])
names(solar_rad_r) <- apply(expand.grid(ssp_scenarios, models), 1, paste, collapse=".")
#aus <- crop(aus, temp_min_r[[1]])

srad_plot <- ggplot() +
  geom_spatraster(data = solar_rad_r) +
  facet_wrap(~ lyr, ncol = 4) +
  scale_fill_grass_c(
    palette = "inferno",
    #labels = scales::label_number(suffix = "º"),
    n.breaks = 20,
    guide = guide_legend(reverse = TRUE)
  ) +
  geom_spatvector(data = aus, fill = NA) +
  theme_void() +
  labs(fill = "", title = expression('Solar radiation of wettest quarter (W m'^-2*'h'^-1*')'))

ggsave(filename=paste0(plot_dir, "/srad.pdf"), plot=srad_plot, width=10, height=6.5) 

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
hours_of_high_RH_r <-  c(hours_of_high_RH_r[[1]], hours_of_high_RH_r[[2]], hours_of_high_RH_r[[3]])
names(hours_of_high_RH_r) <- apply(expand.grid(ssp_scenarios, models), 1, paste, collapse=".")
#aus <- crop(aus, temp_min_r[[1]])

RH_plot <- ggplot() +
  geom_spatraster(data = hours_of_high_RH_r) +
  facet_wrap(~ lyr, ncol = 4) +
  scale_fill_whitebox_c(
    palette = "deep",
    #labels = scales::label_number(suffix = "º"),
    n.breaks = 12,
    guide = guide_legend(reverse = TRUE)
  ) +
  geom_spatvector(data = aus, fill = NA) +
  theme_void() +
  labs(fill = "", title = "Hours of high RH (of wettest quarter)")

ggsave(filename=paste0(plot_dir, "/hours_of_high_RH.pdf"), plot=RH_plot, width=10, height=6.5) 

"----------- Plot infection risk maps ---------------"
output_dir = 'Output/Climate/CMIP6'

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
infection_risk_tmin_r <-  c(infection_risk_tmin_r[[1]], infection_risk_tmin_r[[2]], infection_risk_tmin_r[[3]])
names(infection_risk_tmin_r) <- apply(expand.grid(ssp_scenarios, models), 1, paste, collapse=".")

risk_tmin_plot <- ggplot() +
  geom_spatraster(data = infection_risk_tmin_r) +
  facet_wrap(~ lyr, ncol = 4) +
  scale_fill_whitebox_c(
    palette = "viridi",
    n.breaks = 10,
    guide = guide_legend(reverse = TRUE)
  ) +
  geom_spatvector(data = aus, fill = NA) +
  theme_void() +
  labs(fill = "", title = "Infection risk at min. temperature")

ggsave(filename=paste0(plot_dir, "/risk_tmin.pdf"), plot=risk_tmin_plot, width=10, height=6.5) 

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
infection_risk_tmax_r <-  c(infection_risk_tmax_r[[1]], infection_risk_tmax_r[[2]], infection_risk_tmax_r[[3]])
names(infection_risk_tmax_r) <- apply(expand.grid(ssp_scenarios, models), 1, paste, collapse=".")

risk_tmax_plot <- ggplot() +
  geom_spatraster(data = infection_risk_tmax_r) +
  facet_wrap(~ lyr, ncol = 4) +
  scale_fill_whitebox_c(
    palette = "viridi",
    n.breaks = 10,
    guide = guide_legend(reverse = TRUE)
  ) +
  geom_spatvector(data = aus, fill = NA) +
  theme_void() +
  labs(fill = "", title = "Infection risk at max. temperature")

ggsave(filename=paste0(plot_dir, "/risk_tmax.pdf"), plot=risk_tmax_plot, width=10, height=6.5) 


# climate suitability
threshold_val = 0.2
suitability <- lapply(models, function(m){
  lapply(ssp_scenarios, function(s){
    r <- rast(paste0(file.path(output_dir, m, s),"/suitability_", threshold_val, ".tif"))
    names(r) <- s
    r = mask(r, aus)
    r
  })
})
suitability_r <- list(rast(suitability[[1]]), rast(suitability[[2]]), rast(suitability[[3]]))
suitability_r <-  c(suitability_r[[1]], suitability_r[[2]], suitability_r[[3]])
names(suitability_r) <- apply(expand.grid(ssp_scenarios, models), 1, paste, collapse=".")

suitability_plot <- ggplot() +
  geom_spatraster(data = as.factor(suitability_r)) +
  facet_wrap(~ lyr, ncol = 4) +
  scale_fill_manual(values = c("grey90", "black"), na.value = "transparent") +
  geom_spatvector(data = aus, fill = NA) +
  theme_void() +
  labs(fill = "", title = paste0("Climate suitability at risk = ", threshold_val))

ggsave(filename=paste0(plot_dir, "/suitability_", threshold_val,".pdf"), plot=suitability_plot, width=10, height=6.5) 
