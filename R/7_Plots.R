library(terra)
library(geodata)
target_dir = 'Data/Spatial/Climate/'
nvis <- rast('Data/Spatial/NVIS_MVG_RAS_AA.tif')

# get country boundary
aus <- gadm(country="AUS", level=0, path=tempdir(), resolution=2)
aus <- project(aus, crs(nvis))
aus <- crop(aus, nvis)

"----------- Plot climate input layers ---------------"
# load input layers (if done)
tmin_wet <- rast(paste0(target_dir, 'tminwet.tif')) # or tminminwet.tif
tmax_wet <- rast(paste0(target_dir, 'tmaxwet.tif')) # or tmaxmaxwet.tif
srad_wet <- rast(paste0(target_dir, 'input_srad.tif')) 
hours_of_high_RH <- rast(paste0(target_dir, 'input_h_high_RH.tif'))

# plot climate variables
png(filename="Plots/climate_vars.png", width=1000, height=850, units="px", res = 72*2, pointsize = 12)
par(mfrow = c(2, 2), mar = c(1, 1, 1, 1))
plot(tmin_wet, main = expression("Min. temperature (°C)"), col = viridis::inferno(41)[1:28]); polys(aus, border="grey15", lwd=1, lty=1, alpha=0.5)
plot(tmax_wet, main = expression("Max. temperature (°C)"), col = viridis::inferno(41)); polys(aus, border="grey15", lwd=1, lty=1, alpha=0.5)
plot(srad_wet, main = expression('Solar radiation (W m'^-2*'h'^-1*')'), col = viridis::inferno(50)); polys(aus, border="grey15", lwd=1, lty=1, alpha=0.5)
plot(hours_of_high_RH, main = expression("Hours of high RH"), col = c("#FFFFFF", rev(viridis::mako(48))[1:24])); polys(aus, border="grey15", lwd=1, lty=1, alpha=0.5)
dev.off()

"----------- Plot infection risk maps ---------------"
infection_risk_tmin <- rast(paste0(target_dir, 'infection_risk_tmin.tif'))
infection_risk_tmax <- rast(paste0(target_dir, 'infection_risk_tmax.tif'))

png(filename="Plots/climate_suit.png", width=1000, height=850, units="px", res = 72*2, pointsize = 12)
par(mfrow = c(2, 2), mar = c(1, 1, 1, 1))
plot(infection_risk_tmin, main=expression('Infection risk at min. temp'), col = viridis::viridis(20)); polys(aus, border="grey15", lwd=1, lty=1, alpha=0.5)
plot(infection_risk_tmax, main=expression('Infection risk at max. temp'), col = viridis::viridis(20)); polys(aus, border="grey15", lwd=1, lty=1, alpha=0.5)
threshold_val = 0.2 
clim_mask <- ifel(infection_risk_tmin > threshold_val | 
                    infection_risk_tmax > threshold_val, 1, 0) # 1: suitable, 0: unsuitable
plot(clim_mask, main=expression('Suitable climate at 0.2'), col = c("#FFFFFF", "#000000")); polys(aus, border="grey15", lwd=1, lty=1, alpha=0.5)
dev.off()

"----------- Plot impact w and wo climate mask ---------------" # WIP
clim_mask_resampled <- rast(paste0(target_dir, 'clim_mask_resamp.tif'))

# # create value layer based on MVG classification
# output_table = read.csv(file = "Output/Damage_estimates.csv")
# rclmat <- matrix(c(1:32, output_table$Value), ncol = 2, byrow = F) # is-to-become
# value_r <- classify(mvg_r, rclmat)

# create asset reduction % layer based on MVG classification
impact_perc_median = as.numeric(unname(sapply(output_table2$Impact_perc, function(x) strsplit(x, split = " ")[[1]][1] )))
rclmat <- matrix(c(c(1:32,99), c(impact_perc_median, NA)), ncol = 2, byrow = F) # is-to-become
impact_r <- classify(nvis, rclmat)
plot(impact_r, col = c("#FFFFFF", heat.colors(22, rev = T)))

# create biomass reduction layer based on MVG classification
bioloss_median = as.numeric(unname(sapply(output_table2$Biomass_seq_lost, function(x) strsplit(x, split = " ")[[1]][1] )))
rclmat <- matrix(c(c(1:32,99), c(bioloss_median, NA)), ncol = 2, byrow = F) # is-to-become
bioloss_r <- classify(nvis, rclmat)
plot(bioloss_r, col = c("#FFFFFF", heat.colors(22, rev = T)))

# create new impact raster with impacts only where climate is suitable
impact_r_masked <- ifel(clim_mask_resampled == 1, impact_r, 0) 
bioloss_r_masked <- ifel(clim_mask_resampled == 1, bioloss_r, 0)

# plot damages
png(filename="Plots/damages.png", width=1000, height=850, units="px", res = 72*2, pointsize = 12)
par(mfrow = c(2, 2), mar = c(1, 1, 1, 1))
plot(impact_r, main = expression("Asset reduction (%)"), col = c("#FFFFFF", heat.colors(22, rev = T))); polys(aus, border="grey15", lwd=1, lty=1, alpha=0.5)
plot(bioloss_r, main = expression('Biomass seq. lost (t ha'^-1*')'), col = c("#FFFFFF", heat.colors(22, rev = T))); polys(aus, border="grey15", lwd=1, lty=1, alpha=0.5)
plot(impact_r_masked, main = expression("Asset reduction (%)"), col = c("#FFFFFF", heat.colors(22, rev = T))); polys(aus, border="grey15", lwd=1, lty=1, alpha=0.5)
plot(bioloss_r_masked, main = expression('Biomass seq. lost (t ha'^-1*')'), col = c("#FFFFFF", heat.colors(22, rev = T))); polys(aus, border="grey15", lwd=1, lty=1, alpha=0.5)
dev.off()
