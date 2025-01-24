library(sf)
library(terra)
library(geodata)
nvis <- rast('Data/Spatial/NVIS_MVG_RAS_AA.tif')

"----------- Download, reproject and crop climate layers ---------------"
target_dir = 'Data/Spatial/Climate/'
#options(timeout = max(60*45, getOption("timeout"))) # change timeout to 60*15 instead of 60 seconds

# maximum temperature (°C) -- done
tmax <- worldclim_country('tmax', country = "AUS", res = 0.5, path = tempdir())
tmax <- project(tmax, crs(nvis))
tmax <- crop(tmax, nvis)
writeRaster(tmax, paste0(target_dir, 'tmax.tif'), overwrite=TRUE) # not using saveRDS as readRDS will result in error

# minimum temperature (°C) -- done
tmin <- worldclim_country('tmin', country = "AUS", res = 0.5, path = tempdir())
tmin <- project(tmin, crs(nvis))
tmin <- crop(tmin, nvis)
writeRaster(tmin, paste0(target_dir, 'tmin.tif'), overwrite=TRUE)

# total precipitation (mm) -- done
prec <- worldclim_country('prec', country="AUS", res = 0.5, path = tempdir())
prec <- project(prec, crs(nvis))
prec <- crop(prec, nvis)
writeRaster(prec, paste0(target_dir, 'prec.tif'), overwrite=TRUE)

# vapour pressure (kPa) -- done
vapr <- worldclim_global(var="vapr", res=0.5, path= tempdir()) # not avail just for AUS; need to get global
nvis_ll <- project(nvis, crs(vapr)) # project to same crs as srad (cuz projecting global dataset would take a lot longer)
vapr_ll <- crop(vapr, nvis_ll) # crop in latlong projection
vapr <- project(vapr_ll, crs(nvis)) # then project back to nvis crs
vapr <- crop(vapr, nvis) # crop again cuz extents not equal
#plot(vapr[[1]]); ext(vapr)
writeRaster(vapr, paste0(target_dir, 'vapr.tif'), overwrite=TRUE)

# solar radiation	(kJ m-2 day-1) -- done
srad <- worldclim_country(var = "srad", country = "AUS", res = 0.5, path = tempdir())
srad < -project(srad, crs(nvis))
srad <- crop(srad, nvis)
writeRaster(srad, paste0(target_dir, 'srad.tif'), overwrite=TRUE)

"----------- Process climate layers into input layers we need ---------------"
# i.e., climate variables for the wettest quarter
# do this in native 30s resolution rather than the nvis 1 ha resolution

# load cropped climate layers (if done)
tmax <- rast(paste0(target_dir, 'tmax.tif')) 
tmin <- rast(paste0(target_dir, 'tmin.tif')) 
prec <- rast(paste0(target_dir, 'prec.tif')) 
vapr <- rast(paste0(target_dir, 'vapr.tif')) 
srad <- rast(paste0(target_dir, 'srad.tif')) 

# first determine wettest quarter (for each cell) using precip layer

# -- wettest quarter calculation
wet <- roll(prec, n=3, "sum", "from", na.rm=TRUE, circular=TRUE) # 0.5 min!!
wetqrt <- app(wet, which.max)  # wettest quarter for each cell

# (I) solar radiation: convert units to W m-2 h-1 
srad <- srad * 0.27778/24 # 1 kJ = 0.27778 Watt hours
sradwet <- selectRange(srad, wetqrt)
#writeRaster(sradwet, paste0(target_dir, 'input_srad.tif'), overwrite=TRUE)

# (II) min, max temperature for the wettest quarter
# mean of min temperature
tmpmin <- roll(tmin, n=3, "mean", "from", na.rm=TRUE, circular=TRUE) 
tminwet <- selectRange(tmpmin, wetqrt)
#writeRaster(tminwet, paste0(target_dir, 'tminwet.tif'), overwrite=TRUE)

# mean of max temperature
tmpmax <- roll(tmax, n=3, "mean", "from", na.rm=TRUE, circular=TRUE) 
tmaxwet <- selectRange(tmpmax, wetqrt)
#writeRaster(tmaxwet, paste0(target_dir, 'tmaxwet.tif'), overwrite=TRUE)

# min of min temperature
tmpmin2 <- roll(tmin, n=3, "min", "from", na.rm=TRUE, circular=TRUE) 
tminminwet <- selectRange(tmpmin2, wetqrt)
#writeRaster(tminminwet, paste0(target_dir, 'tminminwet.tif'), overwrite=TRUE)

# max of max temperature
tmpmax2 <- roll(tmax, n=3, "max", "from", na.rm=TRUE, circular=TRUE) 
tmaxmaxwet <- selectRange(tmpmax2, wetqrt)
#writeRaster(tmaxmaxwet, paste0(target_dir, 'tmaxmaxwet.tif'), overwrite=TRUE)

# (III) hours of high relative humidity (RH > 85%) during the wettest quarter
# -- first generate hourly temperature values from tmin and tmax (Eq 7 in https://doi.org/10.1177/0143624407078642)
generate_hourly_temp <- function(tmin, tmax, hour_tmin, hour_tmax, hour = 1:24){
  hourly_temp = (tmax - tmin)/2 * (cos(((hour - hour_tmax)*pi)/(hour_tmax - hour_tmin)) + 1) + tmin
  return(hourly_temp)
} # generate_hourly_temp(15, 30, 6, 15)
hourly_temp_wettest1 <- generate_hourly_temp(tminwet, tmaxwet, 6, 15, hour = 1:24)

# -- then calculate hourly saturation vapr pressure for hourly temperature values
calc_saturation_vp <- function(temperature) {
  0.611 * exp((17.27 * temperature) / (temperature + 237.3))
}
hourly_satr_vp <- calc_saturation_vp(hourly_temp_wettest1)

# -- calculate hourly RH from hourly saturation vapr pressure and vapr (actual vapr pressure)
# first resample vapr layer as extents/reso don't match
s <- rast(nrows=nrow(wetqrt), ncols=ncol(wetqrt),
          xmin=terra::ext(wetqrt)[1], xmax=terra::ext(wetqrt)[2],
          ymin=terra::ext(wetqrt)[3], ymax=terra::ext(wetqrt)[4])
vapr_resampled <- terra::resample(vapr, s, method="bilinear")

# then extract vapr pressure for wettest quarter
vaprwet <- selectRange(vapr_resampled, wetqrt)
#writeRaster(vaprwet, paste0(target_dir, 'vaprwet.tif'), overwrite=TRUE)

# finally, get the hourly RH
hourly_RH <- vaprwet/hourly_satr_vp

# -- get number of hours with RH >85%
hours_of_high_RH <- app(hourly_RH, function(x) length(which(x > 0.85)))
#writeRaster(hours_of_high_RH, paste0(target_dir, 'input_h_high_RH.tif'), overwrite=TRUE)

"----------- Create climate suitability layer ---------------"
# load input layers (if done)
tmin_wet <- rast(paste0(target_dir, 'tminwet.tif')) # or tminminwet.tif
tmax_wet <- rast(paste0(target_dir, 'tmaxwet.tif')) # or tmaxmaxwet.tif
srad_wet <- rast(paste0(target_dir, 'input_srad.tif')) 
hours_of_high_RH <- rast(paste0(target_dir, 'input_h_high_RH.tif'))

# (I) calculate infection risk based on processed climate input layers
#     using the formula provided by the NZ paper

# -- parameter values
C_P1 = -0.000095
C_P2 = 0.0012
C_P3 = 0.0635
C_P4 = 0.0047
M_P1 = 0.0011
M_P2 = 0.0917
M_P3 = -5.3193
M_P4 = 66.288
B = 0.39
temp_min = tmin_wet
temp_max = tmax_wet
solar_rad = srad_wet
hours_of_high_RH = hours_of_high_RH

# -- infection risk for temp_min
C_min = C_P1 * temp_min**3 + C_P2 * temp_min**2 + C_P3 * temp_min + C_P4
M_min = M_P1 * temp_min**3 + M_P2 * temp_min**2 + M_P3 * temp_min + M_P4
Y_min = C_min*exp(-exp(-B * (hours_of_high_RH - M_min)))
L = -0.0008*solar_rad + 1
infection_risk_tmin = Y_min * L

# -- infection risk for temp_max
C_max = C_P1 * temp_max**3 + C_P2 * temp_max**2 + C_P3 * temp_max + C_P4
M_max = M_P1 * temp_max**3 + M_P2 * temp_max**2 + M_P3 * temp_max + M_P4
Y_max = C_max*exp(-exp(-B * (hours_of_high_RH - M_max)))
infection_risk_tmax = Y_max * L

#writeRaster(infection_risk_tmin, paste0(target_dir, 'infection_risk_tmin.tif'), overwrite=TRUE)
#writeRaster(infection_risk_tmax, paste0(target_dir, 'infection_risk_tmax.tif'), overwrite=TRUE)

# -- plot
par(mfrow = c(1, 2))
plot(infection_risk_tmin, main='inf. risk. tmin')
plot(infection_risk_tmax, main='inf. risk. tmax') 

# (II) create climate suitability mask based on thresholds
threshold_val = 0.2 
clim_mask <- ifel(infection_risk_tmin > threshold_val | 
                    infection_risk_tmax > threshold_val, 1, 0) # 1: suitable, 0: unsuitable
plot(clim_mask)
#writeRaster(clim_mask, paste0(target_dir, 'clim_mask.tif'), overwrite=TRUE)

"----------- Resample climate suitability layer to resolution as nvis layer ---------------"
# global(nvis==32, sum, na.rm=TRUE) # use this instead to calculate area_ha in 4_Damages.R
clim_mask <- rast(paste0(target_dir, 'clim_mask.tif')) 

# create empty raster with the same extent and resolution as nvis layer
s <- rast(nrows=nrow(nvis), ncols=ncol(nvis),
          xmin=terra::ext(nvis)[1], xmax=terra::ext(nvis)[2],
          ymin=terra::ext(nvis)[3], ymax=terra::ext(nvis)[4])

# resample climate suitability layer to the resolution of nvis layer (or the other way around if faster?)
clim_mask_resampled <- terra::resample(clim_mask, s, method="near")
writeRaster(clim_mask_resampled, paste0(target_dir, 'clim_mask_resamp.tif'), overwrite=TRUE)

"----------- Calculate area within each MVG that is climate suitable ---------------"
# load climate mask layer (if done)
clim_mask_resampled <- rast(paste0(target_dir, 'clim_mask_resamp.tif'))

# create masked nvis layer
nvis_masked = nvis
nvis_masked[clim_mask_resampled==0] = NA # set all cells that are not suitable to NA

# calculate area in ha for each MVG that is climate suitable (takes a while)
area_aff_mvg <- sapply(1:32, function(x) global(nvis_masked==x, sum, na.rm=TRUE)) # num of cells per MVG = area in ha since res = 100m
area_aff_mvg <- do.call(c, area_aff_mvg)
#save(area_aff_mvg, file = "Input/area_aff_mvg_suit_0.2.RData") 
