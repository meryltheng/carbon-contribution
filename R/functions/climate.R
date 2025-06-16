#===========================================
# Download and subset data
# data from: https://www.longpaddock.qld.gov.au/qld-future-climate/data-info/tern-cmip6/
download_subset_data <- function(variable_id = "hurs", ssp_scenario = "ssp245", model = "ACCESS-CM2", 
                                 ccam_ver = "v2112", ocean = TRUE, ensemble = "r1i1p1f2",
                                 temporal_reso = "mon", time_start="2100-01-16T12:00:00Z", time_end="2100-12-16T12:00:00Z",
                                 lat_range = c(-9.1422, -52.350), lon_range = c(112.9211, 159.1092),
                                 output_file) {
  require(httr)
  
  # Construct the NCSS query parameters
  if(ocean == TRUE){ oc = "oc"} else{oc = NULL}
  ncss_url <- paste0("https://thredds.nci.org.au/thredds/ncss/grid/ig45/QldFCP-2/output/CMIP6/DD/AUS-10i/UQ-DEC/", 
                     model, "/", ssp_scenario, "/", ensemble, "/CCAM", oc, "-", ccam_ver, "/v1-r1/", temporal_reso, "/", 
                     variable_id, "/latest/", variable_id, "_AUS-10i_", model, "_", ssp_scenario, 
                     "_", ensemble, "_UQ-DEC_CCAM", oc, "-", ccam_ver, "_v1-r1_",
                     temporal_reso, "_", str_c(unlist(strsplit(time_start, "-"))[1:2], collapse =""), "-", 
                     str_c(unlist(strsplit(time_end, "-"))[1:2], collapse =""), 
                     ".nc?var=", variable_id, 
                     "&north=", lat_range[1], "&west=", lon_range[1], "&east=", lon_range[2], "&south=", lat_range[2],
                     "&horizStride=1&time_start=", time_start, "&time_end=", time_end, "&&&accept=netcdf3")
  
  print(ncss_url)
  
  # Send the HTTP request to NCSS service
  response <- GET(ncss_url)
  
  if (status_code(response) == 200) {
    # Save the NetCDF data to a file
    writeBin(content(response, "raw"), output_file)
    print(paste("Data successfully saved to", output_file))
  } else {
    print(paste("Failed to retrieve data"))
  }
}

#===========================================
# Generate hourly temperature values from tmin and tmax 
# (Eq 7 in https://doi.org/10.1177/0143624407078642)
generate_hourly_temp <- function(tmin, tmax, hour_tmin, hour_tmax, hour = 1:24){
  hourly_temp = (tmax - tmin)/2 * (cos(((hour - hour_tmax)*pi)/(hour_tmax - hour_tmin)) + 1) + tmin
  return(hourly_temp)
} # generate_hourly_temp(15, 30, 6, 15)

#===========================================
# Calculate saturation vapour pressure from temperature
calc_saturation_vp <- function(temperature) {
  0.611 * exp((17.27 * temperature) / (temperature + 237.3))
}

#===========================================
# Process raw climate data into model inputs (and save)
# (Converts temperature variables from K to C, calculates the wettest quarter 
# and extracts rolling means for each climate variable for wettest quarter, 
# generates hourly temperature values, and backs out relative humidity for 
# every hour.)
prep_clim_vars <- function(model, ssp_scenario, 
                           data_dir = 'Data/Spatial/Climate/CMIP6', 
                           output_dir = 'Input/Climate/CMIP6'){
  require(terra)
  
  # raw data directory path
  data_path = file.path(data_dir, model, ssp_scenario)
  
  # read all variables
  pr <- rast(paste0(data_path, "/pr.nc"))
  tas <- rast(paste0(data_path, "/tas.nc"))
  tasmin <- rast(paste0(data_path, "/tasmin.nc"))
  tasmax <- rast(paste0(data_path, "/tasmax.nc"))
  hurs <- rast(paste0(data_path, "/hurs.nc"))
  rsds <- rast(paste0(data_path, "/rsds.nc"))
  
  # convert temperature variables from K to C
  tas <- tas - 273.15
  tasmin <- tasmin - 273.15
  tasmax <- tasmax - 273.15
  
  # wettest quarter calculation
  wet <- roll(pr, n=3, "sum", "from", na.rm=TRUE, circular=TRUE) # quarterly rolling sum for precip
  wetqrt <- app(wet, which.max)  # wettest quarter for each cell 
  
  # calculate quarterly rolling means for all variables
  tas <- roll(tas, n=3, "mean", "from", na.rm=TRUE, circular=TRUE) 
  tasmin <- roll(tasmin, n=3, "mean", "from", na.rm=TRUE, circular=TRUE) 
  tasmax <- roll(tasmax, n=3, "mean", "from", na.rm=TRUE, circular=TRUE) 
  hurs <- roll(hurs, n=3, "mean", "from", na.rm=TRUE, circular=TRUE) 
  rsds <- roll(rsds, n=3, "mean", "from", na.rm=TRUE, circular=TRUE) 
  
  # extract quarterly rolling means (for all variables) for the wettest quarter
  tas_wet <- selectRange(tas, wetqrt)
  tasmin_wet <- selectRange(tasmin, wetqrt)
  tasmax_wet <- selectRange(tasmax, wetqrt)
  hurs_wet <- selectRange(hurs, wetqrt)
  rsds_wet <- selectRange(rsds, wetqrt)
  
  # approximate hourly hurs (relative humidity)
  tas_hour_wet <- generate_hourly_temp(tasmin_wet, tasmax_wet, 6, 15, hour = 1:24) # generate hourly temperature values
  esat_hour_wet <- calc_saturation_vp(tas_hour_wet) # hourly saturation vap. pressure 
  esat_day_wet <- calc_saturation_vp(tas_wet) # daily mean saturation vap. pressure (from mean temp `tas`)
  e_day_wet <- hurs / 100 * esat_day_wet # estimate the actual vapor pressure from hurs (mean daily)
  hurs_hour_wet <- e_day_wet/esat_hour_wet # get hourly RH
  hours_of_high_RH <- app(hurs_hour_wet, function(x) length(which(x > 0.85)))
  
  # output directory path
  output_path = file.path(output_dir, model, ssp_scenario)
  
  # create folder to contain processed data
  ifelse(!dir.exists(file.path(output_dir, model)), # model-level dir
         dir.create(file.path(output_dir, model)), 
         FALSE) 
  ifelse(!dir.exists(file.path(output_dir, model, ssp_scenario)), # ssp-level dir
         dir.create(file.path(output_dir, model, ssp_scenario)), 
         FALSE) 
  
  # save outputs
  writeRaster(tasmin_wet, paste0(output_path, "/temp_min.tif"), overwrite=TRUE)
  writeRaster(tasmax_wet, paste0(output_path, "/temp_max.tif"), overwrite=TRUE)
  writeRaster(rsds_wet, paste0(output_path, "/solar_rad.tif"), overwrite=TRUE)
  writeRaster(hours_of_high_RH, paste0(output_path, "/hours_of_high_RH.tif"), overwrite=TRUE)
  
  message(paste(m, s, "prep clim inputs done.")) 
}

#===========================================
# Model Myrtle Rust infection pressure from prepared climate inputs (and save)
calc_inf_pressure <- function(model, ssp_scenario, 
                              input_dir = 'Input/Climate/CMIP6', 
                              output_dir = 'Output/Climate/CMIP6'){ #temp_min, temp_max, solar_rad, hours_of_high_RH
  require(terra)
  
  # input directory path
  input_path = file.path(input_dir, model, ssp_scenario)
  
  # read all input vars
  temp_min = rast(paste0(input_path, "/temp_min.tif"))
  temp_max = rast(paste0(input_path, "/temp_max.tif"))
  solar_rad = rast(paste0(input_path, "/solar_rad.tif"))
  hours_of_high_RH = rast(paste0(input_path, "/hours_of_high_RH.tif"))
  
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
  
  # output directory path
  output_path = file.path(output_dir, model, ssp_scenario)
  
  # create folder to contain processed data
  ifelse(!dir.exists(file.path(output_dir, model)), # model-level dir
         dir.create(file.path(output_dir, model)), 
         FALSE) 
  ifelse(!dir.exists(file.path(output_dir, model, ssp_scenario)), # ssp-level dir
         dir.create(file.path(output_dir, model, ssp_scenario)), 
         FALSE) 
  
  # save outputs
  writeRaster(infection_risk_tmin, paste0(output_path, "/infection_risk_tmin.tif"), overwrite=TRUE)
  writeRaster(infection_risk_tmax, paste0(output_path, "/infection_risk_tmax.tif"), overwrite=TRUE)
  #return(list(infection_risk_tmin, infection_risk_tmax))
  
  message(paste(m, s, "inf risk calc done.")) 
}

#===========================================
# Generate climate suitability mask based on a specfied threshold of 
# infection risk
generate_suitability <- function(model = m, ssp_scenario = s, threshold_val = 0.2,
                                 target_dir = 'Output/Climate/CMIP6'){
  require(terra)
  
  # load the infection risk map
  riskmap_path = file.path(target_dir, model, ssp_scenario)
  infection_risk_tmin = rast(paste0(riskmap_path, '/infection_risk_tmin.tif'))
  infection_risk_tmax = rast(paste0(riskmap_path, '/infection_risk_tmax.tif'))
  
  # create climate suitability mask based on thresholds
  clim_mask <- ifel(infection_risk_tmin > threshold_val | 
                      infection_risk_tmax > threshold_val, 1, 0) # 1: suitable, 0: unsuitable
  
  # save output
  writeRaster(clim_mask, paste0(target_dir, "/suitability_", threshold_val, ".tif"), overwrite=TRUE)
  message(paste(m, s, "suitability map for threshold", threshold_val, "done.")) 
}