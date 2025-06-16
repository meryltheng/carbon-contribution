library(httr)
library(terra)
library(dplyr)
source("R/functions/climate.R")

'------ Download the NetCDF subset ------'
target_dir = 'Data/Spatial/Climate/CMIP6'

# models
climate_models <- data.frame(model = c("ACCESS-CM2", "CNRM-CM6-1-HR", "MPI-ESM1-2-LR"),
                             ccam_ver = c("v2112", "v2112", "v2105"),
                             ocean = c(TRUE, FALSE, FALSE),
                             ensemble = c("r2i1p1f1", "r1i1p1f2", "r9i1p1f1"))

# variables and ssp scenarios
variable_ids <- c("tas", "tasmin", "tasmax", "pr", "rsds", "hurs")
ssp_scenarios = c("historical", "ssp126", "ssp245", "ssp370")

# constants
temporal_reso = "mon"
lat_range = c(-9.1, -44.25)
lon_range = c(112.8, 155)

for (m in 1:nrow(climate_models)){
  
  model = climate_models$model[m]
  
  for(s in ssp_scenarios){
    
    if(s == "historical"){
      time_start="2014-01-16T12:00:00Z" # for the year 2014
      time_end="2014-12-16T12:00:00Z"
    } else{
      time_start="2100-01-16T12:00:00Z" # for the year 2100
      time_end="2100-12-16T12:00:00Z"
    }
    
    # create folder to contain data
    ifelse(!dir.exists(file.path(target_dir, model)), # model-level dir
           dir.create(file.path(target_dir, model)), 
           FALSE) 
    ifelse(!dir.exists(file.path(target_dir, model, s)), # ssp-level dir
           dir.create(file.path(target_dir, model, s)), 
           FALSE) 
    
    for(i in variable_ids){
      # Construct the output filename based on the URL (this can be customized)
      output_file <- paste0(file.path(target_dir, model, s), "/", i, ".nc") # paste0(target_dir, i, "_", s, "_2100", ".nc")
      
      # Download and subset data for each URL
      download_subset_data(variable_id = i, ssp_scenario = s, model = model, ccam_ver = climate_models$ccam_ver[m],
                           ocean = climate_models$ocean[m], ensemble = climate_models$ensemble[m],
                           temporal_reso = temporal_reso, time_start=time_start, time_end=time_end,
                           lat_range = lat_range, lon_range = lon_range,
                           output_file = output_file)
    }
  }
}

"----------- Climate suitability calculation -------------"
# define models and ssp_scenarios to loop over
models = c("ACCESS-CM2", "CNRM-CM6-1-HR", "MPI-ESM1-2-LR")
ssp_scenarios = c("historical", "ssp126", "ssp245", "ssp370")

for(m in models){
  for(s in ssp_scenarios){
    
    # process and save climate layers for model input
    prep_clim_vars(model = m, ssp_scenario = s, 
                   data_dir = 'Data/Spatial/Climate/CMIP6', 
                   output_dir = 'Input/Climate/CMIP6')
    
    # calculate infection pressure from processes climate inputs
    calc_inf_pressure(model = m, ssp_scenario = s, 
                      input_dir = 'Input/Climate/CMIP6', 
                      output_dir = 'Output/Climate/CMIP6')
  }
}

for(m in models){
  for(s in ssp_scenarios){
    # generate climate suitability mask based on threshold value
    # 0.2
    generate_suitability(model = m, ssp_scenario = s, threshold_val = 0.2,
                         target_dir = 'Output/Climate/CMIP6')
    # 0.5
    generate_suitability(model = m, ssp_scenario = s, threshold_val = 0.5,
                         target_dir = 'Output/Climate/CMIP6')
  }
}

