#===========================================
# Calculate area within each MVG that is climate suitable
calc_mvg_area_aff <- function(nvis_layer, # mvg 
                              clim_mask){ # climate suitability mask
  require(terra)
  
  # project climate layer to crs of nvis layer
  clim_mask <- project(clim_mask, crs(nvis_layer), method = "near") 
  
  # crop climate layer to nvis layer to match extents
  clim_mask <- crop(clim_mask, nvis_layer) 
  
  # create empty raster with the same extent and resolution as nvis layer
  s <- rast(nrows=nrow(nvis_layer), ncols=ncol(nvis_layer),
            xmin=terra::ext(nvis_layer)[1], xmax=terra::ext(nvis_layer)[2],
            ymin=terra::ext(nvis_layer)[3], ymax=terra::ext(nvis_layer)[4])
  
  # Resample climate suitability layer to resolution as nvis layer
  message('climate mask resampling in progress...')
  clim_mask_resampled <- terra::resample(clim_mask, s, method="near")
  
  # create masked nvis layer
  nvis_masked = nvis_layer
  nvis_masked[clim_mask_resampled==0] = NA # set all cells that are not suitable to NA
  
  # calculate area in ha for each MVG that is climate suitable (takes a while)
  message('climate suitable area calculation in progress...')
  #area_aff_mvg <- sapply(1:32, function(x) global(nvis_masked==x, sum, na.rm=TRUE)) # num of cells per MVG = area in ha since res = 100m
  area_aff_mvg <- freq(nvis_masked)
  area_aff_mvg <- area_aff_mvg$count[1:32]
  
  return(area_aff_mvg)
}

#===========================================
# Compute per ha reduction proportion (%) by MVG from ABC samples
# (Output -- do we want just a single damage value AND the damage range (for ABC samples) for each climate scenario?)
calc_damages <- function(area_aff_mvg = area_aff_mvg,
                         biomass_library = biomass_seq_mvg, # biomass sequestration data (per ha? species/genus level?)
                         impact_dir = 'Python/MR_samples',
                         mvg_table = mvg_info # general MVG info (name, Value_ha, Area_total)
){
  require(tidyverse)
  require(sf)
  damage_samples <- list.files(path = impact_dir, pattern = "*.csv", full.names = T) %>%
    map(function(x){
      # aggregate species impacts by genus
      impact_by_genus <- read.csv(x) %>%
        group_by(Genus) %>%
        summarise(Impact = mean(Damage, na.rm = T)) %>%
        filter(!is.na(Impact)) %>%
        rename(genus = Genus)
      # calculate impact/reduction to biomass sequestered by MVG
      impact_mvg <- left_join(biomass_library, impact_by_genus, by = 'genus') %>% # Assign based on genus match 
        mutate(Impact = if_else(is.na(Impact), 0, Impact, NA)) %>% # if no info, impact = 0
        mutate(Impact_biomass = t_biomass_seq_ha * Impact) %>% # calc absolute impact to biomass per ha (kg)
        group_by(mvgValue) %>%
        summarise(Impact_biomass = sum(Impact_biomass), # summed biomass loss (per ha) across all genera (kg)
                  t_biomass_seq_ha = sum(t_biomass_seq_ha), # biomass seq per ha
                  Impact_prop = sum(Impact_biomass) / sum(t_biomass_seq_ha) ) %>% # impact proportion by MVG
        st_drop_geometry() %>%
        rename(MVG = mvgValue) %>%
        complete(MVG = 1:32, fill = list(Impact_prop = 0)) 
      # compute damages $$
      damage_dat <- mvg_info %>%
        # scenario-specific
        mutate(Area_aff = area_aff_mvg) %>%
        # impact-sample-specific
        mutate(Impact_perc = impact_mvg$Impact_prop,
               C_lost_ha = impact_mvg$Impact_biomass * 0.5) |> # total C seq. lost per ha (kg)
        mutate(C_lost_total = C_lost_ha * Area_aff, # total C seq. lost (kg)
               Value_lost_total = Impact_perc * Value_ha * Area_aff) # damages
      
      
      return(damage_dat)
    })
  return(damage_samples)
}

#===========================================
# Extract summary stats for variable of interest
summarise_var <- function(df = damages, var_name = "Impact_perc", conversion = "*100", rounding = 1){
  var_samples <- lapply(df, function(x) return(eval(parse(text=paste0("x$",var_name)))) )
  var_samples <- do.call(cbind, var_samples)
  var_summ <- apply(var_samples, 1, function(x){
    median_est <- round(eval(parse(text=paste0("median(x, na.rm = T)", conversion))), rounding)
    low_est <- round(eval(parse(text=paste0("quantile(x, 0.025, na.rm = T)", conversion))), rounding)
    upp_est <- round(eval(parse(text=paste0("quantile(x, 0.975, na.rm = T)", conversion))), rounding)
    paste0(median_est, " [", low_est, ",", upp_est, "]")
  })
  return(var_summ)
}