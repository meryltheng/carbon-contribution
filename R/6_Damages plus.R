# Damages by MR accounting for the following:
# (i) climate suitability mask and 
# (ii) confidence intervals (from ABC samples of species impact)

library(tidyverse)
library(sf)
library(terra)
library(ggplot2)
library(reshape2)

# Define function to Calculate area within each MVG that is climate suitable
calc_mvg_area_aff <- function(nvis_layer, # mvg 
                              clim_mask){ # climate suitability mask
  
  # project climate layer to crs of nvis layer
  clim_mask <- project(clim_mask, crs(nvis_layer)) 
  
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

# Define function to Compute per ha reduction proportion (%) by MVG from ABC samples
# Output -- do we want just a single damage value AND the damage range (for ABC samples) for each climate scenario?
calc_damages <- function(area_aff_mvg = area_aff_mvg,
                         biomass_library = biomass_seq_mvg, # biomass sequestration data (per ha? species/genus level?)
                         impact_dir = 'Python/MR_samples',
                         mvg_table = mvg_info # general MVG info (name, Value_ha, Area_total)
){
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

"----------- Load and prepare inputs and data ---------------"
target_dir = 'Python/MR_samples' # impact samples directory
output_dir = 'Output/Climate/CMIP6' # climate mask directory
load("Input/biomass_seq.RData") # load processed data for biomass sequestered
output_table <- read.csv(file = "Output/Damage_estimates.csv")  
nvis <- rast('Data/Spatial/NVIS_MVG_RAS_AA.tif') # nvis mvg layer

# prepare global mvg info
mvg_info <- data.frame(MVG = output_table$MVG,
                       Area_total = output_table$Area, # Total area in ha
                       Value_ha = output_table$Value # Value per ha
) |>
  mutate(Value_total = Value_ha * Area_total) # Total value (Value_ha*Area_total)

biomass_seq_ha <- biomass_seq_mvg |>
  group_by(mvgValue) |>
  summarise(t_biomass_seq_ha = sum(t_biomass_seq_ha)) |> # biomass seq per ha
  st_drop_geometry() %>%
  rename(MVG = mvgValue) %>%
  complete(MVG = 1:32, fill = list(t_biomass_seq_ha = 0)) |>
  select(t_biomass_seq_ha) |>
  as.vector() |>
  unname() |>
  unlist()

mvg_info <- mvg_info |>
  mutate(C_ha = biomass_seq_ha * 0.5) |> # C sequestered per ha (kg)
  mutate(C_total = C_ha * Area_total) # Total C seq. (C_ha*Area_total))
# mvg_info %>%
#   # scenario-specific
#   mutate(Area_aff = area_aff_mvg) %>%
#   # impact-sample-specific
#   mutate(Impact_perc,
#          C_lost_ha,
#          C_lost_total,
#          Value_lost_total) # damages

# define models and ssp_scenarios to loop over
models = c("ACCESS-CM2", "CNRM-CM6-1-HR", "MPI-ESM1-2-LR")
ssp_scenarios = c("historical", "ssp126", "ssp245", "ssp370")
threshold_vals = c(0.2, 0.5)

"----------- Compute per ha impacts (%) and damages (value lost) by MVG from ABC samples ---------------"
# create output storage objects
damages <- list() # to contain mvg-level info
damage_total <- expand.grid( # to contain overall $$ damage
  model = models,
  ssp_scenario = ssp_scenarios,
  threshold_val = threshold_vals,
  sample_id = 1:1000
) |>
  mutate(Damage = NA)

# execute across all scenarios
for(m in models){
  damages[[m]] <- list()
  
  for(s in ssp_scenarios){
    damages[[m]][[s]] <- list()
    
    for(threshold in threshold_vals){
      print(paste('Damages for', m, s, threshold, 'starting...'))
      
      # read climate suitability mask for threshold value = 0.2 or 0.5
      clim_mask = rast(paste0(file.path(output_dir, m, s),"/suitability_", threshold, ".tif"))
      
      # calculate mvg area affected (climate suitable)
      area_aff_mvg <- calc_mvg_area_aff(nvis_layer = nvis, clim_mask)
      
      # calculate damages for each impact sample (from ABC)
      message('calculating damages for each ABC impact sample')
      dmg_table <- calc_damages(area_aff_mvg = area_aff_mvg,
                                biomass_library = biomass_seq_mvg, # biomass sequestration data (per ha? species/genus level?)
                                impact_dir = 'Python/MR_samples',
                                mvg_table = mvg_info # general MVG info (name, Value_ha, Area_total)
      )
      
      # append results
      damages[[m]][[s]][[paste(threshold)]] <- dmg_table
      Value_lost_total <- unlist(lapply(dmg_table, function(x) sum(x$Value_lost_total, na.rm=T)))
      damage_total[damage_total$model==m & 
                     damage_total$ssp_scenario==s & 
                     damage_total$threshold_val==threshold, ]$Damage <- Value_lost_total
      
      # clear mem
      rm(dmg_table, clim_mask, area_aff_mvg, Value_lost_total); gc()
      
      print(paste('Damages for', m, s, threshold, 'done.'))
    }
    
  }
}

# save
write.csv(damage_total, file = "Output/Damage_plus_plot.csv", row.names = F)
save(damages, file = "Output/Damage_plus_data.RData")
#write.csv(damage_total, file = "Output/Damage_plus_table.csv", row.names = F)

# plot
library(ggpubr)
threshold_cols <- c("grey75", "grey10")

damage_total$Damage <- damage_total$Damage/1000000 # in $m
damage_total$threshold_val <- as.factor(damage_total$threshold_val)

ggerrorplot(damage_total, x = "ssp_scenario", y = "Damage", color = "threshold_val", facet.by = 'model',
         palette = threshold_cols, desc_stat = "mean_sd", 
         size = 0.3, width = 0.5, add.params = list(size=0.1),
         position = position_dodge(0.3)) + 
  ggtitle('') + xlab('SSP scenario') + ylab('Damages ($m)') + 
  theme(strip.background = element_blank(), 
        strip.text.x = element_text(size=11))

ggviolin(damage_total, x = "ssp_scenario", y = "Damage"/1000000, fill = "threshold_val", facet.by = 'model',
         palette = threshold_cols, size = 0.3, add = "mean_sd", width = 0.5, add.params = list(size=0.1)) + 
  ggtitle('') + xlab('SSP scenario') + ylab('Damages ($m)') + 
  theme(strip.background = element_blank(), 
        strip.text.x = element_text(size=11))
