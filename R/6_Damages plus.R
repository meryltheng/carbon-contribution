# Damages by MR accounting for the following:
# (i) climate suitability mask and 
# (ii) confidence intervals (from ABC samples of species impact)

library(tidyverse)
library(sf)
library(terra)
library(ggplot2)
library(reshape2)
source("R/functions/damages.R")

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

"----------- Prepare output tables for paper ---------------"
# load sample results if processed
load("Output/Damage_plus_data.RData")
damage_total <- read.csv("Output/Damage_plus_plot.csv")

fixed_vars <- damages[[1]][[1]][[1]][[1]][,-(7:11)] %>%
  mutate(Value_total = Value_total/1000000, # convert to $m
         C_ha = C_ha/1000, # convert to tons
         C_total = C_total/1000) # convert to tons

# summarise damages by MVG across all scenarios
output_tables <- list()

for(m in models){
  output_tables[[m]] <- list()
  
  for(s in ssp_scenarios){
    output_tables[[m]][[s]] <- list()
    
    for(threshold in threshold_vals){
      df_scenario = damages[[m]][[s]][[paste(threshold)]]
      
      out <- cbind(fixed_vars,
                   data.frame(Area_aff = df_scenario[[1]]$Area_aff,
                              Impact_perc = summarise_var(df = df_scenario, var_name = "Impact_perc", conversion = "*100", rounding = 1),
                              C_lost_ha  = summarise_var(df = df_scenario, var_name = "C_lost_ha", conversion = "/1000", rounding = 2), # in tons
                              C_lost_total = summarise_var(df = df_scenario, var_name = "C_lost_total", conversion = "/1000", rounding = 0), # in tons
                              Value_lost_total = summarise_var(df = df_scenario, var_name = "Value_lost_total", conversion = "/1000000", rounding = 0)) # in $mil
      )
      
      output_tables[[m]][[s]][[paste(threshold)]] <- out
      rm(out, df_scenario); gc()
    }
  }
}

save(output_tables, file = "Output/Damage_plus_tables.RData")

# outout list to excel file
library(openxlsx)

wb <- createWorkbook()
for(m in models){
  for(s in ssp_scenarios){
    for(threshold in threshold_vals){
      addWorksheet(wb, sheetName = paste0(m, "_", s, "_", threshold))
      writeData(wb = wb, 
                sheet = paste0(m, "_", s, "_", threshold), 
                x = output_tables[[m]][[s]][[paste(threshold)]])
    }
  }
}

saveWorkbook(wb, file = "Output/Damage_plus_tables.xlsx", overwrite = TRUE)

# summarise total damages by model, ssp_scenario, threshold_val
damage_summ <- damage_total %>%
  group_by(model, ssp_scenario, threshold_val) %>%
  summarise(median_damage = median(Damage)/1000000,
            lo_damage = quantile(Damage, 0.025)/1000000,
            upp_damage = quantile(Damage, 0.975)/1000000) 

total_val = sum(fixed_vars$Value_total) # in mil

damage_summ_table <- damage_summ %>%
  mutate(median_damage_perc = median_damage/total_val*100,
         lo_damage_perc = lo_damage/total_val*100,
         upp_damage_perc = upp_damage/total_val*100) %>%
  mutate(total_damages = paste0(round(median_damage, 2), " [", round(lo_damage, 2), ",", round(upp_damage, 2), "]"), # in $mil
         reduction = paste0(round(median_damage_perc, 2), " [", round(lo_damage_perc, 2), ",", round(upp_damage_perc, 2), "]")) %>% # in %
  select(-median_damage, -lo_damage, -upp_damage, -median_damage_perc, -lo_damage_perc, -upp_damage_perc) 

write.csv(damage_summ_table, file = "Output/Damage_plus_summary.csv", row.names = F)

# C lost total
C_total = sum(fixed_vars$C_total)/1000000 # in mil tons
C_total # 1062.292M tons

# append total damages to summary table
damage_summ_table$total_C_lost <- NA
damage_summ_table$reduction_c <- NA

for(m in models){
  for(s in ssp_scenarios){
    for(threshold in threshold_vals){
      df_scenario = damages[[m]][[s]][[paste(threshold)]]
      
      C_lost <- unlist(lapply(df_scenario, function(x) return(sum(x$C_lost_total, na.rm=T))))
      C_lost_summ <- paste0(round(median(C_lost)/1000/1000000, 1), " [", # in mil tons
                            round(quantile(C_lost, 0.025, na.rm=T)/1000/1000000, 1), ",", 
                            round(quantile(C_lost, 0.975, na.rm=T)/1000/1000000, 1), "]")
      
      C_lost_perc <- (C_lost/1000/1000000)/C_total
      C_lost_perc_summ <- paste0(round(median(C_lost_perc)*100, 1), " [", # in %
                                 round(quantile(C_lost_perc, 0.025)*100, 1), ",", 
                                 round(quantile(C_lost_perc, 0.975)*100, 1), "]")
      
      damage_summ_table[damage_summ_table$model==m & 
                          damage_summ_table$ssp_scenario==s & 
                          damage_summ_table$threshold_val==threshold, ]$total_C_lost <- C_lost_summ
      
      damage_summ_table[damage_summ_table$model==m & 
                          damage_summ_table$ssp_scenario==s & 
                          damage_summ_table$threshold_val==threshold, ]$reduction_c <- C_lost_perc_summ
      # clear mem
      rm(df_scenario, C_lost, C_lost_summ, C_lost_perc, C_lost_perc_summ); gc()
    }
  }
}

write.csv(damage_summ_table, file = "Output/Damage_plus_summary.csv", row.names = F)

