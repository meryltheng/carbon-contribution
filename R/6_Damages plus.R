# Damages by MR accounting for the following:
# (i) climate suitability mask and 
# (ii) confidence intervals (from ABC samples of species impact)

library(tidyverse)
library(sf)
library(ggplot2)
library(reshape2)

"----------- Compute per ha reduction proportion (%) by MVG from ABC samples ---------------"
target_dir = 'Python/MR_samples'
load("Input/biomass_seq.RData") # load processed data for biomass sequestered
load("Input/area_aff_mvg_suit_0.2.RData") # load area affected by MR for each MVG (based on climate suitability)
output_table <- read.csv(file = "Output/Damage_estimates.csv")  

damage_samples <- list.files(path = target_dir, pattern = "*.csv", full.names = T) %>%
  map(function(x){
    # aggregate species impacts by genus
    impact_by_genus <- read.csv(x) %>%
      group_by(Genus) %>%
      summarise(Impact = mean(Damage, na.rm = T)) %>%
      filter(!is.na(Impact)) %>%
      rename(genus = Genus)
    # calculate impact/reduction to biomass sequestered by MVG
    impact_mvg <- left_join(biomass_seq_mvg, impact_by_genus, by = 'genus') %>% # Assign based on genus match 
      mutate(Impact = if_else(is.na(Impact), 0, Impact, NA)) %>% # if no info, impact = 0
      mutate(Impact_biomass = t_biomass_seq_ha * Impact) %>% # calc absolute impact to C-seq
      group_by(mvgValue) %>%
      summarise(Impact_biomass = sum(Impact_biomass), 
                t_biomass_seq_ha = sum(t_biomass_seq_ha),
                Impact_prop = sum(Impact_biomass) / sum(t_biomass_seq_ha) ) %>% # impact proportion by MVG
      st_drop_geometry() %>%
      rename(MVG = mvgValue) %>%
      complete(MVG = 1:32, fill = list(Impact_prop = 0)) 
    # compute damages $$
    damage_dat <- data.frame(MVG = output_table$MVG,
                             Biomass_seq = impact_mvg$t_biomass_seq_ha, # total seq. per ha
                             Biomass_lost = impact_mvg$Impact_biomass,  # total seq. loss per ha
                             Impact_perc = impact_mvg$Impact_prop, # Biomass_lost/Biomass_seq
                             Value = output_table$Value,
                             Area_total = output_table$Area,
                             Area_affe = area_aff_mvg) %>%
      mutate(Value_lost = Impact_perc * Value * Area_affe,
             Total_value = Value * Area_total)  
    return(damage_dat)
  })


"----------- Summary stats ---------------"
# biomass lost by MVG
bioloss_samples <- lapply(damage_samples, function(x) return(x$Biomass_lost) )
bioloss_samples <- do.call(cbind, bioloss_samples)
bioloss_summ <- apply(bioloss_samples, 1, function(x){
  median_est <- round(median(x, na.rm = T)/1000, 2) # in tons
  low_est <- round(quantile(x, 0.025, na.rm = T)/1000, 2)
  upp_est <- round(quantile(x, 0.925, na.rm = T)/1000, 2)
  paste0(median_est, " [", low_est, ",", upp_est, "]")
})

# impact % by MVG
impact_samples <- lapply(damage_samples, function(x) return(x$Impact_perc) )
impact_samples <- do.call(cbind, impact_samples)
impact_summ <- apply(impact_samples, 1, function(x){
  median_est <- round(median(x, na.rm = T)*100, 1)
  low_est <- round(quantile(x, 0.025, na.rm = T)*100, 1)
  upp_est <- round(quantile(x, 0.925, na.rm = T)*100, 1)
  paste0(median_est, " [", low_est, ",", upp_est, "]")
})

# impact_samples_rshp <- melt(impact_samples)
# names(impact_samples_rshp) <- c("MVG", "Sample", "Impact")
# ggplot(impact_samples_rshp, aes(x = Impact)) +
#   geom_density(alpha = 0.7) + 
#   facet_wrap(~MVG)

# value lost by MVG
value_lost_samples <- lapply(damage_samples, function(x) return(x$Value_lost) )
value_lost_samples <- do.call(cbind, value_lost_samples)
value_lost_summ <- apply(value_lost_samples, 1, function(x){
  median_est <- round(median(x, na.rm = T)/1000000)
  low_est <- round(quantile(x, 0.025, na.rm = T)/1000000)
  upp_est <- round(quantile(x, 0.925, na.rm = T)/1000000)
  paste0(median_est, " [", low_est, ",", upp_est, "]")
}) # in $ mil

# value lost TOTAL
val_lost_total <- colSums(value_lost_samples)/1000000 
val_lost_est <- c(round(median(val_lost_total)), round(quantile(val_lost_total,0.025)), round(quantile(val_lost_total,0.925)))
( DAM_R = val_lost_est / sum(output_table$Total_value, na.rm = T) * 100 ) # %

"----------- Prepare output table for paper ---------------"
output_table2 <- data.frame(MVG = output_table$MVG,
                            Biomass_seq_lost = bioloss_summ, # tons
                           Impact_perc = impact_summ,
                           Value = output_table$Value,
                           Area_total = output_table$Area,
                           Value_total = output_table$Total_value,
                           Area_aff = area_aff_mvg,
                           Value_lost = value_lost_summ) 

write.csv(output_table2, file = "Output/Damage_estimates_plus.csv", row.names = F) 
