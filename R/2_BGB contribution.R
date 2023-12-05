library(data.table) # for fast reads of large datasets
library(dplyr)
library(sf)

"----------- Load processed data for AGB ---------------"
load("Input/agb_estimates.RData") # if processed

"----------- Assign record to functional plant type (for BGB) ---------------"
# Quick and dirty way to assign each record to a functional plant type for AGB calculations later

# a. First, assign based on https://doi.org/10.1016/j.foreco.2018.08.043 (Fig. 3)
all_data_bgb <- agb_estimates %>%
  mutate(func_type_bgb = case_when(
    func_type_agb == "F Shrub" ~ "FShrub&Ac",
    func_type_agb == "F Multi" & genus == "Acacia" ~ "FShrub&Ac",
    func_type_agb == "F Euc" ~ "FTree",
    func_type_agb == "F Multi" & genus == "Eucalyptus" ~ "FMallee",
    func_type_agb == "F Other-H" ~ "FTree",
    func_type_agb == "F Other-H" & species == "Pinus pinaster" ~ "FTree", # not present in current dataset
    func_type_agb == "F Other-L" & species == "Pinus radiata"~ "FRadiata" # not present in current dataset
  ))

# subset the no match records
no_match <- all_data_bgb %>%
  filter(is.na(func_type_bgb)) # THERE ARE NONE WOOHOO

"----------- Estimate BGB based on allometric eqns per functional type ---------------" 
# estimate BGB based on eqns given by https://doi.org/10.1016/j.foreco.2018.08.043
biomass_estimates <- all_data_bgb %>%
  mutate(bgb_estimate = case_when(
    func_type_bgb == "FShrub&Ac" ~ exp(-3.553 + 2.185 * log(diameter) * 1.160),
    func_type_bgb == "FMallee" ~ exp(-2.946 + 2.302 * log(diameter) * 1.116),
    func_type_bgb == "FTree" ~ exp(-2.682 + 2.212 * log(diameter) * 1.096),
    func_type_bgb == "FRadiata" ~ exp(-3.740 + 2.299 * log(diameter) * 1.053)
  )) %>%
  mutate(t_biomass_estimate = agb_estimate + bgb_estimate)

"----------- Biomass (AGB + BGB) per ha per genus per MVG ---------------" 
#load("Input/agb_estimates.RData") # if processed

biomass_est_std <- biomass_estimates %>%
  group_by(site, genus) %>%
  summarise(geometry = unique(geometry),
            site = unique(site),
            sitearea_ha = unique(sitearea_ha),
            genus = unique(genus),
            mvgValue = unique(mvgValue),
            mvsValue = unique(mvsValue),
            t_biomass_est = sum(t_biomass_estimate, na.rm = T), # total biomass estimate
            t_biomass_per_ha = sum(t_biomass_estimate, na.rm = T)/unique(sitearea_ha)) # standardised to per ha

biomass_est_mvg <- biomass_est_std %>%
  group_by(mvgValue) %>%
  mutate(n_sites_per_mvg = length(unique(site))) %>%
  group_by(mvgValue, genus) %>%
  summarise(n_sites_per_mvg = unique(n_sites_per_mvg),
            genus = unique(genus),
            mvgValue = unique(mvgValue),
            t_biomass_per_ha = mean(t_biomass_per_ha, na.rm = T)) # standardised to per ha


#save(biomass_estimates, biomass_est_std, biomass_est_mvg, file = "Input/biomass_estimates.RData")

"----------- Visualise! ---------------"
library(ggplot2)
load("Input/biomass_estimates.RData") # if saved

ggplot(biomass_est_mvg, aes(x = as.factor(mvgValue), y = t_biomass_per_ha, fill = genus)) +
  geom_bar(stat = 'identity',color='black') + # position='fill' if want proportion
  labs(x = "MVG", y = "Biomass contribution per ha (kg)") +
  #geom_text(aes(label = n_sites_per_mvg), vjust=0) + # how do I just have one label per mvg
  theme_bw()
