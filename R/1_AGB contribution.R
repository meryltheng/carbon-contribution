library(data.table) # for fast reads of large datasets
library(dplyr)
library(sf)

"----------- Load processed data and compile ---------------"
load("Input/data_clean.RData")

# merge standardised datasets
all_data <- rbind(bpl_data_input,
                  calperum_data_input,
                  robcreek_data_input)

# remove records without genus/species information
all_data <- all_data %>%
  filter(!is.na(genus), species != "Tree") 

"----------- Assign record to functional plant type ---------------"
# Quick and dirty way to assign each record to a functional plant type for AGB calculations later

# Load reference (https://doi.org/10.25901/5b95acd389f07)
file_name = "Data/Australian Tree Biomass Library/AGB_allometrics_database-v4.csv"
agb_data = fread(file_name)

# a. First, assign based on species match 
# create species-functional type reference
agb_ref_sp <- agb_data %>%
  mutate(species = paste(Genus, Species)) %>%
  group_by(species) %>%
  slice(n()) %>%
  select(genus = Genus, species = species, 
         func_type_agb = `Plant Funcational Type`) # func_type_sub = `Plant Funcational Type Sub-category`

# match by species 
all_data_agb <- left_join(all_data, agb_ref_sp, by ='species') %>%
  mutate(func_type_agb = if_else(genus.x == "Eucalyptus" & measurement =="D10", "F Multi", func_type_agb, NA)) %>%
  rename("genus" = "genus.x")

# b. If no match, assign based on genus match 
# subset the no match records
no_match <- all_data_agb %>%
  filter(is.na(func_type_agb)) %>%
  select(-c(func_type_agb, genus.y))

# subset the matches for compilation later
agb_data_1 <- all_data_agb %>%
  filter(!is.na(func_type_agb)) %>%
  select(-c(genus.y))

# create genus-functional type reference
agb_ref_gen <- agb_ref_sp %>%
  group_by(genus) %>%
  # only create genus ref if func_type for all species within genus is same
  mutate(unique = if_else(length(unique(func_type_agb))==1, 1, 0, NA)) %>%
  filter(unique == 1) %>%
  slice(n()) %>%
  select(genus = genus, 
         func_type_agb = func_type_agb) # func_type_sub = func_type_sub

# match by genus 
no_match <- left_join(no_match, agb_ref_gen, by ='genus')

# c. If STILL no match, we fix the common genera (Eucalypts and Acacia) first
still_no_match <- no_match[which(is.na(no_match$func_type_agb)),] # subset the no match records
agb_data_2 <- no_match[which(!is.na(no_match$func_type_agb)),] # subset the matches for compilation later

# assign common genera (Eucalypts and Acacia) to known functional types
still_no_match <- still_no_match %>%
  mutate(func_type_agb = if_else(genus == "Eucalyptus" & measurement %in% c("DBH", "D130"), "F Euc", func_type_agb, NA)) %>%
  mutate(func_type_agb = if_else(genus %in% c("Acacia", "Eucalyptus") & measurement == "D10", "F Multi", func_type_agb, NA)) # assumption that all remaining Acacias are trees not shrubs

# d. remaining no matches (many rainforest species, I think) will be assigned into some group based on knowledge and assumptions (REVISIT at some point)
stillll_no_match <- still_no_match[which(is.na(still_no_match$func_type)),] # subset the no match records
agb_data_3 <- still_no_match[which(!is.na(still_no_match$func_type)),] # subset the matches for compilation later

# assign the rest based on info and assumptions
agb_data_4 <- stillll_no_match %>%
  mutate(func_type_agb = if_else(genus == "Acacia" & measurement %in% c("DBH", "D130"), "F Euc", func_type_agb, NA)) %>% # Acacia wood density greatly variable; I assume here it is similar to Euc
  mutate(func_type_agb = if_else(measurement == "D10", "F Shrub", func_type_agb, NA)) %>% # assume all remaining D10s are shrubs
  mutate(func_type_agb = if_else(genus %in% c("Pinus", "Araucaria", "Agathis") & measurement %in% c("DBH", "D130"), "F Other-L", func_type_agb, NA)) %>% # low density conifers
  mutate(func_type_agb = if_else(is.na(func_type_agb) & measurement %in% c("DBH", "D130"), "F Euc", func_type_agb, NA)) # assume all remaining D130s are similar to Eucs

# e. COMPILE ALL
agb_estimates <- rbind(agb_data_1,
                   agb_data_2,
                   agb_data_3,
                   agb_data_4)

rm(all_data, agb_data_1, agb_data_2, agb_data_3, agb_data_4, no_match, still_no_match, stillll_no_match, bpl_data_input, calperum_data_input, robcreek_data_input)

"----------- Estimate AGB based on allometric eqns per functional type ---------------" 
# estimate AGB based on eqns given by https://doi.org/10.1111/gcb.13201 
agb_estimates <- agb_estimates %>%
  mutate(agb_estimate = case_when(
    func_type_agb == "F Shrub" ~ exp(-3.007 + 2.428 * log(diameter) * 1.128),
    func_type_agb == "F Multi" ~ exp(-2.757 + 2.474 * log(diameter) * 1.079),
    func_type_agb == "F Euc" ~ exp(-2.016 + 2.375 * log(diameter) * 1.067),
    func_type_agb == "F Other-H" ~ exp(-1.693 + 2.220 * log(diameter) * 1.044),
    func_type_agb == "F Other-L" ~ exp(-2.573 + 2.460 * log(diameter) * 1.018)
  ))

save(agb_estimates, file = "Input/agb_estimates.RData")

"----------- AGB per ha per genus per MVG ---------------" 
load("Input/agb_estimates.RData") # if processed

agb_est_std <- agb_estimates %>%
  group_by(site, genus) %>%
  summarise(geometry = unique(geometry),
            site = unique(site),
            sitearea_ha = unique(sitearea_ha),
            genus = unique(genus),
            mvgValue = unique(mvgValue),
            mvsValue = unique(mvsValue),
            t_agb_est = sum(agb_estimate, na.rm = T), # total agb estimate
            t_agb_per_ha = sum(agb_estimate, na.rm = T)/unique(sitearea_ha)) # standardised to per ha

agb_est_mvg <- agb_est_std %>%
  group_by(mvgValue) %>%
  mutate(n_sites_per_mvg = length(unique(site))) %>%
  group_by(mvgValue, genus) %>%
  summarise(n_sites_per_mvg = unique(n_sites_per_mvg),
            genus = unique(genus),
            mvgValue = unique(mvgValue),
            t_agb_per_ha = mean(t_agb_per_ha, na.rm = T)) # standardised to per ha

# check if n_sites per MVG is correct
# agb_data_std %>% group_by(mvgValue) %>%
#   summarise(MVG = unique(mvgValue), n_sites = length(unique(site)))

"----------- Visualise! ---------------"
library(ggplot2)

ggplot(agb_est_mvg, aes(x = as.factor(mvgValue), y = t_agb_per_ha, fill = genus)) +
  geom_bar(stat = 'identity',color='black') + # position='fill' if want proportion
  labs(x = "MVG", y = "AGB contribution per ha (kg)") +
  #geom_text(aes(label = n_sites_per_mvg), vjust=0) + # how do I just have one label per mvg
  theme_bw()
  