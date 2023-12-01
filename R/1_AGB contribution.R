library(data.table) # for fast reads of large datasets
library(dplyr)

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

# Load reference (https://portal.tern.org.au/metadata/TERN/aef6e89d-6fde-42dc-b77b-fc512acb72de)
file_name = "Data/Australian Tree Biomass Library/AGB_Allometric_database-_v2.csv"
allo_data = fread(file_name)

# a. First, assign based on species match 
# create species-functional type reference
allo_ref_sp <- allo_data %>%
  mutate(species = paste(Genus, Species)) %>%
  group_by(species) %>%
  slice(n()) %>%
  select(genus = Genus, species = species, 
         func_type = `Plant Funcational Type`, func_type_sub = `Plant Funcational Type Sub-category`)

# match by species 
all_data_allo <- left_join(all_data, allo_ref_sp, by ='species') %>%
  mutate(func_type = if_else(genus.x == "Eucalyptus" & measurement =="D10", "F Multi", func_type, NA)) %>%
  rename("genus" = "genus.x")

# b. If no match, assign based on genus match 
# subset the no match records
no_match <- all_data_allo %>%
  filter(is.na(func_type)) %>%
  select(-c(func_type, func_type_sub, genus.y))

# subset the matches for compilation later
allo_data_1 <- all_data_allo %>%
  filter(!is.na(func_type)) %>%
  select(-c(genus.y))

# create genus-functional type reference
allo_ref_gen <- allo_ref_sp %>%
  group_by(genus) %>%
  # only create genus ref if func_type for all species within genus is same
  mutate(unique = if_else(length(unique(func_type))==1, 1, 0, NA)) %>%
  filter(unique == 1) %>%
  slice(n()) %>%
  select(genus = genus, 
         func_type = func_type, func_type_sub = func_type_sub)

# match by genus 
no_match <- left_join(no_match, allo_ref_gen, by ='genus')

# c. If STILL no match, we fix the common genera (Eucalypts and Acacia) first
still_no_match <- no_match[which(is.na(no_match$func_type)),] # subset the no match records
allo_data_2 <- no_match[which(!is.na(no_match$func_type)),] # subset the matches for compilation later

# assign common genera (Eucalypts and Acacia) to known functional types
still_no_match <- still_no_match %>%
  mutate(func_type = if_else(genus == "Eucalyptus" & measurement %in% c("DBH", "D130"), "F Euc", func_type, NA)) %>%
  mutate(func_type = if_else(genus == "Eucalyptus" & measurement %in% c("DBH", "D130"), "F Euc", func_type, NA)) %>%
  mutate(func_type = if_else(genus == "Acacia" & measurement == "D10", "F Multi", func_type, NA)) 

# d. remaining no matches (many rainforest species, I think) will be assigned into some group based on knowledge and assumptions (REVISIT at some point)
stillll_no_match <- still_no_match[which(is.na(still_no_match$func_type)),] # subset the no match records
allo_data_3 <- still_no_match[which(!is.na(still_no_match$func_type)),] # subset the matches for compilation later

# assign the rest based on info and assumptions
allo_data_4 <- stillll_no_match %>%
  mutate(func_type = if_else(genus == "Acacia" & measurement %in% c("DBH", "D130"), "F Euc", func_type, NA)) %>%
  mutate(func_type = if_else(measurement == "D10", "F Shrub", func_type, NA)) %>% # assume all remaining D10s are shrubs
  mutate(func_type = if_else(genus %in% c("Pinus", "Araucaria", "Agathis") & measurement %in% c("DBH", "D130"), "F Other-L", func_type, NA)) %>% # low density conifers
  mutate(func_type = if_else(is.na(func_type) & measurement %in% c("DBH", "D130"), "F Euc", func_type, NA))

# e. COMPILE ALL
allo_data <- rbind(allo_data_1,
                   allo_data_2,
                   allo_data_3,
                   allo_data_4)

rm(all_data, allo_data_1, allo_data_2, allo_data_3, allo_data_4, no_match, still_no_match, stillll_no_match, bpl_data_input, calperum_data_input, robcreek_data_input)

"----------- Estimate AGB based on allometric eqns per functional type ---------------" 
# estimate AGB based on eqns given by https://doi.org/10.1111/gcb.13201 
allo_data <- allo_data %>%
  mutate(agb_estim = case_when(
    func_type == "F Shrub" ~ exp(-3.007 + 2.428 * log(diameter) * 1.128),
    func_type == "F Multi" ~ exp(-2.757 + 2.474 * log(diameter) * 1.079),
    func_type == "F Euc" ~ exp(-2.016 + 2.375 * log(diameter) * 1.067),
    func_type == "F Other-H" ~ exp(-1.693 + 2.220 * log(diameter) * 1.044),
    func_type == "F Other-L" ~ exp(-2.573 + 2.460 * log(diameter) * 1.018)
  ))

save(allo_data, file = "Input/allometric_data.RData")

"----------- AGB per ha per genus per MVG ---------------" 
load("Input/allometric_data.RData") # if processed

allo_data_std <- allo_data %>%
  group_by(site, genus) %>%
  summarise(geometry = unique(geometry),
            site = unique(site),
            sitearea_ha = unique(sitearea_ha),
            genus = unique(genus),
            mvgValue = unique(mvgValue),
            mvsValue = unique(mvsValue),
            t_agb_est = sum(agb_estim, na.rm = T), # total agb estimate
            t_agb_per_ha = sum(agb_estim, na.rm = T)/unique(sitearea_ha)) # standardised to per ha

allo_data_mvg <- allo_data_std %>%
  group_by(mvgValue) %>%
  mutate(n_sites_per_mvg = length(unique(site))) %>%
  group_by(mvgValue, genus) %>%
  summarise(n_sites_per_mvg = unique(n_sites_per_mvg),
            genus = unique(genus),
            mvgValue = unique(mvgValue),
            t_agb_per_ha = mean(t_agb_per_ha, na.rm = T)) # standardised to per ha

# check if n_sites per MVG is correct
# allo_data_std %>% group_by(mvgValue) %>%
#   summarise(MVG = unique(mvgValue), n_sites = length(unique(site)))

"----------- Visualise! ---------------"
library(ggplot2)

ggplot(allo_data_mvg, aes(x = as.factor(mvgValue), y = t_agb_per_ha, fill = genus)) +
  geom_bar(stat = 'identity',color='black') + # position='fill' if want proportion
  labs(x = "MVG", y = "AGB contribution per ha (kg)") +
  #geom_text(aes(label = n_sites_per_mvg), vjust=0) + # how do I just have one label per mvg
  theme_bw()
  