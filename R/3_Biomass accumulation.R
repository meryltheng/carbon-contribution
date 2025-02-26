library(data.table) # for fast reads of large datasets
library(dplyr)
library(tidyverse)
library(sf)

"----------- Load processed data for Biomass ---------------"
load("Input/biomass_estimates.RData") # if processed

"----------- Assign growth rates based on species, then genus ---------------"
# Load reference (https://doi.org/10.25901/5b95acd389f07)
file_name = "Data/Tree growth/Ngugi2014-Growth rates of Eucalyptus and other Australian native tree species.csv"
growth_data = fread(file_name)

growth_data_spp <- growth_data %>%
  group_by(species) %>%
  filter(row_number()==1)  %>%
  rename(diameter_increment_yr = Overall_mu)

# a. First, assign based on species match 
biomass_growth <- left_join(biomass_estimates, growth_data_spp[,c("species", "diameter_increment_yr")], by = 'species') 

# b. If no match, assign based on genus match 
# subset the no match records
no_match <- biomass_growth %>%
  filter(is.na(diameter_increment_yr)) %>%
  dplyr::select(-c(diameter_increment_yr))

# subset the matches for compilation later
biomass_growth_1 <- biomass_growth %>%
  filter(!is.na(diameter_increment_yr)) 

# create genus-functional type reference
growth_data_genus <- growth_data_spp %>%
  group_by(genus) %>%
  summarise(diameter_increment_yr = mean(diameter_increment_yr)) # apply the genus mean growth rate on unassigned species records


# match by genus 
no_match <- left_join(no_match, growth_data_genus, by ='genus')

# c. If STILL no match...
still_no_match <- no_match[which(is.na(no_match$diameter_increment_yr)),] # subset the no match records
biomass_growth_2 <- no_match[which(!is.na(no_match$diameter_increment_yr)),] # subset the matches for compilation later

# look for growth data for missing genera
still_no_match <- still_no_match %>%
  mutate(diameter_increment_yr = if_else(genus == "Atherosperma" & measurement %in% c("DBH", "D130"), 0.1, diameter_increment_yr, NA)) %>% # https://www.jstor.org/stable/2679956 (based on Fig 7a midpoint estimate)
  mutate(diameter_increment_yr = if_else(genus == "Phyllocladus" & measurement %in% c("DBH", "D130"), 0.12, diameter_increment_yr, NA)) %>% # https://www.anbg.gov.au/gardens/visiting/exploring/walks/conifers/phyllocladus-asplen.html
  mutate(diameter_increment_yr = if_else(species == "Argyrodendron peralatum" & measurement %in% c("DBH", "D130"), 0.225, diameter_increment_yr, NA)) %>% # Argyrodendron peralatum https://doi.org/10.1016/0378-1127(89)90110-2
  mutate(diameter_increment_yr = if_else(species == "Argyrodendron trifoliolatum" & measurement %in% c("DBH", "D130"), 0.135, diameter_increment_yr, NA)) # Argyrodendron trifoliolatum/polyandrum

# can't find for: Acmena smithii, Austrobuxus swainii, Archirhodomyrtus beckleri, Nothofagus cunninghamii, Nematolepis squamea, Olearia argophylla, etc.

# d. If STILLLL no match...
stillll_no_match <- still_no_match[which(is.na(still_no_match$diameter_increment_yr)),] # subset the no match records
biomass_growth_3 <- still_no_match[which(!is.na(still_no_match$diameter_increment_yr)),] # subset the matches for compilation later

# assume conservative growth rates for all species (with DBH, D130 measurement) without growth rate info
stillll_no_match <- stillll_no_match %>%
  mutate(diameter_increment_yr = if_else(measurement %in% c("DBH", "D130"), 0.1, diameter_increment_yr, NA)) # assume slow growth of 0.1 cm/y

# e. COMPILE ALL
biomass_growth_est <- rbind(biomass_growth_1,
                            biomass_growth_2,
                            biomass_growth_3)

"----------- Dealing with D10s ---------------"
# Can't find diameter growth rate based on D10 measurement
# First, find proportion of biomass contributed by all plants with D10 measurement
biomass_estimates %>% 
  group_by(measurement) %>%
  summarise(biomass = sum(t_biomass_estimate, na.rm=T)) # ~12% biomass contributed by plants with D10 measurement

# ***Discard plants measured by D10***
biomass_growth_est<- biomass_growth_est %>%
  filter(measurement %in% c("DBH", "D130"))

"----------- Estimate biomass (AGB and BGB) sequestered after 1 year of diameter growth ---------------" 
biomass_seq <- biomass_growth_est %>%
  mutate(diameter_1y = diameter + diameter_increment_yr) %>%
  mutate(agb_estimate_1y = case_when(
    func_type_agb == "F Shrub" ~ exp(-3.007 + 2.428 * log(diameter_1y)) * 1.128,
    func_type_agb == "F Multi" ~ exp(-2.757 + 2.474 * log(diameter_1y)) * 1.079,
    func_type_agb == "F Euc" ~ exp(-2.016 + 2.375 * log(diameter_1y)) * 1.067,
    func_type_agb == "F Other-H" ~ exp(-1.693 + 2.220 * log(diameter_1y)) * 1.044,
    func_type_agb == "F Other-L" ~ exp(-2.573 + 2.460 * log(diameter_1y)) * 1.018
  )) %>%
  mutate(bgb_estimate_1y = case_when(
    func_type_bgb == "FShrub&Ac" ~ exp(-3.553 + 2.185 * log(diameter_1y)) * 1.160,
    func_type_bgb == "FMallee" ~ exp(-2.946 + 2.302 * log(diameter_1y)) * 1.116,
    func_type_bgb == "FTree" ~ exp(-2.682 + 2.212 * log(diameter_1y)) * 1.096,
    func_type_bgb == "FRadiata" ~ exp(-3.740 + 2.299 * log(diameter_1y)) * 1.053
  )) %>%
  mutate(t_biomass_estimate_1y = agb_estimate_1y + bgb_estimate_1y) %>%
  mutate(biomass_sequestered = t_biomass_estimate_1y - t_biomass_estimate) 

"----------- Biomass sequestered per ha per year per genus per MVG ---------------" 
biomass_seq_std <- biomass_seq %>%
  group_by(site, genus) %>%
  summarise(geometry = unique(geometry),
            site = unique(site),
            sitearea_ha = unique(sitearea_ha),
            genus = unique(genus),
            mvgValue = unique(mvgValue),
            mvsValue = unique(mvsValue),
            t_biomass_seq = sum(biomass_sequestered, na.rm = T), # total biomass estimate
            t_biomass_seq_ha = sum(biomass_sequestered, na.rm = T)/unique(sitearea_ha)) # standardised to per ha

biomass_seq_mvg <- biomass_seq_std %>%
  group_by(mvgValue) %>%
  mutate(n_sites_per_mvg = length(unique(site))) %>%
  group_by(mvgValue, genus) %>%
  summarise(n_sites_per_mvg = unique(n_sites_per_mvg),
            genus = unique(genus),
            mvgValue = unique(mvgValue),
            t_biomass_seq_ha = mean(t_biomass_seq_ha, na.rm = T)) # standardised to per ha

# class genera with > 10% contribution under "other"
biomass_seq_mvg0 <- biomass_seq_mvg %>% 
  group_by(mvgValue) %>%
  mutate(seq_perc = t_biomass_seq_ha/sum(t_biomass_seq_ha),
         genus_new = if_else(seq_perc > 0.1, genus, "z_Other", NA)) %>%
  group_by(mvgValue, genus_new) %>%
  summarise(n_sites_per_mvg = unique(n_sites_per_mvg),
            genus = unique(genus_new),
            mvgValue = unique(mvgValue),
            t_biomass_seq_ha = sum(t_biomass_seq_ha, na.rm = T)) 

#save(biomass_seq, biomass_seq_std, biomass_seq_mvg, biomass_seq_mvg0, file = "Input/biomass_seq.RData")

"----------- Visualise! ---------------"
library(ggplot2)
#load("Input/biomass_seq.RData") # if saved

# plotting variables
x.labels = c("1 Rainforests and Vine Thickets", #1
             "2 Eucalypt Tall Open Forests", #2
             "3 Eucalypt Open Forests", #3
             "4 Eucalypt Low Open Forests", #4
             "5 Eucalypt Woodlands", #5 
             "6 Acacia Forests and Woodlands", #6
             "7 Callitris Forests and Woodlands", #7
             "8 Casuarina Forests and Woodlands", #8
             "9 Melaleuca Forests and Woodlands", #9
             "10 Other Forests and Woodlands", #10
             "11 Eucalypt Open Woodlands", #11
             "12 Tropical Eucalypt Woodlands/Grasslands", #12
             "13 Acacia Open Woodlands", #13
             "14 Mallee Woodlands and Shrublands", #14
             "15 Low Closed Forests and Tall Closed Shrublands", #15
             "16 Acacia Shrublands", #16
             "17 Other Shrublands", #17
             "18 Heathlands", #18
             "19 Tussock Grasslands", #19
             "20 Hummock Grasslands", #20
             "21 Other Grasslands, Herblands, Sedgelands and Rushlands", #21
             "22 Chenopod Shrublands, Samphire Shrublands and Forblands", #22
             "23 Mangroves", #23
             "24 Inland Aquatic", #24
             "25 Cleared, non-native veg., buildings", #25
             "26 Unclassfied native veg.", #26
             "27 Naturally bare - sand, rock, claypan, mudflat", #27
             "28 Sea and estuaries", #28
             "29 Regrowth, modified native vegetation", #29
             "30 Unclassified forest", #30
             "31 Other Open Woodlands", #31
             "32 Mallee Open Woodlands and Sparse Mallee Shrublands") #32

# PLOT
C_seq_plot <- ggplot(biomass_seq_mvg0, aes(x = as.factor(mvgValue), y = t_biomass_seq_ha/1000*0.5, fill = genus)) +
  geom_bar(stat = 'identity', color='black') + # position='fill' if want proportion
  labs(x = "MVG", y = "Carbon sequestered per ha per year (t)") +
  scale_fill_manual(values = c("dodgerblue2", "#E31A1C", "gold2", "darkorange4", "green4", "orchid1", "black"),
                    name = "Genus") +
  #geom_text(aes(label = n_sites_per_mvg), vjust=0) + # how do I just have one label per mvg
  scale_x_discrete(labels = str_wrap(x.labels[unique(biomass_seq_mvg0$mvgValue)], 7)) +
  theme_bw() +
  theme(axis.text.x=element_text(size=8),
        axis.title.y = element_text(size = 12))


ggsave(filename=paste0("Plots/C_seq_ha_mvg.pdf"), plot=C_seq_plot, width=9, height=5) 


