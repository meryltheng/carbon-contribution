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
    func_type_bgb == "FShrub&Ac" ~ exp(-3.553 + 2.185 * log(diameter)) * 1.160,
    func_type_bgb == "FMallee" ~ exp(-2.946 + 2.302 * log(diameter)) * 1.116,
    func_type_bgb == "FTree" ~ exp(-2.682 + 2.212 * log(diameter)) * 1.096,
    func_type_bgb == "FRadiata" ~ exp(-3.740 + 2.299 * log(diameter)) * 1.053
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
            t_biomass_per_ha = mean(t_biomass_per_ha, na.rm = T))# %>% # standardised to per ha
  #group_by(mvgValue) %>%
 # summarise(biomass_perc = )

# class genera with > 10% contribution under "other"
biomass_est_mvg0 <- biomass_est_mvg %>% 
  group_by(mvgValue) %>%
  mutate(biomass_perc = t_biomass_per_ha/sum(t_biomass_per_ha),
         genus_new = if_else(biomass_perc > 0.1, genus, "z_Other", NA)) %>%
  group_by(mvgValue, genus_new) %>%
  summarise(n_sites_per_mvg = unique(n_sites_per_mvg),
            genus = unique(genus_new),
            mvgValue = unique(mvgValue),
            t_biomass_per_ha = sum(t_biomass_per_ha, na.rm = T)) 

#save(biomass_estimates, biomass_est_std, biomass_est_mvg, file = "Input/biomass_estimates.RData")

"----------- Visualise! ---------------"
library(ggplot2)
load("Input/biomass_estimates.RData") # if saved

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
# ggplot(biomass_est_mvg, aes(x = as.factor(mvgValue), y = t_biomass_per_ha/1000, fill = genus)) +
#   geom_bar(stat = 'identity',color='black') + # position='fill' if want proportion
#   labs(x = "MVG", y = "Biomass contribution per ha (t)") +
#   #geom_text(aes(label = n_sites_per_mvg), vjust=0) + # how do I just have one label per mvg
#   theme_bw()

# with genera > 10% contribution under "other"
ggplot(biomass_est_mvg0, aes(x = as.factor(mvgValue), y = t_biomass_per_ha/1000, fill = genus)) +
  geom_bar(stat = 'identity',color='black') + # position='fill' if want proportion
  labs(x = "MVG", y = "Biomass contribution per ha (t)") +
  scale_fill_manual(values = c("dodgerblue2", "#E31A1C", "gold2", "green4", "#FF7F00", "black"),
                    name = "Genus") +
  #geom_text(aes(label = n_sites_per_mvg), vjust=0) + # how do I just have one label per mvg
  scale_x_discrete(labels = str_wrap(x.labels[unique(biomass_est_mvg0$mvgValue)], 7)) +
  theme(axis.text.x=element_text(size=7)) +
  theme_bw() 

# flipped plot
ggplot(biomass_est_mvg0, aes(x = t_biomass_per_ha/1000, y = as.factor(mvgValue), fill = genus)) +
  geom_bar(stat = 'identity',color='black') + # position='fill' if want proportion
  labs(y = "MVG", x = "Biomass contribution per ha (t)") +
  scale_fill_manual(values = c("dodgerblue2", "#E31A1C", "gold2", "green4", "#FF7F00", "black"),
                    name = "Genus") +
  #geom_text(aes(label = n_sites_per_mvg), vjust=0) + # how do I just have one label per mvg
  scale_y_discrete(labels = str_wrap(x.labels[unique(biomass_est_mvg0$mvgValue)], 25)) +
  theme(axis.text.y=element_text(size=6)) +
  theme_bw() 
