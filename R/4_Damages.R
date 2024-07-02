library(data.table) # for fast reads of large datasets
library(dplyr)
library(tidyverse)
library(sf)
library(raster)

"----------- Load processed data for Biomass ---------------"
load("Input/biomass_seq.RData") # if processed

"----------- Create ref of MR impacts at genus-level ---------------"
# Load Myrtle Rust impact info
file_name = "Data/MR_impact.csv"
impact_data = read.csv(file_name)

# Average impact % across all records for the genus
impact_data_ref <- impact_data %>%
  group_by(Genus) %>%
  summarise(Impact = mean(Damage, na.rm = T)) %>%
  filter(!is.na(Impact)) %>%
  rename(genus = Genus)

"----------- Per ha reduction proportion (%) according to MVG  ---------------"
impact_mvg <- left_join(biomass_seq_mvg, impact_data_ref, by = 'genus') %>% # Assign based on genus match 
  mutate(Impact = if_else(is.na(Impact), 0, Impact, NA)) %>% # if no info, impact = 0
  mutate(Impact_biomass = t_biomass_seq_ha * Impact) %>% # calc absolute impact to C-seq
  group_by(mvgValue) %>%
  summarise(Impact_prop = sum(Impact_biomass) / sum(t_biomass_seq_ha) ) %>% # impact proportion by MVG
  st_drop_geometry() %>%
  rename(MVG = mvgValue) 

"----------- Calculate area of each MVG  ---------------"
# read file
raster_file = 'Data/Spatial/NVIS_MVG_RAS_AA.tif'
MVG <- raster(raster_file)

# calculate in batches due to large file size (otherwise will hit mem limits)
n_cells = MVG@ncols * MVG@nrows
cells_per_batch = 1E7
n_batch <- ceiling(n_cells / (cells_per_batch))

area_store <- list()

# this takes a while...
for (i in 1:n_batch){
  if(i == n_batch){ # last batch 
    vals <- MVG[(cells_per_batch*(i-1)+1):n_cells]
  } else{
    vals <- MVG[(cells_per_batch*(i-1)+1):(cells_per_batch*i)]
  }
  area_store[[i]] <- sapply(1:32, function(x) length(subset(vals, vals == x)) * res(MVG)[1]^2); rm(vals); gc()
  cat("round = ", i, " out of ", n_batch, "\r"); flush.console()
}

# compile batch calculations to obtain area of each MVG
area_dat <- do.call(rbind, area_store)
MVG_area <- colSums(area_dat)/(100*100) # in ha

#save(MVG_area, file = "Input/MVG_area.RData") 

"----------- Calculate damage $$  ---------------"
load(file = "Input/MVG_area.RData") 

# Assign MVGs with no data 0 impact
impact_mvg <- impact_mvg %>%
  complete(MVG = 1:32, fill = list(Impact_prop = 0)) 

# load sequestration value data 
value_data = read.csv("Input/LUP_val_sequestration.csv")
MVG_value_2015 <- as.vector(t(value_data[1,1:32]))
#discount_rate <- 0.03 # discount rate of 3% for environmental assets (Dodd et al. 2020)
#MVG_value_2023 <- MVG_value_2015 * (1 + discount_rate)^(2023-2015) # update to 2023 values

# calculate damage $$
damage_dat <- data.frame(MVG = impact_mvg$MVG,
                         Impact_perc = impact_mvg$Impact_prop,
                         Value = MVG_value_2015,
                         Area = MVG_area) %>%
  mutate(Damage_value = Impact_perc * Value * Area,
         Total_value = Value * Area) 

( DAM_R = sum(damage_dat$Damage_value, na.rm = T) / sum(damage_dat$Total_value, na.rm = T) )

"----------- Prepare output table for report ---------------"
MVG_names = c("1 Rainforests and Vine Thickets", #1
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

output_table <- data.frame(MVG = MVG_names,
                           Impact_perc = round(impact_mvg$Impact_prop*100,1),
                           Value = round(MVG_value_2015,1),
                           Area = MVG_area) %>%
  mutate(Damage_value = round(Impact_perc/100 * Value * Area/1000000), # in millions
         Total_value = round(Value * Area/1000000)) # in millions

write.csv(output_table, file = "Output/Damage_estimates.csv", row.names = F)  