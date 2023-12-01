library(data.table) # for fast reads of large datasets
library(raster)
library(sf)
library(dplyr)

"----------- Read NVIS data (MVG and MVS) ---------------"
# read NVIS raster layers (already projected to Australian Albers)
# read MVG layer
raster_file = 'Data/Spatial/NVIS_MVG_RAS_AA.tif'
MVG <- raster(raster_file)

# read MVS layer
raster_file = 'Data/Spatial/NVIS_MVS_RAS_AA.tif'
MVS <- raster(raster_file)

# Define commonly required projections (for crs transformation later)
LatLong <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
Albers  <- '+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'

"----------- Biomass Plot Library data ---------------"
# https://portal.tern.org.au/metadata/TERN/fc4a7249-ebb2-4ada-8e06-b552bfb297a3

# read and subset to 'recent' data (2013 onwards) 
file_name = 'Data/Biomass Plot Library/biolib_treelist.csv' # tree-level data
all_data <- fread(file_name)
bpl_data <- all_data[obs_time >= "2013-01-01"]
length(unique(bpl_data$tree)) # all records of unique trees (might wanna double-check)
rm(all_data); gc()

# add genus (other cols are non-impt placeholders)
bpl_data[, c("genus", "sp", "spp", "blah") := tstrsplit(species, " ", fixed=TRUE)]
bpl_data$family <- NA

# match record locations to MVG/MVS class
bpl_data_sf <- sf::st_as_sf(bpl_data, coords = c("longitude", "latitude"), crs = LatLong) # GDA94
# reproject to Australian Albers
bpl_data_AA <- sf::st_transform(bpl_data_sf,
                                 crs = Albers)

# append MVG to site info
mvgValue = extract(MVG, sf::st_coordinates(bpl_data_AA))
bpl_data_AA <- cbind(bpl_data_AA, mvgValue)

# append MVG to site info
mvsValue = extract(MVS, sf::st_coordinates(bpl_data_AA))
bpl_data_AA <- cbind(bpl_data_AA, mvsValue)

# reformat to standardise
bpl_data_input <- bpl_data_AA %>%
  select(geometry = geometry,
         #longitude = longitude,
         #latitude = latitude,
         site = site,
         sitearea_ha = sitearea_ha, # total area of site (sum of all subplots)
         obs_time = obs_time,
         plantID = tree,
         family = family,
         genus = genus,
         species = species,
         measurement = measurement, # method used to measure diameter
         diameter = diameter,
         height = ht, # tree height; not alws avail
         agb_est = agb_drymass, # estimated; we won't use this, will use our own calculations for standardisation (note: BPL estimates look to be in the wrong scale not in kg)
         mvgValue = mvgValue,
         mvsValue = mvsValue
         )

# save(bpl_data_AA, file = "bpl_data_AA.RData")

"----------- Robson Creek Rainforest data ---------------"
# https://portal.tern.org.au/metadata/TERN/5a864ec1-f780-4104-a242-7ff00e1f4962

# read and subset to 'recent' data (2013 onwards) 
file_name = 'Data/Robson Creek Rainforest/Robson_Creek_diameter_height_biomass_data.csv' # tree-level data
all_data <- fread(file_name)
robcreek_data <- all_data %>%
  mutate(obs_time = as.Date(phenomenonTime, "%d/%m/%Y")) %>%
  filter(obs_time >= "2013-01-01") %>%
  arrange(obs_time) %>% 
  group_by(plantId) %>%  # one plant has multuple observations
  slice(n()) # last row (most recent observation of the plant)
rm(all_data); gc()

# match record locations to MVG/MVS class
robcreek_data_sf <- sf::st_as_sf(robcreek_data, coords = c("longitude", "latitude"), crs = LatLong) # GDA94
# reproject to Australian Albers
robcreek_data_AA <- sf::st_transform(robcreek_data_sf,
                                crs = Albers)

# append MVG to site info
mvgValue = extract(MVG, sf::st_coordinates(robcreek_data_AA))
robcreek_data_AA <- cbind(robcreek_data_AA, mvgValue)

# append MVG to site info
mvsValue = extract(MVS, sf::st_coordinates(robcreek_data_AA))
robcreek_data_AA <- cbind(robcreek_data_AA, mvsValue)

# reformat to standardise
robcreek_data_input <- robcreek_data_AA %>%
  select(geometry = geometry,
         #longitude = longitude,
         #latitude = latitude,
         site = siteId,
         sitearea_ha = plotLength_metres, 
         obs_time = obs_time,
         plantID = plantId,
         family = family,
         genus = genus,
         species = scientificName,
         measurement = stemDiameterPointOfMeasurement_metres, # method used to measure diameter
         diameter = stemDiameter_centimetres,
         height = stemHeight_metres, # tree height; not alws avail
         agb_est = aboveGroundBiomass_kilograms, # estimated; we won't use this, will use our own calculations for standardisation
         mvgValue = mvgValue,
         mvsValue = mvsValue
  ) %>%
  mutate(site = rep("Robertson Creek", nrow(robcreek_data_AA)), 
         sitearea_ha = rep(1, nrow(robcreek_data_AA)), # 1ha
         measurement = rep("D130", nrow(robcreek_data_AA)), # diameter at 130cm
         )

"----------- Calperum Mallee data ---------------"
# https://portal.tern.org.au/metadata/TERN/c8d054b2-9dfc-4328-8251-f0fd9276083d

# read and subset to 'recent' data (2013 onwards) 
file_name = 'Data/Calperum Mallee/Calperum_Mallee_diameter_height_biomass_data.csv' # tree-level data
all_data <- fread(file_name)
calperum_data <- all_data %>%
  mutate(obs_time = as.Date(phenomenonTime, "%d/%m/%Y")) %>%
  filter(obs_time >= "2013-01-01") %>%
  arrange(obs_time) %>% 
  group_by(plantId) %>%  # one plant has multuple observations
  slice(n()) # last row (most recent observation of the plant); lol 85 plants
rm(all_data); gc()

# match record locations to MVG/MVS class
calperum_data_sf <- sf::st_as_sf(calperum_data, coords = c("longitude", "latitude"), crs = LatLong) # GDA94
# reproject to Australian Albers
calperum_data_AA <- sf::st_transform(calperum_data_sf,
                                     crs = Albers)

# append MVG to site info
mvgValue = extract(MVG, sf::st_coordinates(calperum_data_AA))
calperum_data_AA <- cbind(calperum_data_AA, mvgValue)

# append MVG to site info
mvsValue = extract(MVS, sf::st_coordinates(calperum_data_AA))
calperum_data_AA <- cbind(calperum_data_AA, mvsValue)

# reformat to standardise
calperum_data_input <- calperum_data_AA %>%
  select(geometry = geometry,
         #longitude = longitude,
         #latitude = latitude,
         site = siteId,
         sitearea_ha = plotLength_metres, 
         obs_time = obs_time,
         plantID = plantId,
         family = family,
         genus = genus,
         species = scientificName,
         measurement = stemDiameterPointOfMeasurement_metres, # method used to measure diameter
         diameter = stemDiameter_centimetres,
         height = stemHeight_metres, # tree height; not alws avail
         agb_est = aboveGroundBiomass_kilograms, # estimated; we won't use this, will use our own calculations for standardisation
         mvgValue = mvgValue,
         mvsValue = mvsValue
  ) %>%
  mutate(site = rep("Calperum Mallee", nrow(calperum_data_AA)), 
         sitearea_ha = rep(1, nrow(calperum_data_AA)), # 1ha
         measurement = rep("D10", nrow(calperum_data_AA)), # diameter at 10cm
  )

"----------- Compile and save all data ---------------"
save(bpl_data_input,
     robcreek_data_input,
     calperum_data_input,
     file = "Input/data_clean.RData")