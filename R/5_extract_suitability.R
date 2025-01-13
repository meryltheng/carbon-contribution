library(sf)
library(terra)
library(geodata)
nvis<-rast('Data/Spatial/NVIS_MVG_RAS_AA.tif')

bioclim<-worldclim_country('bio', country="AUS",res=0.5, path='Data/Spatial/')
bioclim<-project(bioclim, crs(nvis))
bioclim<-crop(bioclim,nvis)
saveRDS(bioclim, 'Data/Spatial/bioclim.rds')


tmax<-worldclim_country('tmax',country="AUS",res=0.5, path='Data/Spatial/')
tmax<-project(tmax, crs(nvis))
tmax<-crop(tmax,nvis)
saveRDS(tmax, 'Data/Spatial/tmax.rds')

tmin<-worldclim_country('tmin',country="AUS",res=0.5, path='Data/Spatial/')
tmin<-project(tmin, crs(nvis))
tmin<-crop(tmin,nvis)
saveRDS(tmin, 'Data/Spatial/tmin.rds')

prec<-worldclim_country('prec',country="AUS",res=0.5, path='Data/Spatial/')
prec<-project(prec, crs(nvis))
prec<-crop(prec,nvis)
saveRDS(prec, 'Data/Spatial/prec.rds')

vapr<-worldclim_global(var="vapr",res=0.5, path='Data/Spatial/')
vapr<-crop(vapr,project(nvis, crs(vapr)))
vapr<-project(vapr, crs(nvis))
saveRDS(vapr, 'Data/Spatial/vapr.rds')


### Get 3-month wettest period for each cell

