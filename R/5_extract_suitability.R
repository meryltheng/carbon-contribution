library(sf)
library(terra)
library(geodata)
library(dismo)
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


## Stolen from BIOVARS package
window <- function(x)  { 
  lng <- length(x)
  x <- c(x,  x[1:3])
  m <- matrix(ncol=3, nrow=lng)
  for (i in 1:3) { m[,i] <- x[i:(lng+i-1)] }
  return(apply(m, MARGIN=1, FUN=sum))
}
window_min <- function(x)  { 
  lng <- length(x)
  x <- c(x,  x[1:3])
  m <- matrix(ncol=3, nrow=lng)
  for (i in 1:3) { m[,i] <- x[i:(lng+i-1)] }
  return(apply(m, MARGIN=1, FUN=min))
}
window_max <- function(x)  { 
  lng <- length(x)
  x <- c(x,  x[1:3])
  m <- matrix(ncol=3, nrow=lng)
  for (i in 1:3) { m[,i] <- x[i:(lng+i-1)] }
  return(apply(m, MARGIN=1, FUN=max))
}
prec<-as.matrix(prec)
# precip by quarter (3 months)		
wet <- t(apply(prec, 1, window))
#wettest quarter
wetqrt <- cbind(1:nrow(prec), as.integer(apply(wet, 1, which.max)))
#  Min Temperature of Wettest Quarter 
tmpmin<-as.matrix(tmin)
tmpmin <- t(apply(tmpmin, 1, window))/3 
tminwet<- tmpmin[wetqrt]
bioclim$tminwet<-tminminwet
#  Max Temperature of Wettest Quarter 
tmpmax<-as.matrix(tmax)
tmpmax <- t(apply(tmpmax, 1, window))/3
tmaxwet<- tmpmax[wetqrt]
bioclim$tmaxwet<-tmaxwet

# Min  Min Temperature of Wettest Quarter 
tmpmin2<-as.matrix(tmin)
tmpminmin <- t(apply(tmpmin2, 1, window_min)) 
tminminwet<- tmpminmin[wetqrt]
bioclim$tminminwet<-tminminwet
#Max Max Temperature of Wettest Quarter 
tmpmax2<-as.matrix(tmax)
tmpmaxmax <- t(apply(tmpmax2, 1, window_max))
tmaxmaxwet<- tmpmaxmax[wetqrt]
bioclim$tmaxmaxwet<-tmaxmaxwet

writeRaster(bioclim, 'Data/Spatial/bioclim_newvars.tif', overwrite=TRUE)

##subset solar radiation to correct quarter

solar<-Rast('Data/Spatial/srad.tif')
solar$cell<-over