library(raster)

P4S.latlon <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

#### 1. Define extent and import study area shapefile ####
# Import shapefile
# https://download-directory.github.io/
countries <- rgdal::readOGR("D:/Github/BGE-SDM/GISDATA/Study Area SHP")
#dev.off()
proj4string(countries) <- P4S.latlon
plot(countries)
class(countries)
str(countries); countries@bbox

#width= Unit is meter if x has a longitude/latitude CRS
countriesBuffer <- buffer(countries, width=500000, dissolve=TRUE) # 5 Degree buffer
plot(countriesBuffer, add=T, col='red')
plot(countries, add=T)
class(countriesBuffer)
extent(countries);extent(countriesBuffer)
#xmin: extent(countriesBuffer)[1]

#Define extent
extent <- extent(-43, 82, 25, 81) ##NEED TO BE REARRANGE FOR EAST SIDE

# Download precipitation/Bio/+++ data from WorldClim
global.clim <- raster::getData("worldclim", var="bio", res=5, download=T, path="RDATA")

files.present.bil <- list.files('D:/Github/BGE-SDM/RDATA/wc5/', pattern="[.]bil$", full.names=T) # alternatives for pattern (c|C)(e|E)(l|L)$
files.present.bil

# Loop for cropping with extent and write as ascii
pb <- txtProgressBar(min = 0, max = length(files.present.bil), style = 3)

for(i in files.present.bil)  {
  raster <- raster(i)
  raster <- crop(raster, extent)
  writeRaster(raster,
              filename  = paste("D:/Github/BGE-SDM/RDATA/clipped/", "fe_buffer_", basename(i), sep = ""),
              format    = 'ascii',
              NAflag    = -9999,
              overwrite = T)
  setTxtProgressBar(pb, i)
}
close(pb)

# read worldclim
files.present <- list.files('D:/Github/BGE-SDM/RDATA/clipped/', pattern="[.]asc$", full.names=T) # alternatives for pattern (c|C)(e|E)(l|L)$
files.present[1:19]
