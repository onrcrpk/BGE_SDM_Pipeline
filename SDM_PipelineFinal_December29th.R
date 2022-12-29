rm(list = ls(all=T))

#### Install required libraries ####
library(rgbif)
library(sp)
library(raster) 
library(mapr)
library(dismo)
library(rgeos)
library(ggplot2)
library(countrycode)
library(CoordinateCleaner)
library(dplyr)
library(ggplot2)



### nullModel function ####
nullModel <- function (x, n = 10, rep = 99) # n = number of records
{
  e <- list()
  #stopifnot(n < nrow(x))
  for (r in 1:rep) {
    z <- sample(nrow(x), n) # draw random presences
    pres <- x[z, ] # df with randomly drawn presences
    absc <- x[-z, ] # df with remaining points set as absences
    #rs10k <- sample(1:nrow(absc), 10000) # random select 10k rows
    #absc <- absc[rs10k,]
    rs1k <- sample(1:nrow(absc), 1000) # random select 1k rows
    absc <- absc[rs1k,]
    d <- rbind(pres, absc)
    v <- c(rep(1, nrow(pres)), rep(0, nrow(absc)))
    m.null <- maxent(d,v, args = c("noproduct", "nothreshold", "nohinge"))
    
    e[[r]] <- evaluate(pres, absc, m.null)
    cat("-")
    if (r%%50 == 0) 
      cat(" ", r, "\n")
    flush.console()
  }
  if (r%%50 != 0) {
    cat(" ", r, "\n")
  }
  else {
    cat("\n")
  }
  e
}
#<environment: namespace:dismo>



### Start ####
#Create empty lists.

#Part 1
diptera.total <- 0
diptera.gap <- 0
no.key <- character()
missing.gbif <- character()
raw.data.more <- character()
raw.data.less <- character()
no.country.code <- character()
no.coordinate <- character()
flagged.gbif <- character()
final.results <- data.frame(species=character(0), num.records=numeric(0))
#Part 2
less.than.five <- character()
missing.unique.species <- character()
species.less.auc <- character()
species.modelled <- character()
maxent.results <- data.frame(me.species=character(0), me.records=numeric(0), me.auc =numeric(0), null.auc = numeric(0))


#read study area shapefile.
countries <- rgdal::readOGR("D:/Github/BGE-SDMv2/GISDATA/Study Area SHP")
#assign a CRS.
P4S.latlon <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

#the species in the GAP analysis.
Diptera.files <- list.files("D:/Github/BGE_sdm/T4.1_Gap_Analysis/Diptera csv files", full.names = T, pattern = '[.]csv$')
Diptera.files

## Download precipitation/Bio/+++ data from WorldClim
##global.clim <- raster::getData("worldclim", var="bio", res=5, download=T, path="RDATA")

# read worldclim
files.present <- list.files('D:/Github/BGE_sdm/RDATA/clipped/', pattern="[.]asc$", full.names=T) # alternatives for pattern (c|C)(e|E)(l|L)$
files.present[1:19]

#make a stack with 19 bioclim variables
present.stack <- stack(files.present[1:19])
head(present.stack)

#transform present.stack in SpatialPointsDataFrame
present.df <- as.data.frame(present.stack, xy=T)
coordinates(present.df) <- ~x+y
gridded(present.df) <- T

# assign projection
present.df@proj4string <- P4S.latlon

# Add grid.index value
present.df$grid.index <- present.df@grid.index
str(present.df@data)

#plot Bio01 
image(present.df, 'fe_buffer_bio01') #Annual Mean Temperature

#### Step 1: Download and clean the GBIF data ####

#set the progress bar
pb <- txtProgressBar(min = 0, max = length(Diptera.files), style = 3)

#start the loop until the length of the file.
for(i in 1:length(Diptera.files)){

  #read each of the files in the loop.
  x <- read.csv(Diptera.files[i], header = TRUE, sep = ",")
  #record numbers of all species
  diptera.total <- diptera.total + dim(x)[1]

  #if there is no gap list in the file go next
  if(is.na(sum(x$status == "gap")) == TRUE){ 
    
    next
    
    } else {

    #select the species in the gap list.
    x <- x[x$status == "gap",]
    
    #record number of gap lists species
    diptera.gap <- diptera.gap + dim(x)[1]

    #start the loop for the each gap list species.
    for(j in 1:dim(x)[1]){

      #get the specieskey from the species name.
      speciesKey <- name_backbone(x$species[j])$speciesKey

      #SOME SPECIES CAN NOT GET A SPECIESKEY BUT THEN IT WAS CONTINUE TO DOWNLOAD RANDOM DATA!!!
      if ( is.null(speciesKey) == TRUE){
        
        no.key <- c(no.key, x$species[j])
        
      } else {
        
      #download the GBIF data.
      data.GBIF <- occ_data(taxonKey = speciesKey, limit=10000) # , hasCoordinate = TRUE, hasGeospatialIssue = FALSE

      # assign data to a variable
      data.GBIF <- data.GBIF$data
      
      #check the whether they have data in the GBIF or not.
      if(is.null(data.GBIF) == TRUE ){
        
        #if there no data in the gbif. add species name in the missing.gbif list.
        missing.gbif <- c(missing.gbif, x$species[j])
        
      } else {
      if(dim(data.GBIF)[1] < 5){ 
        
        raw.data.less <- c(raw.data.less, x$species[j])
        
      } else {
      
        raw.data.more <- c(raw.data.more, x$species[j])
        
        if(is.null(data.GBIF$countryCode) == TRUE ){
          
          no.country.code <- c(no.country.code, x$species[j])
          
        } else {
          
          if(is.null(data.GBIF$decimalLatitude) == TRUE | is.null(data.GBIF$decimalLongitude) == TRUE  ){
      
            no.coordinate <- c(no.coordinate, x$species[j])
            
          }else{
             
          #clean NAs in the coordinate.  
          m <- !is.na(data.GBIF$decimalLatitude)
          data.GBIF <- data.GBIF[m,]
          n <- !is.na(data.GBIF$decimalLongitude)
          data.GBIF <- data.GBIF[n,]

          #To use the country - coordinate mismatch test we need to convert the country from ISO2 to ISO3 format.
          #convert country code from ISO2c to ISO3c
          data.GBIF$countryCode <-  countrycode(data.GBIF$countryCode, origin =  'iso2c', destination = 'iso3c')
              
          #https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13152
          #flag problems
          flags <- clean_coordinates(x = data.GBIF,
                                    lon = "decimalLongitude",
                                    lat = "decimalLatitude",
                                    countries = "countryCode",
                                    species = "species",
                                    tests = c("centroids","countries","duplicates","gbif","seas", "equal", "institutions","zeros"))# "capitals", 
          
        #Paramethers for cleaning data:
        #capitals tests a radius around adm-0 capitals.
        #centroids tests a radius around country centroids
        #countries tests if coordinates are from the country indicated in the country column.
        #duplicates tests for duplicate records. 
        #equal tests for equal absolute longitude and latitude.
        #institutions tests a radius around known biodiversity institutions from instiutions. 
        #zeros tests for plain zeros, equal latitude and longitude and a radius around the point 0/0
        #seas tests if coordinates fall into the ocean.
        #gbif tests a one-degree radius around the GBIF headquarters in Copenhagen, Denmark.
        
        #Exclude problematic records
        dat_cl <- data.GBIF[flags$.summary,]
        
        #filter basis of records
        dat_cl <- filter(dat_cl, basisOfRecord == "HUMAN_OBSERVATION" |
                           basisOfRecord == "OBSERVATION" |
                           basisOfRecord == "PRESERVED_SPECIMEN")
        
        #check whether there is data after deleting flagged records.
        
        if(dim(dat_cl)[1] < 5 ){
          
          #record flagged gbif.
          flagged.gbif <- c(flagged.gbif, x$species[j])
          
        } else {
          
          #Prepare the data
          feGBIF <- dat_cl[, c('scientificName', 'taxonomicStatus', 'taxonKey', 'speciesKey', 'decimalLatitude', 'decimalLongitude')]
          
          #save the species data.
          write.csv(feGBIF, paste0('D:/Github/BGE_sdm/T4.1_Gap_Analysis/R/speciesFilesDiptera/', x[j, c("species")], '.csv'))
          
          #species <- read.csv(species.files[j], header = TRUE)
          
          # transform species file in SpatialPointsDataFrame
          coordinates(feGBIF) <- ~decimalLongitude+decimalLatitude
          feGBIF@proj4string <- P4S.latlon #assign CRS
          
          # Get abiotic data and remove duplicates

          # Get climate variables + grid.index
          species.abiotic <- over(feGBIF, present.df) 
          # Link species col and climate data
          feGBIF <- cbind(feGBIF, species.abiotic) 
          #check and remove duplicates
          duplicates <- duplicated(feGBIF@data[,c("scientificName", "grid.index")]) # Duplicates on grid.index
          species.unique <- feGBIF[!duplicates,] # remove duplicates
          species.unique <- na.omit(species.unique)

          #records species name and number of records to final results df.
          species <- x$species[j]
          num.records <- dim(species.unique@data)[1]
          newdf <- data.frame(species,num.records)
          final.results <- rbind(final.results, newdf)
          #write.csv(final.results, 'D:/Github/BGE_sdm/Output/diptera.final.results.csv', row.names=F)
                }
              }
            }
          }
        }
      }
    }
  }
setTxtProgressBar(pb, i)  
}
close(pb)



#check results
diptera.total
diptera.gap
no.key
missing.gbif 
raw.data.less
raw.data.more 
no.country.code 
no.coordinate 
flagged.gbif 




### Create empty mask layer
mask <- raster(files.present[1])
mask <- !is.na(mask) # all values to 1
mask[mask == 0] <- NA # zero values to NA
writeRaster(mask, filename  = "D:/Github/BGE_sdm/RDATA/mask/mask.asc", format = 'ascii', NAflag = -9999, overwrite = T)
#plot(mask, col='red')

#read gap-species
species.files <- list.files('D:/Github/BGE_sdm/T4.1_Gap_Analysis/R/Diptera 28 aralik/', pattern="[.]csv$", full.names=T) # alternatives for pattern (c|C)(e|E)(l|L)$
species.files

#### Part 2 ####

# Loop through all species presence files
pb <- txtProgressBar(min = 0, max = length(species.files), style = 3)

for(i in 1:length(species.files))  {
  
  # read species file i
  species <- read.csv(species.files[i], header = TRUE)
  
  if(dim(species)[1] < 5){
    
    less.than.four <- c(less.than.four, basename(species.files[i]))
    
  } else {
    
    # transform species file in SpatialPointsDataFrame
    coordinates(species) <- ~decimalLongitude+decimalLatitude
    species@proj4string <- P4S.latlon #assign CRS
    
    # Get abiotic data and remove duplicates
    # Get climate variables + grid.index
    species.abiotic <- over(species, present.df) 
    # Link species col and climate data
    species <- cbind(species, species.abiotic) 
    species <- cbind(species@coords, species@data)
    species <- na.omit(species)
    #check and remove duplicates
    duplicates <- duplicated(species[,c("scientificName", "grid.index")]) # Duplicates on grid.index
    species.unique <- species[!duplicates,] # remove duplicates
    species.unique <- na.omit(species.unique)
    
    if (dim(species.unique)[1] < 5) {
      
      #if unique points are less than 5; add species name to missing.species list.
      missing.unique.species <- c(missing.unique.species, basename(species.files[i]))
      
    } else {
      #read empty mask layer
      mask <- raster('D:/Github/BGE_sdm/RDATA/mask/mask.asc')
      coordinates(species.unique) <- ~decimalLongitude+decimalLatitude
      species.unique@proj4string <- P4S.latlon #assign CRS
      
      # add mask to present.df
      present.df@data$mask <- as.data.frame(stack('D:/Github/BGE_sdm/RDATA/mask/mask.asc')) # Add mask layer to spdf

      #prepare data to maxent.
      species.unique.maxent <- as.data.frame(cbind(species.unique@data[, c("scientificName")], species.unique@coords))
      names(species.unique.maxent) <- c("species", "lon", "lat")
      
      # Convert to coll.locs spatial Points Data Frame
      species.unique.maxent$lon <- as.numeric(species.unique.maxent$lon)
      species.unique.maxent$lat <- as.numeric(species.unique.maxent$lat)
      coordinates(species.unique.maxent) <- ~lon+lat
      proj4string(species.unique.maxent) <- P4S.latlon
      #plot(countries); points(species.unique.maxent, pch=19, cex=0.5, col='green')
      
      #### Step 3: Create a buffer around occurences: 500 km ####
      
      #create a buffer around occurences
      s <- circles(species.unique.maxent, d=500000, lonlat=TRUE) 
      # dissolve polygons
      pol <- gUnaryUnion(s@polygons) 
      # extract cell numbers for the circles
      v <- extract(mask, s@polygons, cellnumbers=T)
      # use rbind to combine the elements in list v
      v <- do.call(rbind, v)
      # remove ocean cells
      v <- unique(na.omit(v))
      # to display the results
      m <- mask
      m[] <- NA # empty mask
      m[as.vector(v[,1])] <- 1
      
      # Write mask for buffered areas
      writeRaster(m, filename  = "D:/Github/BGE_sdm/RDATA/mask/mask.500km.asc", format = 'ascii', NAflag = -9999, overwrite = T)
    
      # copy present.df
      present.species.df <- present.df 
      
      # read mask layer
      mask.buffer.df <- as.data.frame(stack('D:/Github/BGE_sdm/RDATA/mask/mask.500km.asc'), xy=T)
      
      # replace mask with mask.500km
      present.species.df$mask <- mask.buffer.df[,'layer'] 
      
      #omit nonvalues.
      present.species.df <- na.omit(present.species.df@data)
      
      #### Step 4: Select uncorrelated variables using VIF ####
      if( dim(present.species.df)[1] < 10000){
        
        # sample all area background points for the VIF
        s <- sample(1:dim(present.species.df)[1])
        
      } else {
        
        # sample 10k background points for the VIF
        s <- sample(1:(dim(present.species.df)[1]), 10000, replace=F)
      }
      
      sample.df <- present.species.df[s,]
      ### VIF ###
      usdm::vif(sample.df[,c(1:19)])
      # VIF threshold 10
      usdm::vifstep(sample.df[,c(1:19)], th=10) 
      # threshold=10
      keep.dat <- usdm::vifstep(sample.df[,c(1:19)], th=10) 
      keep.dat <- keep.dat@results$Variables

      # Add mask layer
      keep.dat <- c(keep.dat, 'mask') 
      sample.df.keep <- sample.df[, (colnames(sample.df) %in% keep.dat)]
      # Species dataframe for keep.dat
      # retrieve species records
      species.unique.df <- species.unique[, (colnames(species.unique@data) %in% keep.dat)]
      #Add mask column
      species.unique.df$mask <- 1
      
      #### Step 5: Run the Maxent Model ####
      ### Create directory for Maxent
      mainDirMaxent <- "D:/Github/BGE_sdm/maxentOutput/"
      #writeRaster(mask, filename  = "D:/Github/BGE_sdm/RDATA/clipped/mask.asc", format = 'ascii', NAflag = -9999, overwrite = T) # Add mask to abiotic data
      
      #SampleWithData (swd) dataframe
      swd <- rbind(species.unique.df@data, sample.df.keep)
      
      # presence/absence (pa) vector
      pa <- c(rep(1, nrow(species.unique.df@data)), rep(0, nrow(sample.df.keep))); length(pa) 
      
      #Run the Maxent Model
      me <- maxent(swd, pa, args = c("noproduct", "nothreshold", "nohinge", "noextrapolate", "outputformat=logistic", "jackknife", "responsecurves", "applyThresholdRule=10 percentile training presence", "projectionlayers=D:/Github/BGE_sdm/RDATA/clipped", "redoifexists"), path=file.path(mainDirMaxent))
      
      #change the name of the output file
      file.rename("D:/Github/BGE_sdm/maxentOutput/species_clipped_thresholded.asc", paste0("D:/Github/BGE_sdm/maxentOutput/", basename(species.files[i]), "_clipped_thresholded.asc"))

      maxentResults <- read.csv(paste(file.path(mainDirMaxent), '/', 'maxentResults.csv', sep=""))
      maxent.auc <- maxentResults$Training.AUC
      # get number of records
      vector <- maxentResults$X.Training.samples 
      vector <- as.vector(sort(vector))
      # define df similar to training data
      x <- (sample.df.keep) 
      
      #### Step 6: Run null-model from source ####
      nm <- nullModel(x, n = vector, rep = 99)
      nm # shows the evaluations of the 'rep' null models created
      auc <- sapply(nm, function(x){slot(x,'auc')})# get just the auc values of  the null models
      auc <- auc[order(auc, decreasing = TRUE)]
      write.csv(auc, paste(file.path(mainDirMaxent), '/',basename(species.files[i]), 'nm_auc.csv', sep=""))
      maxentResults$nm <- auc[5] # Add null-model value to maxentResults
      write.csv(maxentResults, file=paste(file.path(mainDirMaxent), '/',basename(species.files[i]), 'maxentResults.csv', sep=""))
      
      #records AUCs of maxent and null model's 5th order.    
      me.species <- basename(species.files[i])
      me.records <- dim(species.unique@data)[1]
      me.auc <- maxent.auc
      null.auc <- auc[5]
      newdf <- data.frame(me.species,me.records,me.auc,null.auc)
      maxent.results <- rbind(maxent.results, newdf)
      #write.csv(maxent.results, 'D:/Github/BGE_sdm/Output/maxent.results.diptera.28.csv', row.names=F)
      
      setTxtProgressBar(pb, i)
      

      #check and record if maxent auc is higher than null model
      if (maxent.auc < auc[5]){

        species.less.auc <- c(species.less.auc, basename(species.files[i]))
        
      }else{
      
        species.modelled <- c(species.modelled, basename(species.files[i]))
        
      }
    }
  }
}
close(pb)


#read maxent projected probability maps
species.diversity.list <- list.files("D:/Github/BGE_sdm/maxentOutput/", pattern="_clipped_thresholded[.]asc$", full.names=T) # alternatives for pattern (c|C)(e|E)(l|L)$
species.diversity.list

#transform species.diversity.df in SpatialPointsDataFrame
species.diversity.df <- stack(species.diversity.list)
species.diversity.df <- as.data.frame(species.diversity.df, xy=T)
coordinates(species.diversity.df) <- ~x+y
gridded(species.diversity.df) <- T
species.diversity.df@proj4string <- P4S.latlon

head(species.diversity.df)
str(species.diversity.df)
class(species.diversity.df)
species.diversity.df@coords

species.diversity.df@data$diversity <- rowSums(species.diversity.df@data)
image(species.diversity.df, 'diversity') 