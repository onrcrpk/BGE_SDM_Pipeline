# read CSV file
data <- read.csv("D:/deneme/maxent.results.liptera.csv")

# loop through rows of data frame
for (i in 1:nrow(data)) {
  # check value of "status" column
  if (data$status[i] == 1) {
    
    file.rename(paste0("D:/deneme/",data$me.species[i],"_clipped_thresholded.asc"), paste0("D:/deneme/valid",data$me.species[i],"_clipped_thresholded.asc"))
    

    }
  }
