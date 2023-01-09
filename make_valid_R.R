library(dplyr)


# read CSV file
df <- read.csv("D:/Github/BGE_SDM_Pipeline/Output/DipteraMaxentResults.csv")


# Create a new column called "sign" that indicates whether A > B
df <- df %>% mutate(sign = ifelse(me.auc > null.auc, 1, 0))

# loop through rows of data frame
for (i in 1:nrow(df)) {
  # check value of "sing" column
  if (df$sign[i] == 1) {
    
    file.rename(paste0("D:/deneme/", df$me.species[i],"_clipped_thresholded.asc"), paste0("D:/deneme/valid_",df$me.species[i],"_clipped_thresholded.asc"))
    

    }
  }
