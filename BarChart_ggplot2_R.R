# Load the necessary packages
library(ggplot2)
library(tidyverse)
library(dplyr)


# Read in the data and filter out rows with less than 5 records
df <- read.csv("D:/Github/BGE_SDM_Pipeline/Output/HymenopteraMaxentResults.csv") %>%
  filter(me.records >= 5)


#group
df$group <- cut(df$me.records, breaks = c(0,5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 20, 25, 50, 100, 250, 500, Inf), 
                labels = c("5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16-20", "21-25", "26-50", "51-100", "101-250", "251-500", "501+"))

# Use group_by() and summarize() to get the count of each group
df <- df %>%
  group_by(group) %>%
  summarize(n = n())

# Create the bar plot using ggplot2
ggplot(df, aes(x = group, y = n)) +
  geom_col(fill = "#FFB52E", color = "black" , position = "dodge", width = 0.71) +
  geom_text(aes(label = n), vjust = -1, size = 3.5) +
  scale_y_continuous(limits = c(0,250))+
  labs(x = "Number of Unique Cells", y = "Count", title = "Hymenoptera", cex.title = 2,
       cex.labs = 11, cex.axis = 11,  cex.names = 11, las = 2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x =element_text(size=15),
        axis.title.y =element_text(size=15))+
  theme(plot.title = element_text(size = 18))
      

