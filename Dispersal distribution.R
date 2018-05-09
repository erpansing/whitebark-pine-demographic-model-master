##-----------------------------------------------------------##
##                                                           ##
##           Seed dispersal distance distribution            ##
##                    Libby Pansing                          ##
##-----------------------------------------------------------##
library(dplyr)
library(ggplot2)

dispersal <- read.csv("/Users/elizabethpansing/Documents/Test/Field-et-al.-Model/WBP Demographic Model/Lorenz Cache Distance Data.csv")

ggplot(data = dispersal, aes(x = Distance_median)) +
  geom_density(fill = "#01665e") +
  labs(x = "Median dispersal distance", y = "Density") +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.title = element_text(size = 18, face = "bold"))+
  labs(x = "Median dispersal distance (km)", y = "Density")


hist(dispersal$Distance_median, breaks = 10)
