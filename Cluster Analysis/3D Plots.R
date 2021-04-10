library(readxl)
library(dplyr)
library(tidyverse)
library(rgl) # for 3d plots 
library(ggrepel)
library(plotly)

setwd('/Users/alexis/IEHS Dropbox/Rager Lab/Alexis_Payton/1_Compartment_Analysis/Expt1.5_Baseline Cluster Analyses/Output')

data = read_csv("020421NLF_PCA.csv") %>% rename(Protein = X1)
head(data)

plot_ly(data, x = ~PC1, y = ~PC2, z = ~PC3) %>%
  #consider adding the raw values
  add_markers() %>% 
  layout(scene = list(xaxis = list(title = 'PC1'), yaxis = list(title = 'PC2'), zaxis = list(title = 'PC3')))
