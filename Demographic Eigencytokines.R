library(data.table)
library(ggplot2)
library(factoextra)
library(dplyr)
library(tidyverse)
library(gridExtra)
library(rgl) # for 3d plots 
library(plyr)
library(janitor)
library(ltm)
require(foreign)
require(nnet)
require(reshape2)
library(effects)
library(missForest)
library(tidyr)
library(readxl)

setwd("/Users/alexis/IEHS Dropbox/Rager Lab/VAlexis_Payton/1_Compartment_Analysis/Expt1.3_Cluster Distribution Analyses/Input")

Output = ("/Users/alexis/IEHS Dropbox/Rager Lab/Alexis_Payton/1_Compartment_Analysis/Expt1.3_Cluster Distribution Analyses/Output")
Wilcoxon_Input = ("/Users/alexis/IEHS Dropbox/Rager Lab/Alexis_Payton/1_Compartment_Analysis/Expt1.3_Cluster Distribution Analyses/Input")
# Upload master data files 
cytokines <- read_excel("CytokineData_102920.xlsx")
subjects <- data.frame(read_excel("SubjectInfo_102920.xlsx", sheet = 2))
malesubjects <- subjects %>% 
  filter(grepl("M", Sex)) # 24 subjects
maleIDs <- malesubjects$SubjectID
malecytokines <- cytokines %>% 
  filter(SubjectID %in% maleIDs)

# Separating the cytokine data into compartment dfs
malecytokines <- malecytokines %>% 
  group_by(Compartment) %>% 
  group_split
ELF <- malecytokines[[1]]
NLF <- malecytokines[[2]]
Serum <- malecytokines[[3]]
Sputum <- malecytokines[[4]]

# reshaping data 
ELF <- reshape2::dcast(ELF, SubjectID ~ Protein, value.var="Conc_pslog2") %>% 
  column_to_rownames("SubjectID") 
NLF <- reshape2::dcast(NLF, SubjectID ~ Protein, value.var="Conc_pslog2") %>% 
  column_to_rownames("SubjectID") 
Serum <- reshape2::dcast(Serum, SubjectID ~ Protein, value.var="Conc_pslog2") %>% 
  column_to_rownames("SubjectID") 
Sputum <- reshape2::dcast(Sputum, SubjectID ~ Protein, value.var="Conc_pslog2") %>% 
  column_to_rownames("SubjectID") 

#background filter eliminating any cytokines that are not expressed in a compartment  
NLF$I309 <- NULL
Sputum$I309 <- NULL 

#Z-score normalizing across cytokines within each compartment (using pseudolog2 of concentrations)
# the scale function operates across columns 
ELF_scaled <- ELF %>% 
  scale() %>% 
  as.data.frame()
NLF_scaled <- NLF %>% 
  scale() %>% 
  as.data.frame()
Serum_scaled <- Serum %>% 
  scale() %>% 
  as.data.frame()
Sputum_scaled <- Sputum %>% 
  scale() %>% 
  as.data.frame()

#Check that we now have a mean of ~0 and SD of 1 for each cytokine
colMeans(ELF_scaled)
apply(ELF_scaled, 2, sd)
colMeans(NLF_scaled)
apply(NLF_scaled, 2, sd)
colMeans(Serum_scaled)
apply(Serum_scaled, 2, sd)
colMeans(Sputum_scaled)
apply(Sputum_scaled, 2, sd)

#################################################################################################
#################################################################################################
#### Clustering cytokine data using k-means based on subject data - by compartment
#################################################################################################
#################################################################################################

# transforming the data
tELF <- as.data.frame(t(ELF))
tNLF <- as.data.frame(t(NLF))
tSerum <- as.data.frame(t(Serum))
tSputum <- as.data.frame(t(Sputum))

##########################################
# Estimate the optimal number of clusters
##########################################
# ELF
fviz_nbclust(tELF, kmeans, method = "wss") +
  labs(subtitle = "Elbow method") # suggests 4
fviz_nbclust(tELF, kmeans, method = "silhouette") + 
  labs(subtitle="Silhouette method") # suggests 2 

# NLF
fviz_nbclust(tNLF, kmeans, method = "wss") +
  labs(subtitle = "Elbow method") # suggests 3
fviz_nbclust(tNLF, kmeans, method = "silhouette") + 
  labs(subtitle="Silhouette method") # suggests 2

# Serum
fviz_nbclust(tSerum, kmeans, method = "wss") +
  labs(subtitle = "Elbow method") # suggests 2
fviz_nbclust(tSerum, kmeans, method = "silhouette") + 
  labs(subtitle="Silhouette method") # suggests 2 

# Sputum
fviz_nbclust(tSputum, kmeans, method = "wss") +
  labs(subtitle = "Elbow method") # suggests 4
fviz_nbclust(tSputum, kmeans, method = "silhouette") + 
  labs(subtitle="Silhouette method") # suggests 2

###########################################
# Testing different numbers of clusters (k)
##########################################

#########
#  ELF  #
#########
tELF_Cluster_k2 <- kmeans(tELF, centers = 2, nstart=25)
tELF_Cluster_k3 <- kmeans(tELF, centers = 3, nstart=25)
tELF_Cluster_k4 <- kmeans(tELF, centers = 4, nstart=25)
tELF_Cluster_k5 <- kmeans(tELF, centers = 5, nstart=25)
tELF_Cluster_k6 <- kmeans(tELF, centers = 6, nstart=25)
tELF_Cluster_k7 <- kmeans(tELF, centers = 7, nstart=25)
tELF_Cluster_k8 <- kmeans(tELF, centers = 8, nstart=25)
tELF_Cluster_k9 <- kmeans(tELF, centers = 9, nstart=25)
tELF_Cluster_k10 <- kmeans(tELF, centers = 10, nstart=25)

# Plots to compare
te2 <- fviz_cluster(tELF_Cluster_k2, data = tELF) + ggtitle("k = 2")
te3 <- fviz_cluster(tELF_Cluster_k3, data = tELF) + ggtitle("ELF, k = 3")
te4 <- fviz_cluster(tELF_Cluster_k4, data = tELF) + ggtitle("ELF, k = 4")
te5 <- fviz_cluster(tELF_Cluster_k5, data = tELF) + ggtitle("ELF, k = 5")
te6 <- fviz_cluster(tELF_Cluster_k6, data = tELF) + ggtitle("ELF, k = 6")
te7 <- fviz_cluster(tELF_Cluster_k7, data = tELF) + ggtitle("k = 7")
te8 <- fviz_cluster(tELF_Cluster_k8, data = tELF) + ggtitle("k = 8")
te9 <- fviz_cluster(tELF_Cluster_k9, data = tELF) + ggtitle("k = 9")
te10 <- fviz_cluster(tELF_Cluster_k10, data = tELF) + ggtitle("k = 10")

grid.arrange(te2, te3, te4, te5, te6, te7, te8, te9, te10, nrow = 3, top="Cytokine Clusters - ELF")
tELF_subs_plots <- arrangeGrob(te2, te3, te4, te5, te6, te7, te8, te9, te10, nrow=3, top="Cytokine Clusters - ELF") 
#ggsave(paste0(Output,"/", "ELF_cytokines.png"), tELF_subs_plots, width=30, height=24.47)

#########
#  NLF  #
#########
tNLF_Cluster_k2 <- kmeans(tNLF, centers = 2, nstart=25)
tNLF_Cluster_k3 <- kmeans(tNLF, centers = 3, nstart=25)
tNLF_Cluster_k4 <- kmeans(tNLF, centers = 4, nstart=25)
tNLF_Cluster_k5 <- kmeans(tNLF, centers = 5, nstart=25)
tNLF_Cluster_k6 <- kmeans(tNLF, centers = 6, nstart=25)
tNLF_Cluster_k7 <- kmeans(tNLF, centers = 7, nstart=25)
tNLF_Cluster_k8 <- kmeans(tNLF, centers = 8, nstart=25)
tNLF_Cluster_k9 <- kmeans(tNLF, centers = 9, nstart=25)
tNLF_Cluster_k10 <- kmeans(tNLF, centers = 10, nstart=25)

# Plots to compare
tn2 <- fviz_cluster(tNLF_Cluster_k2, data = tNLF) + ggtitle("k = 2")
tn3 <- fviz_cluster(tNLF_Cluster_k3, data = tNLF) + ggtitle("NLF, k = 3")
tn4 <- fviz_cluster(tNLF_Cluster_k4, data = tNLF) + ggtitle("NLF, k = 4")
tn5 <- fviz_cluster(tNLF_Cluster_k5, data = tNLF) + ggtitle("NLF, k = 5")
tn6 <- fviz_cluster(tNLF_Cluster_k6, data = tNLF) + ggtitle("NLF, k = 6")
tn7 <- fviz_cluster(tNLF_Cluster_k7, data = tNLF) + ggtitle("k = 7")
tn8 <- fviz_cluster(tNLF_Cluster_k8, data = tNLF) + ggtitle("k = 8")
tn9 <- fviz_cluster(tNLF_Cluster_k9, data = tNLF) + ggtitle("k = 9")
tn10 <- fviz_cluster(tNLF_Cluster_k10, data = tNLF) + ggtitle("k = 10")

grid.arrange(tn2, tn3, tn4, tn5, tn6, tn7, tn8, tn9, tn10, nrow = 3, top="Cytokine Clusters - NLF")
tNLF_subs_plots <- arrangeGrob(tn2, tn3, tn4, tn5, tn6, tn7, tn8, tn9, tn10, nrow=3, top="Cytokine Clusters - NLF") 
#ggsave(paste0(Output,"/", "NLF_cytokines.png"), tNLF_subs_plots, width=30, height=24.47)

###########
#  Serum  #
###########
tSerum_Cluster_k2 <- kmeans(tSerum, centers = 2, nstart=25)
tSerum_Cluster_k3 <- kmeans(tSerum, centers = 3, nstart=25)
tSerum_Cluster_k4 <- kmeans(tSerum, centers = 4, nstart=25)
tSerum_Cluster_k5 <- kmeans(tSerum, centers = 5, nstart=25)
tSerum_Cluster_k6 <- kmeans(tSerum, centers = 6, nstart=25)
tSerum_Cluster_k7 <- kmeans(tSerum, centers = 7, nstart=25)
tSerum_Cluster_k8 <- kmeans(tSerum, centers = 8, nstart=25)
tSerum_Cluster_k9 <- kmeans(tSerum, centers = 9, nstart=25)
tSerum_Cluster_k10 <- kmeans(tSerum, centers = 10, nstart=25)

# Plots to compare
tse2 <- fviz_cluster(tSerum_Cluster_k2, data = tSerum) + ggtitle("k = 2")
tse3 <- fviz_cluster(tSerum_Cluster_k3, data = tSerum) + ggtitle("Serum, k = 3")
tse4 <- fviz_cluster(tSerum_Cluster_k4, data = tSerum) + ggtitle("Serum, k = 4")
tse5 <- fviz_cluster(tSerum_Cluster_k5, data = tSerum) + ggtitle("Serum, k = 5")
tse6 <- fviz_cluster(tSerum_Cluster_k6, data = tSerum) + ggtitle("Serum, k = 6")
tse7 <- fviz_cluster(tSerum_Cluster_k7, data = tSerum) + ggtitle("k = 7")
tse8 <- fviz_cluster(tSerum_Cluster_k8, data = tSerum) + ggtitle("k = 8")
tse9 <- fviz_cluster(tSerum_Cluster_k9, data = tSerum) + ggtitle("k = 9")
tse10 <- fviz_cluster(tSerum_Cluster_k10, data = tSerum) + ggtitle("k = 10")

grid.arrange(tse2, tse3, tse4, tse5, tse6, tse7, tse8, tse9, tse10, nrow = 3, top="Cytokine Clusters - Serum")
tSerum_subs_plots <- arrangeGrob(tse2, tse3, tse4, tse5, tse6, tse7, tse8, tse9, tse10, nrow=3, top="Cytokine Clusters - Serum") 
#ggsave(paste0(Output,"/", "Serum_cytokines.png"), tSerum_subs_plots, width=30, height=24.47)

############
#  Sputum  #
############
tSputum_Cluster_k2 <- kmeans(tSputum, centers = 2, nstart=25)
tSputum_Cluster_k3 <- kmeans(tSputum, centers = 3, nstart=25)
tSputum_Cluster_k4 <- kmeans(tSputum, centers = 4, nstart=25)
tSputum_Cluster_k5 <- kmeans(tSputum, centers = 5, nstart=25)
tSputum_Cluster_k6 <- kmeans(tSputum, centers = 6, nstart=25)
tSputum_Cluster_k7 <- kmeans(tSputum, centers = 7, nstart=25)
tSputum_Cluster_k8 <- kmeans(tSputum, centers = 8, nstart=25)
tSputum_Cluster_k9 <- kmeans(tSputum, centers = 9, nstart=25)
tSputum_Cluster_k10 <- kmeans(tSputum, centers = 10, nstart=25)

# Plots to compare
tsp2 <- fviz_cluster(tSputum_Cluster_k2, data = tSputum) + ggtitle("k = 2")
tsp3 <- fviz_cluster(tSputum_Cluster_k3, data = tSputum) + ggtitle("Sputum, k = 3")
tsp4 <- fviz_cluster(tSputum_Cluster_k4, data = tSputum) + ggtitle("Sputum, k = 4")
tsp5 <- fviz_cluster(tSputum_Cluster_k5, data = tSputum) + ggtitle("Sputum, k = 5")
tsp6 <- fviz_cluster(tSputum_Cluster_k6, data = tSputum) + ggtitle("Sputum, k = 6")
tsp7 <- fviz_cluster(tSputum_Cluster_k7, data = tSputum) + ggtitle("k = 7")
tsp8 <- fviz_cluster(tSputum_Cluster_k8, data = tSputum) + ggtitle("k = 8")
tsp9 <- fviz_cluster(tSputum_Cluster_k9, data = tSputum) + ggtitle("k = 9")
tsp10 <- fviz_cluster(tSputum_Cluster_k10, data = tSputum) + ggtitle("k = 10")

grid.arrange(tsp2, tsp3, tsp4, tsp5, tsp6, tsp7, tsp8, tsp9, tsp10, nrow = 3, top="Cytokine Clusters - Sputum")
tSputum_subs_plots <- arrangeGrob(tsp2, tsp3, tsp4, tsp5, tsp6, tsp7, tsp8, tsp9, tsp10, nrow=3, top="Cytokine Clusters - Sputum") 
#ggsave(paste0(Output,"/", "Sputum_cytokines.png"), tSputum_subs_plots, width=30, height=24.47)

###########################################
# Exporting k=3,4,5,6 of all compartments
###########################################
grid.arrange(te3, tn3, tse3, tsp3,nrow = 2, top="Cytokines K=3 - MALES")
k3_cytokines <- arrangeGrob(te3, tn3, tse3, tsp3, nrow=2, top="Cytokines K=3 - MALES") 
ggsave(paste0(Output_Folder,"/", "k3_cytokines.png"), k3_cytokines, width=15, height=12.235)

grid.arrange(te4, tn4, tse4, tsp4,nrow = 2, top="Cytokines K=4")
k4_cytokines <- arrangeGrob(te4, tn4, tse4, tsp4, nrow=2, top="Cytokines K=4") 
ggsave(paste0(Output_Folder,"/", "k4_cytokines.png"), k4_cytokines, width=15, height=12.235)

grid.arrange(te5, tn5, tse5, tsp5,nrow = 2, top="Cytokines K=5")
k5_cytokines <- arrangeGrob(te5, tn5, tse5, tsp5, nrow=2, top="Cytokines K=5") 
ggsave(paste0(Output_Folder,"/", "k5_cytokines.png"), k5_cytokines, width=15, height=12.235)

grid.arrange(te6, tn6, tse6, tsp6,nrow = 2, top="Cytokines K=6")
k6_cytokines <- arrangeGrob(te6, tn6, tse6, tsp6, nrow=2, top="Cytokines K=6") 
ggsave(paste0(Output_Folder,"/", "k6_cytokines.png"), k6_cytokines, width=15, height=12.235)



#################################################################################################
#################################################################################################
#### Exporting final cytokine cluster assignments (k=3)
#################################################################################################
#################################################################################################
ELF_k3 <- as.data.frame(tELF_Cluster_k3$cluster) 
colnames(ELF_k3)[1] <- "Cluster"
NLF_k3 <- as.data.frame(tNLF_Cluster_k3$cluster) 
colnames(NLF_k3)[1] <- "Cluster"
Serum_k3 <- as.data.frame(tSerum_Cluster_k3$cluster) 
colnames(Serum_k3)[1] <- "Cluster"
Sputum_k3 <- as.data.frame(tSputum_Cluster_k3$cluster) 
colnames(Sputum_k3)[1] <- "Cluster"

#consider removing later
ELF_cyto = ELF
NLF_cyto = NLF
Sputum_cyto = Sputum
Serum_cyto = Serum

ELF_clus = ELF_k3
NLF_clus = NLF_k3
Sputum_clus = Sputum_k3
Serum_clus = Serum_k3

ELF_cyto <- as.data.frame(t(ELF_cyto)) 
NLF_cyto <- as.data.frame(t(NLF_cyto)) 
Serum_cyto <- as.data.frame(t(Serum_cyto)) 
Sputum_cyto <- as.data.frame(t(Sputum_cyto)) 


#renaming first column, grouping and splitting by "Cluster" column
ELF_clus <- ELF_clus %>% 
  rownames_to_column(var = 'Cytokine') %>% 
  group_by(Cluster) %>% 
  group_split
NLF_clus <- NLF_clus %>% 
  rownames_to_column(var = 'Cytokine') %>% 
  group_by(Cluster) %>% 
  group_split
Serum_clus <- Serum_clus %>% 
  rownames_to_column(var = 'Cytokine') %>% 
  group_by(Cluster) %>% 
  group_split
Sputum_clus <- Sputum_clus %>% 
  rownames_to_column(var = 'Cytokine') %>% 
  group_by(Cluster) %>% 
  group_split

#making dfs for each cluster for PCA analysis 
ELF_1 <- ELF_clus[[1]]
ELF_2 <- ELF_clus[[2]]
ELF_3 <- ELF_clus[[3]]

NLF_1 <- NLF_clus[[1]]
NLF_2 <- NLF_clus[[2]]
NLF_3 <- NLF_clus[[3]]

Serum_1 <- Serum_clus[[1]]
Serum_2 <- Serum_clus[[2]]
Serum_3 <- Serum_clus[[3]]

Sputum_1 <- Sputum_clus[[1]]
Sputum_2 <- Sputum_clus[[2]]
Sputum_3 <- Sputum_clus[[3]]


#making df with subjects' cytokine concentration data for each cluster 
ELF_1 <- ELF_cyto %>% 
  rownames_to_column("Cytokine") %>% 
  filter(Cytokine %in% ELF_1$Cytokine) %>% 
  column_to_rownames(var="Cytokine")
ELF_2 <- ELF_cyto %>% 
  rownames_to_column("Cytokine") %>% 
  filter(Cytokine %in% ELF_2$Cytokine) %>% 
  column_to_rownames(var="Cytokine")
ELF_3 <- ELF_cyto %>% 
  rownames_to_column("Cytokine") %>% 
  filter(Cytokine %in% ELF_3$Cytokine) %>% 
  column_to_rownames(var="Cytokine")

NLF_1 <- NLF_cyto %>% 
  rownames_to_column("Cytokine") %>% 
  filter(Cytokine %in% NLF_1$Cytokine) %>% 
  column_to_rownames(var="Cytokine")
NLF_2 <- NLF_cyto %>% 
  rownames_to_column("Cytokine") %>% 
  filter(Cytokine %in% NLF_2$Cytokine) %>% 
  column_to_rownames(var="Cytokine")
NLF_3 <- NLF_cyto %>% 
  rownames_to_column("Cytokine") %>% 
  filter(Cytokine %in% NLF_3$Cytokine) %>% 
  column_to_rownames(var="Cytokine")

Serum_1 <- Serum_cyto %>% 
  rownames_to_column("Cytokine") %>% 
  filter(Cytokine %in% Serum_1$Cytokine) %>% 
  column_to_rownames(var="Cytokine")
Serum_2 <- Serum_cyto %>% 
  rownames_to_column("Cytokine") %>% 
  filter(Cytokine %in% Serum_2$Cytokine) %>% 
  column_to_rownames(var="Cytokine")
Serum_3 <- Serum_cyto %>% 
  rownames_to_column("Cytokine") %>% 
  filter(Cytokine %in% Serum_3$Cytokine) %>% 
  column_to_rownames(var="Cytokine")

Sputum_1 <- Sputum_cyto %>% 
  rownames_to_column("Cytokine") %>% 
  filter(Cytokine %in% Sputum_1$Cytokine) %>% 
  column_to_rownames(var="Cytokine")
Sputum_2 <- Sputum_cyto %>% 
  rownames_to_column("Cytokine") %>% 
  filter(Cytokine %in% Sputum_2$Cytokine) %>% 
  column_to_rownames(var="Cytokine")
Sputum_3 <- Sputum_cyto %>% 
  rownames_to_column("Cytokine") %>% 
  filter(Cytokine %in% Sputum_3$Cytokine) %>% 
  column_to_rownames(var="Cytokine")

#PCA on each cluster, eigenvectors are in rotation -- for some reason had to convert everything to numeric  
pca_ELF_1 <- ELF_1 %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>% 
  prcomp()
pca_ELF_2 <- ELF_2 %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>% 
  prcomp()
pca_ELF_3 <- ELF_3 %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>%   
  prcomp()

pca_NLF_1 <- NLF_1 %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>% 
  prcomp()
pca_NLF_2 <- NLF_2 %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>% 
  prcomp()
pca_NLF_3 <- NLF_3 %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>% 
  prcomp()

pca_Serum_1 <- Serum_1 %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>% 
  prcomp()
pca_Serum_2 <- Serum_2 %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>% 
  prcomp()
pca_Serum_3 <- Serum_3 %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>% 
  prcomp()

pca_Sputum_1 <- Sputum_1 %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>% 
  prcomp()
pca_Sputum_2 <- Sputum_2 %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>% 
  prcomp()
pca_Sputum_3 <- Sputum_3 %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>% 
  prcomp()

#eigenvector dfs of first principal component 
eigencytokines_ELF_1 <- data.frame(pca_ELF_1$rotation[,"PC1"])
colnames(eigencytokines_ELF_1)[1] <- "Cluster1"
eigencytokines_ELF_2 <- data.frame(pca_ELF_2$rotation[,"PC1"])
colnames(eigencytokines_ELF_2)[1] <- "Cluster2"
eigencytokines_ELF_3 <- data.frame(pca_ELF_3$rotation[,"PC1"])
colnames(eigencytokines_ELF_3)[1] <- "Cluster3"

eigencytokines_NLF_1 <- data.frame(pca_NLF_1$rotation[,"PC1"])
colnames(eigencytokines_NLF_1)[1] <- "Cluster1"
eigencytokines_NLF_2 <- data.frame(pca_NLF_2$rotation[,"PC1"])
colnames(eigencytokines_NLF_2)[1] <- "Cluster2"
eigencytokines_NLF_3 <- data.frame(pca_NLF_3$rotation[,"PC1"])
colnames(eigencytokines_NLF_3)[1] <- "Cluster3"

eigencytokines_Serum_1 <- data.frame(pca_Serum_1$rotation[,"PC1"])
colnames(eigencytokines_Serum_1)[1] <- "Cluster1"
eigencytokines_Serum_2 <- data.frame(pca_Serum_2$rotation[,"PC1"])
colnames(eigencytokines_Serum_2)[1] <- "Cluster2"
eigencytokines_Serum_3 <- data.frame(pca_Serum_3$rotation[,"PC1"])
colnames(eigencytokines_Serum_3)[1] <- "Cluster3"

eigencytokines_Sputum_1 <- data.frame(pca_Sputum_1$rotation[,"PC1"])
colnames(eigencytokines_Sputum_1)[1] <- "Cluster1"
eigencytokines_Sputum_2 <- data.frame(pca_Sputum_2$rotation[,"PC1"])
colnames(eigencytokines_Sputum_2)[1] <- "Cluster2"
eigencytokines_Sputum_3 <- data.frame(pca_Sputum_3$rotation[,"PC1"])
colnames(eigencytokines_Sputum_3)[1] <- "Cluster3"


#collapse all eigencytokine dfs
eigencytokines_ELF <- cbind(eigencytokines_ELF_1, eigencytokines_ELF_2, eigencytokines_ELF_3)
eigencytokines_NLF <- cbind(eigencytokines_NLF_1, eigencytokines_NLF_2, eigencytokines_NLF_3)
eigencytokines_Serum <- cbind(eigencytokines_Serum_1, eigencytokines_Serum_2, eigencytokines_Serum_3)
eigencytokines_Sputum <- cbind(eigencytokines_Sputum_1, eigencytokines_Sputum_2, eigencytokines_Sputum_3)

#scale all eigencytokine dfs
eigencytokines_ELF_scaled <- as.data.frame(scale(eigencytokines_ELF))
eigencytokines_NLF_scaled <- as.data.frame(scale(eigencytokines_NLF))
eigencytokines_Serum_scaled <- as.data.frame(scale(eigencytokines_Serum))
eigencytokines_Sputum_scaled <- as.data.frame(scale(eigencytokines_Sputum))

#export all eigencytokine dfs
write.csv(eigencytokines_ELF, paste0(Wilcoxon_Input,"/", "eigencytokines_ELF_male.csv"), row.names=TRUE)
write.csv(eigencytokines_NLF, paste0(Wilcoxon_Input,"/", "eigencytokines_NLF_male.csv"), row.names=TRUE)
write.csv(eigencytokines_Serum, paste0(Wilcoxon_Input,"/", "eigencytokines_Serum_male.csv"), row.names=TRUE)
write.csv(eigencytokines_Sputum, paste0(Wilcoxon_Input,"/", "eigencytokines_Sputum_male.csv"), row.names=TRUE)

write.csv(eigencytokines_ELF_scaled, paste0(Wilcoxon_Input,"/", "eigencytokines_ELF_scaled_male.csv"), row.names=TRUE)
write.csv(eigencytokines_NLF_scaled, paste0(Wilcoxon_Input,"/", "eigencytokines_NLF_scaled_male.csv"), row.names=TRUE)
write.csv(eigencytokines_Serum_scaled, paste0(Wilcoxon_Input,"/", "eigencytokines_Serum_scaled_male.csv"), row.names=TRUE)
write.csv(eigencytokines_Sputum_scaled, paste0(Wilcoxon_Input,"/", "eigencytokines_Sputum_scaled_male.csv"), row.names=TRUE)

###############################################################################
#now getting eigen vectors for females so do everything over
###############################################################################
rm(list=ls())

setwd("/Users/alexis/IEHS Dropbox/Rager Lab/Alexis_Payton/1_Compartment_Analysis/Expt1.3_Cluster Distribution Analyses/Input")
getwd()
Output = ("/Users/alexis/IEHS Dropbox/Rager Lab/Alexis_Payton/1_Compartment_Analysis/Expt1.3_Cluster Distribution Analyses/Output")
Wilcoxon_Input = ("/Users/alexis/IEHS Dropbox/Rager Lab/Alexis_Payton/1_Compartment_Analysis/Expt1.3_Cluster Distribution Analyses/Input")

cytokines <- data.frame(read_excel("CytokineData_102920.xlsx", sheet = 2))
subjects <- data.frame(read_excel("SubjectInfo_102920.xlsx", sheet = 2))
femalesubjects <- subjects %>% 
  filter(grepl("F", Sex)) # 20 subjects
femaleIDs <- femalesubjects$SubjectID
femalecytokines <- cytokines %>% 
  filter(SubjectID %in% femaleIDs)

# Separating the cytokine data into compartment dfs
femalecytokines <- femalecytokines %>% 
  group_by(Compartment) %>% 
  group_split
ELF <- femalecytokines[[1]]
NLF <- femalecytokines[[2]]
Serum <- femalecytokines[[3]]
Sputum <- femalecytokines[[4]]

# reshaping data 
ELF <- reshape2::dcast(ELF, SubjectID ~ Protein, value.var="Conc_pslog2") %>% 
  column_to_rownames("SubjectID") 
NLF <- reshape2::dcast(NLF, SubjectID ~ Protein, value.var="Conc_pslog2") %>% 
  column_to_rownames("SubjectID") 
Serum <- reshape2::dcast(Serum, SubjectID ~ Protein, value.var="Conc_pslog2") %>% 
  column_to_rownames("SubjectID") 
Sputum <- reshape2::dcast(Sputum, SubjectID ~ Protein, value.var="Conc_pslog2") %>% 
  column_to_rownames("SubjectID") 

#background filter eliminating any cytokines that are not expressed in a compartment  
NLF$I309 <- NULL
Sputum$I309 <- NULL 

#Z-score normalizing across cytokines within each compartment (using pseudolog2 of concentrations)
# the scale function operates across columns 
ELF_scaled <- ELF %>% 
  scale() %>% 
  as.data.frame()
NLF_scaled <- NLF %>% 
  scale() %>% 
  as.data.frame()
Serum_scaled <- Serum %>% 
  scale() %>% 
  as.data.frame()
Sputum_scaled <- Sputum %>% 
  scale() %>% 
  as.data.frame()

#Check that we now have a mean of ~0 and SD of 1 for each cytokine
colMeans(ELF_scaled)
apply(ELF_scaled, 2, sd)
colMeans(NLF_scaled)
apply(NLF_scaled, 2, sd)
colMeans(Serum_scaled)
apply(Serum_scaled, 2, sd)
colMeans(Sputum_scaled)
apply(Sputum_scaled, 2, sd)

#################################################################################################
#################################################################################################
#### Clustering cytokine data using k-means based on subject data - by compartment
#################################################################################################
#################################################################################################

# transforming the data
tELF <- as.data.frame(t(ELF))
tNLF <- as.data.frame(t(NLF))
tSerum <- as.data.frame(t(Serum))
tSputum <- as.data.frame(t(Sputum))

##########################################
# Estimate the optimal number of clusters
##########################################
# ELF
fviz_nbclust(tELF, kmeans, method = "wss") +
  labs(subtitle = "Elbow method") # suggests 4
fviz_nbclust(tELF, kmeans, method = "silhouette") + 
  labs(subtitle="Silhouette method") # suggests 2 

# NLF
fviz_nbclust(tNLF, kmeans, method = "wss") +
  labs(subtitle = "Elbow method") # suggests 3
fviz_nbclust(tNLF, kmeans, method = "silhouette") + 
  labs(subtitle="Silhouette method") # suggests 2

# Serum
fviz_nbclust(tSerum, kmeans, method = "wss") +
  labs(subtitle = "Elbow method") # suggests 2
fviz_nbclust(tSerum, kmeans, method = "silhouette") + 
  labs(subtitle="Silhouette method") # suggests 2 

# Sputum
fviz_nbclust(tSputum, kmeans, method = "wss") +
  labs(subtitle = "Elbow method") # suggests 4
fviz_nbclust(tSputum, kmeans, method = "silhouette") + 
  labs(subtitle="Silhouette method") # suggests 2

###########################################
# Testing different numbers of clusters (k)
##########################################

#########
#  ELF  #
#########
tELF_Cluster_k2 <- kmeans(tELF, centers = 2, nstart=25)
tELF_Cluster_k3 <- kmeans(tELF, centers = 3, nstart=25)
tELF_Cluster_k4 <- kmeans(tELF, centers = 4, nstart=25)
tELF_Cluster_k5 <- kmeans(tELF, centers = 5, nstart=25)
tELF_Cluster_k6 <- kmeans(tELF, centers = 6, nstart=25)
tELF_Cluster_k7 <- kmeans(tELF, centers = 7, nstart=25)
tELF_Cluster_k8 <- kmeans(tELF, centers = 8, nstart=25)
tELF_Cluster_k9 <- kmeans(tELF, centers = 9, nstart=25)
tELF_Cluster_k10 <- kmeans(tELF, centers = 10, nstart=25)

# Plots to compare
te2 <- fviz_cluster(tELF_Cluster_k2, data = tELF) + ggtitle("k = 2")
te3 <- fviz_cluster(tELF_Cluster_k3, data = tELF) + ggtitle("ELF, k = 3")
te4 <- fviz_cluster(tELF_Cluster_k4, data = tELF) + ggtitle("ELF, k = 4")
te5 <- fviz_cluster(tELF_Cluster_k5, data = tELF) + ggtitle("ELF, k = 5")
te6 <- fviz_cluster(tELF_Cluster_k6, data = tELF) + ggtitle("ELF, k = 6")
te7 <- fviz_cluster(tELF_Cluster_k7, data = tELF) + ggtitle("k = 7")
te8 <- fviz_cluster(tELF_Cluster_k8, data = tELF) + ggtitle("k = 8")
te9 <- fviz_cluster(tELF_Cluster_k9, data = tELF) + ggtitle("k = 9")
te10 <- fviz_cluster(tELF_Cluster_k10, data = tELF) + ggtitle("k = 10")

grid.arrange(te2, te3, te4, te5, te6, te7, te8, te9, te10, nrow = 3, top="Cytokine Clusters - ELF")
tELF_subs_plots <- arrangeGrob(te2, te3, te4, te5, te6, te7, te8, te9, te10, nrow=3, top="Cytokine Clusters - ELF") 
ggsave(paste0(Output,"/", "ELF_cytokines_female.png"), tELF_subs_plots, width=30, height=24.47)

#########
#  NLF  #
#########
tNLF_Cluster_k2 <- kmeans(tNLF, centers = 2, nstart=25)
tNLF_Cluster_k3 <- kmeans(tNLF, centers = 3, nstart=25)
tNLF_Cluster_k4 <- kmeans(tNLF, centers = 4, nstart=25)
tNLF_Cluster_k5 <- kmeans(tNLF, centers = 5, nstart=25)
tNLF_Cluster_k6 <- kmeans(tNLF, centers = 6, nstart=25)
tNLF_Cluster_k7 <- kmeans(tNLF, centers = 7, nstart=25)
tNLF_Cluster_k8 <- kmeans(tNLF, centers = 8, nstart=25)
tNLF_Cluster_k9 <- kmeans(tNLF, centers = 9, nstart=25)
tNLF_Cluster_k10 <- kmeans(tNLF, centers = 10, nstart=25)

# Plots to compare
tn2 <- fviz_cluster(tNLF_Cluster_k2, data = tNLF) + ggtitle("k = 2")
tn3 <- fviz_cluster(tNLF_Cluster_k3, data = tNLF) + ggtitle("NLF, k = 3")
tn4 <- fviz_cluster(tNLF_Cluster_k4, data = tNLF) + ggtitle("NLF, k = 4")
tn5 <- fviz_cluster(tNLF_Cluster_k5, data = tNLF) + ggtitle("NLF, k = 5")
tn6 <- fviz_cluster(tNLF_Cluster_k6, data = tNLF) + ggtitle("NLF, k = 6")
tn7 <- fviz_cluster(tNLF_Cluster_k7, data = tNLF) + ggtitle("k = 7")
tn8 <- fviz_cluster(tNLF_Cluster_k8, data = tNLF) + ggtitle("k = 8")
tn9 <- fviz_cluster(tNLF_Cluster_k9, data = tNLF) + ggtitle("k = 9")
tn10 <- fviz_cluster(tNLF_Cluster_k10, data = tNLF) + ggtitle("k = 10")

grid.arrange(tn2, tn3, tn4, tn5, tn6, tn7, tn8, tn9, tn10, nrow = 3, top="Cytokine Clusters - NLF")
tNLF_subs_plots <- arrangeGrob(tn2, tn3, tn4, tn5, tn6, tn7, tn8, tn9, tn10, nrow=3, top="Cytokine Clusters - NLF") 
ggsave(paste0(Output,"/", "NLF_cytokines_female.png"), tNLF_subs_plots, width=30, height=24.47)

###########
#  Serum  #
###########
tSerum_Cluster_k2 <- kmeans(tSerum, centers = 2, nstart=25)
tSerum_Cluster_k3 <- kmeans(tSerum, centers = 3, nstart=25)
tSerum_Cluster_k4 <- kmeans(tSerum, centers = 4, nstart=25)
tSerum_Cluster_k5 <- kmeans(tSerum, centers = 5, nstart=25)
tSerum_Cluster_k6 <- kmeans(tSerum, centers = 6, nstart=25)
tSerum_Cluster_k7 <- kmeans(tSerum, centers = 7, nstart=25)
tSerum_Cluster_k8 <- kmeans(tSerum, centers = 8, nstart=25)
tSerum_Cluster_k9 <- kmeans(tSerum, centers = 9, nstart=25)
tSerum_Cluster_k10 <- kmeans(tSerum, centers = 10, nstart=25)

# Plots to compare
tse2 <- fviz_cluster(tSerum_Cluster_k2, data = tSerum) + ggtitle("k = 2")
tse3 <- fviz_cluster(tSerum_Cluster_k3, data = tSerum) + ggtitle("Serum, k = 3")
tse4 <- fviz_cluster(tSerum_Cluster_k4, data = tSerum) + ggtitle("Serum, k = 4")
tse5 <- fviz_cluster(tSerum_Cluster_k5, data = tSerum) + ggtitle("Serum, k = 5")
tse6 <- fviz_cluster(tSerum_Cluster_k6, data = tSerum) + ggtitle("Serum, k = 6")
tse7 <- fviz_cluster(tSerum_Cluster_k7, data = tSerum) + ggtitle("k = 7")
tse8 <- fviz_cluster(tSerum_Cluster_k8, data = tSerum) + ggtitle("k = 8")
tse9 <- fviz_cluster(tSerum_Cluster_k9, data = tSerum) + ggtitle("k = 9")
tse10 <- fviz_cluster(tSerum_Cluster_k10, data = tSerum) + ggtitle("k = 10")

grid.arrange(tse2, tse3, tse4, tse5, tse6, tse7, tse8, tse9, tse10, nrow = 3, top="Cytokine Clusters - Serum")
tSerum_subs_plots <- arrangeGrob(tse2, tse3, tse4, tse5, tse6, tse7, tse8, tse9, tse10, nrow=3, top="Cytokine Clusters - Serum") 
ggsave(paste0(Output,"/", "Serum_cytokines_female.png"), tSerum_subs_plots, width=30, height=24.47)

############
#  Sputum  #
############
tSputum_Cluster_k2 <- kmeans(tSputum, centers = 2, nstart=25)
tSputum_Cluster_k3 <- kmeans(tSputum, centers = 3, nstart=25)
tSputum_Cluster_k4 <- kmeans(tSputum, centers = 4, nstart=25)
tSputum_Cluster_k5 <- kmeans(tSputum, centers = 5, nstart=25)
tSputum_Cluster_k6 <- kmeans(tSputum, centers = 6, nstart=25)
tSputum_Cluster_k7 <- kmeans(tSputum, centers = 7, nstart=25)
tSputum_Cluster_k8 <- kmeans(tSputum, centers = 8, nstart=25)
tSputum_Cluster_k9 <- kmeans(tSputum, centers = 9, nstart=25)
tSputum_Cluster_k10 <- kmeans(tSputum, centers = 10, nstart=25)

# Plots to compare
tsp2 <- fviz_cluster(tSputum_Cluster_k2, data = tSputum) + ggtitle("k = 2")
tsp3 <- fviz_cluster(tSputum_Cluster_k3, data = tSputum) + ggtitle("Sputum, k = 3")
tsp4 <- fviz_cluster(tSputum_Cluster_k4, data = tSputum) + ggtitle("Sputum, k = 4")
tsp5 <- fviz_cluster(tSputum_Cluster_k5, data = tSputum) + ggtitle("Sputum, k = 5")
tsp6 <- fviz_cluster(tSputum_Cluster_k6, data = tSputum) + ggtitle("Sputum, k = 6")
tsp7 <- fviz_cluster(tSputum_Cluster_k7, data = tSputum) + ggtitle("k = 7")
tsp8 <- fviz_cluster(tSputum_Cluster_k8, data = tSputum) + ggtitle("k = 8")
tsp9 <- fviz_cluster(tSputum_Cluster_k9, data = tSputum) + ggtitle("k = 9")
tsp10 <- fviz_cluster(tSputum_Cluster_k10, data = tSputum) + ggtitle("k = 10")

grid.arrange(tsp2, tsp3, tsp4, tsp5, tsp6, tsp7, tsp8, tsp9, tsp10, nrow = 3, top="Cytokine Clusters - Sputum")
tSputum_subs_plots <- arrangeGrob(tsp2, tsp3, tsp4, tsp5, tsp6, tsp7, tsp8, tsp9, tsp10, nrow=3, top="Cytokine Clusters - Sputum") 
ggsave(paste0(Output,"/", "Sputum_cytokines_female.png"), tSputum_subs_plots, width=30, height=24.47)

###########################################
# Exporting k=3,4,5,6 of all compartments
###########################################
grid.arrange(te3, tn3, tse3, tsp3,nrow = 2, top="Cytokines K=3 - FEMALES")
k3_cytokines <- arrangeGrob(te3, tn3, tse3, tsp3, nrow=2, top="Cytokines K=3 - MALES") 
ggsave(paste0(Output,"/", "k3_cytokines.png"), k3_cytokines, width=15, height=12.235)

grid.arrange(te4, tn4, tse4, tsp4,nrow = 2, top="Cytokines K=4")
k4_cytokines <- arrangeGrob(te4, tn4, tse4, tsp4, nrow=2, top="Cytokines K=4") 
ggsave(paste0(Output,"/", "k4_cytokines.png"), k4_cytokines, width=15, height=12.235)

grid.arrange(te5, tn5, tse5, tsp5,nrow = 2, top="Cytokines K=5")
k5_cytokines <- arrangeGrob(te5, tn5, tse5, tsp5, nrow=2, top="Cytokines K=5") 
ggsave(paste0(Output,"/", "k5_cytokines.png"), k5_cytokines, width=15, height=12.235)

grid.arrange(te6, tn6, tse6, tsp6,nrow = 2, top="Cytokines K=6")
k6_cytokines <- arrangeGrob(te6, tn6, tse6, tsp6, nrow=2, top="Cytokines K=6") 
ggsave(paste0(Output,"/", "k6_cytokines.png"), k6_cytokines, width=15, height=12.235)



#################################################################################################
#################################################################################################
#### Exporting final cytokine cluster assignments (k=3)
#################################################################################################
#################################################################################################
ELF_k3 <- as.data.frame(tELF_Cluster_k3$cluster) 
colnames(ELF_k3)[1] <- "Cluster"
NLF_k3 <- as.data.frame(tNLF_Cluster_k3$cluster) 
colnames(NLF_k3)[1] <- "Cluster"
Serum_k3 <- as.data.frame(tSerum_Cluster_k3$cluster) 
colnames(Serum_k3)[1] <- "Cluster"
Sputum_k3 <- as.data.frame(tSputum_Cluster_k3$cluster) 
colnames(Sputum_k3)[1] <- "Cluster"

write.csv(ELF_k3, paste0(Output,"/", "_ELF_Clusters_female.csv"), row.names = FALSE)
write.csv(NLF_k3, paste0(Output,"/", "_NLF_Clusters_female.csv"), row.names = FALSE)
write.csv(Sputum_k3, paste0(Output,"/", "_Sputum_Clusters_female.csv"), row.names = FALSE)
write.csv(Serum_k3, paste0(Output,"/", "_Serum_Clusters_female.csv"), row.names = FALSE)

#consider removing later
ELF_cyto = ELF
NLF_cyto = NLF
Sputum_cyto = Sputum
Serum_cyto = Serum

ELF_clus = ELF_k3
NLF_clus = NLF_k3
Sputum_clus = Sputum_k3
Serum_clus = Serum_k3

ELF_cyto <- as.data.frame(t(ELF_cyto)) 
NLF_cyto <- as.data.frame(t(NLF_cyto)) 
Serum_cyto <- as.data.frame(t(Serum_cyto)) 
Sputum_cyto <- as.data.frame(t(Sputum_cyto)) 


#renaming first column, grouping and splitting by "Cluster" column
ELF_clus <- ELF_clus %>% 
  rownames_to_column(var = 'Cytokine') %>% 
  group_by(Cluster) %>% 
  group_split
NLF_clus <- NLF_clus %>% 
  rownames_to_column(var = 'Cytokine') %>% 
  group_by(Cluster) %>% 
  group_split
Serum_clus <- Serum_clus %>% 
  rownames_to_column(var = 'Cytokine') %>% 
  group_by(Cluster) %>% 
  group_split
Sputum_clus <- Sputum_clus %>% 
  rownames_to_column(var = 'Cytokine') %>% 
  group_by(Cluster) %>% 
  group_split

#making dfs for each cluster for PCA analysis 
ELF_1 <- ELF_clus[[1]]
ELF_2 <- ELF_clus[[2]]
ELF_3 <- ELF_clus[[3]]

NLF_1 <- NLF_clus[[1]]
NLF_2 <- NLF_clus[[2]]
NLF_3 <- NLF_clus[[3]]

Serum_1 <- Serum_clus[[1]]
Serum_2 <- Serum_clus[[2]]
Serum_3 <- Serum_clus[[3]]

Sputum_1 <- Sputum_clus[[1]]
Sputum_2 <- Sputum_clus[[2]]
Sputum_3 <- Sputum_clus[[3]]


#making df with subjects' cytokine concentration data for each cluster 
ELF_1 <- ELF_cyto %>% 
  rownames_to_column("Cytokine") %>% 
  filter(Cytokine %in% ELF_1$Cytokine) %>% 
  column_to_rownames(var="Cytokine")
ELF_2 <- ELF_cyto %>% 
  rownames_to_column("Cytokine") %>% 
  filter(Cytokine %in% ELF_2$Cytokine) %>% 
  column_to_rownames(var="Cytokine")
ELF_3 <- ELF_cyto %>% 
  rownames_to_column("Cytokine") %>% 
  filter(Cytokine %in% ELF_3$Cytokine) %>% 
  column_to_rownames(var="Cytokine")

NLF_1 <- NLF_cyto %>% 
  rownames_to_column("Cytokine") %>% 
  filter(Cytokine %in% NLF_1$Cytokine) %>% 
  column_to_rownames(var="Cytokine")
NLF_2 <- NLF_cyto %>% 
  rownames_to_column("Cytokine") %>% 
  filter(Cytokine %in% NLF_2$Cytokine) %>% 
  column_to_rownames(var="Cytokine")
NLF_3 <- NLF_cyto %>% 
  rownames_to_column("Cytokine") %>% 
  filter(Cytokine %in% NLF_3$Cytokine) %>% 
  column_to_rownames(var="Cytokine")

Serum_1 <- Serum_cyto %>% 
  rownames_to_column("Cytokine") %>% 
  filter(Cytokine %in% Serum_1$Cytokine) %>% 
  column_to_rownames(var="Cytokine")
Serum_2 <- Serum_cyto %>% 
  rownames_to_column("Cytokine") %>% 
  filter(Cytokine %in% Serum_2$Cytokine) %>% 
  column_to_rownames(var="Cytokine")
Serum_3 <- Serum_cyto %>% 
  rownames_to_column("Cytokine") %>% 
  filter(Cytokine %in% Serum_3$Cytokine) %>% 
  column_to_rownames(var="Cytokine")

Sputum_1 <- Sputum_cyto %>% 
  rownames_to_column("Cytokine") %>% 
  filter(Cytokine %in% Sputum_1$Cytokine) %>% 
  column_to_rownames(var="Cytokine")
Sputum_2 <- Sputum_cyto %>% 
  rownames_to_column("Cytokine") %>% 
  filter(Cytokine %in% Sputum_2$Cytokine) %>% 
  column_to_rownames(var="Cytokine")
Sputum_3 <- Sputum_cyto %>% 
  rownames_to_column("Cytokine") %>% 
  filter(Cytokine %in% Sputum_3$Cytokine) %>% 
  column_to_rownames(var="Cytokine")

#PCA on each cluster, eigenvectors are in rotation -- for some reason had to convert everything to numeric  
pca_ELF_1 <- ELF_1 %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>% 
  prcomp()
pca_ELF_2 <- ELF_2 %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>% 
  prcomp()
pca_ELF_3 <- ELF_3 %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>%   
  prcomp()

pca_NLF_1 <- NLF_1 %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>% 
  prcomp()
pca_NLF_2 <- NLF_2 %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>% 
  prcomp()
pca_NLF_3 <- NLF_3 %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>% 
  prcomp()

pca_Serum_1 <- Serum_1 %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>% 
  prcomp()
pca_Serum_2 <- Serum_2 %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>% 
  prcomp()
pca_Serum_3 <- Serum_3 %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>% 
  prcomp()

pca_Sputum_1 <- Sputum_1 %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>% 
  prcomp()
pca_Sputum_2 <- Sputum_2 %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>% 
  prcomp()
pca_Sputum_3 <- Sputum_3 %>% 
  lapply(as.numeric) %>% 
  as.data.frame() %>% 
  prcomp()

#eigenvector dfs of first principal component 
eigencytokines_ELF_1 <- data.frame(pca_ELF_1$rotation[,"PC1"])
colnames(eigencytokines_ELF_1)[1] <- "Cluster1"
eigencytokines_ELF_2 <- data.frame(pca_ELF_2$rotation[,"PC1"])
colnames(eigencytokines_ELF_2)[1] <- "Cluster2"
eigencytokines_ELF_3 <- data.frame(pca_ELF_3$rotation[,"PC1"])
colnames(eigencytokines_ELF_3)[1] <- "Cluster3"

eigencytokines_NLF_1 <- data.frame(pca_NLF_1$rotation[,"PC1"])
colnames(eigencytokines_NLF_1)[1] <- "Cluster1"
eigencytokines_NLF_2 <- data.frame(pca_NLF_2$rotation[,"PC1"])
colnames(eigencytokines_NLF_2)[1] <- "Cluster2"
eigencytokines_NLF_3 <- data.frame(pca_NLF_3$rotation[,"PC1"])
colnames(eigencytokines_NLF_3)[1] <- "Cluster3"

eigencytokines_Serum_1 <- data.frame(pca_Serum_1$rotation[,"PC1"])
colnames(eigencytokines_Serum_1)[1] <- "Cluster1"
eigencytokines_Serum_2 <- data.frame(pca_Serum_2$rotation[,"PC1"])
colnames(eigencytokines_Serum_2)[1] <- "Cluster2"
eigencytokines_Serum_3 <- data.frame(pca_Serum_3$rotation[,"PC1"])
colnames(eigencytokines_Serum_3)[1] <- "Cluster3"

eigencytokines_Sputum_1 <- data.frame(pca_Sputum_1$rotation[,"PC1"])
colnames(eigencytokines_Sputum_1)[1] <- "Cluster1"
eigencytokines_Sputum_2 <- data.frame(pca_Sputum_2$rotation[,"PC1"])
colnames(eigencytokines_Sputum_2)[1] <- "Cluster2"
eigencytokines_Sputum_3 <- data.frame(pca_Sputum_3$rotation[,"PC1"])
colnames(eigencytokines_Sputum_3)[1] <- "Cluster3"


#collapse all eigencytokine dfs
eigencytokines_ELF <- cbind(eigencytokines_ELF_1, eigencytokines_ELF_2, eigencytokines_ELF_3)
eigencytokines_NLF <- cbind(eigencytokines_NLF_1, eigencytokines_NLF_2, eigencytokines_NLF_3)
eigencytokines_Serum <- cbind(eigencytokines_Serum_1, eigencytokines_Serum_2, eigencytokines_Serum_3)
eigencytokines_Sputum <- cbind(eigencytokines_Sputum_1, eigencytokines_Sputum_2, eigencytokines_Sputum_3)

#scale all eigencytokine dfs
eigencytokines_ELF_scaled <- as.data.frame(scale(eigencytokines_ELF))
eigencytokines_NLF_scaled <- as.data.frame(scale(eigencytokines_NLF))
eigencytokines_Serum_scaled <- as.data.frame(scale(eigencytokines_Serum))
eigencytokines_Sputum_scaled <- as.data.frame(scale(eigencytokines_Sputum))

#export all eigencytokine dfs
write.csv(eigencytokines_ELF, paste0(Wilcoxon_Input,"/", "eigencytokines_ELF_female.csv"), row.names=TRUE)
write.csv(eigencytokines_NLF, paste0(Wilcoxon_Input,"/", "eigencytokines_NLF_female.csv"), row.names=TRUE)
write.csv(eigencytokines_Serum, paste0(Wilcoxon_Input,"/", "eigencytokines_Serum_female.csv"), row.names=TRUE)
write.csv(eigencytokines_Sputum, paste0(Wilcoxon_Input,"/", "eigencytokines_Sputum_female.csv"), row.names=TRUE)

write.csv(eigencytokines_ELF_scaled, paste0(Wilcoxon_Input,"/", "eigencytokines_ELF_scaled_female.csv"), row.names=TRUE)
write.csv(eigencytokines_NLF_scaled, paste0(Wilcoxon_Input,"/", "eigencytokines_NLF_scaled_female.csv"), row.names=TRUE)
write.csv(eigencytokines_Serum_scaled, paste0(Wilcoxon_Input,"/", "eigencytokines_Serum_scaled_female.csv"), row.names=TRUE)
write.csv(eigencytokines_Sputum_scaled, paste0(Wilcoxon_Input,"/", "eigencytokines_Sputum_scaled_female.csv"), row.names=TRUE)
