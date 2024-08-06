# WELCOME TO THE BEHAV3D TUMOR PROFILER (R script edition) Here you can freely
# and in your own environment execute the desired functions of the BEHAV3D Tumor
# Profiler pipeline. Remember to change any desired parameters to your own
# preferences Please set the working directories where required 
# You must execute the pre-processing before executing the modules.
# This script contains the 3 Modules in BEHAV3D Tumor Profiler:
# 1) Heterogeneity Module (line 770)
# 2) Large-scale phenotyping Module (line 1430)
# 3) Small-scale phenotyping Module (line 1883)
# 
# Each of this modules can be run independently from one another and may have common outputs
# The script is designed to output the results of each module separately, to avoid confusion
# For any further explanations, please refer to the github page (https://github.com/AlievaRios/BEHAV3D_Tumor_Profiler).
# THANK YOU FOR USING BEHAV3D TUMOR PROFILER!!!

library(plyr)
library(readr)
library(dplyr)
library(ggplot2)
setwd("C:/Users/Usuario/Documents/Emilio/IVM/data")

################################################################################
################################              ##################################
################################  Checkpoint  ##################################
################################              ##################################
################################################################################

# If the RDS files are in the working directory, it will skip the data input and initial pre-processing

if ( ! file.exists("master_cor.rds")){
  ## function to import the stats of the tracked tumor cells data from csv files extracted by Imaris
  read_plus <- function(flnm) {
    read_csv(flnm, skip = 3, col_types = cols()) %>% 
      mutate(filename = flnm)
  }
  ##import stats of interest
  ## set directory where csv files are located
  
  working_directory <- "C:/stats"
  setwd(working_directory)
  # import volume per organoid object
  pat = "*Displacement"
  files <- list.files(path = working_directory, pattern = pat, full.names = T, recursive = TRUE)
  files<-files[seq(1, length(files), 3)]
  disp_csv <- ldply(files, read_plus)
  # import distance to tumor
  pat = "*Distance_tumor"
  files <- list.files(path = working_directory, pattern = pat, full.names = T, recursive = TRUE)
  dist_tumor_csv <- ldply(files, read_plus)
  # import area per organoid object
  pat = "*Displacement_Delta_Length"
  files <- list.files(path = working_directory, pattern = pat, full.names = T, recursive = TRUE)
  dipl_delta_csv <- ldply(files, read_plus)
  
  # import area per organoid object
  pat = "*Displacement_Length"
  files <- list.files(path = working_directory, pattern = pat, full.names = T, recursive = TRUE)
  dipl_length_csv <- ldply(files, read_plus)
  # import position per organoid object
  pat = "*Position"
  files <- list.files(path = working_directory, pattern = pat, full.names = T, recursive = TRUE)
  pos_csv <- ldply(files, read_plus)
  
  # import speed
  pat = "*Speed"
  files <- list.files(path = working_directory, pattern = pat, full.names = T, recursive = TRUE)
  speed_csv <- ldply(files, read_plus)
  
  # import Time
  pat = "*Time"
  files <- list.files(path = working_directory, pattern = pat, full.names = T, recursive = TRUE)
  time_csv <- ldply(files, read_plus)
  
  ### ## make each ID unique (IDs imported from different imaging files can be repeated):
  category <- as.factor(dist_tumor_csv$filename)
  ranks <- rank(-table(category), ties.method="first")
  ranks <- as.data.frame(ranks)
  ranks$filename <- row.names(ranks)
  dist_tumor_csv <- left_join(dist_tumor_csv, ranks) 
  dist_tumor_csv$ID2 <- with(dist_tumor_csv, interaction(ID, ranks))
  
  category <- as.factor(pos_csv$filename)
  ranks <- rank(-table(category), ties.method="first")
  ranks <- as.data.frame(ranks)
  ranks$filename <- row.names(ranks)
  pos_csv <- left_join(pos_csv, ranks) 
  pos_csv$ID2 <- with(pos_csv, interaction(ID, ranks))
  
  category <- as.factor(time_csv$filename)
  ranks <- rank(-table(category), ties.method="first")
  ranks <- as.data.frame(ranks)
  ranks$filename <- row.names(ranks)
  time_csv <- left_join(time_csv, ranks) 
  time_csv$ID2 <- with(time_csv, interaction(ID, ranks))
  
  category <- as.factor(speed_csv$filename)
  ranks <- rank(-table(category), ties.method="first")
  ranks <- as.data.frame(ranks)
  ranks$filename <- row.names(ranks)
  speed_csv <- left_join(speed_csv, ranks) 
  speed_csv$ID2 <- with(speed_csv, interaction(ID, ranks))
  
  category <- as.factor(dipl_length_csv$filename)  ##this measures how much distnace did a cell move over, a cell can mover over a big distance and stay at the same position
  ranks <- rank(-table(category), ties.method="first")
  ranks <- as.data.frame(ranks)
  ranks$filename <- row.names(ranks)
  dipl_length_csv <- left_join(dipl_length_csv, ranks) 
  dipl_length_csv$ID2 <- with(dipl_length_csv, interaction(ID, ranks))
  
  category <- as.factor(dipl_delta_csv$filename)  ## this variable measures how much distance did the cell move from its origin
  ranks <- rank(-table(category), ties.method="first")
  ranks <- as.data.frame(ranks)
  ranks$filename <- row.names(ranks)
  dipl_delta_csv <- left_join(dipl_delta_csv, ranks) 
  dipl_delta_csv$ID2 <- with(dipl_delta_csv, interaction(ID, ranks))
  
  category <- as.factor(disp_csv$filename)  ## a bit similar to the previous ones, this measures the area that a cell covers overtime
  ranks <- rank(-table(category), ties.method="first")
  ranks <- as.data.frame(ranks)
  ranks$filename <- row.names(ranks)
  disp_csv <- left_join(disp_csv, ranks) 
  disp_csv$ID2 <- with(disp_csv, interaction(ID, ranks))
  
  
  ###Import BLOOD VESSELS stats (does not exist for all tracked positions so needs to be processed separtely.)
  working_directory2 <- "C:/BV_stats"
  setwd(working_directory2)
  
  files <- list.files(path = working_directory2, full.names = T, recursive = TRUE)
  BV_df <- ldply(files, read_plus)
  
  
  
  detach(package:plyr)
  
  ##create dataframe containing the statistics of interest
  master_1 <- left_join(dist_tumor_csv[,c(1,13)],time_csv[,c(1,5,6,8,10,11)], 
                        by.x="ID2", by.y="ID2")
  master_1 <- left_join(master_1,disp_csv[,c(1,11)], 
                        by.x="ID2", by.y="ID2")
  master_1 <- left_join(master_1,dipl_delta_csv[,c(1,11)], 
                        by.x="ID2", by.y="ID2")
  master_1 <- left_join(master_1,dipl_length_csv[,c(1,11)], 
                        by.x="ID2", by.y="ID2")
  master_1 <- left_join(master_1,speed_csv[,c(1,11)], 
                        by.x="ID2", by.y="ID2")
  master <- left_join(master_1,pos_csv[,c(1,2,3,14)], 
                      by.x="ID2", by.y="ID2")
  
  colnames(master) <- c("dist_tumor","ID2","Time","Tracked","TrackID", "filename","ranks","disp2","disp_d","disp_l","speed","x","y","z")

  # Extract mouse data (ID) using regexp (mother+sex+day)
  master$mouse <- master$filename
  master$mouse <- gsub("^.*[_/](\\d*[MF]\\d+).*", "\\1", master$mouse, perl=TRUE)
  
  # Extract position metadata
  master$position <- master$filename
  master$position <- gsub(".*_(CL\\d+_\\d+)_.*", "\\1", master$position, perl=TRUE)
  master$position<-with(master, interaction(mouse, position))
  
  
  # Extract class metadata
  master$class <- master$position
  master$class <- gsub(".*(CL\\d+)_.*", "\\1", master$class, perl=TRUE)
  
  master$filename<-NULL
  ## Round the time
  master$Time<-round(master$Time, digits=2)
  
  
  ##check that the imported tracks look ok
  g<-ggplot(master)+geom_path(aes(x, y, group=TrackID), size=0.2,
       arrow = arrow(ends = "last",type = "open",length = unit(0.01, "inches")))+
      theme_bw()+theme(aspect.ratio = 1)+facet_wrap(class~position,scales = "free")
  
  g
  
  ###Process blood vessels data:
  BV_df<-BV_df[c(1,6,7,8,9,11)]
  colnames(BV_df) <- c("dist_BV","Time","Tracked","TrackID", "ID","filename")
  
  ## rename mouse accroding to filename:
  BV_df$mouse<-BV_df$filename
  BV_df$mouse <- gsub("^.*[_/](\\d*[MF]\\d+).*", "\\1", BV_df$mouse, perl=TRUE)
  
  BV_df$position <- BV_df$filename
  BV_df$position <- gsub(".*_(CL\\d+_\\d+)_.*", "\\1", BV_df$position, perl=TRUE)
  BV_df$position <- with(BV_df, interaction(mouse, position))
  
  BV_df$class <- BV_df$position
  BV_df$class <- gsub(".*(CL\\d+)_.*", "\\1", BV_df$class, perl=TRUE)
  
  
  ## create ranks column based on master ranks:
  ranks_master<-subset(master, position %in% BV_df$position)
  ranks_master<-ranks_master[!duplicated(ranks_master$position),c(6,15)]
  colnames(ranks_master)[1]<-"ranks_bv"
  BV_df<- left_join(BV_df,ranks_master)
  
  ## keep only tracked cells and with a TrackID:
  BV_df<-subset(BV_df, Tracked=="Tracked")
  ## remove any NA:
  BV_df<-na.omit(BV_df)
  
  ## create unique Track_ID for Blood vessels dataframe:
  BV_df$Track2 <- with(BV_df, interaction(ranks_bv, TrackID))
  BV_df$Track2 <- gsub(".", '', BV_df$Track2, fixed = T)
  BV_df$Track2 <- as.numeric(as.character(BV_df$Track2))
  library(ggplot2)
  ##sumarize per TrackID the average, min and max distance to blood vessels
  ggplot(BV_df,aes(dist_BV))+geom_histogram()+facet_wrap(.~ position)
  BV_df<-na.omit(BV_df)
  BV_df$BV_contact<-ifelse(BV_df$dist_BV<4, 0, 1)
  BV_df_sum<-BV_df%>%group_by(Track2)%>%summarise(BV_sd=sd(dist_BV)/sqrt(length(dist_BV)),BV_mean=mean(dist_BV), BV_min=min(dist_BV), BV_max=max(dist_BV), BV_contact=mean(BV_contact))
  
  ### Import the SR101 data.
  library(plyr)
  library(readr)
  library(dplyr)
  ## function to import organoid data from csv files extracted by Imaris
  read_plus <- function(flnm) {
    read_csv(flnm, skip = 3, col_types = cols()) %>% 
      mutate(filename = flnm)
  }
  ##import Sr101 cells for the populations that have them
  
  working_directory <- "C:/SR101_stats"
  setwd(working_directory)
  # import volume per organoid object
  
  # import position per organoid object
  pat = "*Position"
  files <- list.files(path = working_directory, pattern = pat, full.names = T, recursive = TRUE)
  pos_csv <- ldply(files, read_plus)
  
  
  # import Time
  pat = "*Time"
  files <- list.files(path = working_directory, pattern = pat, full.names = T, recursive = TRUE)
  time_csv <- ldply(files, read_plus)
  
  ### ## make each TrackID unique (TrackIDs imported from different imaging files can be repeated):
  
  category <- as.factor(pos_csv$filename)
  ranks <- rank(-table(category), ties.method="first")
  ranks <- as.data.frame(ranks)
  ranks$filename <- row.names(ranks)
  pos_csv <- left_join(pos_csv, ranks) 
  pos_csv$ID2 <- with(pos_csv, interaction(ID, ranks))
  
  category <- as.factor(time_csv$filename)
  ranks <- rank(-table(category), ties.method="first")
  ranks <- as.data.frame(ranks)
  ranks$filename <- row.names(ranks)
  time_csv <- left_join(time_csv, ranks) 
  time_csv$ID2 <- with(time_csv, interaction(ID, ranks))
  
  detach(package:plyr)
  
  ##create dataframe containing the statistics of interest
  master_SR101 <- left_join(pos_csv[,c(1,2,3,10,12)],time_csv[,c(1,9)], 
                            by.x="ID2", by.y="ID2")
  
  
  colnames(master_SR101) <- c("x","y","z","filename","ID2","Time")
  
  ## rename mouse accroding to filename:
  # Extract mouse data (ID) using regexp (mother+sex+day)
  master_SR101$mouse <- master_SR101$filename
  master_SR101$mouse <- gsub("^.*[_/](\\d*[MF]\\d+).*", "\\1", master_SR101$mouse, perl=TRUE)
  
  # Extract position metadata
  master_SR101$position <- master_SR101$filename
  master_SR101$position <- gsub(".*_(CL\\d+_\\d+)_.*", "\\1", master_SR101$position, perl=TRUE)
  master_SR101$position<-with(master_SR101, interaction(mouse, position))
  
  # Extract class metadata
  master_SR101$class <- master_SR101$position
  master_SR101$class <- gsub(".*(CL\\d+)_.*", "\\1", master_SR101$class, perl=TRUE)
  
  master_SR101$filename<-NULL
  
  ## Convert imported time to hours:
  master_SR101$Time<-round(master_SR101$Time, digits=2)
  

  
  
  ### Import the MG data.
  library(plyr)
  library(readr)
  library(dplyr)
  ## function to import organoid data from csv files extracted by Imaris
  read_plus <- function(flnm) {
    read_csv(flnm, skip = 3) %>% 
      mutate(filename = flnm)
  }
  ##import MG cells for the populations that have them
  
  working_directory <- "C:/MG_stats"
  setwd(working_directory)
  # import volume per organoid object
  
  # import position per organoid object
  pat = "*Position"
  files <- list.files(path = working_directory, pattern = pat, full.names = T, recursive = TRUE)
  pos_csv <- ldply(files, read_plus)
  
  
  # import Time
  pat = "*Time"
  files <- list.files(path = working_directory, pattern = pat, full.names = T, recursive = TRUE)
  time_csv <- ldply(files, read_plus)
  
  ### ## make each TrackID unique (TrackIDs imported from different imaging files can be repeated):
  
  category <- as.factor(pos_csv$filename)
  ranks <- rank(-table(category), ties.method="first")
  ranks <- as.data.frame(ranks)
  ranks$filename <- row.names(ranks)
  pos_csv <- left_join(pos_csv, ranks) 
  pos_csv$ID2 <- with(pos_csv, interaction(ID, ranks))
  
  category <- as.factor(time_csv$filename)
  ranks <- rank(-table(category), ties.method="first")
  ranks <- as.data.frame(ranks)
  ranks$filename <- row.names(ranks)
  time_csv <- left_join(time_csv, ranks) 
  time_csv$ID2 <- with(time_csv, interaction(ID, ranks))
  
  detach(package:plyr)
  
  ##create dataframe containing the statistics of interest
  master_MG <- left_join(pos_csv[,c(1,2,3,10,12)],time_csv[,c(1,9)], 
                         by.x="ID2", by.y="ID2")
  
  
  colnames(master_MG) <- c("x","y","z","filename","ID2","Time")
  
  ## rename mouse accroding to filename:
  # Extract mouse data (ID) using regexp (mother+sex+day)
  
  master_MG$mouse<-master_MG$filename
  master_MG$mouse <- gsub("^.*[_/](\\d*[MF]\\d+).*", "\\1", master_MG$mouse, perl=TRUE)
  
  # Extract position metadata
  master_MG$position <- master_MG$filename
  master_MG$position <- gsub(".*_(CL\\d+_\\d+)_.*", "\\1", master_MG$position, perl=TRUE)
  master_MG$position<-with(master_MG, interaction(mouse, position))
  
  # Extract class metadata
  master_MG$class <- master_MG$position
  master_MG$class <- gsub(".*(CL\\d+)_.*", "\\1", master_MG$class, perl=TRUE)
  
  master_MG$filename<-NULL
  
  ## Convert imported time to hours:
  master_MG$Time<-round(master_MG$Time, digits=2)
  
  
## Pre-processing 
  
  library(dplyr)
  library(ggplot2)
  
  ### test that tracks imported are fine:
  master_test<-master%>%group_by(position)%>%mutate(x=x-min(x), y=y-min(y))
  g<-ggplot(master)+geom_path(aes(x, y, group=TrackID, colour=dist_tumor), arrow = arrow(ends = "last", type = "closed",length = unit(0.005, "inches")))+theme_bw()
  g<-g + theme(legend.position = "none")+facet_wrap(.~position, scales = "free")
  g
  
  ## To calculate interactions between cells 
  library(readr)
  library(spatstat)
  library(sp)
  library(gtools)
  library(reshape2)
  library(zoo)
  List2=list()
  for (m in unique(master$position)){
    distance_M4<-master[which(master$position==m),]
    
    List = list()
    for (i in unique(distance_M4$Time)){
      distancei <- distance_M4[distance_M4$Time==i,]
      distancei2<-as.data.frame(distancei[,c(3,11,12,13)])
      coordin<-ppx(distancei2, coord.type=c("t","s","s", "s"))
      dist1<- nndist(coordin)
      dist2<- nndist(coordin, k=2)
      dist3<- nndist(coordin, k=3)
      dist4<- nndist(coordin, k=4)
      dist5<- nndist(coordin, k=5)
      dist6<- nndist(coordin, k=6)
      dist7<- nndist(coordin, k=7)
      dist8<- nndist(coordin, k=8)
      dist9<- nndist(coordin, k=9)
      dist10<- nndist(coordin, k=10)
      dist<-cbind(dist1,dist2,dist3, dist4,dist5, dist6, dist7, dist8, dist9, dist10)  ### calculate distance to the three nearest neighbors
      dist<-data.frame(rowMeans(dist[,c(1:3)]),rowMeans(dist))
      distanceDFi_dist<-cbind(distancei, dist)
      List[[length(List)+1]] <-distanceDFi_dist ## store to list object
    }
    master_dist <- do.call(rbind, List) ## convert List to dta
    List2[[length(List2)+1]] <-master_dist
  }
  master_distance<-do.call(rbind, List2)
  colnames(master_distance)[c(17,18)] <- c("dist_3_neigh", "dist_10_neigh")
  
  ###remove the untracked tracks:
  master_distance<-subset(master_distance, Tracked=="Tracked")
  
  ## create a unique track_ID for master:
  master_distance$Track2 <- with(master_distance, interaction(ranks, TrackID))
  master_distance$Track2 <- gsub(".", '', master_distance$Track2, fixed = T)
  master_distance$Track2 <- as.numeric(as.character(master_distance$Track2))
  
  
  ##check that the imported tracks look ok
  g<-ggplot(master_distance)+geom_path(aes(x, y, group=Track2), size=0.2,
                              arrow = arrow(ends = "last",type = "open",length = unit(0.01, "inches")))+
    theme_bw()+theme(aspect.ratio = 1)+facet_wrap(class~position,scales = "free")
  
  g
  
  
  #### For the SR101 and MG datasets separately find what is the cell distance to these cells
  ##For SR101
  master_distance_SR101<-subset(master_distance, mouse %in% unique(master_SR101$mouse))
  library(RANN)
  ### Calculate distance from Tumor cells to SR101.
  ### Considering that olig2 doens;t move much, let's project the position of the birghtest timepoint from each movie to all timepoints.
  ### Identify birghtest timepoint by the number of SR101 objects:
  brightest_SR101<-master_SR101%>%group_by(Time, position)%>%summarise(n=n())%>%group_by( position)%>%filter(n==max(n))
  master_SR101<-subset(master_SR101, Time %in%brightest_SR101$Time)
  
  List2<-list()
  for (m in unique(master_distance_SR101$position)){
    master_pos <-master_distance_SR101%>%filter(position==m)%>%group_by(Time)%>%arrange(Time)
    master_SR101_pos <-master_SR101%>%filter(position==m)%>%group_by(Time)%>%arrange(Time)
    List1 = list() ## create list for ID of the contacting objects
    for(t in unique(master_pos$Time)){
      master_pos_t <-master_pos%>%filter(Time==t)
      ### calculate distance to 5 neigrest neighorbing SR101 in a radius of 30um
      cells_radius<-nn2(data=master_SR101_pos[c("x","y","z")], query = master_pos_t[c("x","y","z")],k=30, treetype = 'bd',searchtype =  "radius", radius = 30)
      cells_min<-nn2(data=master_SR101_pos[c("x","y","z")], query = master_pos_t[c("x","y","z")],k=1, treetype = 'bd',searchtype =  "standard")
      
      dist_table = data.frame(master_pos_t[c("Time","Track2","mouse","position","class")],cells_min[["nn.dists"]],cells_radius[["nn.dists"]])
      List1[[length(List1)+1]] <-dist_table  
      rm(temp)
    }
    cell_radius_t <- do.call(rbind, List1)
    List2[[length(List2)+1]] <-cell_radius_t 
    rm(cell_radius_t)
  }
  
  master_dist_SR101 <- do.call(rbind, List2)
  ### convert the distant cells to 0
  
  temp2<-master_dist_SR101[c("X1","X2","X3","X4","X5", "X6", "X7", "X8", "X9", "X10", "X11","X12","X13","X14","X15","X16","X17","X18","X19","X20", "X21","X22","X23","X24","X25","X26","X27","X28","X29","X30")]
  temp2[temp2 >1.340781e+153]<-0
  
  ### Calculate number of contacts:
  master_dist_SR101$n_SR101<-rowSums(temp2 != 0)
  ### renamen min_SR101:
  colnames(master_dist_SR101)[colnames(master_dist_SR101) == "cells_min...nn.dists..."] <- "min_SR101"
  
  ### bind information on SR101 distnaces to the main dataframe:
  master_dist_SR101[c("X1","X2","X3","X4","X5", "X6", "X7", "X8", "X9", "X10", "X11","X12","X13","X14","X15","X16","X17","X18","X19","X20", "X21","X22","X23","X24","X25","X26","X27","X28","X29","X30")]<-NULL
  ## For each Track2 find what is the minimal distance to SR101 and the max n of SR101 neiighbors in 30um diameter:
  master_distance_SR101<-master_dist_SR101%>%group_by(Track2, mouse, position, class)%>%summarize(n_SR101=max(n_SR101), min_SR101 =min(min_SR101))
  ### plot per position types the values of 
  p <- ggplot(master_distance_SR101)+geom_jitter(aes(x=class, y=n_SR101))+
    geom_bar(aes(x=class , y=n_SR101),alpha=0.5,stat = "summary") +ggtitle("number of neighboring SR101 cells per position type")+theme_classic()+
    theme(aspect.ratio=0.7)
  
  p
  
  p <- ggplot(master_distance_SR101)+geom_jitter(aes(x=class, y=min_SR101))+
    geom_boxplot(aes(x=class , y=min_SR101),alpha=0.5, outlier.colour = NA) +ggtitle("min distance to SR101 cells per position type")+theme_classic()+
    theme(aspect.ratio=0.7)
  
  p

  
  
  ##For MG
  master_distance_MG<-subset(master_distance, mouse %in% unique(master_MG$mouse))
  library(RANN)
  ### Calculate distance from Tumor cells to MG.
  ### Considering that olig2 doens;t move much, let's project the position of the birghtest timepoint from each movie to all timepoints.
  ### Identify birghtest timepoint by the number of MG objects:
  brightest_MG<-master_MG%>%group_by(Time, position)%>%summarise(n=n())%>%group_by( position)%>%filter(n==max(n))
  master_MG<-subset(master_MG, Time %in%brightest_MG$Time)
  
  List2<-list()
  for (m in unique(master_distance_MG$position)){
    master_pos <-master_distance_MG%>%filter(position==m)%>%group_by(Time)%>%arrange(Time)
    master_MG_pos <-master_MG%>%filter(position==m)%>%group_by(Time)%>%arrange(Time)
    List1 = list() ## create list for ID of the contacting objects
    for(t in unique(master_pos$Time)){
      master_pos_t <-master_pos%>%filter(Time==t)
      ### calculate distance to 5 neigrest neighorbing MG in a radius of 30um
      cells_radius<-nn2(data=master_MG_pos[c("x","y","z")], query = master_pos_t[c("x","y","z")],k=30, treetype = 'bd',searchtype =  "radius", radius = 30)
      cells_min<-nn2(data=master_MG_pos[c("x","y","z")], query = master_pos_t[c("x","y","z")],k=1, treetype = 'bd',searchtype =  "standard")
      
      dist_table = data.frame(master_pos_t[c("Time","Track2","mouse","position","class")],cells_min[["nn.dists"]],cells_radius[["nn.dists"]])
      List1[[length(List1)+1]] <-dist_table  
      rm(temp)
    }
    cell_radius_t <- do.call(rbind, List1)
    List2[[length(List2)+1]] <-cell_radius_t 
    rm(cell_radius_t)
  }
  
  master_dist_MG <- do.call(rbind, List2)
  ### convert the distant cells to 0
  
  temp2<-master_dist_MG[c("X1","X2","X3","X4","X5", "X6", "X7", "X8", "X9", "X10", "X11","X12","X13","X14","X15","X16","X17","X18","X19","X20", "X21","X22","X23","X24","X25","X26","X27","X28","X29","X30")]
  temp2[temp2 >1.340781e+153]<-0
  
  ### Calculate number of contacts:
  master_dist_MG$n_MG<-rowSums(temp2 != 0)
  ### renamen min_MG:
  colnames(master_dist_MG)[colnames(master_dist_MG) == "cells_min...nn.dists..."] <- "min_MG"
  
  ### bind information on MG distnaces to the main dataframe:
  master_dist_MG[c("X1","X2","X3","X4","X5", "X6", "X7", "X8", "X9", "X10", "X11","X12","X13","X14","X15","X16","X17","X18","X19","X20", "X21","X22","X23","X24","X25","X26","X27","X28","X29","X30")]<-NULL
  ## For each Track2 find what is the minimal distance to MG and the max n of MG neiighbors in 30um diameter:
  master_distance_MG<-master_dist_MG%>%group_by(Track2, mouse, position, class)%>%summarize(n_MG=max(n_MG), min_MG =min(min_MG))
  ### plot per position types the values of 
  p <- ggplot(master_distance_MG)+geom_jitter(aes(x=class, y=n_MG))+
    geom_bar(aes(x=class , y=n_MG),alpha=0.5,stat = "summary") +ggtitle("number of neighboring MG cells per position type")+theme_classic()+
    theme(aspect.ratio=0.7)
  
  p
  
  p <- ggplot(master_distance_MG)+geom_jitter(aes(x=class, y=min_MG))+
    geom_boxplot(aes(x=class , y=min_MG),alpha=0.5, outlier.colour = NA) +ggtitle("min distance to MG cells per position type")+theme_classic()+
    theme(aspect.ratio=0.7)
  
  p
  
  
  
  
  ### create empyt dataframe with combination of TIme of itnerest:
  max(master_distance$Time)
  length(seq(from=0, to=5.6, by=0.2))
  empty_df = master_distance[1:29,] ### create a dataframe with the same size as the sequence of Time of interest:
  empty_df[!is.na(empty_df)] <- NA ## make it empyt
  empty_df$Time<-seq(from=0, to=5.6, by=0.2)
  empty_df$Track2<-0.5
  
  
  ### now merge this fake dataframe with my dataframe of interest:
  master_distance<-rbind(master_distance, empty_df)
  library(tidyr)
  ### complete all the datasets so that they have the same timepoints:
  master_distance<-master_distance%>%complete(Track2,Time)
  ### remove the fake track that is not necessary anymore:
  
  master_distance<-subset(master_distance, !Track2==0.5)
  
  
  ##round the number, otherwise they are not equal
  master_distance$Time<-round(master_distance$Time,2)
  ##if there is any track2 and time that was duplicated remove the NA
  # Assuming your dataframe is named your_dataframe
  master_distance <- master_distance %>%
    group_by(Time, Track2) %>%
    dplyr::slice(if (n() == 1) 1 else which(!is.na(dist_tumor)))
  
  ### create a for loop to refill all the NA with interpolated values:
  ### create a list of the columns names that need to be refilled
  column_names<-names(master_distance)
  column_names<-subset(column_names,!column_names %in%c("TrackID","Track2", "ranks","ID","ID2","Time","position","mouse","class","pos", "Tracked","speed"))
  
  ### This is to test only ---------------
  ## create a first dataset with refilled values for speed:
  time_series<-acast(master_distance, Time ~ Track2, value.var='speed',fun.aggregate = mean)
  ## rownames timepoints:
  row.names(time_series)<-unique(master_distance$Time)
  ## get rid of NA
  time_series_zoo<-zoo(time_series, row.names(time_series))
  time_series_zoo<-na.approx(time_series_zoo) ## replace by interpolated value
  time_series<-as.matrix(time_series_zoo)
  time_series2<-melt(time_series)
  data<-time_series2[complete.cases(time_series2), ] 
  colnames(data)<-c("Time", "Track2", "speed")
  ### ----------
  for (i in column_names){
    time_series<-acast(master_distance, Time ~ Track2, value.var=i,fun.aggregate = mean)
    row.names(time_series)<-unique(master_distance$Time)
    ## get rid of NA
    time_series_zoo<-zoo(time_series,row.names(time_series))
    time_series_zoo<-na.approx(time_series_zoo) ## replace by last value
    time_series<-as.matrix(time_series_zoo)
    time_series2<-melt(time_series)
    new<-time_series2[complete.cases(time_series2), ] 
    data[ , ncol(data) + 1] <- new[3]                  # Append new column
    colnames(data)[ncol(data)] <- paste0(i)
    
  }
  ### check that the data is correct by plotting the evolution of a certain variable overtime -----
  p1 <-ggplot(data, aes(Time, disp2, group=Track2, color = as.factor(Track2))) + 
    geom_smooth(method = "loess",size = 0.5, se = F, alpha=0.3, span=1)+
    theme_bw() +
    theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=5), axis.title.y = element_text(size = 10), axis.title.x = element_text(size = 15), legend.text=element_text(size= 10))+
    ggtitle("disp2")+theme(aspect.ratio=1,legend.position = "none")
  
  p1
  ### --------

  
  
  ### add metadata info
  
  master_cor <- data
  master_metadata<- master_distance[c("position","mouse","class", "Track2")]
  master_metadata<-na.omit(master_metadata)
  master_metadata<-master_metadata[!duplicated(master_metadata$Track2),]
 
  library(dplyr)
  #detach(package:plyr, unload = TRUE)
  master_cor<- left_join(master_cor ,master_metadata)

  ### Filter timepoints of interest with same time interval:
  time_of_interest<-as.character(seq(from=0, to=5.6, by=0.2)) ### I don't know why but subsetting from a numeric vector doesn't work!!! So I convert to character
  master_cor$Time<-as.character(master_cor$Time)
  master_cor<-subset(master_cor, Time%in%time_of_interest)
  master_cor$Time<-as.numeric(master_cor$Time) ### now back to numeric
  
  
  
  ##check that the imported tracks look ok
  g<-ggplot(master_cor)+geom_path(aes(x, y, group=Track2), size=0.2,
                                       arrow = arrow(ends = "last",type = "open",length = unit(0.01, "inches")))+
    theme_bw()+theme(aspect.ratio = 1)+facet_wrap(class~position,scales = "free")
  
  g
  
  setwd("C:/results")
  
  saveRDS(master_distance, "C:/results/master_distance.rds")
  saveRDS(master_cor,"master_cor.rds")
  saveRDS(master_distance_MG,"master_distance_MG.rds")
  saveRDS(master_distance_SR101,"master_distance_SR101.rds")
  saveRDS(BV_df_sum,"BV_df_sum.rds")
}else{
  master_cor<-readRDS("master_cor.rds")
  master_distance<-readRDS("master_distance.rds")
  master_distance_MG<-readRDS("master_distance_MG.rds")
  master_distance_SR101<-readRDS("master_distance_SR101.rds")
  BV_df_sum<-readRDS("BV_df_sum.rds")
  print("data already imported")
}

#################################################################################
#############################                  ##################################
#############################  PRE-PROCESSING  ##################################
#############################                  ##################################
#################################################################################

### sumarize the length of the tracks:
## create relative time 
master_cor2<-master_cor %>% 
  group_by(Track2) %>%arrange(Time)%>%mutate(Time2 = Time - dplyr::first(Time))
ggplot(master_cor2,aes(Time2))+geom_histogram()+facet_wrap(.~ position)

max_time<-master_cor2%>%group_by(position)%>%summarise(max_Time=max(Time2))
min_time <- min(max_time$max_Time) #this is the minimal time that any position has
hist(max_time$max_Time)
master_cor2<-master_cor2 %>% 
  group_by(Track2)%>%filter(Time2<min_time)  ##Minimal time for a position is 2.6
## keep only the tracks that are of the same length  (max number of timepoints (13 in this case), this is how they were defined)
# Calculate the most frequent n for each group
subtrack_length<-master_cor2%>%group_by(Track2, position)%>%
  summarise(n=n_distinct(Time2), first_t=min(Time2), last_t=max(Time2))%>%
  ungroup()
hist(subtrack_length$n)
freq_n <- as.numeric(which.max(table(subtrack_length$n))) # Most frequent n
subtrack_length<-subtrack_length%>%filter(n==freq_n)  # we need to automate this to get the most frequent n
#Now we only keep these subtracks that have the same length:
master_cor2<-subset(master_cor2, Track2%in%subtrack_length$Track2)


library(scales)
### create a parameterd for distnace to tumor core
master_cor2<-master_cor2 %>% 
  group_by(position)%>%mutate(dist_tumor_1 = scales::rescale(dist_tumor, to=c(0,100)))

### normalize the distance to tumor factor per position position
master_cor2<-master_cor2 %>% 
  group_by(Track2) %>%arrange(Time2)%>%mutate(dist_tumor = dist_tumor - dplyr::first(dist_tumor))%>%ungroup()

#### check tracks
master_M4_test<-master_cor2%>%group_by(position, mouse)%>%mutate(x=x-min(x), y=y-min(y))
g<-ggplot(master_M4_test)+geom_path(aes(x, y, group=Track2, colour=dist_tumor), arrow = arrow(ends = "last", type = "closed",length = unit(0.005, "inches")))+theme_bw()
g<-g + theme(legend.position = "none")+scale_colour_gradient2(low = "blue", mid = "grey" , high = "red") +facet_wrap(.~position)
g

#number of unqiue tracks that are processed
length(unique(master_cor2$Track2))

### Now we need to do some feature engineering, since we make a sliding window analysis we need to requantify the values of displacement, speed, etc from scratch


## scaling
library(parallel)
library(dplyr)
library(dtwclust)
library(stats)
library(scales)
#detach(package:spatstat, unload = TRUE)
###normalize the data:
master_norm<-master_cor2%>%ungroup()
column_names2<-names(master_cor2)
## select what data to normalize
column_names2<-subset(column_names2,!column_names2 %in%c("Time","Track2","x","y","z","position","mouse","class","Time2", "dist_tumor_1", "dist_3_neigh", "dist_10_neigh", "iteration", "Track2"))
master_norm<-as.data.frame(master_norm)


# transform skwed variables
library(e1071)

# Calculate skewness for each variable
skew_values <- apply(master_norm[c(column_names2)], 2, skewness)

# Identify variables with skewness greater than a threshold (e.g., 0.5)
skewed_variables <- names(skew_values[abs(skew_values) > 2])
# Apply logarithmic transformation to skewed variables
master_norm[ , skewed_variables] <- log1p(master_norm[ , skewed_variables])


### perfrom PCA on the factors I want to use
master_scaled<-master_norm[c(column_names2)]%>%mutate(across(everything(), scale))%>%ungroup() ##first scale
##now to give more weight to the varibale of dist_tumor so the ability of cells to leave
#master_scaled[,c("dist_tumor")]<-master_scaled[,c("dist_tumor")]*2

master_pca<-prcomp(master_norm[c(column_names2)], scale = T)

master_pca2<-master_pca[["x"]] ##select only PC
master_pca3<-data.frame(master_norm[,c(2)],master_pca2[,c(1:3)]) ### take only first 3 components since they explain most variation

colnames(master_pca3)[1]<- "Track2"

pcaCharts <- function(x) {
  x.var <- x$sdev ^ 2
  x.pvar <- x.var/sum(x.var)
  print("proportions of variance:")
  print(x.pvar)
  
  par(mfrow=c(2,2))
  plot(x.pvar,xlab="Principal component", ylab="Proportion of variance explained", ylim=c(0,1), type='b')
  plot(cumsum(x.pvar),xlab="Principal component", ylab="Cumulative Proportion of variance explained", ylim=c(0,1), type='b')
  screeplot(x, main="Absolute Variance", xlab="Principal Component")
  screeplot(x,type="l", main="Absolute Variance")
  par(mfrow=c(1,1))
}

pcaCharts(master_pca)

################################################################################
#######################                         ################################
#######################   Heterogeneity Module  ################################
#######################                         ################################
################################################################################


## Create the output directory 

folder <- "C:/results/heterogeneity_module"

if (file.exists(folder)) {
  cat("The folder already exists")
} else {
  dir.create(folder)
}


###MULTIVARIATE
list_master_pca <- split(master_pca3[,-c(1)],master_pca3$Track2) ## split into list


setwd("C:/results")

if ( ! file.exists("matrix_distmat.rds")){
  ##parallel working
  # load parallel
  library(parallel)
  # create multi-process workers
  workers <- makeCluster(detectCores()-2)
  # load dtwclust in each one, and make them use 1 thread per worker
  invisible(clusterEvalQ(workers, {
    library(dtwclust)
    RcppParallel::setThreadOptions(1L)
  }))
  # register your workers, e.g. with doParallel
  require(doParallel)
  registerDoParallel(workers)
  
  ###MULTIVARIATE
  set.seed(123)
  distmat <- proxy::dist(list_master_pca, method = "dtw")
  saveRDS(distmat, file = "distmat.rds")
  
  matrix_distmat<-as.matrix(distmat)
  saveRDS(matrix_distmat, file = "matrix_distmat.rds")
  

}else{
  matrix_distmat <- readRDS("matrix_distmat.rds")
  }


Track2<-as.numeric(names(list_master_pca))

## UMAP and clustering

library(umap)
set.seed(123)
# Adjust parameters according to the experiment characteristics
umap_dist<- umap(matrix_distmat,n_components=2,input="dist",init = "random", 
                 n_neighbors=7,min_dist=0.5, spread=6,n_epochs=1000 , local_connectivity=1) 

umap_1 <- as.data.frame(umap_dist$`layout`) 
#scatter3Drgl(umap_1$V1, umap_1$V2, umap_1$V3)

plot(umap_dist$`layout`, # Plot the result
     col=rgb(0,0,0,alpha=0.1), 
     pch=19,
     asp=0.4, xlab="UMAP 1", ylab="UMAP 2",
     main="Raw UMAP")
umap_1 <- as.data.frame(umap_dist$`layout`) 
Track2_umap<-cbind(Track2,umap_1)
positiontype<-master_norm[,c("Track2", "position","class","mouse")]
positiontype<- positiontype[!duplicated(positiontype$Track2),]
umap_2 <- left_join(Track2_umap ,positiontype)

set.seed(123)
n_cluster=7
hc <- hclust(dist(umap_dist$`layout`, method = "euclidean"), method = "ward.D2")
hierarchical_clusters <- cutree(hc, k = n_cluster) 
umap_3 <- cbind(hierarchical_clusters, umap_2)

colnames(umap_3)[1]<- "cluster2"
##plot

# Make a color palette that adjusts to the cluster direction
# Generate a continuous palette with 100 colors
continuous_palette <- colorRampPalette(c("cyan4", "darkturquoise","darkorchid4","darkorchid1", "deeppink4","deeppink1", "goldenrod3", "gold"))(50)

# Select 8 colors from the continuous palette
mypalette_1<- continuous_palette[seq(1, length(continuous_palette), length.out = n_cluster)]

ggplot(umap_3, aes(x=V1, y=V2, color=as.factor(cluster2))) +  
  geom_point(size=2, alpha=0.8) + labs(color="cluster")+
  xlab("") + ylab("")+ 
  ggtitle("umap Cluster") +scale_color_manual(values=mypalette_1)+
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)


res_path<-"C:/results/heterogeneity_module" 

# shuffle the dataframe by rows
pdf(paste0(res_path, "/UMAP_large_Scale_positions.pdf"))


ggplot(umap_3, aes(x=V1, y=V2, color=as.factor(class))) +  
  geom_point(size=0.5, alpha=0.8) + labs(color="Environment_cluster")+
  xlab("") + ylab("")+ 
  ggtitle("umap Cluster") +scale_color_manual(values=c("cyan","gold","red"))+
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)+facet_wrap(mouse~position)



ggplot(umap_3, aes(x=V1, y=V2, color=as.factor(class))) +  
  geom_point(size=2, alpha=0.6) + labs(color="Tumor region phenotype")+
  xlab("") + ylab("")+ 
  ggtitle("distribution of large scale regions")+scale_color_manual(values=c("cyan","gold","red"))+
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)
dev.off()





### calculate average stats per track
master_class <- left_join(master_cor2 ,umap_3[c(1:4)], by = c("Track2" = "Track2"))
master_class<-master_class %>% 
  filter(cluster2 != 0)
master_class$cluster2<-as.numeric(master_class$cluster2)
master_class_sum<-master_class%>%group_by(position,class,mouse, Track2)%>%
  arrange(Time)%>%summarize(cluster2=mean(cluster2), V1=mean(V1), V2=mean(V2),
                            dist_tumor=mean(dist_tumor), dist_3_neigh=mean(dist_3_neigh),
                            dist_10_neigh=mean(dist_10_neigh), speed=mean(speed), disp2=mean(disp2), 
                            disp_d=mean(disp_d),disp_l=mean(disp_l), movement=last(dist_tumor), 
                            distance_to_tumor=last(dist_tumor_1),raw_dist_tumor=last(dist_tumor), 
                            Time=first(Time))%>%ungroup()


## Join the information on the SR101 and MG
master_class_sum <- left_join(master_class_sum,master_distance_MG[c("Track2","n_MG","min_MG")] ,by = c("Track2" = "Track2"))
master_class_sum <- left_join(master_class_sum,master_distance_SR101[c("Track2","n_SR101","min_SR101")] ,by = c("Track2" = "Track2"))
master_class_sum <- left_join(master_class_sum,BV_df_sum ,by = c("Track2" = "Track2"))

mid <- 0
master_class_sum$movement2<-ifelse(master_class_sum$movement>0,"away from tumor", ifelse(master_class_sum$movement==0, "no movement", "towards the tumor"))
master_class_sum<-master_class_sum[!is.na(master_class_sum$cluster2),]

saveRDS(master_class_sum,"master_class_sum.rds")
write.csv(master_class_sum,"master_class_sum.csv")

pdf(paste0(res_path, "/UMAP_direction.pdf"))
ggplot(master_class_sum, mapping=aes(x=V1, y=V2, color=movement)) +  
  geom_point(size=2, alpha=0.8) +labs(color="Direction")+
  xlab("") + ylab("") +
  ggtitle("Direction of movement relative to tumor core") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)+theme(aspect.ratio=1)+
  scale_color_gradient2(midpoint = mid, low = "blue",mid="grey80",
                        high = "red3")

dev.off()

library(dplyr)
library(scales)
library(viridis)
library(plotly)
library(pheatmap)

###Create heatmap
sum_all<-master_class_sum%>%group_by(cluster2)%>%
  summarise(displacement2 = median(disp2), 
            speed = median(speed),displacement_delta=median(disp_d),displacement_length=median(disp_l),max_speed=quantile(speed, 0.75), persistance=(displacement_length/displacement_delta),nearest_10_cell= median(dist_10_neigh),dist_tumor_core=median(distance_to_tumor),
            n_SR101=mean(n_SR101, na.rm = T),
            min_SR101=mean(min_SR101, na.rm = T),n_MG=mean(n_MG, na.rm = T),min_MG=mean(min_MG, na.rm = T),min_BV_dist=mean(BV_min, na.rm = T),sd_BV_dist=mean(BV_sd, na.rm = T),mean_BV_dist=mean(BV_mean, na.rm = T),
            movement=median(movement) )

Cluster_movement<-sum_all[,c(1,length(names(sum_all)))]
Cluster_movement <-Cluster_movement[order(Cluster_movement$movement),]

stdize = function(x, ...) {(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
sum_all[-c(1)] <- lapply(sum_all[-c(1)], stdize, na.rm = T)



df<-sum_all[,-c(1)]
df<-df[order(match(rownames(df), Cluster_movement$cluster2)), , drop = FALSE]%>%ungroup()
df<-as.data.frame(df)
rownames(df)<-Cluster_movement$cluster2

heat_m<-pheatmap(df, clustering_method = "complete", cluster_cols=F,cluster_rows=F,cellwidth=20, cellheight=20,
                 treeheight_row = 0, treeheight_col = 0,fontsize = 8,na_col = "grey50", main="Features",
                 angle_col = 90,color = colorRampPalette(c("deepskyblue1", "grey95", "deeppink2"))(50))

pdf(paste0(res_path, "cluster_heatmap_dynamic_features_2.pdf"))
df_cl_features<-df[,c("displacement2",     "speed", "max_speed", "displacement_delta" ,  "displacement_length", "persistance", "movement")]
rownames(df_cl_features)<-rownames(df)
pheatmap(df_cl_features, clustering_method = "complete", cluster_cols=F,cluster_rows=F,cellwidth=20, cellheight=20,
                 treeheight_row = 0, treeheight_col = 0,fontsize = 8,na_col = "grey50", main="Features used for clustering",
                 angle_col = 90,color = colorRampPalette(c("deepskyblue1", "grey95", "deeppink2"))(50))


dev.off()

print(Cluster_movement)  #3 per cluster we can see the direction

## Movement within the clusters in UMAP

## plot clusters giving them a color according to direction
# Convert df to a data frame
# Add cluster2 as a column
df1<-as.data.frame(df)

df1 <- mutate(df1, cluster2 = row.names(df))

# Order df1 by movement and reorder cluster2 accordingly
df1 <- df1 %>%
  arrange(movement) %>%
  mutate(cluster2 = factor(cluster2, levels = unique(cluster2)))
##set first order clusters based on the heatmap:



# Create a palette based on the movement values
pdf(paste0(res_path, "UMAP_cluster_and_other_features.pdf"))


ggplot(master_class_sum, aes(x = V1, y = V2, color = as.factor(cluster2))) +  
  geom_point(size = 2, alpha = 1) + 
  labs(color = "Cluster") +
  xlab("") + 
  ylab("") +
  scale_color_manual(values = setNames(mypalette_1, levels(df1$cluster2))) +  # Match palette to cluster2 levels  ggtitle("UMAP Cluster") +
  theme_light(base_size = 20) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        aspect.ratio = 1)

ggplot(umap_3, aes(x=V1, y=V2, color=as.factor(cluster2))) +  
  geom_point(size=1, alpha=0.8) + labs(color="cluster")+
  xlab("") + ylab("")+ 
  ggtitle("umap Cluster among mice") +scale_color_manual(values = setNames(mypalette_1, levels(df1$cluster2))) +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)+facet_wrap(mouse~.)

ggplot(master_class_sum, mapping=aes(x=V1, y=V2, color=distance_to_tumor)) +  
  geom_point(size=2, alpha=0.6) +
  xlab("") + ylab("") +
  ggtitle("UMAP localization distribution") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)+
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(4, "Spectral")))


ggplot(master_class_sum, mapping=aes(x=V1, y=V2, color=movement)) +  
  geom_point(size=2, alpha=0.6) +
  xlab("") + ylab("") +
  ggtitle("UMAP raw_movement") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)+
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(4, "Spectral")))


ggplot(master_class_sum, mapping=aes(x=V1, y=V2, color=speed)) +  
  geom_point(size=2, alpha=0.6) +
  xlab("") + ylab("") +
  ggtitle("UMAP speed distribution") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)+
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(4, "Spectral")))


ggplot(master_class_sum, mapping=aes(x=V1, y=V2, color=disp2)) +  
  geom_point(size=2, alpha=0.6) +
  xlab("") + ylab("") +
  ggtitle("UMAP disp2distribution") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)+
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(4, "Spectral")))


ggplot(master_class_sum, mapping=aes(x=V1, y=V2, color=disp_l)) +  
  geom_point(size=2, alpha=0.6) +
  xlab("") + ylab("") +
  ggtitle("UMAP disp_length") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)+
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(4, "Spectral")))


dev.off()

### plot enrichment on edge and core:
Number_cell_exp<-umap_3%>%group_by(position,class,mouse)%>%
  summarise(total_cell = n())
Percentage_clus<-left_join(Number_cell_exp,umap_3)
Percentage_clus <- Percentage_clus%>%group_by(cluster2,position,class,mouse)%>%
  summarise(total_cell = mean(total_cell), num_cluster=n())%>%mutate(percentage=num_cluster*100/total_cell)%>%ungroup()
Percentage_clus_2 <- Percentage_clus%>%group_by(cluster2, position,class,mouse)%>%summarise(se_percentage=sd(percentage)/sqrt(length(percentage)),percentage = mean(percentage))
### Plot circular map
Percentage_clus_2$cluster2<-as.factor(Percentage_clus_2$cluster2)
Percentage_clus_2$cluster2<- factor(Percentage_clus_2$cluster2, levels = levels(df1$cluster2))



Per2<-ggplot(Percentage_clus_2, aes(fill=as.factor(cluster2), y=percentage, x="")) + 
  geom_bar( stat="identity", position="fill", width = 0.5)+ coord_flip()+ scale_y_reverse()
Per2 <- Per2 + facet_grid(.~mouse)
Per2<-Per2+theme_void()+  scale_fill_manual(values=mypalette_1)+theme(strip.text.x = element_text(size = 8, angle = 90))
Per2


Per4<-ggplot(Percentage_clus_2, aes(fill=as.factor(mouse), y=percentage, x="")) + 
  geom_bar( stat="identity", position="fill", width = 0.5)+ coord_flip()+ scale_y_reverse()
Per4 <- Per4 + facet_grid(cluster2~.)
Per4<-Per4+theme_void()+theme(strip.text.x = element_text(size = 8, angle = 90))
Per4


pdf(paste0(res_path, "Pie chart.pdf"))
Per2
Per4
dev.off()


# Differences based on clusters
# Analyze the differences between behavioral clusters


### plot differences of the values per cluster
master_class_sum<-as.data.frame(master_class_sum)
master_class_sum$cluster2<-as.character(master_class_sum$cluster2)
df1$mean_movement<-df1$movement
master_class_sum2<-left_join(master_class_sum, df1[,c("cluster2","mean_movement")])
master_class_sum2$cluster2 = factor(master_class_sum2$cluster2, levels=df1$cluster2)


#### try an anova between clusters based speed 
p_speed <- ggplot(master_class_sum2,aes(x=as.factor(cluster2) , y=speed, group=cluster2,fill=cluster2)) +geom_jitter(width=0.2, alpha=0.5) +
  geom_boxplot(outlier.colour = NA, alpha=0.5) +ggtitle("speed per cluster")+theme_classic()+
  theme(aspect.ratio=0.7)+ scale_fill_manual(values=mypalette_1)+coord_cartesian(ylim = c(0,20))

p_speed

a_speed <- aov(speed ~ cluster2, data=master_class_sum2) 
summary(a_speed)
TukeyHSD(a_speed)
### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_speed), file=paste0(res_path,"Tukey_speed.txt"))
capture.output(summary(a_speed), file=paste0(res_path,"aov_speed.txt"))

#### try an anova between clusters based direction
p_direction <- ggplot(master_class_sum2,aes(x=as.factor(cluster2) , y=movement, group=cluster2,fill=cluster2)) +geom_jitter(width=0.2, alpha=0.5) +
  geom_boxplot(outlier.colour = NA, alpha=0.5) +ggtitle("speed per cluster")+theme_classic()+
  theme(aspect.ratio=0.7)+ scale_fill_manual(values=mypalette_1)+coord_cartesian(ylim = c(-15,20))

p_direction

a_direction <- aov(movement ~ cluster2, data=master_class_sum2) 
summary(a_direction)
TukeyHSD(a_direction)
### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_direction), file=paste0(res_path,"Tukey_direction.txt"))
capture.output(summary(a_direction), file=paste0(res_path,"aov_direction.txt"))



#### try an anova between clusters based on distance to dist_3_neigh
p_dist3_neigth <- ggplot(master_class_sum2,aes(x=as.factor(cluster2) , y=dist_3_neigh, group=cluster2, fill=cluster2)) +geom_jitter(width=0.2, alpha=0.5) +
  geom_boxplot(alpha=0.5, outlier.colour = NA)+ggtitle("dist 3 neigh per cluster")+theme_classic()+
  theme(aspect.ratio=0.7)+ scale_fill_manual(values=mypalette_1)+coord_cartesian(ylim = c(10,70))

p_dist3_neigth

a_dist_3_neigh <- aov(dist_3_neigh ~ cluster2, data=master_class_sum2) 
summary(a_dist_3_neigh)
TukeyHSD(a_dist_3_neigh)

capture.output(TukeyHSD(a_dist_3_neigh), file=paste0(res_path,"Tukey_dist_3_neigh.txt"))
capture.output(summary(a_dist_3_neigh), file=paste0(res_path,"aov_dist_3_neigh.txt"))

#### try an anova between clusters based on distance to dist_10_neigh
p_dist_10_neigh <- ggplot(master_class_sum2,aes(x=as.factor(cluster2) , y=dist_10_neigh, group=cluster2,fill=cluster2)) +geom_jitter(width=0.2, alpha=0.5) +
  geom_boxplot( outlier.colour = NA, alpha=0.5) +ggtitle("dist 10 neigh per cluster")+theme_classic()+
  theme(aspect.ratio=0.7)+ scale_fill_manual(values=mypalette_1)+coord_cartesian(ylim = c(20,90))

p_dist_10_neigh

a_dist_10_neigh <- aov(dist_10_neigh ~ cluster2, data=master_class_sum2) 
summary(a_dist_10_neigh)
TukeyHSD(a_dist_10_neigh)

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_dist_10_neigh), file=paste0(res_path,"Tukey_dist_10_neigh.txt"))
capture.output(summary(a_dist_10_neigh), file=paste0(res_path,"aov_dist_10_neigh.txt"))


#### try an anova between clusters based on dist to MG
p_min_MG <- ggplot(subset(master_class_sum2, !is.na(min_MG)) ,aes(x=as.factor(cluster2) , y=min_MG, group=cluster2, fill=cluster2))  +geom_jitter(width=0.3, alpha=0.5) +
  geom_boxplot( outlier.colour = NA, alpha=0.5) +ggtitle("distance to closest MG")+theme_classic()+
  theme(aspect.ratio=0.7)+ scale_fill_manual(values=mypalette_1)+coord_cartesian(ylim = c(2,45))

p_min_MG


a_min_MG <- aov(min_MG ~ cluster2, data= subset(master_class_sum2, !is.na(min_MG))) 
summary(a_min_MG)
TukeyHSD(a_min_MG)


### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_min_MG),  file=paste0(res_path,"Tukey_dist_min_MG.txt"))
capture.output(summary(a_min_MG),  file=paste0(res_path,"aov_dist_min_MG.txt"))


#### try an anova between clusters based on dist to SR101
p_min_SR101 <- ggplot(subset(master_class_sum2, !is.na(min_SR101)) ,aes(x=as.factor(cluster2) , y=min_SR101, group=cluster2, fill=cluster2))  +geom_jitter(width=0.3, alpha=0.5) +
  geom_boxplot(outlier.colour = NA, alpha=0.5) +ggtitle("distance to closest SR101")+theme_classic()+
  theme(aspect.ratio=0.7)+ scale_fill_manual(values=mypalette_1)+coord_cartesian(ylim = c(0,50))

p_min_SR101


a_min_SR101 <- aov(min_SR101 ~ cluster2, data= subset(master_class_sum2, !is.na(min_SR101))) 
summary(a_min_SR101)
TukeyHSD(a_min_SR101)



### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_min_SR101), file=paste0(res_path,"Tukey_dist_min_SR101.txt"))
capture.output(summary(a_min_SR101), file=paste0(res_path,"aov_dist_min_SR101.txt"))



#### try an anova between clusters based on n of  MG
p_n_MG <- ggplot(subset(master_class_sum2, !is.na(n_MG)) ,aes(x=as.factor(cluster2) , y=n_MG, group=cluster2, fill=cluster2))  +geom_jitter(width=0.3, alpha=0.5) +
  geom_boxplot( outlier.colour = NA, alpha=0.5) +ggtitle("n MG")+theme_classic()+
  theme(aspect.ratio=0.7)+ scale_fill_manual(values=mypalette_1)+coord_cartesian(ylim = c(0,15))

p_n_MG


a_n_MG <- aov(n_MG ~ cluster2, data= subset(master_class_sum2, !is.na(min_MG))) 
summary(a_n_MG)
TukeyHSD(a_n_MG)
### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_n_MG), file=paste0(res_path,"Tukey_dist_n_MG.txt") )
capture.output(summary(a_n_MG), file=paste0(res_path,"aov_dist_n_MG.txt") )




#### try an anova between clusters based on n of  SR101
p_nSR101 <- ggplot(subset(master_class_sum2, !is.na(n_SR101)) ,aes(x=as.factor(cluster2) , y=n_SR101, group=cluster2, fill=cluster2))  +geom_jitter(width=0.3, alpha=0.5) +
  geom_boxplot(outlier.colour = NA, alpha=0.5) +ggtitle("n SR101")+theme_classic()+
  theme(aspect.ratio=0.7)+ scale_fill_manual(values=mypalette_1)+coord_cartesian(ylim = c(0,20))

p_nSR101



a_nSR101 <- aov(n_SR101 ~ cluster2, data= subset(master_class_sum2, !is.na(n_SR101))) 
summary(a_nSR101)
TukeyHSD(a_nSR101)
### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_nSR101), file=paste0(res_path,"Tukey_dist_n_SR101.txt") )
capture.output(summary(a_nSR101), file=paste0(res_path,"aov_dist_n_SR101.txt") )



#### try an anova between clusters based on mean distance to BV

p_BV_mean <- ggplot(subset(master_class_sum2, !is.na(BV_mean)) ,aes(x=as.factor(cluster2) , y=BV_mean, group=cluster2, fill=cluster2))  +geom_jitter(width=0.3, alpha=0.5) +
  geom_boxplot(outlier.colour = NA, alpha=0.5) +ggtitle("BV_mean")+theme_classic()+
  theme(aspect.ratio=0.7)+ scale_fill_manual(values=mypalette_1)+coord_cartesian(ylim = c(0,20))

p_BV_mean



a_BV_mean <- aov(BV_mean ~ cluster2, data= subset(master_class_sum2, !is.na(BV_mean))) 
summary(a_BV_mean)
TukeyHSD(a_BV_mean)

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_BV_mean),file=paste0(res_path,"Tukey_dist_mean_BV.txt"))
capture.output(summary(a_BV_mean),file=paste0(res_path,"aov_dist_mean_BV.txt"))


#### try an anova between clusters based on min distance to BV

p_BV_min <- ggplot(subset(master_class_sum2, !is.na(BV_min)) ,aes(x=as.factor(cluster2) , y=BV_min, group=cluster2, fill=cluster2))  +geom_jitter(width=0.2, alpha=0.5) +
  geom_boxplot( outlier.colour = NA, alpha=0.5) +ggtitle("BV_min")+theme_classic()+
  theme(aspect.ratio=0.7)+ scale_fill_manual(values=mypalette_1)+coord_cartesian(ylim = c(0,20))

p_BV_min



a_BV_min <- aov(BV_min ~ cluster2, data= subset(master_class_sum2, !is.na(BV_min))) 
summary(a_BV_min)
TukeyHSD(a_BV_min)


### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_BV_min), file=paste0(res_path,"Tukey_dist_min_BV.txt"))
capture.output(summary(a_BV_min), file=paste0(res_path,"aov_dist_min_BV.txt"))



#### try an anova between clusters based on contact to BV

p_BV_contact <- ggplot(subset(master_class_sum2, !is.na(BV_min)) ,aes(x=as.factor(cluster2) , y=BV_contact, group=cluster2, fill=cluster2))  +geom_jitter(width=0.3, alpha=0.5) +
  geom_boxplot(outlier.colour = NA, alpha=0.5) +ggtitle("BV_contact")+theme_classic()+ scale_fill_manual(values=mypalette_1)+
  theme(aspect.ratio=0.7)

p_BV_contact



a_BV_contact <- aov(BV_contact ~ cluster2, data= subset(master_class_sum2, !is.na(BV_contact))) 
summary(a_BV_contact)
TukeyHSD(a_BV_contact)

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_BV_contact), file=paste0(res_path,"Tukey_dist_contact_BV.txt"))
capture.output(summary(a_BV_contact), file=paste0(res_path,"aov_dist_contact_BV.txt"))



#### try an anova between clusters based on sd distance to BV

p_BV_sd <- ggplot(subset(master_class_sum2, !is.na(BV_min)) ,aes(x=as.factor(cluster2) , y=BV_sd, group=cluster2, fill=cluster2))  +geom_jitter(width=0.3, alpha=0.5) +
  geom_boxplot( outlier.colour = NA, alpha=0.5) +ggtitle("Standart deviation distance")+theme_classic()+
  theme(aspect.ratio=0.7)+ scale_fill_manual(values=mypalette_1)+coord_cartesian(ylim = c(0,1.7))

p_BV_sd



a_BV_sd <- aov(BV_sd ~ cluster2, data= subset(master_class_sum2, !is.na(BV_min))) 
TukeyHSD(a_BV_sd)


### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_BV_sd), file=paste0(res_path,"Tukey_dist_sd_BV.txt"))
capture.output(summary(a_BV_sd), file=paste0(res_path,"aov_dist_sd_BV.txt"))



pdf(paste0(res_path,"per_cluster_features_comparison.pdf" )) ##adjust path

p_speed
p_dist3_neigth
p_dist_10_neigh
p_min_MG
p_min_SR101
p_n_MG
p_nSR101
p_BV_mean
p_BV_min
p_BV_contact
p_BV_sd

dev.off()



### export cell tracks per cluster

### extract track numbers for backprojection
backprojection_location<-paste0(res_path,"backprojection/" )

# Determine position of interest
position_interest<-"2430F11.CL3_1"

which( colnames(master_class)=="cluster2" ) ## column number
names_cell<-master_class[!duplicated(master_class$mouse),]
names_cell$position
Tracks_1<-master_class[which(master_class$position==position_interest),]
Tracks_1<-Tracks_1[!duplicated(Tracks_1$Track2),c(2,which( colnames(master_class)=="cluster2" ))]
##Join with information on unique TrackID for that experiment
# Generated during the data import
master_distance<-readRDS("C:/results/master_distance.rds")
master_distance<-master_distance[c("TrackID", "Track2")]%>% distinct(TrackID, .keep_all = TRUE)
Tracks_1<-left_join(Tracks_1,master_distance)
Tracks_1$Track2<-NULL
Track_1_list<-split(Tracks_1,Tracks_1$cluster2)

write(paste(as.character(Track_1_list), sep="' '", collapse=", "), paste0(backprojection_location, position_interest,".txt"))


###Plot cell trajectories:

## bind cluster information data to original dataset:

which( colnames(master_class)=="cluster2" ) ## column number
which( colnames(master_class)=="position" )
Tracks_1<-master_class
###plot all the tracks
Tracks_1<-Tracks_1[!duplicated(Tracks_1$Track2),c("cluster2", "Track2")]

Tracks_1$cluster2<-as.factor(Tracks_1$cluster2)
master_traject<-left_join(master_cor2,Tracks_1)

master_traject<-na.omit(master_traject)
master_traject<-subset(master_traject, cluster2!=0)
master_traject<-master_traject%>%group_by(position, mouse)%>%mutate(x_new=x-min(x), y_new=y-min(y))

master_traject$cluster2 = factor(master_traject$cluster2, levels=df1$cluster2)

g_all_track<-ggplot(master_traject)+geom_path(aes(x_new, y_new, group=Track2, 
                                        color=as.factor(cluster2)), size=0.2,
                                    arrow = arrow(ends = "last",type = "open",length = unit(0.01, "inches")))+scale_color_manual(values=mypalette_1 )+
  theme_bw()+theme(aspect.ratio = 1)+facet_wrap(class~position)
g_all_track


g_track_int1<-ggplot(subset(master_traject, position%in%position_interest))+geom_path(aes(x_new, y_new, group=Track2, color=as.factor(cluster2)), size=0.5,arrow = arrow(ends = "last",type = "open",length = unit(0.00, "inches")))+
  scale_color_manual(values=mypalette_1)+ggtitle(paste0(position_interest))+
  theme_bw()+theme(aspect.ratio = 1)
g_track_int1

g_track_int2<-ggplot(subset(master_traject, position%in%position_interest))+geom_path(aes(x_new, y_new, group=Track2), size=0.5,arrow = arrow(ends = "last",type = "open",length = unit(0.00, "inches")))+
  scale_color_manual(values=mypalette_1)+ggtitle(paste0(position_interest))+
  theme_bw()+theme(aspect.ratio = 1)
g_track_int2


pdf(paste0(res_path, "backprojection_pos.pdf"))
g_all_track
g_track_int1
g_track_int2
dev.off()


################################################################################
##################                                   ###########################
##################   Large-scale phenotyping Module  ###########################
##################                                   ###########################
################################################################################

# Create a directory for the results
folder <- "/content/results/large_scale_phenotyping_module"

if (file.exists(folder)) {
  cat("The '/content/results/large_scale_phenotyping_module' folder already exists")
} else {
  dir.create(folder)
}

###MULTIVARIATE
list_master_pca <- split(master_pca3[,-c(1)],master_pca3$Track2) ## split into list


setwd("C:/results")

if ( ! file.exists("matrix_distmat.rds")){
  ##parallel working
  # load parallel
  library(parallel)
  # create multi-process workers
  workers <- makeCluster(detectCores()-2)
  # load dtwclust in each one, and make them use 1 thread per worker
  invisible(clusterEvalQ(workers, {
    library(dtwclust)
    RcppParallel::setThreadOptions(1L)
  }))
  # register your workers, e.g. with doParallel
  require(doParallel)
  registerDoParallel(workers)
  
  ###MULTIVARIATE
  set.seed(123)
  distmat <- proxy::dist(list_master_pca, method = "dtw")
  saveRDS(distmat, file = "distmat.rds")
  
  matrix_distmat<-as.matrix(distmat)
  saveRDS(matrix_distmat, file = "matrix_distmat.rds")
  
  
}else{
  matrix_distmat <- readRDS("matrix_distmat.rds")
}


Track2<-as.numeric(names(list_master_pca))


library(umap)
set.seed(123)
# Modify the parameters according to the characteristics of the experiment
umap_dist<- umap(matrix_distmat,n_components=2,input="dist",init = "random", 
                 n_neighbors=7,min_dist=0.5, spread=6,n_epochs=1000 , local_connectivity=1)  

umap_1 <- as.data.frame(umap_dist$`layout`) 

plot(umap_dist$`layout`, # Plot the result
     col=rgb(0,0,0,alpha=0.1), 
     pch=19,
     asp=0.4, xlab="UMAP 1", ylab="UMAP 2",
     main="Raw UMAP")
umap_1 <- as.data.frame(umap_dist$`layout`) 

Track2_umap<-cbind(Track2,umap_1)
positiontype<-master_norm[,c("Track2", "position","class","mouse")]
positiontype<- positiontype[!duplicated(positiontype$Track2),]
umap_2 <- left_join(Track2_umap ,positiontype)

# Clustering
set.seed(123)
n_cluster=7
hc <- hclust(dist(umap_dist$`layout`, method = "euclidean"), method = "ward.D2")
hierarchical_clusters <- cutree(hc, k = n_cluster) 
umap_3 <- cbind(hierarchical_clusters, umap_2)

colnames(umap_3)[1]<- "cluster2"


res_path<-"C:/results/large_scale_phenotyping_module" 

# shuffle the dataframe by rows
pdf(paste0(res_path, "/UMAP_large_Scale_positions.pdf"))


ggplot(umap_3, aes(x=V1, y=V2, color=as.factor(class))) +  
  geom_point(size=0.5, alpha=0.8) + labs(color="Environment_cluster")+
  xlab("") + ylab("")+ 
  ggtitle("umap Cluster") +scale_color_manual(values=c("cyan","gold","red"))+
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)+facet_wrap(mouse~position)



ggplot(umap_3, aes(x=V1, y=V2, color=as.factor(class))) +  
  geom_point(size=2, alpha=0.6) + labs(color="Tumor region phenotype")+
  xlab("") + ylab("")+ 
  ggtitle("distribution of large scale regions")+scale_color_manual(values=c("cyan","gold","red"))+
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)
dev.off()





### calculate average stats per track
master_class <- left_join(master_cor2 ,umap_3[c(1:4)], by = c("Track2" = "Track2"))
master_class<-master_class %>% 
  filter(cluster2 != 0)
master_class$cluster2<-as.numeric(master_class$cluster2)
master_class_sum<-master_class%>%group_by(position,class,mouse, Track2)%>%
  arrange(Time)%>%summarize(cluster2=mean(cluster2), V1=mean(V1), V2=mean(V2),
                            dist_tumor=mean(dist_tumor), dist_3_neigh=mean(dist_3_neigh),
                            dist_10_neigh=mean(dist_10_neigh), speed=mean(speed), disp2=mean(disp2), 
                            disp_d=mean(disp_d),disp_l=mean(disp_l), movement=last(dist_tumor), 
                            distance_to_tumor=last(dist_tumor_1),raw_dist_tumor=last(dist_tumor), 
                            Time=first(Time))%>%ungroup()


## Join the information on the SR101 and MG
master_class_sum <- left_join(master_class_sum,master_distance_MG[c("Track2","n_MG","min_MG")] ,by = c("Track2" = "Track2"))
master_class_sum <- left_join(master_class_sum,master_distance_SR101[c("Track2","n_SR101","min_SR101")] ,by = c("Track2" = "Track2"))
master_class_sum <- left_join(master_class_sum,BV_df_sum ,by = c("Track2" = "Track2"))

mid <- 0
master_class_sum$movement2<-ifelse(master_class_sum$movement>0,"away from tumor", ifelse(master_class_sum$movement==0, "no movement", "towards the tumor"))
master_class_sum<-master_class_sum[!is.na(master_class_sum$cluster2),]

saveRDS(master_class_sum,"master_class_sum.rds")
write.csv(master_class_sum,"master_class_sum.csv")


library(dplyr)
library(scales)
library(viridis)
library(plotly)

###Create heatmap
sum_all<-master_class_sum%>%group_by(cluster2)%>%
  summarise(displacement2 = median(disp2), 
            speed = median(speed),displacement_delta=median(disp_d),displacement_length=median(disp_l),max_speed=quantile(speed, 0.75), persistance=(displacement_length/displacement_delta),nearest_10_cell= median(dist_10_neigh),dist_tumor_core=median(distance_to_tumor),
            n_SR101=mean(n_SR101, na.rm = T),
            min_SR101=mean(min_SR101, na.rm = T),n_MG=mean(n_MG, na.rm = T),min_MG=mean(min_MG, na.rm = T),min_BV_dist=mean(BV_min, na.rm = T),sd_BV_dist=mean(BV_sd, na.rm = T),mean_BV_dist=mean(BV_mean, na.rm = T),
            movement=median(movement) )

Cluster_movement<-sum_all[,c(1,length(names(sum_all)))]
Cluster_movement <-Cluster_movement[order(Cluster_movement$movement),]

stdize = function(x, ...) {(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
sum_all[-c(1)] <- lapply(sum_all[-c(1)], stdize, na.rm = T)



df<-sum_all[,-c(1)]
df<-df[order(match(rownames(df), Cluster_movement$cluster2)), , drop = FALSE]%>%ungroup()
df<-as.data.frame(df)
rownames(df)<-Cluster_movement$cluster2

## plot clusters giving them a color according to direction
# Convert df to a data frame
# Add cluster2 as a column
df1<-as.data.frame(df)

df1 <- mutate(df1, cluster2 = row.names(df))

# Order df1 by movement and reorder cluster2 accordingly
df1 <- df1 %>%
  arrange(movement) %>%
  mutate(cluster2 = factor(cluster2, levels = unique(cluster2)))

# Make a color palette that adjusts to the cluster direction
# Generate a continuous palette with 100 colors
continuous_palette <- colorRampPalette(c("cyan4", "darkturquoise","darkorchid4","darkorchid1", "deeppink4","deeppink1", "goldenrod3", "gold"))(50)

# Select 8 colors from the continuous palette
mypalette_1<- continuous_palette[seq(1, length(continuous_palette), length.out = n_cluster)]


### plot enrichment on edge and core:
Number_cell_exp<-umap_3%>%group_by(position,class,mouse)%>%
  summarise(total_cell = n())
Percentage_clus<-left_join(Number_cell_exp,umap_3)
Percentage_clus <- Percentage_clus%>%group_by(cluster2,position,class,mouse)%>%
  summarise(total_cell = mean(total_cell), num_cluster=n())%>%mutate(percentage=num_cluster*100/total_cell)%>%ungroup()
Percentage_clus_2 <- Percentage_clus%>%group_by(cluster2, position,class,mouse)%>%summarise(se_percentage=sd(percentage)/sqrt(length(percentage)),percentage = mean(percentage))
### Plot circular map
Percentage_clus_2$cluster2<-as.factor(Percentage_clus_2$cluster2)
Percentage_clus_2$cluster2<- factor(Percentage_clus_2$cluster2, levels = levels(df1$cluster2))
Per1<-ggplot(Percentage_clus_2, aes(fill=as.factor(cluster2), y=percentage, x="")) + 
  geom_bar( stat="identity", position="fill", width = 0.5)+ coord_flip()+ scale_y_reverse()
Per1 <- Per1 + facet_grid(class~mouse)
Per1<-Per1+theme_void()+ scale_fill_manual(values=mypalette_1)+theme(strip.text.x = element_text(size = 8, angle = 90))
Per1


Per3<-ggplot(Percentage_clus_2, aes(fill=as.factor(cluster2), y=percentage, x="")) + 
  geom_bar( stat="identity", position="fill", width = 0.5)+ coord_flip()+ scale_y_reverse()
Per3 <- Per3 + facet_grid(class~.)
Per3<-Per3+theme_void()+theme(strip.text.x = element_text(size = 8, angle = 90))+ scale_fill_manual(values=mypalette_1)
Per3



pdf(paste0(res_path, "Pie chart.pdf"))
Per1
Per3
dev.off()


# Additional Statistics
library(emmeans)

# Create an empty dataframe to store results
results_df <- data.frame(cluster2 = character(), p_value = numeric(),contrast=character(), pairwise_p_values=numeric())

# Create an empty list to store plots
plots_list <- list()

# Iterate over unique clusters
for (cluster in unique(Percentage_clus_2$cluster2)) {
  # Subset data for the current cluster
  cluster_data <- subset(Percentage_clus_2, cluster2 == cluster)
  
  # Calculate mean percentage for each mouse and subtract from original percentage
  normalized_data <- cluster_data %>%
    group_by(mouse) %>%
    mutate(normalized_percentage = (percentage - mean(percentage)) / sd(percentage)) %>%
    ungroup()
  
  # Fit model with normalized percentage values
  model <- lm(normalized_percentage ~ class, data = normalized_data)
  
  # Extract estimates and p-value from model summary
  model_summary <- summary(model)
  p_value <- anova(model)$`Pr(>F)`[1]
  # Conduct pairwise comparisons using emmeans
  pairwise_comp <- emmeans(model, pairwise ~ class, adjust = "tukey")
  
  # Extract pairwise comparisons and p-values
  pairwise_p_values <- summary(pairwise_comp)[["contrasts"]][["p.value"]]
  contrast<-summary(pairwise_comp)[["contrasts"]][["contrast"]]
  # Append cluster name, estimates, and p-value to results dataframe
  results_df <- rbind(results_df, data.frame(cluster2 = cluster, p_value = p_value, contrast=contrast,pairwise_p_values=pairwise_p_values))
  # Create boxplot of normalized percentages by class
  plot <- ggplot(normalized_data, aes(x = class, y = normalized_percentage)) +
    geom_boxplot(width=0.5) +geom_jitter(width = 0.3)+
    labs(title = paste("Cluster:", cluster)) +
    theme_bw()+theme(aspect.ratio = 1)
  # Store plot in the list
  plots_list[[cluster]] <- plot
  
}

# Display the results dataframe
print(results_df)
write.csv(results_df,paste0(res_path, "large_scale_per_cl.csv"))

pdf(paste0(res_path, "large_scale_per_cl.pdf"))
plots_list
dev.off()

####Differences based on environmental clusters

### plot differences of the values per cluster
master_class_sum<-as.data.frame(master_class_sum)
master_class_sum$cluster2<-as.character(master_class_sum$cluster2)
df1$mean_movement<-df1$movement
master_class_sum2<-left_join(master_class_sum, df1[,c("cluster2","mean_movement")])
master_class_sum2$cluster2 = factor(master_class_sum2$cluster2, levels=df1$cluster2)


#### try an anova between clusters based on distance to speed
p_class_speed <- ggplot(master_class_sum2,aes(x=as.factor(class) , y=speed, group=class,fill=class)) +
  geom_violin(alpha=1, outlier.colour = NA)+stat_summary(fun = median, geom = "crossbar",  size = 1, color = "grey25")+scale_fill_manual(values=c("cyan","gold","red"))+ggtitle("speed per env cluster")+theme_classic()+
  theme(aspect.ratio=1)+coord_cartesian(ylim = c(0,20))

p_class_speed


a_class_speed <- aov(speed ~ class, data=master_class_sum2) 
summary(a_class_speed)
TukeyHSD(a_class_speed)

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_class_speed), file=paste0(res_path,"/Tukey_speed_cytomap_cl.txt"))




#### try an anova between clusters based on distance to squared displacement
p_class_disp2 <- ggplot(master_class_sum2,aes(x=as.factor(class) , y=disp2, group=class,fill=class)) +
  geom_violin(alpha=1, outlier.colour = NA)+stat_summary(fun = median, geom = "crossbar",  size = 1, color = "grey25")+scale_fill_manual(values=c("cyan","gold","red"))+ggtitle("speed per env cluster")+theme_classic()+
  theme(aspect.ratio=1)+coord_cartesian(ylim = c(0,180))

p_class_disp2


a_class_disp2 <- aov(disp2 ~ class, data=master_class_sum2) 
summary(a_class_disp2)
TukeyHSD(a_class_disp2)

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_class_disp2), file=paste0(res_path,"/Tukey_disp2_cytomap_cl.txt"))



#### try an anova between clusters based on raw tumor movement
p_class_movement <- ggplot(master_class_sum2,aes(x=as.factor(class) , y=movement, group=class,fill=class)) +
  geom_violin(alpha=1, outlier.colour = NA)+stat_summary(fun = median, geom = "crossbar",  size = 1, color = "grey25")+scale_fill_manual(values=c("cyan","gold","red")) +ggtitle("movement direction per env cluster")+theme_classic()+
  theme(aspect.ratio=1)+scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(4, "Spectral")))+coord_cartesian(ylim = c(-10,17))

p_class_movement


a_class_movement <- aov(movement ~ class, data=master_class_sum2) 
summary(a_class_movement)
TukeyHSD(a_class_movement)

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_class_movement), file=paste0(res_path,"/Tukey_move_dir_cytomap_cl.txt"))





pdf(paste0(res_path, "environmental_cluster_stats.pdf"))
p_class_speed
p_class_disp2
p_class_movement
dev.off()




################################################################################
##################                                   ###########################
##################   Small-scale phenotyping Module  ###########################
##################                                   ###########################
################################################################################

# Create a directory for the results
folder <- "/content/results/small_scale_phenotyping_module"

if (file.exists(folder)) {
  cat("The '/content/results/small_scale_phenotyping_module' folder already exists")
} else {
  dir.create(folder)
}


###MULTIVARIATE
list_master_pca <- split(master_pca3[,-c(1)],master_pca3$Track2) ## split into list

setwd("C:/results")

if ( ! file.exists("matrix_distmat.rds")){
  ##parallel working
  # load parallel
  library(parallel)
  # create multi-process workers
  workers <- makeCluster(detectCores()-2)
  # load dtwclust in each one, and make them use 1 thread per worker
  invisible(clusterEvalQ(workers, {
    library(dtwclust)
    RcppParallel::setThreadOptions(1L)
  }))
  # register your workers, e.g. with doParallel
  require(doParallel)
  registerDoParallel(workers)
  
  ###MULTIVARIATE
  set.seed(123)
  distmat <- proxy::dist(list_master_pca, method = "dtw")
  saveRDS(distmat, file = "distmat.rds")
  
  matrix_distmat<-as.matrix(distmat)
  saveRDS(matrix_distmat, file = "matrix_distmat.rds")
  
  
}else{
  matrix_distmat <- readRDS("matrix_distmat.rds")
}


Track2<-as.numeric(names(list_master_pca))

library(dplyr)
library(scales)
library(pheatmap)
library(umap)
set.seed(123)

# Adjust parameters according to the characteristics of the experiment
umap_dist<- umap(matrix_distmat,n_components=2,input="dist",init = "random", 
                 n_neighbors=7,min_dist=0.5, spread=6,n_epochs=1000 , local_connectivity=1) 

umap_1 <- as.data.frame(umap_dist$`layout`) 

Track2_umap<-cbind(Track2,umap_1)
positiontype<-master_norm[,c("Track2", "position","class","mouse")]
positiontype<- positiontype[!duplicated(positiontype$Track2),]
umap_2 <- left_join(Track2_umap ,positiontype)

# Clustering
set.seed(123)
n_cluster=7
hc <- hclust(dist(umap_dist$`layout`, method = "euclidean"), method = "ward.D2")
hierarchical_clusters <- cutree(hc, k = n_cluster) 
umap_3 <- cbind(hierarchical_clusters, umap_2)
colnames(umap_3)[1]<- "cluster2"


### calculate average stats per track
master_class <- left_join(master_cor2 ,umap_3[c(1:4)], by = c("Track2" = "Track2"))
master_class<-master_class %>% 
  filter(cluster2 != 0)
master_class$cluster2<-as.numeric(master_class$cluster2)
master_class_sum<-master_class%>%group_by(position,class,mouse, Track2)%>%
  arrange(Time)%>%summarize(cluster2=mean(cluster2), V1=mean(V1), V2=mean(V2),
                            dist_tumor=mean(dist_tumor), dist_3_neigh=mean(dist_3_neigh),
                            dist_10_neigh=mean(dist_10_neigh), speed=mean(speed), disp2=mean(disp2), 
                            disp_d=mean(disp_d),disp_l=mean(disp_l), movement=last(dist_tumor), 
                            distance_to_tumor=last(dist_tumor_1),raw_dist_tumor=last(dist_tumor), 
                            Time=first(Time))%>%ungroup()


## Join the information on the SR101 and MG
master_class_sum <- left_join(master_class_sum,master_distance_MG[c("Track2","n_MG","min_MG")] ,by = c("Track2" = "Track2"))
master_class_sum <- left_join(master_class_sum,master_distance_SR101[c("Track2","n_SR101","min_SR101")] ,by = c("Track2" = "Track2"))
master_class_sum <- left_join(master_class_sum,BV_df_sum ,by = c("Track2" = "Track2"))

mid <- 0
master_class_sum$movement2<-ifelse(master_class_sum$movement>0,"away from tumor", ifelse(master_class_sum$movement==0, "no movement", "towards the tumor"))
master_class_sum<-master_class_sum[!is.na(master_class_sum$cluster2),]

saveRDS(master_class_sum,"master_class_sum.rds")
write.csv(master_class_sum,"master_class_sum.csv")


###Create heatmap
sum_all<-master_class_sum%>%group_by(cluster2)%>%
  summarise(displacement2 = median(disp2), 
            speed = median(speed),displacement_delta=median(disp_d),displacement_length=median(disp_l),max_speed=quantile(speed, 0.75), persistance=(displacement_length/displacement_delta),nearest_10_cell= median(dist_10_neigh),dist_tumor_core=median(distance_to_tumor),
            n_SR101=mean(n_SR101, na.rm = T),
            min_SR101=mean(min_SR101, na.rm = T),n_MG=mean(n_MG, na.rm = T),min_MG=mean(min_MG, na.rm = T),min_BV_dist=mean(BV_min, na.rm = T),sd_BV_dist=mean(BV_sd, na.rm = T),mean_BV_dist=mean(BV_mean, na.rm = T),
            movement=median(movement) )

Cluster_movement<-sum_all[,c(1,length(names(sum_all)))]
Cluster_movement <-Cluster_movement[order(Cluster_movement$movement),]

stdize = function(x, ...) {(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
sum_all[-c(1)] <- lapply(sum_all[-c(1)], stdize, na.rm = T)

library(viridis)
library(plotly)

df<-sum_all[,-c(1)]
df<-df[order(match(rownames(df), Cluster_movement$cluster2)), , drop = FALSE]%>%ungroup()
df<-as.data.frame(df)
rownames(df)<-Cluster_movement$cluster2

heat_m<-pheatmap(df, clustering_method = "complete", cluster_cols=F,cluster_rows=F,cellwidth=20, cellheight=20,
                 treeheight_row = 0, treeheight_col = 0,fontsize = 8,na_col = "grey50", main=" All features",
                 angle_col = 90,color = colorRampPalette(c("deepskyblue1", "grey95", "deeppink2"))(50))

df_other_features<-df[,!colnames(df)%in%c("displacement2", "speed", "max_speed", "displacement_delta" ,  "displacement_length", "persistance", "movement")]
rownames(df_other_features)<-rownames(df)

pdf(paste0(res_path, "cluster_heatmap_dynamic_features_2.pdf"))

heat_m

pheatmap(df_other_features, clustering_method = "complete", cluster_cols=F,cluster_rows=F,cellwidth=20, cellheight=20,
         treeheight_row = 0, treeheight_col = 0,fontsize = 8,na_col = "grey50",main="Small-scale TME features",
         angle_col = 90,color = colorRampPalette(c("deepskyblue1", "grey95", "deeppink2"))(50))

dev.off()


## plot clusters giving them a color according to direction
library(scales)
# Convert df to a data frame
# Add cluster2 as a column
df1<-as.data.frame(df)

df1 <- mutate(df1, cluster2 = row.names(df))

# Order df1 by movement and reorder cluster2 accordingly
df1 <- df1 %>%
  arrange(movement) %>%
  mutate(cluster2 = factor(cluster2, levels = unique(cluster2)))
##set first order clusters based on the heatmap:



# Create a palette based on the movement values
pdf(paste0(res_path, "UMAP_cluster_and_other_features.pdf"))


ggplot(master_class_sum, mapping=aes(x=V1, y=V2, color=n_SR101)) +  
  geom_point(size=2, alpha=0.6) +
  xlab("") + ylab("") +
  ggtitle("UMAP number of SR101 in 30um radius distribution") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)+
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(4, "Spectral")),na.value="NA")


ggplot(master_class_sum, mapping=aes(x=V1, y=V2, color=min_SR101)) +  
  geom_point(size=2, alpha=0.6) +
  xlab("") + ylab("") +
  ggtitle("UMAP min distance to SR101 in 30um radius distribution") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)+
  scale_color_gradientn(colours = RColorBrewer::brewer.pal(4, "Spectral"),na.value="NA", trans="pseudo_log")


ggplot(master_class_sum, mapping=aes(x=V1, y=V2, color=min_MG)) +  
  geom_point(size=2, alpha=0.6) +
  xlab("") + ylab("") +
  ggtitle("UMAP min distance to MG in 30um radius distribution") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)+
  scale_color_gradientn(colours = RColorBrewer::brewer.pal(4, "Spectral"),na.value="NA", trans="pseudo_log")


ggplot(master_class_sum, mapping=aes(x=V1, y=V2, color=n_MG)) +  
  geom_point(size=2, alpha=0.6) +
  xlab("") + ylab("") +
  ggtitle("UMAP number of MG in 30um radius distribution") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)+
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(4, "Spectral")),na.value="NA")



ggplot(master_class_sum, mapping=aes(x=V1, y=V2, color=BV_mean)) +  
  geom_point(size=2, alpha=0.6) +
  xlab("") + ylab("") +
  ggtitle("UMAP mean distance to Blood vessels") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)+
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(4, "Spectral")),na.value="NA", trans="pseudo_log")



ggplot(master_class_sum, mapping=aes(x=V1, y=V2, color=BV_min)) +  
  geom_point(size=2, alpha=0.6) +
  xlab("") + ylab("") +
  ggtitle("UMAP mean contact with Blood vessels") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)+
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(4, "Spectral")),na.value="NA")



ggplot(master_class_sum, mapping=aes(x=V1, y=V2, color=dist_10_neigh)) +  
  geom_point(size=2, alpha=0.6) +
  xlab("") + ylab("") +
  ggtitle("UMAP big neighbours distribution") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)+
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(4, "Spectral")),trans="pseudo_log")


ggplot(master_class_sum, mapping=aes(x=V1, y=V2, color=dist_3_neigh)) +  
  geom_point(size=2, alpha=0.6) +
  xlab("") + ylab("") +
  ggtitle("UMAP close neighbours distribution") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)+
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(4, "Spectral")),trans="pseudo_log")

dev.off()



## Differences based on clusters
### plot differences of the values per cluster
master_class_sum<-as.data.frame(master_class_sum)
master_class_sum$cluster2<-as.character(master_class_sum$cluster2)
df1$mean_movement<-df1$movement
master_class_sum2<-left_join(master_class_sum, df1[,c("cluster2","mean_movement")])
master_class_sum2$cluster2 = factor(master_class_sum2$cluster2, levels=df1$cluster2)

# Differences per behavioral cluster in small-scale TME features

#### try an anova between clusters based on distance to dist_3_neigh
p_dist3_neigth <- ggplot(master_class_sum2,aes(x=as.factor(cluster2) , y=dist_3_neigh, group=cluster2, fill=cluster2)) +geom_jitter(width=0.2, alpha=0.5) +
  geom_boxplot(alpha=0.5, outlier.colour = NA)+ggtitle("dist 3 neigh per cluster")+theme_classic()+
  theme(aspect.ratio=0.7)+ scale_fill_manual(values=mypalette_1)+coord_cartesian(ylim = c(10,70))

p_dist3_neigth

a_dist_3_neigh <- aov(dist_3_neigh ~ cluster2, data=master_class_sum2) 
summary(a_dist_3_neigh)
TukeyHSD(a_dist_3_neigh)

capture.output(TukeyHSD(a_dist_3_neigh), file=paste0(res_path,"Tukey_dist_3_neigh.txt"))
capture.output(summary(a_dist_3_neigh), file=paste0(res_path,"aov_dist_3_neigh.txt"))

#### try an anova between clusters based on distance to dist_10_neigh
p_dist_10_neigh <- ggplot(master_class_sum2,aes(x=as.factor(cluster2) , y=dist_10_neigh, group=cluster2,fill=cluster2)) +geom_jitter(width=0.2, alpha=0.5) +
  geom_boxplot( outlier.colour = NA, alpha=0.5) +ggtitle("dist 10 neigh per cluster")+theme_classic()+
  theme(aspect.ratio=0.7)+ scale_fill_manual(values=mypalette_1)+coord_cartesian(ylim = c(20,90))

p_dist_10_neigh

a_dist_10_neigh <- aov(dist_10_neigh ~ cluster2, data=master_class_sum2) 
summary(a_dist_10_neigh)
TukeyHSD(a_dist_10_neigh)

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_dist_10_neigh), file=paste0(res_path,"Tukey_dist_10_neigh.txt"))
capture.output(summary(a_dist_10_neigh), file=paste0(res_path,"aov_dist_10_neigh.txt"))


#### try an anova between clusters based on dist to MG
p_min_MG <- ggplot(subset(master_class_sum2, !is.na(min_MG)) ,aes(x=as.factor(cluster2) , y=min_MG, group=cluster2, fill=cluster2))  +geom_jitter(width=0.3, alpha=0.5) +
  geom_boxplot( outlier.colour = NA, alpha=0.5) +ggtitle("distance to closest MG")+theme_classic()+
  theme(aspect.ratio=0.7)+ scale_fill_manual(values=mypalette_1)+coord_cartesian(ylim = c(2,45))

p_min_MG


a_min_MG <- aov(min_MG ~ cluster2, data= subset(master_class_sum2, !is.na(min_MG))) 
summary(a_min_MG)
TukeyHSD(a_min_MG)


### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_min_MG),  file=paste0(res_path,"Tukey_dist_min_MG.txt"))
capture.output(summary(a_min_MG),  file=paste0(res_path,"aov_dist_min_MG.txt"))


#### try an anova between clusters based on dist to SR101
p_min_SR101 <- ggplot(subset(master_class_sum2, !is.na(min_SR101)) ,aes(x=as.factor(cluster2) , y=min_SR101, group=cluster2, fill=cluster2))  +geom_jitter(width=0.3, alpha=0.5) +
  geom_boxplot(outlier.colour = NA, alpha=0.5) +ggtitle("distance to closest SR101")+theme_classic()+
  theme(aspect.ratio=0.7)+ scale_fill_manual(values=mypalette_1)+coord_cartesian(ylim = c(0,50))

p_min_SR101


a_min_SR101 <- aov(min_SR101 ~ cluster2, data= subset(master_class_sum2, !is.na(min_SR101))) 
summary(a_min_SR101)
TukeyHSD(a_min_SR101)



### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_min_SR101), file=paste0(res_path,"Tukey_dist_min_SR101.txt"))
capture.output(summary(a_min_SR101), file=paste0(res_path,"aov_dist_min_SR101.txt"))



#### try an anova between clusters based on n of  MG
p_n_MG <- ggplot(subset(master_class_sum2, !is.na(n_MG)) ,aes(x=as.factor(cluster2) , y=n_MG, group=cluster2, fill=cluster2))  +geom_jitter(width=0.3, alpha=0.5) +
  geom_boxplot( outlier.colour = NA, alpha=0.5) +ggtitle("n MG")+theme_classic()+
  theme(aspect.ratio=0.7)+ scale_fill_manual(values=mypalette_1)+coord_cartesian(ylim = c(0,15))

p_n_MG


a_n_MG <- aov(n_MG ~ cluster2, data= subset(master_class_sum2, !is.na(min_MG))) 
summary(a_n_MG)
TukeyHSD(a_n_MG)
### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_n_MG), file=paste0(res_path,"Tukey_dist_n_MG.txt") )
capture.output(summary(a_n_MG), file=paste0(res_path,"aov_dist_n_MG.txt") )




#### try an anova between clusters based on n of  SR101
p_nSR101 <- ggplot(subset(master_class_sum2, !is.na(n_SR101)) ,aes(x=as.factor(cluster2) , y=n_SR101, group=cluster2, fill=cluster2))  +geom_jitter(width=0.3, alpha=0.5) +
  geom_boxplot(outlier.colour = NA, alpha=0.5) +ggtitle("n SR101")+theme_classic()+
  theme(aspect.ratio=0.7)+ scale_fill_manual(values=mypalette_1)+coord_cartesian(ylim = c(0,20))

p_nSR101



a_nSR101 <- aov(n_SR101 ~ cluster2, data= subset(master_class_sum2, !is.na(n_SR101))) 
summary(a_nSR101)
TukeyHSD(a_nSR101)
### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_nSR101), file=paste0(res_path,"Tukey_dist_n_SR101.txt") )
capture.output(summary(a_nSR101), file=paste0(res_path,"aov_dist_n_SR101.txt") )



#### try an anova between clusters based on mean distance to BV

p_BV_mean <- ggplot(subset(master_class_sum2, !is.na(BV_mean)) ,aes(x=as.factor(cluster2) , y=BV_mean, group=cluster2, fill=cluster2))  +geom_jitter(width=0.3, alpha=0.5) +
  geom_boxplot(outlier.colour = NA, alpha=0.5) +ggtitle("BV_mean")+theme_classic()+
  theme(aspect.ratio=0.7)+ scale_fill_manual(values=mypalette_1)+coord_cartesian(ylim = c(0,20))

p_BV_mean



a_BV_mean <- aov(BV_mean ~ cluster2, data= subset(master_class_sum2, !is.na(BV_mean))) 
summary(a_BV_mean)
TukeyHSD(a_BV_mean)

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_BV_mean),file=paste0(res_path,"Tukey_dist_mean_BV.txt"))
capture.output(summary(a_BV_mean),file=paste0(res_path,"aov_dist_mean_BV.txt"))


#### try an anova between clusters based on min distance to BV

p_BV_min <- ggplot(subset(master_class_sum2, !is.na(BV_min)) ,aes(x=as.factor(cluster2) , y=BV_min, group=cluster2, fill=cluster2))  +geom_jitter(width=0.2, alpha=0.5) +
  geom_boxplot( outlier.colour = NA, alpha=0.5) +ggtitle("BV_min")+theme_classic()+
  theme(aspect.ratio=0.7)+ scale_fill_manual(values=mypalette_1)+coord_cartesian(ylim = c(0,20))

p_BV_min



a_BV_min <- aov(BV_min ~ cluster2, data= subset(master_class_sum2, !is.na(BV_min))) 
summary(a_BV_min)
TukeyHSD(a_BV_min)


### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_BV_min), file=paste0(res_path,"Tukey_dist_min_BV.txt"))
capture.output(summary(a_BV_min), file=paste0(res_path,"aov_dist_min_BV.txt"))



#### try an anova between clusters based on contact to BV

p_BV_contact <- ggplot(subset(master_class_sum2, !is.na(BV_min)) ,aes(x=as.factor(cluster2) , y=BV_contact, group=cluster2, fill=cluster2))  +geom_jitter(width=0.3, alpha=0.5) +
  geom_boxplot(outlier.colour = NA, alpha=0.5) +ggtitle("BV_contact")+theme_classic()+ scale_fill_manual(values=mypalette_1)+
  theme(aspect.ratio=0.7)

p_BV_contact



a_BV_contact <- aov(BV_contact ~ cluster2, data= subset(master_class_sum2, !is.na(BV_contact))) 
summary(a_BV_contact)
TukeyHSD(a_BV_contact)

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_BV_contact), file=paste0(res_path,"Tukey_dist_contact_BV.txt"))
capture.output(summary(a_BV_contact), file=paste0(res_path,"aov_dist_contact_BV.txt"))



#### try an anova between clusters based on sd distance to BV

p_BV_sd <- ggplot(subset(master_class_sum2, !is.na(BV_min)) ,aes(x=as.factor(cluster2) , y=BV_sd, group=cluster2, fill=cluster2))  +geom_jitter(width=0.3, alpha=0.5) +
  geom_boxplot( outlier.colour = NA, alpha=0.5) +ggtitle("Standart deviation distance")+theme_classic()+
  theme(aspect.ratio=0.7)+ scale_fill_manual(values=mypalette_1)+coord_cartesian(ylim = c(0,1.7))

p_BV_sd



a_BV_sd <- aov(BV_sd ~ cluster2, data= subset(master_class_sum2, !is.na(BV_min))) 
TukeyHSD(a_BV_sd)


### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_BV_sd), file=paste0(res_path,"Tukey_dist_sd_BV.txt"))
capture.output(summary(a_BV_sd), file=paste0(res_path,"aov_dist_sd_BV.txt"))



pdf(paste0(res_path,"per_cluster_features_comparison.pdf" )) ##adjust path

p_dist3_neigth
p_dist_10_neigh
p_min_MG
p_min_SR101
p_n_MG
p_nSR101
p_BV_mean
p_BV_min
p_BV_contact
p_BV_sd

dev.off()


####Differences based on environmental clusters with small-scale TME features

#### try an anova between clusters based on distance to dist_10_neigh
p_class_dist_10neigh <- ggplot(master_class_sum2,aes(x=as.factor(class) , y=dist_10_neigh, fill=class)) +
  geom_violin(alpha=1, outlier.colour = NA)+stat_summary(fun = median, geom = "crossbar",  size = 1, color = "grey25")+scale_fill_manual(values=c("cyan","gold","red"))  +ggtitle("dist 10 neigh per cluster")+theme_classic()+
  theme(aspect.ratio=1)+coord_cartesian(ylim = c(20,70))

p_class_dist_10neigh
a_class_dist_10neigh <- aov(dist_10_neigh ~ class, data=master_class_sum2) 
summary(a_class_dist_10neigh)
TukeyHSD(a_class_dist_10neigh)

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_class_dist_10neigh), file=paste0(res_path,"/Tukeydist10neigth_cytomap_cl.txt"))


#### try an anova between clusters based on dist to MG
p_class_minMG <- ggplot(subset(master_class_sum2, !is.na(min_MG)) ,aes(x=as.factor(class) , y=min_MG, fill=class))  +
  geom_violin(alpha=1, outlier.colour = NA)+stat_summary(fun = median, geom = "crossbar",  size = 1, color = "grey25")+scale_fill_manual(values=c("cyan","gold","red"))  +ggtitle("distance to closest MG")+theme_classic()+
  theme(aspect.ratio=1)+coord_cartesian(ylim = c(0,50))

p_class_minMG


a_class_minMG <- aov(min_MG ~ class, data= subset(master_class_sum2, !is.na(min_MG))) 
summary(a_class_minMG)
TukeyHSD(a_class_minMG)

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_class_minMG), file=paste0(res_path,"/Tukey_minMG_cytomap_cl.txt"))



#### try an anova between clusters based on dist to SR101
p_classmin_SR101 <- ggplot(subset(master_class_sum2, !is.na(min_SR101)) ,aes(x=as.factor(class) , y=min_SR101, fill=class))  +
  geom_violin(alpha=1, outlier.colour = NA) +stat_summary(fun = median, geom = "crossbar",  size = 1, color = "grey25")+scale_fill_manual(values=c("cyan","gold","red"))+ggtitle("distance to closest SR101")+theme_classic()+
  theme(aspect.ratio=1)+coord_cartesian(ylim = c(0,50))

p_classmin_SR101 


a_classmin_SR101  <- aov(min_SR101 ~ class, data= subset(master_class_sum2, !is.na(min_SR101))) 
summary(a_classmin_SR101 )
TukeyHSD(a_classmin_SR101 )

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_classmin_SR101 ), file=paste0(res_path,"/Tukey_minSR101_cytomap_cl.txt"))



#### try an anova between clusters based on n of  MG
p_class_n_MG <- ggplot(subset(master_class_sum2, !is.na(min_MG)),aes(x=as.factor(class) , y=n_MG, fill=class)) +
  geom_violin(alpha=1, outlier.colour = NA) +stat_summary(fun = median, geom = "crossbar",  size = 1, color = "grey25")+scale_fill_manual(values=c("cyan","gold","red"))+ggtitle("n of MG per cluster")+theme_classic()+
  theme(aspect.ratio=1)

p_class_n_MG

a_class_n_MG <- aov(n_MG ~ class, data= subset(master_class_sum2, !is.na(min_MG))) 
summary(a_class_n_MG)
TukeyHSD(a_class_n_MG)

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_class_n_MG), file=paste0(res_path,"/Tukey_nMG_cytomap_cl.txt"))



#### try an anova between clusters based on n of  SR101

p_classn_SR101 <- ggplot(subset(master_class_sum2, !is.na(min_SR101)),aes(x=as.factor(class) , y=n_SR101, fill=class)) +
  geom_violin(alpha=1, outlier.colour = NA)+stat_summary(fun = median, geom = "crossbar",  size = 1, color = "grey25")+scale_fill_manual(values=c("cyan","gold","red"))+ggtitle("n SR101 per cluster")+theme_classic()+
  theme(aspect.ratio=1)
p_classn_SR101

a_classn_SR101 <- aov(n_SR101 ~ class, data= subset(master_class_sum2, !is.na(min_SR101))) 
summary(a_classn_SR101)
TukeyHSD(a_classn_SR101)

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_classn_SR101), file=paste0(res_path,"/Tukey_nSR101_cytomap_cl.txt"))


#### try an anova between clusters based on BV distance min
p_class_BV_min <- ggplot(subset(master_class_sum2, !is.na(BV_min)) ,aes(x=as.factor(class) , y=BV_min, fill=class))  +
  geom_violin(alpha=1, outlier.colour = NA)+stat_summary(fun = median, geom = "crossbar",  size = 1, color = "grey25")+ scale_fill_manual(values=c("cyan","gold","red"))  +ggtitle("min distance to BV")+theme_classic()+
  theme(aspect.ratio=1)+coord_cartesian(ylim = c(0,18))

p_class_BV_min


a_class_BV_min <- aov(BV_min ~ class, data= subset(master_class_sum2, !is.na(BV_min))) 
summary(a_class_BV_min)
TukeyHSD(a_class_BV_min)

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_class_BV_min), file=paste0(res_path,"/Tukey_minBV_cytomap_cl.txt"))


#### try an anova between clusters based on BV distance mean
p_class_BV_mean <- ggplot(subset(master_class_sum2, !is.na(BV_mean)) ,aes(x=as.factor(class) , y=BV_mean, fill=class))  +
  geom_violin(alpha=1, outlier.colour = NA)+stat_summary(fun = median, geom = "crossbar",  size = 1, color = "grey25")+ scale_fill_manual(values=c("cyan","gold","red"))  +ggtitle("mean distance to BV")+theme_classic()+
  theme(aspect.ratio=1)+coord_cartesian(ylim = c(0,18))

p_class_BV_mean


a_class_BV_mean <- aov(BV_mean ~ class, data= subset(master_class_sum2, !is.na(BV_mean))) 
summary(a_class_BV_mean)
TukeyHSD(a_class_BV_mean)

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a_class_BV_mean), file=paste0(res_path,"/Tukey_meanBV_cytomap_cl.txt"))


##relation between position relative to the tumor and amount of MG:

pcor_MG_position <- ggplot(subset(master_class_sum2, !is.na(min_MG)),aes(x=distance_to_tumor , y=n_MG)) +geom_jitter()+
  geom_point(alpha=0.5,stat = "summary") +geom_smooth(method = "lm")+ggtitle("scatter plot movement vs min_MG")+theme_classic()+
  theme(aspect.ratio=1)

pcor_MG_position

cor_18<-cor.test(master_class_sum2$n_MG,master_class_sum2$distance_to_tumor,  method = "pearson", use = "pairwise.complete.obs")
cor_18

### save as txt the output of this comparison. Change the path for your own data
capture.output(cor_18, file=paste0(res_path,"/correlation_nMG_to_location_at_tumor_border.txt"))


pdf(paste0(res_path, "environmental_cluster_stats.pdf"))
p_class_dist_10neigh
p_class_minMG
p_classmin_SR101
p_class_n_MG
p_classn_SR101
p_class_BV_min
p_class_BV_mean
pcor_MG_position 
dev.off()
