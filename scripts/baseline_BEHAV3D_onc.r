

library(plyr)
library(readr)
library(dplyr)
library(ggplot2)
setwd("D:/R/Projects/IVM")
if ( ! file.exists("master_cor.rds")){
  ## function to import the stats of the tracked tumor cells data from csv files extracted by Imaris
  read_plus <- function(flnm) {
    read_csv(flnm, skip = 3, col_types = cols()) %>% 
      mutate(filename = flnm)
  }
  ##import stats of interest
  ## set directory where csv files are located  
  
  
  working_directory <- "E:/DATA/Intravital DATA/MA18_MA22/common"
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
  
  category <- as.factor(dipl_length_csv$filename)
  ranks <- rank(-table(category), ties.method="first")
  ranks <- as.data.frame(ranks)
  ranks$filename <- row.names(ranks)
  dipl_length_csv <- left_join(dipl_length_csv, ranks) 
  dipl_length_csv$ID2 <- with(dipl_length_csv, interaction(ID, ranks))
  
  category <- as.factor(dipl_delta_csv$filename)
  ranks <- rank(-table(category), ties.method="first")
  ranks <- as.data.frame(ranks)
  ranks$filename <- row.names(ranks)
  dipl_delta_csv <- left_join(dipl_delta_csv, ranks) 
  dipl_delta_csv$ID2 <- with(dipl_delta_csv, interaction(ID, ranks))
  
  category <- as.factor(disp_csv$filename)
  ranks <- rank(-table(category), ties.method="first")
  ranks <- as.data.frame(ranks)
  ranks$filename <- row.names(ranks)
  disp_csv <- left_join(disp_csv, ranks) 
  disp_csv$ID2 <- with(disp_csv, interaction(ID, ranks))
  
  
  ###Import BLOOD VESSELS stats (does not exist for all tracked positions so needs to be processed separtely.)
  working_directory2 <- "E:/DATA/Intravital DATA/MA18_MA22/BV_stats"
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
  ## rename mouse accroding to filename:
  master$mouse<-master$filename
  master$mouse <- gsub(".*F11.*", "2430F11", master$mouse, perl=TRUE)
  master$mouse <- gsub(".*F18.*", "2430F18", master$mouse, perl=TRUE)
  master$mouse <- gsub(".*F12.*", "2430F12", master$mouse, perl=TRUE)
  master$mouse <- gsub(".*F16.*", "2430F16", master$mouse, perl=TRUE)
  master$mouse <- gsub(".*F17.*", "2430F17", master$mouse, perl=TRUE)
  master$mouse <- gsub(".*F13.*", "2430F13", master$mouse, perl=TRUE)
  
  master$position <- master$filename
  master$position <- gsub(".*CL1_1.*", "CL1_1", master$position, perl=TRUE)
  master$position <- gsub(".*CL1_2.*", "CL1_2", master$position, perl=TRUE)
  master$position <- gsub(".*CL2_1.*", "CL2_1", master$position, perl=TRUE)
  master$position <- gsub(".*CL2_2.*", "CL2_2", master$position, perl=TRUE)
  master$position <- gsub(".*CL3_1.*", "CL3_1", master$position, perl=TRUE)
  master$position <- gsub(".*CL3_2.*", "CL3_2", master$position, perl=TRUE)
  master$position<-with(master, interaction(mouse, position))
  
  master$class <- master$position
  master$class <- gsub(".*CL3.*", "CL3", master$class, perl=TRUE)
  master$class <- gsub(".*CL2.*", "CL2", master$class, perl=TRUE)
  master$class <- gsub(".*CL1.*", "CL1", master$class, perl=TRUE)
  
  
  master$filename<-NULL
  ## Round the time
  master$Time<-round(master$Time, digits=2)
  
  
  
  ###Process blood vessels data:
  BV_df<-BV_df[c(1,6,7,8,9,11)]
  colnames(BV_df) <- c("dist_BV","Time","Tracked","TrackID", "ID","filename")
  ## rename mouse accroding to filename:
  BV_df$mouse<-BV_df$filename
  BV_df$mouse <- gsub(".*F11.*", "2430F11", BV_df$mouse, perl=TRUE)
  BV_df$mouse <- gsub(".*F16.*", "2430F16", BV_df$mouse, perl=TRUE)
  BV_df$mouse <- gsub(".*F13.*", "2430F13", BV_df$mouse, perl=TRUE)
  
  BV_df$position <- BV_df$filename
  BV_df$position <- gsub(".*CL1_1.*", "CL1_1", BV_df$position, perl=TRUE)
  BV_df$position <- gsub(".*CL1_2.*", "CL1_2", BV_df$position, perl=TRUE)
  BV_df$position <- gsub(".*CL2_1.*", "CL2_1", BV_df$position, perl=TRUE)
  BV_df$position <- gsub(".*CL2_2.*", "CL2_2", BV_df$position, perl=TRUE)
  BV_df$position <- gsub(".*CL3_1.*", "CL3_1", BV_df$position, perl=TRUE)
  BV_df$position <- gsub(".*CL3_2.*", "CL3_2", BV_df$position, perl=TRUE)
  BV_df$position<-with(BV_df, interaction(mouse, position))
  
  BV_df$class <- BV_df$position
  BV_df$class <- gsub(".*CL3.*", "CL3", BV_df$class, perl=TRUE)
  BV_df$class <- gsub(".*CL2.*", "CL2", BV_df$class, perl=TRUE)
  BV_df$class <- gsub(".*CL1.*", "CL1", BV_df$class, perl=TRUE)
  
  
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
  List_bv<-list()
  test<-BV_df%>%arrange(Time)
  #measure the min and max number of timepoints
  number_timepoints<-test%>%group_by(Track2, position)%>%summarise(n=n())%>%group_by(position)%>%summarise( max_n=max(n))
  for(i in seq(from=1, to=max(number_timepoints$max_n), by=3)){  ## by needs to be a smaller number, otherwise we miss to many events. Perhaps 1 would be too heavy but we can try
    j=i+7  ### define here the length of the subtracks (here 30 timepoints)
    test3<-test%>%group_by(Track2)%>%slice(i:j)  ### slice the sequence in the definded range
    test3$iteration<-i ## create a variable for each subsequnce
    test3$subsequence<-with(test3, interaction(iteration,Track2,  sep=""))  ## create subsequence names
    List_bv[[length(List_bv)+1]] <-test3
    rm(test3)
  }
  
  BV_df <- do.call("rbind", List_bv)%>%ungroup()
  BV_df<-data.frame(BV_df)
  BV_df$subsequence<-as.numeric(as.character(BV_df$subsequence))
  
  
  BV_df_sum<-BV_df%>%group_by(subsequence, Track2, position)%>%summarise(BV_sd=sd(dist_BV)/sqrt(length(dist_BV)),BV_mean=mean(dist_BV), BV_min=min(dist_BV), BV_max=max(dist_BV), BV_contact=mean(BV_contact))
  
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
  
  working_directory <- "E:/DATA/Intravital DATA/MA18_MA22/SR101/common"
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
  master_SR101$mouse<-master_SR101$filename
  master_SR101$mouse <- gsub(".*F12.*", "2430F12", master_SR101$mouse, perl=TRUE)
  master_SR101$mouse <- gsub(".*F16.*", "2430F16", master_SR101$mouse, perl=TRUE)
  master_SR101$mouse <- gsub(".*F13.*", "2430F13", master_SR101$mouse, perl=TRUE)
  
  master_SR101$position <- master_SR101$filename
  master_SR101$position <- gsub(".*CL1_1.*", "CL1_1", master_SR101$position, perl=TRUE)
  master_SR101$position <- gsub(".*CL1_2.*", "CL1_2", master_SR101$position, perl=TRUE)
  master_SR101$position <- gsub(".*CL2_1.*", "CL2_1", master_SR101$position, perl=TRUE)
  master_SR101$position <- gsub(".*CL2_2.*", "CL2_2", master_SR101$position, perl=TRUE)
  master_SR101$position <- gsub(".*CL3_1.*", "CL3_1", master_SR101$position, perl=TRUE)
  master_SR101$position <- gsub(".*CL3_2.*", "CL3_2", master_SR101$position, perl=TRUE)
  master_SR101$position<-with(master_SR101, interaction(mouse, position))
  
  master_SR101$class <- master_SR101$position
  master_SR101$class <- gsub(".*CL3.*", "CL3", master_SR101$class, perl=TRUE)
  master_SR101$class <- gsub(".*CL2.*", "CL2", master_SR101$class, perl=TRUE)
  master_SR101$class <- gsub(".*CL1.*", "CL1", master_SR101$class, perl=TRUE)
  
  
  master_SR101$filename<-NULL
  
  ## The imported Time is not correct for F12 and F16, convert to hours:
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
  
  working_directory <- "E:/DATA/Intravital DATA/MA18_MA22/MG_cells/common_MG"
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
  master_MG$mouse<-master_MG$filename
  master_MG$mouse <- gsub(".*F11.*", "2430F11", master_MG$mouse, perl=TRUE)
  master_MG$mouse <- gsub(".*F18.*", "2430F18", master_MG$mouse, perl=TRUE)
  master_MG$mouse <- gsub(".*F17.*", "2430F17", master_MG$mouse, perl=TRUE)
  master_MG$position <- master_MG$filename
  master_MG$position <- gsub(".*CL1_1.*", "CL1_1", master_MG$position, perl=TRUE)
  master_MG$position <- gsub(".*CL1_2.*", "CL1_2", master_MG$position, perl=TRUE)
  master_MG$position <- gsub(".*CL2_1.*", "CL2_1", master_MG$position, perl=TRUE)
  master_MG$position <- gsub(".*CL2_2.*", "CL2_2", master_MG$position, perl=TRUE)
  master_MG$position <- gsub(".*CL3_1.*", "CL3_1", master_MG$position, perl=TRUE)
  master_MG$position <- gsub(".*CL3_2.*", "CL3_2", master_MG$position, perl=TRUE)
  master_MG$position<-with(master_MG, interaction(mouse, position))
  
  master_MG$class <- master_MG$position
  master_MG$class <- gsub(".*CL3.*", "CL3", master_MG$class, perl=TRUE)
  master_MG$class <- gsub(".*CL2.*", "CL2", master_MG$class, perl=TRUE)
  master_MG$class <- gsub(".*CL1.*", "CL1", master_MG$class, perl=TRUE)
  
  
  master_MG$filename<-NULL
  
  ## The imported Time is not correct for F12 and F16, convert to hours:
  master_MG$Time<-round(master_MG$Time, digits=2)
  
  
  
  
  
  library(dplyr)
  library(ggplot2)
  
  ### test that tracks imported are fine:
  master_test<-master%>%group_by(position)%>%mutate(x=x-min(x), y=y-min(y))
  g<-ggplot(master)+geom_path(aes(x, y, group=TrackID, colour=dist_tumor), arrow = arrow(ends = "last", type = "closed",length = unit(0.005, "inches")))+theme_bw()
  g<-g + theme(legend.position = "none")+facet_wrap(.~position, scales = "free")
  g
  
  ## For calculate interactions between cells 
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
  
  ##To keep all the information we create subtracks of longuer tracks, this way we can compare changes in behavior as well
  ### generate random subsequences of the same length
  List_control<-list()
  test<-master_distance%>%arrange(Time)
  #measure the min and max number of timepoints
  number_timepoints<-test%>%group_by(Track2, position)%>%summarise(n=n())%>%group_by(position)%>%summarise( max_n=max(n))
  for(i in seq(from=1, to=max(number_timepoints$max_n), by=3)){  ## by needs to be a smaller number, otherwise we miss to many events. Perhaps 1 would be too heavy but we can try
    j=i+7  ### define here the length of the subtracks (here 30 timepoints)
    test3<-test%>%group_by(Track2)%>%slice(i:j)  ### slice the sequence in the definded range
    test3$iteration<-i ## create a variable for each subsequnce
    test3$subsequence<-with(test3, interaction(iteration,Track2,  sep=""))  ## create subsequence names
    List_control[[length(List_control)+1]] <-test3
    rm(test3)
  }
  
  master_distance_subtracks <- do.call("rbind", List_control)%>%ungroup()
  master_distance_subtracks<-data.frame(master_distance_subtracks)
  master_distance_subtracks$subsequence<-as.numeric(as.character(master_distance_subtracks$subsequence))

  #### For the SR101 and MG datasets separately find what is the cell distance to these cells
  ##For SR101
  master_distance_SR101<-subset(master_distance_subtracks, mouse %in% unique(master_SR101$mouse))
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
      cells_radius<-nn2(data=master_SR101_pos[c("x","y","z")], query = master_pos_t[c("x","y","z")],k=30, treetype = 'bd',searchtype =  "radius", radius = 40)
      cells_min<-nn2(data=master_SR101_pos[c("x","y","z")], query = master_pos_t[c("x","y","z")],k=1, treetype = 'bd',searchtype =  "standard")
      
      dist_table = data.frame(master_pos_t[c("Time","Track2","mouse","position","class", "subsequence")],cells_min[["nn.dists"]],cells_radius[["nn.dists"]])
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
  master_distance_SR101<-master_dist_SR101%>%group_by(Track2,subsequence, mouse, position, class)%>%summarize(n_SR101=max(n_SR101), min_SR101 =min(min_SR101))
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
  master_distance_MG<-subset(master_distance_subtracks, mouse %in% unique(master_MG$mouse))
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
      cells_radius<-nn2(data=master_MG_pos[c("x","y","z")], query = master_pos_t[c("x","y","z")],k=30, treetype = 'bd',searchtype =  "radius", radius = 50)
      cells_min<-nn2(data=master_MG_pos[c("x","y","z")], query = master_pos_t[c("x","y","z")],k=1, treetype = 'bd',searchtype =  "standard")
      
      dist_table = data.frame(master_pos_t[c("Time","Track2","mouse","position","class", "subsequence")],cells_min[["nn.dists"]],cells_radius[["nn.dists"]])
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
  master_distance_MG<-master_dist_MG%>%group_by(Track2, mouse, position, class, subsequence)%>%summarize(n_MG=max(n_MG), min_MG =min(min_MG))
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
  max(master_distance_subtracks$Time)
  
  empty_df = master_distance_subtracks[1:length(seq(from=0, to=5.6, by=0.2)),] ### create a dataframe with the same size as the sequence of Time of interest:
  empty_df[!is.na(empty_df)] <- NA ## make it empyt
  empty_df$Time<-seq(from=0, to=5.6, by=0.2)
  empty_df$subsequence<-0.5
  
  
  ### now merge this fake dataframe with my dataframe of interest:
  master_distance<-rbind(master_distance_subtracks, empty_df)
  library(tidyr)
  ### complete all the datasets so that they have the same timepoints:
  master_distance<-master_distance%>%complete(subsequence,Time)
  ### remove the fake track that is not necessary anymore:
  
  master_distance<-subset(master_distance, !subsequence==0.5)
  
  master_distance$Time<-as.numeric(master_distance$Time)
  ##round the number, otherwise they are not equal
  master_distance$Time<-round(master_distance$Time,2)
  ##if there is any track2 and time that was duplicated remove the NA
  master_distance <- master_distance %>%
    group_by(Time, subsequence) %>%
    dplyr::slice(if (n() == 1) 1 else which(!is.na(dist_tumor)))
  
  ### create a for loop to refill all the NA with interpolated values:
  ### create a list of the columns names that need to be refilled
  column_names<-names(master_distance)
  column_names<-subset(column_names,!column_names %in%c("TrackID","Track2", "subsequence","ranks","ID","ID2","Time","position","mouse","class","pos", "Tracked","speed"))
  
  ### This is to test only ---------------
  ## create a first dataset with refilled values for speed:
  time_series<-reshape2::acast(master_distance, Time ~ subsequence, value.var='speed',fun.aggregate = mean)
  ## rownames timepoints:
  row.names(time_series)<-unique(master_distance$Time)
  ## get rid of NA
  time_series_zoo<-zoo(time_series, row.names(time_series))
  time_series_zoo<-na.approx(time_series_zoo) ## replace by interpolated value
  time_series<-as.matrix(time_series_zoo)
  time_series2<-melt(time_series)
  data<-time_series2[complete.cases(time_series2), ] 
  colnames(data)<-c("Time", "subsequence", "speed")
  ### ----------
  for (i in column_names){
    time_series<-acast(master_distance, Time ~ subsequence, value.var=i,fun.aggregate = mean)
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
  p1 <-ggplot(data, aes(Time, disp2, group=subsequence, color = as.factor(subsequence))) + 
    geom_smooth(method = "loess",size = 0.5, se = F, alpha=0.3, span=1)+
    theme_bw() +
    theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=5), axis.title.y = element_text(size = 10), axis.title.x = element_text(size = 15), legend.text=element_text(size= 10))+
    ggtitle("disp2")+theme(aspect.ratio=1,legend.position = "none")
  
  p1
  ### --------
  
  
  ### add metadata info
  
  master_cor <- data
  master_position<- master_distance[c("subsequence" ,"position")]
  master_position<-na.omit(master_position)
  master_position<-master_position[!duplicated(master_position$subsequence),]
  master_mouse<- master_distance[c("subsequence" ,"mouse")]
  master_mouse<-na.omit(master_mouse)
  master_mouse<-master_mouse[!duplicated(master_mouse$subsequence),]
  master_pos<- master_distance[c("subsequence" ,"class")]
  master_pos<-na.omit(master_pos)
  master_pos<-master_pos[!duplicated(master_pos$subsequence),]
  library(dplyr)
  #detach(package:plyr, unload = TRUE)
  master_cor<- left_join(master_cor ,master_position)
  master_cor<- left_join(master_cor ,master_mouse)
  master_cor<- left_join(master_cor ,master_pos)
  ### Filter timepoints of interest with same time interval:
  time_of_interest<-as.character(seq(from=0, to=5.6, by=0.2)) ### I don't know why but subsetting from a numeric vector doesn't work!!! So I convert to character
  master_cor$Time<-as.character(master_cor$Time)
  master_cor<-subset(master_cor, Time%in%time_of_interest)
  master_cor$Time<-as.numeric(master_cor$Time) ### now back to numeric
  setwd("D:/R/Projects/IVM")
  
  saveRDS(master_cor,"master_cor.rds")
  saveRDS(master_distance_MG,"master_distance_MG.rds")
  saveRDS(master_distance_SR101,"master_distance_SR101.rds")
  saveRDS(BV_df_sum,"BV_df_sum.rds")
}else{
  master_cor<-readRDS("master_cor.rds")
  master_distance_MG<-readRDS("master_distance_MG.rds")
  master_distance_SR101<-readRDS("master_distance_SR101.rds")
  BV_df_sum<-readRDS("BV_df_sum.rds")
  print("data already imported")
}



### sumarize the length of the tracks:
## create relative time 
master_cor2<-master_cor %>% 
  group_by(subsequence) %>%arrange(Time)%>%mutate(Time2 = Time - dplyr::first(Time))
ggplot(master_cor2,aes(Time2))+geom_histogram()+facet_wrap(.~ position)

max_time<-master_cor2%>%group_by(position)%>%summarise(max_Time=max(Time2))
print(min(max_time$max_Time)) #this is the minimal time that any position has
hist(max_time$max_Time)
master_cor2<-master_cor2 %>% 
  group_by(subsequence)%>%filter(Time2<2.8)  ##Minimal time for a position is 
library(scales)
### create a parameterd for distnace to tumor core
master_cor2$dist_tumor_1<-master_cor2$dist_tumor

master_cor2<-master_cor2 %>% 
  group_by(position)%>%mutate(dist_tumor_1 = scales::rescale(dist_tumor_1, to=c(0,100)))
### normalize the distance to tumor factor per position position
master_cor2<-master_cor2 %>% 
  group_by(subsequence) %>%arrange(Time2)%>%mutate(dist_tumor = dist_tumor - dplyr::first(dist_tumor))%>%ungroup()

#### check tracks
master_M4_test<-master_cor2%>%group_by(position, mouse)%>%mutate(x=x-min(x), y=y-min(y))
g<-ggplot(master_M4_test)+geom_path(aes(x, y, group=subsequence, colour=dist_tumor), arrow = arrow(ends = "last", type = "closed",length = unit(0.005, "inches")))+theme_bw()
g<-g + theme(legend.position = "none")+scale_colour_gradient2(low = "blue", mid = "grey" , high = "red") +facet_wrap(.~position)
g


### plot per position types the values of 
p <- ggplot(master_cor2)+geom_jitter(aes(x=class, y=disp2))+
  geom_bar(aes(x=class , y=disp_d),alpha=0.5,stat = "summary") +ggtitle("disp2 per position type")+theme_classic()+
  theme(aspect.ratio=0.7)

p




### plot per position types the values of 
p <- ggplot(master_cor2)+
  geom_bar(aes(x=class , y=disp_l),alpha=0.5,stat = "summary") +ggtitle("disp2 per position type")+theme_classic()+
  theme(aspect.ratio=0.7)

p



## PCA and scaling
library(parallel)
library(dplyr)
library(dtwclust)
library(stats)
library(scales)
#detach(package:spatstat, unload = TRUE)
###normalize the data:
master_norm<-master_cor2%>%ungroup()
master_norm$raw_dist_tumor<-master_norm$dist_tumor
column_names2<-names(master_cor2)
## select what data to normalize
column_names2<-subset(column_names2,!column_names2 %in%c("Time","Track2","subsequence","x","y","z","position","mouse","class","Time2", "dist_tumor_1", "dist_3_neigh", "dist_10_neigh", "iteration"))
master_norm<-as.data.frame(master_norm)


### perfrom PCA on the factors I want to use. Change this. I do not want to normalize per mouse.
master_scaled<-master_norm[c(column_names2, "mouse")]%>%group_by(mouse)%>%mutate(across(everything(), scale))%>%ungroup() ##first scale
##now to give more weight to the varibale of dist_tumor so the ability of cells to leave
#master_scaled[,c("dist_tumor")]<-master_scaled[,c("dist_tumor")]*2

master_pca<-prcomp(master_scaled[,-c(6)], scale = F)

#master_pca<-prcomp(master_norm[column_names2], scale = T)
summary(master_pca)

pcaCharts <- function(x) {
  x.var <- x$sdev ^ 2
  x.pvar <- x.var/sum(x.var)
  print("proportions of variance:")
  print(x.pvar)
  
  par(mfrow=c(2,2))
  plot(x.pvar,xlab="Principal component", ylab="Proportion of variance explained", ylim=c(0,1), type='b')
  plot(cumsum(x.pvar),xlab="Principal component", ylab="Cumulative Proportion of variance explained", ylim=c(0,1), type='b')
  screeplot(x)
  screeplot(x,type="l")
  par(mfrow=c(1,1))
}

pcaCharts(master_pca)
master_pca2<-master_pca[["x"]] ##select only PC
master_pca3<-data.frame(master_norm[,c(2)],master_pca2[,c(1:3)]) ### take only first 4 components since they explain most variation
colnames(master_pca3)[1]<- "subsequence"
#rescale_2 <- function(x) (x-min(x))/(max(x) - min(x)) * 100


###MULTIVARIATE
list_master_pca <- split(master_pca3[,-c(1)],master_pca3$subsequence) ## split into list
setwd("D:/R/Projects/IVM")

if ( ! file.exists("matrix_distmat_20240125_n_6_IVM_subsequence_pca3.rds")){
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
  distmat <- proxy::dist(list_master_pca, method = "dtw")
  matrix_distmat<-as.matrix(distmat)
  saveRDS(matrix_distmat, file = "matrix_distmat_20240125_n_6_IVM_subsequence__no_pca.rds")
  
  
}else{
  matrix_distmat <- readRDS("matrix_distmat_20240125_n_6_IVM_subsequence_pca3.rds")
}

#matrix_distmat <- readRDS("matrix_distmat_20240125_n_6_IVM_subsequence_pca3")
subsequence<-as.numeric(names(list_master_pca))

#Try WGCNA module detection
library(WGCNA)
# Convert similarity_matrix to a data frame if it's not already
similarity_matrix <- as.data.frame(matrix_distmat)

# Perform network construction and module detection
net <- blockwiseModules(similarity_matrix, power = 1, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, verbose = 3)

plotDendroAndColors(net$dendrograms[[1]], colors = labels2colors(net$colors))





# Try Leiden clustering
install.packages(c("igraph", "leiden"))
library(igraph)
library(leiden)

# Assuming 'similarity_matrix' is your similarity matrix
# Convert the similarity matrix to a graph object
graph <- graph_from_adjacency_matrix(matrix_distmat, mode = "undirected")

# Perform Leiden clustering
leiden_clusters <- leiden(graph, weights = E(graph)$weight)

# Access cluster assignments
cluster_assignments <- leiden_clusters$membership


library(umap)
set.seed(123)
umap_dist<- umap(matrix_distmat,n_components=2,input="dist",init = "random", 
                 n_neighbors=20,min_dist=0.09, spread=1.5)  
#3D plot if necessary:
#library(plot3Drgl)
umap_1 <- as.data.frame(umap_dist$`layout`) 
#scatter3Drgl(umap_1$V1, umap_1$V2, umap_1$V3)

plot(umap_dist$`layout`, # Plot the result
     col=rgb(0,0,0,alpha=0.1), 
     pch=19,
     asp=0.4)
umap_1 <- as.data.frame(umap_dist$`layout`) 
Track2_umap<-cbind(subsequence,umap_1)
positiontype<-master_norm[,c("subsequence", "position")]
positiontype<- positiontype[!duplicated(positiontype$subsequence),]
postype<-master_norm[,c("subsequence", "class")]
postype<- postype[!duplicated(postype$subsequence),]
mousetype<-master_norm[,c("subsequence","mouse")]
mousetype<- mousetype[!duplicated(mousetype$subsequence),]
umap_2 <- left_join(Track2_umap ,positiontype)
umap_2 <- left_join(umap_2 ,mousetype)
umap_2 <- left_join(umap_2 ,postype)

library(mclust)
set.seed(123)
mc.norm <- Mclust(umap_dist$`layout`,15, start =123) ## slower but better clustering
#km.norm <- kmeans(umap_dist$`layout`,8, nstart = 100)
#library(dbscan)
#cl <- dbscan(umap_dist$`layout`, eps= 0.35, minPts =15) 
## minPts is the min amount of values per cluster
#umap_3 <- cbind(cl$cluster, umap_2)
umap_3 <- cbind(mc.norm$classification, umap_2)
#umap_3 <- cbind(km.norm$cluster, umap_2)

#hc <- hclust(dist(umap_dist$layout))
#hierarchical_clusters <- cutree(hc, k = 12)  
#umap_3 <- cbind(hierarchical_clusters, umap_2)


colnames(umap_3)[1]<- "cluster2"


umap_3$cluster2<-net$colors
##plot
#umap_3<-subset(umap_3, cluster2!=0)
ggplot(umap_3, aes(x=V1, y=V2, color=as.factor(cluster2))) +  
  geom_point(size=2, alpha=0.6) + labs(color="cluster")+
  xlab("") + ylab("")+ 
  ggtitle("umap Cluster") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)

#scatter3Drgl(umap_3$V1, umap_3$V2, umap_3$V3, colvar = umap_3$cluster2)



ggplot(umap_3, aes(x=V1, y=V2, color=as.factor(class))) +  
  geom_point(size=2, alpha=0.6) + labs(color="Environment_cluster")+
  xlab("") + ylab("")+ 
  ggtitle("umap Cluster") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)+facet_wrap(mouse~position)

ggplot(umap_3, aes(x=V1, y=V2, color=as.factor(class))) +  
  geom_point(size=2, alpha=0.6) + labs(color="Environment_cluster")+
  xlab("") + ylab("")+ 
  ggtitle("umap Cluster") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)+facet_wrap(.~mouse)


ggplot(umap_3, aes(x=V1, y=V2, color=as.factor(class))) +  
  geom_point(size=2, alpha=0.6) + labs(color="Tumor region phenotype")+
  xlab("") + ylab("")+ 
  ggtitle("umap Cluster") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)+scale_color_manual(values=c("cyan","yellow","red"))




### calculate average stats per track
master_class <- left_join(master_norm ,umap_3[c(1:4)], by.x = c("subsequence"), by.y = c("subsequence"))

master_class$cluster2<-as.numeric(master_class$cluster2)
master_class_sum<-master_class%>%group_by(subsequence,position,class,mouse)%>%arrange(Time)%>%summarize(cluster2=mean(cluster2), V1=mean(V1), V2=mean(V2),dist_tumor=mean(dist_tumor), 
                                                                                                   dist_3_neigh=mean(dist_3_neigh),dist_10_neigh=mean(dist_10_neigh), speed=mean(speed), disp2=mean(disp2), 
                                                                                                   disp_d=mean(disp_d),disp_l=mean(disp_l), movement=last(dist_tumor), distance_to_tumor=last(dist_tumor_1),raw_dist_tumor=last(raw_dist_tumor))%>%ungroup()


## Join the information on the SR101 and MG
master_class_sum <- left_join(master_class_sum,master_distance_MG[c("subsequence","n_MG","min_MG")] ,by.x = c("subsequence"), by.y = c("subsequence"))
master_class_sum <- left_join(master_class_sum,master_distance_SR101[c("subsequence","n_SR101","min_SR101")] ,by.x = c("subsequence"), by.y = c("subsequence"))
master_class_sum <- left_join(master_class_sum,BV_df_sum ,by.x = c("subsequence"), by.y = c("subsequence"))

mid <- 0
master_class_sum$movement2<-ifelse(master_class_sum$movement>0,"away from tumor", ifelse(master_class_sum$movement==0, "no movement", "towards the tumor"))
p1<-ggplot(master_class_sum, mapping=aes(x=V1, y=V2, color=movement)) +  
  geom_point(size=2, alpha=0.8) +labs(color="Direction")+
  xlab("") + ylab("") +
  ggtitle("direction of movement relative to tumor core") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)+theme(aspect.ratio=1)+
  scale_color_gradient2(midpoint = mid, low = "blue",mid="grey90",
                        high = "red3")
p1

ggplot(umap_3, aes(x=V1, y=V2, color=as.factor(cluster2))) +  
  geom_point(size=2, alpha=0.6) + labs(color="cluster")+
  xlab("") + ylab("")+ 
  ggtitle("umap Cluster") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)



#scatter3Drgl(master_class_sum$V1, master_class_sum$V2, master_class_sum$V3, colvar = master_class_sum$movement)
master_class_sum<-master_class_sum[!is.na(master_class_sum$cluster2),]



library(dplyr)
library(scales)

###Full version
sum_all<-master_class_sum%>%group_by(cluster2)%>%
  summarise(nearest_3_cell= median(dist_3_neigh),nearest_10_cell= median(dist_10_neigh),displacement2 = mean(disp2), 
            speed = mean(speed), dist_tumor_core=median(distance_to_tumor),disp_d=mean(disp_d),
            disp_l=mean(disp_l), movement=median(movement), n_SR101=mean(n_SR101, na.rm = T),
            min_SR101=mean(min_SR101, na.rm = T),n_MG=mean(n_MG, na.rm = T),min_MG=mean(min_MG, na.rm = T),min_BV_dist=mean(BV_min, na.rm = T),sd_BV_dist=mean(BV_sd, na.rm = T),mean_BV_dist=mean(BV_mean, na.rm = T),max_BV_dist=mean(BV_max, na.rm = T),BV_contact=mean(BV_contact, na.rm = T))



### Version for paper
sum_all<-master_class_sum%>%group_by(cluster2)%>%
  summarise(displacement2 = mean(disp2), 
            speed = mean(speed),displacement_delta=mean(disp_d),displacement_length=mean(disp_l), nearest_10_cell= median(dist_10_neigh),dist_tumor_core=median(distance_to_tumor),
            n_SR101=mean(n_SR101, na.rm = T),
            min_SR101=mean(min_SR101, na.rm = T),n_MG=mean(n_MG, na.rm = T),min_MG=mean(min_MG, na.rm = T),min_BV_dist=mean(BV_min, na.rm = T),sd_BV_dist=mean(BV_sd, na.rm = T),mean_BV_dist=mean(BV_mean, na.rm = T),
            movement=median(movement))



#sum_all[is.na(sum_all)]<-0
Cluster_movement<-sum_all[,c(1,15)]
Cluster_movement <-Cluster_movement[order(Cluster_movement$movement),]

stdize = function(x, ...) {(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

sum_all[-c(1)] <- lapply(sum_all[-c(1)], stdize, na.rm = T)

library(viridis)
library(plotly)

df<-sum_all[,-c(1)]

df<-df[order(match(rownames(df), Cluster_movement$cluster2)), , drop = FALSE]%>%ungroup()

rownames(df)<-Cluster_movement$cluster2
library(pheatmap)
heat_m<-pheatmap(df, clustering_method = "complete", cluster_cols=F,cluster_rows=F,cellwidth=20, cellheight=20,
                 treeheight_row = 0, treeheight_col = 0,fontsize = 8,na_col = "grey50",
                 angle_col = 90,color = colorRampPalette(c("green2", "grey90", "pink"))(50))

Cluster_movement





## plot clusters giving them a color according to direction
library(scales)
df1<-as.data.frame(df)
mypal <- col_numeric(palette = c("blue", "green2", "red"), domain = 0:100)
df1$cluster2<-row.names(df)
df1 <-df1[order(df1$movement),] ##order values from low to gih
my_palette <- mypal(unique(df1$movement)) ## create palette based on colors


library(scales)
pal<-rescale(Cluster_movement$movement, to=c(1,100))

my_palette2 <- mypal(pal) ## create palette based on colors

color_names <- c("blue", "cadetblue", "darkolivegreen", "mediumseagreen", "darkolivegreen",
                       "olivedrab", "limegreen", "green3", "greenyellow", "chartreuse",
                       "green", "red")

##order clusters to correspond to the palette
umap_3$cluster2 = factor(umap_3$cluster2, levels=df1$cluster2)
p2<-ggplot(umap_3, aes(x=V1, y=V2, color=as.factor(cluster2))) +  
  geom_point(size=2, alpha=1) + labs(color="cluster")+
  xlab("") + ylab("")+ scale_color_manual(values=my_palette2)+
  ggtitle("umap Cluster ") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)


p2
ggplot(master_class_sum, mapping=aes(x=V1, y=V2, color=distance_to_tumor)) +  
  geom_point(size=2, alpha=0.6) +
  xlab("") + ylab("") +
  ggtitle("UMAP localization distribution") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)+
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(4, "Spectral")))


ggplot(master_class_sum, mapping=aes(x=V1, y=V2, color=raw_dist_tumor)) +  
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



ggplot(master_class_sum, mapping=aes(x=V1, y=V2, color=BV_contact)) +  
  geom_point(size=2, alpha=0.6) +
  xlab("") + ylab("") +
  ggtitle("UMAP mean contact with Blood vessels") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)+
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(4, "Spectral")),na.value="NA")




### convert the highest nearest enighdor for coloring purposes
#master_class_sum$dist_10_neigh<-ifelse(master_class_sum$dist_10_neigh>3,3,master_class_sum$dist_10_neigh)


ggplot(master_class_sum, mapping=aes(x=V1, y=V2, color=dist_10_neigh)) +  
  geom_point(size=2, alpha=0.6) +
  xlab("") + ylab("") +
  ggtitle("UMAP close neighbours distribution") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())+theme(aspect.ratio=1)+
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(4, "Spectral")))




### plot enrichment on edge and core:
Number_cell_exp<-umap_3%>%group_by(position,class,mouse)%>%
  summarise(total_cell = n())
Percentage_clus<-left_join(Number_cell_exp,umap_3)
Percentage_clus <- Percentage_clus%>%group_by(cluster2,position,class,mouse)%>%
  summarise(total_cell = mean(total_cell), num_cluster=n())%>%mutate(percentage=num_cluster*100/total_cell)%>%ungroup()
Percentage_clus_2 <- Percentage_clus%>%group_by(cluster2, position,class,mouse)%>%summarise(se_percentage=sd(percentage)/sqrt(length(percentage)),percentage = mean(percentage))
### Plot circular map

Per1<-ggplot(Percentage_clus_2, aes(fill=as.factor(cluster2), y=percentage, x="")) + 
  geom_bar( stat="identity", position="fill", width = 0.5)+ coord_flip()+ scale_y_reverse()
Per1 <- Per1 + facet_grid(class~mouse)
Per1<-Per1+theme_void()+ scale_fill_manual(values=my_palette2)+theme(strip.text.x = element_text(size = 8, angle = 90))
Per1



Per2<-ggplot(Percentage_clus_2, aes(fill=as.factor(cluster2), y=percentage, x="")) + 
  geom_bar( stat="identity", position="fill", width = 0.5)+ coord_flip()+ scale_y_reverse()
Per2 <- Per2 + facet_grid(.~mouse)
Per2<-Per2+theme_void()+ scale_fill_manual(values=my_palette2)+theme(strip.text.x = element_text(size = 8, angle = 90))
Per2





Per3<-ggplot(Percentage_clus_2, aes(fill=as.factor(cluster2), y=percentage, x="")) + 
  geom_bar( stat="identity", position="fill", width = 0.5)+ coord_flip()+ scale_y_reverse()
Per3 <- Per3 + facet_grid(class~.)
Per3<-Per3+theme_void()+theme(strip.text.x = element_text(size = 8, angle = 90))+ scale_fill_manual(values=my_palette2)
Per3

## Difference between clusters

library(lme4)
library(emmeans)

# Assuming 'cluster2' and 'class' are factor variables in Percentage_clus_2
Percentage_clus_2$cluster2 <- as.factor(Percentage_clus_2$cluster2)
Percentage_clus_2$class <- as.factor(Percentage_clus_2$class)

# Create an empty list to store pairwise comparison results
pairwise_list <- list()

# Get unique levels of cluster2
cluster2_levels <- levels(Percentage_clus_2$cluster2)

# Loop through each level of cluster2
for (cluster_level in cluster2_levels) {
  # Subset the data for the current cluster level
  subset_data <- subset(Percentage_clus_2, cluster2 == cluster_level)
  
  # Fit the mixed-effects model for the current cluster level
  model <- lmer(percentage ~ class + (1|mouse), data = subset_data)
  
  # Obtain estimated marginal means (EMMs) for pairwise comparisons
  emm_model <- emmeans(model, "class")
  
  # Perform pairwise comparisons
  pairwise_comparisons <- pairs(emm_model, text.panel=NULL,upper.panel = NULL)
  
  # Store the results in the list
  pairwise_list[[as.character(cluster_level)]] <- summary(pairwise_comparisons)
}



#### cFisher test to see if all clusters are equally distributed in Class:
Percentage_clus_3 <- Percentage_clus%>%group_by(cluster2,class)%>%summarise(se_percentage=sd(percentage)/sqrt(length(percentage)),percentage = mean(percentage))


library(reshape2)
table1<-acast(Percentage_clus_3, cluster2~class, value.var = "percentage")
result_chisq<-chisq.test(table1)  ## Chisq test, that indicates if there are differences between categories
result_chisq
chisq.test(table1[,c(1,3)])
## Choran=Armitage posthoc test:
library(rcompanion)

table2<-umap_3[c(1,7)]
table2$cluster2<-as.factor(table2$cluster2)
table2<-table(table2)
PT = pairwiseNominalIndependence(table2, fisher=F, gtest=F,
                                 compare="column", method = "bonferroni")  ### pairwise comparison

PT

capture.output(PT, file="E:/SURF_2/Intravital imaging/DATA_IVM_OBS_paper/Results_MA18_MA21_MA22_IVM/percent_chisq_independence.txt")

### save as txt the output of this comparison. Change the path for your own data
capture.output(result_chisq, file="E:/SURF_2/Intravital imaging/DATA_IVM_OBS_paper/Results_MA18_MA21_MA22_IVM/percent_chisq.txt")

### plot position relative to tumor per cluster
master_class_sum<-as.data.frame(master_class_sum)
master_class_sum$cluster2<-as.character(master_class_sum$cluster2)
df1$mean_movement<-df1$movement
master_class_sum2<-left_join(master_class_sum, df1[,c("cluster2","mean_movement")])
master_class_sum2$cluster2 = factor(master_class_sum2$cluster2, levels=df1$cluster2)
#### try an anova between clusters based speed 
p3 <- ggplot(master_class_sum2,aes(x=as.factor(cluster2) , y=speed, group=cluster2,fill=cluster2)) +
  geom_boxplot(outlier.colour = NA) +ggtitle("speed per cluster")+theme_classic()+
  theme(aspect.ratio=0.7)+scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(4, "Spectral")))+coord_cartesian(ylim = c(0,20))

p3

a3 <- aov(speed ~ cluster2, data=master_class_sum2) 
summary(a3)
TukeyHSD(a3)
### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a3), file="E:/SURF_2/Intravital imaging/DATA_IVM_OBS_paper/Results_MA18_MA21_MA22_IVM/Tukey_speed.txt")


#### try an anova between clusters based on distance to dist_3_neigh
p <- ggplot(master_class_sum2,aes(x=as.factor(cluster2) , y=dist_3_neigh, group=cluster2)) +
  geom_boxplot(aes(colour=mean_movement),alpha=0.5, outlier.colour = NA) +ggtitle("dist 3 neigh per cluster")+theme_classic()+
  theme(aspect.ratio=0.7)+scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(4, "Spectral")))+coord_cartesian(ylim = c(10,60))

p

a1 <- aov(dist_3_neigh ~ cluster2, data=master_class_sum2) 
summary(a1)
TukeyHSD(a1)


#### try an anova between clusters based on distance to dist_10_neigh
p4 <- ggplot(master_class_sum2,aes(x=as.factor(cluster2) , y=dist_10_neigh, group=cluster2,fill=cluster2)) +
  geom_boxplot( outlier.colour = NA) +ggtitle("dist 10 neigh per cluster")+theme_classic()+
  theme(aspect.ratio=0.7)+scale_fill_manual(values = color_names)+coord_cartesian(ylim = c(20,70))

p4

a4 <- aov(dist_10_neigh ~ cluster2, data=master_class_sum2) 
summary(a4)
TukeyHSD(a4)

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a4), file="E:/SURF_2/Intravital imaging/DATA_IVM_OBS_paper/Results_MA18_MA21_MA22_IVM/Tukey_dist_10_nei.txt")


#### try an anova between clusters based on dist to MG
p5 <- ggplot(subset(master_class_sum2, !is.na(min_MG)) ,aes(x=as.factor(cluster2) , y=min_MG, group=cluster2, fill=cluster2))  +
  geom_boxplot( outlier.colour = NA) +ggtitle("distance to closest MG")+theme_classic()+
  theme(aspect.ratio=0.7)+scale_fill_manual(values = color_names)+coord_cartesian(ylim = c(2,45))

p5


a5 <- aov(min_MG ~ cluster2, data= subset(master_class_sum2, !is.na(min_MG))) 
summary(a5)
TukeyHSD(a5)


### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a5), file="E:/SURF_2/Intravital imaging/DATA_IVM_OBS_paper/Results_MA18_MA21_MA22_IVM/Tukey_min_MG_dist.txt")



#### try an anova between clusters based on dist to SR101
p6 <- ggplot(subset(master_class_sum2, !is.na(min_SR101)) ,aes(x=as.factor(cluster2) , y=min_SR101, group=cluster2, fill=cluster2))  +
  geom_boxplot(outlier.colour = NA) +ggtitle("distance to closest SR101")+theme_classic()+
  theme(aspect.ratio=0.7)+scale_fill_manual(values = color_names)+coord_cartesian(ylim = c(0,50))

p6


a6 <- aov(min_SR101 ~ cluster2, data= subset(master_class_sum2, !is.na(min_SR101))) 
summary(a6)
TukeyHSD(a6)


### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a6), file="E:/SURF_2/Intravital imaging/DATA_IVM_OBS_paper/Results_MA18_MA21_MA22_IVM/Tukey_min_SR101_dist.txt")



#### try an anova between clusters based on n of  MG
p <- ggplot(subset(master_class_sum2, !is.na(min_MG)),aes(x=as.factor(cluster2) , y=n_MG, group=cluster2, fill=cluster2)) +
  geom_bar(alpha=0.5,stat = "summary") +ggtitle("n of MG per cluster")+theme_classic()+
  theme(aspect.ratio=0.7)+scale_fill_manual(values = color_names)

p


p7 <- ggplot(subset(master_class_sum2, !is.na(min_MG)) ,aes(x=as.factor(cluster2) , y=n_MG, group=cluster2, fill=cluster2))  +
  geom_boxplot( outlier.colour = NA) +ggtitle("n MG")+theme_classic()+
  theme(aspect.ratio=0.7)+scale_fill_manual(values = color_names)+coord_cartesian(ylim = c(0,30))

p7


a7 <- aov(n_MG ~ cluster2, data= subset(master_class_sum2, !is.na(min_MG))) 
summary(a7)
TukeyHSD(a7)
### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a7), file="E:/SURF_2/Intravital imaging/DATA_IVM_OBS_paper/Results_MA18_MA21_MA22_IVM/Tukey_n_MG.txt")




#### try an anova between clusters based on n of  SR101

p <- ggplot(subset(master_class_sum2, !is.na(min_SR101)),aes(x=as.factor(cluster2) , y=n_SR101, group=cluster2)) +
  geom_bar(alpha=0.5,stat = "summary") +ggtitle("n SR101 per cluster")+theme_classic()+
  theme(aspect.ratio=0.7)+scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(4, "Spectral")))

p

p8 <- ggplot(subset(master_class_sum2, !is.na(min_SR101)) ,aes(x=as.factor(cluster2) , y=n_SR101, group=cluster2, fill=cluster2))  +
  geom_boxplot(outlier.colour = NA) +ggtitle("n SR101")+theme_classic()+
  theme(aspect.ratio=0.7)+scale_fill_manual(values = color_names)+coord_cartesian(ylim = c(0,20))

p8



a8 <- aov(n_SR101 ~ cluster2, data= subset(master_class_sum2, !is.na(min_SR101))) 
summary(a8)
TukeyHSD(a8)
### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a8), file="E:/SURF_2/Intravital imaging/DATA_IVM_OBS_paper/Results_MA18_MA21_MA22_IVM/Tukey_n_SR101.txt")



#### try an anova between clusters based on mean distance to BV

p9 <- ggplot(subset(master_class_sum2, !is.na(BV_mean)) ,aes(x=as.factor(cluster2) , y=BV_mean, group=cluster2, fill=cluster2))  +
  geom_boxplot(outlier.colour = NA) +ggtitle("BV_mean")+theme_classic()+
  theme(aspect.ratio=0.7)+scale_fill_manual(values = color_names)+coord_cartesian(ylim = c(0,20))

p9



a9 <- aov(BV_mean ~ cluster2, data= subset(master_class_sum2, !is.na(BV_mean))) 
summary(a9)
TukeyHSD(a9)

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a9), file="E:/SURF_2/Intravital imaging/DATA_IVM_OBS_paper/Results_MA18_MA21_MA22_IVM/Tukey_mean_dist_BV.txt")


#### try an anova between clusters based on min distance to BV

p10 <- ggplot(subset(master_class_sum2, !is.na(BV_min)) ,aes(x=as.factor(cluster2) , y=BV_min, group=cluster2, fill=cluster2))  +
  geom_boxplot( outlier.colour = NA) +ggtitle("BV_min")+theme_classic()+
  theme(aspect.ratio=0.7)+scale_fill_manual(values = color_names)+coord_cartesian(ylim = c(0,20))

p10



a10 <- aov(BV_min ~ cluster2, data= subset(master_class_sum2, !is.na(BV_min))) 
summary(a10)
TukeyHSD(a10)


### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a10), file="E:/SURF_2/Intravital imaging/DATA_IVM_OBS_paper/Results_MA18_MA21_MA22_IVM/Tukey_min_dist_BV.txt")



#### try an anova between clusters based on contact to BV

p <- ggplot(subset(master_class_sum2, !is.na(BV_min)) ,aes(x=as.factor(cluster2) , y=BV_contact, group=cluster2, fill=cluster2))  +
  geom_boxplot(outlier.colour = NA) +ggtitle("BV_contact")+theme_classic()+
  theme(aspect.ratio=0.7)+scale_fill_manual(values = color_names)

p



a1 <- aov(BV_contact ~ cluster2, data= subset(master_class_sum2, !is.na(BV_contact))) 
summary(a1)
TukeyHSD(a1)




#### try an anova between clusters based on sd distance to BV

p10.2 <- ggplot(subset(master_class_sum2, !is.na(BV_min)) ,aes(x=as.factor(cluster2) , y=BV_sd, group=cluster2, fill=cluster2))  +
  geom_boxplot( outlier.colour = NA) +ggtitle("Standart deviation distance")+theme_classic()+
  theme(aspect.ratio=0.7)+scale_fill_manual(values = color_names)+coord_cartesian(ylim = c(0,1.7))

p10.2



a10.2 <- aov(BV_sd ~ cluster2, data= subset(master_class_sum2, !is.na(BV_min))) 
TukeyHSD(a10.2)


### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a10.2), file="E:/SURF_2/Intravital imaging/DATA_IVM_OBS_paper/Results_MA18_MA21_MA22_IVM/Tukey_sd_dist_BV.txt")




####################################################################################################
####Differences based on environmental clusters

#### try an anova between clusters based on distance to speed
p11 <- ggplot(master_class_sum2,aes(x=as.factor(class) , y=speed, group=class,fill=class)) +
  geom_boxplot(alpha=0.5, outlier.colour = NA) +scale_fill_manual(values=c("cyan","gold","red"))+ggtitle("speed per env cluster")+theme_classic()+
  theme(aspect.ratio=0.7)+coord_cartesian(ylim = c(0,17))

p11


a11 <- aov(speed ~ class, data=master_class_sum2) 
summary(a11)
TukeyHSD(a11)

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a11), file="E:/SURF_2/Intravital imaging/DATA_IVM_OBS_paper/Results_MA18_MA21_MA22_IVM/Tukey_speed_cytomap_cl.txt")



#### try an anova between clusters based on raw tumor movement
p12 <- ggplot(master_class_sum2,aes(x=as.factor(class) , y=raw_dist_tumor, group=class,fill=class)) +
  geom_boxplot(alpha=0.5, outlier.colour = NA)+scale_fill_manual(values=c("cyan","gold","red")) +ggtitle("movement direction per env cluster")+theme_classic()+
  theme(aspect.ratio=0.7)+scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(4, "Spectral")))+coord_cartesian(ylim = c(-10,17))

p12


a12 <- aov(raw_dist_tumor ~ class, data=master_class_sum2) 
summary(a12)
TukeyHSD(a12)

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a12), file="E:/SURF_2/Intravital imaging/DATA_IVM_OBS_paper/Results_MA18_MA21_MA22_IVM/Tukey_move_dir_cytomap_cl.txt")



#### try an anova between clusters based on distance to dist_10_neigh
p <- ggplot(master_class_sum2,aes(x=as.factor(class) , y=dist_10_neigh, fill=class)) +
  geom_boxplot(alpha=0.5, outlier.colour = NA)+scale_fill_manual(values=c("cyan","gold","red"))  +ggtitle("dist 10 neigh per cluster")+theme_classic()+
  theme(aspect.ratio=0.7)+coord_cartesian(ylim = c(20,70))

p
a1 <- aov(dist_10_neigh ~ class, data=master_class_sum2) 
summary(a1)
TukeyHSD(a1)



#### try an anova between clusters based on dist to MG
p13 <- ggplot(subset(master_class_sum2, !is.na(min_MG)) ,aes(x=as.factor(class) , y=min_MG, fill=class))  +
  geom_boxplot(alpha=0.5, outlier.colour = NA)+scale_fill_manual(values=c("cyan","gold","red"))  +ggtitle("distance to closest MG")+theme_classic()+
  theme(aspect.ratio=0.7)+coord_cartesian(ylim = c(0,50))

p13


a13 <- aov(min_MG ~ class, data= subset(master_class_sum2, !is.na(min_MG))) 
summary(a13)
TukeyHSD(a13)

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a13), file="E:/SURF_2/Intravital imaging/DATA_IVM_OBS_paper/Results_MA18_MA21_MA22_IVM/Tukey_distMG_cytomap_cl.txt")



#### try an anova between clusters based on dist to SR101
p14 <- ggplot(subset(master_class_sum2, !is.na(min_SR101)) ,aes(x=as.factor(class) , y=min_SR101, fill=class))  +
  geom_boxplot(alpha=0.5, outlier.colour = NA) +scale_fill_manual(values=c("cyan","gold","red"))+ggtitle("distance to closest SR101")+theme_classic()+
  theme(aspect.ratio=0.7)+coord_cartesian(ylim = c(0,50))

p14


a14 <- aov(min_SR101 ~ class, data= subset(master_class_sum2, !is.na(min_SR101))) 
summary(a14)
TukeyHSD(a14)

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a14), file="E:/SURF_2/Intravital imaging/DATA_IVM_OBS_paper/Results_MA18_MA21_MA22_IVM/Tukey_distSR101_cytomap_cl.txt")



#### try an anova between clusters based on n of  MG
p15 <- ggplot(subset(master_class_sum2, !is.na(min_MG)),aes(x=as.factor(class) , y=n_MG, fill=class)) +
  geom_boxplot(alpha=0.5, outlier.colour = NA)  +scale_fill_manual(values=c("cyan","gold","red"))+ggtitle("n of MG per cluster")+theme_classic()+
  theme(aspect.ratio=0.7)

p15

a15 <- aov(n_MG ~ class, data= subset(master_class_sum2, !is.na(min_MG))) 
summary(a15)
TukeyHSD(a15)

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a15), file="E:/SURF_2/Intravital imaging/DATA_IVM_OBS_paper/Results_MA18_MA21_MA22_IVM/Tukey_nMG_cytomap_cl.txt")



#### try an anova between clusters based on n of  SR101

p16 <- ggplot(subset(master_class_sum2, !is.na(min_SR101)),aes(x=as.factor(class) , y=n_SR101, fill=class)) +
  geom_boxplot(alpha=0.5,outlier.colour = NA) +scale_fill_manual(values=c("cyan","gold","red"))+ggtitle("n SR101 per cluster")+theme_classic()+
  theme(aspect.ratio=0.7)
p16

a16 <- aov(n_SR101 ~ class, data= subset(master_class_sum2, !is.na(min_SR101))) 
summary(a16)
TukeyHSD(a16)

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a16), file="E:/SURF_2/Intravital imaging/DATA_IVM_OBS_paper/Results_MA18_MA21_MA22_IVM/Tukey_nSR101_cytomap_cl.txt")


#### try an anova between clusters based on BV distance min
p17 <- ggplot(subset(master_class_sum2, !is.na(BV_min)) ,aes(x=as.factor(class) , y=BV_min, fill=class))  +
  geom_boxplot(alpha=0.5, outlier.colour = NA)+scale_fill_manual(values=c("cyan","gold","red"))  +ggtitle("min distance to BV")+theme_classic()+
  theme(aspect.ratio=0.7)+coord_cartesian(ylim = c(0,18))

p17


a17 <- aov(BV_min ~ class, data= subset(master_class_sum2, !is.na(BV_min))) 
summary(a17)
TukeyHSD(a17)

### save as txt the output of this comparison. Change the path for your own data
capture.output(TukeyHSD(a17), file="E:/SURF_2/Intravital imaging/DATA_IVM_OBS_paper/Results_MA18_MA21_MA22_IVM/Tukey_min_distBV_cytomap_cl.txt")




#### try an anova between clusters based on BV distance
p <- ggplot(subset(master_class_sum2, !is.na(BV_contact)) ,aes(x=as.factor(class) , y=BV_contact, fill=class))  +
  geom_boxplot(alpha=0.5, outlier.colour = NA)+scale_fill_manual(values=c("cyan","gold","red"))  +ggtitle("mean distance to BV")+theme_classic()+
  theme(aspect.ratio=0.7)+coord_cartesian(ylim = c(0,1.5))

p


a1 <- aov(BV_contact ~ class, data= subset(master_class_sum2, !is.na(BV_min))) 
summary(a1)
TukeyHSD(a1)


################################################################################################
##relation between position relative to the tumor and amount of MG:

p18 <- ggplot(subset(master_class_sum2, !is.na(min_MG)),aes(x=distance_to_tumor , y=n_MG)) +geom_jitter()+
  geom_point(alpha=0.5,stat = "summary") +geom_smooth(method = "lm")+ggtitle("scatter plot movement vs min_MG")+theme_classic()+
  theme(aspect.ratio=0.7)

p18

cor_18<-cor.test(master_class_sum2$n_MG,master_class_sum2$distance_to_tumor,  method = "pearson", use = "pairwise.complete.obs")
cor_18

### save as txt the output of this comparison. Change the path for your own data
capture.output(cor_18, file="E:/SURF_2/Intravital imaging/DATA_IVM_OBS_paper/Results_MA18_MA21_MA22_IVM/Pearson_cor_nMG_dist_tum.txt")





pdf("E:/SURF_2/Intravital imaging/DATA_IVM_OBS_paper/Results_MA18_MA21_MA22_IVM/results_MA_18_21_22_20210701.pdf" ) ##adjust path
p1
p2
p3
p4
p5
p6
p7
p8
p9
p10
p10.2
p11
p12
p13
p14
p15
p16
p17
p18
Per3

dev.off()

pdf("E:/SURF_2/Intravital imaging/DATA_IVM_OBS_paper/Results_MA18_MA21_MA22_IVM/MA_18_21_22_heamap.pdf" ) ##adjust path
heat_m

dev.off()



########################################################################################################
### export cell tracks per cluster

### extrack track numbers for backprojection
which( colnames(master_class)=="cluster2" ) ## column number
names_cell<-master_class[!duplicated(master_class$mouse),]
names_cell$position
Tracks_1<-master_class[which(master_class$position=="2430F11.CL3_1"),]
Tracks_1<-Tracks_1[!duplicated(Tracks_1$Track2),c(2,19)]

master_1<-master[which(master$position=="2430F11.CL3_1"),] ### based on the ranks number identify how to convert to original TrackID


Tracks_1$TrackID<-Tracks_1$Track2 - (min(master_1$ranks)*10000000000)  ##reconvert to TrackID
Track_1_list<-split(Tracks_1,Tracks_1$cluster2)

write(paste(as.character(Track_1_list), sep="' '", collapse=", "), "E:/SURF_2/Intravital imaging/DATA_IVM_OBS_paper/IVM cluster backprojection/2430F11_CL3_1/2430F11_CL3_1.txt")





###Plot cell tracjectories:
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

g<-ggplot(master_traject)+geom_path(aes(x_new, y_new, group=Track2, 
                                        color=as.factor(cluster2)), size=0.2,
                                    arrow = arrow(ends = "last",type = "open",length = unit(0.01, "inches")))+scale_color_manual(values=my_palette2)+
  theme_bw()+theme(aspect.ratio = 1)+facet_wrap(class~position)
g

g<-ggplot(master_traject)+geom_path(aes(x_new, y_new, group=Track2, color=as.factor(cluster2)), size=0.2,arrow = arrow(ends = "last",type = "open",length = unit(0.05, "inches")))+scale_color_manual(values=my_palette2)+
  theme_bw()+theme(aspect.ratio = 1)
g


### calculate mean values per cell:



### plot tracks based on speed only:
master_traject$cluster_sp<-ifelse(master_traject$cluster2%in%c(2,5), "fast",ifelse(master_traject$cluster2%in%c(7,6,3,9), "slow", "medium") )
g<-ggplot(master_traject)+geom_path(aes(x_new, y_new, group=Track2, 
                                        color=as.factor(cluster_sp)), size=0.5,
                                    arrow = arrow(ends = "last",type = "open",length = unit(0.01, "inches")))+
  theme_bw()+theme(aspect.ratio = 1)+
  coord_cartesian(xlim=c(0,1200), ylim=c(0,900))+
  facet_wrap(position~.)
g

### plot the whole tumor
g<-ggplot(master_traject)+geom_path(aes(x, y, group=Track2, 
                                        color=as.factor(cluster2)), size=0.5,
                                    arrow = arrow(ends = "last",type = "open",length = unit(0.01, "inches")))+
  theme_bw()+theme(aspect.ratio = 1)+scale_color_manual(values=my_palette)
g
## normalized to start

master_traject<-master_traject%>%group_by(Track2)%>%arrange(Time)%>%mutate(x_norm=x_new-first(x_new), y_norm=y_new-first(y_new))

g<-ggplot(master_traject)+geom_path(aes(x_norm, y_norm, group=Track2, color=as.factor(cluster2)), 
                                    arrow = arrow(ends = "last", type = "closed",length = unit(0.05, "inches")))+
  theme_bw()+scale_color_manual(values=my_palette)+
  facet_wrap(position~.)
g


#### Spatial heterogeneity of behavior:
### calculate for each cell to which cluster does it neighbor beloong to and then check if it is enrihed per cluster 
colnames(Cluster_movement)[2]<-"direction"
Cluster_movement$cluster2<-as.factor(Cluster_movement$cluster2)
master_traject<-left_join(master_traject,Cluster_movement)

master_traject_mean<-master_traject%>%group_by(Track2, position,mouse, pos)%>%
  summarise_all(mean)%>%
  ungroup()
column_names3<-names(master_traject_mean)
column_names3<-subset(column_names3,!column_names3 %in%c("Time","Track2","position","mouse","pos","Time2"))

master_traject_mean[column_names3] <- scale(master_traject_mean[column_names3],center = TRUE, scale = TRUE)

### 

library(sp)
coordinates(master_traject_mean) <- ~x + y  #effectively convert the data into a spatial data frame

### Calculate spatial autocorrelation
library(fields)
### create a coordinate frame for the x and w
master_traject_mean_c2<-master_traject_mean
coords_2 <- as.data.frame(cbind(master_traject_mean_c2$x, master_traject_mean_c2$y))
w2 <- fields:::rdist(coords_2)
library(ape) ## calculate the Moran spatial autocorrelation
moran<-Moran.I(x = master_traject_mean_c2$dist_tumor, w = w2)

moran
List2=list()
for (m in unique(master_traject_mean$position)){
  master_traject_mean_c<-master_traject_mean[which(master_traject_mean$position==m),]
  coords <- as.data.frame(cbind(master_traject_mean_c$x, master_traject_mean$y))
  w <- fields:::rdist(coords)
  moran<-as.data.frame(Moran.I(x = master_traject_mean_c$speed, w = w))
  List2[[length(List2)+1]] <-moran
}
moran_per_position<-do.call(rbind, List2)
moran_per_position$position<-unique(master_traject_mean$position)



List2=list()
for (m in unique(master_traject_mean$position)){
  master_traject_mean_c<-master_traject_mean[which(master_traject_mean$position==m),]
  coords <- as.data.frame(cbind(master_traject_mean_c$x_new, master_traject_mean_c$y_new))
  w <- fields:::rdist(coords)
  moran<-as.data.frame(Moran.I(x = master_traject_mean_c$cluster_sp, w = w))
  List2[[length(List2)+1]] <-moran
}
moran_per_position<-do.call(rbind, List2)
moran_per_position$position<-unique(master_traject_mean$position)


#### measure for the whole tumor (all position combine)
master_traject_mean2<-master_traject_mean
coords <- as.data.frame(cbind(master_traject_mean2$x, master_traject_mean2$y))
w <- fields:::rdist(coords)
library(ape) ## calculate the Moran spatial autocorrelation
moran_speed<-as.data.frame(Moran.I(x = master_traject_mean2$speed, w = w))
moran_dist_tumor<-as.data.frame(Moran.I(x = master_traject_mean2$dist_tumor, w = w))
moran_cluster<-as.data.frame(Moran.I(x = master_traject_mean2$direction, w = w))

moran_speed
#capture.output(moran_speed, file="E:/R/Intravital/Intravital/exports/20191025_Exp1F4_4pos/moran_speed.txt")
moran_dist_tumor
#capture.output(moran_dist_tumor, file="E:/R/Intravital/Intravital/exports/20191025_Exp1F4_4pos/moran_direction.txt")
moran_cluster
#capture.output(moran_cluster, file="E:/R/Intravital/Intravital/exports/20191025_Exp1F4_4pos/moran_cluster.txt")

library(sp)
bubble(master_traject_mean, "speed", maxsize = 1.5)
bubble(master_traject_mean, "dist_tumor", maxsize = 1.5)
bubble(master_traject_mean, "direction", maxsize = 1)
bubble(master_traject_mean, "dist_BV", maxsize = 1)
bubble(master_traject_mean, "dist_MG", maxsize = 1)
bubble(master_traject_mean, "nearest_cell", maxsize = 1)
bubble(master_traject_mean, "acceleration", maxsize = 1)
bubble(master_traject_mean, "dist_BV", maxsize = 1)
bubble(master_traject_mean, "dist_MG", maxsize = 1)


