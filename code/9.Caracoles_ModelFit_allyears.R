# This script will run the empirical model fits for each focal species and each
#       environmental covariate. A separate script will then make the empirical
#       figures for the manuscript
# https://mc-stan.org/docs/2_29/functions-reference/index.html#overview

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 1. SET UP: Import data, create df with abundances ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#rm(list = ls())
setwd("~/Eco_Bayesian/Complexity_caracoles")
#---- 1.1. Import packages ----
#installyears.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
#system.file("~/home/lbuche/R/x86_64-pc-linux-gnu-library/4.1/rstan",package='rstan')
options(mc.cores = parallyearsel::detectCores())

library(rstan)
library("codetools")
library(HDInterval)
library(matrixStats)
library(tidyverse)
library(stringr)
library(dplyr)
library(ggplot2) 

#---- 1.2. IMPORT DATA -----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#rm(list = ls()) # remove environment to not use the objects created before

home.dic <- "/Users/lisabuche/Documents/Projects/Complexity_caracoles/"
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 2. Run analysis ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

args <- commandArgs(trailingOnly = TRUE)

#Extract argument from slurms

run.prelimfit <- 1
run.finalfit <- 1
focal.levels <- c("HOMA","CHFU","LEMA","CETE")
year <-  c("2019",'2020','2021')
complexity.level.int <- c(1,2,3)
summary.interactions <- NULL

for( focal in focal.levels){ 
    for( complexity.level in complexity.level.int){
      
      summary.interactions.n <- data.frame()
      complexity.plant <-c("class","family","code.plant")[complexity.level]
      complexity.animal <- c("group","family","species")[complexity.level]
      FocalPrefix <- focal
      print(paste(focal,complexity.plant,complexity.animal,year))
      
      #---- 1.2. Import the raw data ----
      
      #setwd("~/Eco_Bayesian/Complexity_caracoles")
      home.dic <- "/home/lbuche/Eco_Bayesian/Complexity_caracoles/"
      project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/"
      
      #rm(list = ls())
      abundance.plant <- read.csv(paste0(home.dic,"data/abundance.csv"), sep=",")
      
      competition.plant <- read.csv(paste0(home.dic,"data/competition.csv"), sep=",")
      
      plant.class <- read.csv(paste0(home.dic,"data/plant_code.csv"),sep=",")
      
      herbivorypredator <- read.csv(paste0(home.dic,"data/herbivorypredator.csv"), sep=",")
      
      floral_visitor <- read.csv(paste0(home.dic,"data/floral_visitor.csv"), sep=",")
      
      groupinglevels <- read.csv(paste0(home.dic,"data/groupinglevels.csv"), sep=",")
      
      # Run UploadData.R 2. and 3. or uncomment the following lines 
      plant.class <- read.csv(paste0(home.dic,"data/plant_code.csv"), sep=",")
      
      SpData <- read.csv(paste0(home.dic,"data/SpData.csv")) 
      
      #---- 2.EXTRACT MATRIXES OF ABUNDANCES -----
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #---- 1.4. Import df of abundances----    
      # view complexity levels
      # determine the levels of complexity for plant 
      complexitylevel.plant <- levels(as.factor(plant.class[,complexity.plant]))
      if (complexity.plant != "code.plant"){complexity.minimize <- T}else{complexity.minimize <- F}
      
      FocalPrefix <- focal # "LEMA" or "HOMA" or "CHFU"
      
      # Load in the data and subset out the current focal species and complexity levels
      #SpData <- read.csv("/home/lisavm/Simulations/data/competition.csv")
      
      SpData$seed <- round(SpData$seed)
      
      SpData <- na.omit(SpData) 
      
      SpData <- SpData[which(SpData$year %in% as.numeric(year)),]
      
      SpData <- dplyr::select(SpData,allyears_of(c("day","month", "year","plot","subplot" ,"focal",
                                              "fruit","seed",complexitylevel.plant,FocalPrefix)))
      
      SpDataFocal <- SpData[which(SpData$focal == FocalPrefix),]
      if(levels(as.factor(SpDataFocal$focal)) != FocalPrefix){
        print("WRONG FOCAL")}
      
      
      if (complexity.minimize == T){
        levels.of.focal <- plant.class[which(plant.class$code.plant==FocalPrefix),complexity.plant]
        SpDataFocal[,levels.of.focal] <- SpDataFocal[,levels.of.focal] - SpDataFocal[,FocalPrefix]
      }else{
        complexitylevel.plant <- complexitylevel.plant[which(!complexitylevel.plant %in% focal )]
      }
      
      #head(SpDataFocal)
      
      # Standardized function for species abundance ATTENTION MAYBE need to be log for HTL
      standard.abudance <- function(x){
        x <- (x*5)/max(x,na.rm=T) 
        return(x)}
      
      SpDataFocal[,c(focal[1],complexitylevel.plant)] <- standard.abudance(SpDataFocal[,c(focal,complexitylevel.plant)]) 
      
      if(!max(SpDataFocal[,c(focal,complexitylevel.plant)],na.rm=T) ==5){print("Standardisation of plants UNcorrect")}
      
      
      #import higher trophic levels data 
      source(paste0(home.dic,"code/Group_FvH.R"))
      
      # add H presences and floral visitor visits
      SpData_H <- get(paste0("herbivore","_",complexity.animal))
      
      SpData_H <- SpData_H[which(SpData_H$year %in% as.integer(year) & 
                                   SpData_H$code.plant == FocalPrefix),]
      
      
      SpData_H <- SpDataFocal %>%
        dplyr::select(c("day","month","year","plot","subplot","focal")) %>%
        left_join(SpData_H)
      
      SpData_H[is.na(SpData_H)] <- 0
      
      SpData_FV <- get(paste0("floral_visitor","_",complexity.animal))
      SpData_FV <- SpData_FV[which(SpData_FV$year  %in% as.integer(year)),]
      
      SpData_FV <- SpData_FV[which(SpData_FV$year %in% as.integer(year) & 
                                     SpData_FV$code.plant == FocalPrefix),]
      
      
      SpData_FV <- SpDataFocal %>%
        dplyr::select(c("day","month","year","plot","subplot","focal")) %>%
        left_join(SpData_FV)
      

      SpData_FV[is.na(SpData_FV)] <- 0
      
      # determine the levels of complexity  animals
      complexitylevel.H <- names(SpData_H)[which(!names(SpData_H) %in% c(complexitylevel.plant,"focal","seed",
                                                                         "fruit","day.x",FocalPrefix,
                                                                         "day","month","year","plot",
                                                                         "subplot","plant",
                                                                         "code.plant"))]
      if(ncol(SpData_H[,complexitylevel.H])<1){
        print("wrong level of herbivory")
      }
      SpData_H[,complexitylevel.H] <- standard.abudance(SpData_H[,complexitylevel.H])
      SpData_H[is.na(SpData_H)] <- 0
      if(!max(SpData_H[,complexitylevel.H],na.rm=T) ==5){print("Standardisation of Herbivory UNcorrect")}
      if(sum(colSums(SpData_H[,complexitylevel.H]),na.rm=T)==0){print("No herbivory")}
      
      complexitylevel.FV <- names(SpData_FV)[which(!names(SpData_FV) %in% c(complexitylevel.plant,"seed",
                                                                            "fruit","day.x","focal",FocalPrefix,
                                                                            "day","month","year","plot",
                                                                            "subplot","plant",
                                                                            "code.plant"))]
      if(ncol(SpData_FV[,complexitylevel.FV])<1){
        print("wrong level of FV")
      }
      
      SpData_FV[,complexitylevel.FV] <- standard.abudance(SpData_FV[,complexitylevel.FV])
      SpData_FV[is.na(SpData_FV)] <- 0
      if(!max(SpData_FV[,complexitylevel.FV],na.rm=T) ==5){print("Standardisation of Floral visitors UNcorrect")}
      if(sum(colSums(SpData_FV[,complexitylevel.FV]),na.rm=T)==0){print("No floral visits")}
      
      # Next continue to extract the data needed to run the model. 
      N <- as.integer(nrow(SpDataFocal))
      Fecundity <- SpDataFocal$seed
      plot <- as.integer(factor(as.factor(SpDataFocal$plot), levels = c("1","2","3","4","5","6","7","8","9")))
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #---- 2. ABUNDANCE MATRIX ----
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #---- 2.1. Interaction (direct) matrix of plant with COMP ----
      
      # Now calculate the total number of plant species to use for the model, discounting
      #       any species columns with 0 abundance. Save a vector of the species names
      #       corresponding to each column for easy matching later.
      allyearsSpNames <- names(SpDataFocal)[!names(SpDataFocal) %in% c("day","month","year","plot",
                                                                  "subplot", "focal","fruit","seed")]
      
      allyearsSpAbunds <- SpDataFocal %>% 
        dplyr::select(allyears_of(allyearsSpNames))
      
      SpTotals <- colSums(allyearsSpAbunds)
      SpToKeep <- SpTotals > 0
      
      S <- sum(SpToKeep)
      SpMatrix <- matrix(NA, nrow = N, ncol = S)
      i <- 1
      for(s in 1:ncol(allyearsSpAbunds)){
        if(SpToKeep[s] == 1){
          SpMatrix[,i] <- allyearsSpAbunds[,s]
          i <- i + 1}else{next}
      }
      #SpMatrix <-round((SpMatrix/max(SpMatrix))*100) #scale allyears the interaction between 0 and 100
      #if(max(SpMatrix) == 100){print("scale SpMatrix_plant correct")}
      
      SpNames <- allyearsSpNames[SpToKeep]
      
      #assign(paste0("SpNames_",FocalPrefix),
      #     SpNames)
      Intra <- ifelse(SpNames == FocalPrefix, 1, 0)
      if(sum(Intra)==0){
        Intra <- ifelse(SpNames == Focal, 1, 0)
      }
      # creation of a matrix of S by S of the interaction jk in HOIs_ijk for plants
      matrix_HOIs_plant <- list()
      for (i in 1:N){
        matrix_i <- matrix(nrow=S,ncol=S)
        for (n in 1:S) {
          for (m in 1:S) {
            if (m <= n){
              matrix_i[n,m] = 0
            }else{
              matrix_i[n,m] = SpMatrix[i,n]* SpMatrix[i,m]
            } 
          }
        }
        matrix_i[is.na(matrix_i)] <- 0
        matrix_HOIs_plant[[i]] <-  matrix_i
      }
      
      SpNames_plant_HOIs <- c(apply(combn(SpNames,2),2,paste,collapse="_")) # only when j =/= of k
      
      #---- 2.2. Interaction (direct) matrix of herbivores with FOCAL ----
      
      allyearsSpAbunds_H <- SpData_H %>% 
        dplyr::select(allyears_of(complexitylevel.H))
      
      SpTotals_H <- colSums(allyearsSpAbunds_H, na.rm = T)
      SpToKeep_H  <- SpTotals_H  > 0
      H <- sum(SpToKeep_H)
      if(H==0){RemoveH = 1}else{RemoveH = 0}
      
      SpMatrix_H  <- matrix(NA, nrow = N, ncol = H)
      i <- 1
      for(h in 1:length(SpToKeep_H)){
        if(SpToKeep_H[h] == 1){
          SpMatrix_H[,i] <- allyearsSpAbunds_H[,h]
          i <- i + 1}else{next}
      }
      #SpMatrix_H <-(SpMatrix_H/max(SpMatrix_H,na.rm = T))*100 #scale allyears the interaction between 0 and 100
      #if(max(SpMatrix_H,na.rm = T) == 100){print("scale SpMatrix_H correct")}
      SpMatrix_H[is.na(SpMatrix_H)] <- 0
      SpNames_H <- names(SpToKeep_H)[SpToKeep_H]
      
      if(dim(SpMatrix_H)[2]==0){
        SpMatrix_H  <- matrix(NA, nrow = N, ncol = 1)
        SpMatrix_H[,1] <- rep(0,times=nrow(SpMatrix_H))
      }
      
      #---- 2.3. Interaction (direct) matrix of floral visitors with FOCAL ----
      
      allyearsSpAbunds_FV<- SpData_FV %>% 
        dplyr::select(allyears_of(complexitylevel.FV))
      
      SpTotals_FV <- colSums(allyearsSpAbunds_FV, na.rm = T)
      SpToKeep_FV <- SpTotals_FV  > 0
      
      FV <- sum(SpToKeep_FV)
      if(FV==0){RemoveFV = 1}else{RemoveFV = 0}
      
      
      SpMatrix_FV  <- matrix(NA, nrow = N, ncol = FV)
      i <- 1
      for(s in 1:ncol(allyearsSpAbunds_FV)){
        if(SpToKeep_FV[s] == 1){
          SpMatrix_FV[,i] <- allyearsSpAbunds_FV[,s]
          i <- i + 1}else{next}
      }
      
      if(dim(SpMatrix_FV)[2]==0){
        SpMatrix_FV  <- matrix(NA, nrow = N, ncol = 1)
        SpMatrix_FV[,1] <- rep(0,times=nrow(SpMatrix_FV))
      }
      
      #SpMatrix_FV <-(SpMatrix_FV/max(SpMatrix_FV,na.rm = T))*100 #scale allyears the interaction between 0 and 100
      #if(max(SpMatrix_FV,na.rm = T) == 100){print("scale SpMatrix_FV correct")}
      SpMatrix_FV[is.na(SpMatrix_FV)] <- 0
      SpNames_FV <- names(SpToKeep_FV)[SpToKeep_FV]
      if(length(SpNames_FV)==0){SpNames_FV <- "nopoll"}
      
      #---- 2.4. Interaction (HOIs) matrix of floral visitors with COMPETITORS ----
      
      SpData_FV_comp <- dplyr::right_join(SpData_FV,SpDataFocal,
                                          multiple = "allyears",
                                          #relationship ="many-to-many",
                                          by=c("day","month","year","plot","subplot","focal")) %>%
        unique() %>%
        gather( allyears_of(c(complexitylevel.plant,focal)),
                key = "plant.comp", value = "abundance.plant.comp") %>%
        gather(allyears_of(complexitylevel.FV),
               key = "FV", value = "visits")
      
      SpData_FV_comp$abundance.plant.comp <-as.numeric(SpData_FV_comp$abundance.plant.comp)
      SpData_FV_comp$comp.abund_visits <- SpData_FV_comp$visits*SpData_FV_comp$abundance.plant.comp
      SpData_FV_comp <-  unite(SpData_FV_comp,
                               plant.comp, FV,
                               col = "plant_FV", sep = "_") 
      
      SpData_FV_comp <- aggregate(comp.abund_visits ~ plant_FV + focal + fruit + 
                                    seed + subplot 
                                  + plot + year + month + day, 
                                  SpData_FV_comp, sum,na.action=na.omit)
      
      SpData_FV_comp <-  spread(SpData_FV_comp ,plant_FV, comp.abund_visits)
      SpData_FV_comp <- unique(SpData_FV_comp)
      
      SpData_FV_comp <- SpDataFocal %>%
        dplyr::select(c("day","month","year","plot","subplot","seed","fruit")) %>%
        left_join(SpData_FV_comp,by=c("day","month","year","plot","subplot","seed","fruit"),
                  multiple = "allyears",
                  copy=FALSE)
      
      if( !nrow(SpData_FV_comp) == N){
        print("problem in the numbr of observation, they don't match, 
        try to add  unique(SpData_FV_comp) before the previous line.")
      }
      
      
      allyearsSpAbunds_FV_comp <- SpData_FV_comp  %>% 
        dplyr::select(allyears_of(names(SpData_FV_comp)[!names(SpData_FV_comp) %in% c('subplot', "plot", "year","month",
                                                                                 "day","focal","fruit","seed",
                                                                                 focal, complexitylevel.plant)]))
      allyears.interaction.FV_comp <- as.vector(t(outer(c(complexitylevel.plant,focal), complexitylevel.FV, paste, sep="_")))
      missing.interaction.FV_comp  <- allyears.interaction.FV_comp [which(!allyears.interaction.FV_comp  %in% names(allyearsSpAbunds_FV_comp))]
      
      
      for (n in missing.interaction.FV_comp){
        allyearsSpAbunds_FV_comp[,n] <- NA
      }
      
      allyearsSpAbunds_FV_comp <- dplyr::select(allyearsSpAbunds_FV_comp,
                                           allyears_of(allyears.interaction.FV_comp)) # make the matrix in specific order
      
      
      
      #view(allyearsSpAbunds_FV_comp)
      # this part is strickly to know the potential of interaction at the before the sparsity approach
      
      SpTotals_FV_comp <- colSums(allyearsSpAbunds_FV_comp, na.rm = T)
      SpToKeep_FV_comp <- SpTotals_FV_comp  > 0
      
      FV_comp <- sum(SpToKeep_FV_comp)
      SpMatrix_FV_comp  <- matrix(NA, nrow = N, ncol = FV_comp )
      
      i <- 1
      for(s in 1:ncol(allyearsSpAbunds_FV_comp)){
        if(SpToKeep_FV_comp[s] == 1){
          SpMatrix_FV_comp[,i] <- allyearsSpAbunds_FV_comp[,s]
          i <- i + 1}else{next}
      }
      #SpMatrix_FV_comp <-(SpMatrix_FV_comp/max(SpMatrix_FV_comp,na.rm = T))*100 #scale allyears the interaction between 0 and 100
      #if(max(SpMatrix_FV_comp,na.rm = T) == 100){print("scale SpMatrix_FV_comp correct")}
      
      SpNames_FV_comp <- names(SpToKeep_FV_comp)[SpToKeep_FV_comp]
      
      vector_FV_comp <- vector() # make a vector of the position of each competitor plant species
      # important to be able to divide in the bayesian approach  
      for (n in 1:(length(complexitylevel.plant)+1)){
        vector_FV_comp[n] <- max(grep(c(complexitylevel.plant,focal)[n], SpNames_FV_comp))
        if(vector_FV_comp[n] < 0){vector_FV_comp[n] <- NA} 
      }
      vector_FV_comp <- vector_FV_comp[!is.na(vector_FV_comp)]
      
      # creation matrix_HOIs_ijf
      matrix_HOIs_ijf <- list()
      FV_comp <- ncol(SpMatrix_FV_comp) 
      for (i in 1:N){
        matrix_i <- matrix(nrow=S,ncol=FV)  
        n = 1
        m=1
        if(FV_comp ==0){
          matrix_i <- matrix(nrow=S,ncol=1, data=rep(0,times=S))}else{
            for (fv in 1:FV_comp) {
              if(fv == vector_FV_comp[n] + 1){
                n = n+1
                m = 1 }
              matrix_i[n,m] = SpMatrix_FV_comp[i,fv]
              m = m+1 
            }
          }
        matrix_i[is.na(matrix_i)]<- 0
        matrix_HOIs_ijf[[i]]  <- matrix_i
        if(sum(rowSums(is.na(matrix_i)))==T){
          print("One or multiple row(s) in matrix_HOIs_ijf contain only NAs")
        }
      }
      
      
      #---- 2.5. Interaction (HOIs) matrix of herbivores with COMPETITORS ----
      
      
      SpData_H_comp <- left_join(SpDataFocal,
                                 SpData_H,
                                 multiple = "allyears",
                                 #relationship ="many-to-many",
                                 by=c("day","month","year","plot","subplot","focal")) %>%
        unique() %>%
        gather( allyears_of(c(complexitylevel.plant,focal)),
                key = "plant.comp", value = "abundance.plant.comp") %>%
        gather(allyears_of(complexitylevel.H),
               key = "H", value = "presence")
      
      
      
      SpData_H_comp$abundance.plant.comp <-as.numeric(SpData_H_comp$abundance.plant.comp)
      SpData_H_comp$comp.abund_presence <- SpData_H_comp$presence*SpData_H_comp$abundance.plant.comp
      SpData_H_comp <-  unite(SpData_H_comp,
                              plant.comp, H,
                              col = "plant_H", sep = "_") 
      
      SpData_H_comp <- aggregate(comp.abund_presence ~ plant_H + focal + fruit + 
                                   seed + subplot 
                                 + plot + year + month + day, 
                                 SpData_H_comp, sum,na.action=na.omit)
      
      SpData_H_comp <-  spread(SpData_H_comp ,plant_H, comp.abund_presence)
      
      SpData_H_comp <- SpDataFocal %>%
        dplyr::select(c("day","month","year","plot","subplot","seed","fruit")) %>%
        left_join(SpData_H_comp,by=c("day","month","year","plot","subplot","seed","fruit"),
                  multiple = "allyears",
                  SpData_H_comp)
      
      if( !nrow(SpData_H_comp) == N){
        print("problem in the numbr of observation, they don't match, 
        try to add  unique(SpData_H_comp) before the previous line.")
      }
      
      allyearsSpAbunds_H_comp <- SpData_H_comp  %>% 
        dplyr::select(allyears_of(names(SpData_H_comp)[!names(SpData_H_comp) %in% c('subplot', "plot", "year","month",
                                                                               "day","focal","fruit","seed",
                                                                               focal, complexitylevel.plant)]))
      allyears.interaction.H_comp <- as.vector(t(outer(c(complexitylevel.plant,focal), complexitylevel.H, paste, sep="_")))
      missing.interaction.H_comp <-allyears.interaction.H_comp[which(!allyears.interaction.H_comp %in% names(allyearsSpAbunds_H_comp))]
      
      
      for (n in missing.interaction.H_comp){
        allyearsSpAbunds_H_comp[,n] <- NA
      }
      
      allyearsSpAbunds_H_comp <- dplyr::select(allyearsSpAbunds_H_comp,
                                          allyears_of(allyears.interaction.H_comp)) # make the matrix in specific order
      
      
      
      #view(allyearsSpAbunds_H_comp)
      # this part is strickly to know the potential of interaction at the before the sparsity approach
      
      SpTotals_H_comp <- colSums(allyearsSpAbunds_H_comp, na.rm = T)
      SpToKeep_H_comp <- SpTotals_H_comp  > 0
      
      H_comp <- sum(SpToKeep_H_comp)
      SpMatrix_H_comp  <- matrix(NA, nrow = N, ncol = H_comp)
      i <- 1
      for(s in 1:ncol(allyearsSpAbunds_H_comp)){
        if(SpToKeep_H_comp[s] == 1){
          SpMatrix_H_comp[,i] <- allyearsSpAbunds_H_comp[,s]
          i <- i + 1}else{next}
      }
      #SpMatrix_H_comp <-(SpMatrix_H_comp/max(SpMatrix_H_comp,na.rm = T))*5 #scale allyears the interaction between 0 and 5
      #if(max(SpMatrix_H_comp,na.rm = T) == 100){print("scale SpMatrix_H_comp correct")}
      
      SpNames_H_comp <- names(SpToKeep_H_comp)[SpToKeep_H_comp]
      
      vector_H_comp <- vector() # make a vector of the position of each competitor plant species
      # important to be able to divide in the bayesian approach  
      for (n in 1:(length(complexitylevel.plant)+1)){
        vector_H_comp[n] <- max(grep(c(complexitylevel.plant,focal)[n], SpNames_H_comp))
        if(vector_H_comp[n] < 0){vector_H_comp[n] <- NA} 
      }
      vector_H_comp <- vector_H_comp[!is.na(vector_H_comp)]
      
      # creation of  matrix_HOIs_ijh
      
      
      matrix_HOIs_ijh <- list()
      H_comp <- ncol(SpMatrix_H_comp)
      for (i in 1:N){
        matrix_i <- matrix(nrow=S,ncol=H)  
        n = 1
        m=1
        if(H_comp ==0){
          matrix_i <- matrix(nrow=S,ncol=1, data=rep(0,times=S)) }else{
            for (h_comp in 1:H_comp) {
              if(h_comp == vector_H_comp[n] + 1){
                n = n+1
                m = 1 }
              matrix_i[n,m] = SpMatrix_H_comp[i,h_comp]
              m = m+1 
            }
          }
        matrix_i[is.na(matrix_i)] <- 0
        matrix_HOIs_ijh[[i]] <-  matrix_i
        if(sum(rowSums(is.na(matrix_i)))==T){
          print("One or multiple row(s) in matrix_HOIs_ijf contain only NAs")
        }
      }
      
      
      
      
      #---- 2.6. Return object ----
      run_estimation <- 1
      # rfemove interaction of Herbivore and floral visitors vector 
      if(sum(colSums(SpMatrix_H))==0 & sum(colSums(SpMatrix_FV))==0){
        RemoveFvH <- 1
      }else{RemoveFvH <- 0
      H <- length(names(SpTotals_H[SpToKeep_H]))
      FV <- length(names(SpTotals_FV[SpToKeep_FV]))
      }
      
      
      summary.interactions.n <- tibble(focal= focal,year="allyears",
                                       complexity_plant=complexity.plant,
                                       n.competitors_plant=length(SpNames)-1, # remove focal in count
                                       competitors_plant=paste(collapse=",",
                                                               SpNames[which(!SpNames %in% focal)]),
                                       n.HOIs_plant=length(SpNames_plant_HOIs),
                                       HOIs_plant=paste(collapse=",",SpNames_plant_HOIs))
      
      if(H > 0){
        summary.interactions.n <- add_column(summary.interactions.n ,
                                             n.interactors_H=length(SpNames_H),
                                             interactors_H=paste(collapse=",",SpNames_H), n.HOIs_H=length(SpNames_H_comp),
                                             HOIs_H=paste(collapse=",",SpNames_H_comp))
      }
      
      if(FV >0){
        summary.interactions.n <- add_column(summary.interactions.n ,
                                             n.interactors_FV=length(SpNames_FV),
                                             interactors_FV=paste(collapse=",",SpNames_FV),
                                             n.HOIs_FV=length(SpNames_FV_comp),
                                             HOIs_FV=paste(SpNames_FV_comp,collapse=","))
      }
      if(FV ==0 & RemoveFV ==1){ # need to be higher than o for the model, but not estimated cause removeFV ==1
        FV <- 1
      }
      
      if(H ==0 & RemoveH ==1){ # need to be higher than o for the model, but not estimated cause removeFV ==1
        H <- 1
      }
      assign('summary.interactions.n',  summary.interactions.n,1)
      
      
      # Set the parameters defining the regularized horseshoe prior, as described in
      #       the "Incorporating sparsity-inducing priors" section of the manuscript.
      # Upper bound intrinsic fecundity
      U <- ceiling(log(mean(Fecundity)))
      
      tau0 <- 1
      slab_scale <- sqrt(2)
      slab_df <- 4
      
      DataVec <- list(N=N, 
                      SpNames = SpNames,
                      S=S,
                      U=U,
                      RemoveFvH= RemoveFvH,
                      RemoveH = RemoveH,
                      RemoveFV= RemoveFV,
                      Fecundity=Fecundity, 
                      SpMatrix=SpMatrix,
                      matrix_HOIs_plant=matrix_HOIs_plant,
                      Intra=Intra, tau0=tau0, 
                      slab_scale=slab_scale, slab_df=slab_df,
                      H=H,
                      matrix_HOIs_ijh=matrix_HOIs_ijh,
                      SpMatrix_H=SpMatrix_H,FV=FV,
                      SpMatrix_FV=SpMatrix_FV,
                      matrix_HOIs_ijf=matrix_HOIs_ijf)
      #str(DataVec) #contains 21 elements
      
      #---- 3. PRELIMINARY FIT -----
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #---- 3.1. Set up summary interactions df and parameters ---- 
      #Create summary data frame to see allyears the potential interactions
      
      
      # Set the parameters defining the regularized horseshoe prior, as described in
      #       the "Incorporating sparsity-inducing priors" section of the manuscript.
      tau0 <- 1
      slab_scale <- sqrt(2)
      slab_df <- 4
      
      #---- 4.2. Run  preliminary fit ----
      
      # Now run a preliminary fit of the model to assess parameter shrinkage
      print("preliminary fit beginning")
      options(mc.cores=parallyearsel::detectCores())
      
      rstan_options(auto_write = TRUE) 
      options(mc.cores = parallyearsel::detectCores()) # to use the core at disposition 
      
      if(run.prelimfit == T){
        PrelimFit <- stan(paste0(home.dic,"code/Short_Caracoles_BH_FH_Preliminary.stan"), 
                          data = DataVec,
                          init = "random", 
                          control =list(max_treedepth=15),
                          warmup  = 1000, 
                          iter = 2000,
                          #init_r=2,
                          chains = 8,
                          seed= 1616)
        
        save(file=paste0(project.dic,"results/stan/PrelimFit",
                         FocalPrefix,"_all_",
                         complexity.plant,"_",
                         complexity.animal,".rds"),
             PrelimFit)
      }else{
        load(paste0(project.dic,"results/stan/PrelimFit",FocalPrefix,"_all_",
                    complexity.plant,"_",complexity.animal,".rds"))
      }
      PrelimPosteriors <- rstan::extract(PrelimFit)
      print("preliminary fit done")
      #---- 3.4. Extract the inclusion indice from the Preliminary Data----
      #### Determine which parameters warrant inclusion in the final
      #       model (i.e. the data pulled their posteriors away from 0). The final model
      #       will then be run with only these species-specific parameters, but without
      #       the regularized horseshoe priors.
        SpNames.focal <- c(unlist(strsplit(summary.interactions.n$competitors_plant,",")),focal)
        if(!any(c(SpNames.focal %in% focal) == DataVec$Intra)) {
          print("wrong place of intra")
        } # check otherwise reorder! 
        SpNames <- c(unlist(strsplit(summary.interactions.n$competitors_plant,",")))
        
        Inclusion_ij <- matrix(data = 0, nrow = 1, ncol = length(SpNames.focal),
                               dimnames=list(c(""),
                                             SpNames.focal))
        beta_Inclusion_plant <- matrix(data = 0,nrow = summary.interactions.n$n.competitors_plant+1, 
                                       ncol = summary.interactions.n$n.competitors_plant+1,
                                       dimnames=list(SpNames.focal,
                                                     SpNames.focal))
        
        
        if(DataVec$RemoveH ==0){
          SpNames_H <- unlist(strsplit(summary.interactions.n$interactors_H,","))
          Inclusion_H <- matrix(data = 0,nrow = 1, 
                                ncol = summary.interactions.n$n.interactors_H,
                                dimnames=list(c(""),
                                              SpNames_H))
          beta_Inclusion_H <- matrix(data = 0, nrow = summary.interactions.n$n.competitors_plant +1, 
                                     ncol = summary.interactions.n$n.interactors_H,
                                     dimnames=list(SpNames.focal,
                                                   unlist(strsplit(summary.interactions.n$interactors_H,","))))
        }else{Inclusion_H<- matrix(data = 0,nrow = 1, ncol = 1,
                                   dimnames=list(c(""),
                                                 "NoH"))
        beta_Inclusion_H <- matrix(data = 0,nrow = summary.interactions.n$n.competitors_plant +1, 
                                   ncol =1,
                                   dimnames=list(SpNames.focal,
                                                 "NoH"))
        }
        
        if(DataVec$RemoveFV ==0){
          SpNames_FV <- unlist(strsplit(summary.interactions.n$interactors_FV,","))
          Inclusion_FV<- matrix(data = 0,nrow = 1, ncol = summary.interactions.n$n.interactors_FV,
                                dimnames=list(c(""),
                                              SpNames_FV))
          beta_Inclusion_FV <- matrix(data = 0, nrow = summary.interactions.n$n.competitors_plant+1, 
                                      ncol = summary.interactions.n$n.interactors_FV,
                                      dimnames=list(c(unlist(strsplit(summary.interactions.n$competitors_plant,",")),focal),
                                                    unlist(strsplit(summary.interactions.n$interactors_FV,","))))
          
        }else{Inclusion_FV<- matrix(data = 0,nrow = 1, ncol = 1,
                                    dimnames=list(c(""),
                                                  "NoPoll"))
        beta_Inclusion_FV <- matrix(data = 0, nrow = summary.interactions.n$n.competitors_plant+1, 
                                    ncol = 1,
                                    dimnames=list(c(unlist(strsplit(summary.interactions.n$competitors_plant,",")),focal),
                                                  "NoPoll"))
        }
        
        IntLevel <- 0.2 #0.5 usuallyearsy, 0.75 for Waitzia, shade
        
        #str(PrelimPosteriors$beta_hat_ijk)
        #str(PrelimPosteriors$alpha_hat_ij)
        #For plants as second actor
        for(s in 1:length(SpNames.focal)){
          # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
          Ints_ij <- HDInterval::hdi(PrelimPosteriors$alpha_hat_ij[,s], credMass = IntLevel)
          if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
            Inclusion_ij[1,s] <- 1
          }
        }
        for(s in 1:length(SpNames.focal)){  
          for(m in 1:DataVec$S){
            beta_Ints_ijk <- HDInterval::hdi(PrelimPosteriors$beta_plant_hat_ijk[,s,m], credMass = IntLevel)
            
            if(beta_Ints_ijk[1] > 0 | beta_Ints_ijk[2] < 0){
              beta_Inclusion_plant[s,m] <- 1
            }
          }
          if(DataVec$RemoveH == 0){
            for(m in 1:length(SpNames_H)){
              beta_Ints_ijh <- HDInterval::hdi(PrelimPosteriors$beta_H_hat_ijh[,s,m], credMass = IntLevel)
              
              if(is.na(beta_Ints_ijh[1]) | is.na(beta_Ints_ijh[2])) next
              if(beta_Ints_ijh[1] > 0 | beta_Ints_ijh[2] < 0){
                beta_Inclusion_H[s,m] <- 1
              }
            }
          }
          if(DataVec$RemoveFV ==0){for(m in 1:length(SpNames_FV)){
            beta_Ints_ijf <- HDInterval::hdi(PrelimPosteriors$beta_FV_hat_ijf[,s,m], credMass = IntLevel)
            if(is.na(beta_Ints_ijf[1]) | is.na(beta_Ints_ijf[2])) next
            
            if(beta_Ints_ijf[1] > 0 | beta_Ints_ijf[2] < 0){
              beta_Inclusion_FV[s,m] <- 1
            }
          }
          }
        }
        
        #For H as second actor
        if(DataVec$RemoveH ==0){
          for(s in 1:length(SpNames_H)){
            # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
            Ints_ih <- HDInterval::hdi(PrelimPosteriors$gamma_H_hat_ih[,s], credMass = IntLevel)
            if(is.na(Ints_ih[1]) | is.na(Ints_ih[2])) next
            
            if(Ints_ih[1] > 0 | Ints_ih[2] < 0){
              Inclusion_H[1,s] <- 1
            }
          }
        }
        
        #For pollinator as second actor
        if(DataVec$RemoveFV ==0){
          for(s in 1:length(SpNames_FV)){
            # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
            Ints_if <- HDInterval::hdi(PrelimPosteriors$gamma_FV_hat_if[,s], credMass = IntLevel)
            if(is.na(Ints_if[1]) | is.na(Ints_if[2])) next
            
            if(Ints_if[1] > 0 | Ints_if[2] < 0){
              Inclusion_FV[1,s] <- 1
            }
          }
        }
        
        
        summary.interactions.n$n.competitors_plant_inclus <- sum(Inclusion_ij)
        if (sum(Inclusion_ij) > 0){
          summary.interactions.n$competitors_plant_inclus <- paste(c(colnames(Inclusion_ij))[which(Inclusion_ij>=1)],collapse=",")
        }
        
        if(DataVec$RemoveFV ==0){
          summary.interactions.n$n.competitors_FV_inclus <- sum(Inclusion_FV)
          if (sum(Inclusion_FV)>0){
            summary.interactions.n$competitors_FV_inclus <- paste(c(colnames(Inclusion_FV))[which(Inclusion_FV >=1)],collapse=",")
          }
        }
        if(DataVec$RemoveH ==0){
          summary.interactions.n$n.competitors_H_inclus <- sum(Inclusion_H)
          if (sum(Inclusion_H)>0){
            summary.interactions.n$competitors_H_inclus <- paste(c(colnames(Inclusion_H))[which(Inclusion_H >=1)],collapse=",")
          }
        }
        summary.interactions.n <- as_tibble(summary.interactions.n)
        HOIs.potential <- c("plant","H","FV")
        if(DataVec$RemoveH == 1){HOIs.potential <- c("plant","FV")}
        if(DataVec$RemoveFV == 1){HOIs.potential <- c("plant","H")}
        for ( x in HOIs.potential ){
          matrix_x <- get(paste0("beta_Inclusion_",x))
          summary.interactions.n[,paste("n.HOIs",x,"inclus",sep="_")] = sum(matrix_x)
          
          if (sum(matrix_x  > 0)){
            summary.interactions.n[,paste0("HOIs_",x,"_inclus")]  = paste(paste0(rownames(  matrix_x )[which(  matrix_x >=1, arr.in=TRUE)[,1]],sep="_",
                                                                                 colnames(matrix_x )[which(  matrix_x >=1, arr.in=TRUE)[,2]]),collapse=",")
          } 
        }
        
        Inclusion_allyears <- list(interaction = summary.interactions.n,
                              Inclusion_ij=Inclusion_ij,
                              Inclusion_FV=Inclusion_FV,
                              Inclusion_H=Inclusion_H,
                              beta_Inclusion_plant=beta_Inclusion_plant,
                              beta_Inclusion_FV=beta_Inclusion_FV,
                              beta_Inclusion_H=beta_Inclusion_H)
        Inclusion_allyears <- append(Inclusion_allyears,DataVec)
        
        save(Inclusion_allyears,
             file=paste0(home.dic,"results/Inclusion",FocalPrefix,"_all_",complexity.plant,"_",complexity.animal,".RData"))
        
        summary.interactions <- bind_rows(summary.interactions,summary.interactions.n)
        
        write.csv(summary.interactions,paste0(home.dic,"results/summary.interactions.csv"))

        DataVec.final <- list(Inclusion_ij = Inclusion_ij,
                              Inclusion_FV=Inclusion_FV,
                              Inclusion_H=Inclusion_H, 
                              beta_Inclusion_plant=beta_Inclusion_plant,
                              beta_Inclusion_FV=beta_Inclusion_FV,
                              beta_Inclusion_H=beta_Inclusion_H,
                              run_estimation = 1)
        
        assign("DataVec.final",DataVec.final,1)

      
      #---- 3.3. Preliminary fit posterior check and behavior checks---- 
      ##### Diagnostic plots
      pdf(paste0(home.dic,"figure/Prelimfit_",paste(FocalPrefix,"allyearsyears",complexity.plant,complexity.animal,sep="_"),".pdf"))
      title= paste0("Diagnostic plots for Prelimfit of ",focal, " in allyears years \nwhen plant complexity is ",complexity.plant,
                    " \nand HTL complexity is ",complexity.animal)
      # Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
      source("~/Eco_Bayesian/Test_simulation/code/stan_modelcheck_rem.R") # callyears the functions to check diagnistic plots
      # check the distribution of Rhats and effective sample sizes 
      stan_post_pred_check(PrelimPosteriors,"F_hat",DataVec$Fecundity, max(DataVec$Fecundity)+20)
      
      # N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
      hist(summary(PrelimFit)$summary[,"Rhat"], 
           main = paste("Prelim fit: Histogram of Rhat for",
                        FocalPrefix, " in allyears years \nwhen plant complexity is ",complexity.plant,
                        " \nand HTL complexity is ",complexity.animal))
      hist(summary(PrelimFit)$summary[,"n_eff"],
           main = paste("Prelim fit: Histogram of n_eff for",
                        FocalPrefix, " in allyears years \nwhen plant complexity is ",complexity.plant,
                        " \nand HTL complexity is ",complexity.animal))
      
      # Next check the correlation among key model parameters and identify any
      #       divergent transitions
      par <- c('lambdas','disp_dev','alpha_generic','alpha_intra',"gamma_FV_generic","gamma_H_generic")
      if(DataVec$RemoveFV == 1){par <- c('lambdas','disp_dev','alpha_generic',
                                         'alpha_intra',"gamma_H_generic")
      }
      if(DataVec$RemoveH == 1){par <- c('lambdas','disp_dev',
                                        'alpha_generic','alpha_intra',"gamma_FV_generic")
      }
      if(DataVec$RemoveH == 1 & DataVec$RemoveFV == 1){par <- c('lambdas','disp_dev',
                                                                'alpha_generic','alpha_intra')
      }
      
      stan_trace(PrelimFit,pars=par,inc_warmup = TRUE)
      stan_dens(PrelimFit,pars=par)
      stan_plot(PrelimFit,pars=par)
      
      dev.off()
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #---- 4. FINAL FIT----
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #---- 4.1. Set up parameters ---- 
      
      DataVec.final <- append(DataVec,DataVec.final)
      
      #FinalFit <- stan(file =  "/home/lisavm/Simulations/code/Caracoles_BH_Final.stan", data = DataVec, iter = 1000, chains = 2)
      
      
      #---- 4.2. Run final fit ---- 
      print("final fit begins")
      rstan_options(auto_write = TRUE) 
      options(mc.cores=parallyearsel::detectCores())
      set.seed(1616)
      std.error <- function(x) sd(x)/sqrt(length(x))
      
      list.init <- function(...)list(lambdas= array(as.numeric(rnorm(1,mean = mean(PrelimPosteriors$lambdas),
                                                                     sd= std.error(PrelimPosteriors$lambdas))), dim = 1),
                                     alpha_generic_tilde = array(as.numeric(rnorm(1,mean = mean(PrelimPosteriors$alpha_generic_tilde),
                                                                                  sd= std.error(PrelimPosteriors$alpha_generic_tilde))), dim = 1),
                                     alpha_intra_tilde = array(as.numeric(rnorm(1,mean = mean(PrelimPosteriors$alpha_intra_tilde),
                                                                                sd= std.error(PrelimPosteriors$alpha_intra_tilde))), dim = 1),
                                     disp_dev = array(as.numeric(rnorm(1,mean = mean(PrelimPosteriors$disp_dev),
                                                                       sd= std.error(PrelimPosteriors$disp_dev))), dim = 1))
      
      if(run.finalfit == T){
        options(mc.cores = parallyearsel::detectCores()) # to use the core at disposition 
        FinalFit <- stan(file = paste0(home.dic,"code/Short_Caracoles_BH_Final.stan") ,
                         #fit= PrelimFit, 
                         data = DataVec.final,
                         init =list.init, # allyears initial values are 0 
                         control=list(max_treedepth=15),
                         warmup = 500,
                         iter = 1000, 
                         init_r = 2,
                         chains = 8,
                         seed= 1616) 
        
        save(file=paste0(project.dic,"results/FinalFit",
                         FocalPrefix,"_all_",complexity.plant,
                         "_",complexity.animal,".rds"),
             FinalFit)
      }else{
        load(paste0(project.dic,"results/FinalFit",FocalPrefix,"_all_",complexity.plant,"_",complexity.animal,".rds")) 
      }
      FinalPosteriors <- rstan::extract(FinalFit)
      
      print("final fit is done")
      #---- 4.3. Final fit posterior check and behavior checks---- 
      
      ##### Diagnostic plots and post prediction 
      pdf(paste0(home.dic,"figure/FinalFit_",paste(FocalPrefix,"allyearsyears",complexity.plant,complexity.animal,sep="_"),".pdf"))
      # Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
      source("~/Eco_Bayesian/Test_simulation/code/stan_modelcheck_rem.R") # callyears the functions to check diagnistic plots
      # check the distribution of Rhats and effective sample sizes 
      ##### Posterior check
      
      stan_post_pred_check(FinalPosteriors,"F_hat",
                           DataVec$Fecundity,max(DataVec$Fecundity)+20) 
      
      # N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
      hist(summary(FinalFit)$summary[,"Rhat"],
           main = paste("Finat Fit: Histogram of Rhat for",
                        FocalPrefix, " in allyears years \nwhen plant complexity is ",complexity.plant,
                        " \nand HTL complexity is ",complexity.animal))
      hist(summary(FinalFit)$summary[,"n_eff"],
           main = paste("Finat Fit: Histogram of Rhat for",
                        FocalPrefix, " in allyears years \nwhen plant complexity is ",complexity.plant,
                        " \nand HTL complexity is ",complexity.animal))
      
      # plot the corresponding graphs
      par <- c('lambdas','disp_dev','alpha_generic','alpha_intra',"gamma_FV_generic","gamma_H_generic")
      if(DataVec$RemoveFV == 1){par <- c('lambdas','disp_dev','alpha_generic',
                                         'alpha_intra',"gamma_H_generic")
      }
      if(DataVec$RemoveH == 1){par <- c('lambdas','disp_dev',
                                        'alpha_generic','alpha_intra',"gamma_FV_generic")
      }
      if(DataVec$RemoveH == 1 & DataVec$RemoveFV == 1){par <- c('lambdas','disp_dev',
                                                                'alpha_generic','alpha_intra')
      }
      stan_trace(FinalFit , pars=par,
                 inc_warmup = TRUE)
      stan_dens(FinalFit , pars=par)
      stan_plot(FinalFit , pars=par)
      # Next check the correlation among key model parameters and identify any
      dev.off()
      
      #---- 4.4. Final fit analysis----
      
      #---- Generic parameters---
      df_parameter_n <-NULL
      generic_name <- c("lambdas","alpha_intra","alpha_generic","gamma_H_generic","gamma_FV_generic")
      generic_name <- generic_name[generic_name %in%  names(FinalPosteriors)]
      for ( n in generic_name ){
        df_parameter_hat_n <- FinalPosteriors[[n]] %>% 
          as.data.frame()
        names( df_parameter_hat_n) <-  n 
        
        df_parameter_n <- bind_cols(df_parameter_n,df_parameter_hat_n)
      }
      
      df_parameter_n$focal <- focal
      df_parameter_n$year <- "allyears"
      df_parameter_n$complexity.plant <- complexity.plant
      df_parameter_n$complexity.animal <- complexity.animal
      
      write.csv(df_parameter_n,
                paste0(project.dic,"results/parameters/Parameters_",focal,"_all_",complexity.plant,
                       complexity.animal,"_FinalFit.csv"))
      
      # Species specific parameters
      
      
      load(paste0(home.dic,"results/Inclusion",focal,"_all_",complexity.plant,"_",complexity.animal,".RData"))
      assign(paste0("Inclusion",focal,"_all_",complexity.plant,"_",complexity.animal),Inclusion_allyears)
      
      rm(Inclusion_allyears)
      
      df_inclusion_nat <- as.data.frame(get(paste0("Inclusion",focal,"_all_",complexity.plant,"_",complexity.animal))$interaction) 
      
      
      n.inclus <- c("n.competitors_plant_inclus","n.competitors_FV_inclus","n.competitors_H_inclus",
                    "n.HOIs_plant_inclus","n.HOIs_H_inclus",
                    "n.HOIs_FV_inclus")
      inclus <- c("competitors_plant_inclus","competitors_FV_inclus","competitors_H_inclus",
                  "HOIs_plant_inclus","HOIs_H_inclus",
                  "HOIs_FV_inclus")
      inclus.short <- grep(pattern = '.*inclus$', x = colnames(df_inclusion_nat), value = T)
      inclus.short <- inclus.short[which(!inclus.short %in% n.inclus)]
      inclus.not.short <- inclus[which(!inclus %in% inclus.short)]
      df_inclusion_nat[,inclus.not.short] <-NA
      
      n.potential <- c("n.competitors_plant","n.interactors_FV","n.interactors_H",
                       "n.HOIs_plant","n.HOIs_H",
                       "n.HOIs_FV")
      
      potential <- c("competitors_plant","interactors_FV","interactors_H",
                     "HOIs_plant","HOIs_H",
                     "HOIs_FV")
      if(DataVec$RemoveFV == 1){
        n.inclus <- c("n.competitors_plant_inclus","n.competitors_H_inclus",
                      "n.HOIs_plant_inclus","n.HOIs_H_inclus")
        inclus <- c("competitors_plant_inclus","competitors_H_inclus",
                    "HOIs_plant_inclus","HOIs_H_inclus")
        inclus.short <- grep(pattern = '.*inclus$', x = colnames(df_inclusion_nat), value = T)
        inclus.short <- inclus.short[which(!inclus.short %in% n.inclus)]
        inclus.not.short <- inclus[which(!inclus %in% inclus.short)]
        df_inclusion_nat[,inclus.not.short] <-NA
        
        n.potential <- c("n.competitors_plant","n.interactors_H",
                         "n.HOIs_plant","n.HOIs_H")
        
        potential <- c("competitors_plant","interactors_H",
                       "HOIs_plant","HOIs_H")
        
      }
      if(DataVec$RemoveH == 1){
        n.inclus <- c("n.competitors_plant_inclus","n.competitors_FV_inclus","n.HOIs_plant_inclus","n.HOIs_FV_inclus")
        inclus <- c("competitors_plant_inclus","competitors_FV_inclus","HOIs_plant_inclus", "HOIs_FV_inclus")
        inclus.short <- grep(pattern = '.*inclus$', x = colnames(df_inclusion_nat), value = T)
        inclus.short <- inclus.short[which(!inclus.short %in% n.inclus)]
        inclus.not.short <- inclus[which(!inclus %in% inclus.short)]
        df_inclusion_nat[,inclus.not.short] <-NA
        
        n.potential <- c("n.competitors_plant","n.interactors_FV", "n.HOIs_plant","n.HOIs_FV")
        
        potential <- c("competitors_plant","interactors_FV","HOIs_plant","HOIs_FV")
      }
      if(DataVec$RemoveFvH == 1){
        n.inclus <- c("n.competitors_plant_inclus","n.HOIs_plant_inclus")
        inclus <- c("competitors_plant_inclus","HOIs_plant_inclus")
        inclus.short <- grep(pattern = '.*inclus$', x = colnames(df_inclusion_nat), value = T)
        inclus.short <- inclus.short[which(!inclus.short %in% n.inclus)]
        inclus.not.short <- inclus[which(!inclus %in% inclus.short)]
        df_inclusion_nat[,inclus.not.short] <-NA
        
        n.potential <- c("n.competitors_plant","n.HOIs_plant")
        
        potential <- c("competitors_plant","HOIs_plant")
        
      }
      
      extract.vert <- function(mat,vec,keyname,value,vecname){
        mat.2 <- mat %>% 
          gather(any_of(vec), 
                 key = keyname , value = value) %>%
          mutate(keyname = rep(vecname, each = nrow(mat)))
        return(mat.2)
      }
      
      if (sum(colSums(df_inclusion_nat[,grep(pattern = "^n.*\\inclus$", x = colnames(df_inclusion_nat), value = T)],na.rm=T)) > 0) {
        
        df_inclusion_nat <- df_inclusion_nat %>% 
          gather(any_of(potential), 
                 key="interaction", value= "potential.interactor") %>%
          #mutate(interaction = rep(inclus, times = nrow(df_inclusion_nat))) %>%
          dplyr::select(-any_of(c(inclus,n.inclus,n.potential))) %>%
          mutate(interactor = extract.vert(df_inclusion_nat,inclus,
                                           "interaction","identity",inclus)[,"value"],
                 n.number = extract.vert(df_inclusion_nat,n.inclus,
                                         "interaction.number","n.number",n.inclus)[,"value"],
                 n.potential.interactor = extract.vert(df_inclusion_nat,n.potential ,
                                                       "interaction.pot","n.identity.potential",n.potential)[,"value"]) %>%
          mutate(interaction = dplyr::case_when(interaction == "competitors_plant" ~ "alpha_hat_ij",
                                                interaction == "interactors_FV" ~ "gamma_FV_hat_if",
                                                interaction== "interactors_H" ~ "gamma_H_hat_ih",
                                                interaction== "HOIs_plant" ~ "beta_plant_hat_ijk",
                                                interaction== "HOIs_FV" ~ "beta_FV_hat_ijf",
                                                interaction== "HOIs_H" ~ "beta_H_hat_ijh"))
        df_inclusion_nat_relevant <- df_inclusion_nat %>% 
          dplyr::filter(n.number != 0)
        
        # function to extract specific value of a HOIs matrix 
        extract.matrix.HOIs.plant <- function(n.neighbours,position.n) {
          vec.position <- 1:sum(1:n.neighbours-1)
          data.position <- c(rep(0,times=vec.position[1]),vec.position[1]:c(n.neighbours-1))
          for(n in c(2:n.neighbours)){
            v <- c(rep(0,times=vec.position[n]),
                   data.position[length(data.position)]+1:vec.position[length(vec.position)])[1:n.neighbours]
            if(n == n.neighbours){v <- c(rep(0,times=vec.position[n]))}
            data.position <- c(data.position,v)
          }
          mat.position <- matrix(nrow=n.neighbours,ncol=n.neighbours, byrow = T,
                                 data= data.position)
          position.new <- c()
          for( n in 1:length(position.n)){
            
            vec.position <- 1:n.neighbours
            
            position.new[n]<- paste0(vec.position[apply(mat.position, 1, function(r) any(r %in% c(position.n[n])))],",",# row
                                     vec.position[apply(mat.position, 2, function(r) any(r %in% c(position.n[n])))]) # column
          }
          return(position.new)
        }
        
        extract.matrix.HOIs.HTL <- function(n.plants,n.neighbours,position.n) {
          mat.position <- matrix(ncol=n.neighbours,nrow=n.plants,
                                 data=c(1:(n.plants*n.neighbours)),
                                 byrow = T)
          
          position.new <- c()
          for( n in 1:length(position.n)){
            
            row.position <- 1:n.plants
            col.position <- 1:n.neighbours
            
            position.new[n]<- paste0(row.position[apply(mat.position, 1, function(r) any(r %in% c(position.n[n])))],",",# row
                                     col.position[apply(mat.position, 2, function(r) any(r %in% c(position.n[n])))]) # column
          }
          return(position.new)
        }
        
        if(nrow(df_inclusion_nat_relevant) != 0 ){
          df_parameter_hat <- NULL
          for( i in 1:nrow(df_inclusion_nat_relevant)){
            position.n <- match(stringr::str_split_1(df_inclusion_nat_relevant[i,"interactor"],","),
                                stringr::str_split_1(df_inclusion_nat_relevant[i,"potential.interactor"],","))
            name.allyears <- stringr::str_split_1(df_inclusion_nat_relevant[i,"potential.interactor"],",")
            if(df_inclusion_nat_relevant[i,"interaction"] == c("beta_plant_hat_ijk")){
              n.plants <- DataVec$S
              n.neighbours <- DataVec$S
              position.name <- position.n
              position.n <- extract.matrix.HOIs.plant(n.neighbours, position.n)
            }
            if(df_inclusion_nat_relevant[i,"interaction"] == c("beta_FV_hat_ijf")){
              name.allyears <-paste0(rep(c(stringr::str_split_1(df_inclusion_nat[which(df_inclusion_nat$interaction == "alpha_hat_ij"),
                                                                            "potential.interactor"],","),
                                      focal),each=DataVec$FV),
                                "_",
                                stringr::str_split_1(df_inclusion_nat[which(df_inclusion_nat$interaction == "gamma_FV_hat_if"),
                                                                      "potential.interactor"],","))
              position.n <- match(stringr::str_split_1(df_inclusion_nat_relevant[i,"interactor"],","),
                                  name.allyears)
              name.allyears <- stringr::str_split_1(df_inclusion_nat_relevant[i,"potential.interactor"],",")
              n.plants <- DataVec$S
              n.neighbours <- DataVec$FV
              position.name <- position.n
              position.n <- extract.matrix.HOIs.HTL(n.plants,n.neighbours, position.n)
            }
            if(df_inclusion_nat_relevant[i,"interaction"] == c("beta_H_hat_ijh")){
              name.allyears <- paste0(rep(c(stringr::str_split_1(df_inclusion_nat[which(df_inclusion_nat$interaction == "alpha_hat_ij"),
                                                                             "potential.interactor"],","),
                                       focal),each=DataVec$H),
                                 "_",
                                 stringr::str_split_1(df_inclusion_nat[which(df_inclusion_nat$interaction == "gamma_H_hat_ih"),
                                                                       "potential.interactor"],","))
              
              position.n <- match(stringr::str_split_1(df_inclusion_nat_relevant[i,"interactor"],","),
                                  name.allyears)
              n.plants <- DataVec$S
              n.neighbours <- DataVec$H
              position.name <- position.n
              position.n <- extract.matrix.HOIs.HTL(n.plants,n.neighbours, position.n)
            }
            df_parameter_hat_n <- NULL
            for( n in position.n){
              n.2 <- which(position.n == n)
              if(is.numeric(n)){
                df_parameter_hat_n <- FinalPosteriors[[df_inclusion_nat_relevant[i,"interaction"]]][,n] %>% 
                  as.data.frame()
                names(df_parameter_hat_n) <- paste0(df_inclusion_nat_relevant[i,"interaction"],"_",name.allyears[n])
                
              }else{
                n.list <- stringr::str_split_1(n,",")
                df_parameter_hat_n <- FinalPosteriors[[df_inclusion_nat_relevant[i,"interaction"]]][,as.numeric(n.list[1]),as.numeric(n.list[2])] %>% 
                  as.data.frame()
                names(df_parameter_hat_n) <- paste0(df_inclusion_nat_relevant[i,"interaction"],"_", name.allyears[position.name[n.2]])
              }
              df_parameter_hat <- bind_cols(df_parameter_hat, df_parameter_hat_n) 
              
            }
            
          }
          
          df_parameter_hat$focal <- focal
          df_parameter_hat$year <- "allyears"
          df_parameter_hat$complexity.plant <- complexity.plant
          df_parameter_hat$complexity.animal <- complexity.animal
          write.csv(df_parameter_hat,
                    paste0(project.dic,"results/parameters/Parameters_hat_",focal,"_all_",complexity.plant,
                           complexity.animal,"_FinalFit.csv"))
    
      
      print(paste("allyears done for",FocalPrefix,"allyears years",
                  complexity.plant,complexity.animal,sep=" "))
        }
      }
    
    }
}
write.csv(summary.interactions,paste0(home.dic,"results/summary.interactions.csv"))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 3. Analysis: Inclusion df, figures----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 3.0 CSV for inclusion----

df_inclusion_nat_allyears <- data.frame()
for( focal in c("LEMA","CHFU","HOMA","CETE")){ # "CHFU","HOMA",
    for( complexity.level in c(1:3)){
      
      complexity.animal <- c("group","family","species")[complexity.level]
      complexity.plant <-c("class","family","code.plant")[complexity.level]
      
      load(paste0(home.dic,"results/Inclusion",focal,"_all_",complexity.plant,"_",complexity.animal,".RData"))
      
      assign(paste0("Inclusion",focal,"_all_",complexity.plant,"_",complexity.animal),Inclusion_all)
      rm(Inclusion_all)
      
      df_inclusion_n <- as.data.frame(get(paste0("Inclusion",focal,"_all_",complexity.plant,"_",complexity.animal))$interaction) 
      
      df_inclusion_n$focal <- focal
      df_inclusion_n$year <- as.numeric(year)
      df_inclusion_n$complexity.animal <- complexity.animal
      df_inclusion_nat_allyears <- bind_rows(df_inclusion_nat_allyears,    df_inclusion_n)
      
    }
}
write.csv(df_inclusion_nat_allyears,
          file = paste0(home.dic,"results/Chapt1_Inclusion_parameters_allyears.csv"))

#n.inclus <- grep(pattern = '^n.*inclus$', x = colnames(df_inclusion_nat), value = T)

n.inclus <- c("n.competitors_plant_inclus","n.competitors_FV_inclus","n.competitors_H_inclus",
              "n.HOIs_plant_inclus","n.HOIs_H_inclus",
              "n.HOIs_FV_inclus")
inclus <- c("competitors_plant_inclus","competitors_FV_inclus","competitors_H_inclus",
            "HOIs_plant_inclus","HOIs_H_inclus",
            "HOIs_FV_inclus")
inclus.short <- grep(pattern = '.*inclus$', x = colnames(df_inclusion_nat_allyears), value = T)
inclus.short <- inclus.short[which(!inclus.short %in% n.inclus)]
inclus.not.short <- inclus[which(!inclus %in% inclus.short)]
df_inclusion_nat[,inclus.not.short] <-NA
n.potential <- c("n.competitors_plant","n.interactors_FV","n.interactors_H",
                 "n.HOIs_plant","n.HOIs_H",
                 "n.HOIs_FV")

potential <- c("competitors_plant","interactors_FV","interactors_H",
               "HOIs_plant","HOIs_H",
               "HOIs_FV")
df_inclusion_nat_all_short <- df_inclusion_nat_allyears[which(rowSums(df_inclusion_nat_allyears[,n.inclus],na.rm = T)>0),]
view(df_inclusion_nat_all_short)
write.csv(df_inclusion_nat_all_short,
          file = paste0(home.dic,"results/Chapt1_Inclusion_parameters_short_allyears.csv"))


extract.vert <- function(mat,vec,keyname,value,vecname){
  mat.2 <- mat %>% 
    gather(all_of(vec), 
           key = keyname , value = value) %>%
    mutate(keyname = rep(vecname, each = nrow(mat)))
  return(mat.2)
}

df_inclusion_nat_vertical_allyears <- df_inclusion_nat_allyears %>% 
  gather(all_of(potential), 
         key="interaction", value= "potential.interactor") %>%
  #mutate(interaction = rep(inclus, times = nrow(df_inclusion_nat))) %>%
  select(-all_of(c(inclus,n.inclus,n.potential))) %>%
  mutate(interactor = extract.vert(df_inclusion_nat_allyears,inclus,
                                   "interaction","identity",inclus)[,"value"],
         n.number = extract.vert(df_inclusion_nat_allyears,n.inclus,
                                 "interaction.number","n.number",n.inclus)[,"value"],
         n.potential.interactor = extract.vert(df_inclusion_nat_allyears,n.potential ,
                                               "interaction.pot","n.identity.potential",n.potential)[,"value"]) %>%
  mutate(interaction = case_when(interaction == "competitors_plant" ~ "Plant - plant",
                                 interaction == "HOIs_FV" ~ "Plant - plant - floral visitor",
                                 interaction== "HOIs_H" ~ "Plant - plant - herbivore",
                                 interaction== "HOIs_plant" ~ "Plant - plant - plant",
                                 interaction== "interactors_FV" ~ "Plant - floral visitor",
                                 interaction== "interactors_H" ~ "Plant - herbivore"))



#---- 3.2 Figure of pairwise interactions distributions and Intra Vs Inter ----
#extraction generic parameters
df_param_allyears <- NULL

for( focal in c("CETE","LEMA","HOMA","CHFU")){ # "CHFU","HOMA","CETE"
    for( complexity.level in c(1:3)){
      
      complexity.animal <- c("group","family","species")[complexity.level]
      complexity.plant <-c("class","family","code.plant")[complexity.level]
      
      df_alpha_n <- read.csv(paste0(project.dic,"results/parameters/Parameters_",
                                    focal,"_all_",complexity.plant,
                                    complexity.animal,"_FinalFit.csv"))
      
      if(is.null(df_param_allyears)){df_param_allyears  <- df_alpha_n
      }else{
        df_param_allyears <- full_join(df_param_allyears,df_alpha_n)
      }
  }
}

# species specific parameters

inclusion.interaction.type_allyears <- unique(df_inclusion_nat_vertical_allyears[,c("interaction","interactor")]) %>%
  separate(interactor, into=paste0("interactor", 1:max(df_inclusion_nat_vertical_allyears$n.number,na.rm=T)), sep=",")  %>%
  gather(paste0("interactor", 1:max(df_inclusion_nat_vertical_allyears$n.number,na.rm=T)),key="to.remove",value="interactor") %>%
  filter(!is.na(interactor)) %>%
  select(interaction,interactor)


df_param_hat_allyears <- NULL
for ( i in 1:nrow(df_inclusion_nat_all_short)){
  focal <- df_inclusion_nat_all_short[i,'focal']
  
  complexity.plant <- df_inclusion_nat_all_short[i,'complexity_plant']
  complexity.animal <- df_inclusion_nat_all_short[i,'complexity.animal']
  
  df_param_hat_n <- read.csv(paste0(project.dic,"results/parameters/Parameters_hat_",
                                    focal,"_all_",complexity.plant,
                                    complexity.animal,"_FinalFit.csv"))
  name.species.specific <- names(df_param_hat_n)[which(!names(df_param_hat_n) %in% c("X","focal","year",
                                                                                     "complexity.plant",
                                                                                     "complexity.animal" ))]
  
  name.species.specific_n <- sapply(  strsplit(name.species.specific, '_'), `[`, 5)
  for(n in which(is.na(name.species.specific_n))){
    name.species.specific_n[n] <- sapply(  strsplit(name.species.specific, '_'), `[`, 4)[n]
    if(!is.na(sapply(strsplit(name.species.specific, '_'), `[`, 6)[n])){
      name.species.specific_n[n] <- paste0(sapply(  strsplit(name.species.specific, '_'), `[`, 4)[n],"_",
                                           sapply(strsplit(name.species.specific, '_'), `[`, 5)[n])
      
    }
  }
  
  for (n in which(!is.na(sapply(strsplit(name.species.specific, '_'), `[`, 6)))){ # in case of HOIs we need to put bacj togetehr the two species together in the name of interactors
    name.species.specific_n[n] <- paste0(sapply(  strsplit(name.species.specific, '_'), `[`, 5)[n],"_",
                                         sapply(strsplit(name.species.specific, '_'), `[`, 6)[n])
  }
  


gr <- colorRampPalette(c("#009E73"))(200)                      
re <- colorRampPalette(c("#E69F00"))(200)

df_param_vert_allyears <- df_param_allyears %>%
  gather(alpha_intra,alpha_generic,
         gamma_H_generic,gamma_FV_generic,
         key="parameter",value="estimate") %>%
  filter(!is.na(estimate)) %>%
  mutate(parameter = case_when(parameter=="gamma_H_generic" ~ "Plant - herbivore",
                               parameter=="gamma_FV_generic" ~ "Plant - floral visitor",
                               parameter=="alpha_intra" ~ "Intraspecific",
                               parameter=="alpha_generic" ~ "Plant - plant"))

write.csv(df_param_allyears,
          file = paste0(home.dic,"results/Chapt1_Parameters_values_allyears.csv"))


df_param_vert_allyears <- df_param_vert_allyears %>%
  mutate(parameter = factor(parameter,
                            levels = rev(c("Plant - plant", "Intraspecific",
                                           "Plant - floral visitor","Plant - herbivore"
                            )))) %>%
  mutate(complexity.plant  = factor(complexity.plant,
                                    levels = c("class","family","code.plant")),
         complexity.plant  = case_when(complexity.plant=="class"~ "1.Order",
                                       complexity.plant=="family"~ "2.Family",
                                       complexity.plant=="code.plant"~ "3.Species")) %>%
  mutate(focal.name = case_when(focal == "CHFU" ~"Chamaemelum \nfuscatum",
                                focal == "CETE" ~"Centaurium \ntenuiflorum",
                                focal == "HOMA" ~"Hordeum \nmarinum",
                                focal == "LEMA" ~"Leontodon \nmaroccanus")) %>%
  mutate(parameter = fct_relevel(parameter,rev(c("Plant - plant", "Intraspecific",
                                 "Plant - floral visitor","Plant - herbivore"))))

plot.alphas_allyears <- df_param_vert_allyears %>%
  filter(!(focal =="HOMA" & parameter== "Plant - floral visitor")) %>%
  ggplot() +  
  geom_density_ridges_gradient(aes(x=estimate, y=parameter,
                                   fill= after_stat(x)),
                               scale = 1) + 
  facet_grid(as.factor(focal.name)~ complexity.plant) +
  #scale_x_continuous(limits=c(-0.5,0.5)) + 
  #scale_y_discrete(labels=c("Generic Inter-specific",
  #                          "Intra-specific")) + 
  labs(y="Interactions", x= "Estimated distribution",
       color = "Species-specific \n parameters") + 
  scale_color_manual(values=safe_colorblind_palette) +
  scale_fill_gradientn(
    colours=c(re,"white",gr),limits = c(-2, 2),
    #midpoint = 0,
    #space = "Lab",
    na.value = "grey50",
    guide = "none",
    aesthetics = "fill"
  ) +
  theme_minimal() +
  theme(strip.placement = "outside",
        legend.key.size = unit(1, 'cm'),
        title =element_text(size=16),
        axis.text.x= element_text(size=16),
        axis.text.y= element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        strip.text = element_text(size=16),
        panel.border = element_rect(color = "black", 
                                    fill = NA))

plot.alphas_allyears 
ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/Parameters_distribution_species_allyears.pdf"),
       dpi="retina",
       width = 21,
       height = 16,
       units = c("cm"),
       plot.alphas
)

