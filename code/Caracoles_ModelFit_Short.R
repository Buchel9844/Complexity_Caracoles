# This script will run the empirical model fits for each focal species and each
#       environmental covariate. A separate script will then make the empirical
#       figures for the manuscript
# https://mc-stan.org/docs/2_29/functions-reference/index.html#overview

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 1. SET UP: Import data, create df with abundances ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
#---- 1.1. Import the Data ----
setwd("~/Eco_Bayesian/Complexity_caracoles")
abundance.plant <- read.csv("data/abundance.csv", sep=",")

competition.plant <- read.csv("data/competition.csv", sep=",")

plant.class <- read.csv("data/plant_code.csv", sep=",")

herbivorypredator <- read.csv("data/herbivorypredator.csv", sep=",")

floral_visitor <- read.csv("data/floral_visitor.csv", sep=",")

# Run UploadData.R 2. and 3. or uncomment the following lines 
plant.class <- read.csv("data/plant_code.csv", sep=",")


#---- 1.2. Import packages ----
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
library(rstan)
#install.packages("HDInterval")
library("HDInterval")
#install.packages("tidyverse")
library(tidyverse)
#install.packages("dplyr")
library(dplyr)
options(mc.cores = parallel::detectCores()) # to use the core at disposition 

rstan_options(auto_write = TRUE) 

#---- 1.3. Set parameters ----
summary.interactions <- data.frame()
#years <- "2020"
focal <- "LEMA"
FocalPrefix<- "LEMA"
complexity.plant  <- "class"
complexity.animal <- "group"
#"2017",'2018',
for (years in c("2019","2020","2021")){
  #for ( focal in c("CHFU","LEMA",'HOMA',"CETE")){
  # for ( complexity.plant in c("code.plant","family","class","rareORabundant")){   
  #for ( complexity.animal in c("species","group","family")){ }
  # }
  # }
  
  
  #focal <- "CHFU"
  #complexity  <- "family"
  #---- 1.4. Import df of abundances----    
  # view complexity levels 
  summary.interactions.n <- data.frame()
  print(paste(focal,complexity.plant,complexity.animal))
  
  # determine the levels of complexity for plant 
  complexitylevel.plant <- levels(as.factor(plant.class[,complexity.plant]))
  
  
  if (complexity.plant != "code.plant"){complexity.minimize <- T}else{complexity.minimize <- F}
  
  FocalPrefix <- focal # "LEMA" or "HOMA" or "CHFU"
  
  # Load in the data and subset out the current focal species and complexity levels
  #SpData <- read.csv("/home/lisavm/Simulations/data/competition.csv")
  SpData <- read.csv("data/SpData.csv")
  
  SpData$seed <- round(SpData$seed)
  
  SpData <- na.omit(SpData) 
  SpData <- subset(SpData, year == years)
  
  SpData <- dplyr::select(SpData,all_of(c("day","month", "year","plot","subplot" ,"focal",
                                          "fruit","seed",complexitylevel.plant,FocalPrefix)))
  
  SpDataFocal <- subset(SpData, focal == FocalPrefix )
  #FocalPrefix <- "CHFU" # "LEMA" or "HOMA" or "CHFU"
  
  
  ggplot(SpDataFocal, aes(x=seed)) + geom_density() 
  
  #head(SpDataFocal)
  
  
  
  if (complexity.minimize == T){
    levels.of.focal <- plant.class[which(plant.class$code.plant==FocalPrefix),complexity.plant]
    SpDataFocal[,levels.of.focal] <- SpDataFocal[,levels.of.focal] - SpDataFocal[,FocalPrefix]
  }
  
  head(SpDataFocal)
  
  # Standardized function for species abundance ATTENTION MAYBE need to be log for HTL
  standard.abudance <- function(x){
    x <- (x*5)/max(x,na.rm=T) 
    return(x)
  }
  
  SpDataFocal[,c(focal,complexitylevel.plant)] <- standard.abudance(SpDataFocal[,c(focal,complexitylevel.plant)]) 
  
  if(!max(SpDataFocal[,c(focal,complexitylevel.plant)],na.rm=T) ==5){print("Standardisation of plants UNcorrect")}
  
  
  #import higher trophic levels data 
  source("code/Group_FvH.R")
  # add H presences and floral visitor visits
  SpData_H <- subset(get(paste0("herbivore","_",complexity.animal)),year %in% years & 
                       code.plant == focal)
  
  SpData_H <- left_join(subset(SpDataFocal,select=c("day","month","year","plot","subplot")),
                        SpData_H)
  
  SpData_FV <- subset(get(paste0("floral_visitor","_",complexity.animal)),year %in% years & 
                        code.plant == focal)
  SpData_FV <- left_join(subset(SpDataFocal,select=c("day","month","year","plot","subplot")),
                         SpData_FV)
  
  
  # determine the levels of complexity  animals
  complexitylevel.H <- names(SpData_H)[which(!names(SpData_H) %in% c(complexitylevel.plant,"focal","seed",
                                                                     "fruit","day.x",FocalPrefix,
                                                                     "day","month","year","plot",
                                                                     "subplot","plant",
                                                                     "code.plant"))]
  
  SpData_H[,complexitylevel.H] <- standard.abudance(SpData_H[,complexitylevel.H])
  if(!max(SpData_H[,complexitylevel.H],na.rm=T) ==5){print("Standardisation of Herbivory UNcorrect")}
  
  complexitylevel.FV <- names(SpData_FV)[which(!names(SpData_FV) %in% c(complexitylevel.plant,"seed",
                                                                        "fruit","day.x","focal",FocalPrefix,
                                                                        "day","month","year","plot",
                                                                        "subplot","plant",
                                                                        "code.plant"))]
  
  SpData_FV[,complexitylevel.FV] <- standard.abudance(SpData_FV[,complexitylevel.FV])
  if(!max(SpData_FV[,complexitylevel.FV],na.rm=T) ==5){print("Standardisation of Floral visitors UNcorrect")}
  
  
  # Next continue to extract the data needed to run the model. 
  N <- as.integer(nrow(SpDataFocal))
  Fecundity <- SpDataFocal$seed
  plot <- as.integer(factor(as.factor(SpDataFocal$plot), levels = c("1","2","3","4","5","6","7","8","9")))
  
  
  #---- 1.5. Pre-analysis of the data----    
  
  head(SpData_FV)
  
  
  ggplot(gather(SpDataFocal,c(FocalPrefix,complexitylevel.plant),
                key="plant.complexity",value="abundance"),
         aes(y=seed,x=abundance, color=plant.complexity)) + geom_line() + theme_bw()
  
  #ggplot(gather(SpData_FV,all_of(complexitylevel.FV),
  #              key="FV.complexity",value="abundance"),
  #       aes(y=Fecundity,x=abundance, color=FV.complexity)) + geom_line() + theme_bw()
  
  #ggplot(gather(SpData_H,all_of(complexitylevel.H),
  #             key="H.complexity",value="abundance"),
  #     aes(y=seed, x=abundance, color=H.complexity)) + geom_line() + theme_bw()
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #---- 2. ABUNDANCE MATRIX ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #---- 2.1. Interaction (direct) matrix of plant with COMP ----
  
  # Now calculate the total number of plant species to use for the model, discounting
  #       any species columns with 0 abundance. Save a vector of the species names
  #       corresponding to each column for easy matching later.
  AllSpNames <- names(SpDataFocal)[!names(SpDataFocal) %in% c("day","month","year","plot",
                                                              "subplot", "focal","fruit","seed")]
  
  AllSpAbunds <- SpDataFocal %>% 
    dplyr::select(all_of(AllSpNames))
  
  SpTotals <- colSums(AllSpAbunds)
  SpToKeep <- SpTotals > 0
  
  S <- sum(SpToKeep)
  SpMatrix <- matrix(NA, nrow = N, ncol = S)
  i <- 1
  for(s in 1:ncol(AllSpAbunds)){
    if(SpToKeep[s] == 1){
      SpMatrix[,i] <- AllSpAbunds[,s]
      i <- i + 1
    }else{next}
  }
  #SpMatrix <-round((SpMatrix/max(SpMatrix))*100) #scale all the interaction between 0 and 100
  #if(max(SpMatrix) == 100){print("scale SpMatrix_plant correct")}
  
  SpNames <- AllSpNames[SpToKeep]
  
  #assign(paste0("SpNames_",FocalPrefix),
  #     SpNames)
  Intra <- ifelse(SpNames == FocalPrefix, 1, 0)
  
  # creation of a matrix of S by S of the interaction jk in HOIs_ijk for plants
  matrix_HOIs_plant <- list()
  for (i in 1:N){
    matrix_i <- matrix(nrow=S,ncol=S)
    for (n in 1:S) {
      for (m in 1:S) {
        if (m <= n){
          matrix_i[n,m] = 0
        }
        else{
          matrix_i[n,m] = SpMatrix[i,n]* SpMatrix[i,m]
        } 
      }
    }
    matrix_i[is.na(matrix_i)] <- 0
    matrix_HOIs_plant[[i]] <-  matrix_i
  }
  
  #---- 2.2. Interaction (direct) matrix of herbivores with FOCAL ----
  
  AllSpAbunds_H <- SpData_H %>% 
    dplyr::select(all_of(complexitylevel.H))
  
  SpTotals_H <- colSums(AllSpAbunds_H, na.rm = T)
  SpToKeep_H  <- SpTotals_H  > 0
  H <- sum(SpToKeep_H)
  
  SpMatrix_H  <- matrix(NA, nrow = N, ncol = H)
  i <- 1
  for(h in 1:length(SpToKeep_H)){
    if(SpToKeep_H[h] == 1){
      SpMatrix_H[,i] <- AllSpAbunds_H[,h]
      i <- i + 1
    }else{next}
  }
  #SpMatrix_H <-(SpMatrix_H/max(SpMatrix_H,na.rm = T))*100 #scale all the interaction between 0 and 100
  #if(max(SpMatrix_H,na.rm = T) == 100){print("scale SpMatrix_H correct")}
  SpMatrix_H[is.na(SpMatrix_H)] <- 0
  SpNames_H <- names(SpToKeep_H)[SpToKeep_H]
  
  #---- 2.3. Interaction (direct) matrix of floral visitors with FOCAL ----
  
  AllSpAbunds_FV<- SpData_FV %>% 
    dplyr::select(all_of(complexitylevel.FV))
  
  SpTotals_FV <- colSums(AllSpAbunds_FV, na.rm = T)
  SpToKeep_FV <- SpTotals_FV  > 0
  
  FV <- sum(SpToKeep_FV)
  SpMatrix_FV  <- matrix(NA, nrow = N, ncol = FV)
  i <- 1
  for(s in 1:ncol(AllSpAbunds_FV)){
    if(SpToKeep_FV[s] == 1){
      SpMatrix_FV[,i] <- AllSpAbunds_FV[,s]
      i <- i + 1
    }else{next}
  }
  #SpMatrix_FV <-(SpMatrix_FV/max(SpMatrix_FV,na.rm = T))*100 #scale all the interaction between 0 and 100
  #if(max(SpMatrix_FV,na.rm = T) == 100){print("scale SpMatrix_FV correct")}
  SpMatrix_FV[is.na(SpMatrix_FV)] <- 0
  SpNames_FV <- names(SpToKeep_FV)[SpToKeep_FV]
  
  #---- 2.4. Interaction (HOIs) matrix of floral visitors with COMPETITORS ----
  
  SpData_FV_comp <- left_join(SpDataFocal,
                              SpData_FV,
                              by=c("day","month","year","plot","subplot")) %>%
    unique() %>%
    gather( all_of(c(complexitylevel.plant,focal)),
            key = "plant.comp", value = "abundance.plant.comp") %>%
    gather(all_of(complexitylevel.FV),
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
  
  SpData_FV_comp <- left_join(subset(SpDataFocal,select=c("day","month","year","plot","subplot","seed","fruit")),
                              SpData_FV_comp,by=c("day","month","year","plot","subplot","seed","fruit"),
                              copy=FALSE)
  if( !nrow(SpData_FV_comp) == N){
    print("problem in the numbr of observation, they don't match, 
        try to add  unique(SpData_FV_comp) before the previous line.")
  }
  
  
  AllSpAbunds_FV_comp <- SpData_FV_comp  %>% 
    dplyr::select(all_of(names(SpData_FV_comp)[!names(SpData_FV_comp) %in% c('subplot', "plot", "year","month",
                                                                             "day","focal","fruit","seed",
                                                                             focal, complexitylevel.plant)]))
  all.interaction.FV_comp <- as.vector(t(outer(c(complexitylevel.plant,focal), complexitylevel.FV, paste, sep="_")))
  missing.interaction.FV_comp  <-all.interaction.FV_comp [which(!all.interaction.FV_comp  %in% names(AllSpAbunds_FV_comp))]
  
  
  for (n in missing.interaction.FV_comp){
    AllSpAbunds_FV_comp[,n] <- NA
  }
  
  AllSpAbunds_FV_comp <- dplyr::select(AllSpAbunds_FV_comp,
                                       all_of(all.interaction.FV_comp)) # make the matrix in specific order
  
  
  
  #view(AllSpAbunds_FV_comp)
  # this part is strickly to know the potential of interaction at the before the sparsity approach
  
  SpTotals_FV_comp <- colSums(AllSpAbunds_FV_comp, na.rm = T)
  SpToKeep_FV_comp <- SpTotals_FV_comp  > 0
  
  FV_comp <- sum(SpToKeep_FV_comp)
  SpMatrix_FV_comp  <- matrix(NA, nrow = N, ncol = FV_comp )
  
  i <- 1
  for(s in 1:ncol(AllSpAbunds_FV_comp)){
    if(SpToKeep_FV_comp[s] == 1){
      SpMatrix_FV_comp[,i] <- AllSpAbunds_FV_comp[,s]
      i <- i + 1
    }else{next}
  }
  #SpMatrix_FV_comp <-(SpMatrix_FV_comp/max(SpMatrix_FV_comp,na.rm = T))*100 #scale all the interaction between 0 and 100
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
    for (fv in 1:FV_comp) {
      if(fv == vector_FV_comp[n] + 1){
        n = n+1
        m = 1 }
      matrix_i[n,m] = SpMatrix_FV_comp[i,fv]
      m = m+1 
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
                             by=c("day","month","year","plot","subplot")) %>%
    unique() %>%
    gather( all_of(c(complexitylevel.plant,focal)),
            key = "plant.comp", value = "abundance.plant.comp") %>%
    gather(all_of(complexitylevel.H),
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
  
  SpData_H_comp <- left_join(subset(SpDataFocal,select=c("day","month","year","plot","subplot","seed","fruit")),
                             SpData_H_comp,by=c("day","month","year","plot","subplot","seed","fruit"),
                             SpData_H_comp)
  if( !nrow(SpData_H_comp) == N){
    print("problem in the numbr of observation, they don't match, 
        try to add  unique(SpData_H_comp) before the previous line.")
  }
  
  AllSpAbunds_H_comp <- SpData_H_comp  %>% 
    dplyr::select(all_of(names(SpData_H_comp)[!names(SpData_H_comp) %in% c('subplot', "plot", "year","month",
                                                                           "day","focal","fruit","seed",
                                                                           focal, complexitylevel.plant)]))
  all.interaction.H_comp <- as.vector(t(outer(c(complexitylevel.plant,focal), complexitylevel.H, paste, sep="_")))
  missing.interaction.H_comp <-all.interaction.H_comp[which(!all.interaction.H_comp %in% names(AllSpAbunds_H_comp))]
  
  
  for (n in missing.interaction.H_comp){
    AllSpAbunds_H_comp[,n] <- NA
  }
  
  AllSpAbunds_H_comp <- dplyr::select(AllSpAbunds_H_comp,
                                      all_of(all.interaction.H_comp)) # make the matrix in specific order
  
  
  
  #view(AllSpAbunds_H_comp)
  # this part is strickly to know the potential of interaction at the before the sparsity approach
  
  SpTotals_H_comp <- colSums(AllSpAbunds_H_comp, na.rm = T)
  SpToKeep_H_comp <- SpTotals_H_comp  > 0
  
  FV_comp <- sum(SpToKeep_H_comp)
  SpMatrix_H_comp  <- matrix(NA, nrow = N, ncol = FV_comp)
  i <- 1
  for(s in 1:ncol(AllSpAbunds_H_comp)){
    if(SpToKeep_H_comp[s] == 1){
      SpMatrix_H_comp[,i] <- AllSpAbunds_H_comp[,s]
      i <- i + 1
    }else{next}
  }
  #SpMatrix_H_comp <-(SpMatrix_H_comp/max(SpMatrix_H_comp,na.rm = T))*5 #scale all the interaction between 0 and 5
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
    for (h_comp in 1:H_comp) {
      if(h_comp == vector_H_comp[n] + 1){
        n = n+1
        m = 1 }
      matrix_i[n,m] = SpMatrix_H_comp[i,h_comp]
      m = m+1 
    }
    matrix_i[is.na(matrix_i)] <- 0
    matrix_HOIs_ijh[[i]] <-  matrix_i
    if(sum(rowSums(is.na(matrix_i)))==T){
      print("One or multiple row(s) in matrix_HOIs_ijf contain only NAs")
    }
  }
  
  
  
  #---- 2.6. Interaction (HOIs) matrix of herbivores with herbivores ----
  SpData_2H <- left_join(SpData_H,
                         SpData_H,
                         by=c("day","month","year","plot","subplot","code.plant"),
                         suffix=c(".h1",".h2")) %>%
    unique() %>%
    gather( all_of(paste0(complexitylevel.H,".h1")),
            key = "herbivore.1", value = "presence.h.1") %>%
    gather(all_of(paste0(complexitylevel.H,".h2")),
           key = "herbivore.2", value = "presence.h.2")
  
  SpData_2H$presence.h <- SpData_2H$presence.h.1*SpData_2H$presence.h.2
  SpData_2H <-  unite(SpData_2H,
                      herbivore.1, herbivore.2,
                      col = "herbivore_herbivore", sep = "_") 
  
  SpData_2H <- aggregate(presence.h ~ herbivore_herbivore + code.plant + subplot 
                         + plot + year + month + day, 
                         SpData_2H, sum,na.action=na.omit)
  
  SpData_2H <-  spread(SpData_2H ,herbivore_herbivore, presence.h  )
  
  SpData_2H <- left_join(subset(SpDataFocal,select=c("day","month","year","plot","subplot")),
                         SpData_2H)
  
  
  AllSpAbunds_2H <- SpData_2H  %>% 
    dplyr::select(all_of(names(SpData_2H)[!names(SpData_2H) %in% c('subplot', "plot", "year","month",
                                                                   "day","focal","fruit","seed",
                                                                   focal, complexitylevel.plant)]))
  # check all posible interactions of h with h and see if they are included in the dataframe
  all.interaction_2H <- as.vector(t(outer(paste0(complexitylevel.H,".h1"),
                                          paste0(complexitylevel.H,".h2"),
                                          paste, sep="_")))
  
  missing.interaction_2H <-all.interaction_2H[which(!all.interaction_2H %in% 
                                                      names(AllSpAbunds_2H))]
  
  for (n in missing.interaction_2H){ # if not present, include them
    AllSpAbunds_2H[,n] <- NA
  }
  
  AllSpAbunds_2H <- dplyr::select(AllSpAbunds_2H,
                                  all_of(all.interaction_2H)) # make the matrix in specific order
  
  
  
  #view(AllSpAbunds_2H)
  # this part is strickly to know the potential of interaction at the before the sparsity approach
  
  SpTotals_2H <- colSums(AllSpAbunds_2H, na.rm = T)
  SpToKeep_2H <- SpTotals_2H  > 0
  
  H_H  <- sum(SpToKeep_2H)
  SpMatrix_2H  <- matrix(NA, nrow = N, ncol = H_H )
  i <- 1
  for(s in 1:ncol(AllSpAbunds_2H)){
    if(SpToKeep_2H[s] == 1){
      SpMatrix_2H[,i] <- AllSpAbunds_2H[,s]
      i <- i + 1
    }else{next}
  }
  #SpMatrix_2H <-(SpMatrix_2H/max(SpMatrix_2H,na.rm = T))*100 #scale all the interaction between 0 and 100
  #if(max(SpMatrix_2H,na.rm = T) == 100){print("scale SpMatrix_2H correct")}
  
  SpNames_2H <- names(SpToKeep_2H)[SpToKeep_2H]
  
  vector_2H <- vector() # make a vector of the position of each competitor plant species
  # important to be able to divide in the bayesian approach  
  for (n in 1:(length(complexitylevel.H))){
    hname <- paste0("^",paste0(complexitylevel.H,".h1")[n],"_")
    if(length(grep(hname, SpNames_2H)) ==0){vector_2H[n] <- NA
    next }
    vector_2H[n] <- max(grep(hname, SpNames_2H))
    if(vector_2H[n] < 0){vector_2H[n] <- NA} 
  }
  
  vector_2H <- vector_2H[!is.na(vector_2H)]
  
  # creation of matrix_HOIs_ihh
  H_H <- ncol(SpMatrix_2H)
  matrix_HOIs_ihh <- list()
  for (i in 1:N){
    matrix_i <- matrix(nrow=H,ncol=H)  
    n = 1
    m=1
    for (h_h in 1:H_H) {
      if(h_h == vector_2H[n] + 1){
        n = n+1
        m = 1 }
      matrix_i[n,m] = SpMatrix_2H[i,h_h]
      m = m+1 
    }
    matrix_i[is.na(matrix_i)] <- 0
    matrix_i <- round( matrix_i)
    matrix_HOIs_ihh[[i]]  <- matrix_i
    if(sum(rowSums(is.na(matrix_i)))==T){
      print("One or multiple row(s) in matrix_HOIs_ijf contain only NAs")
    }
  }
  
  #---- 2.7. Interaction (HOIs) matrix of floral visitors with floral visitors ----
  SpData_2FV <- left_join(SpData_FV,
                          SpData_FV,
                          by=c("day","month","year","plot","subplot","code.plant"),
                          suffix=c(".fv1",".fv2")) %>%
    unique() %>%
    gather( all_of(paste0(complexitylevel.FV,".fv1")),
            key = "floral.vis.1", value = "abundance.fv.1") %>%
    gather(all_of(paste0(complexitylevel.FV,".fv2")),
           key = "floral.vis.2", value = "abundance.fv.2")
  
  SpData_2FV$abundance.fv <- SpData_2FV$abundance.fv.1*SpData_2FV$abundance.fv.2
  SpData_2FV <-  unite(SpData_2FV,
                       floral.vis.1, floral.vis.2,
                       col = "Fv_Fv", sep = "_") 
  
  SpData_2FV <- aggregate(abundance.fv ~ Fv_Fv + code.plant + subplot 
                          + plot + year + month + day, 
                          SpData_2FV, sum,na.action=na.omit)
  
  SpData_2FV <-  spread(SpData_2FV ,Fv_Fv, abundance.fv  )
  
  SpData_2FV <- left_join(subset(SpDataFocal,select=c("day","month","year","plot","subplot")),
                          SpData_2FV)
  
  
  AllSpAbunds_2FV <- SpData_2FV  %>% 
    dplyr::select(all_of(names(SpData_2FV)[!names(SpData_2FV) %in% c('subplot', "plot", "year","month",
                                                                     "day","focal","fruit","seed",
                                                                     focal, complexitylevel.plant)]))
  # check all posible interactions of Fv with FV and see if they are included in the dataframe
  all.interaction.2FV <- as.vector(t(outer(paste0(complexitylevel.FV,".fv1"),
                                           paste0(complexitylevel.FV,".fv2"),
                                           paste, sep="_")))
  
  missing.interaction.2FV <-all.interaction.2FV[which(!all.interaction.2FV %in% 
                                                        names(AllSpAbunds_2FV))]
  
  for (n in missing.interaction.2FV){ # if not present, include them
    AllSpAbunds_2FV[,n] <- NA
  }
  
  AllSpAbunds_2FV <- dplyr::select(AllSpAbunds_2FV,
                                   all_of(all.interaction.2FV)) # make the matrix in specific order
  
  
  
  #view(AllSpAbunds_2FV)
  # this part is strickly to know the potential of interaction at the before the sparsity approach
  
  SpTotals_2FV <- colSums(AllSpAbunds_2FV, na.rm = T)
  SpToKeep_2FV <- SpTotals_2FV  > 0
  
  FV_FV  <- sum(SpToKeep_2FV)
  SpMatrix_2FV  <- matrix(NA, nrow = N, ncol = FV_FV )
  i <- 1
  for(s in 1:ncol(AllSpAbunds_2FV)){
    if(SpToKeep_2FV[s] == 1){
      SpMatrix_2FV[,i] <- AllSpAbunds_2FV[,s]
      i <- i + 1
    }else{next}
  }
  #SpMatrix_2FV <-(SpMatrix_2FV/max(SpMatrix_2FV,na.rm = T))*100 #scale all the interaction between 0 and 100
  #if(max(SpMatrix_2FV,na.rm = T) == 100){print("scale SpMatrix_2FV correct")}
  
  SpNames_2FV <- names(SpToKeep_2FV)[SpToKeep_2FV]
  
  vector_2FV <- vector() # make a vector of the position of each competitor plant species
  # important to be able to divide in the bayesian approach  
  for (n in 1:(length(complexitylevel.FV))){
    fvname <- paste0("^",paste0(complexitylevel.FV,".fv1")[n],"_")
    if(length(grep(fvname, SpNames_2FV)) ==0){vector_2FV[n] <- NA
    next }
    vector_2FV[n] <- max(grep(fvname, SpNames_2FV))
    if(vector_2FV[n] < 0){vector_2FV[n] <- NA} 
  }
  vector_2FV <- vector_2FV[!is.na(vector_2FV)]
  
  # creation of matrix_HOIs_iff
  
  FV_FV <- ncol(SpMatrix_2FV)
  matrix_HOIs_iff <- list()
  for (i in 1:N){
    matrix_i <- matrix(nrow=FV,ncol=FV)  
    n = 1
    m=1
    for (f_f in 1:FV_FV) {
      if(f_f == vector_2FV[n] + 1){
        n = n+1
        m = 1 }
      matrix_i[n,m] = SpMatrix_2FV[i,f_f]
      m = m+1 
    }
    matrix_i[is.na(matrix_i)] <- 0
    matrix_i <- round( matrix_i)
    matrix_HOIs_iff[[i]]  <- matrix_i
    if(sum(rowSums(is.na(matrix_i)))==T){
      print("One or multiple row(s) in matrix_HOIs_ijf contain only NAs")
    }
  }
  
  
  #---- 2.8. Interaction (HOIs) matrix of floral visitors with herbivores----
  common.sp.FvH <- complexitylevel.FV[complexitylevel.FV %in% complexitylevel.H]
  
  SpData_FvH <- left_join(SpData_FV,
                          SpData_H,
                          by=c("day","month","year","plot","subplot","code.plant"),
                          suffix=c(".fv",".h")) %>%
    unique() %>%
    tidyr::gather( all_of(c(complexitylevel.FV[!complexitylevel.FV %in% common.sp.FvH],
                            paste0(common.sp.FvH,".fv"))),
                   key = "floral.vis", value = "abundance.fv") %>%
    tidyr::gather(all_of(c(complexitylevel.H[!complexitylevel.H %in% common.sp.FvH],
                           paste0(common.sp.FvH,".h"))),
                  key = "herbivore", value = "presence.h")
  # to avoid multiplying apple and pear scale everything on a scale from 0 to 100
  # we will scale the multiplication again later on
  SpData_FvH$presence.h <- SpData_FvH$presence.h/max(SpData_FvH$presence.h,na.rm = T)*100
  SpData_FvH$abundance.fv <- SpData_FvH$abundance.fv/max(SpData_FvH$abundance.fv,na.rm = T)*100
  
  
  SpData_FvH$presence.h.abundance.fv <- SpData_FvH$presence.h*SpData_FvH$abundance.fv
  SpData_FvH <-  unite(SpData_FvH,
                       floral.vis, herbivore,
                       col = "FvH", sep = "_") 
  
  SpData_FvH <- aggregate(presence.h.abundance.fv ~FvH + code.plant + subplot 
                          + plot + year + month + day, 
                          SpData_FvH, sum,na.action=na.omit)
  
  SpData_FvH <-  spread(SpData_FvH ,FvH, presence.h.abundance.fv )
  
  SpData_FvH <- left_join(subset(SpDataFocal,select=c("day","month","year","plot","subplot")),
                          SpData_FvH)
  
  
  AllSpAbunds_FvH <- SpData_FvH  %>% 
    dplyr::select(all_of(names(SpData_FvH)[!names(SpData_FvH) %in% c('subplot', "plot", "year","month",
                                                                     "day","focal","fruit","seed",
                                                                     focal, complexitylevel.plant)]))
  all.interaction.FvH <- as.vector(t(outer(c(complexitylevel.FV[!complexitylevel.FV %in% common.sp.FvH],
                                             "beetle.fv","bug.fv"),
                                           c(complexitylevel.H[!complexitylevel.H %in% common.sp.FvH],
                                             "beetle.h","bug.h"),paste, sep="_")))
  missing.interaction.FvH <-all.interaction.FvH[which(!all.interaction.FvH %in% 
                                                        names(AllSpAbunds_FvH))]
  
  
  for (n in missing.interaction.FvH){
    AllSpAbunds_FvH[,n] <- NA
  }
  
  AllSpAbunds_FvH <- dplyr::select(AllSpAbunds_FvH,
                                   all_of(all.interaction.FvH)) # make the matrix in specific order
  
  
  
  #view(AllSpAbunds_FvH)
  # this part is strickly to know the potential of interaction at the before the sparsity approach
  
  SpTotals_FvH <- colSums(AllSpAbunds_FvH, na.rm = T)
  SpToKeep_FvH <- SpTotals_FvH  > 0
  
  FvH  <- sum(SpToKeep_FvH)
  SpMatrix_FvH  <- matrix(NA, nrow = N, ncol = FvH )
  i <- 1
  for(s in 1:ncol(AllSpAbunds_FvH)){
    if(SpToKeep_FvH[s] == 1){
      SpMatrix_FvH[,i] <- AllSpAbunds_FvH[,s]
      i <- i + 1
    }else{next}
  }
  #SpMatrix_FvH <-(SpMatrix_FvH/max(SpMatrix_FvH,na.rm = T))*100 #scale all the interaction between 0 and 100
  #if(max(SpMatrix_FvH,na.rm = T) == 100){print("scale SpMatrix_FvH correct")}
  
  SpNames_FvH <- names(SpToKeep_FvH)[SpToKeep_FvH]
  
  vector_FvH <- vector() # make a vector of the position of each competitor plant species
  # important to be able to divide in the bayesian approach  
  for (n in 1:(length(complexitylevel.FV))){
    fvname <- paste0("^",c(complexitylevel.FV[!complexitylevel.FV %in% common.sp.FvH],
                           "beetle.fv","bug.fv")[n],"_")
    if(length(grep(fvname, SpNames_FvH)) ==0){vector_FvH[n] <- NA
    next }
    vector_FvH[n] <- max(grep(fvname, SpNames_FvH))
    if(vector_FvH[n] < 0){vector_FvH[n] <- NA} 
  }
  vector_FvH <- vector_FvH[!is.na(vector_FvH)]
  
  # creation of matrix_HOIs_ifh
  FV_H <- ncol(SpMatrix_FvH)
  matrix_HOIs_ifh <- list()
  for (i in 1:N){
    matrix_i <- matrix(nrow=FV,ncol=H)  
    n = 1;
    m=1
    for (Fvh in 1:FV_H) {
      if(Fvh == vector_FvH[n] + 1){
        n = n+1;
        m = 1 }
      matrix_i[n,m] = SpMatrix_FvH[i,Fvh];
      m = m+1 ;
    }
    matrix_i[is.na(matrix_i)] <- 0
    matrix_i <- round( matrix_i)
    matrix_HOIs_ifh[[i]] <-   matrix_i
    if(sum(rowSums(is.na(matrix_i)))==T){
      print("One or multiple row(s) in matrix_HOIs_ijf contain only NAs")
    }
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #---- 3'. SHORT PRELIMINARY FIT----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #---- 3'.1. Set up summary interactions df and parameters ---- 
  # Set the parameters defining the regularized horseshoe prior, as described in
  #       the "Incorporating sparsity-inducing priors" section of the manuscript.
  tau0 <- 1
  slab_scale <- sqrt(2)
  slab_df <- 4
  
  # dimension of poll and herb direct interaction 
  FV <- ncol(SpMatrix_FV)
  H <- ncol(SpMatrix_H)
  run_estimation <- 1
  
  
  DataVec <- c("N", "S", "H","FV",
               "Fecundity", "SpMatrix",
               "SpMatrix_H","SpMatrix_FV",
               "Intra", "tau0", "slab_scale", "slab_df")
  
  #---- 3'.2. Run  preliminary fit ----
  
  # Now run a preliminary fit of the model to assess parameter shrinkage
  print("preliminary fit beginning")
  options(mc.cores=parallel::detectCores())
  #install.packages("codetools")
  library("codetools")
  
  
  Short_PrelimFit <- stan(file = "code/Short_Caracoles_BH_FH_Preliminary.stan", 
                          data = DataVec,
                          #init = 0 , # all initial values are 0 
                          control=list(max_treedepth=12),
                          iter = 1000, 
                          #options(mc.cores=parallel::detectCores()), 
                          chains = 2)
  
  
  save(file=paste0("results/Short_PrelimFit",focal,"_",years,"_",complexity.plant,"_",complexity.animal,".rds"),
       Short_PrelimFit)
  
  load(paste0("results/Short_PrelimFit",focal,"_",years,"_",complexity.plant,"_",complexity.animal,".rds"))
  Short_PrelimPosteriors <- rstan::extract(Short_PrelimFit)
  print("preliminary fit done")
  
  #---- 3'.3. Preliminary fit posterior check and behavior checks---- 
  ##### Diagnostic plots
  pdf(paste0("figure/Short_PrelimFit_",paste(focal,years,complexity.plant,complexity.animal,sep="_"),".pdf"))
  title= paste0("Diagnostic plots for Short_PrelimFit of ",focal, " in ", years," \nwhen plant complexity is ",complexity.plant,
                " \nand HTL complexity is ",complexity.animal)
  # Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
  source("Test_simulation/code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
  # check the distribution of Rhats and effective sample sizes 
  stan_post_pred_check(Short_PrelimPosteriors,"F_hat",Fecundity,1000) 
  
  # N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
  hist(summary(Short_PrelimFit)$summary[,"Rhat"], 
       main = paste("Short_PrelimFit fit: Histogram of Rhat for",
                    focal, " in ", years," \nwhen plant complexity is ",complexity.plant,
                    " \nand HTL complexity is ",complexity.animal))
  hist(summary(Short_PrelimFit)$summary[,"n_eff"],
       main = paste("Short_PrelimFit fit: Histogram of n_eff for",
                    focal, " in ", years," \nwhen plant complexity is ",complexity.plant,
                    " \nand HTL complexity is ",complexity.animal))
  
  # Next check the correlation among key model parameters and identify any
  #       divergent transitions
  pairs(Short_PrelimFit, pars = c("lambdas", "alpha_intra","alpha_generic_tilde",
                                  "gamma_H","gamma_FV"))
  
  
  # plot the traceplot and distribution graphs
  stan_model_check(Short_PrelimFit,
                   param =c('lambdas','alpha_generic_tilde',"alpha_hat_ij_tilde",
                            'disp_dev'))
  
  
  dev.off()
  
  #---- 3'.4. Preliminary fit analysis---- 
  #### Determine which parameters warrant inclusion in the final
  #       model (i.e. the data pulled their posteriors away from 0). The final model
  #       will then be run with only these species-specific parameters, but without
  #       the regularized horseshoe priors.
  Inclusion_ij <- matrix(data = 0, nrow = 1, ncol = length(SpNames),
                         dimnames=list(c(""),c(names(SpTotals[SpTotals!=SpToKeep]))))
  Inclusion_FV<- matrix(data = 0,nrow = 1, ncol = length(SpNames_FV),
                        dimnames=list(c(""),c(SpNames_FV)))
  Inclusion_H <- matrix(data = 0,nrow = 1, 
                        ncol = length(SpNames_H),
                        dimnames=list(c(""),c(SpNames_H)))
  
  IntLevel <- 0.5#0.5 usually, 0.75 for Waitzia, shade
  
  #str(PrelimPosteriors$beta_hat_ijk)
  #str(PrelimPosteriors$alpha_hat_ij)
  #For plants as second actor
  for(s in 1:length(SpNames)){
    # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
    Ints_ij <- HDInterval::hdi(Short_PrelimPosteriors$alpha_hat_ij[,s], credMass = IntLevel)
    if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
      Inclusion_ij[1,s] <- 1
    }
  }
  
  #For H as second actor
  for(s in 1:length(SpNames_H)){
    # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
    Ints_ih <- HDInterval::hdi(Short_PrelimPosteriors$gamma_H_hat_ih[,s], credMass = IntLevel)
    if(Ints_ih[1] > 0 | Ints_ih[2] < 0){
      Inclusion_H[1,s] <- 1
    }
  }
  
  for(s in 1:length(SpNames_FV)){
    # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
    Ints_if <- HDInterval::hdi(PrelimPosteriors$gamma_FV_hat_if[,s], credMass = IntLevel)
    if(Ints_if[1] > 0 | Ints_if[2] < 0){
      Inclusion_FV[1,s] <- 1
    }
  }
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #---- 4'. SHORT FINAL FIT----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #---- 4'.1. Set up parameters ---- 
  run_estimation = 1
  DataVec <- c("N", "S", "H","FV","Intra","run_estimation",
               "Fecundity", "SpMatrix",
               "SpMatrix_H","SpMatrix_FV",
               "Inclusion_ij","Inclusion_FV","Inclusion_H"
  )
  #FinalFit <- stan(file =  "/home/lisavm/Simulations/code/Caracoles_BH_Final.stan", data = DataVec, iter = 1000, chains = 2)
  
  
  #---- 4.2. Run final fit ---- 
  print("final fit begins")
  
  Short_FinalFit <- stan(file = "code/Short_Caracoles_BH_Final.stan", 
                         data = DataVec,
                         #init=0, # all initial values are 0 
                         control=list(max_treedepth=12),
                         iter = 1000, 
                         chains = 2)
  
  save(file=paste0("results/Short_FinalFit",focal,"_",years,"_",complexity.plant,"_",complexity.animal,".rds"),
       Short_FinalFit)
  load(paste0("results/Short_FinalFit",focal,"_",years,"_",complexity.plant,"_",complexity.animal,".rds")) 
  
  FinalPosteriors <- rstan::extract(Short_FinalFit)
  print("final fit is done")
  #---- 4.3. Final fit posterior check and behavior checks---- 
  
  ##### Diagnostic plots and post prediction 
  pdf(paste0("figure/Short_FinalFit_",paste(focal,years,complexity.plant,complexity.animal,sep="_"),".pdf"))
  # Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
  source("Test_simulation/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
  # check the distribution of Rhats and effective sample sizes 
  ##### Posterior check
  
  stan_post_pred_check(FinalPosteriors,"F_hat",Fecundity) 
  
  # N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
  hist(summary(Short_FinalFit)$summary[,"Rhat"],
       main = paste("Finat Fit: Histogram of Rhat for",
                    focal, " in ", years," \nwhen plant complexity is ",complexity.plant,
                    " \nand HTL complexity is ",complexity.animal))
  hist(summary(Short_FinalFit)$summary[,"n_eff"],
       main = paste("Finat Fit: Histogram of Rhat for",
                    focal, " in ", years," \nwhen plant complexity is ",complexity.plant,
                    " \nand HTL complexity is ",complexity.animal))
  
  # plot the corresponding graphs
  stan_model_check(Short_FinalFit,
                   param =c('lambdas','alpha_generic_tilde','gamma_H_generic_tilde','gamma_FV_generic_tilde',
                            'disp_dev'))
  
  # Next check the correlation among key model parameters and identify any
  pairs(Short_FinalFit, pars = c("lambdas", "alpha_generic", "alpha_intra",
                                 "gamma_FV","gamma_H"))
  
  
  dev.off()
  
  #---- 4.4. Final fit analysis---
  
  #FileName <- paste("/home/lisavm/Simulations/", FocalPrefix, "_", "_Short_FinalFit.rdata", sep = "")
  FileName <- paste("results/", years,"_",complexity.plant,"_",FocalPrefix, "_Short_FinalFit.rdata", sep = "")
  
  save(Short_FinalFit, SpNames, N,S,FV,H, Fecundity,years,complexity.plant,
       FocalPrefix, 
       SpMatrix,SpMatrix_H,SpMatrix_FV, 
       matrix_HOIs_plant,matrix_HOIs_ijh, matrix_HOIs_ijf, 
       matrix_HOIs_ihh, matrix_HOIs_ifh,matrix_HOIs_iff,
       Intra,Inclusion_ij,Inclusion_FV,Inclusion_H,beta_Inclusion_plant,
       beta_Inclusion_H,beta_Inclusion_FV,
       beta_Inclusion_2H,beta_Inclusion_2FV,beta_Inclusion_FvH,
       tau0, slab_scale, slab_df,file = FileName)
  str(Short_FinalFit)
  
  
  
}
