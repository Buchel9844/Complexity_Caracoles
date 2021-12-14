# This script will run the empirical model fits for each focal species and each
#       environmental covariate. A separate script will then make the empirical
#       figures for the manuscript

rm(list = ls())

#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
library(rstan)
#install.packages(HDInterval)
library("HDInterval")
library("tidyverse")
library(dplyr)
options(mc.cores = parallel::detectCores())
#rstan_options(auto_write = TRUE)
summary.interactions <- data.frame()
# Run UploadData.R 2. and 3. or uncomment the following lines 
plant.class <- read.csv("/Users/lisabuche/Code/Project/HOI_Caracoles/data/plant_code.csv", sep=",")
years <- "2020"
focal <- "CHFU"
complexity.plant  <- "class"
complexity.animal <- "group"
source("/Users/lisabuche/Code/Project/HOI_Caracoles/code/Group_FvH.R")
#focal <- "CHFU"
#complexity  <- "family"
for (years in c("2017","2018","2019","2020","2021")){
for (focal in c("MESU","LEMA","HOMA","CHFU")){
for(complexity.plant  in c("code.plant","family","class")){ #add "code.plant" ,make sure it is the name of a column of plant.class
for(complexity.animal  in c("species","family","group")){ 

#--- setting ups----    
  # view complexity levels 
  summary.interactions.n <- data.frame()
  print(focal)
  print(complexity.plant)
  # determine the levels of complexity for plant 
  complexitylevel.plant <- levels(as.factor(plant.class[,complexity.plant]))
  
  print(complexity.animal)
  if (complexity.plant != "code.plant"){complexity.minimize <- T}else{complexity.minimize <- F}
  
FocalPrefix <-focal # "LEMA" or "HOMA" or "CHFU"

# Load in the data and subset out the current focal species and complexity levels
#SpData <- read.csv("/home/lisavm/Simulations/data/competition.csv")
SpData<- read.csv("/Users/lisabuche/Code/Project/HOI_Caracoles/data/competition.csv")
SpData$seed <- round(SpData$seed)

SpData <- na.omit(SpData) 
SpData <- subset(SpData, year == years)
str(SpData)
SpData<- select(SpData,all_of(c("day","month", "year","plot","subplot" ,"focal",
                               "fruit","seed",complexitylevel.plant,FocalPrefix)))

SpDataFocal <- subset(SpData, focal == FocalPrefix )
#head(SpData)


if (complexity.minimize == T){
  levels.of.focal <- plant.class[which(plant.class$code.plant==FocalPrefix),complexity.plant]
  SpDataFocal[,levels.of.focal] <- SpDataFocal[,levels.of.focal] - SpDataFocal[,FocalPrefix]
}

head(SpDataFocal)
# add H presences and floral visitor visits
SpData_H <- subset(get(paste0("herbivore","_",complexity.animal)),year %in% years & 
                            code.plant == focal)
SpData_H <- left_join(subset(SpDataFocal,select=c("day","month","year","plot","subplot")),
                             SpData_H)

SpData_FV <- subset(get(paste0("floral_visitor","_",complexity.animal)),year %in% years & 
                            code.plant == focal)
SpData_FV <- left_join(subset(SpDataFocal,select=c("day","month","year","plot","subplot")),
                            SpData_FV)

view(SpDataFocal)
# determine the levels of complexity  animals
complexitylevel.H <- names(SpData_H)[which(!names(SpData_H) %in% c(complexitylevel.plant,"focal","seed",
                                                                                         "fruit","day.x",FocalPrefix,
                                                                                         "day","month","year","plot",
                                                                                         "subplot","plant",
                                                                                         "code.plant"))]
complexitylevel.FV <- names(SpData_FV)[which(!names(SpData_FV) %in% c(complexitylevel.plant,"seed",
                                                                                              "fruit","day.x","focal",FocalPrefix,
                                                                                              "day","month","year","plot",
                                                                                              "subplot","plant",
                                                                                              "code.plant"))]



# Next continue to extract the data needed to run the model. 
N <- as.integer(nrow(SpDataFocal))
Fecundity <- as.integer(SpDataFocal$seed)  
plot <- as.integer(factor(as.factor(SpDataFocal$plot), levels = c("1","2","3","4","5","6","7","8","9")))

#---- Interaction (direct) matrix of plant with COMP ----

# Now calculate the total number of plant species to use for the model, discounting
#       any species columns with 0 abundance. Save a vector of the species names
#       corresponding to each column for easy matching later.
AllSpNames <- names(SpDataFocal)[!names(SpDataFocal) %in% c("day","month","year","plot",
                                                            "subplot", "focal","fruit","seed")]

AllSpAbunds <- SpDataFocal %>% 
  select(all_of(AllSpNames))

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
SpMatrix <-(SpMatrix/max(SpMatrix))*100 #scale all the interaction between 0 and 100
if(max(SpMatrix) == 100){print("scale SpMatrix_plant correct")}

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

#---- Interaction (direct) matrix of herbivores with FOCAL ----
 
AllSpAbunds_H <- SpData_H %>% 
  select(all_of(complexitylevel.H))

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
SpMatrix_H <-(SpMatrix_H/max(SpMatrix_H,na.rm = T))*100 #scale all the interaction between 0 and 100
if(max(SpMatrix_H,na.rm = T) == 100){print("scale SpMatrix_H correct")}
SpMatrix_H[is.na(SpMatrix_H)] <- 0
SpNames_H <- names(SpToKeep_H)[SpToKeep_H]

#---- Interaction (direct) matrix of floral visitors with FOCAL ----

AllSpAbunds_FV<- SpData_FV %>% 
  select(all_of(complexitylevel.FV))

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
SpMatrix_FV <-(SpMatrix_FV/max(SpMatrix_FV,na.rm = T))*100 #scale all the interaction between 0 and 100
if(max(SpMatrix_FV,na.rm = T) == 100){print("scale SpMatrix_FV correct")}
SpMatrix_FV[is.na(SpMatrix_FV)] <- 0
SpNames_FV <- names(SpToKeep_FV)[SpToKeep_FV]

#---- Interaction (HOIs) matrix of floral visitors with COMPETITORS ----

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

SpData_FV_comp <- left_join(subset(SpDataFocal,select=c("day","month","year","plot","subplot")),
                                  SpData_FV_comp)

SpData_FV_comp <- unique(SpData_FV_comp)

AllSpAbunds_FV_comp <- SpData_FV_comp  %>% 
  select(all_of(names(SpData_FV_comp)[!names(SpData_FV_comp) %in% c('subplot', "plot", "year","month",
                                                                                "day","focal","fruit","seed",
                                                                                focal, complexitylevel.plant)]))
all.interaction.FV_comp <- as.vector(t(outer(c(complexitylevel.plant,focal), complexitylevel.FV, paste, sep="_")))
missing.interaction.FV_comp  <-all.interaction.FV_comp [which(!all.interaction.FV_comp  %in% names(AllSpAbunds_FV_comp))]


for (n in missing.interaction.FV_comp){
  AllSpAbunds_FV_comp[,n] <- NA
}

AllSpAbunds_FV_comp <- select(AllSpAbunds_FV_comp,
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
SpMatrix_FV_comp <-(SpMatrix_FV_comp/max(SpMatrix_FV_comp,na.rm = T))*100 #scale all the interaction between 0 and 100
if(max(SpMatrix_FV_comp,na.rm = T) == 100){print("scale SpMatrix_FV_comp correct")}

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


#---- Interaction (HOIs) matrix of herbivores with COMPETITORS ----


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

SpData_H_comp <- left_join(subset(SpDataFocal,select=c("day","month","year","plot","subplot")),
                                  SpData_H_comp)

SpData_H_comp <- unique(SpData_H_comp)

AllSpAbunds_H_comp <- SpData_H_comp  %>% 
  select(all_of(names(SpData_H_comp)[!names(SpData_H_comp) %in% c('subplot', "plot", "year","month",
                                                                                "day","focal","fruit","seed",
                                                                                focal, complexitylevel.plant)]))
all.interaction.H_comp <- as.vector(t(outer(c(complexitylevel.plant,focal), complexitylevel.H, paste, sep="_")))
missing.interaction.H_comp <-all.interaction.H_comp[which(!all.interaction.H_comp %in% names(AllSpAbunds_H_comp))]


for (n in missing.interaction.H_comp){
  AllSpAbunds_H_comp[,n] <- NA
}

AllSpAbunds_H_comp <- select(AllSpAbunds_H_comp,
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
SpMatrix_H_comp <-(SpMatrix_H_comp/max(SpMatrix_H_comp,na.rm = T))*100 #scale all the interaction between 0 and 100
if(max(SpMatrix_H_comp,na.rm = T) == 100){print("scale SpMatrix_H_comp correct")}

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



#---- Interaction (HOIs) matrix of Herbivores with Herbivores ----
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
  select(all_of(names(SpData_2H)[!names(SpData_2H) %in% c('subplot', "plot", "year","month",
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

AllSpAbunds_2H <- select(AllSpAbunds_2H,
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
SpMatrix_2H <-(SpMatrix_2H/max(SpMatrix_2H,na.rm = T))*100 #scale all the interaction between 0 and 100
if(max(SpMatrix_2H,na.rm = T) == 100){print("scale SpMatrix_2H correct")}

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
  matrix_HOIs_ihh[[i]]  <- matrix_i
  if(sum(rowSums(is.na(matrix_i)))==T){
    print("One or multiple row(s) in matrix_HOIs_ijf contain only NAs")
  }
}

#---- Interaction (HOIs) matrix of floral visitors with floral visitors ----
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
  select(all_of(names(SpData_2FV)[!names(SpData_2FV) %in% c('subplot', "plot", "year","month",
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

AllSpAbunds_2FV <- select(AllSpAbunds_2FV,
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
SpMatrix_2FV <-(SpMatrix_2FV/max(SpMatrix_2FV,na.rm = T))*100 #scale all the interaction between 0 and 100
if(max(SpMatrix_2FV,na.rm = T) == 100){print("scale SpMatrix_2FV correct")}

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
  
  matrix_HOIs_iff[[i]]  <- matrix_i
  if(sum(rowSums(is.na(matrix_i)))==T){
    print("One or multiple row(s) in matrix_HOIs_ijf contain only NAs")
  }
}


#---- Interaction (HOIs) matrix of floral visitors with herbivores----
common.sp.FvH <- complexitylevel.FV[complexitylevel.FV %in% complexitylevel.H]

SpData_FvH <- left_join(SpData_FV,
                                SpData_H,
                                by=c("day","month","year","plot","subplot","code.plant"),
                                suffix=c(".fv",".h")) %>%
  unique() %>%
  gather( all_of(c(complexitylevel.FV[!complexitylevel.FV %in% common.sp.FvH],
                   paste0(common.sp.FvH,".fv"))),
          key = "floral.vis", value = "abundance.fv") %>%
  gather(all_of(c(complexitylevel.H[!complexitylevel.H %in% common.sp.FvH],
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
  select(all_of(names(SpData_FvH)[!names(SpData_FvH) %in% c('subplot', "plot", "year","month",
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

AllSpAbunds_FvH <- select(AllSpAbunds_FvH,
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
SpMatrix_FvH <-(SpMatrix_FvH/max(SpMatrix_FvH,na.rm = T))*100 #scale all the interaction between 0 and 100
if(max(SpMatrix_FvH,na.rm = T) == 100){print("scale SpMatrix_FvH correct")}

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

  matrix_HOIs_ifh[[i]] <-   matrix_i
  if(sum(rowSums(is.na(matrix_i)))==T){
    print("One or multiple row(s) in matrix_HOIs_ijf contain only NAs")
  }
}


#---- Preliminary fit ---- 
# create summary data frame to see all the potential interactions 
summary.interactions.n <- tibble(focal= focal,year=years,
                                                              complexity_plant=complexity.plant,
                                                              n.competitors_plant=length(SpNames)-1,
                                                              competitors_plant=paste(collapse=",",SpNames[which(!SpNames %in% focal)]),
                                                              n.competitors_H=length(SpNames_H),
                                                              competitors_H=paste(collapse=",",SpNames_H),
                                                              n.competitors_FVitor=length(SpNames_FV),
                                                              competitors_FVitor=paste(collapse=",",SpNames_FV),
                                                              n.HOIs_H=length(SpNames_H_comp),
                                                              HOIs_H=paste(collapse=",",SpNames_H_comp),
                                                              n.HOIs_FVitor=length(SpNames_FV_comp),
                                                              HOIs_FVitor=paste(SpNames_FV_comp,collapse=","),
                                                              n.HOIs_2FV=length(SpNames_2FV),
                                                              HOIs_2FV=paste(collapse=",",SpNames_2FV),
                                                              n.HOIs_2h=length(SpNames_2H),
                                                              HOIs_2H=paste(collapse=",",SpNames_2H),
                                                              n.HOIs_FvH=length(SpNames_FvH),
                                                              HOIs_FvH=paste(collapse=",",SpNames_FvH)
                                                              )



# Set the parameters defining the regularized horseshoe prior, as described in
#       the "Incorporating sparsity-inducing priors" section of the manuscript.
tau0 <- 1
slab_scale <- sqrt(2)
slab_df <- 4

# dimension of poll and herb direct interaction 
FV <- ncol(SpMatrix_FV)
H <- ncol(SpMatrix_H)

DataVec <- c("N", "S", "H","FV",
             "Fecundity", "SpMatrix",
             "SpMatrix_H","SpMatrix_FV",
             "matrix_HOIs_plant","matrix_HOIs_ijh",
             "matrix_HOIs_ijf", "matrix_HOIs_ihh",
             "matrix_HOIs_ifh","matrix_HOIs_iff",
             "Intra", "tau0", "slab_scale", "slab_df")

# Now run a perliminary fit of the model to assess parameter shrinkage
print("preliminary fit beginning")

PrelimFit <- stan(file = "/Users/lisabuche/Code/Project/HOI_Caracoles/code/Caracoles_BH_FH_Preliminary.stan", 
                  data = DataVec,
                  iter = 100, 
                  chains = 3)
PrelimPosteriors <- rstan::extract(PrelimFit)
print("prelimi nary fit done")
#---- Preliminary fit analysis---- 
pdf(paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/figure/beta_Inclusion_",years,"_",complexity,"_",FocalPrefix,".pdf"))
##### Diagnostic plots
# First check the distribution of Rhats and effective sample sizes 
# N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
hist(summary(PrelimFit)$summary[,"Rhat"])
hist(summary(PrelimFit)$summary[,"n_eff"])
# Next check the correlation among key model parameters and identify any
pairs(PrelimFit, pars = c("lambdas", "alpha_generic", "alpha_intra","beta_generic"))
# Finally, check for autocorrelation in the posteriors of key model parameters
# N.B.  ACF =  autocorrelation function = 
# coefficient of correlation between two values in a time series
# how to interpret ACF : https://medium.com/analytics-vidhya/interpreting-acf-or-auto-correlation-plot-d12e9051cd14
for (i in 1:nyear){
  print(acf(PrelimPosteriors$lambdas[,i,1]))
  print(acf(PrelimPosteriors$lambdas[,i,2]))
}
print(acf(PrelimPosteriors$alpha_generic[,1]))
print(acf(PrelimPosteriors$alpha_intra[,1]))
print(acf(PrelimPosteriors$beta_generic[,1]))
dev.off()
# Briefly, the scale is from -1 to 1 because it is the correlation coefficient.
# From the graph we can see the lags do not have significant effect
# (within the bounds - cannot tell them from being zero). 
# The ACF function says if the current value depends consistently on previous 
# values (the lags)
str(PrelimPosteriors)
#### If the diagnostic plots don't reveal any problems wiht the model fit, now
#       move on to determining which parameters warrant inclusion in the final
#       model (i.e. the data pulled their posteriors away from 0). The final model
#       will then be run with only these species-specific parameters, but without
#       the regularized horseshoe priors.
Inclusion_ij <- matrix(data = 0, nrow = 1, ncol = length(SpNames),
                       dimnames=list(c(""),c(names(SpTotals[SpTotals!=SpToKeep]))))
Inclusion_FV<- matrix(data = 0,nrow = 1, ncol = length(SpNames_FV),
                             dimnames=list(c(""),c(names(SpTotals_FV[SpTotals_FV!=SpToKeep_FV]))))
Inclusion_H<- matrix(data = 0,nrow = 1, 
                     ncol = length(c(names(SpTotals_H[SpTotals_H!=SpToKeep_H]))),
                             dimnames=list(c(""),c(names(SpTotals_H[SpTotals_H!=SpToKeep_H]))))


beta_Inclusion_plant <- matrix(data = 0,nrow = length(SpNames), 
                               ncol = length(SpNames),
                               dimnames=list(c(names(SpTotals[SpTotals!=SpToKeep])),
                                             c(names(SpTotals[SpTotals!=SpToKeep]))))

beta_Inclusion_FV <- matrix(data = 0,nrow = length(SpNames), 
                            ncol = length(names(SpTotals_FV[SpTotals_FV!=SpToKeep_FV])),
                                   dimnames=list(c(names(SpTotals[SpTotals!=SpToKeep])),
                                                 c(names(SpTotals_FV[SpTotals_FV!=SpToKeep_FV]))))

beta_Inclusion_H<- matrix(data = 0,nrow = length(SpNames), 
                          ncol = length(SpTotals_H[SpTotals_H!=SpToKeep_H]),
                                  dimnames=list(c(names(SpTotals[SpTotals!=SpToKeep])),
                                                c(names(SpTotals_H[SpTotals_H!=SpToKeep_H]))))

beta_Inclusion_2FV<- matrix(data = 0,nrow = length(names(SpTotals_FV[SpTotals_FV!=SpToKeep_FV])),
                            ncol = length(names(SpTotals_FV[SpTotals_FV!=SpToKeep_FV])),
                              dimnames=list(c(names(SpTotals_FV[SpTotals_FV!=SpToKeep_FV])),
                                            c(names(SpTotals_FV[SpTotals_FV!=SpToKeep_FV]))))

beta_Inclusion_2H<- matrix(data = 0,nrow = length(c(names(SpTotals_H[SpTotals_H!=SpToKeep_H]))), 
                           ncol = length(c(names(SpTotals_H[SpTotals_H!=SpToKeep_H]))),
                                    dimnames=list(c(names(SpTotals_H[SpTotals_H!=SpToKeep_H])),
                                                  c(names(SpTotals_H[SpTotals_H!=SpToKeep_H]))))

beta_Inclusion_FvH<- matrix(data = 0,nrow = length(names(SpTotals_FV[SpTotals_FV!=SpToKeep_FV])),
                            ncol = length(c(names(SpTotals_H[SpTotals_H!=SpToKeep_H]))),
                            dimnames=list(c(names(SpTotals_FV[SpTotals_FV!=SpToKeep_FV])),
                                          c(names(SpTotals_H[SpTotals_H!=SpToKeep_H]))))


          
IntLevel <- 0.5 #0.5 usually, 0.75 for Waitzia, shade
is.list(PrelimPosteriors$beta_hat_ij)
#str(PrelimPosteriors$beta_hat_ijk)
#str(PrelimPosteriors$alpha_hat_ij)
#For plants as second actor
for(s in 1:length(SpNames)){
    # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
    Ints_ij <- HDInterval::hdi(PrelimPosteriors$alpha_hat_ij[,s], credMass = IntLevel)
    if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
      Inclusion_ij[1,s] <- 1
    }
    for(m in 1:S){
      beta_Ints_ijk <- HDInterval::hdi(PrelimPosteriors$beta_plant_hat_ijk[,s,m], credMass = IntLevel)
      
      if(beta_Ints_ijk[1] > 0 | beta_Ints_ijk[2] < 0){
        beta_Inclusion_plant[s,m] <- 1
      }
    }
      for(m in 1:length(SpNames_H)){
        beta_Ints_ijh <- HDInterval::hdi(PrelimPosteriors$beta_H_hat_ijh[,s,m], credMass = IntLevel)
        
        if(beta_Ints_ijh[1] > 0 | beta_Ints_ijh[2] < 0){
          beta_Inclusion_H[s,m] <- 1
        }
      }
      for(m in 1:length(SpNames_FV)){
        beta_Ints_ijf <- HDInterval::hdi(PrelimPosteriors$beta_FV_hat_ijf[,s,m], credMass = IntLevel)
        
        if(beta_Ints_ijf[1] > 0 | beta_Ints_ijf[2] < 0){
          beta_Inclusion_FV[s,m] <- 1
        }
      }
  }

#For H as second actor
for(s in 1:length(SpNames_H)){
  # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
  Ints_ih <- HDInterval::hdi(PrelimPosteriors$gamma_H_hat_ih[,s], credMass = IntLevel)
  if(Ints_ih[1] > 0 | Ints_ih[2] < 0){
    Inclusion_H[1,s] <- 1
  }
  for(m in 1:length(SpNames_H)){
    beta_Ints_ihh <- HDInterval::hdi(PrelimPosteriors$beta_2H_hat_ihh[,s,m], credMass = IntLevel)
    
    if(beta_Ints_ihh[1] > 0 | beta_Ints_ihh[2] < 0){
      beta_Inclusion_2H[s,m] <- 1
    }
  }
}

#For pollinator as second actor
for(s in 1:length(SpNames_FV)){
  # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
  Ints_if <- HDInterval::hdi(PrelimPosteriors$gamma_FV_hat_if[,s], credMass = IntLevel)
  if(Ints_if[1] > 0 | Ints_if[2] < 0){
    Inclusion_FV[1,s] <- 1
  }
  for(m in 1:length(SpNames_FV)){
    beta_Ints_iff <- HDInterval::hdi(PrelimPosteriors$beta_2FV_hat_iff[,s,m], credMass = IntLevel)
    
    if(beta_Ints_iff[1] > 0 | beta_Ints_iff[2] < 0){
      beta_Inclusion_2FV[s,m] <- 1
    }
  }
  for(m in 1:length(SpNames_H)){
    beta_Ints_ifh <- HDInterval::hdi(PrelimPosteriors$beta_FvH_hat_ifh[,s,m], credMass = IntLevel)
    
    if(beta_Ints_ifh[1] > 0 | beta_Ints_ifh[2] < 0){
      beta_Inclusion_FvH[s,m] <- 1
    }
  }
}


summary.interactions.n$n.competitors_plant_inclus <- sum(Inclusion_ij)
if (sum(Inclusion_ij) >0){
  summary.interactions.n$competitors_plant_inclus <- c(colnames(Inclusion_ij))[Inclusion_ij]
}
write.csv(Inclusion_ij,
          paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/results/Inclusion_ij_",years,"_",complexity.plant,
                 "_",complexity.animal,"_",FocalPrefix,".csv"))

summary.interactions.n$n.competitors_FV_inclus <- sum(Inclusion_FV)
if (sum(Inclusion_FV)>0){
  summary.interactions.n$competitors_FV_inclus <- c(colnames(Inclusion_FV))[which(Inclusion_FV >=1)]
}
write.csv(Inclusion_FV,
          paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/results/Inclusion_FV_",years,"_",complexity.plant,
                 "_",complexity.animal,"_",FocalPrefix,".csv"))


summary.interactions.n$n.competitors_H_inclus <- sum(Inclusion_H)
if (sum(Inclusion_H)>0){
  summary.interactions.n$competitors_H_inclus <- c(colnames(Inclusion_H))[which(Inclusion_H >=1)]
}
write.csv(Inclusion_H,
          paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/results/Inclusion_H_",years,"_",complexity.plant,
                 "_",complexity.animal,"_",FocalPrefix,".csv"))


summary.interactions.n <- as_tibble(summary.interactions.n)
for ( x in c("plant","H","FV","2FV","2H","FvH")){
  matrix_x <- get(paste0("beta_Inclusion_",x))
  summary.interactions.n[,paste("n.HOIs",x,"inclus",sep="_")] = sum(matrix_x)
  
  if (sum(matrix_x  > 0)){
    summary.interactions.n[,paste0("HOIs_",x,"_inclus")]  = paste(paste0(rownames(  matrix_x )[which(  matrix_x >=1, arr.in=TRUE)[,1]],sep="_",
                                                            colnames(matrix_x )[which(  matrix_x >=1, arr.in=TRUE)[,2]]),collapse=",")
  } 
    write.csv(matrix_x,
            paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/results/beta_Inclusion_",x,"_",years,"_",complexity.plant,
                   "_",complexity.animal,"_",FocalPrefix,".csv"))
}





  }# for complexity 
}}
  summary.interactions <- bind_rows(summary.interactions,summary.interactions.n)
}

write.csv(summary.interactions,"/Users/lisabuche/Code/Project/HOI_Caracoles/results/summary.interactions.csv")


#---- Final fit ---- 

DataVec <- c("N", "S", "Fecundity", "year", "SpMatrix","Intra",
             "Inclusion_ij", "beta_Inclusion")
#FinalFit <- stan(file = "/home/lisavm/Simulations/code/Caracoles_BH_Final.stan", data = DataVec, iter = 1000, chains = 2)
FinalFit <- stan(file = "/Users/lisabuche/Code/Project/HOI_Caracoles/code/Caracoles_BH_Final.stan", 
                 data = DataVec,
                 iter = 100, 
                 chains = 2)

FinalPosteriors <- extract(FinalFit)

# Diagnostic figures

return(acf(FinalPosteriors$lambdas[,1]))
return(acf(FinalPosteriors$lambdas[,2]))
return(acf(FinalPosteriors$alpha_generic[,1]))
return(acf(FinalPosteriors$alpha_intra[,1]))
return(acf(FinalPosteriors$beta_generic[,1]))

#FileName <- paste("/home/lisavm/Simulations/", FocalPrefix, "_", "_FinalFit.rdata", sep = "")
FileName <- paste("/Users/lisabuche/Code/Project/HOI_Caracoles/code/", years,"_",complexity,"_",FocalPrefix, "_FinalFit.rdata", sep = "")

save(FinalFit, SpNames, N, S, Fecundity, year, SpMatrix, Inclusion_ij,
     beta_Inclusion_ij,
     tau0, slab_scale, slab_df, Intra, file = FileName)
str(FinalFit)

