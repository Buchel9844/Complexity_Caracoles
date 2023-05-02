# extraction of matrixes

extract.matrix <- function(focal= c("LEMA","HOMA","CHFU","CETE"),
                           year=c('2018',"2019","2020","2021"),
                           complexity.plant =c("class","family","code.plant"),
                           complexity.animal= c("group","family","species"),
                           plant.class,SpData ){
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
  SpData <- subset(SpData, year == years)
  
  SpData <- dplyr::select(SpData,all_of(c("day","month", "year","plot","subplot" ,"focal",
                                          "fruit","seed",complexitylevel.plant,FocalPrefix)))
  
  SpDataFocal <- subset(SpData, focal == FocalPrefix )
  
  #FocalPrefix <- "CHFU" # "LEMA" or "HOMA" or "CHFU"
  
  
  #ggplot(SpDataFocal, aes(x=seed)) + geom_density() 
  
  #head(SpDataFocal)
  
  if (complexity.minimize == T){
    levels.of.focal <- plant.class[which(plant.class$code.plant==FocalPrefix),complexity.plant]
    SpDataFocal[,levels.of.focal] <- SpDataFocal[,levels.of.focal] - SpDataFocal[,FocalPrefix]}else{
      SpDataFocal <- rename(SpDataFocal, c(focal_neighbour = focal[1]))
      complexitylevel.plant <- complexitylevel.plant[which(!complexitylevel.plant %in% focal )]
      focal <- "focal_neighbour"}
  
  #head(SpDataFocal)
  
  # Standardized function for species abundance ATTENTION MAYBE need to be log for HTL
  standard.abudance <- function(x){
    x <- (x*5)/max(x,na.rm=T) 
    return(x)}
  
  SpDataFocal[,c(focal[1],complexitylevel.plant)] <- standard.abudance(SpDataFocal[,c(focal,complexitylevel.plant)]) 
  
  if(!max(SpDataFocal[,c(focal,complexitylevel.plant)],na.rm=T) ==5){print("Standardisation of plants UNcorrect")}
  
  
  #import higher trophic levels data 
  source("code/Group_FvH.R")
  
  # add H presences and floral visitor visits
  SpData_H <- subset(get(paste0("herbivore","_",complexity.animal)),year %in% years & 
                       code.plant == FocalPrefix)
  
  SpData_H <- left_join(subset(SpDataFocal,select=c("day","month","year","plot","subplot")),
                        multiple = "all",
                        SpData_H)
  SpData_H[is.na(SpData_H)] <- 0
  
  SpData_FV <- subset(get(paste0("floral_visitor","_",complexity.animal)),year %in% years & 
                        code.plant == FocalPrefix)
  SpData_FV <- left_join(subset(SpDataFocal,select=c("day","month","year","plot","subplot")),
                         SpData_FV)
  SpData_FV[is.na(SpData_FV)] <- 0
  
  # determine the levels of complexity  animals
  complexitylevel.H <- names(SpData_H)[which(!names(SpData_H) %in% c(complexitylevel.plant,"focal","seed",
                                                                     "fruit","day.x",FocalPrefix,
                                                                     "day","month","year","plot",
                                                                     "subplot","plant",
                                                                     "code.plant"))]
  
  SpData_H[,complexitylevel.H] <- standard.abudance(SpData_H[,complexitylevel.H])
  SpData_H[is.na(SpData_H)] <- 0
  if(!max(SpData_H[,complexitylevel.H],na.rm=T) ==5){print("Standardisation of Herbivory UNcorrect")}
  if(sum(colSums(SpData_H[,complexitylevel.H]),na.rm=T)==0){print("No herbivory")}
  
  complexitylevel.FV <- names(SpData_FV)[which(!names(SpData_FV) %in% c(complexitylevel.plant,"seed",
                                                                        "fruit","day.x","focal",FocalPrefix,
                                                                        "day","month","year","plot",
                                                                        "subplot","plant",
                                                                        "code.plant"))]
  
  SpData_FV[,complexitylevel.FV] <- standard.abudance(SpData_FV[,complexitylevel.FV])
  SpData_FV[is.na(SpData_FV)] <- 0
  if(!max(SpData_FV[,complexitylevel.FV],na.rm=T) ==5){print("Standardisation of Floral visitors UNcorrect")}
  if(sum(colSums(SpData_FV[,complexitylevel.FV]),na.rm=T)==0){print("No floral visits")}
  
  # Next continue to extract the data needed to run the model. 
  N <- as.integer(nrow(SpDataFocal))
  Fecundity <- SpDataFocal$seed
  plot <- as.integer(factor(as.factor(SpDataFocal$plot), levels = c("1","2","3","4","5","6","7","8","9")))
  
  
  #---- 1.5. Pre-analysis of the data----    
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
      i <- i + 1}else{next}
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
        }else{
          matrix_i[n,m] = SpMatrix[i,n]* SpMatrix[i,m]
        } 
      }
    }
    matrix_i[is.na(matrix_i)] <- 0
    matrix_HOIs_plant[[i]] <-  matrix_i
  }
  
  SpNames_plant_HOIs <- c(apply(combn(SpNames,2),2,paste,collapse="_"),
                          paste0(SpNames,sep="_",SpNames))
  
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
      i <- i + 1}else{next}
  }
  #SpMatrix_H <-(SpMatrix_H/max(SpMatrix_H,na.rm = T))*100 #scale all the interaction between 0 and 100
  #if(max(SpMatrix_H,na.rm = T) == 100){print("scale SpMatrix_H correct")}
  SpMatrix_H[is.na(SpMatrix_H)] <- 0
  SpNames_H <- names(SpToKeep_H)[SpToKeep_H]
  
  if(dim(SpMatrix_H)[2]==0){
    SpMatrix_H  <- matrix(NA, nrow = N, ncol = 1)
    SpMatrix_H[,1] <- rep(0,times=nrow(SpMatrix_H))
  }
  
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
      i <- i + 1}else{next}
  }
  
  if(dim(SpMatrix_FV)[2]==0){
    SpMatrix_FV  <- matrix(NA, nrow = N, ncol = 1)
    SpMatrix_FV[,1] <- rep(0,times=nrow(SpMatrix_FV))
  }
  
  #SpMatrix_FV <-(SpMatrix_FV/max(SpMatrix_FV,na.rm = T))*100 #scale all the interaction between 0 and 100
  #if(max(SpMatrix_FV,na.rm = T) == 100){print("scale SpMatrix_FV correct")}
  SpMatrix_FV[is.na(SpMatrix_FV)] <- 0
  SpNames_FV <- names(SpToKeep_FV)[SpToKeep_FV]
  if(length(SpNames_FV)==0){SpNames_FV <- "nopoll"}
  
  #---- 2.4. Interaction (HOIs) matrix of floral visitors with COMPETITORS ----
  
  SpData_FV_comp <- left_join(SpDataFocal,
                              SpData_FV,
                              multiple = "all",
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
                              multiple = "all",
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
      i <- i + 1}else{next}
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
                             multiple = "all",
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
                             multiple = "all",
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
  
  H_comp <- sum(SpToKeep_H_comp)
  SpMatrix_H_comp  <- matrix(NA, nrow = N, ncol = H_comp)
  i <- 1
  for(s in 1:ncol(AllSpAbunds_H_comp)){
    if(SpToKeep_H_comp[s] == 1){
      SpMatrix_H_comp[,i] <- AllSpAbunds_H_comp[,s]
      i <- i + 1}else{next}
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
  
  
  
  #---- 2.6. Interaction (HOIs) matrix of herbivores with herbivores ----
  SpData_2H <- left_join(SpData_H,
                         SpData_H,
                         multiple = "all",
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
                         multiple = "all",
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
    if(H_H ==0){
      matrix_i <- matrix(nrow=1,ncol=1, data=rep(0,times=S)) 
    }else{
      for (h_h in 1:H_H) {
        if(h_h == vector_2H[n] + 1){
          n = n+1
          m = 1 }
        matrix_i[n,m] = SpMatrix_2H[i,h_h]
        m = m+1 
      }
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
                          multiple = "all",
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
                          multiple = "all",
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
    if(FV_FV ==0){
      matrix_i <- matrix(nrow=1,ncol=1, data=rep(0,times=1)) 
    }else{
      for (f_f in 1:FV_FV) {
        if(!is.na(vector_2FV[n]) & f_f == vector_2FV[n] + 1){
          n = n+1
          m = 1 }
        matrix_i[n,m] = SpMatrix_2FV[i,f_f]
        m = m+1 
      }
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
                          multiple = "all",
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
  if(max(SpData_FvH$abundance.fv,na.rm = T)<=0 | max(SpData_FvH$presence.h,na.rm = T) <= 0){
    SpData_FvH$presence.h <- rep(0,times=length(SpData_FvH$presence.h))
    SpData_FvH$abundance.fv <-rep(0,times=length(SpData_FvH$abundance.fv))
  }else{
    SpData_FvH$presence.h <- SpData_FvH$presence.h/max(SpData_FvH$presence.h,na.rm = T)*100
    SpData_FvH$abundance.fv <- SpData_FvH$abundance.fv/max(SpData_FvH$abundance.fv,na.rm = T)*100
  }
  
  SpData_FvH$presence.h.abundance.fv <- SpData_FvH$presence.h*SpData_FvH$abundance.fv
  SpData_FvH <-  unite(SpData_FvH,
                       floral.vis, herbivore,
                       col = "FvH", sep = "_") 
  
  SpData_FvH <- aggregate(presence.h.abundance.fv ~ FvH + code.plant  + subplot +
                            + plot + year + month + day, 
                          SpData_FvH, sum,na.action=na.omit)
  
  SpData_FvH <-  spread(SpData_FvH ,FvH, presence.h.abundance.fv )
  
  SpData_FvH <- left_join(subset(SpDataFocal,select=c("day","month","year","plot","subplot")),
                          multiple = "all",
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
  FvH_Fv <- length(vector_FvH) # extract the number of floral visitor present in FV - H interactions
  FvH_h <- length(which(c(complexitylevel.H,"beetle.h","bug.h") %in% unlist(strsplit(SpNames_FvH, "_")))) # same for herbivores
  
  
  
  if(length(names(SpTotals_FV[SpToKeep_FV]))==0){FV = 1}else{FV <- length(names(SpTotals_FV[SpToKeep_FV]))}
  if(length(names(SpTotals_H[SpToKeep_H]))==0){H = 1}else{H <- length(names(SpTotals_H[SpToKeep_H]))}
  
  # creation of matrix_HOIs_ifh
  FV_H <- ncol(SpMatrix_FvH)
  matrix_HOIs_ifh <- list()
  for (i in 1:N){
    matrix_i <- matrix(nrow=FvH_Fv,ncol=FvH_h)  
    n = 1;
    m=1
    if( FV_H==0){
      matrix_i <- matrix(nrow=FV,ncol=H,data=rep(0,times=H)) 
    }else{
      if(FV==0){
        matrix_i <- matrix(nrow=1,ncol=H,data=rep(0,times=H)) 
      }else{
        for (Fvh in 1:FV_H) {
          if(Fvh == vector_FvH[n] + 1){
            n = n+1;
            m = 1 }
          matrix_i[n,m] = SpMatrix_FvH[i,Fvh];
          m = m+1 ;
        }
      }
    }
    matrix_i[is.na(matrix_i)] <- 0
    matrix_i <- round( matrix_i)
    matrix_HOIs_ifh[[i]] <-   matrix_i
    if(sum(rowSums(is.na(matrix_i)))==T){
      print("One or multiple row(s) in matrix_HOIs_ijf contain only NAs")
    }
  }
  
  
  #---- 3. Return object ----
  
  summary.interactions.n <- tibble(focal= focal,year=as.integer(years),
                                   complexity_plant=complexity.plant,
                                   n.competitors_plant=length(SpNames)-1, # remove focal in count
                                   competitors_plant=paste(collapse=",",
                                                           SpNames[which(!SpNames %in% focal)]),
                                   n.interactors_H=length(SpNames_H),
                                   interactors_H=paste(collapse=",",SpNames_H),
                                   n.interactors_FV=length(SpNames_FV),
                                   interactors_FV=paste(collapse=",",SpNames_FV),
                                   n.HOIs_plant=length(SpNames_plant_HOIs),
                                   HOIs_plant=paste(collapse=",",SpNames_plant_HOIs),
                                   n.HOIs_H=length(SpNames_H_comp),
                                   HOIs_H=paste(collapse=",",SpNames_H_comp),
                                   n.HOIs_FV=length(SpNames_FV_comp),
                                   HOIs_FV=paste(SpNames_FV_comp,collapse=","),
                                   n.HOIs_2FV=length(SpNames_2FV),
                                   HOIs_2FV=paste(collapse=",",SpNames_2FV),
                                   n.HOIs_2h=length(SpNames_2H),
                                   HOIs_2H=paste(collapse=",",SpNames_2H),
                                   n.HOIs_FvH=length(SpNames_FvH),
                                   HOIs_FvH=paste(collapse=",",SpNames_FvH))
  
  assign('summary.interactions.n',  summary.interactions.n,1)
  
  # dimension of poll and herb direct interaction 
  if(length(names(SpTotals_FV[SpToKeep_FV]))==0){FV = 1 
  FvH_Fv = 1
  }else{FV <- length(names(SpTotals_FV[SpToKeep_FV]))}
  
  if(length(names(SpTotals_H[SpToKeep_H]))==0){H = 1 
  FvH_H = 1
  }else{HExtract <- length(names(SpTotals_H[SpToKeep_H]))}
  
  
  run_estimation <- 1
  # rfemove interaction of Herbivore and floral visitors vector 
  if(sum(colSums(SpMatrix_H))==0|sum(colSums(SpMatrix_FV))==0){
    RemoveFvH <- 1
  }else{RemoveFvH <- 0}
  # Set the parameters defining the regularized horseshoe prior, as described in
  #       the "Incorporating sparsity-inducing priors" section of the manuscript.
  
  tau0 <- 1
  slab_scale <- sqrt(2)
  slab_df <- 4
  
  DataVec <- list(N=N, 
                  SpNames = SpNames,
                  S=S, H=H,FV=FV,FvH_h =FvH_h,
                  FvH_Fv=FvH_Fv,
                  RemoveFvH= RemoveFvH,
                  Fecundity=Fecundity, 
                  SpMatrix=SpMatrix,
                  SpMatrix_H=SpMatrix_H,
                  SpMatrix_FV=SpMatrix_FV,
                  matrix_HOIs_plant=matrix_HOIs_plant,
                  matrix_HOIs_ijh=matrix_HOIs_ijh,
                  matrix_HOIs_ijf=matrix_HOIs_ijf, 
                  matrix_HOIs_ihh = matrix_HOIs_ihh,
                  matrix_HOIs_ifh=matrix_HOIs_ifh,
                  matrix_HOIs_iff=matrix_HOIs_iff,
                  Intra=Intra, tau0=tau0, slab_scale=slab_scale, slab_df=slab_df
  )
  
  assign('DataVec',DataVec,1)
  assign("complexitylevel.H",complexitylevel.H,1)
  assign("complexitylevel.FV",complexitylevel.FV,1)
  
}
