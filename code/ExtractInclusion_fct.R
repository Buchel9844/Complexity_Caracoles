#---- 3.4. Preliminary fit analysis---- 
#### Determine which parameters warrant inclusion in the final
#       model (i.e. the data pulled their posteriors away from 0). The final model
#       will then be run with only these species-specific parameters, but without
#       the regularized horseshoe priors.
extract.inclusion <- function(summary.interactions.n, DataVec){
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

IntLevel <- 0.2 #0.5 usually, 0.75 for Waitzia, shade

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

Inclusion_all <- list(interaction = summary.interactions.n,
                      Inclusion_ij=Inclusion_ij,
                      Inclusion_FV=Inclusion_FV,
                      Inclusion_H=Inclusion_H,
                      beta_Inclusion_plant=beta_Inclusion_plant,
                      beta_Inclusion_FV=beta_Inclusion_FV,
                      beta_Inclusion_H=beta_Inclusion_H)
Inclusion_all <- append(Inclusion_all,DataVec)

save(Inclusion_all,
     file=paste0(home.dic,"results/Inclusion",FocalPrefix,"_",year,"_",complexity.plant,"_",complexity.animal,".RData"))

summary.interactions <- bind_rows(summary.interactions,summary.interactions.n)

write.csv(summary.interactions,paste0(home.dic,"results/summary.interactions.csv"))
assign("summary.interactions",summary.interactions,1)

DataVec.final <- list(Inclusion_ij = Inclusion_ij,
                      Inclusion_FV=Inclusion_FV,
                      Inclusion_H=Inclusion_H, 
                      beta_Inclusion_plant=beta_Inclusion_plant,
                      beta_Inclusion_FV=beta_Inclusion_FV,
                      beta_Inclusion_H=beta_Inclusion_H,
                      run_estimation = 1)

assign("DataVec.final",DataVec.final,1)
}
