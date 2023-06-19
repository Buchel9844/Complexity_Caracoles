#---- 3.4. Preliminary fit analysis---- 
#### Determine which parameters warrant inclusion in the final
#       model (i.e. the data pulled their posteriors away from 0). The final model
#       will then be run with only these species-specific parameters, but without
#       the regularized horseshoe priors.
extract.inclusion <- function(summary.interactions.n, DataVec){
  SpNames<- c(unlist(strsplit(summary.interactions.n$competitors_plant,",")),focal)
  if(!any(c(SpNames %in% focal) == DataVec$Intra)) {
    print("wrong place of intra")
  } # check otherwise reorder! 
  
  Inclusion_ij <- matrix(data = 0, nrow = 1, ncol = summary.interactions.n$n.competitors_plant + 1,
                         dimnames=list(c(""),
                                       SpNames))
  
  SpNames_FV <- unlist(strsplit(summary.interactions.n$interactors_FV,","))
  Inclusion_FV<- matrix(data = 0,nrow = 1, ncol = summary.interactions.n$n.interactors_FV,
                        dimnames=list(c(""),
                                      SpNames_FV))
  SpNames_H <- unlist(strsplit(summary.interactions.n$interactors_H,","))
  
  Inclusion_H <- matrix(data = 0,nrow = 1, 
                        ncol = summary.interactions.n$n.interactors_H,
                        dimnames=list(c(""),
                                      SpNames_H))
  beta_Inclusion_plant <- matrix(data = 0,nrow = summary.interactions.n$n.competitors_plant + 1, 
                                 ncol = summary.interactions.n$n.competitors_plant + 1,
                                 dimnames=list(c(unlist(strsplit(summary.interactions.n$competitors_plant,",")),focal),
                                               c(unlist(strsplit(summary.interactions.n$competitors_plant,",")),focal)))
  
  beta_Inclusion_FV <- matrix(data = 0, nrow = summary.interactions.n$n.competitors_plant+1, 
                              ncol = summary.interactions.n$n.interactors_FV,
                              dimnames=list(c(unlist(strsplit(summary.interactions.n$competitors_plant,",")),focal),
                                            unlist(strsplit(summary.interactions.n$interactors_FV,","))))
  beta_Inclusion_H<- matrix(data = 0, nrow = summary.interactions.n$n.competitors_plant +1, 
                            ncol = summary.interactions.n$n.interactors_H,
                            dimnames=list(c(unlist(strsplit(summary.interactions.n$competitors_plant,",")),focal),
                                          unlist(strsplit(summary.interactions.n$interactors_H,","))))
  
  beta_Inclusion_2FV<- matrix(data = 0, nrow = summary.interactions.n$n.interactors_FV, 
                              ncol = summary.interactions.n$n.interactors_FV,
                              dimnames=list(unlist(strsplit(summary.interactions.n$interactors_FV,",")),
                                            unlist(strsplit(summary.interactions.n$interactors_FV,","))))
  beta_Inclusion_2H<-matrix(data = 0, nrow = summary.interactions.n$n.interactors_H, 
                            ncol = summary.interactions.n$n.interactors_H,
                            dimnames=list(unlist(strsplit(summary.interactions.n$interactors_H,",")),
                                          unlist(strsplit(summary.interactions.n$interactors_H,","))))
  
  level.H <- c(complexitylevel.H,"beetle.h","bug.h")
  level.Fv <- c(complexitylevel.FV,"beetle.fv","bug.fv")
  
  beta_Inclusion_FvH<- matrix(data = 0,nrow = DataVec$FvH_Fv,
                              ncol = DataVec$FvH_h,
                              dimnames=list(level.Fv[which(level.Fv %in% unlist(strsplit(unlist(strsplit(summary.interactions.n$HOIs_FvH,",")),"_")))],
                                            level.H[which(level.H %in% unlist(strsplit(unlist(strsplit(summary.interactions.n$HOIs_FvH,",")),"_")))]
                              ))
  
  IntLevel <- 0.2 #0.5 usually, 0.75 for Waitzia, shade
  
  #str(PrelimPosteriors$beta_hat_ijk)
  #str(PrelimPosteriors$alpha_hat_ij)
  #For plants as second actor
  for(s in 1:length(SpNames)){
    # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
    Ints_ij <- HDInterval::hdi(PrelimPosteriors$alpha_hat_ij[,s], credMass = IntLevel)
    if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
      Inclusion_ij[1,s] <- 1
    }
    for(m in 1:DataVec$S){
      beta_Ints_ijk <- HDInterval::hdi(PrelimPosteriors$beta_plant_hat_ijk[,s,m], credMass = IntLevel)
      
      if(beta_Ints_ijk[1] > 0 | beta_Ints_ijk[2] < 0){
        beta_Inclusion_plant[s,m] <- 1
      }
    }
    for(m in 1:length(SpNames_H)){
      beta_Ints_ijh <- HDInterval::hdi(PrelimPosteriors$beta_H_hat_ijh[,s,m], credMass = IntLevel)
      if(is.na(beta_Ints_ijh[1]) | is.na(beta_Ints_ijh[2])) next
      if(beta_Ints_ijh[1] > 0 | beta_Ints_ijh[2] < 0){
        beta_Inclusion_H[s,m] <- 1
      }
    }
    for(m in 1:length(SpNames_FV)){
      beta_Ints_ijf <- HDInterval::hdi(PrelimPosteriors$beta_FV_hat_ijf[,s,m], credMass = IntLevel)
      if(is.na(beta_Ints_ijf[1]) | is.na(beta_Ints_ijf[2])) next
      
      if(beta_Ints_ijf[1] > 0 | beta_Ints_ijf[2] < 0){
        beta_Inclusion_FV[s,m] <- 1
      }
    }
  }
  
  #For H as second actor
  for(s in 1:length(SpNames_H)){
    # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
    Ints_ih <- HDInterval::hdi(PrelimPosteriors$gamma_H_hat_ih[,s], credMass = IntLevel)
    if(is.na(Ints_ih[1]) | is.na(Ints_ih[2])) next
    
    if(Ints_ih[1] > 0 | Ints_ih[2] < 0){
      Inclusion_H[1,s] <- 1
    }
    for(m in 1:length(SpNames_H)){
      beta_Ints_ihh <- HDInterval::hdi(PrelimPosteriors$beta_2H_hat_ihh[,s,m], credMass = IntLevel)
      if(is.na(beta_Ints_ihh[1]) | is.na(beta_Ints_ihh[2])) next
      if(beta_Ints_ihh[1] > 0 | beta_Ints_ihh[2] < 0){
        beta_Inclusion_2H[s,m] <- 1
      }
    }
  }
  
  #For pollinator as second actor
  for(s in 1:length(SpNames_FV)){
    # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
    Ints_if <- HDInterval::hdi(PrelimPosteriors$gamma_FV_hat_if[,s], credMass = IntLevel)
    if(is.na(Ints_if[1]) | is.na(Ints_if[2])) next
    
    if(Ints_if[1] > 0 | Ints_if[2] < 0){
      Inclusion_FV[1,s] <- 1
    }
    for(m in 1:length(SpNames_FV)){
      beta_Ints_iff <- HDInterval::hdi(PrelimPosteriors$beta_2FV_hat_iff[,s,m], credMass = IntLevel)
      if(is.na(beta_Ints_iff[1]) | is.na(beta_Ints_iff[2])) next
      
      if(beta_Ints_iff[1] > 0 | beta_Ints_iff[2] < 0){
        beta_Inclusion_2FV[s,m] <- 1
      }
    }
    
  }
  
  for(s in 1:DataVec$FvH_Fv){
    for(m in 1:DataVec$FvH_h){
      beta_Ints_ifh <- HDInterval::hdi(PrelimPosteriors$beta_FvH_hat_ifh[,s,m], credMass = IntLevel)
      if(is.na(beta_Ints_ifh[1]) | is.na(beta_Ints_ifh[2])) next
      
      if(beta_Ints_ifh[1] > 0 | beta_Ints_ifh[2] < 0){
        beta_Inclusion_FvH[s,m] <- 1
      }
    }
  }
  
  if(dim(beta_Inclusion_2FV)[2]==0){
    beta_Inclusion_FV <- matrix(nrow=nrow(beta_Inclusion_FV),
                                ncol=1,
                                data=rep(0,times=nrow(beta_Inclusion_FV)))
  }
  
  summary.interactions.n$n.competitors_plant_inclus <- sum(Inclusion_ij)
  if (sum(Inclusion_ij) > 0){
    summary.interactions.n$competitors_plant_inclus <- paste(c(colnames(Inclusion_ij))[which(Inclusion_ij>=1)],collapse="_")
  }
  
  
  summary.interactions.n$n.competitors_FV_inclus <- sum(Inclusion_FV)
  if (sum(Inclusion_FV)>0){
    summary.interactions.n$competitors_FV_inclus <- paste(c(colnames(Inclusion_FV))[which(Inclusion_FV >=1)],collapse="_")
  }
  
  
  summary.interactions.n$n.competitors_H_inclus <- sum(Inclusion_H)
  if (sum(Inclusion_H)>0){
    summary.interactions.n$competitors_H_inclus <- paste(c(colnames(Inclusion_H))[which(Inclusion_H >=1)],collapse="_")
  }
  
  summary.interactions.n <- as_tibble(summary.interactions.n)
  for ( x in c("plant","H","FV","2FV","2H","FvH")){
    matrix_x <- get(paste0("beta_Inclusion_",x))
    summary.interactions.n[,paste("n.HOIs",x,"inclus",sep="_")] = sum(matrix_x)
    
    if (sum(matrix_x  > 0)){
      summary.interactions.n[,paste0("HOIs_",x,"_inclus")]  = paste(paste0(rownames(  matrix_x )[which(  matrix_x >=1, arr.in=TRUE)[,1]],sep="_",
                                                                           colnames(matrix_x )[which(  matrix_x >=1, arr.in=TRUE)[,2]]),collapse=",")
    } 
  }
  if(dim(beta_Inclusion_FvH)[1]==0){
    beta_Inclusion_FvH<- matrix(data = 0,nrow = 1,
                                ncol = length(SpNames_H),
                                dimnames=list("nopoll",
                                              SpNames_H))
    
  }
  if(dim(beta_Inclusion_2FV)[1]==0){
    beta_Inclusion_2FV<- matrix(data = 0,nrow = 1,
                                ncol = 1,
                                dimnames=list("nopoll",
                                              "nopoll"))
    
  }
  if(dim(beta_Inclusion_FV)[2]==0){
    beta_Inclusion_FV<- matrix(data = 0,nrow = length(SpNames), 
                               ncol = 1,
                               dimnames=list(c(names(SpTotals[SpTotals!=SpToKeep])),
                                             "nopoll"))
    
  }
  Inclusion_all <- list(interaction = summary.interactions.n,
                        Inclusion_ij=Inclusion_ij,
                        Inclusion_FV=Inclusion_FV,
                        Inclusion_H=Inclusion_H,
                        beta_Inclusion_plant=beta_Inclusion_plant,beta_Inclusion_FV=beta_Inclusion_FV,beta_Inclusion_H=beta_Inclusion_H,
                        beta_Inclusion_2FV=beta_Inclusion_2FV,beta_Inclusion_2H=beta_Inclusion_2H,beta_Inclusion_FvH=beta_Inclusion_FvH)
  Inclusion_all <- append(Inclusion_all,DataVec)
  
  save(Inclusion_all,
       file=paste0(home.dic,"results/Inclusion",FocalPrefix,"_",years,"_",complexity.plant,"_",complexity.animal,".RData"))
  
  summary.interactions <- bind_rows(summary.interactions,summary.interactions.n)
  
  write.csv(summary.interactions,paste0(home.dic,"results/summary.interactions.csv"))
  assign("summary.interactions",summary.interactions,1)
  
  DataVec.final <- list(Inclusion_ij = Inclusion_ij,
                        Inclusion_FV=Inclusion_FV,
                        Inclusion_H=Inclusion_H, 
                        beta_Inclusion_plant=beta_Inclusion_plant,
                        beta_Inclusion_FV=beta_Inclusion_FV,
                        beta_Inclusion_H=beta_Inclusion_H,
                        beta_Inclusion_2FV=beta_Inclusion_2FV,
                        beta_Inclusion_2H = beta_Inclusion_2H,
                        beta_Inclusion_FvH =beta_Inclusion_FvH
  )
  
  assign("DataVec.final",DataVec.final,1)
}
