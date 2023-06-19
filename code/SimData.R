# Function to create a simulated dataset 
simul_data <- function(
    # number of species has to be a multiple of 10
  S = 10,   # number of focal groups / species
  K = 10,   # number of neighbour focals 
  pI = 0.1,   # proportion of interactions which are NOT observed 
  HTL= T, # add higher trophic level 
  FvH = 10,   # number of HTL, here floral visitor or not floral visitor
  pI_HTL = 0.3, # proportion of HTL interactions which are NOT observed 
  ...
) {
  
  # Neighbourhood abundance counts
  #-------------------------------
  # we assume we have a different number of observations for each focal
  S_obs <- c(rpois(.1*S, 300), rpois(.2*S, 80), rpois(.5*S, 50), rpois(.2*S, 35))
  
  # we assume that S = K, and that the number of observations of K is relative to their abundance
  approx_K_obs <- round(jitter(10*(S_obs), amount = 30))
  summedSK <- round(jitter(sapply(S_obs, '*', approx_K_obs)/1000, 20))
  
  # randomly select observations that will be set to 0
  obs_to_rm <- sample(seq_along(1:(S*K)), pI*S*K)
  # remove unobserved interactions
  for(i in 1:(pI*S*K)) {
    summedSK[[obs_to_rm[i]]] <- 0
  }
  
  # set up an empty matrix with columns for each neighbours
  K_Nmat <- matrix(ncol = K)
  for(s in 1:S){ 
    # create a matrix of observations for each focal species
    foo <- matrix(data = 0, nrow = S_obs[s], ncol = K)
    for(j in 1:K) {
      if(summedSK[s, j] > 0)
        for (i in 1:summedSK[s, j]*2){ 
          # randomly select an observation
          cellN <- round(runif(1, 1, S_obs[s]))
          # fill it with a 1 
          foo[cellN, j] <- foo[cellN, j] + 1  
          # and so on until neighbours are all accounted for
        }
    }
    if(identical(summedSK[s, ], colSums(foo)) == F) message('Error in abundance tallies!')
    K_Nmat <- rbind(K_Nmat, foo)
  }
  K_Nmat <- K_Nmat[-1, ]  # remove the NA row
  if (nrow(K_Nmat) != sum(S_obs)) message('Error in total number of observations!')
  if (sum(colSums(K_Nmat) != colSums(summedSK)) > 1) message('Error in neighbour abundances')
  colnames(K_Nmat) <- paste0('K', 1:K)
  Spmax <- max(K_Nmat)
  #K_Nmat <-round((K_Nmat/max(K_Nmat))*100) #scale all the abundance between 0 and 100
  
  # Parameters for within guilt pairwise interaction
  #-------------------------------------------------
  # log of species-specific intrinsic performance
  sim_lambda <- runif(S, 2, 3)  
  
  # 'true' interaction strengths based on random draw
  #a <- rnorm( 8000, 0, .3 )
  #b <- rnorm( 2000, -0.7, .3 )
  #c <- rnorm( 1000,  0.3, .2 )
  #customdist <- c( a, b, c ) ## make a weird dist with Kurtosis and Skew
  #sim_truealpha <- matrix(data = sample(customdist, S*K, replace = F),
  #                       nrow = S, ncol = K)
  
  # 'true' interaction strengths based choosen values
  sim_genericalpha <- -0.3 # I put the same alpha generic for all focal 
  
  sim_alphaspecific <- matrix(nrow=S,ncol=K,
                              sample(rep(c(0,-0.2),times=50),S*K, replace = TRUE)) # I put specific alpha specific, eitehr 0 or 0.5
  
  # Within guilt HOIs abundance
  #-------------------------------
  
  K_Nmat_HOIs<- list()
  
  for(foc in 1:K){
    K_Nmat_species <- list()
    q <- K_Nmat[c((sum(S_obs[0:(foc-1)])+1):sum(S_obs[0:(foc)])),]
    
    for (obs in 1:S_obs[foc]){
      matrix_i <- matrix(nrow=S,ncol=S)
      for (n in 1:S) {
        for (m in 1:S) {
          if (m <= n){
            matrix_i[n,m] = 0
          }
          else{
            matrix_i[n,m] =  (q[obs,n]*q[obs,m])
          } 
        }
      }
      matrix_i[is.na(matrix_i)] <- 0
      K_Nmat_species[[obs]] <- matrix_i
    }
    K_Nmat_HOIs[[foc]] <- K_Nmat_species
  }
  # Parameters for within guilt HOIs
  #-------------------------------------------------
  
  
  
  # 'true' interaction strengths 
  HOIs_a <- rnorm( 8000, -0.1, 0.15)
  HOIs_b <- rnorm( 2000, 0.1, 0.15 )
  HOIs_c <- rnorm( 1000,  -0.05, 0.08 )
  customdist <- c(  HOIs_a, HOIs_b , HOIs_c ) ## make a weird dist with Kurtosis and Skew
  
  sim_HOIs_specific <- list()
  for(foc in c(1:S)){
    sim_HOIs_specific[[foc]] <- matrix(data = sample(rep(c(0,-0.1),times=50), K*K, replace = F),
                                       nrow = K, ncol = K)
  }
  
  sim_HOIs_generic  <- matrix(nrow=S,ncol=K)
  
  for(foc in c(1:S)){
    for(comp in c(1:K)){
      sim_HOIs_generic[foc,comp] <-  -0.1 # mean(sim_HOIs_specific[[foc]][comp,])
    }
  }
  
  
  # Parameters for HTL
  #-------------------
  
  # HTL abundance counts
  #-------------------------------
  if(HTL == T){
    # set up an empty matrix with columns for each HTL
    # we assume we have a different number of observations for each focal: S_obs  comuted above
    
    # we assume that S = HTL, and that the number of observations of HTL is relative to their abundance
    approx_HTL_obs <-   approx_K_obs
    summedSHTL <-   summedSK
    
    # randomly select observations that will be set to 0
    obs_to_rm_HTL <- sample(seq_along(1:(S* FvH)), pI_HTL*S* FvH)
    # remove unobserved interactions
    for(i in 1:(pI_HTL*S* FvH)) {
      summedSHTL[[obs_to_rm_HTL[i]]] <- 0
    }
    
    # set up an empty matrix with columns for each neighbours
    HTL_Nmat <- matrix(ncol =  FvH)
    for(s in 1:S){ 
      # create a matrix of observations for each focal species
      HTL_foo <- matrix(data = 0, nrow = S_obs[s], ncol =  FvH)
      for(h in 1: FvH) {
        if(summedSHTL[s, h] > 0)
          for (i in 1:summedSHTL[s, h]*2){ 
            # randomly select an observation
            HTL_cellN <- round(runif(1, 1, S_obs[s]))
            # fill it with a 1 
            HTL_foo[HTL_cellN, h] <- HTL_foo[HTL_cellN, h] + 1  
            # and so on until neighbours are all accounted for
          }
      }
      if(identical(summedSHTL[s , ], colSums(HTL_foo)) == F) message('Error in abundance tallies!')
      HTL_Nmat <- rbind(HTL_Nmat, HTL_foo)
    }
    HTL_Nmat <- HTL_Nmat[-1, ]  # remove the NA row
    if (nrow(HTL_Nmat) != sum(S_obs)) message('Error in total number of observations!')
    if (sum(colSums(HTL_Nmat) != colSums(summedSHTL)) > 1) message('Error in neighbour abundances')
    colnames(HTL_Nmat) <- paste0('HTL', 1: FvH)
    #HTL_Nmat <-round((HTL_Nmat/Spmax)*100) #scale all the abundance between 0 and 100 based on max interactions
    
  }
  # Parameters for HTL
  #-------------------
  if( HTL == T){
    # 'true' interaction strengths based choosen values
    sim_genericHTL <- - sim_genericalpha # make HTL less than direct interactions
    
    # 'true' HTL interaction strengths (a bit different than above)
    HTL_a <- rnorm( 8000, 0.1, 0.2 )
    HTL_b <- rnorm( 2000, 0.2, .05 )
    HTL_c <- rnorm( 1000,  0.1, 0.15 )
    HTL_customdist <- c( HTL_a, HTL_b, HTL_c ) ## make a weird dist with Kurtosis and Skew
    sim_specificHTL <- matrix(data = sample(rep(c(0,0.),times=50), S* FvH, replace = F),
                              nrow = S, ncol = FvH)
    
  }  
  # Observed seed set
  #------------------
  seeds <- rep(NA, length = sum(S_obs))
  counter <- 1
  
  for(s in 1:S) {
    
    for (i in (counter):(counter + S_obs[s] - 1)) {
      # print(i)
      # multiply neighbour abundances by 'true' alpha values
      if( HTL == T){ # add HTL effect
        beta_compk <- matrix(nrow=K,ncol=K)
        for(k in 1:K){
          beta_compk[k,] <-  (sim_HOIs_generic[s,k] + 
                                sim_HOIs_specific[[s]][k,])*K_Nmat_HOIs[[s]][[i - counter + 1]][k,]
        }
        seeds[i] <-  round(exp(sim_lambda[s] + sum(K_Nmat[i, ] * (sim_genericalpha +  sim_alphaspecific[s,])) +
                                 sum(beta_compk)  +
                                 sum(HTL_Nmat[i, ] * (sim_genericHTL + sim_specificHTL[s,]))))
      }else{
        beta_compk <- matrix(nrow=K,ncol=K)
        for(k in 1:K){
          beta_compk[k,] <- (sim_HOIs_generic[s,k] + sim_HOIs_specific[[s]][k,])*K_Nmat_HOIs[[s]][[i - counter + 1]][k,]
        }
        seeds[i] <- round(exp(sim_lambda[s] + sum(K_Nmat[i, ] * (sim_genericalpha +  sim_alphaspecific[s, ]))
                              + sum(beta_compk)))
      }
    }
    counter <- counter + S_obs[s]
  }
  
  # seeds <- round(log(seeds)) # to include facilitation we took the exponential of the above expression. 
  # To have the "observed" fecundity we need to take its natural logarithm
  # Simulated dataset
  #------------------
  focal <- do.call('c', mapply(rep, 1:S, S_obs, SIMPLIFY = F))    # focal identifier
  
  colnames(sim_alphaspecific) <- paste0('alpha', 1: K)
  row.names(sim_alphaspecific) <- paste0('focal', 1: S)
  if( HTL == T){ # add HTL effect
    colnames(sim_specificHTL) <- paste0('HTL', 1: FvH)
    row.names(sim_specificHTL) <- paste0('focal', 1: S)
    
    simdata <- cbind(focal, seeds, K_Nmat,HTL_Nmat)
    simdata <- as.data.frame(simdata)
    
    return(list(simdata= simdata, sim_lambda = sim_lambda,
                sim_alpha_generic= sim_genericalpha, sim_alpha_specific=sim_alphaspecific,
                sim_HOIs_generic=sim_HOIs_generic,sim_HOIs_specific=sim_HOIs_specific,
                sim_HTL_generic=sim_genericHTL,sim_HTL_specific=sim_specificHTL))
  }else{
    simdata <- cbind(focal, seeds, K_Nmat)
    simdata <- as.data.frame(simdata)
    return(list(simdata= simdata, sim_lambda = sim_lambda,
                sim_alpha_generic= sim_genericalpha, sim_alpha_specific=sim_alphaspecific,
                sim_HOIs_generic=sim_HOIs_generic,sim_HOIs_specific=sim_HOIs_specific))
  }
  
}



