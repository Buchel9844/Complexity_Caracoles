simul_data <- function(
    # number of species has to be a multiple of 10
  S = 1,   # number of focal groups / species
  K = 5,   # number of neighbour focals 
  pI = 0.1,   # proportion of interactions which are NOT observed 
  HTL= T, # add higher trophic levels 
  Pol = 5,   # number of  floral visitor 
  H = 5,   # number of  herbivores 
  pI_HTL = 0.1, # proportion of HTL interactions which are NOT observed 
  
  ...
) {
  # Neighbourhood abundance counts
  #-------------------------------
  # we assume we have a different number of observations for each focal
  S_obs <- rep(200,times=S)
  K_obs <- abs(rnorm(5,90,10))
  # we assume that S = K, and that the number of observations of K is relative to their abundance
  approx_K_obs <- round(jitter(10*(S_obs), amount = 30))
  summedSK <- round(jitter(sapply(S_obs, '*', approx_K_obs)/1000, 20))
  
  # set up an empty matrix with columns for each neighbours
  # create a matrix of observations for each focal species
  K_Nmat <- matrix(ncol = K)
  for(s in 1:S){ 
    # create a matrix of observations for each focal species
    foo <- matrix(data = 0, nrow = S_obs[s], ncol = K)
    for(j in 1:K) {
      for (i in 1:  K_obs[j]){ 
        # randomly select an observation
        cellN <- round(runif(1, 1, S_obs))
        # fill it with a 1 
        foo[cellN, j] <-   foo[cellN, j] + 1  
        # and so on until neighbours are all accounted for
      }
    }
    K_Nmat <- rbind(K_Nmat, foo)
  }
  K_Nmat <- K_Nmat[-1,]
  colSums(K_Nmat)
  
  if (nrow(K_Nmat) != sum(S_obs)) message('Error in total number of observations!')
  colnames(K_Nmat) <- paste0("K",c(1: K))
  Spmax <- max(K_Nmat,na.rm=T)
  if(!Spmax < 15){print("start again spmax too low")}
  
  #K_Nmat <-round((K_Nmat/max(K_Nmat))*100) #scale all the abundance between 0 and 100
  
  # Parameters for within guilt pairwise interaction
  #-------------------------------------------------
  # log of species-specific intrinsic performance
  sim_lambda <- abs(rnorm(S, 2, 2)) 
  
  # 'true' interaction strengths based on random draw
  #a <- rnorm( 8000, 0, .3 )
  #b <- rnorm( 2000, -0.7, .3 )
  #c <- rnorm( 1000,  0.3, .2 )
  #customdist <- c( a, b, c ) ## make a weird dist with Kurtosis and Skew
  #sim_truealpha <- matrix(data = sample(customdist, S*K, replace = F),
  #                       nrow = S, ncol = K)
  
  # 'true' interaction strengths based choosen values
  sim_alpha_generic <- -0.4 # same for intra specific interactions
  sim_alpha_specific <- matrix(nrow=S,ncol=K,
                               sample(rep(c(0,-0.2),times=50),S*K, 
                                      replace = TRUE)) 
  sim_alpha_specific[1] <- 0 # no specific intraspecific interactions
  diag(sim_alpha_specific) <- - abs(diag(sim_alpha_specific))
  # Within guilt HOIs abundance
  #-------------------------------
  
  K_Nmat_HOIs<- list()
  
  for(foc in 1:S){
    K_Nmat_species <- list()
    q <- K_Nmat[c(1:S_obs[foc]),]
    
    for (obs in 1:S_obs[foc]){
      matrix_i <- matrix(nrow=K,ncol=K)
      for (n in 1:K) {
        for (m in 1:K) {
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
  
  sim_HOIs_specific <- list()
  for(foc in c(1:S)){
    sim_HOIs_specific[[foc]] <- matrix(data = sample(rep(c(0,-0.1),times=50), K*K, replace = F),
                                       nrow = K, ncol = K)
  }
  
  sim_HOIs_generic  <- matrix(nrow=S,ncol=K)
  
  for(foc in c(1:S)){
    for(comp in c(1:K)){
      sim_HOIs_generic[foc,comp] <-  0 # mean(sim_HOIs_specific[[foc]][comp,])
    }
  }
  # Pol abundance counts
  #-------------------------------
  Pol_obs <- round(rnorm(Pol,90,10))
  # set up an empty matrix with columns for each neighbours
  Pol_Nmat <- matrix(ncol =  Pol)
  for(s in 1:S){ 
    # create a matrix of observations for each focal species
    Pol_foo <- matrix(data = 0, nrow = S_obs[s], ncol =  Pol)
    for(h in 1: Pol) {
      for (i in 1:Pol_obs[h]){ 
        # randomly select an observation
        Pol_cellN <- round(runif(1, 1, S_obs[s]))
        # fill it with a 1 
        Pol_foo[Pol_cellN, h] <- Pol_foo[Pol_cellN, h] + 1  
        # and so on until neighbours are all accounted for
      }
    }
    Pol_Nmat <- rbind(Pol_Nmat, Pol_foo)
  }
  Pol_Nmat <- Pol_Nmat[-1, ]  # remove the NA row
  if (nrow(Pol_Nmat) != sum(S_obs)) message('Error in total number of observations!')
  if (sum(colSums(Pol_Nmat)) != sum(Pol_obs))  message('Error in neighbour abundances')
  colnames(Pol_Nmat) <- paste0('Pol', 1: Pol)
  #HTL_Nmat <-round((HTL_Nmat/Spmax)*100) #scale all the abundance between 0 and 100 based on max interactions
  
  
  # Parameters for Pol
  #-------------------
  
  # 'true' interaction strengths based choosen values
  sim_generic_Pol <-  abs(sim_alpha_generic) # make HTL less than direct interactions
  
  sim_specific_Pol <- matrix(data = sample(rep(c(0,0.2),times=50), S* Pol, replace = F),
                             nrow = S, ncol = Pol)
  
  
  
  # H abundance counts
  #-------------------------------
  H_obs <- round(rnorm(H,90,10))
  # set up an empty matrix with columns for each neighbours
  H_Nmat <- matrix(ncol =  H)
  for(s in 1:S){ 
    # create a matrix of observations for each focal species
    H_foo <- matrix(data = 0, nrow = S_obs[s], ncol =  H)
    for(h in 1: H) {
      for (i in 1:H_obs[h]){ 
        # randomly select an observation
        H_cellN <- round(runif(1, 1, S_obs[s]))
        # fill it with a 1 
        H_foo[H_cellN, h] <- H_foo[H_cellN, h] + 1  
        # and so on until neighbours are all accounted for
      }
    }
    H_Nmat <- rbind(H_Nmat, H_foo)
  }
  H_Nmat <- H_Nmat[-1, ]  # remove the NA row
  if (nrow(H_Nmat) != sum(S_obs)) message('Error in total number of observations!')
  if (sum(colSums(H_Nmat)) != sum(H_obs))  message('Error in neighbour abundances')
  colnames(H_Nmat) <- paste0('H', 1: H)
  #HTL_Nmat <-round((HTL_Nmat/Spmax)*100) #scale all the abundance between 0 and 100 based on max interactions
  
  
  # Parameters for H
  #-------------------
  
  # 'true' interaction strengths based choosen values
  sim_generic_H <-  abs(sim_alpha_generic) # make HTL less than direct interactions
  
  sim_specific_H <- matrix(data = sample(rep(c(0,0.2),times=50), S* H, replace = F),
                             nrow = S, ncol = H)
  
  # Observed seed set
  #------------------
  seeds <- rep(NA, length = sum(S_obs))
  counter <- 1
  
  for(s in 1:S) {
    
    for (i in (counter):(counter + S_obs[s] - 1)) {
      # print(i)
      # multiply neighbour abundances by 'true' alpha values
      beta_compk <- matrix(nrow=K,ncol=K)
      for(k in 1:K){
        beta_compk[k,] <-  (sim_HOIs_generic[s,k] + 
                              sim_HOIs_specific[[s]][k,])*K_Nmat_HOIs[[s]][[i - counter + 1]][k,]
      }
      seeds[i] <- round(exp(sim_lambda[s] + 
                              sum(K_Nmat[i, ] * (sim_alpha_generic +  sim_alpha_specific[s, ])) + 
                              sum(beta_compk) +
                              sum(Pol_Nmat[i, ] * (sim_generic_Pol + sim_specific_Pol[s,])) +
                              sum(H_Nmat[i, ] * (sim_generic_H + sim_specific_H[s,]))))
    }
    counter <- counter + S_obs[s]
  }
  
  # seeds <- round(log(seeds)) # to include facilitation we took the exponential of the above expression. 
  # To have the "observed" fecundity we need to take its natural logarithm
  # Simulated dataset
  #------------------
  focal <- do.call('c', mapply(rep, c("i"), S_obs, SIMPLIFY = F))    # focal identifier
  
  colnames(sim_alpha_specific) <- paste0('alpha', c(1:K))
  row.names(sim_alpha_specific) <- paste0('focal', c("i"))
  
  simdata <- cbind(focal, seeds, K_Nmat,Pol_Nmat,H_Nmat)
  simdata <- as.data.frame(simdata)
  return(list(simdata= simdata, sim_lambda = sim_lambda,
              sim_alpha_generic = sim_alpha_generic, sim_alpha_specific=sim_alpha_specific,
              sim_generic_Pol =  sim_generic_Pol, sim_specific_Pol = sim_specific_Pol,
              sim_generic_H =  sim_generic_H, sim_specific_H = sim_specific_H,
              sim_HOIs_generic=sim_HOIs_generic,sim_HOIs_specific=sim_HOIs_specific))
}