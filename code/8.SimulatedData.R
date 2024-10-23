library(rstan)
rstan_options(auto_write = TRUE) 
library("HDInterval")
library("tidyverse")
library(dplyr)
library(ggpubr)
library(ggplot2)
library(ggthemes)
library(MASS)  # the glm.neg model requires this packagefor the fits to be performed and/or to converge	
library("knitr")
library("ggpubr")
library(gridExtra)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 1. SET UP: Simulate data, create df with abundances ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#setwd("~/Eco_Bayesian/Complexity_caracoles")
home.dic <- "/home/lbuche/Eco_Bayesian/Complexity_caracoles/"
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/"
inclusion.df <- NULL
Simulated.parameter.df <- read.csv(paste0(home.dic,"results/Simulated.parameter.df.csv"))
Simulated.parameter.df <- Simulated.parameter.df[,-1]

Generic.parameter.df <- read.csv(paste0(home.dic,"results/Generic.parameter.df.csv"))
Generic.parameter.df  <- Generic.parameter.df[,-1]

Specific.parameter.df <- read.csv(paste0(home.dic,"results/Specific.parameter.df.csv"))
Specific.parameter.df  <- Specific.parameter.df[,-1]
init <- max(Generic.parameter.df$sim)
for( sim.i in init:100){
  inclusion.df.n <- NULL
  # rm(list = ls()) # empty envi
  source("/home/lbuche/Eco_Bayesian/Test_simulation/code/simul_data.R")
  simul_data_list <- simul_data( S = 1,   # number of focal groups / species
                                 K = 5,   # number of neighbour of focals 
                                 pI = 0,# proportion of interactions which are NOT observed 
                                 HTL= T, # add higher trophic level 
                                 Pol = 5,   # number of HTL, here floral visitor or not floral visitor
                                 H = 5,   # number of HTL, here floral visitor or not floral visitor
                                 pI_HTL = 0 # proportion of HTL interactions which are NOT observed 
  )
  
  Simulated.parameter.df.n <- data.frame(lambda =simul_data_list$sim_lambda,
                                         alpha_generic = simul_data_list$sim_alpha_generic,
                                         alpha_intra = simul_data_list$sim_alpha_generic,
                                         FV_generic =simul_data_list$sim_generic_Pol,
                                         H_generic = simul_data_list$sim_generic_H,
                                         
                                         sim=sim.i)
  
  Simulated.parameter.df <- bind_rows(Simulated.parameter.df.n,Simulated.parameter.df)
  
  write.csv(Simulated.parameter.df,
            file=paste0(home.dic,"results/Simulated.parameter.df.csv"))
  
  summary.interactions <- data.frame()
  summary.interactions.n <- data.frame()
  
  focal <- "i"
  #
  # subset data set
  simdata <- simul_data_list$simdata
  SpData <- simdata[simdata$focal==focal,]
  SpData[ , c(2:ncol(SpData))] <- apply(SpData[ , c(2:ncol(SpData))], 2,            # Specify own function within apply
                                        function(x) as.numeric(as.character(x)))
  
  SpData$seeds <- as.integer(SpData$seeds)
  SpData$seeds[is.na(SpData$seeds)] <- 0
  SpData <- SpData%>%
    dplyr::filter(seeds < 3000)
  # Extract the data needed to run the model. 
  N <- as.integer(nrow(SpData))
  Fecundity <- SpData$seeds 
  plot(density(Fecundity ))
  U <- floor(log(mean(Fecundity)))
  if(U ==0){
    U <- 1
  }
  FV <- 5
  H <- 5
  K <- 5
  
  # Abundance matrix of plants for pairwise interactions 
  
  # Now calculate the total number of plant species to use for the model, discounting
  #       any species columns with 0 abundance. Save a vector of the species names
  #       corresponding to each column for easy matching later.
  AllSpNames <- names(SpData)[!names(SpData) %in% 
                                c("focal","seeds",paste0('Pol', 1: FV),
                                  paste0('H', 1: H))]
  
  AllSpAbunds <- SpData %>% 
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
  
  # determine the vector intra that designate with neighbour is the focal
  Intra <- ifelse(SpNames == "K1", 1, 0)
  
  # Abundance matrix of plants*plants for HOIs
  
  # creation of a matrix of S by S of the interaction jk in HOIs_ijk for plants
  matrix_HOIs_plant <- list()
  for (i in 1:N){
    matrix_i <- matrix(nrow=K,ncol=K)
    for (n in 1:K) {
      for (m in 1:K) {
        if (m <= n){
          matrix_i[n,m] = 0
        }
        else{
          matrix_i[n,m] = (SpMatrix[i,n]* SpMatrix[i,m]) 
        } 
      }
    }
    matrix_i[is.na(matrix_i)] <- 0
    matrix_HOIs_plant[[i]] <-  matrix_i
  }
  print(head(SpMatrix))
  print(matrix_HOIs_plant[[i]])
  
  # POLINATORS
  # Abundance matrix of HOIs of FV on pairwise interactions is empty as we do not consider them here
  
  FVSpNames <- names(SpData)[!names(SpData) %in% 
                               c("focal","seeds",paste0('K', 1: K),
                                 paste0('H', 1: H))]
  
  FVSpAbunds <- SpData %>% 
    dplyr::select(all_of(FVSpNames))
  
  FVSpTotals <- colSums(FVSpAbunds)
  FVSpToKeep <- FVSpTotals > 0
  
  FV <- sum(FVSpToKeep)
  SpMatrix_FV <- matrix(NA, nrow = N, ncol = FV)
  i <- 1
  for(s in 1:ncol(FVSpAbunds)){
    if(SpToKeep[s] == 1){
      SpMatrix_FV[,i] <- FVSpAbunds[,s]
      i <- i + 1
    }else{next}
  }
  
  FVMatrix <- matrix(0, nrow = K, ncol = FV)
  matrix_HOIs_ijf <- list()
  
  for(s in 1:N){
    matrix_HOIs_ijf[[s]] <-  FVMatrix
  }
  
  FVSpNames <- FVSpNames[FVSpToKeep]
  
  # HERBIVORES
  # Abundance matrix of HOIs of H on pairwise interactions is empty as we do not consider them here
  
  HSpNames <- names(SpData)[!names(SpData) %in% 
                              c("focal","seeds",paste0('K', 1: K),
                                paste0('Pol', 1: FV))]
  
  HSpAbunds <- SpData %>% 
    dplyr::select(all_of(HSpNames))
  
  HSpTotals <- colSums(HSpAbunds)
  HSpToKeep <- HSpTotals > 0
  
  H <- sum(HSpToKeep)
  SpMatrix_H <- matrix(NA, nrow = N, ncol = H)
  i <- 1
  for(s in 1:ncol(HSpAbunds)){
    if(SpToKeep[s] == 1){
      SpMatrix_H[,i] <- HSpAbunds[,s]
      i <- i + 1
    }else{next}
  }
  
  HMatrix <- matrix(0, nrow = K, ncol = H)
  matrix_HOIs_ijh <- list()
  
  for(s in 1:N){
    matrix_HOIs_ijh[[s]] <-  HMatrix
  }
  
  HSpNames <- HSpNames[HSpToKeep]
  
  simul_data_list_parameters <- list(N=N, S=S,FV=FV,H=H,
                                     Intra=Intra,
                                     Fecundity = Fecundity,
                                     SpMatrix =SpMatrix,
                                     SpMatrix_H =SpMatrix_H,
                                     SpMatrix_FV =SpMatrix_FV,
                                     matrix_HOIs_plant = matrix_HOIs_plant,
                                     matrix_HOIs_ijf = matrix_HOIs_ijf,
                                     matrix_HOIs_ijh = matrix_HOIs_ijh,
                                     SpNames=SpNames,
                                     simul_data_list)
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #---- 2. FIT PRELIMINARY FIT ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #---- 2.1. Stan model ----
  print("preliminary test fit for plant pairwise interaction model beginning")
  
  # Set the parameters defining the regularized horseshoe prior, as described in
  #       the "Incorporating sparsity-inducing priors" section of the manuscript.
  tau0 <- 1
  slab_scale <- sqrt(4)
  slab_df <- 4
  
  
  # Set the parameters defining the regularized horseshoe prior, as described in
  #       the "Incorporating sparsity-inducing priors" section of the manuscript.
  
  run_simple <- 0
  run_complex  <- 1
  run_estimation <- 1
  RemoveH <- 0
  RemoveFvH <- 0
  RemoveFV <- 0
  
  DataVec <- list(N=N,
                  S=S,
                  U=U,
                  RemoveFvH= RemoveFvH,
                  RemoveH = RemoveH,
                  RemoveFV= RemoveFV,
                  Fecundity=Fecundity, 
                  SpMatrix=SpMatrix,
                  matrix_HOIs_plant=matrix_HOIs_plant,
                  Intra=Intra, 
                  tau0=tau0, 
                  slab_scale=slab_scale, 
                  slab_df=slab_df,
                  H=H,
                  matrix_HOIs_ijh=matrix_HOIs_ijh,
                  SpMatrix_H=SpMatrix_H,FV=FV,
                  SpMatrix_FV=SpMatrix_FV,
                  matrix_HOIs_ijf=matrix_HOIs_ijf)
  
  library("codetools")
  options(mc.cores = parallel::detectCores())
  
  print("preliminary test fit for plant pairwise interactions model starts")
  
  options(mc.cores = parallel::detectCores())
  TestFit_sparsity <- stan(file = paste0(home.dic,"code/Short_Caracoles_BH_FH_Preliminary.stan"), 
                           data = DataVec,
                           init="random", # all initial values are 0 
                           control=list(max_treedepth=15),
                           warmup = 300,
                           iter = 600, 
                           chains = 3)
  
  print("preliminary test fit for plant pairwise interactions model done for",)
  Post_TestFit_sparsity <- rstan::extract(TestFit_sparsity)
  
  #---- 2.2. Check model ----
  if( sim.i ==1){
    save(file= paste0(home.dic,"results/TestFit_sparsity_complex.rds"),
         TestFit_sparsity)
    load( paste0(home.dic,"results/TestFit_sparsity_complex.rds"))
    Post_TestFit_sparsity <- rstan::extract(TestFit_sparsity)
    
    # Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
    # check the distribution of Rhats and effective sample sizes 
    ##### Posterior check
    source("~/Eco_Bayesian/Test_simulation/code/stan_modelcheck_rem.R") # call the functions to check##### Diagnostic plots and post prediction 
    pdf(paste0(home.dic,"figure/PrelimFit_sparsity_complex.pdf"))
    
    stan_post_pred_check(Post_TestFit_sparsity ,"F_hat",Fecundity,200) 
    
    stan_model_check(TestFit_sparsity,
                     params=c('lambdas','disp_dev','alpha_generic',
                              'alpha_intra',"gamma_FV_generic","gamma_H_generic"))
    
    # N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
    hist(summary(TestFit_sparsity)$summary[,"Rhat"], 
         main = "Prelim fit: Histogram of Rhat \n
     Plant Pairwise interaction model",
     xlab = "Rhat")
    hist(summary(TestFit_sparsity)$summary[,"n_eff"],
         main = "Prelim fit: Histogram of Neff \n
     Plant Pairwise interaction model",
     xlab = "Neff")
    
    stan_trace(TestFit_sparsity, pars=c('lambdas','disp_dev','alpha_generic',
                                        'alpha_intra',"gamma_FV_generic","gamma_H_generic"),
               inc_warmup = TRUE)
    stan_trace(TestFit_sparsity, pars=c("beta_plant_hat_ijk"),
               inc_warmup = TRUE)
    stan_dens(TestFit_sparsity, pars=c('lambdas','disp_dev','alpha_generic',
                                       'alpha_intra',"gamma_FV_generic","gamma_H_generic"))
    stan_plot(TestFit_sparsity, pars=c('lambdas','disp_dev','alpha_generic',
                                       'alpha_intra',"gamma_FV_generic","gamma_H_generic"))
    
    dev.off()
  }
  
  #---- 2.3. Analysis species specific inclusion ---- 
  
  Inclusion_ij <- matrix(data = 0, nrow = 1, 
                         ncol = DataVec$S,
                         dimnames=list(c(""),SpNames))
  
  Inclusion_FV <- matrix(data = 0,nrow = 1, ncol = length(FVSpNames),
                         dimnames=list(c(""),c(FVSpNames)))
  
  Inclusion_H <- matrix(data = 0,nrow = 1, ncol = length(HSpNames),
                        dimnames=list(c(""),c(HSpNames)))
  
  beta_Inclusion_plant <- matrix(data = 0,nrow = length(SpNames), 
                                 ncol = length(SpNames),
                                 dimnames=list(SpNames,
                                               SpNames))
  
  
  IntLevel <- 0.3#0.5 usually, 0.75 for Waitzia, shade
  focal.neighbours <- paste0("K",focal)
  plant.neighbours <- SpNames
  Pol.neighbours <- FVSpNames
  H.neighbours <- HSpNames
  
  for(s in 1:length(plant.neighbours)){
    # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
    Ints_ij <- HDInterval::hdi(Post_TestFit_sparsity$alpha_hat_ij[,s], credMass = IntLevel)
    
    if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
      Inclusion_ij[1,s] <- 1
    }
    for(m in 1:length(plant.neighbours)){
      beta_Ints_ijk <- HDInterval::hdi(Post_TestFit_sparsity$beta_plant_hat_ijk[,s,m], credMass = IntLevel)
      
      if(beta_Ints_ijk[1] > 0 | beta_Ints_ijk[2] < 0){
        beta_Inclusion_plant[s,m] <- 1
      }
    }
  }
  
  #For Pollinator as second actor
  for(s in 1:length(Pol.neighbours)){
    # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
    Ints_if <- HDInterval::hdi(Post_TestFit_sparsity$gamma_FV_hat_if[,s], credMass = IntLevel)
    
    if(Ints_if[1] > 0 | Ints_if[2] < 0){
      Inclusion_FV[1,s] <- 1
    }
  }
  
  #For Herbivore as second actor
  for(s in 1:length(H.neighbours)){
    # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
    Ints_ih <- HDInterval::hdi(Post_TestFit_sparsity$gamma_H_hat_ih[,s], credMass = IntLevel)
    
    if(Ints_ih[1] > 0 | Ints_ih[2] < 0){
      Inclusion_H[1,s] <- 1
    }
  }
  
  Inclusion_ij.df <- Inclusion_ij %>%
    as.data.frame() %>%
    gather(key="neigh","inclusion") %>%
    mutate(parameter="plant.plant")
  
  Inclusion_ip.df <-Inclusion_FV %>%
    as.data.frame() %>%
    gather(key="neigh","inclusion") %>%
    mutate(parameter="plant.pol")
  Inclusion_ih.df <-Inclusion_H %>%
    as.data.frame() %>%
    gather(key="neigh","inclusion") %>%
    mutate(parameter="plant.her")
  Inclusion_ijk.df <- beta_Inclusion_plant  %>%
    as.data.frame() %>%
    rownames_to_column(var="neigh") %>%
    gather(c(paste0("K",1:5)),
           key="neigh.2",value="inclusion")%>%
    mutate(parameter="HOI.plant")
  
  
  inclusion.df.n <- bind_rows(Inclusion_ij.df,Inclusion_ip.df,
                              Inclusion_ih.df,Inclusion_ijk.df) %>%
    mutate(sim = sim.i) %>%
    left_join(simul_data_list$df.param)
  inclusion.df <- bind_rows(inclusion.df,inclusion.df.n)
  
  if(sim.i==1){
    Results_Sparsity_sim <- list(Inclusion_ij=Inclusion_ij,
                                 beta_Inclusion_plant=beta_Inclusion_plant,
                                 Inclusion_FV=Inclusion_FV,
                                 Inclusion_H=Inclusion_H) 
    save(file= paste0(home.dic,"results/TestFit_sparsity_complex.RData"),
         Results_Sparsity_sim )
  }
  Results_Sparsity_sim <- list(Inclusion_ij=Inclusion_ij,
                               beta_Inclusion_plant=beta_Inclusion_plant,
                               Inclusion_FV=Inclusion_FV,
                               Inclusion_H=Inclusion_H,
                               beta_Inclusion_FV = matrix(0, ncol=FV, nrow=S), # no HOIs of pol tested
                               beta_Inclusion_H = matrix(0, ncol=H, nrow=S)) # no HOIs of H tested
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #---- 3. FIT FINAL FIT ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #----3.1. Stan model ----
  
  DataVec.final <- append(DataVec,Results_Sparsity_sim)
  
  run_simple <- 0
  run_complex  <- 1
  run_estimation <- 1
  
  library("codetools")
  options(mc.cores = parallel::detectCores())
  std.error <- function(x) sd(x)/sqrt(length(x))
  
  
  Final_sparsity_model <- stan(file = paste0(home.dic,"code/Short_Caracoles_BH_Final.stan"), 
                               data = DataVec.final,
                               init="random", # all initial values are 0 
                               control=list(max_treedepth=12), #  increase the limit to avoid wrong sampling
                               warmup=300,
                               iter = 600, 
                               chains = 3)
  print("final test fit for sparsity model done")
  
  Post_Final_sparsity_model  <- rstan::extract(Final_sparsity_model )
  
  #---- 2.2. Check model ----
  if( sim.i==1){
    # Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
    # check the distribution of Rhats and effective sample sizes 
    ##### Posterior check
    save(file= paste0(home.dic,"results/TestFit_sparsity_Final_HTLmodel.rds"),
         Final_sparsity_model )
    source("~/Eco_Bayesian/Test_simulation/code/stan_modelcheck_rem.R") # call the functions to check 
    ##### Diagnostic plots and post prediction 
    pdf(paste0(home.dic,"figure/FinalFit_sparsity_complex.pdf"))
    
    stan_post_pred_check(Post_Final_sparsity_model ,"F_hat",Fecundity,300) 
    
    
    # N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
    hist(summary(Final_sparsity_model)$summary[,"Rhat"], 
         main = "Prelim fit: Histogram of Rhat \n
     Plant Pairwise interaction model",
     xlab = "Rhat")
    hist(summary(Final_sparsity_model)$summary[,"n_eff"],
         main = "Prelim fit: Histogram of Neff \n
     Plant Pairwise interaction model",
     xlab = "Neff")
    # plot the corresponding graphs
    #stan_model_check(FinalFit,
    #                param =c('lambdas','c','alpha_initial','alpha_slope','c'))
    # Next check the correlation among key model parameters and identify any
    #pairs(FinalFit, pars = c("lambdas",'alpha_initial','alpha_slope','c'))
    
    # functions from Rstan pacakges
    
    stan_trace(Final_sparsity_model , pars=c('lambdas','disp_dev',
                                             'alpha_generic',
                                             'alpha_intra',
                                             "gamma_FV_generic",
                                             "gamma_H_generic",
                                             "alpha_hat_ij",
                                             "gamma_FV_hat_if",
                                             "gamma_H_hat_ih"),
               inc_warmup = TRUE)
    stan_dens(Final_sparsity_model , pars=c('lambdas','disp_dev',
                                            'alpha_generic',
                                            'alpha_intra',
                                            "gamma_FV_generic",
                                            "gamma_H_generic",
                                            "alpha_hat_ij",
                                            "gamma_FV_hat_if",
                                            "gamma_H_hat_ih"))
    stan_plot(Final_sparsity_model , pars=c('lambdas','disp_dev',
                                            'alpha_generic',
                                            'alpha_intra',
                                            "gamma_FV_generic",
                                            "gamma_H_generic",
                                            "alpha_hat_ij",
                                            "gamma_FV_hat_if",
                                            "gamma_H_hat_ih"))
    
    dev.off()
  }
  
  Generic.parameter.df.n <- data.frame(lambda = median(Post_Final_sparsity_model$lambdas),
                                       alpha_generic = median(Post_Final_sparsity_model$alpha_generic),
                                       alpha_intra = median(Post_Final_sparsity_model$alpha_intra) ,
                                       FV_generic = median(Post_Final_sparsity_model$gamma_FV_generic),
                                       H_generic = median(Post_Final_sparsity_model$gamma_H_generic),
                                       sim=sim.i,
                                       Rhat=max(summary(Final_sparsity_model)$summary[,"Rhat"],na.rm=T),
                                       Neff =max(summary(Final_sparsity_model)$summary[,"n_eff"],na.rm=T))
  Generic.parameter.df <- bind_rows(Generic.parameter.df,Generic.parameter.df.n)
  write.csv(Generic.parameter.df,
            file=paste0(home.dic,"results/Generic.parameter.df.csv"))
  
  Alpha_specific.df <- Post_Final_sparsity_model$alpha_hat_ij %>%
    as.data.frame() %>%
    summarise_all(median) %>%
    gather(key="neigh",value="median.estimated") %>%
    mutate(neigh = paste0("K",1:5)) %>%
    mutate(parameter="plant.plant")
  
  Alpha_ip_specific.df <- Post_Final_sparsity_model$gamma_FV_hat_if%>%
    as.data.frame() %>%
    summarise_all(median) %>%
    gather(key="neigh",value="median.estimated") %>%
    mutate(neigh = FVSpNames) %>%
    mutate(parameter="plant.pol")
  
  Alpha_ih_specific.df <- Post_Final_sparsity_model$gamma_H_hat_ih%>%
    as.data.frame() %>%
    summarise_all(median) %>%
    gather(key="neigh",value="median.estimated") %>%
    mutate(neigh = HSpNames) %>%
    mutate(parameter="plant.her")
  Beta_specific.df <- Post_Final_sparsity_model$beta_plant_hat_ijk %>%
    as.data.frame() %>%
    summarise_all(median) %>%
    gather(key="neigh1.2",value="median.estimated") %>%
    mutate(neigh = paste0("K",rep(1:5,time=5)),
           neigh.2 =  paste0("K",rep(1:5,each=5)),
           parameter="HOI.plant") %>%
    dplyr::select(-"neigh1.2")
  
  Specific.parameter.df.n <- bind_rows(Alpha_specific.df,Alpha_ip_specific.df ,
                                       Alpha_ih_specific.df,Beta_specific.df) %>%
    left_join(inclusion.df.n)
  
  Specific.parameter.df <- bind_rows(Specific.parameter.df,Specific.parameter.df.n)
  
  write.csv(Specific.parameter.df,
            file=paste0(home.dic,"results/Specific.parameter.df.csv"))
  
  
}


Specific.parameter.df <- read.csv(paste0(home.dic,"results/Specific.parameter.df.csv"))
#view(Specific.parameter.df)
summary(Specific.parameter.df)
prob.detection.df <- Specific.parameter.df %>%
  mutate(true.detection = case_when((inclusion==1 & !value==0) ~ 1,
                                    T~0)) %>%
  aggregate(true.detection ~ parameter, sum) %>%
  left_join(Specific.parameter.df %>%
              mutate(should.detected = case_when((!value==0) ~ 1,
                                                 T~0)) %>%
              aggregate(should.detected  ~ parameter, sum)) %>%
  left_join(Specific.parameter.df %>%
              mutate(should.NOT.detected = case_when((inclusion==1 & value==0) ~ 1,
                                                     T~0)) %>%
              aggregate(should.NOT.detected  ~ parameter, sum)) %>%
  left_join(Specific.parameter.df %>%
              mutate(total.param = 1) %>%
              aggregate(total.param  ~ parameter, sum)) %>%
  mutate(probability.of.detection= true.detection /should.detected,
         probability.of.WRONG.detection= should.NOT.detected/total.param)

prob.detection.df
Specific.parameter.plot <- Specific.parameter.df %>%
  filter(inclusion ==1 & !value==0) %>%
  ggplot(aes(x=median.estimated - value,y=parameter)) +
  geom_density_ridges(quantile_lines = TRUE,scale=1,
                      quantiles = c(0.025,0.5, 0.975)) +
  geom_vline(xintercept=0,color="red") +
  theme_clean()
Specific.parameter.plot 


Simulated.parameter.df <- read.csv(paste0(home.dic,"results/Simulated.parameter.df.csv"))
Simulated.parameter.df <- Simulated.parameter.df[,-1]

Generic.parameter.df <- read.csv(paste0(home.dic,"results/Generic.parameter.df.csv"))
Generic.parameter.df  <- Generic.parameter.df[,-1]
#view(Generic.parameter.df)
summary(Generic.parameter.df)
generic.plot <- Generic.parameter.df %>%
  gather(alpha_generic,alpha_intra,FV_generic,
         H_generic, lambda,
         key="parameter",
         value="estimate") %>%
  left_join(Simulated.parameter.df %>%
              gather(alpha_generic,alpha_intra,FV_generic,
                     H_generic,lambda,
                     key="parameter",value="true.value"),
            by=c("sim","parameter"),
            multiple = "all") %>%
  mutate(diff = estimate -true.value) %>%
  filter(!parameter=="lambda") %>%
  ggplot(aes(x=estimate-true.value,y=parameter)) +
  geom_density_ridges(quantile_lines = TRUE,scale=1, 
                      quantiles = c(0.025,0.5, 0.975)) +
  geom_vline(xintercept=0,color="red") +
  theme_clean()
generic.plot

