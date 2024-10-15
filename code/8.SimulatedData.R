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
home.dic <- "/Users/lisabuche/Documents/Projects/Complexity_caracoles/"
project.dic <- "/Users/lisabuche/Documents/Projects/Complexity_caracoles/"


# rm(list = ls()) # empty envi
print("S = 1, K = 5, H = 5, FV = 5, pI = 25% and pI_HTL = 65%")
source(paste0(home.dic,"code/8.1.simul_data_toolbox.R"))
set.seed(1698)
S = 1 # one focal species
K = 5 # 5 plant neighbours
H = 5 # 5 herbivores
FV = 5 # 5 floral visitor

simul_data_list <- simul_data( S = 1,   # number of focal groups / species
                               K = 5,   # number of neighbour of focals 
                               pI = 0.1,# proportion of interactions which are NOT observed 
                               HTL= T, # add higher trophic level 
                               Pol = 5,   # number of HTL, here floral visitor or not floral visitor
                               H = 5,   # number of HTL, here floral visitor or not floral visitor
                               pI_HTL = 0.1 # proportion of HTL interactions which are NOT observed 
)


save(file= paste0(home.dic,"results/TestFit_sparsity_simdata.RData"),
     simul_data_list)

simdata <- simul_data_list$simdata
str(simdata)
ggplot(simdata[which(simdata$focal=="i"),]) + 
  geom_density(aes(x=as.numeric(seeds))) 


summary.interactions <- data.frame()
summary.interactions.n <- data.frame()

focal <- "i"
focal.neighbours <- paste0("K",focal)
plant.neighbours <- paste0("K",1:5)
Pol.neighbours <- paste0("FV",1:5)
H.neighbours <- paste0("H",1:5)

# subset data set
SpData <- simdata[simdata$focal==focal,]
SpData[ , c(2:ncol(SpData))] <- apply(SpData[ , c(2:ncol(SpData))], 2,            # Specify own function within apply
                                      function(x) as.numeric(as.character(x)))

SpData$seeds <- as.integer(SpData$seeds)
SpData$seeds[is.na(SpData$seeds)] <- 0
# Extract the data needed to run the model. 
N <- as.integer(nrow(SpData))
Fecundity <- SpData$seeds 
U <- floor(log(mean(Fecundity)))
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
                Intra=Intra, tau0=tau0, 
                slab_scale=slab_scale, slab_df=slab_df,
                H=H,
                matrix_HOIs_ijh=matrix_HOIs_ijh,
                SpMatrix_H=SpMatrix_H,FV=FV,
                SpMatrix_FV=SpMatrix_FV,
                matrix_HOIs_ijf=matrix_HOIs_ijf)

library("codetools")
options(mc.cores = parallel::detectCores())

print("preliminary test fit for plant pairwise interactions model starts")

options(mc.cores = parallel::detectCores())
TestFit_sparsity <- stan(file = paste0(home.dic,"code/Caracoles_R_Preliminary.stan"), 
                         data = DataVec,
                         init="random", # all initial values are 0 
                         control=list(max_treedepth=15),
                         warmup = 500,
                         iter = 1000, 
                         chains = 4,
                         seed= 1616)

save(file= paste0(home.dic,"results/TestFit_sparsity_complex.rds"),
     TestFit_sparsity)

load( paste0(home.dic,"results/TestFit_sparsity_complex.rds"))

print("preliminary test fit for plant pairwise interactions model done")
Post_TestFit_sparsity <- rstan::extract(TestFit_sparsity)

#---- 2.2. Check model ----

# Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
# check the distribution of Rhats and effective sample sizes 
##### Posterior check
source("~/Eco_Bayesian/Test_simulation/code/stan_modelcheck_rem.R") # call the functions to check 
##### Diagnostic plots and post prediction 
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
stan_dens(TestFit_sparsity, pars=c('lambdas','disp_dev','alpha_generic',
                                   'alpha_intra',"gamma_FV_generic","gamma_H_generic"))
stan_plot(TestFit_sparsity, pars=c('lambdas','disp_dev','alpha_generic',
                                   'alpha_intra',"gamma_FV_generic","gamma_H_generic"))

dev.off()

print(max(summary(TestFit_sparsity)$summary[,"Rhat"], na.rm=T),
      min(summary(TestFit_sparsity)$summary[,"n_eff"], na.rm=T))
#---- 2.3. Analysis species specific inclusion ---- 

Inclusion_ij <- matrix(data = 0, nrow = 1, 
                       ncol = length(plant.neighbours),
                       dimnames=list(c(""),c(plant.neighbours)))

Inclusion_FV <- matrix(data = 0,nrow = 1, ncol = length(Pol.neighbours),
                       dimnames=list(c(""),c(Pol.neighbours)))

Inclusion_H <- matrix(data = 0,nrow = 1, ncol = length(H.neighbours),
                      dimnames=list(c(""),c(H.neighbours)))



beta_Inclusion_plant <- matrix(data = 0,nrow = length(plant.neighbours), 
                               ncol = length(plant.neighbours),
                               dimnames=list(c(plant.neighbours),
                                             c(plant.neighbours)))


IntLevel <- 0.2 #0.5 usually, 0.75 for Waitzia, shade
for(s in 1:length(plant.neighbours)){
  # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
  Ints_ij <- HDInterval::hdi(Post_TestFit_sparsity$alpha_hat_ij[,s], credMass = IntLevel)
  
  if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
    Inclusion_ij[1,s] <- 1
  }
  for(m in 1:S){
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

Results_Sparsity_sim <- list(Inclusion_ij=Inclusion_ij,
                             beta_Inclusion_plant=beta_Inclusion_plant,
                             Inclusion_FV=Inclusion_FV,
                             Inclusion_H=Inclusion_H,
                             beta_Inclusion_FV = matrix(0, ncol=FV, nrow=S), # no HOIs of pol tested
                             beta_Inclusion_H = matrix(0, ncol=H, nrow=S)) # no HOIs of H tested

save(file= paste0(home.dic,"results/TestFit_sparsity_complex.RData"),
     Results_Sparsity_sim )

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
set.seed(1616)
std.error <- function(x) sd(x)/sqrt(length(x))


Final_sparsity_model <- stan(file = paste0(home.dic,"code/Caracoles_R_Final.stan"), 
                             data = DataVec.final,
                             init="random", # all initial values are 0 
                             control=list(max_treedepth=12), #  increase the limit to avoid wrong sampling
                             iter = 1000, 
                             chains = 3)
print("final test fit for sparsity model done")
save(file= paste0(home.dic,"results/TestFit_sparsity_Final_HTLmodel.rds"),
     Final_sparsity_model )

Post_Final_sparsity_model  <- rstan::extract(Final_sparsity_model )

#---- 2.2. Check model ----

# Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
# check the distribution of Rhats and effective sample sizes 
##### Posterior check
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


Results_Sparsity_sim <- list(Inclusion_ij=Inclusion_ij,
                             beta_Inclusion_plant=beta_Inclusion_plant,
                             Inclusion_FV=Inclusion_FV,
                             Inclusion_H=Inclusion_H,
                             lambdas = Post_Final_sparsity_model$lambdas,
                             alpha_generic = Post_Final_sparsity_model$alpha_generic,
                             alpha_intra = Post_Final_sparsity_model$alpha_intra,
                             gamma_FV_generic = Post_Final_sparsity_model$gamma_FV_generic,
                             alpha_hat_ij = Post_Final_sparsity_model$alpha_hat_ij,
                             gamma_FV_hat_if = Post_Final_sparsity_model$gamma_FV_hat_if,
                             gamma_H_hat_ih = Post_Final_sparsity_model$gamma_H_hat_ih) # no HOIs of H tested

save(file= paste0(home.dic,"results/TestFit_sparsity_parameters.RData"),
     Results_Sparsity_sim )


