# This script will run the empirical model fits for each focal species and each
#       environmental covariate. A separate script will then make the empirical
#       figures for the manuscript

rm(list = ls())
remotes::install_version("ggplot2", version = "3.3.5", 
                         repos = "http://cran.us.r-project.org")  
#library(ggplot2)
library("rstan")
#install.packages("HDInterval")
library("HDInterval")
#library(here)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

FocalLetter <- "W" # "W" or "A"
FocalPrefix <- "WAAC" # "Waitzia" or "ARCA"
FocalSpecies <- "Waitzia.acuminata" # "Waitzia.acuminata" or "Arctotheca.calendula"

EnvCov <- "Phos" # "Phos" or "Shade"
EnvCol <- 71  # 72 for Canopy or 71 for Phosphorous

# Load in the data and subset out the current focal species.
SpData <- read.csv("data/water_full_env.csv")
str(SpData)
SpData <- subset(SpData, select = -c(X.NA., Seedcount.extrapolated.integer))
SpData <- na.omit(SpData) 

SpDataFocal <- subset(SpData, Focal.sp.x == FocalLetter)

# Next continue to extract the data needed to run the model. 
N <- as.integer(nrow(SpDataFocal))
Fecundity <- as.integer(SpDataFocal$Number.flowers.total)
reserve <- as.integer(as.factor(SpDataFocal$Reserve.x))
env <- as.vector(scale(SpDataFocal[,EnvCol]))

# Now calculate the total number of species to use for the model, discounting
#       any species columns with 0 abundance. Save a vector of the species names
#       corresponding to each column for easy matching later.
AllSpAbunds <- SpDataFocal[,10:69]
AllSpNames <- names(SpDataFocal[10:69])
SpTotals <- colSums(AllSpAbunds)
SpToKeep <- SpTotals > 0
S <- sum(SpToKeep)
SpMatrix <- matrix(NA, nrow = N, ncol = S)
i <- 1
for(s in 1:ncol(AllSpAbunds)){
  if(SpToKeep[s] == 1){
    SpMatrix[,i] <- AllSpAbunds[,s]
    i <- i + 1
  }
}
SpNames <- AllSpNames[SpToKeep]
Intra <- ifelse(SpNames == FocalSpecies, 1, 0)


# Set the parameters defining the regularized horseshoe prior, as described in
#       the "Incorporating sparsity-inducing priors" section of the manuscript.
tau0 <- 1
slab_scale <- sqrt(2)
slab_df <- 4
DataVec <- c("N", "S", "Fecundity", "reserve", "SpMatrix", "env", 
             "Intra", "tau0", "slab_scale", "slab_df")

# Now run a perliminary fit of the model to assess parameter shrinkage
PrelimFit <- stan(file = "Code/Topher_BH_FH_Preliminary.stan", data = DataVec, iter = 3000, 
                  chains = 3)
PrelimPosteriors <- extract(PrelimFit)

##### Diagnostic plots
# First check the distribution of Rhats and effective sample sizes 
# N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
hist(summary(PrelimFit)$summary[,"Rhat"])
hist(summary(PrelimFit)$summary[,"n_eff"])
# Next check the correlation among key model parameters and identify any
#       c
pairs(PrelimFit, pars = c("lambdas", "alpha_generic", "alpha_intra","beta_generic"))
# Finally, check for autocorrelation in the posteriors of key model parameters
# N.B.  ACF =  autocorrelation function = 
# coefficient of correlation between two values in a time series
# how to interpret ACF : https://medium.com/analytics-vidhya/interpreting-acf-or-auto-correlation-plot-d12e9051cd14
acf(PrelimPosteriors$lambdas[,1,1])
acf(PrelimPosteriors$lambdas[,1,2])
acf(PrelimPosteriors$lambdas[,2,1])
acf(PrelimPosteriors$lambdas[,2,2])
acf(PrelimPosteriors$alpha_generic[,1])
acf(PrelimPosteriors$alpha_generic[,2])
acf(PrelimPosteriors$alpha_intra[,1])
acf(PrelimPosteriors$alpha_intra[,2])
acf(PrelimPosteriors$beta_generic[,1])
acf(PrelimPosteriors$beta_generic[,2])

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
Inclusion_ij <- matrix(data = 0, nrow = 2, ncol = S)
Inclusion_eij <- matrix(data = 0, nrow = 2, ncol = S)
beta_Inclusion_ij <- matrix(data = 0, nrow = S, ncol = S)
beta_Inclusion_eij <- matrix(data = 0, nrow = S, ncol = S)
beta_Inclusion <-list()
IntLevel <- 0.5 #0.5 usually, 0.75 for Waitzia, shade
is.list(PrelimPosteriors$beta_hat_ij)
str(PrelimPosteriors$beta_hat_ijk)
str(PrelimPosteriors$alpha_hat_ij)
for(i in 1:2){
  for(s in 1:S){
    # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
    Ints_ij <- HDInterval::hdi(PrelimPosteriors$alpha_hat_ij[,i,s], credMass = IntLevel)
    Ints_eij <- HDInterval::hdi(PrelimPosteriors$alpha_hat_eij[,i,s], credMass = IntLevel)
    if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
      Inclusion_ij[i,s] <- 1
    }
    if(Ints_eij[1] > 0 | Ints_eij[2] < 0){
      Inclusion_eij[i,s] <- 1
    }
    
    for(m in 1:S){
    beta_Ints_ij <- HDInterval::hdi(PrelimPosteriors$beta_hat_ijk[,i,s,m], credMass = IntLevel)
    beta_Ints_eij <- HDInterval::hdi(PrelimPosteriors$beta_hat_eijk[,i,s,m], credMass = IntLevel)
    
    if(beta_Ints_eij[1] > 0 | beta_Ints_eij[2] < 0){
      beta_Inclusion_eij[s,m] <- 1
    }
    
    if(beta_Ints_ij[1] > 0 | beta_Ints_ij[2] < 0){
      beta_Inclusion_ij[s,m] <- 1
    }
    }
  }
  beta_Inclusion[[i]] <- list(beta_Inclusion_eij,beta_Inclusion_ij)
}

sum(Inclusion_ij)
sum(Inclusion_eij)
sum(beta_Inclusion[[1]][[1]]) + sum(beta_Inclusion[[2]][[1]]) # 0 means that no specific HOIs are relevant
sum(beta_Inclusion[[1]][[2]]) + sum(beta_Inclusion[[2]][[2]]) # 0 means that no specific HOIs are relevant

DataVec <- c("N", "S", "Fecundity", "reserve", "SpMatrix", "env", "Intra",
             "Inclusion_ij", "Inclusion_eij","beta_Inclusion_ij", "beta_Inclusion_eij")
FinalFit <- stan(file = "Code/Topher_BH_Final.stan", data = DataVec, iter = 3000, chains = 3)
FinalPosteriors <- extract(FinalFit)

# Diagnostic figures
hist(summary(FinalFit)$summary[,"Rhat"])
hist(summary(FinalFit)$summary[,"n_eff"])
pairs(FinalFit, pars = c("lambdas", "alpha_generic", "alpha_intra","beta_generic"))
acf(FinalPosteriors$lambdas[,1,1])
acf(FinalPosteriors$lambdas[,1,2])
acf(FinalPosteriors$lambdas[,2,1])
acf(FinalPosteriors$lambdas[,2,2])
acf(FinalPosteriors$alpha_generic[,1])
acf(FinalPosteriors$alpha_generic[,2])
acf(FinalPosteriors$alpha_intra[,1])
acf(FinalPosteriors$alpha_intra[,2])
acf(FinalPosteriors$beta_generic[,1])
acf(FinalPosteriors$beta_generic[,2])

FileName <- paste("Results/", FocalPrefix, "_", EnvCov, "_FinalFit.rdata", sep = "")
save(FinalFit, SpNames, N, S, Fecundity, reserve, SpMatrix, env, Inclusion_ij,
     Inclusion_eij,beta_Inclusion_ij,beta_Inclusion_eij, 
     tau0, slab_scale, slab_df, Intra, file = FileName)
