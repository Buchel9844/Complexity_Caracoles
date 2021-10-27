# This script will run the empirical model fits for each focal species and each
#       environmental covariate. A separate script will then make the empirical
#       figures for the manuscript

rm(list = ls())

#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
library(rstan)
#install.packages(HDInterval)
library("HDInterval")
library("tidyverse")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
#focal <- "CHFU"
#complexity  <- "family"
for ( focal in c("MESU","LEMA","HOMA","CHFU")){
  for(complexity  in c("family","class","rareORabundant")){ #add "code.plant" ,make sure it is the name of a column of plant.class
  # set complexity levels 
  complexity.minimize <- F
  print(focal)
  print(complexity)

FocalPrefix <-focal # "LEMA" or "HOMA" or "CHFU"
#FocalSpecies <- "Hordeum marinum" # "Leontodon maroccanus" or " Hordeum marinum" or "PUPA"
plant.class <- read.csv("data/plant_class.csv", sep=",")
# determine the level of complexity based on any levels of the plant_class dataframe
# code.plant or rareORabundant
# if levels of complexity is not code.plant then uncomment this line 

complexitylevel <- levels(as.factor(plant.class[,complexity])) # change at line 35 too !!!
print(complexitylevel)
# Load in the data and subset out the current focal species.
#SpData <- read.csv("/home/lisavm/Simulations/data/competition.csv")
SpData <- read.csv("/Users/lisabuche/Code/Project/HOI_Caracoles/data/competition.csv")
#view(SpData)

SpData <- na.omit(SpData) 
years <- c("2017","2018","2019")
SpDataFocal <- subset(SpData, focal == FocalPrefix & year %in% years)
#view(SpDataFocal)
SpDataFocal <- select(SpDataFocal,
  all_of(c("day","month", "year","plot","subplot" ,"focal","fruit","seed",complexitylevel,FocalPrefix)))
head(SpDataFocal)

if (complexity.minimize == T){
levels.of.focal <- plant.class[which(plant.class$code.plant==FocalPrefix),complexity]
SpDataFocal[,levels.of.focal] <- SpDataFocal[,levels.of.focal] - SpDataFocal[,FocalPrefix]
}
# Next continue to extract the data needed to run the model. 
N <- as.integer(nrow(SpDataFocal))
Fecundity <- as.integer(SpDataFocal$seed) # X fruit ? 
year <- as.integer(factor(as.factor(SpDataFocal$year), levels = c("2017", "2018", "2019")))
plot <- as.integer(factor(as.factor(SpDataFocal$plot), levels = c("1","2","3","4","5","6","7","8","9")))

# Now calculate the total number of species to use for the model, discounting
#       any species columns with 0 abundance. Save a vector of the species names
#       corresponding to each column for easy matching later.
AllSpNames <- names(SpDataFocal)[!names(SpDataFocal) %in% c("day","month","year","plot",
                                                            "subplot", "focal","fruit","seed")]
library(dplyr)

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
  }
}
SpNames <- AllSpNames[SpToKeep]
#assign(paste0("SpNames_",FocalPrefix),
 #     SpNames)
Intra <- ifelse(SpNames == FocalPrefix, 1, 0)


# Set the parameters defining the regularized horseshoe prior, as described in
#       the "Incorporating sparsity-inducing priors" section of the manuscript.
tau0 <- 1
slab_scale <- sqrt(2)
slab_df <- 4
nyear <- length(years)
DataVec <- c("N", "S", "Fecundity", "year", "SpMatrix",
             "Intra", "tau0", "slab_scale", "slab_df")

# Now run a perliminary fit of the model to assess parameter shrinkage
print("preliminary fit beginning")

PrelimFit <- stan(file = "/Users/lisabuche/Code/Project/HOI_Caracoles/code/Caracoles_BH_FH_Preliminary.stan", 
                  data = DataVec,
                  iter = 100, 
                  chains = 3)
PrelimPosteriors <- rstan::extract(PrelimFit)
print("prelimi nary fit done")
pdf(paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/figure/beta_Inclusion_",complexity,"_",FocalPrefix,".pdf"))
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
Inclusion_ij <- matrix(data = 0, nrow = nlevels(as.factor(year)), ncol = S+1)
beta_Inclusion <-data.frame()
IntLevel <- 0.5 #0.5 usually, 0.75 for Waitzia, shade
is.list(PrelimPosteriors$beta_hat_ij)
#str(PrelimPosteriors$beta_hat_ijk)
#str(PrelimPosteriors$alpha_hat_ij)
for(i in 1:nlevels(as.factor(year))){
  beta_Inclusion_ij <- matrix(data = 0, nrow = S, ncol = S)
  for(s in 1:S){
    # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
    Ints_ij <- HDInterval::hdi(PrelimPosteriors$alpha_hat_ij[,i,s], credMass = IntLevel)
    if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
      Inclusion_ij[i,s] <- 1
    }
    for(m in 1:S){
      beta_Ints_ij <- HDInterval::hdi(PrelimPosteriors$beta_hat_ijk[,i,s,m], credMass = IntLevel)
      
      if(beta_Ints_ij[1] > 0 | beta_Ints_ij[2] < 0){
        beta_Inclusion_ij[s,m] <- 1
      }
    }
  }
  beta_Inclusion_ij <- as.data.frame(beta_Inclusion_ij)
  beta_Inclusion_ij$year <- years[as.integer(levels(as.factor(year)))[i]]
  beta_Inclusion <- rbind(beta_Inclusion,beta_Inclusion_ij)
  Inclusion_ij[i,S+1]<- years[as.integer(levels(as.factor(year)))[i]]
  
}

names(beta_Inclusion) <- c(names(SpTotals[SpTotals!=SpToKeep]),"year")
Inclusion_ij  <- as.data.frame(Inclusion_ij )
names(Inclusion_ij) <- c(names(SpTotals[SpTotals!=SpToKeep]),"year")

#write.csv(beta_Inclusion,
#          paste("/home/lisavm/Simulations/beta_Inclusion_",FocalPrefix,".csv"))
#write.csv(Inclusion_ij,
 #         paste("/home/lisavm/Simulations/Inclusion_ij_",FocalPrefix,".csv"))

write.csv(beta_Inclusion,
          paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/results/beta_Inclusion_",complexity,"_",FocalPrefix,".csv"))
write.csv(Inclusion_ij,
          paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/results/Inclusion_ij_",complexity,"_",FocalPrefix,".csv"))

  }
}
return(sum(Inclusion_ij))
return(sum(beta_Inclusion[,1:nyear])) # 0 means that no specific HOIs are relevant overall years


DataVec <- c("N", "S", "Fecundity", "year", "SpMatrix","Intra",
             "Inclusion_ij", "beta_Inclusion")
FinalFit <- stan(file = "/home/lisavm/Simulations/code/Caracoles_BH_Final.stan", data = DataVec, iter = 1000, chains = 2)
FinalFit <- stan(file = "/Users/lisabuche/Code/Project/HOI_Caracoles/code/Caracoles_BH_Final.stan", 
                 data = DataVec,
                 iter = 100, 
                 chains = 2)

FinalPosteriors <- extract(FinalFit)

# Diagnostic figures
for (i in 1:nyear){
  return(acf(FinalPosteriors$lambdas[,i,1]))
  return(acf(FinalPosteriors$lambdas[,i,2]))
}
return(acf(FinalPosteriors$alpha_generic[,1]))
return(acf(FinalPosteriors$alpha_intra[,1]))
return(acf(FinalPosteriors$beta_generic[,1]))

FileName <- paste("/home/lisavm/Simulations/", FocalPrefix, "_", "_FinalFit.rdata", sep = "")
FileName <- paste("/Users/lisabuche/Code/Project/HOI_Caracoles/code/", FocalPrefix, "_", "_FinalFit.rdata", sep = "")

save(FinalFit, SpNames, N, S, Fecundity, year, SpMatrix, Inclusion_ij,
     beta_Inclusion_ij,
     tau0, slab_scale, slab_df, Intra, file = FileName)
str(FinalFit)

