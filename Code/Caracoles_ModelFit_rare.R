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
rstan_options(auto_write = TRUE)
summary.interactions <- data.frame()
# Run UploadData.R 2. and 3. or uncomment the following lines 
plant.class <- read.csv("data/plant_code.csv", sep=",")
years <- "2020"
focal <- "LEMA"
complexity.plant  <- "class"
complexity.animal <- "group"
#focal <- "CHFU"
#complexity  <- "family"
for (years in c("2017","2018","2019","2020","2021")){
for (focal in c("MESU","LEMA","HOMA","CHFU")){
for(complexity.plant  in c("code.plant","family","class")){ #add "code.plant" ,make sure it is the name of a column of plant.class
for(complexity.animal  in c("species","family","group")){ #add "code.plant" ,make sure it is the name of a column of plant.class
    
  # view complexity levels 
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

SpData <- na.omit(SpData) 
SpData <- subset(SpData, year == years)
SpData<- select(SpData,all_of(c("day","month", "year","plot","subplot" ,"focal",
                               "fruit","seed",complexitylevel.plant,FocalPrefix)))

SpDataFocal <- subset(SpData, focal == FocalPrefix )
#head(SpData)


if (complexity.minimize == T){
  levels.of.focal <- plant.class[which(plant.class$code.plant==FocalPrefix),complexity.plant]
  SpDataFocal[,levels.of.focal] <- SpDataFocal[,levels.of.focal] - SpDataFocal[,FocalPrefix]
}



SpDataherbivore <- subset(get(paste0("herbivore","_",complexity.animal)),year %in% years)
names(SpDataherbivore)[which(names(SpDataherbivore) == "plant")] <- "focal"

SpDataherbivore <- left_join(SpData,SpDataherbivore,by=c("month","year","plot",
                                            "subplot","focal"))

SpDatafloralvis <- subset(get(paste0("floral_visitor","_",complexity.animal)),year %in% years)
names(SpDatafloralvis)[which(names(SpDatafloralvis) == "plant")] <- "focal"

SpDatafloralvis <- left_join(SpData,SpDatafloralvis,by=c("month","year","plot",
                                                              "subplot","focal"))
view(SpDatafloralvis)
# determine the levels of complexity  animals
complexitylevel.herbivore <- names(SpDataherbivore)[which(!names(SpDataherbivore) %in% c(complexitylevel.plant,"focal","seed",
                                                                                         "fruit","day.x",FocalPrefix,
                                                                                         "day","month","year","plot",
                                                                                         "subplot","plant",
                                                                                         "code.plant"))]
complexitylevel.floral.visitor <- names(SpDatafloralvis)[which(!names(SpDatafloralvis) %in% c(complexitylevel.plant,"seed",
                                                                                              "fruit","day.x","focal",FocalPrefix,
                                                                                              "day","month","year","plot",
                                                                                              "subplot","plant",
                                                                                              "code.plant"))]



# Next continue to extract the data needed to run the model. 
N <- as.integer(nrow(SpDataFocal))
Fecundity <- as.integer(SpDataFocal$seed)  
plot <- as.integer(factor(as.factor(SpDataFocal$plot), levels = c("1","2","3","4","5","6","7","8","9")))

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
  }
}
SpNames <- AllSpNames[SpToKeep]

#assign(paste0("SpNames_",FocalPrefix),
 #     SpNames)
Intra <- ifelse(SpNames == FocalPrefix, 1, 0)

# Interaction matrix of herbivores with FOCAL
 
SpDataherbivore_focal <- subset(SpDataherbivore,focal == focal)
  
AllSpAbunds_herbivore <- SpDataherbivore_focal %>% 
  select(all_of(complexitylevel.herbivore))

SpTotals_herbivore <- colSums(AllSpAbunds_herbivore, na.rm = T)
SpToKeep_herbivore  <- SpTotals_herbivore  > 0

H <- sum(SpToKeep_herbivore)
SpMatrix_herbivore  <- matrix(NA, nrow = N, ncol = H)
i <- 1
for(s in 1:ncol(AllSpAbunds_herbivore)){
  if(SpToKeep_herbivore[s] == 1){
    SpMatrix_herbivore[,i] <- AllSpAbunds_herbivore[,s]
    i <- i + 1
  }
}
SpNames_herbivore <- names(SpToKeep_herbivore)[SpToKeep_herbivore]

# Interaction matrix of floral visitors with FOCAL

SpDatafloralvis_focal <- subset(SpDatafloralvis,focal == focal)

AllSpAbunds_floralvis<- SpDatafloralvis_focal %>% 
  select(all_of(complexitylevel.floral.visitor))

SpTotals_floralvis <- colSums(AllSpAbunds_floralvis, na.rm = T)
SpToKeep_floralvis <- SpTotals_floralvis  > 0

FV <- sum(SpToKeep_floralvis)
SpMatrix_floralvis  <- matrix(NA, nrow = N, ncol = FV)
i <- 1
for(s in 1:ncol(AllSpAbunds_floralvis)){
  if(SpToKeep_floralvis[s] == 1){
    SpMatrix_floralvis[,i] <- AllSpAbunds_floralvis[,s]
    i <- i + 1
  }
}
SpNames_floralvis <- names(SpToKeep_floralvis)[SpToKeep_floralvis]

summary.interactions <- bind_rows(summary.interactions,tibble(focal= focal,year=years,
                                                              complexity_plant=complexity.plant,
                                                              n.competitors_plant=length(SpNames)-1,
                                                              competitors_plant=list(SpNames[which(!SpNames %in% focal)]),
                                                              n.competitors_herbivore=length(SpNames_herbivore),
                                                              competitors_herbivore=list(SpNames_herbivore),
                                                              n.competitors_floralvisitor=length(SpNames_floralvis),
                                                              competitors_floralvisitor=list(SpNames_floralvis)))

# Interaction matrix of herbivores with COMPETITORS
SpDatafloralvis_comp <- select(SpDatafloralvis[which(!SpDatafloralvis$focal == focal),],
                               c("month","year","plot",
                                 "subplot","focal",all_of(complexitylevel.floral.visitor)))%>%
  gather( all_of(complexitylevel.floral.visitor),
            key = "floral.visitor", value = "visits") %>%
  unite(focal, floral.visitor,
      col = "plant_floral.visitor", sep = "_") 
SpDatafloralvis_comp <- aggregate(visits ~ plant_floral.visitor + subplot + plot + year + month, 
                                    SpDatafloralvis_comp, sum)
SpDatafloralvis_comp <-  spread(SpDatafloralvis_comp ,plant_floral.visitor, visits)



SpDatafloralvis_comp <- left_join(SpDataFocal,SpDatafloralvis_comp, 
                                  by=c('subplot', "plot", "year","month"))

View(SpDatafloralvis_comp)

AllSpAbunds_floralvis_comp <- SpDatafloralvis_comp  %>% 
  select(all_of(names(SpDatafloralvis_comp)[!names(SpDatafloralvis_comp) %in% c('subplot', "plot", "year","month",
                                                                                "day","focal","fruit","seed",focal, complexitylevel.plant)]))

SpTotals_floralvis_comp <- colSums(AllSpAbunds_floralvis_comp, na.rm = T)
SpToKeep_floralvis_comp <- SpTotals_floralvis_comp  > 0

FV_comp <- sum(SpToKeep_floralvis_comp)
SpMatrix_floralvis_comp  <- matrix(NA, nrow = N, ncol = FV_comp )
i <- 1
for(s in 1:ncol(AllSpAbunds_floralvis_comp)){
  if(SpToKeep_floralvis_comp[s] == 1){
    SpMatrix_floralvis_comp[,i] <- AllSpAbunds_floralvis_comp[,s]
    i <- i + 1
  }
}

SpNames_floralvis <- names(SpToKeep_floralvis)[SpToKeep_floralvis]

# Interaction matrix of floral visitors with COMPETITORS


# Set the parameters defining the regularized horseshoe prior, as described in
#       the "Incorporating sparsity-inducing priors" section of the manuscript.
tau0 <- 1
slab_scale <- sqrt(2)
slab_df <- 4

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
Inclusion_ij <- matrix(data = 0, nrow = 1, ncol = length(SpNames))
beta_Inclusion <- data.frame(nrow = length(SpNames), ncol = length(SpNames))
IntLevel <- 0.5 #0.5 usually, 0.75 for Waitzia, shade
is.list(PrelimPosteriors$beta_hat_ij)
#str(PrelimPosteriors$beta_hat_ijk)
#str(PrelimPosteriors$alpha_hat_ij)
  for(s in 1:length(SpNames)){
    # hdi : Calculate the highest density interval (HDI) for a probability distribution for a given probability mass
    Ints_ij <- HDInterval::hdi(PrelimPosteriors$alpha_hat_ij[,s], credMass = IntLevel)
    if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
      Inclusion_ij[1,s] <- 1
    }
    for(m in 1:S){
      beta_Ints_ij <- HDInterval::hdi(PrelimPosteriors$beta_hat_ijk[,s,m], credMass = IntLevel)
      
      if(beta_Ints_ij[1] > 0 | beta_Ints_ij[2] < 0){
        beta_Inclusion[s,m] <- 1
      }
    }
}
beta_Inclusion <- as.data.frame(beta_Inclusion)
names(beta_Inclusion) <- c(names(SpTotals[SpTotals!=SpToKeep]))
Inclusion_ij  <- as.data.frame(Inclusion_ij )
names(Inclusion_ij) <- c(names(SpTotals[SpTotals!=SpToKeep]))

#write.csv(beta_Inclusion,
#          paste("/home/lisavm/Simulations/beta_Inclusion_",FocalPrefix,".csv"))
#write.csv(Inclusion_ij,
 #         paste("/home/lisavm/Simulations/Inclusion_ij_",FocalPrefix,".csv"))

write.csv(beta_Inclusion,
          paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/results/beta_Inclusion_",years,"_",complexity,"_",FocalPrefix,".csv"))
write.csv(Inclusion_ij,
          paste0("/Users/lisabuche/Code/Project/HOI_Caracoles/results/Inclusion_ij_",years,"_",complexity,"_",FocalPrefix,".csv"))

  }
}}}

return(sum(Inclusion_ij))
return(sum(beta_Inclusion)) # 0 means that no specific HOIs are relevant overall years


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

