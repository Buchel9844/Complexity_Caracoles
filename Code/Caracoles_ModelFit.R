# This script will run the empirical model fits for each focal species and each
#       environmental covariate. A separate script will then make the empirical
#       figures for the manuscript
# https://mc-stan.org/docs/2_29/functions-reference/index.html#overview

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 1. SET UP: Import data, create df with abundances ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#rm(list = ls())

#---- 1.1. Import packages ----
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
#system.file("~/home/lbuche/R/x86_64-pc-linux-gnu-library/4.1/rstan",package='rstan')
#installed.packages('~/home/lbuche/R/x86_64-pc-linux-gnu-library/4.1/rstan') 
library(rstan)
#install.packages("matrixStats", dependency =T)
#install.packages("HDInterval")
#system.file("~/home/lbuche/R/x86_64-pc-linux-gnu-library/4.1/HDInterval",package='HDInterval')
#installed.packages('~/home/lbuche/R/x86_64-pc-linux-gnu-library/4.1/HDInterval')
library(HDInterval)
#library("HDInterval")
#install.packages("tidyverse")
#system.file("~/home/lbuche/R/x86_64-pc-linux-gnu-library/4.1/tidyverse",package='tidyverse')
#installed.packages('~/home/lbuche/R/x86_64-pc-linux-gnu-library/4.1/tidyverse')
library(tidyverse)
#library(tidyverse)
#install.packages("dplyr")
#system.file("~/home/lbuche/R/x86_64-pc-linux-gnu-library/4.1/dplyr",package='dplyr')
#installed.packages('~/home/lbuche/R/x86_64-pc-linux-gnu-library/4.1/dply')
library(dplyr)
#library(dplyr)
#system.file("~/usr/local/easybuild-2019/easybuild/software/mpi/gcc/10.2.0/openmpi/4.0.5/r/4.1.0/lib64/R/library/ggplot2",package='ggplot2')
#installed.packages('~/usr/local/easybuild-2019/easybuild/software/mpi/gcc/10.2.0/openmpi/4.0.5/r/4.1.0/lib64/R/library/ggplot2')
library(ggplot2) 
#library(ggplot2)
options(mc.cores = parallel::detectCores()) # to use the core at disposition 

rstan_options(auto_write = TRUE) 
#---- 1.2. Import the Data ----
#setwd("~/Eco_Bayesian/Complexity_caracoles")
home.dic <- "/home/lbuche/Eco_Bayesian/Complexity_caracoles/"
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/"

#rm(list = ls())
abundance.plant <- read.csv(paste0(home.dic,"data/abundance.csv"), sep=",")

competition.plant <- read.csv(paste0(home.dic,"data/competition.csv"), sep=",")

plant.class <- read.csv(paste0(home.dic,"data/plant_code.csv"),sep=",")

herbivorypredator <- read.csv(paste0(home.dic,"data/herbivorypredator.csv"), sep=",")

floral_visitor <- read.csv(paste0(home.dic,"data/floral_visitor.csv"), sep=",")

# Run UploadData.R 2. and 3. or uncomment the following lines 
plant.class <- read.csv(paste0(home.dic,"data/plant_code.csv"), sep=",")



#---- 1.3. Set parameters ----
summary.interactions <- data.frame()
years <- "2020"
focal <- "CHFU"
FocalPrefix<- "CHFU"
complexity.plant  <- "family"
complexity.animal <- "family"

#Extract argument from slurms
args <- commandArgs(trailingOnly = TRUE)
run.prelimfit <- as.numeric(args[1])
run.finalfit <- as.numeric(args[2])
#"2017",'2018',
#for ( focal in c("LEMA","HOMA","CHFU","CETE")){
#for (years in c('2018',"2019","2020","2021")){
#    for( complexity.level in c(1:2)){
#     complexity.animal <- c("group","family","species")[complexity.level]
#     complexity.plant <-c("class","family","code.plant")[complexity.level]
#  if(c(years == "2018") & focal == "CHFU") next
#if(c(years == "2020") & focal == "HOMA") next

#   }
#  }
#}
#{
#  {
#    {

#---- 1.IMPORT DATA -----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

summary.interactions.n <- data.frame()
print(paste(focal,complexity.plant,complexity.animal,years))
SpData <- read.csv(paste0(home.dic,"data/SpData.csv"))

#---- 2.EXTRACT MATRIXES OF ABUNDANCES -----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source(paste0(home.dic,"code/ExtractMatrix_fct.R"))
extract.matrix(focal,year=years,
               complexity.plant ,
               complexity.animal,
               plant.class,SpData)
#Returns two elemens
#view(summary.interactions.n)
#view(DataVec) #contains 21 elements

#---- 3. PRELIMINARY FIT -----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 3.1. Set up summary interactions df and parameters ---- 
#Create summary data frame to see all the potential interactions


# Set the parameters defining the regularized horseshoe prior, as described in
#       the "Incorporating sparsity-inducing priors" section of the manuscript.
tau0 <- 1
slab_scale <- sqrt(2)
slab_df <- 4

#---- 3.2. Run  preliminary fit ----

# Now run a preliminary fit of the model to assess parameter shrinkage
print("preliminary fit beginning")
options(mc.cores=parallel::detectCores())
#install.packages("codetools")
library("codetools")
library("rstan")
init.list = list(lambdas=rnorm(1,log(min(DataVec$Fecundity)),1))
rstan_options(auto_write = TRUE) 
if(run.prelimfit == T){
  PrelimFit <- stan(paste0(home.dic,"code/Caracoles_BH_FH_Preliminary.stan"), 
                    data = DataVec,
                    init = "random", # all initial values are 0 
                    control=list(max_treedepth=15),
                    warmup  = 500, 
                    iter = 1000,
                    #init_r=2,
                    chains = 8,
                    seed= 1616)
  
  save(file=paste0(project.dic,"results/PrelimFit",FocalPrefix,years,"_",complexity.plant,"_",complexity.animal,".rds"),
       PrelimFit)
}

#load(paste0(project.dic,"results/PrelimFit",FocalPrefix,"_",years,"_",complexity.plant,"_",complexity.animal,".rds"))


PrelimPosteriors <- rstan::extract(PrelimFit)
print("preliminary fit done")
#---- 3.4. Extract the inclusion indice from the Preliminary Data----
source(paste0(home.dic,"code/ExtractInclusion_fct.R"))
extract.inclusion(summary.interactions.n, DataVec)
view(summary.interactions)

#---- 3.3. Preliminary fit posterior check and behavior checks---- 
##### Diagnostic plots
pdf(paste0(home.dic,"figure/Prelimfit_",paste(FocalPrefix,years,complexity.plant,complexity.animal,sep="_"),".pdf"))
title= paste0("Diagnostic plots for Prelimfit of ",focal, " in ", years," \nwhen plant complexity is ",complexity.plant,
              " \nand HTL complexity is ",complexity.animal)
# Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
source("~/Eco_Bayesian/Test_simulation/code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
# check the distribution of Rhats and effective sample sizes 
stan_post_pred_check(PrelimPosteriors,"F_hat",DataVec$Fecundity, max(DataVec$Fecundity)+20)

# N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
hist(summary(PrelimFit)$summary[,"Rhat"], 
     main = paste("Prelim fit: Histogram of Rhat for",
                  FocalPrefix, " in ", years," \nwhen plant complexity is ",complexity.plant,
                  " \nand HTL complexity is ",complexity.animal))
hist(summary(PrelimFit)$summary[,"n_eff"],
     main = paste("Prelim fit: Histogram of n_eff for",
                  FocalPrefix, " in ", years," \nwhen plant complexity is ",complexity.plant,
                  " \nand HTL complexity is ",complexity.animal))

# Next check the correlation among key model parameters and identify any
#       divergent transitions
pairs(PrelimFit, pars = c("lambdas", "alpha_intra","alpha_generic_tilde",
                          "gamma_H","gamma_FV"))

for( i in 1:length(DataVec$SpNames)){
  pairs(PrelimFit, pars = c(paste0("alpha_hat_ij_tilde[",i,"]"), 
                            paste0("beta_plant_generic_tilde[",i,"]")))
  
}

# plot the traceplot and distribution graphs
stan_model_check(PrelimFit,
                 param =c('lambdas','alpha_generic_tilde',"alpha_hat_ij_tilde",
                          'beta_plant_generic_tilde','disp_dev'))


dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 4. FINAL FIT----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 4.1. Set up parameters ---- 

DataVec.final <- append(DataVec,DataVec.final)

#FinalFit <- stan(file =  "/home/lisavm/Simulations/code/Caracoles_BH_Final.stan", data = DataVec, iter = 1000, chains = 2)


#---- 4.2. Run final fit ---- 
print("final fit begins")
rstan_options(auto_write = TRUE) 
options(mc.cores=parallel::detectCores())
list.init <- function(...)list(lambdas= as.numeric(rnorm(2,mean = mean(PrelimPosteriors$lambdas),
                                                         sd= var(PrelimPosteriors$lambdas)^2)),
                               alpha_generic_tilde = as.numeric(rnorm(2,mean = mean(PrelimPosteriors$alpha_generic_tilde),
                                                                      sd= var(PrelimPosteriors$alpha_generic_tilde)^2)),
                               alpha_intra_tilde =as.numeric(rnorm(2,mean = mean(PrelimPosteriors$alpha_intra_tilde),
                                                                   sd= var(PrelimPosteriors$alpha_intra_tilde)^2)),
                               disp_dev =as.numeric(rnorm(2,mean = mean(PrelimPosteriors$disp_dev),
                                                          sd= var(PrelimPosteriors$disp_dev)^2)))

if(run.finalfit == T){
  FinalFit <- stan(file = paste0(home.dic,"code/Caracoles_BH_Final.stan") ,
                   #fit= PrelimFit, 
                   data = DataVec.final,
                   init =list.init, # all initial values are 0 
                   control=list(max_treedepth=15),
                   warmup = 500,
                   iter = 1000, 
                   init_r = 2,
                   chains = 8) 
  
  save(file=paste0(project.dic,"results/FinalFit",FocalPrefix,"_",years,"_",complexity.plant,"_",complexity.animal,".rds"),
       FinalFit)
}
load(paste0(project.dic,"results/FinalFit",FocalPrefix,"_",years,"_",complexity.plant,"_",complexity.animal,".rds")) 

FinalPosteriors <- rstan::extract(FinalFit)
print("final fit is done")
#---- 4.3. Final fit posterior check and behavior checks---- 

##### Diagnostic plots and post prediction 
pdf(paste0("figure/FinalFit_",paste(FocalPrefix,years,complexity.plant,complexity.animal,sep="_"),".pdf"))
# Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
source("~/Eco_Bayesian/Test_simulation/code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
# check the distribution of Rhats and effective sample sizes 
##### Posterior check

stan_post_pred_check(FinalPosteriors,"F_hat",Fecundity,max(Fecundity)+20) 

# N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
hist(summary(FinalFit)$summary[,"Rhat"],
     main = paste("Finat Fit: Histogram of Rhat for",
                  FocalPrefix, " in ", years," \nwhen plant complexity is ",complexity.plant,
                  " \nand HTL complexity is ",complexity.animal))
hist(summary(FinalFit)$summary[,"n_eff"],
     main = paste("Finat Fit: Histogram of Rhat for",
                  FocalPrefix, " in ", years," \nwhen plant complexity is ",complexity.plant,
                  " \nand HTL complexity is ",complexity.animal))

# plot the corresponding graphs
stan_model_check(FinalFit,
                 param =c('lambdas','alpha_generic_tilde','gamma_H_generic_tilde','gamma_FV_generic_tilde',
                          'beta_plant_generic_tilde','disp_dev'))

# Next check the correlation among key model parameters and identify any
pairs(FinalFit, pars = c("lambdas", "alpha_generic", "alpha_intra",
                         "gamma_FV","gamma_H"))

pairs(FinalFit, pars = c("alpha_generic", "alpha_intra","beta_plant_generic"))
pairs(FinalFit, pars = c("gamma_FV", "beta_FV_generic","beta_2FV_generic","beta_FvH_generic"))
pairs(FinalFit, pars = c("gamma_H", "beta_H_generic","beta_2H_generic","beta_FvH_generic"))


dev.off()

#---- 4.4. Final fit analysis---

#FileName <- paste("/home/lisavm/Simulations/", FocalPrefix, "_", "_FinalFit.rdata", sep = "")
FileName <- paste("results/", years,"_",complexity.plant,complexity.animal,"_",FocalPrefix, "_FinalFit.rdata", sep = "")

save(append(FinalFit, DataVec.final),file = FileName)



