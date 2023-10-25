options(mc.cores = parallel::detectCores())
library(rstan)
library("codetools")
library(HDInterval)
library(matrixStats)
library(tidyverse)
library(stringr)
library(dplyr)
library(ggplot2) 
library(shinystan) 
library(Hmisc) 
#### predictive check for all model 
args <- commandArgs(trailingOnly = TRUE)
#args <- c(1,1,1)

#year <- as.numeric(args[4])
#focal <- args[3]
#FocalPrefix <- focal 
complexity.level <- as.numeric(args[3])
complexity.plant <-c("class","family","code.plant")[complexity.level]
complexity.animal <- c("group","family","species")[complexity.level]
#print(paste(focal,complexity.plant,complexity.animal,year))
# FocalPrefix <- focal
#Extract argument from slurms
run.prelimfit <- as.numeric(args[1])
run.finalfit <- as.numeric(args[2])
home.dic <- "/home/lbuche/Eco_Bayesian/Complexity_caracoles/"
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/"
ModelCheck <- NULL
postfecundity <- list()
pdf(paste0(home.dic,"figure/ModelsBehavior.pdf"))
title= paste0("Appendix  : Diagnostic plots for Final fit of all models")

for( focal in c("LEMA","CHFU","HOMA","CETE")){ # "CHFU","HOMA",
  for( year in c("2019",'2020','2021')){
    for( complexity.level in c(1:3)){
      if(complexity.level ==3 & focal=="HOMA") next
    FocalPrefix <- focal
    print(paste(focal,complexity.plant,complexity.animal,year))
    
    #---- 1.2. Import the raw data ----
    
    #setwd("~/Eco_Bayesian/Complexity_caracoles")
    home.dic <- "/home/lbuche/Eco_Bayesian/Complexity_caracoles/"
    project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/"
    
    #rm(list = ls())
    abundance.plant <- read.csv(paste0(home.dic,"data/abundance.csv"), sep=",")
    
    competition.plant <- read.csv(paste0(home.dic,"data/competition.csv"), sep=",")
    
    plant.class <- read.csv(paste0(home.dic,"data/plant_code.csv"),sep=",")
    
    herbivorypredator <- read.csv(paste0(home.dic,"data/herbivorypredator.csv"), sep=",")
    
    floral_visitor <- read.csv(paste0(home.dic,"data/floral_visitor.csv"), sep=",")
    
    groupinglevels <- read.csv(paste0(home.dic,"data/groupinglevels.csv"), sep=",")
    
    # Run UploadData.R 2. and 3. or uncomment the following lines 
    plant.class <- read.csv(paste0(home.dic,"data/plant_code.csv"), sep=",")
    
    SpData <- read.csv(paste0(home.dic,"data/SpData.csv"))
    
    #---- 2.EXTRACT MATRIXES OF ABUNDANCES -----

    source(paste0(home.dic,"code/ExtractMatrix_fct.R"))
    extract.matrix(focal,year=year,
                   complexity.plant ,
                   complexity.animal,
                   plant.class,SpData)
    #Returns two elemens
    #view(summary.interactions.n)
    #str(DataVec) #contains 21 elements
    #---- 3. FINAL FIT----
    #---- 3.1. Set up parameters ---- 
    
    
    load(paste0(project.dic,"results/stan/FinalFit",FocalPrefix,"_",year,"_",complexity.plant,"_",complexity.animal,".rds")) 
    
    FinalPosteriors <- rstan::extract(FinalFit)
    
    print("final fit is done")
    #---- 4.3. Final fit posterior check and behavior checks---- 
    
    ##### Diagnostic plots and post prediction 
    # Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
    source("~/Eco_Bayesian/Test_simulation/code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
    # check the distribution of Rhats and effective sample sizes 
    ##### Posterior check
    
    stan_post_pred_check_plot(FinalPosteriors,"F_hat",
                         DataVec$Fecundity,max(DataVec$Fecundity)+40,
                         title= paste0(focal," in year ",year," for complexity: ",complexity.plant)) 
    
    # N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
    
    
    hist(summary(FinalFit)$summary[,"Rhat"],
         main = paste("Finat Fit: Histogram of Rhat for",
                      FocalPrefix, " in ", year," \nwhen complexity is ",complexity.plant))
    hist(summary(FinalFit)$summary[,"n_eff"],
         main = paste("Finat Fit: Histogram of Rhat for",
                      FocalPrefix, " in ", year," \nwhen complexity is ",complexity.plant))
    
    # plot the corresponding graphs
    par <- c('lambdas','alpha_generic','alpha_intra',"gamma_FV_generic","gamma_H_generic")
    if(DataVec$RemoveFV == 1){par <- c('lambdas','alpha_generic',
                                       'alpha_intra',"gamma_H_generic")
    }
    if(DataVec$RemoveH == 1){par <- c('lambdas',
                                      'alpha_generic','alpha_intra',"gamma_FV_generic")
    }
    if(DataVec$RemoveH == 1 & DataVec$RemoveFV == 1){par <- c('lambdas',
                                                              'alpha_generic','alpha_intra')
    }
    
    trace <- stan_trace(FinalFit , pars=par,inc_warmup = TRUE)
    print(trace)
    dens <- stan_dens(FinalFit , pars=par)
    print(dens)
    plot.stan <- stan_plot(FinalFit , pars=par)
    print( plot.stan)
    parplot <- stan_par(FinalFit, par)
    print(parplot)
    pairsplot <-  pairs(FinalFit,pars=par)
    print(pairsplot)
    
  # Next check the correlation among key model parameters and identify any
    
    # make data frame with rhat 
    mc <- data.frame(focal =FocalPrefix, 
                     year =year,
                     complexity = complexity.plant,
                     Rhat = max(summary(FinalFit)$summary[,"Rhat"],na.rm =T),
                     Neff = min(summary(FinalFit)$summary[,"n_eff"],na.rm = T))
    ModelCheck <- bind_rows(ModelCheck,mc)
    remove(FinalFit)
    
    
   }
 }
}
dev.off()
#---- 4. Shiny App for specific model -----
library(shinystan)

focal = "LEMA"
year = "2020"
complexity.level <- 2
complexity.plant <-c("class","family","code.plant")[complexity.level]
complexity.animal <- c("group","family","species")[complexity.level]

#load(paste0(project.dic,"results/stan/FinalFit",FocalPrefix,"_",year,"_",complexity.plant,"_",complexity.animal,".rds")) 

#my_sso <- launch_shinystan(as.shinystan(FinalFit, model_name = "my_model",warmup = 500))
#launch_shinystan(my_sso)

