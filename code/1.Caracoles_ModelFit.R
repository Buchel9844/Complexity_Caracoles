# This script will run the empirical model fits for each focal species and each
#       environmental covariate. A separate script will then make the empirical
#       figures for the manuscript
# https://mc-stan.org/docs/2_29/functions-reference/index.html#overview

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 1. SET UP: Import data, create df with abundances ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#rm(list = ls())
setwd("~/Eco_Bayesian/Complexity_caracoles")
#---- 1.1. Import packages ----
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
#system.file("~/home/lbuche/R/x86_64-pc-linux-gnu-library/4.1/rstan",package='rstan')
options(mc.cores = parallel::detectCores())

library(rstan)
library("codetools")
library(HDInterval)
library(matrixStats)
library(tidyverse)
library(stringr)
library(dplyr)
library(ggplot2) 

#---- 1.IMPORT DATA -----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#rm(list = ls()) # remove environment to not use the objects created before

#---- 1.1. Import the arguments----
args <- commandArgs(trailingOnly = TRUE)

#Extract argument from slurms
focal.list <- c("CHFU","HOMA","LEMA","CETE")[as.numeric(args[1])]
run.prelimfit <- 1 # as.numeric(args[1])
run.finalfit <- 1 #as.numeric(args[2])
#focal.levels <- c("CETE","HOMA","CHFU","LEMA)
summary.interactions <- data.frame() 
#focal="CHFU"
#year="2019"
#complexity.level=1
for( focal in focal.list){ # "CHFU","HOMA","LEMA,"CETE"
  for( year in c("2019",'2020','2021')){ # '2020','2021'"2018" "2019",'2020','2021'
    for( complexity.level in c(1:3)){ #complexity.level.int
      summary.interactions.n <- data.frame()
      complexity.plant <-c("class","family","code.plant")[complexity.level]
      complexity.animal <- c("group","family","species")[complexity.level]
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
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      source(paste0(home.dic,"code/ExtractMatrix_fct.R"))
      extract.matrix(focal,year=year,
                     complexity.plant ,
                     complexity.animal,
                     plant.class,SpData=SpData)
      
      #Returns two elemens
      #view(summary.interactions.n)
      #str(DataVec) #contains 21 elements
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
      
      rstan_options(auto_write = TRUE) 
      options(mc.cores = parallel::detectCores()) # to use the core at disposition 
      
      if(run.prelimfit == T){
        PrelimFit <- stan(paste0(home.dic,"code/Short_Caracoles_BH_FH_Preliminary.stan"), 
                          data = DataVec,
                          init = "random", 
                          control =list(max_treedepth=15),
                          warmup  = 1000, 
                          iter = 2000,
                          #init_r=2,
                          chains = 8,
                          seed= 1616)
        
        save(file=paste0(project.dic,"results/stan/PrelimFit",
                         FocalPrefix,"_",year,"_",
                         complexity.plant,"_",
                         complexity.animal,".rds"),
             PrelimFit)
      }else{
        load(paste0(project.dic,"results/stan/PrelimFit",FocalPrefix,"_",year,"_",
                    complexity.plant,"_",complexity.animal,".rds"))
      }
      PrelimPosteriors <- rstan::extract(PrelimFit)
      print("preliminary fit done")
      #---- 3.4. Extract the inclusion indice from the Preliminary Data----
      source(paste0(home.dic,"code/ExtractInclusion_fct.R"))
      extract.inclusion(summary.interactions.n, DataVec)
      view(summary.interactions)
      
      #---- 3.3. Preliminary fit posterior check and behavior checks---- 
      ##### Diagnostic plots
      pdf(paste0(home.dic,"figure/Prelimfit_",paste(FocalPrefix,year,complexity.plant,complexity.animal,sep="_"),".pdf"))
      title= paste0("Diagnostic plots for Prelimfit of ",focal, " in ", year," \nwhen plant complexity is ",complexity.plant,
                    " \nand HTL complexity is ",complexity.animal)
      # Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
      source("~/Eco_Bayesian/Test_simulation/code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
      # check the distribution of Rhats and effective sample sizes 
      stan_post_pred_check(PrelimPosteriors,"F_hat",DataVec$Fecundity, max(DataVec$Fecundity)+20)
      
      # N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
      hist(summary(PrelimFit)$summary[,"Rhat"], 
           main = paste("Prelim fit: Histogram of Rhat for",
                        FocalPrefix, " in ", year," \nwhen plant complexity is ",complexity.plant,
                        " \nand HTL complexity is ",complexity.animal))
      hist(summary(PrelimFit)$summary[,"n_eff"],
           main = paste("Prelim fit: Histogram of n_eff for",
                        FocalPrefix, " in ", year," \nwhen plant complexity is ",complexity.plant,
                        " \nand HTL complexity is ",complexity.animal))
      
      # Next check the correlation among key model parameters and identify any
      #       divergent transitions
      par <- c('lambdas','disp_dev','alpha_generic','alpha_intra',"gamma_FV_generic","gamma_H_generic")
      if(DataVec$RemoveFV == 1){par <- c('lambdas','disp_dev','alpha_generic',
                                         'alpha_intra',"gamma_H_generic")
      }
      if(DataVec$RemoveH == 1){par <- c('lambdas','disp_dev',
                                        'alpha_generic','alpha_intra',"gamma_FV_generic")
      }
      if(DataVec$RemoveH == 1 & DataVec$RemoveFV == 1){par <- c('lambdas','disp_dev',
                                                                'alpha_generic','alpha_intra')
      }
      
      stan_trace(PrelimFit,pars=par,inc_warmup = TRUE)
      stan_dens(PrelimFit,pars=par)
      stan_plot(PrelimFit,pars=par)
      
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
      set.seed(1616)
      std.error <- function(x) sd(x)/sqrt(length(x))
      
      list.init <- function(...)list(lambdas= array(as.numeric(rnorm(1,mean = mean(PrelimPosteriors$lambdas),
                                                                     sd= std.error(PrelimPosteriors$lambdas))), dim = 1),
                                     alpha_generic_tilde = array(as.numeric(rnorm(1,mean = mean(PrelimPosteriors$alpha_generic_tilde),
                                                                                  sd= std.error(PrelimPosteriors$alpha_generic_tilde))), dim = 1),
                                     alpha_intra_tilde = array(as.numeric(rnorm(1,mean = mean(PrelimPosteriors$alpha_intra_tilde),
                                                                                sd= std.error(PrelimPosteriors$alpha_intra_tilde))), dim = 1),
                                     disp_dev = array(as.numeric(rnorm(1,mean = mean(PrelimPosteriors$disp_dev),
                                                                       sd= std.error(PrelimPosteriors$disp_dev))), dim = 1))
      
      if(run.finalfit == T){
        options(mc.cores = parallel::detectCores()) # to use the core at disposition 
        FinalFit <- stan(file = paste0(home.dic,"code/Short_Caracoles_BH_Final.stan") ,
                         #fit= PrelimFit, 
                         data = DataVec.final,
                         init =list.init, # all initial values are 0 
                         control=list(max_treedepth=15),
                         warmup = 500,
                         iter = 1000, 
                         init_r = 2,
                         chains = 8,
                         seed= 1616) 
        
        save(file=paste0(project.dic,"results/stan/FinalFit",
                         FocalPrefix,"_",year,"_",complexity.plant,
                         "_",complexity.animal,".rds"),
             FinalFit)
      }else{
        load(paste0(project.dic,"results/stan/FinalFit",FocalPrefix,"_",year,"_",complexity.plant,"_",complexity.animal,".rds")) 
      }
      FinalPosteriors <- rstan::extract(FinalFit)
      library(rstanarm)
      library(rstan)
      library(rlang)
      #install.packages("/home/lbuche/Eco_Bayesian/punim1670/rstan_2.32.6.tar.gz", 
      #                 repos = NULL, type="source",
      #                 dependencies = TRUE, INSTALL_opts = '--no-lock')
      
      
      print("final fit is done")
      #---- 4.3. Final fit posterior check and behavior checks---- 
      
      ##### Diagnostic plots and post prediction 
      pdf(paste0(home.dic,"figure/FinalFit_",paste(FocalPrefix,year,complexity.plant,complexity.animal,sep="_"),".pdf"))
      # Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
      source("~/Eco_Bayesian/Test_simulation/code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
      # check the distribution of Rhats and effective sample sizes 
      ##### Posterior check
      
      stan_post_pred_check(FinalPosteriors,"F_hat",
                           DataVec$Fecundity,max(DataVec$Fecundity)+20) 
      
      # N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
      hist(summary(FinalFit)$summary[,"Rhat"],
           main = paste("Finat Fit: Histogram of Rhat for",
                        FocalPrefix, " in ", year," \nwhen plant complexity is ",complexity.plant,
                        " \nand HTL complexity is ",complexity.animal))
      hist(summary(FinalFit)$summary[,"n_eff"],
           main = paste("Finat Fit: Histogram of Rhat for",
                        FocalPrefix, " in ", year," \nwhen plant complexity is ",complexity.plant,
                        " \nand HTL complexity is ",complexity.animal))
      
      # plot the corresponding graphs
      par <- c('lambdas','disp_dev','alpha_generic','alpha_intra',"gamma_FV_generic","gamma_H_generic")
      if(DataVec$RemoveFV == 1){par <- c('lambdas','disp_dev','alpha_generic',
                                         'alpha_intra',"gamma_H_generic")
      }
      if(DataVec$RemoveH == 1){par <- c('lambdas','disp_dev',
                                        'alpha_generic','alpha_intra',"gamma_FV_generic")
      }
      if(DataVec$RemoveH == 1 & DataVec$RemoveFV == 1){par <- c('lambdas','disp_dev',
                                                                'alpha_generic','alpha_intra')
      }
      stan_trace(FinalFit , pars=par,
                 inc_warmup = TRUE)
      stan_dens(FinalFit , pars=par)
      stan_plot(FinalFit , pars=par)
      # Next check the correlation among key model parameters and identify any
      dev.off()
      
      #---- 4.4. Final fit analysis---
      
      #FileName <- paste("/home/lisavm/Simulations/", FocalPrefix, "_", "_FinalFit.rdata", sep = "")
      #FileName <- paste(home.dic,"results/", year,"_",complexity.plant,
      #                  complexity.animal,"_",FocalPrefix, "_FinalFit.rdata", sep = "")
      
      #save(DataVec.final,file = FileName)
      
      source(paste0(home.dic,"code/ExtractEstimates.R"))
      
      print(paste("ALL done for",FocalPrefix,year,
                  complexity.plant,complexity.animal,sep=" "))
      #  }
    }
  }
}
