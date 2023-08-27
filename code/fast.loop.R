#rm(list = ls())
setwd("~/Eco_Bayesian/Complexity_caracoles")
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

home.dic <- "/home/lbuche/Eco_Bayesian/Complexity_caracoles/"
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/"

for ( focal in c("CETE","LEMA")){ # "HOMA","CHFU","CETE","LEMA
  for (year in c("2018","2019","2020","2021")){ #2018","2019","2020","2021
    if(c(year == "2018") & focal == "CHFU") next
    if(c(year == "2020") & focal == "HOMA") next
    for( complexity.level in c(2)){
      complexity.animal <- c("group","family","species")[complexity.level]
      complexity.plant <-c("class","family","code.plant")[complexity.level]
      FocalPrefix = focal
      print(paste(FocalPrefix,year,
                  complexity.plant,complexity.animal,sep=" "))
      # extract species specific
      source(paste0(home.dic,"code/ExtractEstimates.R"))
    }
  }
}
