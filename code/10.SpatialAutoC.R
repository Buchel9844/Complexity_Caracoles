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

library(sf)
library(spdep)
#toBibtex(citation(package="spdep") )
#---- 1.IMPORT DATA -----
#~~~~~~~~~~~~~~~~~~~~~~~~~~
home.dic <- "/home/lbuche/Eco_Bayesian/Complexity_caracoles/"
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/"
moran.df <- read.csv(paste0(home.dic, "results/Moran.df.csv"))
moran.df <- moran.df[,-1]
moran.df <- NULL
for( focal in c("CETE","HOMA","LEMA","CHFU")){ #"CETE","HOMA","LEMA","CHFU"
  for( year in c("2019",'2020','2021')){ # '2020','2021'"2018" "2019",'2020','2021'
    for( complexity.level in c(1:3)){ #complexity.level.int
      complexity.plant <-c("class","family","code.plant")[complexity.level]
      complexity.animal <- c("group","family","species")[complexity.level]
      FocalPrefix <- focal
      print(paste(focal,complexity.plant,complexity.animal,year))
      
      #---- Extract SpData specific for the focal and year and complexity----
      #rm(list = ls())
      abundance.plant <- read.csv(paste0(home.dic,"data/abundance.csv"), sep=",")
      competition.plant <- read.csv(paste0(home.dic,"data/competition.csv"), sep=",")
      plant.class <- read.csv(paste0(home.dic,"data/plant_code.csv"), sep=",")
      SpData <- read.csv(paste0(home.dic,"data/SpData.csv")) 
      complexitylevel.plant <- levels(as.factor(plant.class[,complexity.plant]))
      if (complexity.plant != "code.plant"){complexity.minimize <- T}else{complexity.minimize <- F}
      
      FocalPrefix <- focal # "LEMA" or "HOMA" or "CHFU"
      
      SpData <-  SpData %>%
        dplyr::select(all_of(c("day","month","year","plot",
                               "subplot","focal","fruit","seed",
                               "abundance",complexity.plant))) %>%
        aggregate(abundance ~ day + month + year + plot + focal+
                    subplot + .[,complexity.plant] +
                    fruit + seed , FUN=sum) %>%
        rename(complexity.plant = ".[, complexity.plant]") %>%
        tidyr::spread(complexity.plant ,abundance) %>%
        left_join(SpData %>% dplyr::select(all_of(c("day","month","year","plot",
                                                    "subplot","focal","fruit","seed",
                                                    "abundance","code.plant"))) %>%
                    dplyr::filter(code.plant ==  FocalPrefix) %>% tidyr::spread(code.plant,abundance) 
        ) %>%
        as.data.frame() %>%
        mutate(seed = round(seed))
      SpData <- SpData[which(SpData$year == as.numeric(year)),]
      
      if(as.numeric(levels(as.factor(SpData$year))) != as.numeric(year)){
        print("WRONG YEAR")}
      SpData <- dplyr::select(SpData,all_of(c("day","month", "year","plot","subplot" ,"focal",
                                              "fruit","seed",complexitylevel.plant,FocalPrefix)))
      
      SpDataFocal <- SpData[which(SpData$focal == FocalPrefix),]
      if(levels(as.factor(SpDataFocal$focal)) != FocalPrefix){
        print("WRONG FOCAL")}
      
      
      if (complexity.minimize == T){
        levels.of.focal <- plant.class[which(plant.class$code.plant==FocalPrefix),complexity.plant]
        SpDataFocal[,levels.of.focal] <- SpDataFocal[,levels.of.focal] - SpDataFocal[,FocalPrefix]
      }else{
        complexitylevel.plant <- complexitylevel.plant[which(!complexitylevel.plant %in% focal )]
      }
      
      #str(SpDataFocal)
      
      #---- Import residual----
      
      load(paste0(project.dic,"results/stan/FinalFit",FocalPrefix,"_",year,"_",complexity.plant,"_",complexity.animal,".rds")) 
      
      FinalPosteriors <- rstan::extract(FinalFit)
      
      mu <- FinalPosteriors$F_hat # matrix with nrow = draws and ncol = observations
      disp_dev <-FinalPosteriors$disp_dev
      phi <- (disp_dev^2)^(-1)
      
      # generating posterior predictions
      seed_pred <- matrix(nrow = dim(mu)[1], ncol = dim(mu)[2])
      for (i in 1:dim(mu)[1]) {     # for each posterior draw
        for (j in 1:dim(mu)[2]) {    # for each observation 
          # draw from the predicted distribution
          seed_pred[i, j] <- rnbinom(1, mu = mu[i, j], size = phi[i])  
          if(!is.na(seed_pred[i, j]) & seed_pred[i, j] > 100000){seed_pred[i, j]  <- 100000}
        }
      }
      
      seed.res.df <- seed_pred %>%
        as.data.frame() %>%
        summarise_all(~ median(.x, na.rm = TRUE)) %>%
        gather(key="obs",value="fitted.fec") %>%
        bind_cols(SpDataFocal %>%
                    dplyr::select(all_of(c("day","year","plot","subplot","focal","seed")))) %>%
        mutate(ID= fitted.fec - seed)
      
      dummy.df <- data.frame(subplot = rep(as.vector(outer(c("A","B","C","D","E","F"),
                                                           c(1:6), paste, sep="")),times=1),
                             X = rep(c(0,1.5,3,4.5,6,7.5),times=6),
                             Y = rep(c(0,1.5,3,4.5,6,7.5),each=6))
      for(n.plot in as.numeric(levels(as.factor( seed.res.df$plot)))){
        seed.res.df.p <- seed.res.df %>%
          filter(plot ==n.plot) %>%
          left_join(  dummy.df) %>%
          aggregate(ID~ X + Y,mean)
        if(nrow(seed.res.df.p)<5) next
        coords <- cbind(seed.res.df.p$X,seed.res.df.p$Y)
        #established who our neighbors are by creating an nb object
        nb <- dnearneigh(coords,0,7.5)
        #assign weights to each neighbor relationship.
        lw <- nb2listw( nb,style="W") # style = "W" indicates that the weights for each spatial unit are standardized to sum to 1 (this is known as row standardization).
        #moran.plot(as.numeric(scale(seed.res.df.p$ID)), listw=lw,
        #           xlab="Standardized Fecundity", 
        #       ylab="Standardized Lagged Fecundity")
        I <- moran(seed.res.df.p$ID, lw, length(nb), Szero(lw))[1]
        nsim= 999
        if(length(nb) <- 10){nsim=99}
        if(I==0) next
        if(I < 0 ){
          MC <- moran.mc(seed.res.df.p$ID, lw, nsim=  nsim, alternative="less")
          # View results (including p-value)
        }
        if(I >0 ){
          MC <- moran.mc(seed.res.df.p$ID, lw, nsim=  nsim, alternative="greater")
          # View results (including p-value)
        }
        
        moran.df.n <- data.frame(Moran.I = I$I,
                                 statistic = MC$statistic,
                                 observedrank = MC$parameter,
                                 alternativeH = MC$alternative,
                                 npermutation = nsim,
                                 pvalue = MC$p.value,
                                 focal=focal,
                                 plot=n.plot,
                                 year=as.numeric(year),
                                 complexity=complexity.plant,
                                 nobs=nrow(seed.res.df.p))
        #https://mgimond.github.io/simple_moransI_example/
        #a positive Moran’s I value suggests clustering whereas a negative Moran’s I value suggests dispersion
        moran.df <- bind_rows(moran.df,moran.df.n)
        write.csv(  moran.df,
                    file=paste0(home.dic, "results/Moran.df.csv"))
      }
    }
  }
}
view(moran.df)
moran.df <- moran.df %>%
  mutate(significance = case_when(pvalue<0.05 ~ "*",
                                  pvalue<0.01 ~ "**",
                                  pvalue<0.001 ~ "***",
                                  T~""))

write.csv(  moran.df,
            file=paste0(home.dic, "results/Moran.df.csv"))
moran.df <- read.csv(paste0(home.dic, "results/Moran.df.csv"))

moran.sum <- moran.df  %>% 
  group_by(complexity,focal) %>%
  summarise(n_positive = sum(statistic>0),
            n_negative= sum(statistic<0),
            n_significant0.01 = sum(pvalue <=0.01),
            n_significant0.05 = sum(pvalue <=0.05), 
            n= length(pvalue)) %>%
  mutate(percentage.sign.0.01 =  n_significant0.01/n,
         percentage.sign.0.05 =  n_significant0.05/n)
moran.sum
write.csv(  moran.sum,
            file=paste0(home.dic, "results/Moran.sum.csv"))
moran.family.df <-  moran.df %>%
  dplyr::filter(complexity =="family")
write.csv(  moran.family.df ,
            file=paste0(home.dic, "results/Moran.family.csv"))


nrow(moran.family.df[which(moran.family.df$pvalue<=0.01),])/nrow(moran.family.df)
nrow(moran.family.df[which(moran.family.df$pvalue<=0.05),])/nrow(moran.family.df)

write.csv(  moran.family.df,
            file=paste0(home.dic, "results/moran.family.df.csv"))


moran.class.df <-  moran.df %>%
  dplyr::filter(complexity =="class")
nrow(moran.class.df[which(moran.class.df$pvalue<=0.01),])/nrow(moran.class.df )
nrow(moran.class.df[which(moran.class.df$pvalue<=0.05),])/nrow(moran.class.df )


moran.species.df <-  moran.df %>%
  dplyr::filter(complexity =="code.plant")
nrow(moran.species.df[which(moran.species.df$pvalue<=0.01),])/nrow(moran.species.df)
nrow(moran.species.df[which(moran.species.df$pvalue<=0.05),])/nrow(moran.species.df)


# useful link for Moran I 
#https://cran.r-project.org/web/packages/spdep/vignettes/nb.html#:~:text=The%20dnearneigh()%20function%20is,argument%20to%20handle%20geographical%20coordinates.
