# this script make ...
# packages
library(ggplot2)
library(tidyverse)
library(ggsci)
library(tidyverse)
#install.packages("ggthemes")
library(ggthemes)
#install.packages("ggridges")
library(ggridges)
#install.packages("wesanderson")
library(wesanderson)
library(ggpubr)
library(cowplot)
library(magick)
library(png)
library(grid)

# define you directory
home.dic <- "/Users/lisabuche/Documents/Projects/Complexity_caracoles/"


##
abund.mean.df <- NULL
Trophic.mean.df <- NULL
for( focal in c("CETE","LEMA","HOMA","CHFU")){ # "CHFU","HOMA","CETE"
  for( year in c("2019",'2020','2021')){
    for( complexity.level in c(1:3)){
      
      complexity.animal <- c("group","family","species")[complexity.level]
      complexity.plant <-c("class","family","code.plant")[complexity.level]
      load(paste0(home.dic,"results/stan/Inclusion",
                  focal,"_", year,"_",complexity.plant,"_",
                  complexity.animal,".RData")
           )

      abund.mean.df.i <- data.frame(focal = focal, year= year, 
                                    complexity.animal = complexity.animal,
                                    complexity.plant =  complexity.plant,
                                    name.sp = c(Inclusion_all$SpNames,
                                                 colnames(Inclusion_all$Inclusion_H),
                                                 colnames(Inclusion_all$Inclusion_FV)),
                                      abund.sp = c(Inclusion_all$SpMatrix %>% colMeans(),
                                                   Inclusion_all$SpMatrix_H %>% colMeans(),
                                                   Inclusion_all$SpMatrix_FV %>% colMeans()),
                                      parameter = c(rep("Plant - plant",length(Inclusion_all$SpNames)-1),
                                                      "Intraspecific",
                                                      rep("Plant - herbivore",length(colnames(Inclusion_all$Inclusion_H))),
                                                      rep("Plant - floral visitor",length(colnames(Inclusion_all$Inclusion_FV))))
                                      )
      Trophic.mean.df.i <- data.frame(focal = focal, year= year, 
                                    complexity.animal = complexity.animal,
                                    complexity.plant =  complexity.plant,
                                    abund.sp.mean = c(mean(rowSums(Inclusion_all$SpMatrix[,!Inclusion_all$Intra])),
                                                 Inclusion_all$SpMatrix[Inclusion_all$Intra] %>% 
                                                   mean(),
                                                 mean(rowSums(Inclusion_all$SpMatrix_H)),
                                                 mean(rowSums(Inclusion_all$SpMatrix_FV))),
                                    abund.sp.median = c(median(rowSums(Inclusion_all$SpMatrix[,!Inclusion_all$Intra])),
                                                 Inclusion_all$SpMatrix[Inclusion_all$Intra] %>% 
                                                   median(),
                                                 median(rowSums(Inclusion_all$SpMatrix_H)),
                                                 median(rowSums(Inclusion_all$SpMatrix_FV))),
                                    parameter = c("Plant - plant",
                                                    "Intraspecific",
                                                    "Plant - herbivore",
                                                    "Plant - floral visitor")
                                    )
      
      abund.mean.df <- bind_rows(   abund.mean.df,   abund.mean.df.i)
      Trophic.mean.df <- bind_rows( Trophic.mean.df,   Trophic.mean.df.i)
      
      remove(Inclusion_all)
    }
  }
}
view(Trophic.mean.df)
write.csv(Trophic.mean.df,
          "results/Trophic.mean.df")

# for specific observation

summary.interactions.short <- read.csv(paste0(home.dic,"results/inclusion/Chapt1_Inclusion_parameters_short.csv"))

Hat.mean.df <- NULL
for( n in 1:nrow( summary.interactions.short)){
      df_hat <- summary.interactions.short[n,]
      
          load(paste0(home.dic,"results/stan/Inclusion",
                  df_hat$focal,"_", df_hat$year,"_",df_hat$complexity_plant,"_",
                  df_hat$complexity.animal,".RData"))
          mean.plant_hat <-NA
          mean.FV_hat <-NA
          mean.H_hat <-NA
   if(df_hat$n.competitors_plant_inclus==1){
    mean.plant_hat <- mean(Inclusion_all$SpMatrix[,df_hat$competitors_plant_inclus == Inclusion_all$SpNames])
    competitors_plant_inclus <- df_hat$competitors_plant_inclus
      }
  if(df_hat$n.competitors_plant_inclus>1){
            mean.plant_hat <- colMeans(Inclusion_all$SpMatrix[,Inclusion_all$SpNames %in% strsplit(df_hat$competitors_plant_inclus,split=",")[[1]]])
            competitors_plant_inclus <- strsplit(df_hat$competitors_plant_inclus,split=",")[[1]]
            }
          
   if(df_hat$n.competitors_FV_inclus==1){
      mean.FV_hat <-mean(Inclusion_all$SpMatrix_FV[,df_hat$competitors_FV_inclus == colnames(Inclusion_all$Inclusion_FV)])
      }
    if(df_hat$n.competitors_H_inclus==1){
      mean.H_hat <-mean(Inclusion_all$SpMatrix_H[,df_hat$competitors_H_inclus == colnames(Inclusion_all$Inclusion_H)])
      }
          
          
      Hat.mean.df.i <- data.frame(focal = df_hat$focal, year= df_hat$year, 
                                      complexity.animal = df_hat$complexity.animal,
                                      complexity.plant =  df_hat$complexity_plant,
                                  
                                      abund.sp.mean = c(mean.plant_hat,
                                                        mean.FV_hat,
                                                        mean.H_hat),
                                      parameter_hat =c(competitors_plant_inclus,
                                                       df_hat$competitors_FV_inclus,
                                                       df_hat$competitors_H_inclus ),
                                      parameter = c(rep("Plant - plant",length(competitors_plant_inclus)),
                                                    "Plant - floral visitor",
                                                    "Plant - herbivore"))
      
      Hat.mean.df <- bind_rows(  Hat.mean.df,  Hat.mean.df.i)
      
}
Hat.mean.df <-  Hat.mean.df %>%
  filter(!is.na(abund.sp.mean))
view(Hat.mean.df)


# unstandard abundance
plant.class <- read.csv(paste0(home.dic,"data/plant_code.csv"), sep=",")
complexitylevel.plant <- levels(as.factor(plant.class[,complexity.plant]))

SpData <- read.csv(paste0(home.dic,"data/SpData.csv"))
SpDataFocal <- NULL
for( n in c("CHFU","LEMA","CETE","HOMA")){
SpDataFocal.n <- SpData[which(SpData$focal ==n),] %>%
    mutate(abundance.all = (rowSums(.[c("Forb","Grass")]) - rowSums(.[c(n)]))) %>%
    aggregate(abundance.all ~ focal + year, median)
SpDataFocal <- bind_rows(SpDataFocal,SpDataFocal.n)
       
}  

SpDataFocal <- SpDataFocal %>%
  filter(year %in% c("2019","2020","2021"))

max(SpDataFocal$abundance.all)

median(SpDataFocal$abundance.all)


