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

      specific.int.names <- c(colnames(Inclusion_all$Inclusion_FV)[Inclusion_all$Inclusion_FV==1],
              colnames(Inclusion_all$Inclusion_H)[Inclusion_all$Inclusion_H==1],
                    colnames(Inclusion_all$Inclusion_ij)[Inclusion_all$Inclusion_ij==1])
      
      abund.mean.df.i <- data.frame(focal = focal, year= year, 
                                    complexity.animal = complexity.animal,
                                    complexity.plant =  complexity.plant,
                                    name.sp = c(Inclusion_all$SpNames,
                                                 colnames(Inclusion_all$Inclusion_H),
                                                 colnames(Inclusion_all$Inclusion_FV)),
                                      abund.sp = c(Inclusion_all$SpMatrix %>% colMeans(),
                                                   Inclusion_all$SpMatrix_H %>% colMeans(),
                                                   Inclusion_all$SpMatrix_FV %>% colMeans()),
                                    abund.sp.min = c(Inclusion_all$SpMatrix %>% apply(.,2,min),
                                                 Inclusion_all$SpMatrix_H %>% apply(.,2,min),
                                                 Inclusion_all$SpMatrix_FV %>%  apply(.,2,min)),
                                    abund.sp.max = c(Inclusion_all$SpMatrix %>%  apply(.,2,max),
                                                 Inclusion_all$SpMatrix_H %>%  apply(.,2,max),
                                                 Inclusion_all$SpMatrix_FV %>%  apply(.,2,max)),
                                      parameter = c(rep("Plant - plant",length(Inclusion_all$SpNames)-1),
                                                      "Intraspecific",
                                                      rep("Plant - herbivore",length(colnames(Inclusion_all$Inclusion_H))),
                                                      rep("Plant - floral visitor",length(colnames(Inclusion_all$Inclusion_FV))))
                                      ) %>%
        mutate(inclus = case_when( name.sp %in% specific.int.names ~ "inclus",
                                   T ~ "non-inclus"))
               
      
      Trophic.mean.df.i <- data.frame(focal = focal, year= year, 
                                    complexity.animal = complexity.animal,
                                    complexity.plant =  complexity.plant,
                                    abund.sp.mean = c(mean(rowSums(Inclusion_all$SpMatrix[,!Inclusion_all$Intra])),
                                                 Inclusion_all$SpMatrix[,Inclusion_all$Intra] %>% 
                                                   mean(),
                                                 mean(rowSums(Inclusion_all$SpMatrix_H)),
                                                 mean(rowSums(Inclusion_all$SpMatrix_FV))),
                                    abund.sp.median = c(median(rowSums(Inclusion_all$SpMatrix[,!Inclusion_all$Intra])),
                                                 Inclusion_all$SpMatrix[,Inclusion_all$Intra] %>% 
                                                   median(),
                                                 median(rowSums(Inclusion_all$SpMatrix_H)),
                                                 median(rowSums(Inclusion_all$SpMatrix_FV))),
                                    abund.sp.max= c(max(Inclusion_all$SpMatrix[,!Inclusion_all$Intra]),
                                                        max(Inclusion_all$SpMatrix[,Inclusion_all$Intra]),
                                                        max(Inclusion_all$SpMatrix_H),
                                                        max(Inclusion_all$SpMatrix_FV)),
                                    abund.sp.min= c(min(Inclusion_all$SpMatrix[,!Inclusion_all$Intra]),
                                                    min(Inclusion_all$SpMatrix[,Inclusion_all$Intra]),
                                                    min(Inclusion_all$SpMatrix_H),
                                                    min(Inclusion_all$SpMatrix_FV)),
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
          "results/Trophic.mean.df.csv")
view(abund.mean.df)
write.csv(abund.mean.df,
          "results/abund.mean.df.csv")

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


