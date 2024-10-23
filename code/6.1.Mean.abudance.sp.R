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
home.dic <- "/home/lbuche/Eco_Bayesian/Complexity_caracoles/"

df_param_all <- read.csv(paste0(home.dic,
                                "results/Chapt1_Parameters_values.csv"))

##
abund.mean.df <- NULL
Trophic.mean.df <- NULL
for( focal.n in c("CETE","LEMA","HOMA","CHFU")){ # "CHFU","HOMA","CETE"
  for( year.n in c("2019",'2020','2021')){
    for( complexity.level in c(1:3)){
      print(paste0(focal.n,year.n,complexity.level))
      
      complexity.animal.n <- c("group","family","species")[complexity.level]
      complexity.plant.n <-c("class","family","code.plant")[complexity.level]
      load(paste0(home.dic,"results/Inclusion",
                  focal.n,"_", year.n,"_",complexity.plant.n,"_",
                  complexity.animal.n,".RData")
      )
      
      df_param_name_hat <- df_param_all %>%
        dplyr::filter(.[,"focal"] == as.character(focal.n)  & 
                        .[,"year"] == as.character(year.n)& 
                        .[,"complexity.plant"] == as.character(complexity.plant.n) &
                        !is.na(estimate_hat )) %>%
        dplyr::select(parameter,parameter_hat,complexity.plant) %>%
        unique()
      
      specific.int.names <-  df_param_name_hat$parameter_hat
      
      
      abund.mean.df.i <- data.frame(focal = focal.n, year= year.n, 
                                    complexity.animal = complexity.animal.n,
                                    complexity.plant =  complexity.plant.n,
                                    parameter_hat = c(Inclusion_all$SpNames,
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
        mutate(inclus = case_when( parameter_hat %in% specific.int.names ~ "inclus",
                                   T ~ "non-inclus"))
      
      if(sum(!specific.int.names %in% levels(as.factor(abund.mean.df.i$parameter_hat)))>0){
        missing_HOIs <- specific.int.names[!specific.int.names %in% levels(as.factor(abund.mean.df.i$parameter_hat))]
        param_hois <- df_param_name_hat$parameter[which(df_param_name_hat$parameter_hat==missing_HOIs)]
        missing_HOIs1 <- strsplit( missing_HOIs, '_')[[1]][1]
        missing_HOIs2 <- strsplit( missing_HOIs, '_')[[1]][2]
        
        abund.mean.df.i.hois <- data.frame(focal = focal.n, year= year.n, 
                                           complexity.animal = complexity.animal.n,
                                           complexity.plant =  complexity.plant.n,
                                           parameter_hat = missing_HOIs,
                                           inclus="inclus",
                                           parameter = param_hois,
                                           abund.sp = abund.mean.df.i$abund.sp[which(abund.mean.df.i$parameter_hat ==missing_HOIs1)]*
                                             abund.mean.df.i$abund.sp[which(abund.mean.df.i$parameter_hat ==missing_HOIs2)],
                                           abund.sp.max = abund.mean.df.i$abund.sp.max[which(abund.mean.df.i$parameter_hat ==missing_HOIs1)]*
                                             abund.mean.df.i$abund.sp.max[which( abund.mean.df.i$parameter_hat ==missing_HOIs2)],
                                           abund.sp.min = abund.mean.df.i$abund.sp.min[which(abund.mean.df.i$parameter_hat ==missing_HOIs1)]*
                                             abund.mean.df.i$abund.sp.min[which(abund.mean.df.i$parameter_hat ==missing_HOIs2)])
        abund.mean.df.i <- bind_rows( abund.mean.df.i, abund.mean.df.i.hois)
      }
      
      Trophic.mean.df.i <- data.frame(focal = focal.n, year= year.n, 
                                      complexity.animal = complexity.animal.n,
                                      complexity.plant =  complexity.plant.n,
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




