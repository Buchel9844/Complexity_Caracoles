################################################################################
# 00 - Import package
#################################################################################
#rm(list = ls(all = TRUE))
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(dplyr)) {install.packages("dplyr"); library(dplyr)}
if(!require(ggthemes)) {install.packages("ggthemes"); library(ggthemes)}
if(!require(vegan)) {install.packages("vegan"); library(vegan)}
if(!require(viridis)) {install.packages("viridis"); library(viridis)}
if(!require(corrplot)) {install.packages("corrplot"); library(corrplot)}
if(!require(RColorBrewer)) {install.packages("RColorBrewer"); library(RColorBrewer)}
if(!require(scales)) {install.packages("scales"); library(scales)}
if(!require(ggsci)) {install.packages("ggsci"); library(ggsci)}

setwd("~/Documents/Code/Project/HOI_Caracoles")
################################################################################
# 01 - Prepare the data
#################################################################################
source()
################################################################################
# 02 - Run the different models 
#################################################################################
source()
# lets figure out who the focal competitors are
id.plant.focal <- levels(as.factor(competition.plant$id.plant))
# lets figure out who the potential competitors are
id.plant.comp <- names(competition.plant)[which(!names(competition.plant) %in% c("year","plot","subplot","id.plant","fruit", "seed"))]
id.int.type <- names(HigherTL.wide.int.type)[which(!names(HigherTL.wide.int.type) %in% c("year","plot","subplot","id.plant"))]

for( plant.spc in id.plant.focal){
  model.of.plant.spc <-  list()
  # A ---- "Model.Null" -----
  for (space in c(TRUE,FALSE)){ for(time in c(TRUE,FALSE)){
  #fit a negative-binomial model without any effect of competition 
  model_of_[[paste(plant.spc,"Model.no.HOIs",sep="_")]] <- fit.fecundity.model(subset(Plant.and.id.inter,id.plant==plant.spc) ,
                                                                               type="negbin",
                                                                               fit.plot= FALSE,
                                                                               fit.year = FALSE,
                                                                               fit.alpha=FALSE, 
                                                                               fit.gama=FALSE,
                                                                               fit.betas.plant.plant = FALSE, 
                                                                               fit.betas.plant.HTL = FALSE, 
                                                                               fit.betas.HTL.HTL = FALSE)
    }   
  }
}

data <- subset(Plant.and.id.inter,id.plant=="BEMA")


################################################################################
# 02 - Run the different models 
#################################################################################

