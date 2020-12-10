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





