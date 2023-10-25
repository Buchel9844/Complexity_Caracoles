#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 0. METADATE: Create meta data file to Summarise our data----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
install.packages("readr")
install.packages("devtools")
library(devtools)
library(readr)
devtools::install_github("ropenscilabs/dataspice")

library(dataspice)

# 0.1. CREATE TEMPLATES
dataspice::create_spice() # creates template from directory

# 0.2. EDIT TEMPLATES
#attributes.csv: This is where most of the user data entry will take place. 
# For each variable, its name, units, and a written description are filled in.
prep_attributes()
edit_attributes()

# access.csv: Includes a row for each file that was read in, 
# and documents the name of each file and its format.
prep_access()
edit_access()

# creators.csv: One row for each creator, and gives their affiliation,
# contact email, ORCID, etc
dataspice::edit_creators() #

#biblio.csv: Citation information about the project, 
# as much or as little data as possible can be included, 
# but if things like bounding box coordinates are not included, 
# then when the website is generated there will not be a bounding box map generated
edit_biblio()



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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 2.IMPORT DATA -----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#rm(list = ls()) # remove environment to not use the objects created before

home.dic <- "/home/lbuche/Eco_Bayesian/Complexity_caracoles/"
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/"

#---- 2.1. Competiton -----
comp.files <- list.files(path = paste0(home.dic,"data/rawdata/"),pattern = "competition_")

full.comp <- NULL

for(i.comp in 1:length(comp.files)){
  my.file <- read.csv2(paste0(home.dic,"data/rawdata/",comp.files[i.comp]),sep=",",
                       header = T,stringsAsFactors = F)
  if (length(names(my.file))==1){
    my.file <- read.csv2(paste0(home.dic,"data/rawdata/",comp.files[i.comp]),sep=";",
                         header = T,stringsAsFactors = F)
  }
  #readr::read_delim(paste("raw_data/",abundance.files[i.abund],sep=""),delim = ";")
  
  # replace records with NA's in number of individuals by 0
  my.file[is.na(my.file)] <- 0
  my.file[my.file$plot==0] <- NA
  my.file <- my.file[complete.cases(my.file),]
  # convert names of columns to lower case
  names(my.file) <- tolower(names(my.file))
  
  # convert date into day|month|year
  if(length(my.file[names(my.file) %in% c("year","month","day")])==0){
    # in some cases this happens
    my.file <- my.file[which(my.file$date != "na"),]
    
    my.file$temp.dates <- as.Date(my.file$date,"%d/%m/%Y")
    my.file$year <- as.numeric(format(my.file$temp.dates, format = "%Y"))
    my.file$month <- as.numeric(format(my.file$temp.dates, format = "%m"))
    my.file$day <- as.numeric(format(my.file$temp.dates, format = "%d"))
  }
  
  # clean up
  

    full.comp <- bind_rows(full.comp,my.file)
  
}
names(full.comp)

full.comp <- full.comp %>%
  dplyr::select(-c(comments, observers)) %>%
  mutate(seed = NA) 


full.comp.plot <- full.comp %>%
  filter(focal %in% c("HOMA","CHFU","LEMA","CETE")) %>%
  aggregate(fruit ~ focal + year, length) %>%
  #filter(year %in% c("2020","2021")) %>%
  ggplot(aes(y=fruit,x=as.factor(year),fill=as.factor(focal))) +
  geom_bar(stat="identity", position="dodge") +
  theme_bw() +
  scale_fill_colorblind() +
  labs(y="number of observations", x="year",fill="focal") 
full.comp.plot

ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/ObservationsFocal.pdf"),
       dpi="retina",
       width = 21,
       height = 16,
       units = c("cm"),
       full.comp.plot)

#---- 2.2. Seed production -----


seed.files <- list.files(path = paste0(home.dic,"data/rawdata/"),pattern = "numb_seed_")

seed.prod <- NULL

for(i.seed in 1:length(seed.files)){
  my.file <- read.csv2(paste0(home.dic,"data/rawdata/",seed.files[i.seed]),sep=",",
                       header = T,stringsAsFactors = F)
  if (length(names(my.file))==1){
    my.file <- read.csv2(paste0(home.dic,"data/rawdata/",seed.files[i.seed]),sep=";",
                         header = T,stringsAsFactors = F)
  }
  #readr::read_delim(paste("raw_data/",abundance.files[i.abund],sep=""),delim = ";")
  
  # replace records with NA's in number of individuals by 0
  
  my.file$seeds.fruit[which(my.file$seeds.fruit ==0)] <-NA
  my.file <- my.file %>% drop_na(seeds.fruit)
  # convert names of columns to lower case
  names(my.file) <- tolower(names(my.file))
  
  # convert date into day|month|year
  if(length(my.file[names(my.file) %in% c("year","month","day")])==0){
    # in some cases this happens
    my.file <- my.file[which(my.file$date != "na"),]
    
    my.file$temp.dates <- as.Date(my.file$date,"%d/%m/%Y")
    my.file$year <- as.numeric(format(my.file$temp.dates, format = "%Y"))
    my.file$month <- as.numeric(format(my.file$temp.dates, format = "%m"))
    my.file$day <- as.numeric(format(my.file$temp.dates, format = "%d"))
  }
  
  # clean up
  
  
  seed.prod <- bind_rows(seed.prod,my.file)
  
}
view(seed.prod)


seed.prod.plot <- seed.prod %>%
  filter(focal %in% c("HOMA","CHFU","LEMA","CETE")) %>%
  #filter(year %in% c("2020","2021")) %>%
  ggplot(aes(seeds.fruit,group=as.factor(year),color=as.factor(year))) +
  geom_density() +
  facet_wrap(focal~., scales="free") + 
  theme_bw() 

ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/seed.prod.plot.pdf"),
       dpi="retina",
       width = 21,
       height = 16,
       units = c("cm"),
       seed.prod.plot)

#---- 2.3. Join seed and interactions -----
seed.prod <- seed.prod %>%
  filter(year %in% c("2020","2021"))

for(i in c(1:nrow(full.comp))){
  year.i <- full.comp[i,"year"] 
  focal.i <- full.comp[i,"focal"]
  seeds.fruit.i <- NA
  if (year.i == "2015"|year.i == "2016"){
    seeds.fruit.i <- seed.2015[which(seed.2015$year==year.i & seed.2015$focal == focal.i),
                               "seeds.fruit"]
    if(is.na(seeds.fruit.i)){seeds.fruit.i <- mean(seed.2015[which(seed.2015$focal == focal.i),
                                                             "seeds.fruit"], na.rm = TRUE)}
    full.comp$seed[i] <- seeds.fruit.i*full.comp[i,"fruit"]
  }
  if (year.i == "2017"|year.i == "2018"){
    seeds.fruit.i <- mean(seed.2015[which(seed.2015$focal == focal.i),
                                    "seeds.fruit"], na.rm = TRUE)
    full.comp$seed[i] <-  seeds.fruit.i*full.comp[i,"fruit"]
  }
  if (year.i == "2019"|year.i == "2020"|year.i == "2021"){
    plot.i <- full.comp[i,"plot"] 
    subplot.i <- full.comp[i,"subplot"] 
    
    seeds.fruit.i <- mean(seed.2019.2021[which(seed.2019.2021$year==year.i &
                                                 seed.2019.2021$focal == focal.i &
                                                 seed.2019.2021$plot == plot.i &
                                                 seed.2019.2021$subplot == subplot.i),
                                         "seeds.fruit"],na.rm=T)
    if(is.na(seeds.fruit.i)){
      seeds.fruit.i <- mean(seed.2019.2021[which(seed.2019.2021$year==year.i &
                                                   seed.2019.2021$focal == focal.i &
                                                   seed.2019.2021$plot == plot.i),
                                           "seeds.fruit"],na.rm=T)
      
    }
    if(is.na(seeds.fruit.i)){
      seeds.fruit.i <- mean(seed.2019.2021[which(seed.2019.2021$year==year.i &
                                                   seed.2019.2021$focal == focal.i),
                                           "seeds.fruit"],na.rm=T)
      
    }
    if(is.na(seeds.fruit.i)){
      seeds.fruit.i <- mean(c(seed.2019.2021[which(seed.2019.2021$focal == focal.i),
                                             "seeds.fruit"],
                              seed.2015[which(seed.2015$focal == focal.i),
                                        "seeds.fruit"]),na.rm=T)
      
    }
    
    full.comp$seed[i] <-  seeds.fruit.i*full.comp[i,"fruit"]
  }
  if (is.na(full.comp$seed[i])){
    seeds.fruit.i <- mean(c(seed.2019.2021[which(seed.2019.2021$focal == focal.i),
                                           "seeds.fruit"], 
                            seed.2015[which(seed.2015$focal == focal.i),
                                      "viable.seeds.fruit"]),na.rm=T)
    full.comp$seed[i] <-  seeds.fruit.i*full.comp[i,"fruit"]
  }
  
}
#check the present of Na in seed production
full.comp.NA <- full.comp[is.na(full.comp$seed),]
levels(as.factor(cbind(full.comp.NA$year,full.comp.NA$focal)))

#remove fruit ==0 
full.comp <- full.comp[-which( full.comp$fruit ==0 ),]

#threshold of seed. production to 1000

full.comp$seed[which( full.comp$seed > 1000 )] <- 1000

#capitalize the plant species names
not.capitalized.name <- c("day","month","year", "plot","subplot","focal","fruit",
                          "seed","observers","comments")
names(full.comp)[!names(full.comp) %in% not.capitalized.name] <-toupper(names(full.comp)[!names(full.comp) %in% not.capitalized.name])

view(full.comp)
# 3 - write the csv for all years and per years
readr::write_delim(full.comp ,file = "UnderTheHood/data/competition.csv",delim = ",")


readr::write_delim(full.comp ,file = "data/competition.csv",delim = ",")
competition <- full.comp 
save(competition ,file = "data/competition.rda")