#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 0. METADATE: Create meta data file to Summarise our data----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
install.packages("readr")
install.packages("devtools")
devtools::install_github("ropenscilabs/dataspice")

library(dataspice)

# 0.1. CREATE TEMPLATES
create_spice(paste0(home.dic,"data/rawdata/")) # creates template from directory

# 0.2. EDIT TEMPLATES
#attributes.csv: This is where most of the user data entry will take place. 
# For each variable, its name, units, and a written description are filled in.
edit_attributes()

# access.csv: Includes a row for each file that was read in, 
# and documents the name of each file and its format.
edit_access()

# creators.csv: One row for each creator, and gives their affiliation,
# contact email, ORCID, etc
edit_creators() #

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
view(full.comp)

full.comp <- full.comp %>%
  dplyr::select(-c(comments, observers)) %>%
  mutate(seed = NA) 


full.comp.plot <- full.comp %>%
  filter(focal %in% c("HOMA","CHFU","LEMA","CETE")) %>%
  filter(year %in% c("2019","2020","2021")) %>%
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
  labs(color="year",x="seeds per fruit") + 
  scale_color_colorblind()+
  facet_wrap(focal~., scales="free") + 
  theme_bw() 
seed.prod.plot
ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/seed.prod.plot.pdf"),
       dpi="retina",
       width = 21,
       height = 16,
       units = c("cm"),
       seed.prod.plot)

#---- 2.3. Join seed and interactions -----

focal.comp <- full.comp %>%
  filter(focal %in% c("HOMA","CHFU","LEMA","CETE")) %>%
  filter(year %in% c("2019","2020","2021"))
view(focal.comp)
set.seed(1616)

for(i in c(1:nrow(focal.comp))){
  year.i <- focal.comp[i,"year"] 
  focal.i <- focal.comp[i,"focal"]
  plot.i <- focal.comp[i,"plot"] 
  subplot.i <- focal.comp[i,"subplot"] 
  seeds.fruit.i <- NA
  
  if(focal.i =="CETE" & year.i  %in% c("2019","2021")){
    year.i = "2020"} # seed sample in 2020 
  if(focal.i == c("HOMA") & year.i %in% c("2019","2021")){
    year.i = c("2020") # seed sample in 2020 and 2019 
  }
  if(focal.i  %in% c("CHFU","LEMA") & year.i  %in% c("2019")){
    year.i = c("2020","2019","2021") # seed sample in all years
  }
  
  seeds.fruit.i <- abs(rnorm(1,mean = mean(seed.prod[which(seed.prod$year %in% year.i &
                                                             seed.prod$focal == focal.i &
                                                             seed.prod$plot == plot.i &
                                                             seed.prod$subplot == subplot.i),
                                                     "seeds.fruit"], na.rm=T),
                             sd = sd(seed.prod[which(seed.prod$year  %in% year.i &
                                                       seed.prod$focal == focal.i &
                                                       seed.prod$plot == plot.i),
                                               "seeds.fruit"], na.rm=T)))
  
  if(is.na(seeds.fruit.i)){
    seeds.fruit.i <- abs(rnorm(1,mean(seed.prod[which(seed.prod$year  %in% year.i &
                                                        seed.prod$focal == focal.i &
                                                        seed.prod$plot == plot.i),
                                                "seeds.fruit"], na.rm=T),
                               sd = sd(seed.prod[which(seed.prod$year  %in% year.i &
                                                         seed.prod$focal == focal.i),
                                                 "seeds.fruit"], na.rm=T)))
    
    
  }
  if(is.na(seeds.fruit.i)){
    seeds.fruit.i <- abs(rnorm(1,mean(seed.prod[which(seed.prod$year %in% year.i &
                                                        seed.prod$focal == focal.i),
                                                "seeds.fruit"], na.rm=T),
                               sd = sd(seed.prod[which(seed.prod$year %in% year.i),
                                                 "seeds.fruit"], na.rm=T)))
    
  }
  
  focal.comp$seed[i] <-  seeds.fruit.i*focal.comp[i,"fruit"]
}

view(full.comp)
focal.comp.seed.plot <- focal.comp %>%
  filter(focal %in% c("HOMA","CHFU","LEMA","CETE")) %>%
  filter(year %in% c("2019","2020","2021")) %>%
  filter(seed < 3000) %>%
  ggplot(aes(seed,group=as.factor(year),color=as.factor(year))) +
  geom_density() +
  labs(color="year",x="total seeds per individual") + 
  scale_color_colorblind() +
  facet_wrap(focal~., scales="free") + 
  theme_bw()  
focal.comp.seed.plot 


ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/full.comp.seed.plot.pdf"),
       dpi="retina",
       width = 21,
       height = 16,
       units = c("cm"),
       focal.comp.seed.plot)


#check the present of Na in seed production
focal.comp.NA <- focal.comp[is.na(focal.comp$seed),]
levels(as.factor(cbind(focal.comp.NA$year,focal.comp.NA$focal)))

#threshold of seed. production to 1000

focal.comp <- focal.comp[which(focal.comp$seed <= 3000),]

#capitalize the plant species names
not.capitalized.name <- c("day","month","year", "plot","subplot","focal","fruit",
                          "seed","observers","comments")
names(focal.comp)[!names(focal.comp) %in% not.capitalized.name] <-toupper(names(focal.comp)[!names(focal.comp) %in% not.capitalized.name])

view(focal.comp)
# 3 - write the csv for all years and per years
readr::write_delim(focal.comp ,file = paste0(home.dic,"data/competition.csv"),
                   delim = ",")
#----Make SpDATA----

competition.plant <- read.csv(paste0(home.dic,"data/competition.csv"), sep=",")

species.neigh <- names(competition.plant )
species.neigh  <- species.neigh[!species.neigh %in% c("day","month","year",
                                                      "plot","subplot","focal","fruit","comments",
                                                      "observers", "site",
                                                      "treatment","seed")]

competition.to.keep <- competition.plant  %>%  
  dplyr::select(all_of(c(species.neigh)))%>%
  mutate_at( species.neigh, as.numeric) %>%
  colSums(na.rm=T)


# change the format of the rows to numeric 
competition.plant[species.neigh] <- sapply(competition.plant[species.neigh],as.numeric)

# change na values to 0
competition.plant[is.na(competition.plant)] <- 0

plant.class <- read.csv(paste0(home.dic,"data/plant_code.csv"),sep=",")


competition.plant_long <- competition.plant %>%
  gather(any_of(plant.class$code.plant), key="code.plant",
         value="abundance") %>%
  #left_join(plant.class, by="code.plant") %>%
  left_join(plant.class %>%
              select(code.plant,family,class),
            by="code.plant")  %>%
  dplyr::select(day,month,year,plot,subplot,focal,fruit,seed,code.plant,family,class,
                abundance) %>%
  aggregate(abundance ~ day + month + year + plot + focal+
              subplot + code.plant + family + class +
              fruit + seed , FUN=sum) 

readr::write_delim(competition.plant_long ,
                   file = paste0(home.dic,"data/SpData.csv"),
                   delim = ",")

