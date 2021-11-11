plot.coords <- read.csv("Caracoles/UnderTheHood/raw_data/abiotic_data/Caracoles_allplotsposition.csv", sep=";")
str(plot.coords)

#---- Import the Data ----
abundance.plant <- read.csv("/Users/lisabuche/Code/Project/Caracoles/data/abundance.csv", sep=",")

competition.plant <- read.csv("/Users/lisabuche/Code/Project/Caracoles/data/competition.csv", sep=",")

plant.class <- read.csv("/Users/lisabuche/Code/Project/Caracoles/data/plant_code.csv", sep=",")

#################################################################################
# Plant abundance to compute rare or abundant, done! 
#################################################################################


#---- visualise abundance ----
abundance.plant.gathered <- aggregate(individuals ~ year + month + day + species,
                                            abundance.plant,mean)
  library(tidyr)
abundance.plant.gathered <- unite(abundance.plant.gathered,year,month,day,
                                        col="date",sep = "-", remove = F)


abundance.plant.gathered$date <- as.Date(abundance.plant.gathered$date)
library(ggplot2)
ggplot(abundance.plant.gathered,aes(x=date,y=individuals,color=species)) +
  geom_line() + 
  xlab("") +
  scale_x_date(date_labels = "%Y %b") +
  scale_y_sqrt("mean abundance per subplot per day ") + 
  theme_bw()

#---- Determine the rare species througt the years -----
library(dplyr)
abundance.sp.years <- full_join(aggregate(individuals ~ year + species, abundance.plant,sum),
       aggregate(individuals ~ year,abundance.plant,sum),
       by = c("year"),
       suffix = c(".sp.year", ".year"))
 head(abundance.sp.years)

 abundance.sp.years$ratio <- abundance.sp.years$individuals.sp.year/abundance.sp.years$individuals.year
 
 abundance.sp.all.years <- aggregate(ratio ~  species, abundance.sp.years,sum)
 
 abundance.sp.years$rare <- F
 abundance.sp.years$rare[which(abundance.sp.years$ratio < 0.01)] <- T


 abundance.rare <- aggregate(rare ~ species,abundance.sp.years,mean)
 
 abundance.rare$rare.overall <- "abundant"
 # give the species which were rare (>1% of the total abundance) over more than one year 
 abundance.rare$rare.overall[which(abundance.rare$rare >= 0.4)] <- "rare"
 
 abundance.rare <- abundance.rare[,-2]
 names(abundance.rare) <-c("code.plant","rareORabundant")

 
 #---- Join dataset with plant.class ----
 
 
plant.class <- right_join( plant.class, abundance.rare,by = c("code.plant")) 
head( plant.class)
View( plant.class)

#################################################################################
# Plant competiton
#################################################################################

#---- Join dataset with competition----

code.plant.comp <- names(competition.plant)[!names(competition.plant) %in%
                                              c("day","month", "year","plot","subplot" ,
                                                "focal","fruit","seed","comments","observers")]
competition.plant.test <- gather(competition.plant,all_of(code.plant.comp),
                                  key = "code.plant", value = "abundance")
head(competition.plant.test)
 
competition.plant.test <- inner_join(competition.plant.test, plant.class, by = c("code.plant"))
head(competition.plant.test)

# merge data frame horyzntally
names(competition.plant)
head(plant.class)
 for ( n in names(plant.class)[c(!names(plant.class) %in% c("number","code.plant","species","Native.sp"))]){
   competition.plant.test.test <- competition.plant.test %>%
     dplyr::select(year,plot,subplot,focal,fruit,abundance,seed,n)
   competition.plant.test.test$abundance <- as.numeric(competition.plant.test.test$abundance)
   competition.plant.test.test <- aggregate(abundance ~ competition.plant.test.test[,n] + year + plot
                                            + subplot+ focal+ fruit + seed,
                                           competition.plant.test.test,sum)
   names(competition.plant.test.test)[1] <- c(n) 
   
   competition.plant.test.test <-spread(competition.plant.test.test, key=n ,value='abundance')
   
   competition.plant <- inner_join(competition.plant,
                                   competition.plant.test.test, 
                                   by = c("year","plot","subplot","focal","fruit","seed"))
   }
 head(competition.plant)
 write.csv( competition.plant,"data/competition.csv")

