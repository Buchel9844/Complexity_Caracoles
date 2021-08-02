#################################################################################
# Plant abundance and competition
#################################################################################
plot.coords <- read.csv("Caracoles/UnderTheHood/raw_data/abiotic_data/Caracoles_allplotsposition.csv", sep=";")
str(plot.coords)

#################################################################################
# Plant abundance and competition
#################################################################################
#---- Import the Data ----
abundance.plant <- read.csv("Caracoles/UnderTheHood/data/to_delete/abundances_ME_wide.csv", sep=";")
str(abundance.plant)
code.plant <- c(names(abundance.plant[,6:ncol(abundance.plant)]))
abundance.plant.gathered <- gather(abundance.plant, all_of(code.plant),
               key = "code.plant", value = "abundance")
head(abundance.plant.gathered)

competition.plant <- read.csv("Caracoles/UnderTheHood/data/to_delete/competition_ME_wide.csv", sep=";")
str(competition.plant)
names(competition.plant)[4] <- "id.plant"


#---- Classify the plant ----
plant.class <- read.csv("data/plant_class.csv", sep=",")

#---- visualise abundance ----
abundance.plant.gathered.short <- aggregate(abundance ~ year+month+day +code.plant,
                                            abundance.plant.gathered,mean)
  
abundance.plant.gathered.short <- unite(abundance.plant.gathered.short,year,month,day,
                                        col="date",sep = "-", remove = F)

head(abundance.plant.gathered)

abundance.plant.gathered.short$date <- as.Date(abundance.plant.gathered.short$date)

ggplot(abundance.plant.gathered.short,aes(x=date,y=abundance,color=code.plant)) +
  geom_line() + 
  xlab("") +
  scale_x_date(date_labels = "%Y %b")

#---- to determine the rare species througt the years -----
abundance.sp.years <- right_join(aggregate(abundance ~ year + code.plant, abundance.plant.gathered,sum),
       aggregate(abundance ~ year,abundance.plant.gathered,sum),
       by = c("year"),
       suffix = c(".sp.year", ".year"))
 head(abundance.sp.years)

 abundance.sp.years$ratio <- abundance.sp.years$abundance.sp.year/abundance.sp.years$abundance.year
 abundance.sp.years$ratio.all.years <- abundance.sp.years$abundance.year/ sum(abundance.plant.gathered$abundance,na.rm=T)
 
 abundance.sp.years$rare <- F
 abundance.sp.years$rare[which(abundance.sp.years$ratio < 0.01)] <- T
 abundance.sp.years$rare.all.years <- F
 abundance.sp.years$rare[which(abundance.sp.years$ratio.all.years < 0.01)] <- T
 
 
 abundance.rare <- aggregate(rare ~ code.plant,abundance.sp.years,mean)
 abundance.rare$rare.overall <- "abundant"
 # give the species which were rare (>1% of the total abundance) over more than one year 
 
 abundance.rare$rare.overall[which(abundance.rare$rare >= 0.4)] <- "rare"
 abundance.rare <- abundance.rare[,-2]
 
 #---- Join dataset ----
 code.plant.comp <- c(names(competition.plant[,7:ncol(competition.plant)]))
 
 plant.class <- right_join(plant.class,abundance.rare, by = c("code.plant"))

  competition.plant.test <- gather(competition.plant,all_of(code.plant.comp),
                                  key = "code.plant", value = "abundance")
 
 head(competition.plant.test)
 
 competition.plant.test <- inner_join(competition.plant.test,plant.class, by = c("code.plant"))
head(competition.plant.test.test)

# merge data frame horyzntally

 for ( n in names(plant.class)[c(!names(plant.class) %in% c("code.plant","species"))]){
   competition.plant.test.test <- competition.plant.test %>%
     dplyr::select(year,plot,subplot,id.plant,fruit,seed,abundance,n)
   competition.plant.test.test <-aggregate(abundance ~ competition.plant.test.test[,n] + year + plot + subplot+ id.plant+ fruit + seed,
                                           competition.plant.test.test,sum)
   
   names(competition.plant.test.test)[1] <- c(n) 
   
   competition.plant.test.test <-spread(competition.plant.test.test, key=n ,value='abundance')
   
   competition.plant <- inner_join(competition.plant,
                                   competition.plant.test.test, 
                                   by = c("year","plot","subplot","id.plant","fruit","seed"))
   }
 

 head(competition.plant)

#################################################################################
# HigherTL & floral visitor
#################################################################################
#---- Import the Data ----
herbivor2019 <- read.csv("Caracoles/UnderTheHood/raw_data/foodweb/foodweb_2019.csv", sep=";")
herbivor2020 <- read.csv("Caracoles/UnderTheHood/raw_data/foodweb/foodweb_2020.csv", sep=",")

herbivor2020$Day <- as.integer(herbivor2020$Day)
names(herbivor2020)[which(!names(herbivor2020) %in% names(herbivor2019))]
herbivor2019[,names(herbivor2020)[which(!names(herbivor2020) %in% names(herbivor2019))]] <- NA
herbivor <- bind_rows(herbivor2019,herbivor2020)
nrow(herbivor[which(is.na(herbivor$Plant)),]) # plant_simpla ?? Maria ! 
herbivor$Sampling  <- "herbivor"

# here remove pollinators from herbivore data.frame ? 

# Upload floral_visitor data set
floral_visitors2016_2019 <- read.csv("Caracoles/UnderTheHood/raw_data/floral_visitor/pollinator_2016_2019.csv", sep=";")
floral_visitors2016_2019 <- floral_visitors2016_2019[,-c(ncol(floral_visitors2016_2019):100)]
floral_visitors2020 <- read.csv("Caracoles/UnderTheHood/raw_data/floral_visitor/pollinator_2020.csv", sep=";")
floral_visitors2016_2019[,names(floral_visitors2020)[which(!names(floral_visitors2020) %in% names(floral_visitors2016_2019))]] <- NA
floral_visitors2020[,names(floral_visitors2016_2019)[which(!names(floral_visitors2016_2019) %in% names(floral_visitors2020))]] <- NA
floral_visitors2020$Plot <- as.character(floral_visitors2020$Plot)

floral_visitors <- bind_rows(floral_visitors2016_2019,floral_visitors2020)
floral_visitors$Sampling <- "floral_visitor"
#names(floral_visitors)[which(!names(floral_visitors) %in% names(herbivor))]

# How to combine both dataset ???
herbivor <- subset(herbivor,select=c("Year","Plot","Subplot","Sampling","Group" ,"Order","ID","Plant","X."))
floral_visitors <- subset(floral_visitors,select=c("Year","Plot","Subplot","Sampling","Group" ,"Order","ID_Simple","Plant","Visits"))

names(herbivor) = names(floral_visitors) <- c("year","plot","subplot","Sampling","group.focal" ,"order.focal","id.focal","id.plant","frequency")
floral_visitors$year <- as.integer(floral_visitors$year)
floral_visitors$plot <- as.integer(floral_visitors$plot)

# bind all HigherTL in the same dataframe
HigherTL <- bind_rows(herbivor,floral_visitors)

# compute the mean of frequency of each herbivor for a specific set of sp1-year-plot-subplot
HigherTL <- HigherTL %>%
  subset(id.focal != "soil") %>%
  group_by(Sampling,id.focal,id.plant,year, plot, subplot) %>%
  dplyr::summarize(Mean = mean(frequency, na.rm=TRUE),
                               group.focal= group.focal ,
                               order.focal=order.focal) %>%
  ungroup() %>%
  drop_na(id.focal,id.plant,year, plot)%>%

  filter(group.focal != "0")
##########################################
#---- construction of different group ----

# Make the group uniform
HigherTL$group.focal[which(HigherTL$group.focal =="Beetles")] <- "Beetle"
HigherTL$group.focal[which(HigherTL$group.focal =="Caterpillars")] <- "Caterpillar"
HigherTL$group.focal[which(HigherTL$group.focal =="Snails")] <- "Snail"
HigherTL$group.focal[which(HigherTL$group.focal =="Bugs")] <- "Bug"

HigherTL.id.group <- levels(as.factor(HigherTL$group.focal))

# divide the groups present in HigherTL.id.group throught the 4 following categories
id.floral.visitors <- c("Bee","Butterfly","Dragonfly","Fly","Neuroptera","Wasp")
  id.herbivor <- c("Bug","Caterpillar","Grasshoppers","Larva","Mite","Slater","Snail","Trips")
  id.herbivor.and.floral.visitors <- c("Ant","Beetle")
id.predator <- c("Mantis","Spider")

# check if all the levels of HigherTL.id.group are included in on of the 4 above categories
HigherTL.id.group[which(! HigherTL.id.group %in% c(id.floral.visitors,id.herbivor,
                                       id.herbivor.and.floral.visitors,id.predator))]
# should be equal to character(0)
# Correct the columns of Sampling based on the above categories
HigherTL$id.inter <- NA
HigherTL$id.inter[which(HigherTL$group.focal %in% id.floral.visitors)]  <- "floral.visitors"
HigherTL$id.inter[which(HigherTL$group.focal %in% id.herbivor)]  <- "herbivor"
HigherTL$id.inter[which(HigherTL$group.focal %in% id.herbivor.and.floral.visitors)]  <- "herbivor.and.floral.visitors"
HigherTL$id.inter[which(HigherTL$group.focal %in% id.predator)]  <- "predator"



# ----functional group---- 
#Data set not ready 

Bee <- c("Andrena","Apis","Osmia","Eucera_sp","Hymenoptera","Lasioglossum_sp","Andrenidae",
          "Andrena cinerea","Apidae", "Eucera sp.","Halictidae", "Lassioglosum immunitum") # social bees and solidar bees
Ants <- c("Formicidae")
Wasp <- c("Chrysididae","Chrysotoxum_sp")
Hymenoptera_rare <- c("Apoidea")
Hymenoptera <- c(Bee, Ants,Hymenoptera_rare)


Butterfly <- c("Euchloe_sp","Lasiocampa_trifolii","Lepidoptera","Hesperiidae",
               "Thymelicus sp.","Lasiocampidae","Lasiocampa trifolii","Pieridae", "Colias croceus")
Lepidoptera <- c(Butterfly)


Flies <- c("Diptera")
Humbleflies <- c("Bombyliidae")
House_flies <- c("Calliphoridae ","Musca_sp","Muscidae",
                 "Sarcophaga_sp")
Small_flies <- c("Dilophus_sp","Bibionidae","Dilophus_sp.",
                 "Empididae","Stratiomyidae","Nemotelus_sp.","Odontomyia_sp.",
                 "Ulidiidae")
Hoverflies <- c("Episyrphus_balteatus","Episyrphus_sp","Eupeodes_corollae","Eupeodes_sp",
                "Lomatia_sp","Nemotelus_sp","Sphaerophoria_scripta","Sphaerophoria_sp",
                "Syrphidae","Chrysotoxum sp.","Eristalis_sp.","Scaeva_sp.","Melanostoma_sp.")
Diptera <- c(Humbleflies,House_flies,Small_flies,Hoverflies)



Small_beetles <- c("Brassicogethes_sp","Curculionidae","Nitidulidae","Bruchidae",
                   "Meligethes sp.","Chrysomelidae") 
Big_beetles <- c("Cerambycidae","Lagorina_sericea","Meloidae")
Beetles <- c("Coleoptera")

Medium_beetles <- c("Coleoptera","Malachius_bipustulatus",
                    "Melyridae","Mordellidae","Oedemeridae",
                    "Psilothrix_viridicoerulea","Anthaxia_semicuprea","Buprestidae",
                    "Cantharidae","Buprestidae","Anthaxia semicuprea","Cantharidae",
                    "Chrysomelidae", "Chrysomelidae","Meloidae",
                    "Lagorina sericea","Malachius_bipustulatus",
                    "Melyridae","Psilothrix viridicoerulea","Mordellidae","Oedemeridae") # ~flower_beetles

Coleoptera <- c(Small_beetles,Big_beetles,Medium_beetles)

#check01: check that all the floral visitors id are inlcuded in the groups
HigherTL.id.funct <- levels(as.factor(HigherTL$id.focal))

#HigherTL_family <- c(Andrena, Apis, Osmia)
#HigherTL.id.funct[which(!HigherTL.id.funct %in% HigherTL_family)]

HigherTL_functional <- c(Bee,Ants,Wasp,Hymenoptera_rare,
                         Butterfly,Flies,Humbleflies,
                         House_flies,Small_flies,Hoverflies,
                         Small_beetles,Big_beetles,Beetles,Medium_beetles)


HigherTL.id.funct[which(!HigherTL.id.funct %in% HigherTL_functional)]

HigherTL_order <- c(Hymenoptera,Lepidoptera,Diptera,Coleoptera)
HigherTL_family[which(HigherTL_family%in% HigherTL.id)]


################################################################################
# Gather all dataframes
#################################################################################
#check:
HigherTL[which(is.na(HigherTL$id.inter)),]
HigherTL <- as.data.frame(HigherTL)

# ----construct dataframe of higherTL based on the different categories---

HigherTL.wide.int.type <- HigherTL %>% 
  subset(select=c("id.inter","id.plant" ,
                   "year","plot","subplot","Mean")) %>% 
  group_by(id.inter,id.plant,year, plot, subplot) %>%
  dplyr::summarize(Mean = mean(Mean, na.rm=TRUE)) %>%
  ungroup() %>%
  spread(id.inter,Mean)

str(HigherTL.wide.int.type)

Plant.and.id.interaction <-  left_join(competition.plant,HigherTL.wide.int.type,
                                 by = c("year", "plot", "subplot", "id.plant"),
                                keep = F)

  
#test <-  left_join(competition.plant,HigherTL_wide,by = c("year", "plot", "subplot", "focal"),
#                   keep = F)


# no HigherTL before 2017? or after 2018 ?
# how to group the herbivor and the floral visitor? 

