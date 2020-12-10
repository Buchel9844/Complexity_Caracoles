################################################################################
# Plant abundance
#################################################################################

abundance.plant <- read.csv("Caracoles/data/abundances_ME_wide.csv", sep=";")
str(abundance.plant)
competition.plant <- read.csv("Caracoles/data/competition_ME_wide.csv", sep=";")
str(competition.plant)

################################################################################
# HigherTL & floral visitor
#################################################################################
#---- Import the Data ----
herbivor2019 <- read.csv("Caracoles/raw_data/Metadata_Herbivorous_2019.csv", sep=";")
herbivor2020 <- read.csv("Caracoles/raw_data/raw_Herbivorous_data_2020.csv", sep=",")
str(herbivor2019)
herbivor2020$Day <- as.integer(herbivor2020$Day)
names(herbivory2020)[which(!names(herbivory2020) %in% names(herbivory2019))]
herbivory2019[,names(herbivory2020)[which(!names(herbivory2020) %in% names(herbivory2019))]] <- NA
herbivor <- bind_rows(herbivor2019,herbivor2020)

levels(as.factor(herbivor$ID))


floral_visitors2016_2019 <- read.csv("Caracoles/raw_data/Metadata_Pollinators_2019_2016_bueno.csv", sep=";")
floral_visitors2020 <- read.csv("Caracoles/raw_data/raw_Pollinators_2020_1.csv", sep=";")
names(floral_visitors2020)[which(!names(floral_visitors2020) %in% names(floral_visitors2016_2019))]

HigherTL <- bind_rows(herbivory,floral_visitors)
str(HigherTL)
# compute the mean of frequency of each herbivor for a specific set of sp1-year-plot-subplot
HigherTL <- HigherTL %>% 
  subset(sp1 != "soil") %>% 
  group_by(sp1,sp2,year, plot, subplot) %>% 
  dplyr::summarize(Mean = mean(frequency, na.rm=TRUE),
                   interaction_type=interaction_type) %>%
  ungroup()
str(HigherTL)

names(HigherTL) <- c("focal","sp2", "year","plot","subplot","Mean","interaction_type")

# for check 00
check00 <- HigherTL[which(HigherTL$focal=="HOMA"&
                             HigherTL$year== 2018 &
                             HigherTL$plot == 1 &
                             HigherTL$subplot =="C3"), c("sp2","Mean")]


# spread sp2 
HigherTL_wide <- HigherTL %>% 
  group_by(focal,year, plot, subplot) %>% 
  spread(sp2,Mean) %>%
  ungroup()

#check that the frequency found in the HigherTL_wide  
# correspond to the ones from the initial data HigherTL (check00)
check00sp2 <- check00$sp2
HigherTL_wide[which(HigherTL_wide$focal=="HOMA"&
                       HigherTL_wide$year== 2018 &
                       HigherTL_wide$plot == 1 &
                       HigherTL_wide$subplot =="C3"),check00sp2]
# compare with:
check00

#---- Regroup different group ----
HigherTL.id <- levels(as.factor(HigherTL$sp2))

# family group
Andrena <- c("Andrena_cinerea_elliptica","Andrena_hesperia","Andrena_humilis","Andrena_sp",
             "Andrenidae")
Apis <- c("Apidae","Apis_mellifera")

Osmia <- c("Osmia_ligurica","Osmia_submicans")

# functional group
Hymenoptera <- c(Bee, Ants,Hymenoptera_rare)
Bees <- c("Andrena","Apis","Osmia","Eucera_sp","Hymenoptera","Lasioglossum_sp","Andrenidae",
          "Andrena cinerea","Apidae", "Eucera sp.","Halictidae", "Lassioglosum immunitum") # social bees and solidar bees
Ants <- c("Formicidae")
Wasp <- c("Chrysididae","Chrysomelidae","Chrysotoxum_sp")
Hymenoptera_rare <- c("Apoidea")

Lepidoptera <- c(Butterfly)
Butterfly <- c("Euchloe_sp","Lasiocampa_trifolii","Lepidoptera","Hesperiidae",
               "Thymelicus sp.","Lasiocampidae","Lasiocampa trifolii","Pieridae", "Colias croceus")

Diptera <- c(Humbleflies,House_flies,Small_flies,Hoverflies)
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


Coleoptera <- c(Small_beetles,Big_beetles,Flower_beetles)
Small_beetles <- c("Brassicogethes_sp","Curculionidae","Nitidulidae","Bruchidae",
                   "Meligethes sp.")
Big_beetles <- c("Cerambycidae","Lagorina_sericea","Meloidae")
beetles <- c("Coleoptera")

Medium_beetles <- c("Coleoptera","Malachius_bipustulatus",
                    "Melyridae","Mordellidae","Oedemeridae",
                    "Psilothrix_viridicoerulea","Anthaxia_semicuprea","Buprestidae",
                    "Cantharidae","Buprestidae","Anthaxia semicuprea","Cantharidae",
                    "Chrysomelidae", "Chrysomelidae","Meloidae",
                    "Lagorina sericea","Malachius_bipustulatus",
                    "Melyridae","Psilothrix viridicoerulea","Mordellidae","Oedemeridae") # ~flower_beetles

# Herbivor or Pollinator
Herbivor <- 
Pollinator <- 
#check01: check that all the floral visitors id are inlcuded in the groups

HigherTL_family <- c(Andrena, Apis, Osmia)
HigherTL.id[which(!HigherTL.id %in% HigherTL_family)]

HigherTL_functional <- c(Bees, Butterfly,Humbleflies,Hymenoptera_rare,Wasp,
                                Small_beetles,House_flies,Big_beetles,Medium_beetles,
                                Small_flies,Hoverflies,Ants)
HigherTL.id[which(!HigherTL.id %in% HigherTL_functional)]

HigherTL_order <- c(Hymenoptera,Lepidoptera,Diptera,Coleoptera)
HigherTL_family[which(HigherTL_family%in% HigherTL.id)]




################################################################################
# Gather all dataframes
#################################################################################

test <-  left_join(competition.plant,HigherTL_wide,by = c("year", "plot", "subplot", "focal"),
                   keep = F)


# no HigherTL before 2017? or after 2018 ?
# how to group the herbivor and the floral visitor? 

