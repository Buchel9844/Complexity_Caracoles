################################################################################
# Groups of floral_visitor and Herbivor
#################################################################################
#---- floral visitor----
FVtoKeep_species <- groupinglevels$grouping[which(groupinglevels$HTL =="floral.visitor" &
                                                    groupinglevels$complexity == "species" & groupinglevels$YesNo == "yes")]
FVtoKeep_family <- groupinglevels$grouping[which(groupinglevels$HTL =="floral.visitor" &
                                                   groupinglevels$complexity == "family" & groupinglevels$YesNo == "yes")]
FVtoKeep_group <- groupinglevels$grouping[which(groupinglevels$HTL =="floral.visitor" &
                                                  groupinglevels$complexity == "functional.groups" & groupinglevels$YesNo == "yes")]


names(floral_visitor)[which(names(floral_visitor)=="plant")] <- "code.plant"

floral_visitor_species <- dplyr::select(floral_visitor,c("day","month","year","plot","subplot",
                                                         "id_final","number_visits","code.plant")) 
floral_visitor_species <-stats::aggregate(number_visits ~ id_final +year +plot +subplot + code.plant,
                                          data= floral_visitor_species,
                                          sum)
floral_visitor_species <-  dplyr::filter(floral_visitor_species,
                                         id_final %in% FVtoKeep_species)


floral_visitor_species <-  spread( floral_visitor_species,id_final, number_visits)
floral_visitor_species[is.na(floral_visitor_species)] <- 0
head(floral_visitor_species)

length(which(!names(floral_visitor_species) %in% c("plot","year","subplot","code.plant")))

floral_visitor_family <- dplyr::select(floral_visitor,c("plot","subplot","year",
                                                        "family","number_visits","code.plant")) 
floral_visitor_family <-stats::aggregate(number_visits ~ family + year + plot + code.plant +subplot ,
                                         data=  floral_visitor_family,
                                         sum)

floral_visitor_family  <- dplyr::filter(floral_visitor_family , 
                                        family %in% FVtoKeep_family)


floral_visitor_family  <-  spread(  floral_visitor_family,family, number_visits)
floral_visitor_family[is.na(floral_visitor_family)] <- 0
head( floral_visitor_family )

length(which(!names( floral_visitor_family) %in% c("plot","year","subplot","code.plant")))


floral_visitor_group <- dplyr::select(floral_visitor,c("plot","subplot","year",
                                                       "group","number_visits","code.plant")) 
floral_visitor_group <-stats::aggregate(number_visits ~ group + year +plot + code.plant+subplot ,
                                        data = floral_visitor_group,
                                        sum)

floral_visitor_group <- dplyr::filter(floral_visitor_group , group %in% FVtoKeep_group)


floral_visitor_group  <-  spread(  floral_visitor_group,group, number_visits)

floral_visitor_group[is.na(floral_visitor_group)] <- 0

head(floral_visitor_group)
length(which(!names( floral_visitor_group) %in% c("plot","year","subplot","code.plant")))

#---- herbivor ----
HtoKeep_species <- groupinglevels$grouping[which(groupinglevels$HTL =="herbivore" &
                                                   groupinglevels$complexity == "species" & groupinglevels$YesNo == "yes")]
HtoKeep_family <- groupinglevels$grouping[which(groupinglevels$HTL =="herbivore" &
                                                  groupinglevels$complexity == "family" & groupinglevels$YesNo == "yes")]
HtoKeep_group <- groupinglevels$grouping[which(groupinglevels$HTL =="herbivore" &
                                                 groupinglevels$complexity == "functional.groups" & groupinglevels$YesNo == "yes")]



levels(as.factor(herbivorypredator$type))
names(herbivorypredator)[which(names(herbivorypredator)=="plant")] <- "code.plant"
herbivore <- herbivorypredator[which(herbivorypredator$type=="H/FV"|herbivorypredator$type=="H"),]


herbivore_species <- dplyr::select( herbivore, c("day","month","year","plot","subplot",
                                                 "id_final","number_animal","code.plant")) 
herbivore_species <-stats::aggregate(number_animal ~ id_final + year +plot + code.plant +subplot ,
                                     data=  herbivore_species,
                                     sum)
herbivore_species <- dplyr::filter(herbivore_species, 
                                   id_final %in% HtoKeep_species)

herbivore_species <-  spread(  herbivore_species,id_final, number_animal)
herbivore_species[is.na( herbivore_species)] <- 0
head( herbivore_species)
length(which(!names(herbivore_species) %in% c("plot","year","subplot","code.plant")))

#groupinglevels <- rbind(groupinglevels,data.frame(complexity="species", HTL="herbivore",
#                                                  grouping = names(herbivore_species)[which(!names(herbivore_species) %in% c("plot","year","subplot","code.plant"))]))


herbivore_family <- dplyr::select( herbivore,c("day","month","year","plot","subplot",
                                               "family","number_animal","code.plant")) 
herbivore_family <-stats::aggregate(number_animal ~ family + year +plot  + code.plant +subplot ,
                                    data=   herbivore_family,
                                    sum)

herbivore_family <- dplyr::filter(herbivore_family, family %in% HtoKeep_family)
herbivore_family  <-  spread(   herbivore_family,family, number_animal)
herbivore_family[is.na( herbivore_family)] <- 0
head(herbivore_family )
length(which(!names(herbivore_family) %in% c("plot","year","subplot","code.plant")))



herbivore_group <- dplyr::select( herbivore,c("day","month","year","plot","subplot",
                                              "group","number_animal","code.plant")) 
herbivore_group <-stats::aggregate(number_animal ~ group + year + plot  + code.plant +subplot ,
                                   data=   herbivore_group,
                                   sum)
herbivore_group <- dplyr::filter(herbivore_group, group %in% HtoKeep_group)

herbivore_group  <-  spread(   herbivore_group,group, number_animal)

herbivore_group[is.na( herbivore_group)] <- 0

head(herbivore_group)
length(which(!names(herbivore_group) %in% c("plot","year","subplot","code.plant")))


