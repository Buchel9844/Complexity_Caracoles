################################################################################
# Groups of floral_visitor and Herbivor
#################################################################################
#---- floral visitor----
#FVtoKeep_species <- groupinglevels$grouping[which(groupinglevels$HTL =="floral.visitor" &
#                                                    groupinglevels$complexity == "species" & (groupinglevels$YesNo == "yes"|groupinglevels$YesNo == "maybe"))]
#FVtoKeep_family <- groupinglevels$grouping[which(groupinglevels$HTL =="floral.visitor" &
#                                                  groupinglevels$complexity == "family" & (groupinglevels$YesNo == "yes"|groupinglevels$YesNo == "maybe"))]
#FVtoKeep_group <- groupinglevels$grouping[which(groupinglevels$HTL =="floral.visitor" &
#                                                 groupinglevels$complexity == "functional.groups" & (groupinglevels$YesNo == "yes"|groupinglevels$YesNo == "maybe"))]

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

test.fv.sp <- floral_visitor_species %>%
  group_by(id_final) %>%
  summarise(mean.sp = mean(number_visits),
            var.sp = var(number_visits)) %>%
  filter(0<round(var.sp,digits=2)) 
#filter(mean.sp<=var.sp) 
FVtoKeep_species <- test.fv.sp$id_final[test.fv.sp$id_final %in% FVtoKeep_species]

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
test.fv.fam <- floral_visitor_family %>%
  group_by(family) %>%
  summarise(mean.sp = mean(number_visits),
            var.sp = var(number_visits)) %>%
  #filter(mean.sp<=var.sp) %>%
  filter(0<round(var.sp,digits=2)) 
FVtoKeep_family <- test.fv.fam$family[test.fv.fam$family %in% FVtoKeep_family]

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

test.fv.group <- floral_visitor_group %>%
  group_by(group) %>%
  summarise(mean.sp = mean(number_visits,na.rm = T),
            var.sp = var(number_visits,na.rm = T)) %>%
  #filter(mean.sp<=var.sp) %>%
  filter(0<round(var.sp,digits=2)) 
FVtoKeep_group <- test.fv.group$group[test.fv.group$group %in% FVtoKeep_group]

floral_visitor_group <- dplyr::filter(floral_visitor_group , 
                                      group %in% FVtoKeep_group)


floral_visitor_group  <-  spread(  floral_visitor_group,group, number_visits)

floral_visitor_group[is.na(floral_visitor_group)] <- 0

head(floral_visitor_group)
length(which(!names( floral_visitor_group) %in% c("plot","year","subplot","code.plant")))

#---- herbivor ----
#HtoKeep_species <- groupinglevels$grouping[which(groupinglevels$HTL =="herbivore" &
#                                                  groupinglevels$complexity == "species" & (groupinglevels$YesNo == "maybe"|groupinglevels$YesNo == "yes"))]
#HtoKeep_family <- groupinglevels$grouping[which(groupinglevels$HTL =="herbivore" &
#                                                 groupinglevels$complexity == "family" & (groupinglevels$YesNo == "maybe"|groupinglevels$YesNo == "yes"))]
#HtoKeep_group <- groupinglevels$grouping[which(groupinglevels$HTL =="herbivore" &
#                                                  groupinglevels$complexity == "functional.groups" & (groupinglevels$YesNo == "maybe"|groupinglevels$YesNo == "yes"))]

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
test.h.sp <- herbivore_species %>%
  group_by(id_final) %>%
  summarise(mean.sp = mean(number_animal),
            var.sp = var(number_animal)) %>%
  #filter(mean.sp<=var.sp) %>%
  filter(0<round(var.sp,digits=2)) 

HtoKeep_species <- test.h.sp$id_final[test.h.sp$id_final %in% HtoKeep_species]

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
test.h.fam <- herbivore_family %>%
  group_by(family) %>%
  summarise(mean.sp = mean(number_animal),
            var.sp = var(number_animal)) %>%
  #filter(mean.sp<=var.sp) %>%
  filter(0<round(var.sp,digits=2)) 

HtoKeep_family <- test.h.fam$family[test.h.fam$family %in% HtoKeep_family]

herbivore_family <- dplyr::filter(herbivore_family, family %in% HtoKeep_family)
herbivore_family  <-  spread(   herbivore_family,family, number_animal)
herbivore_family[is.na( herbivore_family)] <- 0
head(herbivore_family )
length(which(!names(herbivore_family) %in% c("plot","year","subplot","code.plant")))



herbivore_group <- dplyr::select( herbivore,
                                  c("day","month","year","plot","subplot",
                                    "group","number_animal","code.plant")) 
herbivore_group <-stats::aggregate(number_animal ~ group + year + plot  + code.plant +subplot ,
                                   data=   herbivore_group,
                                   sum)
test.h.group <- herbivore_group %>%
  group_by(group) %>%
  summarise(mean.sp = mean(number_animal),
            var.sp = var(number_animal)) %>%
  #filter(mean.sp<=var.sp) %>%
  filter(0<round(var.sp,digits=2)) 

HtoKeep_group <- test.h.group$group[test.h.group$group %in% HtoKeep_group]

herbivore_group <- dplyr::filter(herbivore_group, group %in% HtoKeep_group)

herbivore_group  <-  spread(   herbivore_group,group, number_animal)

herbivore_group[is.na( herbivore_group)] <- 0

head(herbivore_group)
length(which(!names(herbivore_group) %in% c("plot","year","subplot","code.plant")))


H.obs.variation <- bind_rows(test.h.sp %>%
                               dplyr::filter(id_final%in% HtoKeep_species) %>%
                               mutate(grouping="species") %>%
                               rename("ID"=id_final),
                             test.h.fam %>%
                               dplyr::filter(family %in% HtoKeep_family) %>%
                               mutate(grouping="family") %>%
                               rename("ID"=family),
                             test.h.group %>%
                               dplyr::filter(group %in% HtoKeep_group) %>%
                               mutate(grouping="group") %>%
                               rename("ID"=group))
write.csv(H.obs.variation ,
          file=paste0(home.dic,"results/H.obs.variation.csv"))

FV.obs.variation <- bind_rows(test.fv.sp %>%
                                dplyr::filter(id_final %in% all_of(FVtoKeep_species)) %>%
                                mutate(grouping="species") %>%
                                rename("ID"=id_final),
                              test.fv.fam %>%
                                dplyr::filter(family %in% FVtoKeep_family) %>%
                                mutate(grouping="family") %>%
                                rename("ID"=family),
                              test.fv.group %>%
                                dplyr::filter(group %in% FVtoKeep_group) %>%
                                mutate(grouping="group") %>%
                                rename("ID"=group))
write.csv(FV.obs.variation ,
          file=paste0(home.dic,"results/FV.obs.variation.csv"))
