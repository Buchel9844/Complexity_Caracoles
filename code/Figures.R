# This script will make the figures of chap 1, of the natural data part, number- XXX of the manuscript

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

#remotes::install_github("coolbutuseless/ggpattern",force = TRUE,dependencies=T)
library(ggpattern)

home.dic <- "/home/lbuche/Eco_Bayesian/Complexity_caracoles/"
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/"


colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
scales::show_col(colorBlindGrey8)

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
scales::show_col(safe_colorblind_palette)
#---- 3. Natural data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#---- 3.0.0 All potential interactions----

SpData <- read.csv(paste0(home.dic,"data/SpData.csv")) 

complexitylevel.plant.species <- levels(as.factor(plant.class$species))
complexitylevel.plant.family <- levels(as.factor(plant.class$family))
complexitylevel.plant.class <- levels(as.factor(plant.class$class))

groupinglevels <- read.csv(paste0(home.dic,"data/groupinglevels.csv"), sep=",")

FVtoKeep_species <- groupinglevels$grouping[which(groupinglevels$HTL =="floral.visitor" &
                                                    groupinglevels$complexity == "species" & groupinglevels$YesNo == "yes")]
FVtoKeep_family <- groupinglevels$grouping[which(groupinglevels$HTL =="floral.visitor" &
                                                   groupinglevels$complexity == "family" & groupinglevels$YesNo == "yes")]
FVtoKeep_group <- groupinglevels$grouping[which(groupinglevels$HTL =="floral.visitor" &
                                                  groupinglevels$complexity == "functional.groups" & groupinglevels$YesNo == "yes")]
HtoKeep_species <- groupinglevels$grouping[which(groupinglevels$HTL =="herbivore" &
                                                   groupinglevels$complexity == "species" & groupinglevels$YesNo == "yes")]
HtoKeep_family <- groupinglevels$grouping[which(groupinglevels$HTL =="herbivore" &
                                                  groupinglevels$complexity == "family" & groupinglevels$YesNo == "yes")]
HtoKeep_group <- groupinglevels$grouping[which(groupinglevels$HTL =="herbivore" &
                                                 groupinglevels$complexity == "functional.groups" & groupinglevels$YesNo == "yes")]


summary.potential.interaction <- data.frame(levelparametercomplexity=c("high","medium","low"),
                                            grouping=c("species","family","functional group"),
                                            plant = c(paste0("e.g., ",complexitylevel.plant.species[3]),
                                                      paste0("e.g., ",complexitylevel.plant.family[3]),
                                                      paste0("e.g., ",complexitylevel.plant.class[1])),
                                            "number\nplant\ninteraction" = c(length(complexitylevel.plant.species),
                                                                             length(complexitylevel.plant.family),
                                                                             length(complexitylevel.plant.class)),
                                            floralvisitor = c(paste0("e.g., ",FVtoKeep_species[2]),
                                                              paste0("e.g., ",FVtoKeep_family[2]),
                                                              paste0("e.g., ",FVtoKeep_group[2])),
                                            "number\nfloralvisitor\ninteraction" = c(length(FVtoKeep_species),
                                                                                     length(FVtoKeep_family),
                                                                                     length(FVtoKeep_group)),
                                            herbivore = c(paste0("e.g., ",HtoKeep_species[3]),
                                                          paste0("e.g., ",HtoKeep_family[3]),
                                                          paste0("e.g., ",HtoKeep_group[3])),
                                            "number\nherbivore\ninteraction" = c(length(HtoKeep_species),
                                                                                 length(HtoKeep_family),
                                                                                 length(HtoKeep_group)))



write.csv(summary.potential.interaction,
          file = paste0(home.dic,"results/summary.potential.interaction.csv"))


#---- 3.0 CSV for inclusion----

df_inclusion_nat <- data.frame()
for( focal in c("LEMA","CHFU","HOMA","CETE")){ # "CHFU","HOMA",
  for( year in c("2019",'2020','2021')){
    for( complexity.level in c(1:3)){
      
      
      complexity.animal <- c("group","family","species")[complexity.level]
      complexity.plant <-c("class","family","code.plant")[complexity.level]
      
      load(paste0(home.dic,"results/Inclusion",focal,"_",
                  year,"_",complexity.plant,"_",complexity.animal,".RData"))
      
      assign(paste0("Inclusion",focal,"_",
                    year,"_",complexity.plant,"_",complexity.animal),Inclusion_all)
      rm(Inclusion_all)
      
      df_inclusion_n <- as.data.frame(get(paste0("Inclusion",focal,"_",
                                                 year,"_",complexity.plant,"_",complexity.animal))$interaction) 
      
      df_inclusion_n$focal <- focal
      df_inclusion_n$year <- as.numeric(year)
      df_inclusion_n$complexity.animal <- complexity.animal
      df_inclusion_nat <- bind_rows(df_inclusion_nat,    df_inclusion_n)
      
    }
  }
}
write.csv(df_inclusion_nat,
          file = paste0(home.dic,"results/Chapt1_Inclusion_parameters.csv"))

summary.interactions <- read.csv(paste0(home.dic,"results/Chapt1_Inclusion_parameters.csv"))
df_inclusion_nat <- summary.interactions 
#n.inclus <- grep(pattern = '^n.*inclus$', x = colnames(df_inclusion_nat), value = T)

n.inclus <- c("n.competitors_plant_inclus","n.competitors_FV_inclus","n.competitors_H_inclus",
              "n.HOIs_plant_inclus","n.HOIs_H_inclus",
              "n.HOIs_FV_inclus")
inclus <- c("competitors_plant_inclus","competitors_FV_inclus","competitors_H_inclus",
            "HOIs_plant_inclus","HOIs_H_inclus",
            "HOIs_FV_inclus")
inclus.short <- grep(pattern = '.*inclus$', x = colnames(df_inclusion_nat), value = T)
inclus.short <- inclus.short[which(!inclus.short %in% n.inclus)]
inclus.not.short <- inclus[which(!inclus %in% inclus.short)]
df_inclusion_nat[,inclus.not.short] <-NA
n.potential <- c("n.competitors_plant","n.interactors_FV","n.interactors_H",
                 "n.HOIs_plant","n.HOIs_H",
                 "n.HOIs_FV")

potential <- c("competitors_plant","interactors_FV","interactors_H",
               "HOIs_plant","HOIs_H",
               "HOIs_FV")
df_inclusion_nat_short <- df_inclusion_nat[which(rowSums(df_inclusion_nat[,n.inclus],na.rm = T)>0),]
view(df_inclusion_nat_short)
extract.vert <- function(mat,vec,keyname,value,vecname){
  mat.2 <- mat %>% 
    gather(all_of(vec), 
           key = keyname , value = value) %>%
    mutate(keyname = rep(vecname, each = nrow(mat)))
  return(mat.2)
}

df_inclusion_nat_vertical <- df_inclusion_nat %>% 
  gather(all_of(potential), 
         key="interaction", value= "potential.interactor") %>%
  #mutate(interaction = rep(inclus, times = nrow(df_inclusion_nat))) %>%
  select(-all_of(c(inclus,n.inclus,n.potential))) %>%
  mutate(interactor = extract.vert(df_inclusion_nat,inclus,
                                   "interaction","identity",inclus)[,"value"],
         n.number = extract.vert(df_inclusion_nat,n.inclus,
                                 "interaction.number","n.number",n.inclus)[,"value"],
         n.potential.interactor = extract.vert(df_inclusion_nat,n.potential ,
                                               "interaction.pot","n.identity.potential",n.potential)[,"value"]) %>%
  mutate(interaction = case_when(interaction == "competitors_plant" ~ "Plant - plant",
                                 interaction == "HOIs_FV" ~ "Plant - plant - floral visitor",
                                 interaction== "HOIs_H" ~ "Plant - plant - herbivore",
                                 interaction== "HOIs_plant" ~ "Plant - plant - plant",
                                 interaction== "interactors_FV" ~ "Plant - floral visitor",
                                 interaction== "interactors_H" ~ "Plant - herbivore"))


pal <- wes_palette("Zissou1", 100, type = "continuous")

df_inclusion_nat_ratio <- df_inclusion_nat_vertical %>%
  mutate(ratio =  (n.number/(n.potential.interactor-1)) ) %>%
  mutate(year =  as.integer(year)) %>%
  #mutate(interactor =  as.character(interactor)) %>%
  mutate(interaction =factor(interaction))%>%
  mutate(complexity_plant =factor(complexity_plant,
                                  levels=c("class","family","code.plant")))

df_inclusion_nat_heatmap <- ggplot(df_inclusion_nat_ratio,
                                   aes(x=interaction,y=as.factor(year),fill= as.numeric(ratio))) + 
  scale_fill_gradientn(colours =  c("lightgrey",wes_palette("Zissou1", type = "continuous")), #c("lightgrey", pal),
                       na.value = "white",
                       values =seq(0,100,0.01),
                       lim=c(0,100)) + 
  geom_tile()+ ylab("") + facet_grid(focal~ complexity_plant) +
  theme_bw() + theme(axis.text.x = element_text(angle=45))

ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/df_inclusion_plot_ratio_ALL.pdf"),
       dpi="retina",scale=1, 
       df_inclusion_nat_heatmap )

scale_fill_colorblind7 = function(.ColorList = 2L:8L, ...){
  scale_fill_discrete(..., type = colorblind_pal()(8)[.ColorList])
}

# Plot
df_inclusion_nat_vertical$interaction <- factor(df_inclusion_nat_vertical$interaction,
                                                levels = c("Plant - plant", "Plant - herbivore", "Plant - floral visitor",
                                                           "Plant - plant - plant","Plant - plant - herbivore","Plant - plant - floral visitor"))
df_inclusion_nat_vertical.2 <- df_inclusion_nat_vertical %>%
  mutate(complexity_plant = case_when(complexity_plant == "class" ~ "Class",
                                      complexity_plant == "family" ~ "Family",
                                      complexity_plant == "code.plant" ~ "Species")) %>%
  mutate(type.interaction = case_when((interaction =="Plant - plant" | interaction==  "Plant - floral visitor"|
                                         interaction== "Plant - herbivore") ~ "Direct interactions",
                                      (interaction ==  "Plant - plant - floral visitor"| interaction== "Plant - plant - herbivore"|
                                         interaction== "Plant - plant - plant") ~ "HOIs"))


df_inclusion_nat_bar <- ggplot(df_inclusion_nat_vertical.2,
                               aes(x=as.factor(year),y=n.number,
                                   fill = interaction,pattern = type.interaction)) +
  geom_bar_pattern(position = "stack",stat= "identity",
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) +
  scale_pattern_manual(values = c("Direct interactions" = "none",
                                  "HOIs"= "stripe")) + 
  facet_grid(focal~ complexity_plant,scale="free") +
  theme_bw() + #theme(axis.text.x = element_text(angle=45)) + 
  scale_fill_manual(values=c("#117733","#DDCC77","#88CCEE",
                             "#117733","#DDCC77","#88CCEE")) + 
  labs(x="year",y="Number of species specific interactions",
       fill =" interaction") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = c("none","none","none",
                                                             "stripe","stripe","stripe")
         ))) #+ 
# geom_text(aes(label=interactor, y=position.label.2),
#      size=df_inclusion_nat_vertical.2$size.label, angle=df_inclusion_nat_vertical.2$angle.label)
#_pattern(position = "stack",stat= "identity",

df_inclusion_nat_bar 
ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/df_inclusion_plot_bar.pdf"),
       dpi="retina",
       width = 21,
       height = 16,
       units = c("cm"),
       df_inclusion_nat_bar)



#---- 4. Results for interactions values ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#---- 4.1 Figure of pairwise interactions distributions and Intra Vs Inter ----
#extraction generic parameters
df_param <- NULL

for( focal in c("CETE","LEMA","HOMA","CHFU")){ # "CHFU","HOMA","CETE"
  for( year in c("2019",'2020','2021')){
    for( complexity.level in c(1:3)){
      
      complexity.animal <- c("group","family","species")[complexity.level]
      complexity.plant <-c("class","family","code.plant")[complexity.level]
      
      df_alpha_n <- read.csv(paste0(project.dic,"results/parameters/Parameters_",
                                    focal,"_", year,"_",complexity.plant,
                                    complexity.animal,"_FinalFit.csv"))
      
      if(is.null(df_param)){df_param  <- df_alpha_n
      }else{
        df_param <- full_join(df_param,df_alpha_n)
      }
    }
  }
}

# species specific parameters
df_param_hat <- NULL
inclusion.interaction.type <- unique(df_inclusion_nat_vertical[,c("interaction","interactor")])
for ( i in 1:nrow(df_inclusion_nat_short)){
  year <- df_inclusion_nat_short[i,'year']
  focal <- df_inclusion_nat_short[i,'focal']
  
  complexity.plant <- df_inclusion_nat_short[i,'complexity_plant']
  complexity.animal <- df_inclusion_nat_short[i,'complexity.animal']
  
  df_param_hat_n <- read.csv(paste0(project.dic,"results/parameters/Parameters_hat_",
                                    focal,"_", year,"_",complexity.plant,
                                    complexity.animal,"_FinalFit.csv"))
  name.species.specific <- names(df_param_hat_n)[which(!names(df_param_hat_n) %in% c("X","focal","year",
                                                                                     "complexity.plant",
                                                                                     "complexity.animal" ))]
  
  df_param_hat_n <-   df_param_hat_n %>%
    gather(all_of(name.species.specific),
           key="parameter_hat",value="estimate_hat") %>%
    mutate(parameter_hat =   sub(".*\\_", "", parameter_hat ))
  type_n <- inclusion.interaction.type$interaction[which(inclusion.interaction.type$interactor %in% levels(as.factor(df_param_hat_n$parameter_hat)))]
  df_param_hat_n$parameter <- rep(type_n, each= nrow(df_param_hat_n)/length(type_n))
  
  if(is.null(df_param_hat)){
    df_param_hat <- df_param_hat_n
  }else{
    df_param_hat <- bind_rows(df_param_hat,df_param_hat_n)
  }
}
head(df_param_hat)
levels(as.factor(df_param_hat$parameter))

gr <- colorRampPalette(c("#009E73"))(200)                      
re <- colorRampPalette(c("#E69F00"))(200)
df_param_vert <- df_param %>%
  gather(alpha_intra,alpha_generic,
         gamma_H_generic,gamma_FV_generic,
         key="parameter",value="estimate") %>%
  filter(!is.na(estimate)) %>%
  mutate(parameter = case_when(parameter=="gamma_H_generic" ~ "Plant - herbivore",
                               parameter=="gamma_FV_generic" ~ "Plant - floral visitor",
                               parameter=="alpha_intra" ~ "Intraspecific",
                               parameter=="alpha_generic" ~ "Plant - plant"))
dplyr::anti_join(df_param_vert,df_param_hat, by = c("X", "focal", "year", "complexity.plant", 
                                                    "complexity.animal", "parameter"))
df_param_all <- full_join(df_param_vert,df_param_hat,
                          by = c("X", "focal", "year", "complexity.plant", 
                                 "complexity.animal", "parameter"),
                          multiple = "all")

write.csv(df_param_all,
          file = paste0(home.dic,"results/Chapt1_Parameters_values.csv"))


df_param_all <- df_param_all %>%
  mutate(parameter = factor(parameter,
                            levels = rev(c("Plant - plant", "Intraspecific",
                                           "Plant - floral visitor","Plant - herbivore"
                            ))))

plot.alphas <- ggplot(df_param_all[which(df_param_all$complexity.plant=="code.plant"),]) +  
  geom_density_ridges_gradient(aes(x=estimate, y=parameter,
                                   fill= after_stat(x)),
                               scale = 1) + 
  geom_boxplot(aes(x=estimate + estimate_hat, y = parameter,
                   color=parameter_hat),
               width=0.6,outlier.shape = NA) +
  facet_grid(as.factor(focal) ~ as.factor(year)) +
  #scale_x_continuous(limits=c(-0.5,0.5)) + 
  #scale_y_discrete(labels=c("Generic Inter-specific",
  #                          "Intra-specific")) + 
  labs(y="Interactions", x= "Estimated distribution",
       color = "Species-specific \n parameters") + 
  scale_color_manual(values=safe_colorblind_palette) +
  scale_fill_gradientn(
    colours=c(re,"white",gr),limits = c(-2, 2),
    #midpoint = 0,
    #space = "Lab",
    na.value = "grey50",
    guide = "none",
    aesthetics = "fill"
  ) +
  theme_minimal() +
  theme(strip.placement = "outside",
        legend.key.size = unit(1, 'cm'),
        title =element_text(size=16),
        axis.text.x= element_text(size=16),
        axis.text.y= element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        strip.text = element_text(size=16),
        panel.border = element_rect(color = "black", 
                                    fill = NA))

plot.alphas
ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/Parameters_distribution_Species.pdf"),
       dpi="retina",
       width = 21,
       height = 16,
       units = c("cm"),
       plot.alphas
)


plot.alphas.small <- df_param_all %>%
  filter(complexity.plant=="family") %>%
  filter(year==2020) %>%
  mutate(focal.name = case_when(focal == "CHFU" ~"Chamaemelum fuscatum",
                                focal == "CETE" ~"Centaurium  tenuiflorum",
                                focal == "HOMA" ~"Hordeum marinum",
                                focal == "LEMA" ~"Leontodon maroccanus"))  %>%
  mutate(parameter = factor(parameter,
                            levels = rev(c("Plant - plant", "Intraspecific",
                                           "Plant - floral visitor","Plant - herbivore"
                            )))) %>%
  ggplot() +  
  geom_density_ridges_gradient(aes(x=estimate, y=parameter,
                                   fill= after_stat(x)),
                               scale = 1) + 
  geom_boxplot(aes(x=estimate + estimate_hat, y = parameter,
                   color=parameter_hat),
               width=0.3,alpha=0.4,outlier.shape = NA) +
  facet_wrap(.~as.factor(focal.name) , ncol=2, nrow=2) +
  #scale_x_continuous(limits=c(-0.5,0.5)) + 
  #scale_y_discrete(labels=c("Generic Inter-specific",
  #                          "Intra-specific")) + 
  labs(y="Interactions", x= "Estimated distribution",
       color = "Family-specific \n parameters") + 
  scale_color_manual(values=rep(c("black"),times=9)) +
  scale_fill_gradientn(
    colours=c(re,"white",gr),limits = c(-2, 2),
    #midpoint = 0,
    #space = "Lab",
    na.value = "grey50",
    guide = "none",
    aesthetics = "fill"
  ) +
  theme_minimal() +
  theme(strip.placement = "outside",
        legend.key.size = unit(1, 'cm'),
        title =element_text(size=16),
        axis.text.x= element_text(size=16),
        axis.text.y= element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        strip.text = element_text(size=16,face="italic"),
        panel.border = element_rect(color = "black", 
                                    fill = NA))

plot.alphas.small
ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/Parameters_distribution_species_fam2020.pdf"),
       # dpi="retina",
       # width = 21,
       # height = 16,
       # units = c("cm"),
       plot.alphas.small
)

#---- 4.2. Figure of pairwise interactions sign percentage ----
df_param_perc  <- NULL
for( focal in c("CETE","LEMA","HOMA","CHFU")){ # "CHFU","HOMA","CETE"
  for( year in c("2019",'2020','2021')){
    for( complexity.level in c(1:3)){
      complexity.animal <- c("group","family","species")[complexity.level]
      complexity.plant <-c("class","family","code.plant")[complexity.level]
      df_param_perc_n  <- NULL
      df_param_perc_n <-  df_param_all[which(df_param_all$focal==focal  & 
                                               df_param_all$year == year &
                                               df_param_all$complexity.plant == complexity.plant),] %>%
        group_by(parameter) %>%
        mutate(positive.perc = length(which(estimate > 0))/length(estimate),
               neg.perc = length(which(estimate < 0))/length(estimate),
               null.perc = length(which(estimate==0))/length(estimate),
               mean.estimate = mean(estimate),
               median.estimate = median(estimate),
               ten.estimate = quantile(estimate, probs = 0.12, na.rm = T),
               ninety.estimate = quantile(estimate, probs = 0.88, na.rm = T),
               var.estimate=var(estimate)) %>%
        mutate(positive.perc_hat = length(which(estimate_hat > 0))/length(estimate_hat),
               neg.perc_hat = length(which(estimate_hat < 0))/length(estimate_hat),
               null.perc_hat = length(which(estimate_hat==0))/length(estimate_hat)) %>%
        mutate(mean.estimate_hat = mean(estimate_hat, na.rm = T),
               median.estimate_hat = median(estimate_hat, na.rm = T),
               ten.estimate_hat = quantile(estimate_hat, probs = 0.12, na.rm = T),
               ninety.estimate_hat = quantile(estimate_hat, probs = 0.88, na.rm = T),
               var.estimate_hat=var(estimate_hat)) %>%
        ungroup() %>%
        mutate(positive.perc_overall = length(which(estimate_hat > 0))/length(estimate_hat),
               positive.perc_mean = mean(positive.perc),
               positive.perc_var = var(positive.perc),
               neg.perc_mean = mean(neg.perc),
               neg.perc_var = var(neg.perc)) %>%
        select(focal,year,complexity.plant,parameter,
               positive.perc, neg.perc,null.perc,
               positive.perc_hat,neg.perc_hat, null.perc_hat,
               positive.perc_mean, positive.perc_var,
               neg.perc_mean, neg.perc_var,
               mean.estimate,median.estimate,
               ten.estimate,ninety.estimate,var.estimate,
               mean.estimate_hat,median.estimate_hat,
               ten.estimate_hat,ninety.estimate_hat,
               var.estimate_hat,positive.perc_overall ) %>%
        unique()
      
      df_param_perc <- bind_rows(df_param_perc,df_param_perc_n)
      
    }
  }
}
str(df_param_perc)  
view(df_param_perc)
df_param_perc$positive.perc_hat[which(df_param_perc$neg.perc_hat==0 & 
                                        df_param_perc$positive.perc_hat==0)] <- NA

df_param_perc$neg.perc_hat[which(df_param_perc$neg.perc_hat==0 & 
                                   is.na(df_param_perc$positive.perc_hat))] <- NA

df_param_perc <- df_param_perc %>%
  mutate(complexity.plant  = factor(complexity.plant,
                                    levels = c("class","family","code.plant")),
         complexity.plant  = case_when(complexity.plant=="class"~ "1.Order",
                                       complexity.plant=="family"~ "2.Family",
                                       complexity.plant=="code.plant"~ "3.Species"),
         parameter = factor(parameter,
                            levels = c("Plant - plant", "Intraspecific",
                                       "Plant - floral visitor","Plant - herbivore"
                            )))

levels(df_param_perc$parameter)

plot_param_perc <- df_param_perc %>%
  mutate(focal.name = case_when(focal == "CHFU" ~"Chamaemelum fuscatum",
                                focal == "CETE" ~"Centaurium  tenuiflorum",
                                focal == "HOMA" ~"Hordeum marinum",
                                focal == "LEMA" ~"Leontodon maroccanus")) %>%
  mutate(parameter = fct_relevel(parameter,"Plant - plant", "Intraspecific",
                                 "Plant - floral visitor","Plant - herbivore")) %>%
  ggplot(aes(y=as.factor(parameter),x=positive.perc, fill=as.factor(parameter))) + 
  geom_rect(xmin=0.4,xmax=0.6,ymin=-Inf,ymax=Inf,
            fill="grey", color="grey", size=0.5, alpha=0.02) + 
  geom_point(size=6, color="black",
             aes(shape=as.factor(year),fill=as.factor(parameter)), alpha=0.8, na.rm = TRUE) + 
  geom_point(size=4,color="red",
             aes(y=as.factor(parameter), x=positive.perc_hat, shape=as.factor(year)), 
             alpha=0.9) + 
  geom_point(size=6,color="black",fill="black",
             aes(y="Sum all", x=positive.perc_mean,shape=as.factor(year)), 
             alpha=0.9) + 
  geom_linerange(aes(y="Sum all", xmin = positive.perc_mean - positive.perc_var^2 , 
                     xmax = positive.perc_mean +positive.perc_var^2))+
  facet_grid(focal.name  ~complexity.plant)  +
  geom_vline(xintercept=0.5,color="grey") + 
  scale_fill_manual(values=rev(c("#332288","#117733","#DDCC77","#88CCEE"))) +
  scale_shape_manual(values=c(21,22,24)) +
  scale_y_discrete(limits=rev(c("Plant - plant", "Intraspecific",
                                "Plant - herbivore", "Plant - floral visitor","Sum all"))) +
  scale_x_continuous("ratio Positive(+)/Negative(-)",
                     breaks=c(0,0.25,0.5,0.75,1),
                     labels = rev(c("100% + \n 0% -", 
                                    "75% + \n25% -",
                                    "50% + \n50% -",
                                    "25% + \n75% -",
                                    "0% + \n100% -"))) +
  theme_minimal() +
  labs(title = "Percentage of Positive vs Negative direct interactions",
       y="",
       shape="year",
       fill="interaction") +
  guides(color= "none",fill="none") +
  theme(panel.grid.minor = element_blank(),
        panel.spacing.y=unit(0.5, "lines") ,
        panel.spacing.x=unit(4,"lines"),
        legend.key.size = unit(1, 'cm'),
        title =element_text(size=16),
        axis.text.x= element_text(size=16),
        axis.text.y= element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        strip.text = element_text(size=16),
        panel.border = element_rect(color = "black", 
                                    fill = NA))
plot_param_perc

ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/Plot_param_perc.pdf"),
       #dpi="retina",
       #width = 25,
       #height = 16,
       #units = c("cm"),
       plot_param_perc
)

plot_param_perc_small <- df_param_perc %>%
  filter(#year== 2020 & 
    complexity.plant=="2.Family") %>%
  filter(!(focal =="HOMA" & parameter== "Plant - floral visitor")) %>%
  mutate(focal.name = case_when(focal == "CHFU" ~"Chamaemelum \nfuscatum",
                                focal == "CETE" ~"Centaurium  \ntenuiflorum",
                                focal == "HOMA" ~"Hordeum \nmarinum",
                                focal == "LEMA" ~"Leontodon \nmaroccanus")) %>%
  mutate(parameter = fct_relevel(parameter,"Plant - plant", "Intraspecific",
                                 "Plant - floral visitor","Plant - herbivore"),
         ParYear=paste0(parameter,"_",year)) %>%
  mutate(ParYearInt=case_when(ParYear == "Plant - plant_2021" ~17,
                              ParYear == "Plant - plant_2020"~16,
                              ParYear == "Plant - plant_2019"~15,
                              ParYear == "Intraspecific_2021"~12,
                              ParYear ==  "Intraspecific_2020"~11,
                              ParYear ==  "Intraspecific_2019"~10,
                              ParYear == "Plant - herbivore_2021"~7,
                              ParYear == "Plant - herbivore_2020"~ 6,
                              ParYear == "Plant - herbivore_2019"~ 5,
                              ParYear == "Plant - floral visitor_2021"~2,
                              ParYear == "Plant - floral visitor_2020"~1,
                              ParYear == "Plant - floral visitor_2019"~0)) %>%
  ggplot(aes(x=positive.perc,y=ParYearInt,label=year,
             shape=as.factor(focal.name))) +
  geom_vline(xintercept=0.5,color="black") +
  geom_point( color="black",
              aes(fill=positive.perc), 
              size=10,
              alpha=0.8, na.rm = TRUE) +
  geom_point(size=6,color="red",
             aes( x=1-neg.perc_hat,
                  fill=1-neg.perc_hat), 
             alpha=0.9, na.rm = TRUE) + 
  #geom_point(color="black",
  #          aes(y="Sum all", 
  #              x=positive.perc_mean,
  #             fill=positive.perc_mean), 
  #         alpha=0.9) + 
  geom_vline(xintercept=0.5,color="grey") + 
  scale_fill_gradientn(colours=c("#E69F00","white","#009E73"),
                       limits = c(0, 1),
                       guide = "none", aesthetics = "fill") +
  scale_shape_manual(values=c(21,22,23,24)) +
  scale_x_continuous("",
                     breaks=c(0,0.25,0.5,0.75,1),
                     labels = rev(c("100% + \n 0% -", 
                                    "75% + \n25% -",
                                    "50% + \n50% -",
                                    "25% + \n75% -",
                                    "0% + \n100% -"))) +
  scale_y_continuous(breaks =c(16,11,6,1),
                     minor_breaks =c(0,2,5,7,10,12,15,17),
                     labels=c("Plant - plant", "Intraspecific",
                              "Plant - herbivore", "Plant - floral visitor")) +
  theme_minimal() +
  labs(y="",x="",
       shape="Focal plant species",
       fill="Direction of effect") +
  guides(color= "none",fill="none",
         shape=guide_legend(override.aes = list(size = 10,color="black"),
                            nrow = 4,
                            direction="vertical",
                            byrow = TRUE,
                            title.hjust = 0.1)) +
  geom_text(x = 1.05, # Set the position of the text to always be at '14.25'
            hjust = 0,
            size = 7,
            family= "sans",
            fontface="plain") +
  coord_cartesian(xlim = c(0, 1), # This focuses the x-axis on the range of interest
                  clip = 'off') +   # This keeps the labels from disappearing
  theme(legend.position = c(1.28, 0.5),
        panel.grid.minor.x = element_blank(),
        panel.spacing.y=unit(0.5, "lines") ,
        #panel.spacing.x=unit(1.5,"lines"),
        legend.key.size = unit(1, 'cm'),
        axis.title =element_text(size=24),
        axis.text.x= element_text(size=22),
        axis.text.y= element_text(size=24),
        legend.text=element_text(size=22, face="italic"),
        legend.title=element_text(size=24),
        strip.text = element_text(size=20, face="italic"),
        panel.border = element_rect(color = "white", 
                                    fill = NA),
        plot.margin=grid::unit(c(1,100,35,1), "mm"))
plot_param_perc_small

img <- readPNG("~/Eco_Bayesian/Complexity_caracoles/figure/xaxisfig2.png")
g <- rasterGrob(img, interpolate=TRUE)

theme_set(theme_cowplot())
plot_param_perc_small_2 <- ggdraw() +
  draw_image(img,  x = -0.02, y = -0.41, scale = .52) +
  #draw_image(img,  x = 0.2125 , y = -0.4, scale = .35) +
  draw_plot(plot_param_perc_small)    
plot_param_perc_small_2 

ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/Plot_param_perc_fam.pdf"),
       #dpi="retina",
       #width = 25,
       #height = 16,
       #units = c("cm"),
       plot_param_perc_small_2 
)

#---- 4.3. Figure of species Strength ----

df_param_strength  <- NULL
for( focal in c("CETE","LEMA","HOMA","CHFU")){ # "CHFU","HOMA","CETE"
  for( year in c("2019",'2020','2021')){
    for( complexity.level in c(1:3)){
      complexity.animal <- c("group","family","species")[complexity.level]
      complexity.plant <-c("class","family","code.plant")[complexity.level]
      df_param_strength_n <- NULL
      df_param_strength_n <-  df_param_all[which(df_param_all$focal==focal  & 
                                                   df_param_all$year == year &
                                                   df_param_all$complexity.plant == complexity.plant),] %>%
        group_by(parameter) %>%
        mutate(mean.estimate = mean(exp(estimate)),
               median.estimate = median(exp(estimate)),
               ten.estimate = quantile(exp(estimate), probs = 0.12, na.rm = T),
               ninety.estimate = quantile(exp(estimate), probs = 0.88, na.rm = T),
               var.estimate=var(exp(estimate))) %>%
        mutate(mean.estimate_hat = mean(exp(estimate_hat), na.rm = T),
               median.estimate_hat = median(exp(estimate_hat), na.rm = T),
               ten.estimate_hat = quantile(exp(estimate_hat), probs = 0.12, na.rm = T),
               ninety.estimate_hat = quantile(exp(estimate_hat), probs = 0.88, na.rm = T),
               var.estimate_hat=var(exp(estimate_hat))) %>%
        ungroup() %>%
        mutate(mean.estimate.all = mean(exp(estimate)),
               median.estimate.all = median(exp(estimate)),
               ten.estimate.all = quantile(exp(estimate), probs = 0.12, na.rm = T),
               ninety.estimate.all = quantile(exp(estimate), probs = 0.88, na.rm = T),
               var.estimate.all =var(exp(estimate)))%>%
        select(focal,year,complexity.plant,parameter,
               mean.estimate,median.estimate,
               ten.estimate,ninety.estimate,var.estimate,
               mean.estimate_hat,median.estimate_hat,
               ten.estimate_hat,ninety.estimate_hat,
               var.estimate_hat,
               mean.estimate.all,median.estimate.all,
               ten.estimate.all,ninety.estimate.all,var.estimate.all) %>%
        unique()
      
      df_param_strength <- bind_rows(df_param_strength,df_param_strength_n)
      
    }
  }
}
str(df_param_strength)  
plot_param_strength <- df_param_strength %>%
  mutate(complexity.plant  = factor(complexity.plant,
                                    levels = c("class","family","code.plant")),
         complexity.plant  = case_when(complexity.plant=="class"~ "Order",
                                       complexity.plant=="family"~ "Family",
                                       complexity.plant=="code.plant"~ "Species"),
         focal.name = case_when(focal == "CHFU" ~"Chamaemelum \nfuscatum",
                                focal == "CETE" ~"Centaurium  \ntenuiflorum",
                                focal == "HOMA" ~"Hordeum \nmarinum",
                                focal == "LEMA" ~"Leontodon \nmaroccanus")) %>%
  mutate(parameter = fct_relevel(parameter,"Plant - plant", "Intraspecific",
                                 "Plant - floral visitor","Plant - herbivore"),
         complexity.plant = fct_relevel(complexity.plant,"Order",
                                        "Family","Species")) %>%
  ggplot(aes(x=mean.estimate,y=as.factor(parameter))) +
  geom_linerange(aes( xmin = mean.estimate - var.estimate , 
                      xmax = mean.estimate + var.estimate))+ 
  geom_linerange(aes(y="Sum all", xmin = mean.estimate.all - var.estimate.all , 
                     xmax = mean.estimate.all + var.estimate.all))+ 
  geom_point(aes(shape=as.factor(year),fill=mean.estimate),color="black",
             size=7) +
  geom_point(aes(x=mean.estimate_hat,
                 shape=as.factor(year),fill=mean.estimate_hat),color="red",
             size=5) +
  geom_point(aes(y="Sum all", x=mean.estimate.all,
                 shape=as.factor(year),
                 fill=mean.estimate.all),
             size=7,color="black") +
  geom_vline(xintercept=1,color="darkgrey") + 
  scale_fill_gradientn(colours=c("#E69F00","white","#009E73"),
                       limits = c(0.5, 1.5),
                       guide = "none", aesthetics = "fill") +
  scale_shape_manual("year",values=c(21,22,24)) +
  scale_x_continuous(breaks = c(0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5),
                     labels = c("-30%","-20%","-10%","No effect",
                                "+10%","+20%","+30%","+40%","+50%"),
                     limits = c(0.7,1.5)) +
  scale_y_discrete(limits=rev(c("Plant - plant", "Intraspecific",
                                "Plant - herbivore", "Plant - floral visitor","Sum all"))) +
  facet_grid( focal.name ~ complexity.plant) +
  theme_minimal() +
  labs(y="",x="Strength on fecundity") +
  guides(color= "none",fill="none",
         shape=guide_legend(override.aes = list(size = 6) ) ) +
  theme(panel.grid.minor = element_blank(),
        #panel.grid.major.x = element_blank(),
        panel.spacing.y=unit(0.5, "lines") ,
        panel.spacing.x=unit(2.5,"lines"),
        legend.key.size = unit(1, 'cm'),
        title =element_text(size=16),
        axis.text.x= element_text(size=14, angle=66, vjust=0.8),
        axis.text.y= element_text(size=16),
        legend.text=element_text(size=14, face="italic"),
        legend.title=element_text(size=16),
        strip.text = element_text(size=16, face="italic"),
        panel.border = element_rect(color = "black", 
                                    fill = NA),
        plot.margin=grid::unit(c(1,1,1,1), "mm"))
plot_param_strength

ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/plot_param_strength.pdf"),
       #dpi="retina",
       #width = 25,
       #height = 16,
       #units = c("cm"),
       plot_param_strength
)

# SMALL
df_param_strength$mean.estimate_hat
plot_param_strength <- df_param_strength %>%
  filter(complexity.plant=="family") %>%
  filter(!(focal =="HOMA" & parameter== "Plant - floral visitor")) %>%
  mutate(focal.name = case_when(focal == "CHFU" ~"Chamaemelum \nfuscatum",
                                focal == "CETE" ~"Centaurium \ntenuiflorum",
                                focal == "HOMA" ~"Hordeum \nmarinum",
                                focal == "LEMA" ~"Leontodon \nmaroccanus")) %>%
  mutate(parameter = fct_relevel(parameter,"Plant - plant", "Intraspecific",
                                 "Plant - floral visitor","Plant - herbivore"),
         ParYear=paste0(parameter,"_",year)) %>%
  mutate(ParYearInt=case_when(ParYear == "Plant - plant_2021" ~17,
                              ParYear == "Plant - plant_2020"~16,
                              ParYear == "Plant - plant_2019"~15,
                              ParYear == "Intraspecific_2021"~12,
                              ParYear ==  "Intraspecific_2020"~11,
                              ParYear ==  "Intraspecific_2019"~10,
                              ParYear == "Plant - herbivore_2021"~7,
                              ParYear == "Plant - herbivore_2020"~ 6,
                              ParYear == "Plant - herbivore_2019"~ 5,
                              ParYear == "Plant - floral visitor_2021"~2,
                              ParYear == "Plant - floral visitor_2020"~1,
                              ParYear == "Plant - floral visitor_2019"~0)) %>%
  ggplot(aes(x=mean.estimate,y=ParYearInt,label=year)) +
  geom_vline(xintercept=1,color="black") + 
  #geom_linerange(aes( xmin = mean.estimate - sqrt(var.estimate , 
  #                   xmax = mean.estimate + sqrt(var.estimate))+ 
  geom_point(aes(fill=mean.estimate,shape=focal.name),
             color="black",size=10,alpha=0.95) +
  geom_point(aes(x=mean.estimate_hat,fill=mean.estimate_hat,
                 shape=as.factor(focal.name)),na.rm=T,
             color="red",size=5,alpha=0.8) +
  
  scale_fill_gradientn(colours=c("#E69F00","white","#009E73"),
                       limits = c(0.5, 1.5),
                       guide = "none", aesthetics = "fill") +
  scale_shape_manual(values=c(21,22,23,24)) +
  scale_x_continuous(breaks = c(0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5),
                     labels = c("-30%","-20%","-10%","No\neffect",
                                "+10%","+20%","+30%","+40%","+50%"),
                     limits = c(0.7,1.5)) +
  scale_y_continuous(breaks =c(16,11,6,1),
                     minor_breaks =c(0,2,5,7,10,12,15,17),
                     labels=c("Plant - plant", "Intraspecific",
                              "Plant - herbivore", "Plant - floral visitor")) +
  theme_minimal() +
  labs(y="",x="",
       shape="Focal plant species",
       fill="Direction of effect") +
  guides(color= "none",fill="none",
         shape=guide_legend(override.aes = list(size = 10,color="black"),
                            nrow = 4,
                            direction="vertical",
                            byrow = TRUE,
                            title.hjust = 0.1)) +
  geom_text(x = 1.54, # Set the position of the text to always be at '14.25'
            hjust = 0,
            size = 7,
            family= "sans",
            fontface="plain") +
  coord_cartesian(xlim = c(0.7, 1.5), # This focuses the x-axis on the range of interest
                  clip = 'off') +   # This keeps the labels from disappearing
  theme(legend.position = c(1.25, 0.5),
        panel.grid.minor.x = element_blank(),
        panel.spacing.y=unit(0.5, "lines") ,
        #panel.spacing.x=unit(1.5,"lines"),
        legend.key.size = unit(1, 'cm'),
        axis.title =element_text(size=24),
        axis.text.x= element_text(size=22),
        axis.text.y= element_text(size=24),
        legend.text=element_text(size=22, face="italic"),
        legend.title=element_text(size=24),
        strip.text = element_text(size=20, face="italic"),
        panel.border = element_rect(color = "white", 
                                    fill = NA),
        plot.margin=grid::unit(c(1,100,35,1), "mm"))# This widens the right margin
plot_param_strength 
img <- readPNG("~/Eco_Bayesian/Complexity_caracoles/figure/xaxisfig3.png")
g <- rasterGrob(img, interpolate=TRUE)

theme_set(theme_cowplot())
plot_param_strength_2 <- ggdraw() +
  draw_image(img,  x = -0.025, y = -0.41, scale = .53) +
  #draw_image(img,  x = 0.2125 , y = -0.4, scale = .35) +
  draw_plot(plot_param_strength)    
plot_param_strength_2 

ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/plot_param_strength_small.pdf"),
       #dpi="retina",
       #width = 25,
       #height = 16,
       #units = c("cm"),
       plot_param_strength_2
)


var.estimate.plot <- test.full.df %>%
  filter(complexity.plant=="family") %>%
  ggplot(aes(  x = median.estimate, y = sqrt(var.estimate))) +
  geom_point(aes(shape=parameter, 
                 fill=positive.perc),size=7,
             alpha=0.8, stroke = 1.5) +
  geom_point(aes( y = sqrt(var.estimate_hat), x = sqrt(median.estimate_hat),
                  fill=positive.perc_hat),size=7,
             color="black",alpha=0.8, stroke = 1.5,shape=25) +
  scale_shape_manual("Type of\ngeneric interaction",values=c(21,22,23,24,25)) +
  scale_fill_gradientn("sign of interactions",
                       colours=c("#E69F00","white","#009E73"), limits = c(0, 1),
                       #midpoint = 0,
                       #space = "Lab",na.value = "grey50",
                       # guide = "none",
                       aesthetics = "fill",
                       breaks=c(0,0.5,1),
                       labels=c("100% Fac", "50/50", "100% Comp")) +
  labs(x="Median of estimated generic interactions",
       y="Standard deviation of estimated generic interactions") +
  scale_x_continuous(limits= c(-0.4,0.4)) +
  guides(#fill = guide_legend(override.aes = list(color="white")),
    shape = guide_legend(override.aes = list(size=5,alpha=1, stroke = 2))) +
  theme_minimal() +
  theme(legend.key.size = unit(1, 'cm'),
        axis.title.x= element_text(size=16, hjust =0.5),
        axis.title.y= element_text(size=16, hjust = 0.5),
        axis.text= element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        strip.text = element_text(size=16))
var.estimate.plot

ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/var.estimate.plot.pdf"),
       #dpi="retina",
       #width = 25,
       #height = 16,
       #units = c("cm"),
       var.estimate.plot 
)


#---- 4.3 Figure of specific interaction appearring ----
tesaurus.plant <- read.csv(paste0(home.dic,"data/plant_code.csv"),sep=",") %>%
  dplyr::select("family","code.plant","class")%>%
  mutate_all(.funs=tolower) %>%
  rename("species" = "code.plant" )

tesaurus.all <- read.csv(paste0(home.dic,"data/tesaurus_clean.csv")) %>%
  dplyr::select("group","family","species_id") %>%
  mutate_all(.funs=tolower) %>%
  rename("class" ="group" ,"species" = "species_id" ) %>%
  bind_rows(tesaurus.plant)

HTL.inclus <- read.csv(paste0(home.dic,"data/groupinglevels.csv"), sep=",") %>%
  dplyr::filter(YesNo =="yes") %>%
  dplyr::select("complexity","grouping","HTL") %>%
  rename("guild" = "HTL") %>%
  mutate_all(.funs=tolower) %>%
  mutate(complexity = case_when(complexity =="functional.groups" ~ "class",
                                TRUE ~ complexity))

species.tesaurus.inclus <- read.csv(paste0(home.dic,"data/plant_code.csv"),sep=",") %>%
  dplyr::select("family","code.plant","class") %>%
  gather(c(family,code.plant,class),key = "complexity",value="grouping" ) %>%
  mutate(complexity = case_when(complexity =="code.plant" ~ "species",
                                TRUE ~ complexity),
         guild = "plant" ) %>%
  mutate_all(.funs=tolower) %>%
  bind_rows(HTL.inclus) %>%
  unique() 

compl.levels <- levels(as.factor(species.tesaurus.inclus$complexity))
species.tesaurus.inclus.full <- NULL
for( n in c(1:length(compl.levels))){
  n.name <- compl.levels[n]
  species.tesaurus.inclus.n <- species.tesaurus.inclus %>%
    dplyr::filter(complexity ==  n.name) %>%
    dplyr::select(grouping,guild)
  if( n ==1){
    species.tesaurus.inclus.full <- species.tesaurus.inclus.n %>%
      mutate(complexity = n.name) %>%
      mutate(class = grouping)
  }else{
    names(species.tesaurus.inclus.n) <- c( n.name,"guild")
    species.tesaurus.inclus.n <- left_join(species.tesaurus.inclus.n,
                                           unique(tesaurus[,compl.levels[1:n]]),
                                           multiple = "all") %>%
      mutate(complexity = n.name)
    names(species.tesaurus.inclus.n)[1] <-c("grouping")
    species.tesaurus.inclus.full <- bind_rows(species.tesaurus.inclus.full, species.tesaurus.inclus.n)
  }
}

view(species.tesaurus.inclus.full) 

names(species.tesaurus.inclus.full)[1] <-c("interactor")
plant.species <- species.tesaurus.inclus.full$interactor[which(species.tesaurus.inclus.full$guild=="plant")]
FV.species <- species.tesaurus.inclus.full$interactor[which(species.tesaurus.inclus.full$guild=="floral.visitor")]
H.species <- species.tesaurus.inclus.full$interactor[which(species.tesaurus.inclus.full$guild=="herbivore")]


df_inclusion_nat_species <- df_inclusion_nat_vertical.2[,c("complexity_plant","interactor","n.number","type.interaction")] %>%
  dplyr::filter(!is.na(interactor)) %>%
  mutate_all(.funs=tolower) %>%
  mutate(interactor.length = sapply(gregexpr(",",interactor),length)) %>% 
  mutate(interactor.length = case_when(interactor.length > 1 ~interactor.length +1, 
                                       T ~ interactor.length),
         n.number = as.numeric(n.number)/as.numeric(interactor.length) , 
         interactor = strsplit(as.character(interactor), ",")) %>% 
  unnest(interactor) %>% 
  mutate(n.number = 1,
         guild= case_when(interactor %in% tolower(plant.species) ~ "plant",
                          interactor %in% tolower(FV.species) ~ "floral.visitor",
                          interactor %in% tolower(H.species) ~ "herbivore"),
         complexity = case_when(complexity_plant =="code.plant" ~ "species",
                                T ~ complexity_plant),
         n.number  = as.numeric(n.number)) %>%
  aggregate(n.number ~ interactor +  complexity + guild, sum) %>%
  unique() %>%
  full_join(species.tesaurus.inclus.full, by = join_by(interactor,guild,complexity)) %>%
  arrange(n.number)


potential.interactor <- df_inclusion_nat_vertical.2[which(df_inclusion_nat_vertical.2$type.interaction == "Direct interactions"),
                                                    c("potential.interactor")]  %>%
  unique()
potential.interactor <- tolower(unique(str_split_1(paste(potential.interactor, collapse = ","),",")))
potential.interactor <- potential.interactor[which(!potential.interactor =="na")]


df_inclusion_nat_species <- df_inclusion_nat_species[which(!is.na(df_inclusion_nat_species$class)),]
df_inclusion_nat_species <- df_inclusion_nat_species  %>%
  mutate(n.number.bis = case_when(interactor %in% potential.interactor ~ 0,
                                  T~ NA),
         n.number = case_when(n.number.bis == 0 & is.na(n.number) ~ n.number.bis,
                              T ~ n.number)) %>%
  dplyr::filter(!is.na(n.number)) %>%
  mutate(interactor = as.character(interactor),
         class =as.character(class),
         guild = factor(guild, levels=c('plant','herbivore','floral.visitor'))) %>%
  arrange(desc(n.number))

view(df_inclusion_nat_species) 


plot_inclusion_species_order <-  ggplot(df_inclusion_nat_species[which(df_inclusion_nat_species$guild=="plant"),], 
                                        aes(x=fct_inorder(interactor), y = n.number, fill=fct_inorder(class))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c(safe_colorblind_palette,"blue","orange")) +
  facet_wrap(complexity ~ ., nrow=1,
             scales= "free") +
  theme_bw() +
  labs(y="number of relevant \nspecies-specific interaction retained",
       x="interacting species",
       fill="functional group of species") +
  theme(axis.text.x = element_text(angle = 66, hjust=1,size=16),
        legend.position = "top",legend.key.size = unit(1, 'cm'),
        title =element_text(size=16),
        axis.text.y= element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        strip.text = element_text(size=16))

plot_inclusion_species_order  
ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/plot_inclusion_species_order.pdf"),
       plot_inclusion_species_order )



var.estimate.plot <- ggplot(test.full.df, 
                            aes(  x = sqrt(var.abundance), y = sqrt(var.estimate))) +
  #geom_smooth(color="black",alpha=0.1) +
  geom_point(aes(color=parameter, 
                 fill=neg.perc,size=abs(median.estimate)),
             shape=21,alpha=0.8, stroke = 1.5) +
  geom_point(aes( y = sqrt(var.estimate_hat), x = sqrt(var.abundance_hat),
                  fill=neg.perc_hat, size=abs(median.estimate_hat)),
             shape=21, color="#661100",alpha=0.8, stroke = 1.5) +
  scale_size_continuous("strenght of interactions\n(absolute value)",
                        range = c(1, 8),
                        breaks=c(0.1,0.2,0.3),
                        labels = c("< 0.1", "[0.1,0.2]", "> 0.2")) +
  scale_fill_gradient("Sign of interaction", 
                      low = "white",high = "black",
                      breaks=c(0,0.5,1),
                      labels=c("100% Fac", "50/50", "100% Comp")) +
  scale_color_manual(values=c("#DDCC77","#88CCEE","#999933","#117733","#661100")) + 
  labs(x="Standard deviation in the observed abundance of neighbours",
       y="Standard deviation in estimation of interaction") +
  #geom_tile(aes(fill=abs(median.estimate))) +
  #geom_errorbarh(aes(xmin=abs(ten.estimate),xmax=abs(ninety.estimate)  ))+
  guides(#fill = guide_legend(override.aes = list(color="white")),
    colour = guide_legend(override.aes = list(size=5,alpha=1, stroke = 2))) +
  theme_minimal() +
  theme(legend.key.size = unit(1, 'cm'),
        axis.title.x= element_text(size=16, hjust =0.5),
        axis.title.y= element_text(size=16, hjust = 0.5),
        axis.text= element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        strip.text = element_text(size=16))

var.estimate.plot 


#---- 4.4. Figure of gradient neighbors ----
herbivorypredator <- read.csv(paste0(home.dic,"data/herbivorypredator.csv"), sep=",")

floral_visitor <- read.csv(paste0(home.dic,"data/floral_visitor.csv"), sep=",")

groupinglevels <- read.csv(paste0(home.dic,"data/groupinglevels.csv"), sep=",")


HtoKeep_group <- groupinglevels$grouping[which(groupinglevels$HTL =="herbivore" &
                                                 groupinglevels$complexity == "functional.groups" & groupinglevels$YesNo == "yes")]

herbivorabundance <- herbivorypredator %>%
  filter(plant %in% c("HOMA","CHFU","LEMA","CETE")) %>%
  filter(year %in% c("2019","2020","2021")) %>%
  dplyr::filter(group %in% HtoKeep_group) %>%
  aggregate(number_animal ~ subplot + plant + year + plot, sum) %>%
  rename("abundance"=number_animal,
         "focal" = plant ) %>%
  mutate(neighbours = "herbivore")


FVtoKeep_group <- groupinglevels$grouping[which(groupinglevels$HTL =="floral.visitor" &
                                                  groupinglevels$complexity == "functional.groups" & groupinglevels$YesNo == "yes")]

HTL_abundance <- floral_visitor %>%
  filter(plant %in% c("HOMA","CHFU","LEMA","CETE")) %>%
  filter(year %in% c("2019","2020","2021")) %>%
  dplyr::filter(group %in% FVtoKeep_group) %>%
  aggregate(number_visits ~ subplot + plant + year + plot, sum)  %>%
  mutate(neighbours = "floral visitor") %>%
  rename("abundance"=number_visits,
         "focal" = plant) %>%
  bind_rows(herbivorabundance)


SpData <- read.csv(paste0(home.dic,"data/SpData.csv")) 
plant.neigh.df <- SpData %>%
  filter(focal %in% c("CETE","HOMA","CHFU","LEMA"))%>%
  filter(year %in% c("2019","2020","2021"))
plant.neigh.df$species.specific <- NA
for(n in 1:nrow(plant.neigh.df)){
  print(n)
  if(plant.neigh.df$year[n] %in% df_inclusion_nat_short$year[which(df_inclusion_nat_short$focal == plant.neigh.df$focal[n])]){
    interactor <- df_inclusion_nat_vertical$interactor[which(df_inclusion_nat_vertical$year ==plant.neigh.df$year[n] &
                                                               df_inclusion_nat_vertical$focal == plant.neigh.df$focal[n] &
                                                               df_inclusion_nat_vertical$interaction =="Plant - plant")]
    interactor <-  interactor[!is.na( interactor)] %>%
      paste(collapse=",") %>%
      stringr::str_split_1(",")
    
    # keep the highest level of 
    interactor <- plant.class[which(plant.class$class %in% interactor |
                                      plant.class$family %in% interactor |
                                      plant.class$code.plant %in% interactor ),"code.plant"]
    
    
    plant.neigh.df[n,] <- plant.neigh.df[n,] %>%
      mutate(species.specific = select(., all_of(interactor)) %>% 
               rowSums(na.rm = TRUE))
  }else{
    next
  }
}

plant.neigh.density <- plant.neigh.df %>%
  filter(focal %in% c("CETE","HOMA","CHFU","LEMA"))%>%
  filter(year %in% c("2019","2020","2021")) %>%
  mutate(Intraspecific = case_when(focal =="CETE" ~ CETE,
                                   focal =="HOMA" ~ HOMA,
                                   focal =="LEMA" ~ LEMA,
                                   focal =="CHFU" ~ CHFU)) %>%
  mutate(Interspecific = select(., all_of(complexitylevel.plant.class)) %>% 
           rowSums(na.rm = TRUE)) %>%
  gather(c(Intraspecific,Interspecific,species.specific), key="neighbours", value="abundance") %>%
  mutate(neighbours = case_when(neighbours =="Intraspecific" ~ "intra- plant",
                                neighbours =="Interspecific" ~ "inter- plant",
                                neighbours =="species.specific" ~ "plant(s) with\nspecific effect")) %>%
  dplyr::select(all_of(names(HTL_abundance))) %>%
  bind_rows(HTL_abundance) %>%
  filter(!is.na(abundance)) %>%
  ggplot(aes(x=abundance, y=neighbours, fill=neighbours)) +
  geom_density_ridges( #stat = "binline", bins = 20, scale = 0.95, draw_baseline = FALSE,
    scale = 1,
    jittered_points = TRUE,
    position = position_points_jitter(width = 0.005, height = 0),
    point_shape = '|', point_size = 1, point_alpha = 1, alpha = 0.7) +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE)+
  geom_vline(xintercept = 0, color="black", linetype="dashed")+
  scale_fill_manual(values=c("#DDCC77","#88CCEE","#999933","#117733", "#661100"))+ 
  guides(fill="none") +
  scale_x_continuous(limits=c(-5,60)) +
  labs(x="Total number of individual in the neighbourhood",
       y="type of neighbours",
       title="Frequency of abundances observed in the neighboorhood") +
  facet_grid(focal~year, scales="free") +
  theme(legend.position = "top",legend.key.size = unit(1, 'cm'),
        title =element_text(size=16),
        axis.text.y= element_text(size=12),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        strip.text = element_text(size=16))

plant.neigh.density

ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/plant.neigh.density.pdf"),
       dpi="retina",
       #width = 21,
       #height = 16,
       #units = c("cm"),
       plant.neigh.density)

#---- 4.5. Figure of gradient neighbors related to strength ----
head(df_param_perc)
test.df <- plant.neigh.df %>%
  filter(focal %in% c("CETE","HOMA","CHFU","LEMA"))%>%
  filter(year %in% c("2019","2020","2021")) %>%
  mutate(Intraspecific = case_when(focal =="CETE" ~ CETE,
                                   focal =="HOMA" ~ HOMA,
                                   focal =="LEMA" ~ LEMA,
                                   focal =="CHFU" ~ CHFU)) %>%
  mutate(Interspecific = select(., all_of(complexitylevel.plant.class)) %>% 
           rowSums(na.rm = TRUE)) %>%
  gather(c(Intraspecific,Interspecific,species.specific), key="neighbours", value="abundance") %>%
  mutate(neighbours = case_when(neighbours =="Intraspecific" ~ "Intraspecific",
                                neighbours =="Interspecific" ~ "Plant - plant",
                                neighbours =="species.specific" ~ "Interspecific_hat")) %>%
  dplyr::select(all_of(names(HTL_abundance))) %>%
  bind_rows(HTL_abundance) %>%
  filter(!is.na(abundance))  %>%
  group_by(cbind(neighbours, focal, year)) %>%
  mutate(mean.abundance = mean(abundance),
         median.abundance = median(abundance),
         ten.abundance = quantile(abundance, probs = 0.12, na.rm = T),
         ninety.abundance = quantile(abundance, probs = 0.88, na.rm = T),
         var.abundance = var(abundance)) %>%
  ungroup() %>%
  select(year,focal,neighbours,mean.abundance,median.abundance,ten.abundance,ninety.abundance,var.abundance) %>%
  unique()  %>%
  rename("parameter" = neighbours)

view(test.df)

test.df <- bind_rows(test.df[which(test.df$parameter != "Plant - plant" &
                                     test.df$parameter !="Interspecific_hat"),],
                     dplyr::full_join(test.df[which(test.df$parameter == "Plant - plant" ),],
                                      test.df[which(test.df$parameter == "Interspecific_hat"),], 
                                      by= c("year", "focal"),
                                      suffix= c("","_hat"))
) %>%
  mutate(parameter = case_when(parameter == "floral visitor" ~ "Plant - floral visitor",
                               parameter == "herbivore" ~ "Plant - herbivore",
                               T~parameter))

test.full.df <- dplyr::full_join(df_param_perc,
                                 test.df , 
                                 by= c("year", "focal","parameter"),
                                 multiple = "all",
                                 suffix= c("",""))


test.full.df[nrow(test.full.df)+1,] <- NA
test.full.df[nrow(test.full.df)+1,"parameter"] <- "Species-specific"
test.full.df <- test.full.df[-which(is.na(test.full.df$parameter)),] 


var.abundance.estimate.plot <- ggplot(test.full.df, 
                                      aes(  x = sqrt(var.abundance)/mean.abundance, y = sqrt(var.estimate))) +
  #geom_smooth(color="black",alpha=0.1) +
  geom_point(aes(shape=parameter, 
                 fill=positive.perc,
                 size=abs(mean.estimate)),
             alpha=0.8, stroke = 1.5) +
  geom_point(aes( y = sqrt(var.estimate_hat), x = sqrt(var.abundance_hat)/mean.abundance_hat,
                  fill=positive.perc_hat, size=abs(median.estimate_hat)),
             color="red",alpha=0.8, stroke = 1.5,shape=25) +
  scale_shape_manual(values=c(21,22,23,24,25)) +
  scale_size_continuous("strenght of interactions\n(absolute value)",
                        range = c(1, 8),
                        breaks=c(0.1,0.2,0.3),
                        labels = c("< 0.1", "[0.1,0.2]", "> 0.2")) +
  scale_fill_gradientn("sign of interactions",
                       colours=c("#E69F00","white","#009E73"), limits = c(0, 1),
                       #midpoint = 0,
                       #space = "Lab",na.value = "grey50",
                       # guide = "none",
                       aesthetics = "fill",
                       breaks=c(0,0.5,1),
                       labels=rev(c("100% Fac", "50/50", "100% Comp"))) +
  #scale_color_manual(values=c("#DDCC77","#88CCEE","#999933","#117733","#661100")) + 
  labs(x="Standard deviation divided by mean in the observed abundance of neighbours",
       y="Standard deviation in estimation of interaction") +
  #geom_tile(aes(fill=abs(median.estimate))) +
  #geom_errorbarh(aes(xmin=abs(ten.estimate),xmax=abs(ninety.estimate)  ))+
  guides(#fill = guide_legend(override.aes = list(color="white")),
    size = guide_legend(override.aes = list(shape=1, color="black")),
    shape = guide_legend(override.aes = list(size=5,alpha=1, stroke = 2,
                                             color=c("black","black","black","black","red")))) +
  theme_minimal() +
  theme(legend.key.size = unit(1, 'cm'),
        axis.title.x= element_text(size=16, hjust =0.5),
        axis.title.y= element_text(size=16, hjust = 0.5),
        axis.text= element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        strip.text = element_text(size=16))
var.abundance.estimate.plot

ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/var.abundance.estimate.plot.pdf"),
       dpi="retina",
       width = 30,
       height = 21,
       units = c("cm"),
       var.abundance.estimate.plot)


#---- 5.0. Potential figure for manuscripts ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 5.1. Variation intrinsic parameters ----

plot_lambda_violin_std_all <- ggplot(df_param,
                                     aes(y=scale(lambdas), x=as.factor(year),
                                         fill=as.factor(focal))) +
  geom_violin(position = position_dodge(width = 0.3),
              width=3) + 
  theme_bw() + facet_grid(.~complexity.plant) + 
  labs(y="Intrinsic fecundity (scaled)", x="year") 


plot_alpha_intra_violin_all <- ggplot(df_param,
                                      aes(y=alpha_intra, x=as.factor(year),
                                          fill=as.factor(focal)))+
  geom_hline(yintercept = 0, color="grey") +
  geom_violin(position = position_dodge(width = 0.3),
              width=3) + ylim(-1,1) +
  theme_bw() + 
  labs(y="Intra-specific interactions", x="year") +
  facet_grid(.~complexity.plant)


plot_alpha_intra_density <- ggplot(df_param[which(df_param$complexity.plant=="family"),],
                                   aes(y=alpha_intra,fill=as.factor(focal)))+
  geom_density(alpha=0.5) + ylim(-2,2) +
  theme_bw() + labs(y="",x="distributions")

plot_alpha_violin <- ggplot(df_param,
                            aes(y=alpha_generic, x=as.factor(year),
                                fill=as.factor(focal)))+
  geom_hline(yintercept = 0, color="grey") +
  geom_violin(position = position_dodge(width = 0.3),
              width=3) +  ylim(-1,1) +
  theme_bw() +
  labs(y="Generic Inter-specific interactions", x="year") 

plot_alpha_density <- ggplot(df_param[which(df_param$complexity.plant=="family"),],
                             aes(y=alpha,fill=as.factor(focal)))+
  geom_density(alpha=0.5) + ylim(-1,1) +
  theme_bw() + labs(y="",x="distributions")

ggsave( "figures/LambdavsIntravsInter.pdf", 
        ggarrange(plot_lambda_violin_std_all,#plot_lambda_density,
                  plot_alpha_intra_violin_all,#plot_alpha_intra_density,
                  plot_alpha_violin,#plot_alpha_density,
                  ncol=1,nrow=3, common.legend = T,
                  widths=(c(2,1)))
) 
#---- 5.2. Fit of model to observations ----
#---- 5.3. Correlation between alpha_intra and lambda ----
df.corr <- NULL
for( focal in c("CETE","LEMA","CHFU","HOMA")){ # "CHFU","HOMA",
  for( year in c("2019","2020","2021")){
    res <- cor.test(df_param$alpha_intra[which(df_param$focal==focal &
                                                 df_param$year==year)],
                    df_param$lambdas[which(df_param$focal==focal & 
                                             df_param$year==year)],
                    method=c("pearson"))
    if(res$p.value < 0.05){
      significant = "yes"
    }else{significant = "no"}
    
    df.corr <- bind_rows(  df.corr,
                           data.frame( focal = focal, 
                                       year = year,
                                       p.value =    res$p.value, 
                                       corr =   res$estimate,
                                       significant = significant))
  }
}
write.csv(df.corr, file="~/Eco_Bayesian/Complexity_caracoles/results/CORRIntravsLambda")



#---- 6.0. Appendix figures ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
home.dic <- ""
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/"
#---- 6.1. Neighbourhood abundances ----
SpData <- read.csv(paste0(home.dic,"data/SpData.csv"))
plant.class <- read.csv(paste0(home.dic,"data/plant_code.csv"),sep=",")
complexity.plant <-c("class","family","code.plant")[2]

complexitylevel.plant <- levels(as.factor(plant.class[,complexity.plant]))


Neighbourhood_plot <- SpData %>%
  filter(focal %in% c("HOMA","CHFU","LEMA","CETE")) %>%
  filter(year %in% c("2019","2020","2021")) %>%
  dplyr::select(year, plot, focal,complexitylevel.plant) %>%
  gather(complexitylevel.plant, key=family, value=observation) %>%
  aggregate(observation ~ year +plot + focal + family, length) %>%
  ggplot(aes(x=as.factor(plot),y=observation,fill=family )) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values=safe_colorblind_palette) +
  facet_grid(focal~year, scale="free_y") +
  labs(title="Neigbourhood observations",
       y = "number of times a family is observed in the neighbours of the focal", x="plot") + 
  theme_bw()
Neighbourhood_plot

ggsave(paste0("figures/Neighbourhood_plot.pdf"),
       Neighbourhood_plot)

#---- 6.2. HTL abundances ----
#rm(list = ls())

herbivorypredator <- read.csv(paste0(home.dic,"data/herbivorypredator.csv"), sep=",")

floral_visitor <- read.csv(paste0(home.dic,"data/floral_visitor.csv"), sep=",")

groupinglevels <- read.csv(paste0(home.dic,"data/groupinglevels.csv"), sep=",")

head(floral_visitor)

HtoKeep_group <- groupinglevels$grouping[which(groupinglevels$HTL =="herbivore" &
                                                 groupinglevels$complexity == "functional.groups" & groupinglevels$YesNo == "yes")]

herbivorgroup_plot <- herbivorypredator %>%
  filter(plant %in% c("HOMA","CHFU","LEMA","CETE")) %>%
  filter(year %in% c("2019","2020","2021")) %>%
  aggregate(number_animal ~ group + plant + year + plot, length) %>%
  dplyr::filter(group %in% HtoKeep_group) %>%
  rename("numberobservations"=number_animal ) %>%
  ggplot(aes(x=as.factor(plot),y=numberobservations,fill=group )) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values=safe_colorblind_palette) +
  facet_grid(plant~year, scale="free_y") +
  labs(title="Herbivore observations", fill="Herbivore\nfunctional group",
       y = "number of observations", x="plot") + 
  theme_bw()

ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/Herb_plot.pdf"),
       dpi="retina",
       #width = 21,
       #height = 16,
       #units = c("cm"),
       herbivorgroup_plot)


FVtoKeep_group <- groupinglevels$grouping[which(groupinglevels$HTL =="floral.visitor" &
                                                  groupinglevels$complexity == "functional.groups" & groupinglevels$YesNo == "yes")]

HTL_plot <- floral_visitor %>%
  filter(plant %in% c("HOMA","CHFU","LEMA","CETE")) %>%
  filter(year %in% c("2019","2020","2021")) %>%
  dplyr::filter(group %in% FVtoKeep_group) %>%
  aggregate(number_visits ~ group + plant + year + plot, length) %>%
  rename("numberobservations"=number_visits ) %>%
  bind_rows(herbivorgroup) %>%
  ggplot(aes(x=as.factor(plot),y=numberobservations,fill=group )) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values=safe_colorblind_palette) +
  facet_grid(plant~year, scale="free_y") +
  labs(title="Floral visitors and herbivore observations",
       y = "number of observations", x="plot") + 
  theme_bw()

Poll_plot <- floral_visitor %>%
  filter(plant %in% c("HOMA","CHFU","LEMA","CETE")) %>%
  filter(year %in% c("2019","2020","2021")) %>%
  dplyr::filter(group %in% FVtoKeep_group) %>%
  aggregate(number_visits ~ group + plant + year + plot , length) %>%
  rename("numberobservations"=number_visits ) %>%
  ggplot(aes(x=as.factor(plot),y=numberobservations,fill=group )) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values=safe_colorblind_palette) +
  facet_grid(plant~year, scale="free_y") +
  labs(title="Floral visitors observations",
       fill="Floral visitor\nfunctional group",
       y = "number of visits", x="plot") + 
  theme_bw()


ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/Poll_plot.pdf"),
       dpi="retina",
       #width = 21,
       #height = 16,
       #units = c("cm"),
       Poll_plot)


