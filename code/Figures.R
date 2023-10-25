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
#
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
#---- 3.0 CSV for inclusion----

df_inclusion_nat <- data.frame()
for( focal in c("LEMA","CHFU","HOMA","CETE")){ # "CHFU","HOMA",
for( year in c("2019",'2020','2021')){
  for( complexity.level in c(1:3)){
    if(complexity.level ==3 & focal=="HOMA") next
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
      if(complexity.level ==3 & focal=="HOMA") next

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


plot.alphas <- ggplot(df_param_all[which(df_param_all$complexity.plant=="family"),]) +  
  geom_density_ridges_gradient(aes(x=estimate, y=parameter,
                                   fill= after_stat(x)),
                               scale = 1) + 
  geom_boxplot(aes(x=estimate + estimate_hat, y = parameter,
                   color=parameter_hat),
               width=0.3,alpha=0.4,outlier.shape = NA) +
  facet_wrap(as.factor(focal)~as.factor(year), 
             ncol = 3,nrow=4) +
  #scale_x_continuous(limits=c(-0.5,0.5)) + 
  #scale_y_discrete(labels=c("Generic Inter-specific",
  #                          "Intra-specific")) + 
  labs(y="Interactions", x= "Estimated distribution",
       color = "species-specific \n parameters") + 
  scale_color_manual(values=rep(c("black"),times=9)) +
  scale_fill_gradientn(
    colours=c(re,"white",gr),limits = c(-2, 2),
    #midpoint = 0,
    #space = "Lab",
    na.value = "grey50",
    guide = "none",
    aesthetics = "fill"
  ) +
  theme_bw() 

plot.alphas
ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/Parameters_distribution_species.pdf"),
       dpi="retina",
       width = 21,
       height = 16,
       units = c("cm"),
       plot.alphas
)


#---- 4.2 Figure of pairwise interactions sign percentage ----
  df_param_perc  <- NULL
  for( focal in c("CETE","LEMA","HOMA","CHFU")){ # "CHFU","HOMA","CETE"
    for( year in c("2019",'2020','2021')){
      for( complexity.level in c(1:3)){
        if(complexity.level ==3 & focal=="HOMA") next
        
        complexity.animal <- c("group","family","species")[complexity.level]
        complexity.plant <-c("class","family","code.plant")[complexity.level]
        df_param_perc_n  <- NULL
        df_param_perc_n <-  df_param_all[which(df_param_all$focal==focal  & 
                                            df_param_all$year == year &
                                              df_param_all$complexity.plant == complexity.plant),] %>%
          group_by(parameter) %>%
          mutate(positive.perc = length(which(estimate > 0))/length(estimate),
                 neg.perc = length(which(estimate < 0))/length(estimate),
                 null.perc = length(which(estimate==0))/length(estimate))%>%
          mutate(positive.perc_hat = length(which(estimate_hat > 0))/length(estimate_hat),
                 neg.perc_hat = length(which(estimate_hat < 0))/length(estimate_hat),
                 null.perc_hat = length(which(estimate_hat==0))/length(estimate_hat)) %>%
          ungroup() %>%
          mutate(positive.perc_mean = mean(positive.perc),
                 positive.perc_var = var(positive.perc),
                 neg.perc_mean = mean(neg.perc),
                 neg.perc_var = var(neg.perc),) %>%
          select(focal,year,complexity.plant,parameter,
                 positive.perc, neg.perc,null.perc,
                 positive.perc_hat,neg.perc_hat, null.perc_hat,
                 positive.perc_mean, positive.perc_var,
                 neg.perc_mean, neg.perc_var) %>%
          unique()
        
        df_param_perc <- bind_rows(df_param_perc,df_param_perc_n)

      }
    }
  }
str(df_param_perc)
view(df_param_perc)
df_param_perc$neg.perc_hat[which(df_param_perc$neg.perc_hat==0 & 
                                   df_param_perc$positive.perc_hat==0)] <- NA
df_param_perc$complexity.plant <- factor(df_param_perc$complexity.plant,
                                  levels = c("class","family","code.plant"))

df_param_perc$parameter <- factor(df_param_perc$parameter,
                                                levels = c("Plant - floral visitor","Plant - herbivore",
                                                           "Plant - plant", "Intraspecific"))

levels(df_param_perc$parameter)
plot_param_perc <- ggplot(df_param_perc, aes(y=parameter,
                       x=neg.perc, fill=parameter)) + 
  geom_rect(xmin=0.4,xmax=0.6,ymin=-Inf,ymax=Inf,
            fill="grey", color="grey", size=0.5, alpha=0.02) + 
  geom_point(size=4, color="black",
             aes(shape=as.factor(year),fill=parameter), alpha=0.8) + 
  geom_point(size=2,color="red",
             aes(y=parameter, x=neg.perc_hat,fill=parameter, shape=as.factor(year)), 
             alpha=0.9) + 
  geom_point(size=4,color="black",fill="black",
             aes(y="ALL", x=neg.perc_mean,shape=as.factor(year)), 
             alpha=0.9) + 
  geom_linerange(aes(y="ALL", xmin = neg.perc_mean - neg.perc_var^2 , 
                     xmax = neg.perc_mean +neg.perc_var^2))+
  facet_grid(focal ~complexity.plant)  +
  geom_vline(xintercept=0.5,color="grey") + 
  scale_fill_manual(values=rev(c("#332288","#117733","#DDCC77","#88CCEE"))) +
  scale_shape_manual(values=c(21,22,24)) +
  scale_x_continuous("ratio Facilitation/Competition",
                     breaks=c(0,0.25,0.5,0.75,1),
                     labels = c("100% F \n 0% C", 
                                "75% F \n25% C",
                                "50% F \n50% C",
                                "25% F \n75% C",
                                "0% F \n100% C")) +
  theme_bw() +
  labs(title = "Percentage of Facilitation vs Competition in Direct interactions",
       y="",
       shape="year",
       fill="interaction") +
  guides(color= "none",fill="none") +
  theme(panel.grid.minor = element_blank(),
        panel.spacing.y=unit(0.5, "lines") ,
        panel.spacing.x=unit(1.5,"lines"))
plot_param_perc


ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/Plot_param_perc.pdf"),
       dpi="retina",
       width = 25,
       height = 16,
       units = c("cm"),
       plot_param_perc
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
  theme(axis.text.x = element_text(angle = 66, hjust=1),
        legend.position = "top")
  
plot_inclusion_species_order  
ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/plot_inclusion_species_order.pdf"),
        plot_inclusion_species_order )


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

#rm(list = ls())
abundance.plant <- read.csv(paste0(home.dic,"data/abundance.csv"), sep=",")

competition.plant <- read.csv(paste0(home.dic,"data/competition.csv"), sep=",")

plant.class <- read.csv(paste0(home.dic,"data/plant_code.csv"),sep=",")

herbivorypredator <- read.csv(paste0(home.dic,"data/herbivorypredator.csv"), sep=",")

floral_visitor <- read.csv(paste0(home.dic,"data/floral_visitor.csv"), sep=",")

groupinglevels <- read.csv(paste0(home.dic,"data/groupinglevels.csv"), sep=",")

head(floral_visitor)

HtoKeep_group <- groupinglevels$grouping[which(groupinglevels$HTL =="herbivore" &
                                                 groupinglevels$complexity == "functional.groups" & groupinglevels$YesNo == "yes")]

herbivorgroup <- herbivorypredator %>%
  filter(plant %in% c("HOMA","CHFU","LEMA","CETE")) %>%
  filter(year %in% c("2019","2020","2021")) %>%
  aggregate(number_animal ~ group + plant + year + plot, length) %>%
  dplyr::filter(group %in% HtoKeep_group) %>%
  rename("numberobservations"=number_animal )
  

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
  

ggsave(paste0("figures/HTL_plot.pdf"),
       HTL_plot )


