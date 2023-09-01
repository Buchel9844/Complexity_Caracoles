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
#remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)

home.dic <- "/home/lbuche/Eco_Bayesian/Complexity_caracoles/"
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/"


home.dic <- ""
project.dic <- ""
load("results/inclusion/InclusionLEMA_2018_class_group.RData")
assign(paste0("InclusionLEMA_2018_class_group"),Inclusion_all)

view(bind_rows(InclusionLEMA_2021_class_group$interaction,
               InclusionLEMA_2020_class_group$interaction,
               InclusionLEMA_2018_class_group$interaction))

#---- 3. Natural data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 3.0 CSV for inclusion----

df_inclusion_nat <- data.frame()
for( focal in c("CETE","LEMA","HOMA","CHFU")){ # "CHFU","HOMA",
for( year in c("2018","2019",'2020','2021')){
  for( complexity.level in 1){
    if(year == "2018" & focal =="CHFU") next
    if(complexity.level == 2 & focal =="HOMA"  & year=="2020") next
    if((year == "2018"|year == "2019"|year == "2021") & focal =="HOMA") next

    complexity.animal <- c("group","family","species")[complexity.level]
    complexity.plant <-c("class","family","code.plant")[complexity.level]
    
    if(year == "2018" & focal == "CHFU") next
    load(paste0(home.dic,"results/inclusion/Inclusion",focal,"_",
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

summary.interactions <- read.csv(paste0(home.dic,"results/Chapt1_Inclusion_parameters.csv"))
write.csv(df_inclusion_nat,
           file = paste0(home.dic,"results/Chapt1_Inclusion_parameters.csv"))
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
  scale_pattern_manual(values = c("Direct interactions" = "none", "HOIs"= "stripe")) + 
  facet_grid(focal~ complexity_plant,scale="free") +
  theme_bw() + #theme(axis.text.x = element_text(angle=45)) + 
  scale_fill_manual(values=c("#117733","#DDCC77","#88CCEE","#117733","#DDCC77","#88CCEE")) + 
  labs(x="year",y="Number of species specific interactions",
       fill =" interaction") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = c("none","none","none",
                                                             "stripe","stripe","stripe")
                                                 ))) #+ 
 # geom_text(aes(label=interactor, y=position.label.2),
          #      size=df_inclusion_nat_vertical.2$size.label, angle=df_inclusion_nat_vertical.2$angle.label)

df_inclusion_nat_bar 
ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/df_inclusion_plot_bar.pdf"),
       dpi="retina",
       width = 21,
       height = 16,
       units = c("cm"),
       df_inclusion_nat_bar)

#---- 3.1  Figure of pairwise interactions distributions and Intra Vs Inter ----
#extraction generic parameters
df_param <- NULL


  for( focal in c("CETE","LEMA","HOMA","CHFU")){ # "CHFU","HOMA","CETE"
    for( year in c("2018","2019",'2020','2021')){
      for( complexity.level in 1){
          if(year == "2018" & focal =="CHFU") next
          if(complexity.level == 2 & focal =="HOMA"  & year=="2020") next
          
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
for ( i in 1:nrow(df_inclusion_nat_short)){
  year <- df_inclusion_nat_short[i,'year']
  focal <- df_inclusion_nat_short[i,'focal']
  if(focal == "CETE") next
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
  
      if(is.null(df_param_hat)){
      df_param_hat <- df_param_hat_n
      }else{
        df_param_hat <- full_join(df_param_hat,df_param_hat_n)
      }
}

  gr <- colorRampPalette(c("#009E73"))(200)                      
  re <- colorRampPalette(c("#E69F00"))(200)
  df_param_vert <- df_param %>%
    gather(alpha_intra,alpha_generic,
         gamma_H_generic,gamma_FV_generic,
         key="parameter",value="estimate") %>%
    filter(!is.na(estimate)) %>%
    mutate(parameter = case_when(parameter=="gamma_H_generic" ~ "Plant-H",
                                 parameter=="gamma_FV_generic" ~ "Plant-Fv",
                                 parameter=="alpha_intra" ~ "Intraspecific",
                                 parameter=="alpha_generic" ~ "Plant-Plant"))
  df_param_all <- full_join(df_param_vert,df_param_hat)
  
  plot.alphas <- ggplot(  df_param_vert) +  
    geom_density_ridges_gradient(aes(x=estimate, y=parameter,
                                     fill= after_stat(x)),
                                 scale = 1) + 
    #geom_boxplot(aes(x=estimate + estimate_hat, y=4, color=parameter_hat ),
    #             width=0.3,alpha=0.4,outlier.shape = NA) +
    facet_wrap(as.factor(focal)~as.factor(year), 
                ncol = 4,nrow=4) +
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
  ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/Figure1.pdf"),
         dpi="retina",
         width = 21,
         height = 16,
         units = c("cm"),
         plot.alphas
  )


#---- 4. Results for one species ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#---- 4.1 Figure of pairwise interactions ----
  df_param_perc  <- NULL
  for( focal in c("CETE","LEMA","HOMA","CHFU")){ # "CHFU","HOMA","CETE"
    for( year in c("2018","2019",'2020','2021')){
      for( complexity.level in 1){
        if(year == "2018" & focal =="CHFU") next
        if(complexity.level == 2 & focal =="HOMA"  & year=="2020") next
        
        complexity.animal <- c("group","family","species")[complexity.level]
        complexity.plant <-c("class","family","code.plant")[complexity.level]
        df_param_perc_n  <- NULL
        df_param_perc_n <-  df_param_all[which(df_param_all$focal==focal  & 
                                            df_param_all$year == year),] %>%
          group_by(parameter) %>%
          mutate(positive.perc = length(which(estimate>0))/length(estimate),
                 neg.perc = length(which(estimate < 0))/length(estimate),
                 null.perc = length(which(estimate==0))/length(estimate))%>%
          ungroup() %>%
          select(focal,year,complexity.plant,parameter,
                 positive.perc, neg.perc,null.perc ) %>%
          unique()
        
        df_param_perc <- bind_rows(df_param_perc,df_param_perc_n)

      }
    }
  }
head(df_param_perc)

ggplot(df_param_perc, aes(y=positive.perc,
                       x=neg.perc, color=parameter)) + 
  geom_point(size=2,aes(shape=as.factor(year))) + facet_grid(focal ~.) +
  ylim(0,1) + xlim(0,1) + 
  scale_color_colorblind() +
  theme_bw()
  

  
  
  
#---- 5.0. Figure for manuscripts ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
#---- 5.1. Fig.2. Distribution direct interactions ----
  gr <- colorRampPalette(c("#009E73"))(200)                      
  re <- colorRampPalette(c("#E69F00"))(200)

  
  plot.gamma.H <-  ggplot(gather(df,gamma_H,
                                 key="gamma_H",value="estimate")) +  
    geom_density_ridges_gradient(aes(x=estimate, y=gamma_H,
                                     fill= after_stat(x)),
                                 scale = 5) + theme_bw() +
    scale_fill_gradientn(
      colours=c(re,"white",gr),limits = c(-2, 2),
      #midpoint = 0,
      #space = "Lab",
      na.value = "grey50",
      guide = "none",
      aesthetics = "fill"
    ) + xlim(-0.6,0.6)
  
  plot.gamma.FV <-  ggplot(gather(df,gamma_FV,
                                  key="gamma_FV",value="estimate")) +  
    geom_density_ridges_gradient(aes(x=estimate, y=gamma_FV,
                                     fill= after_stat(x)),
                                 scale = 5) + theme_bw() +
    scale_fill_gradientn(
      colours=c(re,"white",gr),limits = c(-2, 2),
      #midpoint = 0,
      #space = "Lab",
      na.value = "grey50",
      guide = "none",
      aesthetics = "fill"
    ) + xlim(-0.6,0.6)
  
  
  ggsave( "figures/Chapt1Fig2.pdf", 
          ggarrange(plot.alpha, plot.gamma.H,plot.gamma.FV,ncol=1,nrow=3, 
                    align = c("v"))
  ) 
  
#---- 5.2. Fig.3. Variation intrinsic parameters ----

  plot_lambda_violin_std_all <- ggplot(df_param,
                                   aes(y=scale(lambdas), x=as.factor(year),
                                       fill=as.factor(focal))) +
    geom_violin(position = position_dodge(width = 0.3),
                width=3) + 
    theme_bw() + guides(fill="none") + facet_grid(.~complexity.plant) + 
    labs(y="Intrinsic fecundity (scaled)", x="year") 
  
  
  plot_lambda_violin_std <- ggplot(df_param[which(df_param$complexity.plant=="family"),],
                               aes(y=scale(lambdas), x=as.factor(year),
                                   fill=as.factor(focal))) +
    geom_violin(position = position_dodge(width = 0.3),
                width=3) + 
    theme_bw() + guides(fill="none") + 
    labs(y="Intrinsic fecundity (scaled)", x="year") 
  
  plot_lambda_density <- ggplot(df_param[which(df_param$complexity.plant=="family"),],
                                aes(y=scale(lambdas),fill=as.factor(focal)))+
    geom_density(alpha=0.5) +
    theme_bw() + labs(y="",x="distributions")
  
  
  ggarrange(plot_lambda_violin,plot_lambda_density, widths=(c(2,1)))
  
  
  plot_alpha_intra_violin_all <- ggplot(df_param,
                                       aes(y=alpha_intra, x=as.factor(year),
                                           fill=as.factor(focal)))+
    geom_hline(yintercept = 0, color="grey") +
    geom_violin(position = position_dodge(width = 0.3),
                width=3) + 
    theme_bw() + guides(fill="none") +  ylim(-2,2) +
    labs(y="Intra-specific interactions", x="year") +
    facet_grid(.~complexity.plant)
  
  plot_alpha_intra_violin <- ggplot(df_param[which(df_param$complexity.plant=="family"),],
                                    aes(y=alpha_intra, x=as.factor(year),
                                        fill=as.factor(focal)))+
    geom_hline(yintercept = 0, color="grey") +
    geom_violin(position = position_dodge(width = 0.3),
                width=3) + 
    theme_bw() + guides(fill="none") +  ylim(-2,2) +
    labs(y="Intra-specific interactions", x="year") 
  
  plot_alpha_intra_density <- ggplot(df_param[which(df_param$complexity.plant=="family"),],
                                     aes(y=alpha_intra,fill=as.factor(focal)))+
    geom_density(alpha=0.5) + ylim(-2,2) +
    theme_bw() + labs(y="",x="distributions")
  
  plot_alpha_violin <- ggplot(df_param[which(df_param$complexity.plant=="family"),],
                              aes(y=alpha_generic, x=as.factor(year),
                                  fill=as.factor(focal)))+
    geom_hline(yintercept = 0, color="grey") +
    geom_violin(position = position_dodge(width = 0.3),
                width=3) +
    theme_bw() + guides(fill="none") + ylim(-2,2) +
    labs(y="Generic Inter-specific interactions", x="year") 
  
  plot_alpha_density <- ggplot(df_param[which(df_param$complexity.plant=="family"),],
                               aes(y=alpha,fill=as.factor(focal)))+
    geom_density(alpha=0.5) + ylim(-2,2) +
    theme_bw() + labs(y="",x="distributions")
  
  ggsave( "figures/Chapt1Fig3.pdf", 
          ggarrange(plot_lambda_violin,plot_lambda_density,
                    plot_alpha_intra_violin,plot_alpha_intra_density,
                    plot_alpha_violin,plot_alpha_density,
                    ncol=2,nrow=3,
                    widths=(c(2,1)))
  ) 
  #---- 6.0. SUPP Figure for manuscripts ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #---- 6.1. Fit of model to observations ----
  #---- 6.2. Correlation between alpha_intra and lambda ----
  df.corr <- NULL
  for( focal in c("CETE","LEMA")){ # "CHFU","HOMA",
      res <- cor.test(df_param$alpha_intra[which(df_param$focal==focal)],
           df_param$lambdas[which(df_param$focal==focal)],
           method=c("pearson"))
      if(res$p.value < 0.05){
        significant = "yes"
      }else{significant = "no"}
      
      df.corr <- bind_rows(  df.corr,
                             data.frame( focal = focal, 
                                         p.value =    res$p.value, 
                                         corr =   res$estimate,
                                         significant = significant))
    }


  