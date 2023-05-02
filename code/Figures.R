# This script will make the figures of chap 1, of the natural data part, number- XXX of the manuscript

library(ggplot2)
library(tidyverse)
library(ggsci)
library(tidyverse)
library(ggthemes)
#install.packages("ggridges ")
library(ggridges)
#install.packages("wesanderson")
library(wesanderson)
library(ggpubr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 1.Simulated data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#---- 1.1. Figure of lambda distribution ----
df_lambda <- data.frame()
focal = "LEMA"
for (complexity.plant in  c("class","family")){
  for (years.n in c("2019","2020")){
    if(complexity.plant == "family" & years.n =="2019") next
    load(paste0(project.dic,"results/FinalFit",
                focal,"_",years.n,"_",complexity.plant,"_group.rds"))
    
    df_lambda_n <- FinalFit %>% 
      as.data.frame() %>% 
      dplyr::select("lambdas[1]")
    names(df_lambda_n) <-c("lambda")
    df_lambda_n$year <- years.n
    df_lambda_n$focal <- focal
    df_lambda_n$complexity.plant <- complexity.plant
    df_lambda_n$complexity.animal <- "group"
    df_lambda_n$y <- 
      df_lambda <- bind_rows(df_lambda,df_lambda_n)
  }
}

color_li <-  colorblind_pal()(8)


plot_lambda <-ggplot(df_lambda, aes(x=lambda, y=cbind(year,complexity.plant)))+
  geom_joy(scale = 2, alpha=0.5) +
  scale_y_discrete(name = NULL,expand=c(0.01, 0)) +
  scale_x_continuous(name = NULL,expand=c(0.01, 0)) +
  geom_vline(xintercept=simul_data_list$sim_lambda[1],color=color_li[1]) +
  #scale_colour_manual(values=alpha(color_li[1:8],0.7)) +
  labs(title ="Intrinsic fecundity distributions 
of focal K1 ", x="Intrinsic fecundity",y="density") +  theme_joy()

#---- 1.2. Figure of plant interactions distribution ----

df_alpha<- data.frame()

for (state in c("simple","alternative","complex","verycomplex","sparsity_Final_HTL")){
  df_alpha_n <- get(paste0("TestFit_",state,"model")) %>% 
    as.data.frame() %>% 
    dplyr::select("alpha_generic[1]")
  names(df_alpha_n) <-c("alpha")
  df_alpha_n$model <- state
  df_alpha <- bind_rows(df_alpha,df_alpha_n)
}
df_alpha$model <- factor(df_alpha$model,
                         levels = c("simple","alternative","complex",
                                    "verycomplex","sparsity_Final_HTL"))


library(ggthemes)

plot_alpha <-ggplot(df_alpha,aes(x=alpha, y=model))+
  geom_joy(scale = 2, alpha=0.5) +
  scale_y_discrete(name = NULL,expand=c(0.01, 0)) +
  scale_x_continuous(name = NULL,expand=c(0.01, 0),lim=c(-2,2)) +
  geom_vline(xintercept=simul_data_list$sim_alpha_generic[1],
             color=color_li[1],size=1) +
  #scale_colour_manual(values=alpha(color_li[1:8],0.7)) +
  labs(title ="Plant generic pairwise interaction
distributions of focal K1 ",
       x="Plant generic pairwise interaction",
       y="density") +
  theme_joy()


#---- 1.3. Figure of pollinator interactions distribution ----

df_gamma <- data.frame()

for (state in c("simple","alternative","complex","verycomplex","sparsity_Final_HTL")){
  df_gamma_n <- get(paste0("TestFit_",state,"model")) %>% 
    as.data.frame() %>% 
    dplyr::select("gamma_generic[1]")
  names(df_gamma_n) <-c("gamma")
  df_gamma_n$model <- state
  df_gamma <- bind_rows(df_gamma,df_gamma_n)
}
df_gamma$model <- factor(df_gamma$model,
                         levels = c("simple","alternative","complex",
                                    "verycomplex","sparsity_Final_HTL"))
plot_gamma <-ggplot(df_gamma,aes(x=gamma, y=model))+
  geom_joy(scale = 2, alpha=0.5) +
  scale_y_discrete(name = NULL,expand=c(0.01, 0)) +
  scale_x_continuous(name = NULL,expand=c(0.01, 0),lim=c(-2,2)) +
  geom_vline(xintercept=simul_data_list$sim_HTL_generic[1],
             color=color_li[1],size=1) +
  #scale_colour_manual(values=alpha(color_li[1:8],0.7)) +
  labs(title ="Pollinator generic pairwise interaction
distributions of focal K1 ",
       x="Pollinator generic pairwise interaction",
       y="density")+
  theme_joy()


#---- 2.3. Figure of HOIs distribution ---

df_beta <- data.frame( observation =paste0("Neighbour ", SpNames),
                       beta=simul_data_list_parameters$Spbeta_generic, model="initial")

for (state in c("simple","alternative","complex","verycomplex","sparsity_Final_HTL")){
  df_beta_n <- get(paste0("TestFit_",state,"model")) %>% 
    as.data.frame() %>% 
    dplyr::select(paste0("beta_plant_generic[",c(1:length(SpNames)),"]"))
  names(df_beta_n) <- paste0("Neighbour ", SpNames)
  df_beta_n <- gather( df_beta_n,key="observation",value="beta")
  df_beta_n$model <- state
  df_beta <- bind_rows(df_beta,df_beta_n)
}
df_beta$model <- factor(df_beta$model,
                        levels = c("initial","simple","alternative","complex",
                                   "verycomplex","sparsity_Final_HTL"))
df_beta$observation <- factor(  df_beta$observation,
                                levels = c(paste0("Neighbour ", SpNames)))

df_beta_init <- data.frame( observation =paste0("Neighbour ", SpNames),
                            beta=simul_data_list_parameters$Spbeta_generic, model="initial")
df_beta_init$observation <- factor(  df_beta_init$observation,
                                     levels = c(paste0("Neighbour ", SpNames)))


plot_beta  <- ggplot() + 
  geom_density(data =df_beta,  
               mapping =aes(x=beta,y=..scaled..,color=model),size=1) +
  geom_vline(data =df_beta_init,
             mapping =aes(xintercept = beta),
             color=color_li[1],size=0.5) +
  facet_wrap(~observation, ncol = 3) +
  theme_bw() + theme(strip.background = element_blank(),
                     panel.grid.minor = element_blank())+
  scale_colour_manual(values=alpha(color_li[2:8],0.7)) + 
  guides(color = "none") +
  labs(title ="Plant generic HOIs distribution of focal K1 ",
       x="Plant generic pairwise interaction",
       y="density")

#---- 2.4. Regroup figure ---
library(ggpubr)
ggsave("~/Eco_Bayesian/Test_simulation/figures/SimData_parameters_distributions.pdf",
       dpi="retina",scale=3, 
       ggpubr::ggarrange(plot_lambda,plot_alpha,plot_gamma,
                         ncol = 3, 
                         labels = c("A", "B","C")
                         
       ) 
)       

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 3. Natural data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 3.0 CSV for inclusion----

df_inclusion_nat <- data.frame()
for( focal in c("CETE","CHFU","HOMA","LEMA")){
  for( year in c("2018","2019",'2020','2021')){
    for( complexity.level in 1:2){
      complexity.animal <- c("group","family","species")[complexity.level]
      complexity.plant <-c("class","family","code.plant")[complexity.level]
      
      if(year == "2018" & focal == "CHFU") next
      if(year == "2020" & focal == "HOMA" & complexity.level == 2) next
      load(paste0(project.dic,"results/Inclusion",focal,"_",
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

summary.interactions <- read.csv("results/Chapt1_Inclusion_parameters.csv")
write.csv(df_inclusion_nat,
          file = "~/Eco_Bayesian/Complexity_caracoles/results/Chapt1_Inclusion_parameters.csv")
df_inclusion_nat <- summary.interactions 
n.inclus <- c("n.competitors_plant_inclus","n.competitors_FV_inclus","n.competitors_H_inclus",
              "n.HOIs_plant_inclus","n.HOIs_H_inclus",
              "n.HOIs_FV_inclus","n.HOIs_2FV_inclus","n.HOIs_2H_inclus","n.HOIs_FvH_inclus")
n.inclus  <- n.inclus[which(n.inclus %in% names(df_inclusion_nat))]

inclus <- c("competitors_plant_inclus","competitors_FV_inclus","competitors_H_inclus",
            "HOIs_plant_inclus","HOIs_H_inclus",
            "HOIs_FV_inclus","HOIs_2FV_inclus","HOIs_2H_inclus","HOIs_FvH_inclus")
inclus.short  <- inclus[which(inclus %in% names(df_inclusion_nat))]
n.potential <- c("n.competitors_plant","n.interactors_FV","n.interactors_H",
                 "n.HOIs_plant","n.HOIs_H",
                 "n.HOIs_FV","n.HOIs_2FV","n.HOIs_2h","n.HOIs_FvH")
n.potential  <- n.potential[which(n.potential %in% names(df_inclusion_nat))]

potential <- c("competitors_plant","interactors_FV","interactors_H",
               "HOIs_plant","HOIs_H",
               "HOIs_FV","HOIs_2FV","HOIs_2H","HOIs_FvH")
potential  <- potential[which(potential %in% names(df_inclusion_nat))]

df_inclusion_nat_short <- df_inclusion_nat[rowSums(df_inclusion_nat[,n.inclus])>0,]

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
  select(-all_of(c(inclus.short,n.inclus,n.potential))) %>%
  mutate(#interactor = extract.vert(df_inclusion_nat,inclus,
    #                       "interaction","identity",inclus)[,"value"],
    n.number = extract.vert(df_inclusion_nat,n.inclus,
                            "interaction","n.number",inclus)[,"value"],
    n.potential.interactor = extract.vert(df_inclusion_nat,n.potential ,
                                          "interaction","n.identity.potential",potential)[,"value"])

pal <- wes_palette("Zissou1", 9, type = "continuous")

df_inclusion_nat_ratio <- df_inclusion_nat_vertical %>%
  mutate(ratio =  n.number/n.potential.interactor)  %>%
  mutate(year =  as.integer(year)) %>%
  #mutate(interactor =  as.character(interactor)) %>%
  mutate(interaction =factor(interaction,levels=potential))

df_inclusion_nat_heatmap <- ggplot(df_inclusion_nat_ratio,
                                   aes(x=interaction,y=as.factor(year),fill=ratio)) + 
  scale_fill_gradientn(colours = c("lightgrey", pal),
                       na.value = "white",
                       values =seq(0,1,0.1),
                       lim=c(0,1)) + 
  geom_tile(color="black")+ ylab("") + facet_grid(focal~ complexity_plant) +
  theme_bw() + theme(axis.text.x = element_text(angle=45))

ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/df_inclusion_plot_ratio_ALL.pdf"),
       dpi="retina",scale=1, 
       df_inclusion_nat_heatmap )

#---- 3.1 Figure of pairwise interactions distribution ----

df_alpha_nat <- data.frame()

for( focal in c("CETE","CHFU","HOMA","LEMA")){
  for( year in c("2018","2019",'2020','2021')){
    for( complexity.level in 1:2){
      complexity.animal <- c("group","family","species")[complexity.level]
      complexity.plant <-c("class","family","code.plant")[complexity.level]
      
      if(year == "2018" & focal == "CHFU") next
      load(paste0(project.dic,"results/FinalFit",focal,"_",
                  year,"_",complexity.plant,"_",complexity.animal,".rds"))
      
      assign(paste0("FinalFit",focal,"_",
                    year,"_",complexity.plant,"_",complexity.animal),FinalFit)
      rm(FinalFit)
      
      df_alpha_n <- get(paste0("FinalFit",focal,"_",
                               year,"_",complexity.plant,"_",complexity.animal)) %>% 
        as.data.frame() %>% 
        dplyr::select("alpha_generic[1]")
      names(df_alpha_n) <-c("alpha")
      df_alpha_n$focal <- focal
      df_alpha_n$year <- as.integer(year)
      df_alpha_n$complexity.plant <- complexity.plant
      df_alpha_n$complexity.animal <- complexity.animal
      df_alpha_nat <- bind_rows(df_alpha_nat,df_alpha_n)
      
    }
  }
}
write.csv(df_alpha_nat,
          file = "~/Eco_Bayesian/Complexity_caracoles/results/Chapt1_df_alpha_nat.csv")
df_alpha_nat <- read.csv("~/Eco_Bayesian/Complexity_caracoles/results/Chapt1_df_alpha_nat.csv")

view(df_alpha_nat)

library(ggthemes)
gr <- colorRampPalette(c("dark green"))(200)                      
re <- colorRampPalette(c( "red2"))(200)

df_alpha_nat <- df_alpha_nat[which(df_alpha_nat$alpha > -1 & df_alpha_nat$alpha < 1),]


plot_alpha_nat <- ggplot(df_alpha_nat,aes(x=alpha, y=as.factor(year),
                                          fill = stat(x),group = year)) +
  geom_density_ridges_gradient(scale = 1.4,alpha=0.5, 
                               rel_min_height = 0.01,
                               quantile_lines = TRUE, quantiles = 2) +
  scale_fill_gradientn(
    colours=c(re,"white",gr),limits = c(-1, 1),
    #midpoint = 0,
    #space = "Lab",
    na.value = "grey50",
    guide = "none",
    aesthetics = "fill"
  ) +
  geom_vline(xintercept = 0, color="darkgrey") +
  scale_y_discrete(expand=c(0.01, 0), #limits=c('2018','2021'),
                   breaks = c('2018','2019','2020','2021')) +
  scale_x_continuous(name = NULL,expand=c(0.01, 0)) +
  facet_wrap( focal~complexity.plant, ncol = 2, scales= "free_x")+
  #scale_colour_manual(values=alpha(color_li[1:8],0.7)) +
  labs(title =paste0("Plant generic pairwise interaction distributions"),
       x="Plant generic pairwise interaction",
       y="density") +
  theme_bw() +
  theme(axis.text=element_text(size=18)) 

ggsave(paste0("~/Eco_Bayesian/Complexity_caracoles/figure/",focal,"_GenericAlpha.pdf"),
       dpi="retina",scale=1, 
       plot_alpha_nat
)
#---- 3.2. Figure of pairwise interactions POINTs ----

df_alphas <- data.frame()

for( focal in c("CETE","CHFU","HOMA","LEMA")){
  for( year in c("2018","2019",'2020','2021')){
    for( complexity.level in 1:2){
      complexity.animal <- c("group","family","species")[complexity.level]
      complexity.plant <-c("class","family","code.plant")[complexity.level]
      
      if(year == "2018" & focal == "CHFU") next
      if(year == "2020" & focal == "CHFU") next
      
      load(paste0(project.dic,"/results/FinalFit",focal,"_",
                  year,"_",complexity.plant,"_",complexity.animal,".rds"))
      
      df_alpha_n <- as.data.frame(FinalFit) %>% 
        dplyr::select("alpha_generic[1]") %>% 
        apply(2, summary) %>% 
        as.data.frame() %>%
        rownames_to_column(var="stats")%>%
        spread(key="stats",value="alpha_generic[1]")
      
      
      
      assign(paste0("FinalFit",focal,"_",
                    year,"_",complexity.plant,"_",complexity.animal),FinalFit)
      rm(FinalFit)
      
      df_alpha_n$focal <- focal
      df_alpha_n$year <- as.numeric(year)
      df_alpha_n$complexity.plant <- complexity.plant
      df_alpha_n$complexity.animal <- complexity.animal
      df_alphas <- bind_rows(df_alphas,df_alpha_n)
    }
  }
}

df_alphas[which(df_alphas$focal == "LEMA"),c(3,4)] <- df_alphas[which(df_alphas$focal == "LEMA"),c("X1st.Qu..1","X3rd.Qu..1")] 
df_alphas <- df_alphas[,-c(1,2)]
df_alphas <- df_alphas[,-c(11,12)]
names(df_alphas) <- c("FirstQ","ThirdQ", "Max","Mean","Median","Min","focal","year",           
                      "complexity.plant","complexity.animal")

write.csv(  df_alphas,
            file = "~/Eco_Bayesian/Complexity_caracoles/results/Chapt1_df_alphas.csv")
df_alphas <- read.csv("~/Eco_Bayesian/Complexity_caracoles/results/Chapt1_df_alphas.csv")

assign(paste0("df_param_",focal),df)

ggplot(df_alphas,aes(x=Mean, y = year, group=complexity.plant)) + 
  geom_point(aes(shape=complexity.plant, 
                 color=cut(Mean, c(-Inf,0, 3))),
             size = 5) +
  geom_errorbarh(aes(xmin=FirstQ,
                     xmax=ThirdQ,
                     color=cut(Mean, c(-Inf,0, 3))),
                 height=0.1,linewidth=0.5) +
  geom_vline(xintercept=0,color="black") +
  scale_color_manual(name = "",
                     values = c("(-Inf,0]" = "red",
                                "(0,3]" = "dark green"),
                     labels = c("competition", "facilitation")) +
  facet_grid(.~focal,scales="free") +
  labs(title ="Generic pairwise interaction distributions", 
       x="Plants generic pairwise interaction") +
  theme_bw() 






#---- 3.3. Figure of intrinsic growth rate  ----
df_lambda_nat <- data.frame()

for( year in c("2019",'2020','2021')){
  for( complexity.plant in c("class","family")){
    
    df_lambda_n <- get(paste0("FinalFit",focal,"_",
                              year,"_",complexity.plant,"_",complexity.animal)) %>% 
      as.data.frame() %>% 
      dplyr::select("lambdas[1]")
    names(df_lambda_n) <-c("lambda")
    df_lambda_n$focal <- focal
    df_lambda_n$year <- year
    df_lambda_n$complexity.plant <- complexity.plant
    df_lambda_n$complexity.animal <- complexity.animal
    
    df_lambda_nat <- bind_rows(df_lambda_nat,df_lambda_n)
  }
}

color_li <-  colorblind_pal()(8)


plot_lambda_nat <- ggplot(df_lambda_nat, aes(x=lambda, y=year))+
  geom_joy(scale = 2, alpha=0.5) +
  scale_y_discrete(name = NULL,expand=c(0.01, 0)) +
  scale_x_continuous(name = NULL,expand=c(0.01, 0)) +
  facet_grid(.~complexity.plant) +
  labs(title =paste0("Intrinsic fecundity distributions of ",focal), 
       x="Intrinsic fecundity",y="density") +  theme_joy() +
  theme(axis.text=element_text(size=10))

#---- 3.4. Figure of pollinator interactions distribution ----

df_gamma_FV_nat <- data.frame()
for( focal in c("CETE","CHFU","HOMA","LEMA")){
  for( year in c("2018","2019",'2020','2021')){
    for( complexity.level in 1:2){
      complexity.animal <- c("group","family","species")[complexity.level]
      complexity.plant <-c("class","family","code.plant")[complexity.level]
      
      df_gamma_FV_n <- get(paste0("FinalFit",focal,"_",
                                  year,"_",complexity.plant,"_",complexity.animal)) %>% 
        as.data.frame() %>% 
        dplyr::select("gamma_FV[1]")
      names(df_gamma_FV_n) <-c("gamma_FV")
      
      df_gamma_FV_n$focal <- focal
      df_gamma_FV_n$year <- year
      df_gamma_FV_n$complexity.plant <- complexity.plant
      df_gamma_FV_n$complexity.animal <- complexity.animal
      
      df_gamma_FV_nat <- bind_rows(df_gamma_FV_nat,df_gamma_FV_n)
    }
  }
}


plot_gamma_FV_nat <-ggplot(df_gamma_FV_nat,aes(x=gamma_FV, y=year))+
  geom_joy(scale = 2, alpha=0.5) +
  scale_y_discrete(name = NULL,expand=c(0.01, 0)) +
  scale_x_continuous(name = NULL,expand=c(0.01, 0)) +
  facet_grid(.~complexity.plant) +
  labs(title =paste0("Pollinator generic pairwise interaction
distributions with ",focal), 
       x="Pollinator generic pairwise interaction",
       y="density") +  theme_joy()+
  theme(axis.text=element_text(size=10))

#---- 3.5. Figure of Herbivore interactions distribution ----

df_gamma_H_nat <- data.frame()
for( year in c("2019",'2020','2021')){
  for( complexity.plant in c("class","family")){
    
    df_gamma_H_n <- get(paste0("FinalFit",focal,"_",
                               year,"_",complexity.plant,"_",complexity.animal)) %>% 
      as.data.frame() %>% 
      dplyr::select("gamma_H[1]")
    names(df_gamma_H_n) <-c("gamma_H")
    
    df_gamma_H_n$focal <- focal
    df_gamma_H_n$year <- year
    df_gamma_H_n$complexity.plant <- complexity.plant
    df_gamma_H_n$complexity.animal <- complexity.animal
    
    df_gamma_H_nat <- bind_rows(df_gamma_H_nat,df_gamma_H_n)
  }
}


plot_gamma_H_nat <-ggplot(df_gamma_H_nat,aes(x=gamma_H, y=year))+
  geom_joy(scale = 2, alpha=0.5) +
  scale_y_discrete(name = NULL,expand=c(0.01, 0)) +
  scale_x_continuous(name = NULL,expand=c(0.01, 0)) +
  facet_grid(.~complexity.plant) +
  labs(title =paste0("Herbivore generic pairwise interaction
distributions with ",focal), 
       x="Pollinator generic pairwise interaction",
       y="density") +  theme_joy()+
  theme(axis.text=element_text(size=10))

#---- 3.6. Regroup figure ----
library(ggpubr)
ggsave("~/Eco_Bayesian/Complexity_caracoles/figure/Coeff_parameters.pdf",
       dpi="retina",scale=4, 
       ggpubr::ggarrange(plot_lambda_nat,plot_alpha_nat,
                         plot_gamma_H_nat, plot_gamma_FV_nat,
                         
                         ncol = 2, nrow=2,
                         labels = c("A", "B","C",'D')
                         
       ) 
)          

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- 4. Results for one species ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#---- 4.1 Figure of pairwise interactions ----

focal <- "LEMA"

df <- data.frame()

for( year in c("2018","2019",'2020','2021')){
  for( complexity.level in 1:2){
    complexity.animal <- c("group","family","species")[complexity.level]
    complexity.plant <-c("class","family","code.plant")[complexity.level]
    
    load(paste0(project.dic,"results/FinalFit",focal,"_",
                year,"_",complexity.plant,"_",complexity.animal,".rds"))
    
    assign(paste0("FinalFit",focal,"_",
                  year,"_",complexity.plant,"_",complexity.animal),FinalFit)
    rm(FinalFit)
    
    df_alpha_n <- get(paste0("FinalFit",focal,"_",
                             year,"_",complexity.plant,"_",complexity.animal)) %>% 
      as.data.frame() %>% 
      dplyr::select("lambdas[1]","alpha_generic[1]","gamma_H[1]","gamma_FV[1]")
    
    names(df_alpha_n) <-c("lambda","alpha","gamma_H","gamma_FV")
    df_alpha_n$focal <- focal
    df_alpha_n$year <- year
    df_alpha_n$complexity.plant <- complexity.plant
    df_alpha_n$complexity.animal <- complexity.animal
    df<- bind_rows(df,df_alpha_n)
    
  }
}

write.csv(df,
          file = paste0("~/Eco_Bayesian/Complexity_caracoles/results/Chapt1_df_alpha",focal,".csv"))
df <- read.csv(paste0("~/Eco_Bayesian/Complexity_caracoles/results/Chapt1_df_alpha",focal,".csv"))

assign(paste0("df_param_",focal),df)




gr <- colorRampPalette(c("dark green"))(200)                      
re <- colorRampPalette(c( "red2"))(200)

df_param_LEMA <- df_param_LEMA[,-1]
df_param_LEMA_vert <- gather(df_param_LEMA,lambda,alpha,gamma_H,gamma_FV, 
                             key="parameter", value="estimation")

df_param_LEMA_vert <- df_param_LEMA_vert[which((!df_param_LEMA_vert$parameter == "lambda" & df_param_LEMA_vert$estimation > -1 & df_param_LEMA_vert$estimation < 1) |  df_param_LEMA_vert$parameter == "lambda"),]

plot_parameters_nat <-   ggplot(df_param_LEMA_vert,aes(x=estimation, y=year,height=after_stat(density),
                                                       fill = after_stat(x),group = year))+
  geom_density_ridges_gradient(scale = 1.5,alpha=0.9,
                               rel_min_height = 0.001,
                               quantile_lines = TRUE, quantiles = 2) +
  scale_fill_gradientn(
    colours=c(re,"white",gr),limits = c(-2, 2),
    #midpoint = 0,
    #space = "Lab",
    na.value = "grey50",
    guide = "none",
    aesthetics = "fill"
  ) +
  scale_y_discrete(name = NULL,expand=c(0.01, 0)) +
  geom_vline(xintercept=0,color="darkgrey") +
  #scale_x_continuous(name = NULL,expand=c(0.01, 0),lim=c(-1,1)) +
  facet_wrap(parameter~ complexity.plant, scales = "free_x", ncol=2) +
  #scale_colour_manual(values=alpha(color_li[1:8],0.7)) +
  labs(title =paste0("Plant generic pairwise interaction distributions of ",focal),
       x="Plant generic pairwise interaction",
       y="density") +
  theme_bw()+
  theme(axis.text=element_text(size=10))

ggsave("~/Eco_Bayesian/Complexity_caracoles/figure/Coeff_parameters_LEMA.pdf",
       dpi="retina",scale=1, plot_parameters_nat 
)    


