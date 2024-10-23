options(mc.cores = parallel::detectCores())
library(rstan)
library(cli)
library("codetools")
library(HDInterval)
library(matrixStats)
library(tidyverse)
library(stringr)
library(dplyr)
library(ggplot2) 
library(shinystan) 
library(Hmisc) 
#### predictive check for all model 
home.dic <- "/home/lbuche/Eco_Bayesian/Complexity_caracoles/"
project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/"
ModelCheck <- NULL
postfecundity <- list()
model.compare <- NULL
pdf(paste0(home.dic,"figure/ModelsBehavior.pdf"))
#title= paste0("Appendix  : Diagnostic plots for Final fit of all models")

for( focal in c("LEMA","CHFU","HOMA","CETE")){ # "CHFU","HOMA","CETE"
  for( year in c("2019",'2020','2021')){
    for( complexity.level in c(1:3)){
      
      complexity.plant <-c("class","family","code.plant")[complexity.level]
      complexity.animal <- c("group","family","species")[complexity.level]
      FocalPrefix <- focal
      print(paste(focal,complexity.plant,complexity.animal,year))
      
      #---- 1.2. Import the raw data ----
      
      #setwd("~/Eco_Bayesian/Complexity_caracoles")
      home.dic <- "/home/lbuche/Eco_Bayesian/Complexity_caracoles/"
      project.dic <- "/data/projects/punim1670/Eco_Bayesian/Complexity_caracoles/"
      
      #rm(list = ls())
      abundance.plant <- read.csv(paste0(home.dic,"data/abundance.csv"), sep=",")
      
      competition.plant <- read.csv(paste0(home.dic,"data/competition.csv"), sep=",")
      
      plant.class <- read.csv(paste0(home.dic,"data/plant_code.csv"),sep=",")
      
      herbivorypredator <- read.csv(paste0(home.dic,"data/herbivorypredator.csv"), sep=",")
      
      floral_visitor <- read.csv(paste0(home.dic,"data/floral_visitor.csv"), sep=",")
      
      groupinglevels <- read.csv(paste0(home.dic,"data/groupinglevels.csv"), sep=",")
      
      # Run UploadData.R 2. and 3. or uncomment the following lines 
      plant.class <- read.csv(paste0(home.dic,"data/plant_code.csv"), sep=",")
      
      SpData <- read.csv(paste0(home.dic,"data/SpData.csv"))
      
      #---- 2.EXTRACT MATRIXES OF ABUNDANCES -----
      
      source(paste0(home.dic,"code/ExtractMatrix_fct.R"))
      extract.matrix(focal,year=year,
                     complexity.plant ,
                     complexity.animal,
                     plant.class,
                     SpData)
      #Returns two elemens
      #view(summary.interactions.n)
      #str(DataVec) #contains 21 elements
      #----- 2. Output of prelminiray fit ----
      
      load(file=paste0(home.dic,"results/Inclusion",FocalPrefix,"_",year,
                       "_",complexity.plant,"_",complexity.animal,".RData"))
      
      
      
      #---- 3. FINAL FIT----
      #---- 3.1. Set up parameters ---- 
      
      load(paste0(project.dic,"results/stan/FinalFit",FocalPrefix,"_",year,"_",complexity.plant,"_",complexity.animal,".rds")) 
      
      FinalPosteriors <- rstan::extract(FinalFit)
      
      
      #---- 4.3. Final fit posterior check and behavior checks---- 
      
      ##### Diagnostic plots and post prediction 
      # Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
      source("~/Eco_Bayesian/Test_simulation/code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots
      # check the distribution of Rhats and effective sample sizes 
      ##### Posterior check
      
      stan_post_pred_check_plot(FinalPosteriors,"F_hat",
                                DataVec$Fecundity,max(DataVec$Fecundity)+20,
                                title= paste0(focal," in year ",year," for complexity: ",complexity.plant)) 
      
      # N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
      
      
      hist(summary(FinalFit)$summary[,"Rhat"],
           main = paste("Finat Fit: Histogram of Rhat for",
                        FocalPrefix, " in ", year," \nwhen complexity is ",complexity.plant))
      hist(summary(FinalFit)$summary[,"n_eff"],
           main = paste("Finat Fit: Histogram of Neff for",
                        FocalPrefix, " in ", year," \nwhen complexity is ",complexity.plant))
      #
      # plot the corresponding graphs
      par <- c('lambdas','alpha_generic','alpha_intra',"gamma_FV_generic","gamma_H_generic")
      if(DataVec$RemoveFV == 1){par <- c('lambdas','alpha_generic',
                                         'alpha_intra',"gamma_H_generic")
      }
      if(DataVec$RemoveH == 1){par <- c('lambdas',
                                        'alpha_generic','alpha_intra',"gamma_FV_generic")
      }
      if(DataVec$RemoveH == 1 & DataVec$RemoveFV == 1){par <- c('lambdas',
                                                                'alpha_generic','alpha_intra')
      }
      
      trace <- stan_trace(FinalFit , pars=par,inc_warmup = TRUE)
      print(trace)
      dens <- stan_dens(FinalFit , pars=par)
      print(dens)
      plot.stan <- stan_plot(FinalFit , pars=par)
      print( plot.stan)
      pairsplot <-  pairs(FinalFit,pars=par)
      print(pairsplot)
      
      # Next check the Rsquarred of fitted data
      
      bayes_R2 <- function(post.draws,
                           var_name){
        mu <- post.draws[[var_name]]
        var_mu <- apply(mu, 1, var)
        
        disp_dev <- post.draws$disp_dev[1]
        phi <- (disp_dev^2)^(-1)
        var_phi <- mu + (mu^2/phi)
        
        var_mu/ (var_mu + var_phi)
      }
      
      bayes_RMDS <- function(post.draws,
                             var_name,
                             obs.var){
        y <- post.draws[[var_name]]
        x <-  obs.var
        n <- length(x)
        rmds <- c()
        for(i in 1:nrow(y)){
          rmds[i] <- sqrt(sum((x - y[i,])^2)/n)
        }
        data.frame(mean.fec=mean(x),
                   sd.fec= sd(x),
                   mean.rmds = mean(rmds),
                   min.rmds = min(rmds),
                   max.rmds =max(rmds),
                   Sim.Bias = (mean(x)-mean(y))^2)
      }
      
      
      Rsquarred <- bayes_R2(FinalPosteriors,
                            var_name = 'F_hat') 
      
      # extract log likelihood 
      loo.fit <- rstan::loo(FinalFit,pars ="F_sim")
      
      loo.fit.df <- as.data.frame(loo.fit$estimates)
      
      save(file=paste0(project.dic,"results/stan/Loo.fit ",
                       FocalPrefix,"_",year,"_",complexity.level,".rds"),
           loo.fit)
      
      # number of parameter: lambda, generic plant plant, generic pollinator, generic herbivore and
      # intrasepcific plant plant = 5
      n.parameter <- length(par) + sum(Inclusion_all$Inclusion_ij,Inclusion_all$Inclusion_FV,
                                       Inclusion_all$Inclusion_H,Inclusion_all$beta_Inclusion_plant,
                                       Inclusion_all$beta_Inclusion_FV,Inclusion_all$beta_Inclusion_H)
      
      # make data frame with rhat 
      mc <- data.frame(focal =FocalPrefix, 
                       year =year,
                       complexity = complexity.plant,
                       n.parameter = n.parameter,
                       n.obs = DataVec$N,
                       Rhat = max(summary(FinalFit)$summary[,"Rhat"],na.rm =T),
                       Neff = min(summary(FinalFit)$summary[,"n_eff"],na.rm = T),
                       Rsquarred=Rsquarred %>% mean,
                       Rsquarred.sd = Rsquarred %>% sd()) %>%
        bind_cols( bayes_RMDS(FinalPosteriors,
                              var_name = 'F_hat',
                              DataVec$Fecundity)) %>%
        mutate(perc.K = sum(loo.fit$diagnostics$pareto_k > min(1-(1/log10(DataVec$N)),0.7))/
                 length(loo.fit$diagnostics$pareto_k),
               p_loo = paste0(round(loo.fit.df["p_loo","Estimate"],digits = 1),"+-",
                              round(loo.fit.df["p_loo","SE"],digits = 1))) %>%
        mutate(p_loo.check = case_when(loo.fit.df["p_loo","Estimate"]-loo.fit.df["p_loo","SE"] < n.parameter ~ "Good specification",
                                       T~"misspecification")) %>%
        mutate(elpd_loo = loo.fit.df["elpd_loo","Estimate"],
               elpd_loo_se = loo.fit.df["elpd_loo","SE"])
      
      ModelCheck <- bind_rows(ModelCheck,mc)
      remove(FinalFit)
      
    }
    load(file=paste0(project.dic,"results/stan/Loo.fit ",
                     FocalPrefix,"_",year,"_",1,".rds"))
    assign("loofit.Class",
           loo.fit) 
    
    load(file=paste0(project.dic,"results/stan/Loo.fit ",
                     FocalPrefix,"_",year,"_",2,".rds"))
    assign("loofit.Family",
           loo.fit) 
    
    load(file=paste0(project.dic,"results/stan/Loo.fit ",
                     FocalPrefix,"_",year,"_",3,".rds"))
    assign("loofit.Species",
           loo.fit) 
    model.compare.n <- loo::loo_compare(loofit.Class,loofit.Family,loofit.Species) %>%
      as.data.frame() %>%
      rownames_to_column(var="model") %>%
      mutate(model= case_when(model=="model1" ~ "Class",
                              model=="model2" ~ "Family",
                              model=="model3" ~ "Species")) %>%
      mutate(Sign.int = se_diff*1.96) %>%
      mutate(significance = case_when(abs(elpd_diff) > abs(Sign.int) ~ "significantly diff",
                                      T~"not significantly diff")) %>%
      mutate(focal=focal,
             year=year)
    model.compare <- bind_rows(  model.compare, model.compare.n)
    
    write.csv(model.compare,
              file=paste0(home.dic,"results/Model.compare.loo.csv"))
    write.csv(ModelCheck, 
              paste0(home.dic ,"results/ModelCheck.csv"))
  }
}
dev.off()



Model.check <- read.csv(paste0(home.dic ,"results/ModelCheck.csv"))

#---- 2.Illustration -----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ggsave(
  plot=ggarrange(ggplot(Model.check, aes(sd.fec,mean.rmds,color=focal)) + 
                   geom_jitter(size=3,width = 3, height = 3,
                               aes(shape=as.factor(year)))+
                   geom_abline()+
                   labs(y="Mean RMDS across chains",
                        x="Standard deviation of observed fecundity",
                        shape="year") +
                   theme_bw(), 
                 ggplot(Model.check, aes(elpd_loo,mean.rmds,color=focal)) +
                   scale_x_reverse() +
                   geom_jitter(size=3,width = 3, height = 3,
                               aes(shape=as.factor(year)))+
                   labs(y="Mean RMDS across chains",
                        x="ELPD_loo",
                        shape="year") +
                   theme_bw(),
                 nrow=1,
                 labels = c("a.","b."),
                 common.legend = T,
                 legend="bottom"
  ),
  filename=paste0(home.dic,"figure/ModelBehaviorRMDS.linear.pdf"))


df_param_all <- read.csv(paste0(home.dic,"results/Chapt1_Parameters_values.csv"))

Trophic.mean.df <- read.csv(paste0(home.dic,"results/Trophic.mean.df.csv"))

abund.mean.df <- read.csv(paste0(home.dic,"results/abund.mean.df.csv"))

param.fec.pred <- NULL
for( focal in c("LEMA","CHFU","HOMA","CETE")){ # "CHFU","HOMA","CHFU","HOMA","CETE"
  for( year in c("2019",'2020','2021')){ # '2020','2021'"2018" "2019",'2020','2021'
    for( complexity.level in c(1:3)){ #complexity.level.int
      
      
      summary.interactions.n <- data.frame()
      complexity.plant <-c("class","family","code.plant")[complexity.level]
      complexity.animal <- c("group","family","species")[complexity.level]
      FocalPrefix <- focal
      print(paste(focal,complexity.plant,complexity.animal,year))
      
      load(file=paste0(home.dic,"results/Inclusion",focal,"_",year,
                       "_",complexity.plant,"_",complexity.animal,".RData"))
      
      
      df_param_n <- df_param_all[which(df_param_all[,"focal"] %in% focal &
                                         df_param_all[,"year"] %in%  year &
                                         df_param_all[,"complexity.plant"] %in%  complexity.plant),]
      
      Trophic.mean.df.n <- Trophic.mean.df[which(Trophic.mean.df[,"focal"] %in% focal &
                                                   Trophic.mean.df[,"year"] %in%  year &
                                                   Trophic.mean.df[,"complexity.plant"] %in%  complexity.plant),] %>%
        mutate(year = as.numeric(year))
      
      abund.mean.df.n <- abund.mean.df[which(abund.mean.df[,"focal"] %in% focal &
                                               abund.mean.df[,"year"] %in%  year &
                                               abund.mean.df[,"complexity.plant"] %in%  complexity.plant),] %>%
        rename("abund.sp.min_hat"="abund.sp.min") %>%
        rename("abund.sp.max_hat"="abund.sp.max") %>%
        rename("abund.sp.mean_hat"="abund.sp") %>%
        dplyr::select("focal","year","complexity.plant","parameter_hat","abund.sp.mean_hat",
                      "abund.sp.min_hat", "abund.sp.max_hat","parameter") %>%
        mutate(year = as.numeric(year))
      
      param.effect <- df_param_n %>%
        left_join(Trophic.mean.df.n,
                  by = c("focal", "year", "complexity.plant",
                         "complexity.animal", "parameter")) %>%
        left_join(abund.mean.df.n,
                  by = c("focal", "year", "complexity.plant", "parameter",
                         "parameter_hat")) %>%
        dplyr::select("focal", "year", "complexity.plant",
                      "lambdas","parameter","estimate",
                      "abund.sp.min","abund.sp.mean","abund.sp.max")
      if(nrow(param.effect)>16000){
        param.effect <- param.effect  %>%
          unique() }
      
      param.effect <- param.effect  %>%
        mutate(n.ite= rep(1:4000,times=4)) %>%
        mutate(effect.min = estimate*abund.sp.min,
               effect.mean = estimate*abund.sp.mean,
               effect.max = estimate*abund.sp.max) %>%
        dplyr::select("n.ite","focal", "year", "complexity.plant",
                      "lambdas","parameter","effect.min","effect.mean","effect.max")
      
      param.effect_hat <- df_param_n %>%
        left_join(Trophic.mean.df.n,
                  by = c("focal", "year", "complexity.plant",
                         "complexity.animal", "parameter")) %>%
        left_join(abund.mean.df.n,
                  by = c("focal", "year", "complexity.plant", "parameter",
                         "parameter_hat")) %>%
        dplyr::select("focal", "year", "complexity.plant",
                      "lambdas","parameter","parameter_hat",
                      "estimate_hat","abund.sp.min_hat",
                      "abund.sp.mean_hat","abund.sp.max_hat") %>%
        dplyr::filter(!is.na(estimate_hat)) 
      if(nrow( param.effect_hat)>0){
        param.effect_hat <- param.effect_hat %>%
          mutate(n.ite= rep(1:4000,times=nrow(.)/4000)) %>%
          mutate(effect.min_hat = estimate_hat*abund.sp.min_hat,
                 effect.mean_hat = estimate_hat*abund.sp.mean_hat,
                 effect.max_hat = estimate_hat*abund.sp.max_hat) %>%
          dplyr::select("n.ite","focal", "year", "complexity.plant",
                        "lambdas","parameter",
                        "effect.min_hat","effect.mean_hat","effect.max_hat") %>%
          aggregate(.~parameter + lambdas + n.ite + focal + year + complexity.plant, sum)
        
        param.effect <- param.effect %>%
          left_join(param.effect_hat)
      }else{
        param.effect <- param.effect %>%
          mutate(effect.mean_hat=0,
                 effect.max_hat=0,
                 effect.min_hat=0)
      }
      for(i in 1:4000){
        param.effect.n <- param.effect[param.effect$n.ite==i,]
        param.fec.pred.n <-  data.frame(param.effect.n[1,c("n.ite","focal","year","complexity.plant")],
                                        effect.total.mean = Inclusion_all$U*param.effect.n$lambdas[1] + sum(param.effect.n$effect.mean + sum(param.effect.n$effect.mean_hat,na.rm = T)),
                                        effect.total.min = Inclusion_all$U*param.effect.n$lambdas[1] + sum(param.effect.n$effect.min + sum(param.effect.n$effect.min_hat,na.rm = T)),
                                        effect.total.max = Inclusion_all$U*param.effect.n$lambdas[1] + sum(param.effect.n$effect.max + sum(param.effect.n$effect.max_hat,na.rm = T)))
        if(param.fec.pred.n$effect.total.min > param.fec.pred.n$effect.total.max ){
          names(param.fec.pred.n)[6:7] <- c("effect.total.max",
                                            "effect.total.min")
        }
        param.fec.pred <- bind_rows(param.fec.pred,param.fec.pred.n)
      }
    } 
  }
}

ModelCheck <- read.csv(paste0(home.dic,"results/ModelCheck.csv")) %>%
  rename("complexity.plant"=complexity)
y.lim <- data.frame(focal=rep(c("CETE","CHFU","HOMA","LEMA"),each=3),
                    year=rep(c(2019,2020,2021),times=4),
                    upp = c(1500,1500,1500,
                            1000,300,750,
                            100,700,700,
                            1000,300,500)
)
library(ggridges)
range.plot.rmds <- param.fec.pred %>%
  inner_join(y.lim) %>%
  filter(exp(effect.total.max) < upp &
           exp(effect.total.min) < upp &
           exp(effect.total.mean) < upp) %>%
  ggplot() +
  geom_density_ridges(aes(y=complexity.plant,x=exp(effect.total.min)),
                      alpha=0.2,fill="blue",scale=1) +
  geom_density_ridges(aes(y=complexity.plant,x=exp(effect.total.max)),
                      alpha=0.2,fill="red",scale=1) +
  geom_density_ridges(aes(y=complexity.plant,x=exp(effect.total.mean)),
                      alpha=0.2,fill="orange",scale=1) +
  #stat_summary(aes(y=complexity.plant,x=exp(effect.total.mean)),
  #            color="orange",
  #           fun = "median", geom = "point", size = 2) +
  #stat_summary(aes(y=complexity.plant,x=exp(effect.total.max)),
  #            fun = "median", geom = "point", 
  #           size = 2,color="red",alpha=0.5) +
  #stat_summary(aes(y=complexity.plant,x=exp(effect.total.min)),
  #         fun = "median", geom = "point",
  #         size = 2,color="blue",alpha=0.5) +
  
  geom_errorbarh(data=ModelCheck,height = .2,
                 aes(y=complexity.plant,
                     xmin=min.rmds,
                     xmax=max.rmds),
                 color="black",size=1,alpha=0.6,
                 position=position_dodge(0.3)) +
  geom_point(data=ModelCheck,size=2,shape=19,color="green",
             aes(y=complexity.plant,
                 x=mean.rmds),
             color="black",
             position=position_dodge(0.3)) +
  theme_bw() +
  labs(y="grouping factor",
       x="Predicted fecundity")+
  #scale_x_continuous(aes(limits=c(0,upp))) +
  facet_wrap(focal ~ year,ncol=3,
             scale="free") +
  guides(shape=guide_legend("RMDS",
                            override.aes = list(size = 10,color="black"),
                            nrow = 1,
                            direction="horizontal",
                            byrow = T,
                            title.hjust = 0.1))+
  theme(strip.placement = "outside",
        legend.key.size = unit(1, 'cm'),
        title =element_text(size=12),
        axis.text.x= element_text(size=16),
        axis.text.y= element_text(size=16),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        strip.text = element_text(size=12),
        panel.border = element_rect(color = "black", 
                                    fill = NA))
range.plot.rmds

ggsave(paste0(home.dic,"figure/range.plot.rmds.pdf"),
       #dpi="retina",
       width = 30,
       height = 20,
       units = c("cm"),
       range.plot.rmds
)
