# 00 - Import package and Data
#################################################################################
#---- install C++ ----
# Source = https://github.com/rmacoslib/r-macos-rtools#how-do-i-use-the-installer
# and https://github.com/stan-dev/rstan/wiki/Configuring-C---Toolchain-for-Mac

# enable some compiler optimizations to improve the estimation speed of the model:
dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, "Makevars")
if (!file.exists(M)) file.create(M)
cat("\nCXX14FLAGS += -O3 -mtune=native -arch x86_64 -ftemplate-depth-256",
    file = M, sep = "\n", append = FALSE)

#---- install Stan ---- 
# Source = https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
# remove any existing RStan via
remove.packages("rstan")
if (file.exists(".RData")) file.remove(".RData")

# Restart R, then run the following
#Sys.setenv(DOWNLOAD_STATIC_LIBV8 = 1) # only necessary for Linux without the nodejs library / headers
install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)

# Verify Installation

example(stan_model, package = "rstan", run.dontrun = TRUE) # run the RStan example/test model:

library("rstan") # observe startup messages
rstan_options(auto_write = TRUE)

#---- install BRMS package ---- 

if(!require(brms)){install.packages("brms"); library(brms)} # for the analysis
if(!require(haven)){install.packages("haven"); library(haven)}# to load the SPSS .sav file
library(tidyverse) # needed for data manipulation.
library(RColorBrewer) # needed for some extra colours in one of the graphs
library(ggplot2)
# to dowload ggmcmc I had to first upload reshape (reshape_0.8.8.tar.gz) and 
# GGally (GGally_2.1.2.tar.gz)
install.packages("/Users/lisabuche/Downloads/GGally_2.1.2.tar.gz",
                 repos=NULL, type="source", dependencies=TRUE) # need ggplot2 >=3.3.4
#require("devtools")
library(reshape)
library(GGally) #you may have to restar R for this to work
install.packages("/Users/lisabuche/Downloads/ggmcmc_1.5.1.1.tar.gz",
                 repos=NULL, type="source", dependencies=TRUE)

library(ggmcmc) # to load the SPSS .sav file and graphical tools with ggplot2
library(ggthemes)
if(!require(ggridges)){install.packages("ggridges"); library(ggridges)}# to load the SPSS .sav file

#################################################################################
# Building a Multilevel Model in BRMS Tutorial: Popularity Data
# 01 - Import Data and first look
#################################################################################


#~~~~ 01 - Import Data and first look----
#---- Downloading the data ----
# Description of the data: the main outcome variable is the pupil popularity,
# a popularity rating on a scale of 1–10 derived by a sociometric procedure.
# Typically, a sociometric procedure asks all pupils in a class to rate all the
# other pupils, and then assigns the average received popularity rating to each
# pupil. Because of the sociometric procedure, group effects as apparent from 
# higher-level variance components are rather strong. There is a second outcome 
# variable: pupil popularity as rated by their teacher, on a scale from 1 to 10.
# The explanatory variables are pupil gender (boy = 0, girl = 1), pupil 
# extraversion (10-point scale), and teacher experience in years.
# Goal of the data: find models and test hypotheses about the relation 
# between the characteristics of pupils in different classes and the popularity of pupils (according to their classmates).

popular2data <- read_sav(file = "https://github.com/MultiLevelAnalysis/Datasets-third-edition-Multilevel-book/blob/master/chapter%202/popularity/SPSS/popular2.sav?raw=true")

popular2data <- select(popular2data, pupil, class, extrav, sex, texp, popular) # we select just the variables we will use
head(popular2data) # we have a look at the first 6 observations

#---- Plotting the Data ----

ggplot(data    = popular2data,
       aes(x   = extrav,
           y   = popular,
           group=class, # to have the 100 regression lines
           col = class))+ #to add the colours for different classes:  
                          #vizualize the nested multilevel structure of the data
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ #to add some random noise for plotting purposes
  geom_smooth(method = lm,
              se     = FALSE, 
              #col    = "black", # if you remove this line, you have 100 regression lines
              size   = .5, 
              alpha  = .8)+ # to add regression line
  theme_minimal()+
  theme(legend.position = "none")+
  scale_color_gradientn(colours = rainbow(100))+
  labs(title    = "Popularity vs. Extraversion",
       subtitle = "add colours for different classes")

# To colour code the extremes, we need to write a small function that calculates the regression lines and adds a collumn indicating which clusters have the most extreme.
f1 <- function(data, x, y, grouping, n.highest = 3, n.lowest = 3){
  groupinglevel <- data[,grouping]
  res           <- data.frame(coef = rep(NA, length(unique(groupinglevel))), group = unique(groupinglevel))
  names(res)    <- c("coef", grouping)
  for(i in 1:length(unique(groupinglevel))){
    data2    <- as.data.frame(data[data[,grouping] == i,])
    res[i,1] <- as.numeric(lm(data2[, y] ~ data2[, x])$coefficients[2])
  }
  top    <- res %>% top_n(n.highest, coef)
  bottom <- res %>% top_n(-n.lowest, coef)
  res    <- res %>% mutate(high_and_low = ifelse(coef %in% top$coef, "top",  ifelse(coef %in% bottom$coef, "bottom", "none")))
  data3  <- left_join(data, res)
  return(data3)
}

f1(data = as.data.frame(popular2data), 
   x    = "extrav",
   y    = "popular",
   grouping = "class",
   n.highest = 3, 
   n.lowest = 3) %>%
  ggplot()+
  geom_point(aes(x     = extrav,
                 y     = popular, 
                 fill  = class, 
                 group = class),
             size     =  1, 
             alpha    = .5, 
             position = "jitter", 
             shape    = 21, 
             col      = "white")+
  geom_smooth(aes(x     = extrav,
                  y     = popular,
                  col   = high_and_low,
                  group = class,
                  size  = as.factor(high_and_low),
                  alpha = as.factor(high_and_low)),
              method = lm,
              se     = FALSE)+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_fill_gradientn(colours = rainbow(100))+
  scale_color_manual(values=c("top"      = "blue",
                              "bottom"   = "red",
                              "none"     = "grey40"))+
  scale_size_manual(values=c("top"       = 1.2,
                             "bottom"   = 1.2,
                             "none"     = .5))+
  scale_alpha_manual(values=c("top"      = 1,
                              "bottom"    = 1,
                              "none"      =.3))+
  labs(title="Linear Relationship Between Popularity and Extraversion for 100 Classes",
       subtitle="The 6 with the most extreme relationship have been highlighted red and blue")
       

#~~~~ 02 - Building models----
#---- Analysing the Data ----
#---- Intercept only model ----
interceptonlymodeltest <- brm(popular ~ 1 + (1 | class), 
                              # “popular” = dependent variable we want to predict.
                              # |x = random effects/slopes with x as the grouping variables
                              # the dependent variable ‘popular’ is predicted by an intercept and a random error term for the intercept
                              data   = popular2data, 
                              warmup = 100, # how many iterations we want to discard per chain (warmup or burnin phase).
                              iter   = 200, # how many iteration we want the MCMC to run.
                              chains = 2,  # how many chains we want to run.
                              inits  = "random", # what our initial values are for the different chains for the parameters of interest. or we can just tell brms that we want random values as initial values.
                              cores  = 2)  #the cores function tells STAN to make use of 2 CPU cores simultaneously instead of just 1.
# many warnings appear, let's check if the model as converged or not
summary(interceptonlymodeltest)
# the convergence is only partial, so we need to change the prior informations

interceptonlymodel <- brm(popular ~ 1 + (1|class),  
                          data = popular2data, 
                          warmup = 1000, 
                          iter = 3000, 
                          cores = 2, 
                          chains = 2, 
                          seed = 123) #the seed for random number generation to make results reproducible
# notes on seed: A random seed specifies the start point when a computer 
# generates a random number sequence. The purpose of the seed is to allow
# the user to "lock" the pseudo-random number generator, to allow replicable
# analysis.


summary(interceptonlymodel)

# see https://www.rensvandeschoot.com/tutorials/brms-started/ to compute the 
# the posterior mean of the residual variance on the class level
# the residual variance on the first level (pupil level)
# the  Intraclass correlation (ICC) 
# function to do it automatically 
hyp <- "sd_class__Intercept^2 / (sd_class__Intercept^2 + sigma^2) = 0"
hypothesis(interceptonlymodel, hyp, class = NULL)
# This function will also indicate if 0 is included in the 95% CCI of the ICC.
# In our example that is not the case which means a multilevel model is warrante

#---- First Level Predictors ----

model1 <- brm(popular ~ 1 + sex + extrav + (1|class),  
              # ex and extrav are first level predictors, 
              # now included as fixed effect
              data = popular2data, 
              warmup = 1000, iter = 3000, 
              cores = 2, chains = 2, 
              seed = 123) 
# no priors yet, BRMS will pick priors that are non or very weakly informative,
# so that their influence on the results will be negligible.
summary(model1)

# visualize the convergence in so-called caterpillar plots.
model1tranformed <- ggs(model1) 
# the ggs function transforms the brms output into a longformat tibble,
# that we can use to make different types of plots.

ggplot(filter(model1tranformed, 
              Parameter %in% c("b_Intercept", "b_extrav", "b_sex")),
       aes(x   = Iteration,
           y   = value, 
           col = as.factor(Chain)))+
  geom_line() +
  geom_vline(xintercept = 1000)+
  facet_grid(Parameter ~ . ,
             scale  = 'free_y',
             switch = 'y')+
  labs(title = "Caterpillar Plots", 
       col   = "Chains")

# To test whether all regression coefficients are different from zero,
# we can look at the Credible Intervals that are listed in the summary output 
# or we can visually represent them in density plots. If we do so, we clearly 
# see that zero is not included in any of the density plots, meaning that we can 
# be reasonably certain the regression coefficients are different from zero

ggplot(filter(model1tranformed,
              #Parameter == "b_Intercept",  # check one parameter at a time  
              Parameter =="b_extrav",
              #Parameter =="b_sex",
              Iteration > 1000),
       aes(x = value))+
  geom_density(fill  = "yellow", 
               alpha = .5)+
  geom_vline(xintercept = 0, 
             col  = "red",
             size = 1)+
  scale_x_continuous(name   = "Value",
                     limits = c(-1, 3)) + # make clearer the range of the parameter
  geom_vline(xintercept = summary(model1)$fixed[3,3:4], # need to change for the diff parameters
             col = "blue",
             linetype = 2) +
  theme_light() +
  labs(title = "Posterior Density of Intercept") # of Regression Coefficient for Sex 

#---- First and Second Level Predictors ----

model2 <- brm(popular ~ 1 + sex + extrav + texp + (1|class),  
              data = popular2data, 
              warmup = 1000, iter = 3000, 
              cores = 2, chains = 2, 
              seed = 123)

# texp (teacher experience) =  a predictor variable on the second level
# see url to calculate the explained variance at level 1 and at level 2.


#---- First and Second Level Predictors with Random Slopes (1) ----
model3 <- brm(popular ~ 1 + sex + extrav + (1 + sex + extrav | class),  
              data = popular2data, 
              warmup = 1000, iter = 3000, 
              cores = 2, chains = 2, 
              seed = 123) #to run the model
summary(model3)
# seeurl for interpretation and why sex should not be a random factor
#---- First and Second Level Predictors with Random Slopes (2) ----
# We continue after omitting the random slope of sex.
model4 <- brm(popular ~ 1 + sex + extrav + texp + (1 + extrav | class),  
              data = popular2data, 
              warmup = 1000, iter = 3000, 
              cores = 2, chains = 2, 
              seed = 123)
summary(model4)

##---- First and Second Level Predictors with Random Slopes and Crosslevel Interaction
# add the cross-level interaction between Extraversion and Teacher experience
model5 <- brm(popular ~ 1 + sex + extrav + texp + extrav:texp + (1 + extrav|class), 
              data  = popular2data, warmup = 1000,
              iter  = 3000, chains = 2, 
              seed  = 123, control = list(adapt_delta = 0.97),
              cores = 2) 
# to reach a usuable number effective samples in the posterior distribution of
# the interaction effect, we need many more iteration. This sampler will take
# quite some time and you might want to run it with a few less iterations.

summary(model5)$fixed




#################################################################################
# Influence of Priors: Popularity Data
#################################################################################
# source: https://www.rensvandeschoot.com/tutorials/brms-priors/
# 4 types of priors
# 1 - With an estimate far off the value we found in the data with uninformative priors 
# with a wide variance 
# 2 - With an estimate close to the value we found in the data with uninformative priors
# with a small variance 
# 3 - With an estimate far off the value we found in the data with uninformative priors 
# with a small variance 
# 4. With an estimate far off the value we found in the data with uninformative priors 
# with a small variance (2).
# In this tutorial we will only focus on priors for the regression coefficients 
# and not on the error and variance terms, 
# since we are most likely to actually have information on the size and direction
# of a certain effect and less (but not completely) unlikely to have prior knowledge 
# on the unexplained variances. You might have to play around a little bit with the
# controls of the brm() function and specifically the adapt_delta and max_treedepth. 
# Thankfully BRMS will tell you when to do so.

#---- Effect of Priors----
# With the get_prior() command we can see which priors we can specify for this model.
get_prior(popular ~ 0 + intercept + sex + extrav + texp + extrav:texp + (1 + extrav | class),
          data = popular2data)
# To place a prior on the fixed intercept, one needs to include 0 + intercept
prior1 <- c(set_prior("normal(-10,100)", class = "b", coef = "extrav"),
            set_prior("normal(10,100)", class = "b", coef = "extrav:texp"),
            set_prior("normal(-5,100)", class = "b", coef = "sex"),
            set_prior("normal(-5,100)", class = "b", coef = "texp"),
            set_prior("normal(10,100)", class = "b", coef = "intercept" ))
model6 <- brm(popular ~ 0 + intercept + sex + extrav + texp + extrav:texp + (1 + extrav|class), 
             data  = popular2data, warmup = 1000,
             iter  = 3000, chains = 2, 
             prior = prior1,
             seed  = 123, control = list(adapt_delta = 0.97),
             cores = 2,
             sample_prior = TRUE) 
# to reach a usuable number effective samples in the posterior distribution of the interaction effect,
# we need many more iteration. This sampler will take quite some time and you might want to run it with
# a few less iterations.

# To see which priors were inserted, use the prior_summary() command