#---- 1. Create fake data set ----

tau0 <- 1
slab_scale <- sqrt(2)
slab_df <- 4
N=100
S = 5
SpMatrix=matrix(0,nrow=100,ncol=5)

matrix_HOIs_plant <- list()
for (i in 1:N){
  matrix_i <- matrix(nrow=S,ncol=S)
  for (n in 1:S) {
    for (m in 1:S) {
      if (m <= n){
        matrix_i[n,m] = 0
      }else{
        matrix_i[n,m] = SpMatrix[i,n]* SpMatrix[i,m]
      } 
    }
  }
  matrix_i[is.na(matrix_i)] <- 0
  matrix_HOIs_plant[[i]] <-  matrix_i
}

H = 5
SpMatrix_H=matrix(0,nrow=100,ncol=5)
matrix_HOIs_ijh <- matrix_HOIs_plant 

FV = 5
SpMatrix_FV=matrix(0,nrow=100,ncol=5)
matrix_HOIs_ijf <- matrix_HOIs_plant 


DataVec <- list(N=100, 
                SpNames = c(),
                S=5,
                U=1,
                RemoveFvH= 0,
                RemoveH = 0,
                RemoveFV= 0,
                Fecundity=rep(0,100), 
                SpMatrix=matrix(0,nrow=100,ncol=5),
                matrix_HOIs_plant=matrix_HOIs_plant,
                Intra=c(1,0,0,0,0), 
                tau0=tau0, 
                slab_scale=slab_scale, 
                slab_df=slab_df,
                H=H,
                matrix_HOIs_ijh=matrix_HOIs_ijh,
                SpMatrix_H=SpMatrix_H,
                FV=FV,
                SpMatrix_FV=SpMatrix_FV,
                matrix_HOIs_ijf=matrix_HOIs_ijf)
#---- 2. Run  preliminary fit ----

# Now run a preliminary fit of the model to assess parameter shrinkage



print("preliminary fit beginning")
options(mc.cores=parallel::detectCores())

rstan_options(auto_write = TRUE) 
options(mc.cores = parallel::detectCores()) # to use the core at disposition 

PrelimFit <- stan(paste0(home.dic,"code/Short_Caracoles_BH_FH_Preliminary.stan"), 
                    data = DataVec,
                    init = "random", 
                    control =list(max_treedepth=15),
                    warmup  = 1000, 
                    iter = 2000,
                    #init_r=2,
                    chains = 4,
                    seed= 1616)
save(file=paste0(project.dic,"results/stan/PrelimFit_NoData.rds"),
     PrelimFit)

PrelimPosteriors <- rstan::extract(PrelimFit)
print("preliminary fit done")


#---- 3. Preliminary fit posterior check and behavior checks---- 
##### Diagnostic plots
pdf(paste0(home.dic,"figure/Prelimfit_Nodata.pdf"))
title= paste0("Diagnostic plots for Prelimfit of Nodata")
# Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
source("~/Eco_Bayesian/Test_simulation/code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots

# N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
hist(summary(PrelimFit)$summary[,"Rhat"], 
     main = paste("Prelim fit: Histogram of Rhat for Nodata"))
hist(summary(PrelimFit)$summary[,"n_eff"],
     main = paste("Prelim fit: Histogram of n_eff for Nodata"))

# Next check the correlation among key model parameters and identify any
#       divergent transitions
par <- c('lambdas','disp_dev','alpha_generic','alpha_intra',"gamma_FV_generic","gamma_H_generic")

stan_trace(PrelimFit,pars=par,inc_warmup = TRUE)
stan_dens(PrelimFit,pars=par)
stan_plot(PrelimFit,pars=par)

dev.off()

#---- 2. Run  final fit ----

DataVec.final <- list(Inclusion_ij = matrix(1,1,5),
                      Inclusion_FV=matrix(1,1,5),
                      Inclusion_H=matrix(1,1,5), 
                      beta_Inclusion_plant=matrix(1,5,5),
                      beta_Inclusion_FV=matrix(1,5,5),
                      beta_Inclusion_H=matrix(1,5,5),
                      run_estimation = 1)

DataVec.final <- append(DataVec,DataVec.final)

options(mc.cores = parallel::detectCores()) # to use the core at disposition 
FinalFit <- stan(file = paste0(home.dic,"code/Short_Caracoles_BH_Final.stan") ,
                 #fit= PrelimFit, 
                 data = DataVec.final,
                 init ="random", # all initial values are 0 
                 control=list(max_treedepth=15),
                 warmup = 500,
                 iter = 1000, 
                 init_r = 2,
                 chains = 4,
                 seed= 1616) 

save(file=paste0(project.dic,"results/stan/FinalFit_NoData.rds"),
     FinalFit)

#---- 3. Preliminary fit posterior check and behavior checks---- 
##### Diagnostic plots
pdf(paste0(home.dic,"figure/FinalFit_Nodata.pdf"))
title= paste0("Diagnostic plots for Final fit of Nodata")
# Internal checks of the behaviour of the Bayes Modelsummary(PrelimFit)
source("~/Eco_Bayesian/Test_simulation/code/stan_modelcheck_rem.R") # call the functions to check diagnistic plots

# N.B. amount by which autocorrelation within the chains increases uncertainty in estimates can be measured
hist(summary(FinalFit)$summary[,"Rhat"], 
     main = paste("Final fit: Histogram of Rhat for Nodata"))
hist(summary(FinalFit)$summary[,"n_eff"],
     main = paste("Final fit: Histogram of n_eff for Nodata"))

# Next check the correlation among key model parameters and identify any
#       divergent transitions
par <- c('lambdas','disp_dev','alpha_generic','alpha_intra',
         "gamma_FV_generic","gamma_H_generic",
         'alpha_hat_ij[1]',"gamma_H_hat_ih[1]",
         "gamma_FV_hat_if[1]","beta_H_hat_ijh[1,1]",
         "beta_FV_hat_ijf[1,1]")

stan_trace(FinalFit,pars=par,inc_warmup = TRUE)
stan_dens(FinalFit,pars=par)
stan_plot(FinalFit,pars=par)

dev.off()
