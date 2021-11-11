
#### Fit_fecundity_model, modified from Mayfield and Stouffer, 2018
# the model 'negbin-mm' used the factor year as random effect 
# gamma
#pollinator
################################################################################################
# Description of the function
##############################################################################################

fit.fecundity.model <- function(
  data,
  type=c("negbin", "negbin-mm", "poisson", "linear","inverse"),
  fit.plot= FALSE,
  fit.year = FALSE,
  fit.alpha=FALSE, 
  fit.gama=FALSE,
  fit.betas.plant.plant = FALSE, 
  fit.betas.plant.HTL = FALSE, 
  fit.betas.HTL.HTL = FALSE
){
  #  ---- Preliminary operations ----  
  # make sure the user has provided an eligible model
  match.arg(type)
  
  # some of the models require differzent packages for the fits to be performed and/or to converge		
  if(type == 'negbin') library(MASS)
  if(type == 'negbin-mm') library(glmmADMB)
  if(type == 'inverse') library(glm2)
  
  # drop extraneous levels that could complicate the model-fitting code
  data$year <- as.factor(data$year)
  data$year <- droplevels(data$year)
  
   # remove species that are never observed to co-occur
  for(sp in c(id.plant.comp,id.int.type)){
  data[which(is.na(data[,sp])),sp] <- 0 # replace NA by 0
  if (sum(data[,sp],na.rm = T)==0){
    data[,sp] <- NULL
  }
  }
  
  # lets figure out the observed competitors 
  competitors.plant <- colnames(data)[colnames(data) %in% id.plant.comp]
  competitors.HTL <- colnames(data)[colnames(data) %in% id.int.type]
  
  ##########################################################################################################
  # Description of the different coefficient based on the ecological network
  ##########################################################################################################
  # A ---- For pairwise competition----
  # all combinaison of plant-plant is possible 
  all.alpha <- competitors.plant
  
  # B ---- For mutualistic effect----
  # Combinaison of all plant-HTL(insect, pollinator,predator) possible 
  all.gama <- competitors.HTL
  
  # C ---- For HOI on mutualistic effects-- two HTL ----
  
  # betas between heterospecific neighbors
  if(length(all.gama)>1){
    betas.HTL.HTL <- combn(all.gama,2)
    betas.HTL.HTL <- apply(betas.HTL.HTL,2,paste,collapse=":")
    
    }else{betas.HTL.HTL <-  c() } # c() need to be replaced by intrabetas.polinator.on.mutu if soft HOI included
  
  
  # D ---- For HOI on Competition----
  ### fit betas between plant and HTL ###
  betas.plant.HTL <- as.vector(outer(all.alpha,all.gama, paste, sep=":"))
 
  ### fit betas between heterospecific neighbors ###
  if(length(all.alpha)>1){
    betas.plant.plant <- combn(all.alpha,2)
    betas.plant.plant <- apply(betas.plant.plant,2,paste,collapse=":")
    # combine all betas together into a single variable
  }else{betas.plant.plant <- c()}
  
  #all.betas <- c()

  ##########################################################################################################
  # Description of the different situations (panels)
  ##########################################################################################################   
  # start with no competition model
  base.model.formula <- "seed ~ 1"
  
  # add years as a fixed effect 
  if((fit.years) > 0 & nlevels(data$year)>1){
    base.model.formula <- paste0(base.model.formula,"+ year")
  }
  
  # add plot as a random effect 
  if((fit.plot) > 0 & nlevels(data$plot)>1){
    base.model.formula <- paste0(base.model.formula,"+ 1/plot")
  }
  
  
  # in case it's the null
  model.formula <- base.model.formula
  
  # A ---- add pairwise coefficients for all competitor species-----
  
  if((fit.alpha) && (length(competitors.plant) > 0)){
    # if there are alpha to add, add them to the model formula
    if(length(all.alpha)>0){
      model.formula <- paste0(model.formula, " + ",paste0(all.alpha, collapse=" + "))
    }
  }
  
  # B ---- add direct (mutualistic) effect between plant and pollinator, gama ----
  
  if((fit.gama) && (length(competitors.plant) > 0)){
    if(length(all.gama) >0){
      model.formula <- paste0(model.formula, " + ", paste0(all.gama, collapse=" + "))
      
    }
  }
  
  # C ---- add HOI on mutualistic effects----
  ### fit.betas.HTL.HTL ###
  if(fit.betas.HTL.HTL && (length(competitors.plant) > 0)){
    # add the betas to the model formula
    if(length(betas.HTL.HTL )>0){
      model.formula <- paste0(model.formula, " + ", paste0( betas.HTL.HTL, collapse=" + "))
    }
  }    
  
  # D ---- add HOI on Competition----
  ### fit.betas.plant.HT ###
  if(fit.betas.plant.HTL  && (length(competitors.plant) > 0)){
    # add the betas to the model formula
    if(length(betas.plant.HTL)>0){
      model.formula <- paste0(model.formula, " + ", paste0( betas.plant.HTL, collapse=" + "))
    }
  }   
  
  ### fit.betas.plant.plant ###
  if(fit.betas.plant.plant  && (length(competitors.plant) > 0)){
    # add the betas to the model formula
    if(length(betas.plant.plant)>0){
      model.formula <- paste0(model.formula, " + ", paste0(betas.plant.plant, collapse=" + "))
    }
  } 
  
  # E ---- Remove interaction equal to 0---- 
  X <- as.data.frame(model.matrix(as.formula(model.formula), data))
  Xnames <- names(X)
   if ( length(Xnames) > 2){
      Xremove <- c()
      if (length(names(which(colSums(X,na.rm=T) == 0 ))) != 0){
          Xremove <- c(names(which(colSums(X,na.rm=T) == 0)),
  "(Intercept)","year2016","year2017","year2018","year2019"  )
      }
      else{ Xremove <- c("(Intercept)","year2017")
      }
      all.coeff <- c()
     all.coeff <- Xnames[which(!Xnames %in%Xremove)]
      if (length(all.coeff)>0) {
          model.formula <- paste0(base.model.formula, " + ", paste0(all.coeff, collapse=" + "))}
  }
  ##########################################################################################################
  # fit the model
  # note that the call depends on the model form hence the many ifs (and they are not nested for clarity and because they don't need to be)
  ########################################################################################################## ##################

  # fit the negative binomial model as described in the main text
  if(type=='negbin'){
    m <- glm.nb(
      as.formula(model.formula),
      data=data,
      na.action=na.fail,
      control=glm.control(maxit=1000)
    )
  }

    # fit a poisson form which is equivalent to the negative binomial without the dispersion parameter
  if(type=='poisson'){
    m <- glm(
      as.formula(model.formula),
      data=data, 
      family='poisson',
      control=glm.control(maxit=1000)
    )
  }
  
  
  
  
  # fit the linear form
  if(type=="linear"){
    m <- glm(
      as.formula(model.formula),
      data=data,
      control=glm.control(maxit=1000),
      ...
    )
  }
  
  # note that the mixed model and inverse forms are more finicky and harder to get to converge
  # the methods used to fit them also require a non-singular model matrix which is an issue for the beta model (that has many coefficients that cannot be inferred)
  
  # if(type=="negbin-mm" || type=="inverse"){
  # identify linearly independent predictors in this model
  # mm <- model.matrix(as.formula(model.formula), data)
  
  # separate out the terms that don't have to do with per capita effects
  # mm.base <- colnames(mm)[which(!colnames(mm) %in% all.alpha & !colnames(mm) %in% all.gama)]
  
  # figure out which columns (i.e. coefficients) to remove to eliminate rank deficiencies using qr decomposition
  # qrmm <- qr(mm)
  
  # separate out the names of these predictors
  #reduced.predictors <- colnames(mm)[qrmm$pivot[seq(qrmm$rank)]]
  
  # remove anything that is included in the base model so that the model formula command below makes sense
  # reduced.predictors <- reduced.predictors[which(!reduced.predictors %in% mm.base)]
  
  # separate these into alpha and betas
  #reduced.alpha <- reduced.predictors[which(!reduced.predictors %in% all.gama)]
  #reduced.betas <- reduced.predictors[which(reduced.predictors %in% all.gama)]
  
  # put them back together in a sorted order (so that coefficients in the the model summaries resemble what R would normally produce)
  # reduced.predictors <- c(sort(reduced.alpha), sort(reduced.betas))
  
  # reconstruct a non-rank-deficient model formula
  #if(length(reduced.predictors) > 0){
  #    model.formula <- paste0(base.model.formula," + ", paste0(reduced.predictors, collapse=" + "))
  #}else{
  #    model.formula <- base.model.formula
  #}
  #}
  # fit the negative binomial mixed model (with year level random effect)
  if(type=='negbin-mm'){
    m <- glmmadmb(
      as.formula(model.formula),
      #random = as.formula(~ (1/"plot")),
      data=data,
      family='nbinom'
    )
  }
  
  # fit the inverse form
  if(type=='inverse'){
    m<-glm2(
      as.formula(model.formula),
      data=data,
      family=gaussian(link="inverse"),
      control=glm.control(maxit=1000),
      ...
    )
  }
  
  return(m)
}

