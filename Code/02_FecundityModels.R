
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
  fit.alphas=FALSE, 
  fit.gama=FALSE,
  fit.betas.polinator.on.comp =FALSE,
  fit.betas.polinator.on.mutu =FALSE,
  fit.betas.plant.on.comp = FALSE, 
  #fit.betas.plant.on.mutu = FALSE,
  sample.alphas=NULL,
  sample.betas=NULL,
  ...
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
  not.numerical.col <- c("focal","seeds","treatment", "year")
  
  # lets figure out who the potential competitors are
  competitors <- colnames(data)[!colnames(data) %in% not.numerical.col]
  
  # remove species that are never observed to co-occur
  for(sp in names(which(colSums(data[,competitors, drop=FALSE],na.rm = T)==0))){
    data[,sp] <- NULL
  }
  
  # lets figure out the observed competitors 
  competitors <- colnames(data)[!colnames(data) %in% not.numerical.col]
  
  ##########################################################################################################
  # Description of the different coefficient based on the ecological network
  ##########################################################################################################
  # A ---- For pairwise competition----
  # all combinaison of plant-plant is possible 
  all.alphas <- competitors[! competitors %in% pollinator.id]
  
  # B ---- For mutualistic effect----
  # Combinaison of all plant-pollinator possible 
  all.gama <- competitors[!competitors %in% plant.id]
  
  # C ---- For HOI on mutualistic effects-- two pollinators ----
  
  # betas between heterospecific neighbors
  if(length(all.gama)>1){
    interbetas.polinator.on.mutu <- combn(all.gama,2)
    interbetas.polinator.on.mutu.2 <- apply(interbetas.polinator.on.mutu,2,paste,collapse=":")
    
    # combine all betas together into a single variable
    all.betas.polinator.on.mutu <- c(interbetas.polinator.on.mutu.2)
  }
  else{all.betas.polinator.on.mutu <-  c() } # c() need to be replaced by intrabetas.polinator.on.mutu if soft HOI included
  
  
  # D ---- For HOI on Competition----
  ### fit.betas.polinator.on.comp ###
  betas.polinator.on.comp <- as.vector(outer(all.alphas,all.gama, paste, sep=":"))
  imp.interaction <- c("Tomato:Lucilia_sericata","Tomato:Osmia_bicornis","Vicia:Lucilia_sericata")
  imp.interaction.no.link <- c("Raphanus:Bombus_terrestris")
  all.betas.polinator.on.comp <- betas.polinator.on.comp[!betas.polinator.on.comp %in% imp.interaction]
  if (network == "no_link"){
    all.betas.polinator.on.comp <- all.betas.polinator.on.comp[!all.betas.polinator.on.comp %in% imp.interaction.no.link ]
  }
  ### fit.betas.plant.on.comp ###
  
  # betas between heterospecific neighbors
  if(length(all.alphas)>1){
    interbetas.plant.on.comp <- combn(all.alphas,2)
    interbetas.plant.on.comp <- apply(interbetas.plant.on.comp,2,paste,collapse=":")
    # combine all betas together into a single variable
    all.betas.plant.on.comp <- c(interbetas.plant.on.comp)
  }
  else{all.betas.plant.on.comp <- c()}
  
  #all.betas <- c()
  #all.betas <-c(all.betas.plant.on.comp,all.betas.polinator.on.comp,all.betas.polinator.on.mutu)
  
  ##########################################################################################################
  # Description of the different situations (panels)
  ##########################################################################################################   
  # start with no competition model
  base.model.formula <- "seeds ~ 1"
  
  # add a fixed effect 
  if(nlevels(data$year)>1){
    base.model.formula <- paste0(base.model.formula,"+ year")
  }
  
  # in case it's the null
  model.formula <- base.model.formula
  
  
  
  # A ---- add pairwise coefficients for all competitor species-----
  
  if((fit.alphas) && (length(competitors) > 0)){
    # if there are alphas to add, add them to the model formula
    if(length(all.alphas)>0){
      model.formula <- paste0(model.formula, " + ",paste0(all.alphas, collapse=" + "))
    }
  }
  
  # B ---- add mutualistic effect between plant and pollinator, gama ----
  
  if((fit.gama) && (length(competitors) > 0)){
    if(length(all.gama) >0){
      model.formula <- paste0(model.formula, " + ", paste0(all.gama, collapse=" + "))
      
    }
  }
  
  # C ---- add HOI on mutualistic effects----
  ### fit.betas.polinator.on.mutu ###
  if(fit.betas.polinator.on.mutu && (length(competitors) > 0)){
    # add the betas to the model formula
    if(length(all.betas.polinator.on.mutu )>0){
      model.formula <- paste0(model.formula, " + ", paste0(all.betas.polinator.on.mutu, collapse=" + "))
    }
  }    
  
  # D ---- add HOI on Competition----
  ### fit.betas.polinator.on.comp ###
  if(fit.betas.polinator.on.comp  && (length(competitors) > 0)){
    # add the betas to the model formula
    if(length(all.betas.polinator.on.comp)>0){
      model.formula <- paste0(model.formula, " + ", paste0(all.betas.polinator.on.comp, collapse=" + "))
    }
  }   
  
  ### fit.betas.plant.on.comp ###
  if(fit.betas.plant.on.comp  && (length(competitors) > 0)){
    # add the betas to the model formula
    if(length(all.betas.plant.on.comp)>0){
      model.formula <- paste0(model.formula, " + ", paste0(all.betas.plant.on.comp, collapse=" + "))
    }
  } 
  
  # E ---- Remove interaction equal to 0---- 
  #X <- as.data.frame(model.matrix(as.formula(model.formula), data))
  #Xnames <- names(X)
  # if ( length(Xnames) > 2){
  #    Xremove <- c()
  #    if (length(names(which(colSums(X,na.rm=T) == 0 ))) != 0){
  #        Xremove <- c(names(which(colSums(X,na.rm=T) == 0)),"(Intercept)","year2017")
  #    }
  #    else{ Xremove <- c("(Intercept)","year2017")
  #    }
  #    all.coeff <- c()
  #   all.coeff <- Xnames[which(!Xnames %in%Xremove)]
  #    if (length(all.coeff)>0) {
  #        model.formula <- paste0(base.model.formula, " + ", paste0(all.coeff, collapse=" + "))}
  #}
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
      control=glm.control(maxit=1000),
      ...
    )
  }
  
  # fit a poisson form which is equivalent to the negative binomial without the dispersion parameter
  if(type=='poisson'){
    m <- glm(
      as.formula(model.formula),
      data=data, 
      family='poisson',
      control=glm.control(maxit=1000),
      ...
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
  # mm.base <- colnames(mm)[which(!colnames(mm) %in% all.alphas & !colnames(mm) %in% all.betas)]
  
  # figure out which columns (i.e. coefficients) to remove to eliminate rank deficiencies using qr decomposition
  # qrmm <- qr(mm)
  
  # separate out the names of these predictors
  #reduced.predictors <- colnames(mm)[qrmm$pivot[seq(qrmm$rank)]]
  
  # remove anything that is included in the base model so that the model formula command below makes sense
  # reduced.predictors <- reduced.predictors[which(!reduced.predictors %in% mm.base)]
  
  # separate these into alphas and betas
  #reduced.alphas <- reduced.predictors[which(!reduced.predictors %in% all.betas)]
  #reduced.betas <- reduced.predictors[which(reduced.predictors %in% all.betas)]
  
  # put them back together in a sorted order (so that coefficients in the the model summaries resemble what R would normally produce)
  # reduced.predictors <- c(sort(reduced.alphas), sort(reduced.betas))
  
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
      random = as.formula(~ (1/"random.variable")),
      data=data,
      family='nbinom',
      ...
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

