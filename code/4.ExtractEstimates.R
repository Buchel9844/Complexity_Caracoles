# Extract estimates
#install.packages("stringr")

#---- Generic parameters---
load(paste0(project.dic,"results/stan/FinalFit",FocalPrefix,"_",year,"_",complexity.plant,"_",complexity.animal,".rds")) 
FinalPosteriors <- rstan::extract(FinalFit)

df_parameter_n <-NULL
generic_name <- c("lambdas","alpha_intra","alpha_generic","gamma_H_generic","gamma_FV_generic")
generic_name <- generic_name[generic_name %in%  names(FinalPosteriors)]
for ( n in generic_name ){
  df_parameter_hat_n <- FinalPosteriors[[n]] %>% 
    as.data.frame()
  names( df_parameter_hat_n) <-  n 
  
  df_parameter_n <- bind_cols(df_parameter_n,df_parameter_hat_n)
}

df_parameter_n$focal <- focal
df_parameter_n$year <- year
df_parameter_n$complexity.plant <- complexity.plant
df_parameter_n$complexity.animal <- complexity.animal

write.csv(df_parameter_n,
          paste0(project.dic,"results/parameters/Parameters_",focal,"_", year,"_",complexity.plant,
                 complexity.animal,"_FinalFit.csv"))

# Species specific parameters


load(paste0(home.dic,"results/Inclusion",focal,"_",
            year,"_",complexity.plant,"_",complexity.animal,".RData"))
assign(paste0("Inclusion",focal,"_",
              year,"_",complexity.plant,"_",complexity.animal),Inclusion_all)

rm(Inclusion_all)

df_inclusion_nat <- as.data.frame(get(paste0("Inclusion",focal,"_",
                                             year,"_",complexity.plant,"_",complexity.animal))$interaction) 


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
if(DataVec$RemoveFV == 1){
  n.inclus <- c("n.competitors_plant_inclus","n.competitors_H_inclus",
                "n.HOIs_plant_inclus","n.HOIs_H_inclus")
  inclus <- c("competitors_plant_inclus","competitors_H_inclus",
              "HOIs_plant_inclus","HOIs_H_inclus")
  inclus.short <- grep(pattern = '.*inclus$', x = colnames(df_inclusion_nat), value = T)
  inclus.short <- inclus.short[which(!inclus.short %in% n.inclus)]
  inclus.not.short <- inclus[which(!inclus %in% inclus.short)]
  df_inclusion_nat[,inclus.not.short] <-NA
  
  n.potential <- c("n.competitors_plant","n.interactors_H",
                   "n.HOIs_plant","n.HOIs_H")
  
  potential <- c("competitors_plant","interactors_H",
                 "HOIs_plant","HOIs_H")
  
}
if(DataVec$RemoveH == 1){
  n.inclus <- c("n.competitors_plant_inclus","n.competitors_FV_inclus","n.HOIs_plant_inclus","n.HOIs_FV_inclus")
  inclus <- c("competitors_plant_inclus","competitors_FV_inclus","HOIs_plant_inclus", "HOIs_FV_inclus")
  inclus.short <- grep(pattern = '.*inclus$', x = colnames(df_inclusion_nat), value = T)
  inclus.short <- inclus.short[which(!inclus.short %in% n.inclus)]
  inclus.not.short <- inclus[which(!inclus %in% inclus.short)]
  df_inclusion_nat[,inclus.not.short] <-NA
  
  n.potential <- c("n.competitors_plant","n.interactors_FV", "n.HOIs_plant","n.HOIs_FV")
  
  potential <- c("competitors_plant","interactors_FV","HOIs_plant","HOIs_FV")
}
if(DataVec$RemoveFvH == 1){
  n.inclus <- c("n.competitors_plant_inclus","n.HOIs_plant_inclus")
  inclus <- c("competitors_plant_inclus","HOIs_plant_inclus")
  inclus.short <- grep(pattern = '.*inclus$', x = colnames(df_inclusion_nat), value = T)
  inclus.short <- inclus.short[which(!inclus.short %in% n.inclus)]
  inclus.not.short <- inclus[which(!inclus %in% inclus.short)]
  df_inclusion_nat[,inclus.not.short] <-NA
  
  n.potential <- c("n.competitors_plant","n.HOIs_plant")
  
  potential <- c("competitors_plant","HOIs_plant")
  
}

extract.vert <- function(mat,vec,keyname,value,vecname){
  mat.2 <- mat %>% 
    gather(any_of(vec), 
           key = keyname , value = value) %>%
    mutate(keyname = rep(vecname, each = nrow(mat)))
  return(mat.2)
}

if (sum(colSums(df_inclusion_nat[,grep(pattern = "^n.*\\inclus$", x = colnames(df_inclusion_nat), value = T)],na.rm=T)) > 0) {
  
  df_inclusion_nat <- df_inclusion_nat %>% 
    gather(any_of(potential), 
           key="interaction", value= "potential.interactor") %>%
    #mutate(interaction = rep(inclus, times = nrow(df_inclusion_nat))) %>%
    dplyr::select(-any_of(c(inclus,n.inclus,n.potential))) %>%
    mutate(interactor = extract.vert(df_inclusion_nat,inclus,
                                     "interaction","identity",inclus)[,"value"],
           n.number = extract.vert(df_inclusion_nat,n.inclus,
                                   "interaction.number","n.number",n.inclus)[,"value"],
           n.potential.interactor = extract.vert(df_inclusion_nat,n.potential ,
                                                 "interaction.pot","n.identity.potential",n.potential)[,"value"]) %>%
    mutate(interaction = dplyr::case_when(interaction == "competitors_plant" ~ "alpha_hat_ij",
                                          interaction == "interactors_FV" ~ "gamma_FV_hat_if",
                                          interaction== "interactors_H" ~ "gamma_H_hat_ih",
                                          interaction== "HOIs_plant" ~ "beta_plant_hat_ijk",
                                          interaction== "HOIs_FV" ~ "beta_FV_hat_ijf",
                                          interaction== "HOIs_H" ~ "beta_H_hat_ijh"))
  df_inclusion_nat_relevant <- df_inclusion_nat %>% 
    dplyr::filter(n.number != 0)
  
  # function to extract specific value of a HOIs matrix 
  extract.matrix.HOIs.plant <- function(n.neighbours,position.n) {
    vec.position <- 1:sum(1:n.neighbours-1)
    data.position <- c(rep(0,times=vec.position[1]),vec.position[1]:c(n.neighbours-1))
    for(n in c(2:n.neighbours)){
      v <- c(rep(0,times=vec.position[n]),
             data.position[length(data.position)]+1:vec.position[length(vec.position)])[1:n.neighbours]
      if(n == n.neighbours){v <- c(rep(0,times=vec.position[n]))}
      data.position <- c(data.position,v)
    }
    mat.position <- matrix(nrow=n.neighbours,ncol=n.neighbours, byrow = T,
                           data= data.position)
    position.new <- c()
    for( n in 1:length(position.n)){
      
      vec.position <- 1:n.neighbours
      
      position.new[n]<- paste0(vec.position[apply(mat.position, 1, function(r) any(r %in% c(position.n[n])))],",",# row
                               vec.position[apply(mat.position, 2, function(r) any(r %in% c(position.n[n])))]) # column
    }
    return(position.new)
  }
  
  extract.matrix.HOIs.HTL <- function(n.plants,n.neighbours,position.n) {
    mat.position <- matrix(ncol=n.neighbours,nrow=n.plants,
                           data=c(1:(n.plants*n.neighbours)),
                           byrow = T)
    
    position.new <- c()
    for( n in 1:length(position.n)){
      
      row.position <- 1:n.plants
      col.position <- 1:n.neighbours
      
      position.new[n]<- paste0(row.position[apply(mat.position, 1, function(r) any(r %in% c(position.n[n])))],",",# row
                               col.position[apply(mat.position, 2, function(r) any(r %in% c(position.n[n])))]) # column
    }
    return(position.new)
  }
  
  if(nrow(df_inclusion_nat_relevant) != 0 ){
    df_parameter_hat <- NULL
    for( i in 1:nrow(df_inclusion_nat_relevant)){
      position.n <- match(stringr::str_split_1(df_inclusion_nat_relevant[i,"interactor"],","),
                          stringr::str_split_1(df_inclusion_nat_relevant[i,"potential.interactor"],","))
      name.all <- stringr::str_split_1(df_inclusion_nat_relevant[i,"potential.interactor"],",")
      if(df_inclusion_nat_relevant[i,"interaction"] == c("beta_plant_hat_ijk")){
        n.plants <- DataVec$S
        n.neighbours <- DataVec$S
        position.name <- position.n
        position.n <- extract.matrix.HOIs.plant(n.neighbours, position.n)
      }
      if(df_inclusion_nat_relevant[i,"interaction"] == c("beta_FV_hat_ijf")){
        name.all <-paste0(rep(c(stringr::str_split_1(df_inclusion_nat[which(df_inclusion_nat$interaction == "alpha_hat_ij"),
                                                             "potential.interactor"],","),
                                focal),each=DataVec$FV),
                          "_",
                          stringr::str_split_1(df_inclusion_nat[which(df_inclusion_nat$interaction == "gamma_FV_hat_if"),
                                                       "potential.interactor"],","))
        position.n <- match(stringr::str_split_1(df_inclusion_nat_relevant[i,"interactor"],","),
                            name.all)
        name.all <- stringr::str_split_1(df_inclusion_nat_relevant[i,"potential.interactor"],",")
        n.plants <- DataVec$S
        n.neighbours <- DataVec$FV
        position.name <- position.n
        position.n <- extract.matrix.HOIs.HTL(n.plants,n.neighbours, position.n)
      }
      if(df_inclusion_nat_relevant[i,"interaction"] == c("beta_H_hat_ijh")){
        name.all <- paste0(rep(c(stringr::str_split_1(df_inclusion_nat[which(df_inclusion_nat$interaction == "alpha_hat_ij"),
                                                              "potential.interactor"],","),
                                 focal),each=DataVec$H),
                           "_",
                           stringr::str_split_1(df_inclusion_nat[which(df_inclusion_nat$interaction == "gamma_H_hat_ih"),
                                                        "potential.interactor"],","))
        
        position.n <- match(stringr::str_split_1(df_inclusion_nat_relevant[i,"interactor"],","),
                            name.all)
        n.plants <- DataVec$S
        n.neighbours <- DataVec$H
        position.name <- position.n
        position.n <- extract.matrix.HOIs.HTL(n.plants,n.neighbours, position.n)
      }
      df_parameter_hat_n <- NULL
      for( n in position.n){
        n.2 <- which(position.n == n)
        if(is.numeric(n)){
          df_parameter_hat_n <- FinalPosteriors[[df_inclusion_nat_relevant[i,"interaction"]]][,n] %>% 
            as.data.frame()
          names(df_parameter_hat_n) <- paste0(df_inclusion_nat_relevant[i,"interaction"],"_",name.all[n])
          
        }else{
          n.list <- stringr::str_split_1(n,",")
          df_parameter_hat_n <- FinalPosteriors[[df_inclusion_nat_relevant[i,"interaction"]]][,as.numeric(n.list[1]),as.numeric(n.list[2])] %>% 
            as.data.frame()
          names(df_parameter_hat_n) <- paste0(df_inclusion_nat_relevant[i,"interaction"],"_", name.all[position.name[n.2]])
        }
        df_parameter_hat <- bind_cols(df_parameter_hat, df_parameter_hat_n) 
        
      }
      
    }
    
    df_parameter_hat$focal <- focal
    df_parameter_hat$year <- year
    df_parameter_hat$complexity.plant <- complexity.plant
    df_parameter_hat$complexity.animal <- complexity.animal
    write.csv(df_parameter_hat,
              paste0(project.dic,"results/parameters/Parameters_hat_",focal,"_", year,"_",complexity.plant,
                     complexity.animal,"_FinalFit.csv"))
  }
}
