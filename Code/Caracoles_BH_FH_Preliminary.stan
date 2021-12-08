// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

data{
  int<lower = 1> N; // Number of plots
  int<lower = 1> S; // Number of plant species
  int<lower = 1> H; // Number of herbivores species
  int<lower = 1> FV; // Number of floral visitors species
  int<lower = 1> H_comp; // Number of HOIs with herbivores species
  int<lower = 1> FV_comp; // Number of HOIs with floral visitors species
  int<lower = 1> FV_FV;// Number of HOIs with 2 floral visitors species
  int<lower = 1> FV_H;// Number of HOIs with floral visitors species and herbivores
  int<lower = 1> H_H;// Number of HOIs with 2 herbivores species
  vector[S] vector_floralvis_comp // Specific vector of order for the HOIs_ijf
  vector[S] vector_herbivore_comp// Specific vector of order for the HOIs_ijh
  vector[FV] vector_FV_herbivore// Specific vector of order for the HOIs_ifh
  vector[H] vector_H_herbivore// Specific vector of order for the HOIs_ihh
  vector[FV] vector_FV_FV// Specific vector of order for the HOIs_iff
  
  int Fecundity[N]; // Fecundity of the focal species in each plot
  //int plot[N];   // Indicator variable for plot
  matrix[N,S] SpMatrix; // Matrix of abundances for each plant species 
  matrix[N,H] SpMatrix_herbivore; // Matrix of abundances for each herbivores species 
  matrix[N,FV] SpMatrix_floralvis; // Matrix of abundances for each floral visitors species

  matrix[N,H_comp] SpMatrix_herbivore_comp; // Matrix of abundances for each herbivores species and competitor
  matrix[N,FV_comp] SpMatrix_floralvis_comp; // Matrix of abundances for each herbivores species and competitor
  
  matrix[N,FV_FV] SpMatrix_FV_FV; // Matrix of abundances for each floral visitors species with each floral visitors species
  matrix[N,H_H] SpMatrix_H_herbivore; // Matrix of abundances for each herbivores species with each herbivores species
  matrix[N,FV_H] SpMatrix_FV_herbivore; // Matrix of abundances for each floral visitors species with each herbivores species

  
  // vector[N] env;   // Environmental values for each plot
  int<lower = 0> Intra[S]; // Indicator boolean variable to identify the focal species (0 for non-focal and 1 for focal). Included for easier calculations
  // The below values define the regularized horseshoe priors used for species-specific parameters
  real tau0; 		// determines the scale of the global shrinkage parameter (tau)
  real slab_scale;	// scale for significant alpha_sp values
  real slab_df;		// effective degrees of freedom for significant alpha_sp values

}

transformed data{
  real slab_scale2 = square(slab_scale);
  real half_slab_df = 0.5*slab_df;
}

parameters{
  vector[1] lambdas;
  vector[1] alpha_generic_tilde; // direct interaction inter plants - plants ; generic
  vector[1] alpha_intra_tilde; // direct interaction intra plants - plants;
  vector[S] alpha_hat_ij_tilde; // direct interaction inter plants - plants ; species -specific term 
  vector<lower = 0>[S] local_shrinkage_ij; // direct interaction inter plants - plants, shrinkage
   
  vector[1] gamma_H_generic_tilde; // direct interaction plants - herbivores ; generic
  vector[H] gamma_H_hat_ih_tilde; // direct interaction plants - herbivores ; species -specific term
  vector<lower = 0>[H] gamma_H_local_shrinkage_ih; // direct interaction plants - herbivores ; shrinkage
  
  vector[1] gamma_FV_generic_tilde; // direct interaction plants - FV ; generic
  vector[FV] gamma_FV_hat_if_tilde; // direct interaction plants - FV ; species -specific term
  vector<lower = 0>[FV] gamma_FV_local_shrinkage_if; // direct interaction plants - FV ; shrinkage
  
  vector[1] beta_plant_generic_tilde; // HOIs plants
  matrix[S,S] beta_plant_hat_ijk_tilde; // HOIs plants
  matrix<lower = 0>[S,S] beta_plant_local_shrinkage_ijk; // HOIs plants
  
  vector[1] beta_H_generic_tilde; // HOIs herbivores
  matrix[S,H] beta_H_hat_ijh_tilde; // HOIs herbivores 
  matrix<lower = 0>[S,H] beta_H_local_shrinkage_ijh; // HOIs herbivores
  
  vector[1] beta_FV_generic_tilde; // HOIs floral visitors
  matrix[S,FV] beta_FV_hat_ijf_tilde; // HOIs floral visitors
  matrix<lower = 0>[S,FV] beta_FV_local_shrinkage_ijf; // HOIs floral visitors
  
  vector[1] beta_2FV_generic_tilde; // HOIs 2 floral visitors
  matrix[FV,FV] beta_2FV_hat_iff_tilde; // HOIs 2 floral visitors
  matrix<lower = 0>[FV,FV] beta_2FV_local_shrinkage_iff; // HOIs 2 floral visitors
  
  vector[1] beta_2H_generic_tilde; // HOIs 2 herbivores
  matrix[H,H] beta_2H_hat_ihh_tilde; // HOIs 2 herbivores 
  matrix<lower = 0>[H,H] beta_2H_local_shrinkage_ihh; // HOIs 2 herbivores
  
  vector[1] beta_FvH_generic_tilde; // HOIs 1 floral visitor & 1 herbivore
  matrix[FV,H] beta_FvH_hat_ifh_tilde; // HOIs 1 floral visitor & 1 herbivore
  matrix<lower = 0>[FV,H] beta_FvH_local_shrinkage_ifh; // HOIs 1 floral visitor & 1 herbivore
  
  
  real<lower = 0> c2_tilde;
  real<lower = 0> tau_tilde;
}

transformed parameters{
  // Calculate the scaled parameters needed for the regularized horeshoe prior here from the normalized (and thus easier to sample)
  // 	counterparts declared in the parameters block
  real c2;
  real tau;
    vector[S] alpha_hat_ij; // direct interaction inter plants
  vector[S] local_shrinkage_ij_tilde; // direct interaction inter plants
    vector[H] gamma_H_hat_ih;  // direct interaction plants - H
  vector[H] gamma_H_local_shrinkage_ih_tilde;  // direct interaction plants - H
    vector[FV] gamma_FV_hat_if;  // direct interaction plants - FV
  vector[FV] gamma_FV_local_shrinkage_if_tilde;  // direct interaction plants - FV

  
  matrix[S,S] beta_plant_hat_ijk; // HOIs plant
  matrix[S,S] beta_plant_local_shrinkage_ijk_tilde; // HOIs plant
  
   matrix[S,H] beta_H_hat_ijh; // HOIs herbivores
  matrix[S,H] beta_H_local_shrinkage_ijh_tilde;// HOIs herbivores
   matrix[S,FV] beta_FV_hat_ijf; // HOIs floral visitors
  matrix[S,FV] beta_FV_local_shrinkage_ijf_tilde; // HOIs floral visitors
  
    matrix[FV,FV] beta_2FV_hat_iff; // HOIs 2 floral visitors
  matrix[FV,FV] beta_2FV_local_shrinkage_iff_tilde; // HOIs 2 floral visitors
  
    matrix[H,H] beta_2H_hat_ihh; // HOIs 2 herbivores
  matrix[H,H] beta_2H_local_shrinkage_ihh_tilde; // HOIs 2 herbivores
 
    matrix[FV,H] beta_FvH_hat_ifh; // HOIs 1 floral visitor & 1 herbivore
  matrix[FV,H] beta_FvH_local_shrinkage_ifh_tilde; // HOIs 1 floral visitor & 1 herbivore
 
  
  vector[1] alpha_generic;  // direct interaction INTER plants
  vector[1] alpha_intra;  // direct interaction INTRA plants
  vector[1] gamma_H;  // direct interaction plants - H
  vector[1] gamma_FV;  // direct interaction plants - FV
  vector[1] beta_plant_generic; //HOIs plant 
  vector[1] beta_H_generic; //HOIs herbivore
  vector[1] beta_FV_generic; //HOIs FV
  vector[1] beta_2FV_generic; //HOIs 2 floral visitors
  vector[1] beta_2H_generic; //HOIs 2 herbivores
  vector[1] beta_FvH_generic; //HOIs 1 floral visitor & 1 herbivore


  tau = tau0*tau_tilde; 	// tau ~ cauchy(0, tau0)
  c2 = slab_scale2*c2_tilde;	// c2 ~ inv_gamma(half_slab_df, half_slab_df*slab_scale2)

  // This calculation follows equation 2.8 in Piironen and Vehtari 2013

    for(s in 1:S){ // direct interaction INTER plants
      local_shrinkage_ij_tilde[s] = sqrt( c2 * square(local_shrinkage_ij[s]) / (c2 + square(tau) * square(local_shrinkage_ij[s])) );
      alpha_hat_ij[s] = tau * local_shrinkage_ij_tilde[s] * alpha_hat_ij_tilde[s];


    for(k in 1:S){ //HOIs plant-plant 
      beta_local_shrinkage_ijk_tilde[s,k] = sqrt( c2 * square(beta_local_shrinkage_ijk[s,k]) / (c2 + square(tau) * square(beta_local_shrinkage_ijk[s,k])) );
      beta_hat_ijk[s,k] = tau * beta_local_shrinkage_ijk_tilde[s,k] * beta_hat_ijk_tilde[s,k];

       }
       for(fv in 1:FV){ //HOIs plant-plant-FV 
      beta_FV_local_shrinkage_ijf_tilde[s,fv] = sqrt( c2 * square(beta_FV_local_shrinkage_ijf[s,fv]) / (c2 + square(tau) * square(beta_FV_local_shrinkage_ijf[s,fv])) );
      beta_FV_hat_ijf[s,fv] = tau * beta_FV_local_shrinkage_ijf_tilde[s,fv] * beta_FV_hat_ijf_tilde[s,fv];

       }
       for(h in 1:H){ //HOIs plant-plant-H
      beta_H_local_shrinkage_ijh_tilde[s,h] = sqrt( c2 * square(beta_H_local_shrinkage_ijh[s,h]) / (c2 + square(tau) * square(beta_H_local_shrinkage_ijh[s,h])) );
      beta_H_hat_ijh[s,h] = tau * beta_H_local_shrinkage_ijh_tilde[s,h] * beta_H_hat_ijh_tilde[s,h];
       }
    }
    
    for(h in 1:H){ # direct interaction plant - H
      gamma_H_local_shrinkage_ih_tilde[h] = sqrt( c2 * square(gamma_H_local_shrinkage_ih[h]) / (c2 + square(tau) * square(gamma_H_local_shrinkage_ih[h])) );
      gamma_H_hat_ih[h] = tau * gamma_H_local_shrinkage_ih_tilde[h] * gamma_H_hat_ih_tilde[h];
      
      for(h2 in 1:H){ //HOIs plant-H-H 
      beta_2H_local_shrinkage_ihh_tilde[h,h2] = sqrt( c2 * square(beta_2H_local_shrinkage_ihh[h,h2]) / (c2 + square(tau) * square(beta_2H_local_shrinkage_ihh[h,h2])) );
      beta_2H_hat_ihh[h,h2] = tau * beta_2H_local_shrinkage_ihh_tilde[h,h2] * beta_2H_hat_ihh_tilde[h,h2];
       }

}

    for(fv in 1:FV){ # direct interaction plant - FV
      gamma_FV_local_shrinkage_if_tilde[fv] = sqrt( c2 * square(gamma_FV_local_shrinkage_if[fv]) / (c2 + square(tau) * square(gamma_FV_local_shrinkage_if[fv])));
      gamma_FV_hat_if[fv] = tau * gamma_FV_local_shrinkage_if_tilde[fv] * gamma_FV_hat_if_tilde[fv];
      
      for(fv2 in 1:FV){ //HOIs plant-FV-FV 
      beta_2FV_local_shrinkage_iff_tilde[fv,fv2] = sqrt( c2 * square(beta_2FV_local_shrinkage_iff[fv,fv2]) / (c2 + square(tau) * square(beta_2FV_local_shrinkage_ifh[fv,fv2])) );
      beta_2FV_hat_iff[fv,fv2] = tau * beta_2FV_local_shrinkage_iff_tilde[fv,fv2] * beta_2FV_hat_iff_tilde[fv,fv2];

       }
       
      for(h in 1:H){ //HOIs plant-FV-H 
      beta_FvH_local_shrinkage_ifh_tilde[fv,h] = sqrt( c2 * square(beta_FvH_local_shrinkage_ifh[fv,h]) / (c2 + square(tau) * square(beta_FvH_local_shrinkage_ifh[fv,h])) );
      beta_FvH_hat_ifh[fv,h] = tau * beta_FvH_local_shrinkage_ifh_tilde[fv,h] * beta_FvH_hat_ifh_tilde[fv,h];

       }
}


  // scale the lambdas,alphas,gammas and HOIs values
  alpha_generic[1] = 3 * alpha_generic_tilde[1] - 6;

  alpha_intra[1] = 3 * alpha_intra_tilde[1] - 6;

  
  gamma_H[1] = 3 * gamma_H_generic_tilde[1] - 6;

  gamma_FV[1] = 3 * gamma_FV_generic_tilde[1] - 6;

  beta_plant_generic[1] = 3 * beta_generic_tilde[1] - 6;

  beta_H_generic[1] = 3 * beta_H_generic_tilde[1] - 6;

  beta_FV_generic[1] = 3 * beta_FV_generic_tilde[1] - 6;

  beta_2FV_generic[1] = 3 * beta_2FV_generic_tilde[1] - 6;

  beta_2H_generic[1] = 3 * beta_2H_generic_tilde[1] - 6;

  beta_FvH_generic[1] = 3 * beta_FvH_generic_tilde[1] - 6;

}

model{
  // Declare objects necessary for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //     a matrix of the species specific alpha values for each species and plot (interaction_effects), and a matrix
  //     of the the alpha*N values for each species.
  vector[N] F_hat;
  vector[N] D_effects;
  vector[N] HOI_effects;
  
  vector[N] lambda_i; #lambda of focal i for obs s
  
  matrix[S,S] matrix_HOIs;
  matrix[N,S] alpha_ij; #alpha of focal i
  matrix[N,H] gamma_H_ih; #gamma herbivores of focal i
  matrix[N,FV] gamma_FV_if; #gamma floral visitors of focal i
  
  # the matrix of interaction of plant-plant-plant HOIs was not created prior to the script
  matrix[S,S] beta_ijk; # HOIs plant 
  matrix[N,S] matrix_beta_ijk; # HOIs plant 
  matrix[1,S] Spmatrixi; # HOIs plant 
  matrix[1,S] Spmatrixi_HTL;

  matrix[S,FV] beta_F_ijf;
    matrix[N,S] matrix_beta_F_ijf # HOIs 1 floral visitor
  matrix[S,H] beta_H_ijh; 
    matrix[N,S] matrix_beta_H_ijh # HOIs 1 herbivore
  
  matrix[FV,FV] beta_2F_iff;
    matrix[N,FV] matrix_beta_2F_iff # HOIs 2 floral visitors
  matrix[H,H] beta_2H_ihh; 
    matrix[N,H] matrix_beta_2H_ihh # HOIs 2 herbivores
  matrix[FV,H] beta_FvH_ifh; 
    matrix[N,FV] matrix_beta_FvH_ifh # HOIs 1 floral visitor and 1 herbivore
  
  // set regular priors
  alpha_generic_tilde ~ normal(0,1);
  alpha_intra_tilde ~ normal(0,1);
  
  gamma_H_generic_tilde ~ normal(0,1);
  gamma_FV_generic_tilde ~ normal(0,1);
  
  beta_generic_tilde ~ normal(0,1);
  beta_H_generic_tilde ~ normal(0,1);
  beta_FV_generic_tilde ~ normal(0,1);
  beta_2FV_generic_tilde ~ normal(0,1);
  beta_2H_generic_tilde ~ normal(0,1);
  beta_FvH_generic_tilde ~ normal(0,1);
  lambdas ~ normal(0, 1);


  // set the hierarchical priors for the Finnish horseshoe (regularized horseshoe) (Piironen and Vehtari 2017)
  // Following the stan implementation from https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html
    tau_tilde ~ cauchy(0,1);
    c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);
    
    alpha_hat_ij_tilde ~ normal(0,1);
    local_shrinkage_ij ~ cauchy(0,1);
    
    gamma_H_hat_ih_tilde ~ normal(0,1);
    gamma_H_local_shrinkage_ih ~ cauchy(0,1);
    
    gamma_FV_hat_if_tilde ~ normal(0,1); 
    gamma_FV_local_shrinkage_if ~ cauchy(0,1);

  
for (s in 1:S){
    beta_hat_ijk_tilde[s,] ~ normal(0,1); 
    beta_local_shrinkage_ijk[s,] ~ cauchy(0,1);
    
        beta_FV_hat_ijf_tilde[s,] ~ normal(0,1);
    beta_FV_local_shrinkage_ijf[s,] ~ cauchy(0,1);
    
    beta_H_hat_ijh_tilde[s,] ~ normal(0,1);
    beta_H_local_shrinkage_ijh[s,] ~ cauchy(0,1);
    

}
for (h in 1:H){
      beta_2H_hat_ihh_tilde[h,] ~ normal(0,1);
    beta_2H_local_shrinkage_ihh[h,] ~ cauchy(0,1);
    }
for (fv in 1:FV){
      beta_2FV_hat_iff_tilde[fv,] ~ normal(0,1);
    beta_2FV_local_shrinkage_iff[fv,] ~ cauchy(0,1);
    
    beta_FvH_hat_ifh_tilde[fv,] ~ normal(0,1);
    beta_FvH_local_shrinkage_ifh[fv,] ~ cauchy(0,1);
}


  // implement the biological model
  for(i in 1:N){ // for the observation N 
    lambda_i[i] = exp(lambdas[1]);
    
    // creation of a matrix of S by S of the interaction jk in HOIs_ijk for plants
    
    Spmatrixi[1,] = SpMatrix[i,]; 
    for (n in 1:S) {
          for (m in 1:S) {
            if (m <= n){
              matrix_HOIs[n,m] = 0;
            }
            else{
        matrix_HOIs[n,m] = Spmatrixi[1,n]* Spmatrixi[1,m];
            } 
        }
        }
        
#---- Interaction (HOIs) matrix of herbivores with COMPETITORS ----

matrix[1,H_comp] Spmatrixi_HTL_H
matrix[S,H] matrix_HOIs_ijh

Spmatrixi_HTL_H <- NULL
matrix_HOIs_ijh <- matrix(ncol=H,nrow=H)
Spmatrixi_HTL_H = SpMatrix_herbivore_comp[i,]; 
n = 1
m=1
for (h_comp in 1:H_comp) {
  if(h_comp == vector_herbivore_comp[n] + 1){
    n = n+1
    m = 1 }
  matrix_HOIs_ijh[n,m] = Spmatrixi_HTL_H[h_comp];
  m = m+1 
}

if(sum(rowSums(is.na(matrix_HOIs_ijh))==FV)){
  print("One or multiple row(s) in matrix_HOIs_ijh contain only NAs")
}
#---- Interaction (HOIs) matrix of floral visitors with COMPETITORS ----
matrix[1,FV_comp] Spmatrixi_HTL_FV
matrix[S,FV] matrix_HOIs_ijf

Spmatrixi_HTL_FV <- NULL
matrix_HOIs_ijf <- matrix(ncol=FV,nrow=S)
Spmatrixi_HTL_FV = SpMatrix_floralvis_comp[i,]; 
n = 1
m=1

for (fv in 1:FV_comp) {
  if(fv == vector_floralvis_comp[n] + 1){
    n = n+1
    m = 1 }
  matrix_HOIs_ijf[n,m] = Spmatrixi_HTL_FV[fv];
  m = m+1 
}

if(sum(rowSums(is.na(matrix_HOIs_ijf))==T)){
  print("One or multiple row(s) in matrix_HOIs_ijf contain only NAs")
}
#---- Interaction (HOIs) matrix of herbivores with herbivores ----

matrix[1,H_H] Spmatrixi_HTL_2H
matrix[H,H] matrix_HOIs_ihh

Spmatrixi_HTL_2H <- NULL
matrix_HOIs_ihh <- matrix(ncol=H,nrow=H)
Spmatrixi_HTL_2H = SpMatrix_H_herbivore[i,]; 
n = 1
m=1

for (h_h in 1:H_H) {
  if(h_h == vector_H_herbivore[n] + 1){
    n = n+1
    m = 1 }
  matrix_HOIs_ihh[n,m] = Spmatrixi_HTL_2H[h_h];
  m = m+1 
}

if(sum(rowSums(is.na(matrix_HOIs_ihh))==FV)){
  print("One or multiple row(s) in matrix_HOIs_ihh contain only NAs")
}

#---- Interaction (HOIs) matrix of floral visitors with floral visitors ----

matrix[1,FV_FV] Spmatrixi_HTL_2FV
matrix[FV,FV] matrix_HOIs_iff

Spmatrixi_HTL_2FV <- NULL
matrix_HOIs_iff <- matrix(ncol=FV,nrow=FV)
Spmatrixi_HTL_2FV = SpMatrix_FV_FV[i,]; 
n = 1
m=1

for (f_f in 1:FV_FV) {
  if(f_f == vector_FV_FV[n] + 1){
    n = n+1
    m = 1 }
  matrix_HOIs_iff[n,m] = Spmatrixi_HTL_2FV[f_f];
  m = m+1 
}
if(sum(rowSums(is.na(matrix_HOIs_iff))==FV)){
  print("One or multiple row(s) in matrix_HOIs_iff contain only NAs")
  }
#---- Interaction (HOIs) matrix of floral visitors with herbivore ----

matrix[1,FV_H] Spmatrixi_HTL_FvH
matrix[S,FV] matrix_HOIs_ifh

Spmatrixi_HTL_FvH <- NULL
matrix_HOIs_ifh <- matrix(ncol=FV,nrow=H)
Spmatrixi_HTL_FvH <- SpMatrix_FV_herbivore[i,]; 
n = 1
m=1
for (Fvh in 1:FV_H) {
  if (n > nrow(matrix_HOIs_ifh)) next
  if(Fvh == vector_FV_herbivore[n] + 1){
    n = n+1
    m = 1 }
  matrix_HOIs_ifh[n,m] = Spmatrixi_HTL_FvH[Fvh];
  m = m+1 
}


if(sum(rowSums(is.na(matrix_HOIs_ifh))==T)){
  print("One or multiple row(s) in matrix_HOIs_ifh contain only NAs")
}

#---- Model ----        
    for(s in 1:S){ // for one competing species j in alpha_ij, here s = species j
        alpha_ij[i,s] = exp((1-Intra[s]) * alpha_generic[1] + Intra[s] * alpha_intra[1] + (1-Intra[s]) * alpha_hat_ij[s]);
        
        for(k in 1:S){ // for all third competing species k in HOIs_ijk, here k = plant species  
        beta_ijk[s,k] = exp(beta_generic[1] + beta_hat_ijk[s,k]) ;
        }
        
        matrix_beta_ijk[i,s] = sum(beta_ijk[s,] .* matrix_HOIs[s,]);
        
        // for one herbivore species h in beta_H_ijh, here h = herbivor species h and j = floral visitor species f
    for(h in 1:H){
        beta_H_ijh[s,h] = exp(beta_H_generic[1] + beta_H_hat_ih[s,h]);
}
         matrix_beta_H_ijh[i,s] = sum(beta_H_ijh[s,] .* matrix_HOIs_ifh[s,]);

        // for one floral visitor species h in beta_H_ijf, here h = herbivor species h and j = floral visitor species f
    for(fv in 1:FV){
        beta_F_ijf[s,fv] = exp(beta_FV_generic[1] + beta_FV_hat_ijf[s,fv]);
}
matrix_beta_F_ijf[i,s] = sum(beta_F_ijf[s,] .* matrix_HOIs_ijf[s,]);
    }
    
    for(h in 1:H){ // for one herbivore species h in gamma_H_ih, here h = species h
        gamma_H_ih[i,h] = exp(gamma_H[1] + gamma_H_hat_ih[h]);
        
        // for two herbivores inbeta_ihh, here h = herbivore species h
    for(h2 in 1:H){
        beta_2H_ihh[h,h2] = exp(beta_2H_generic[1] + beta_2H_hat_ihh[h,h2]);
}
matrix_beta_2H_ihh[i,h] = sum(beta_2H_ihh[h,] .* matrix_HOIs_ihh[h,]);

}
    for(fv in 1:FV){ // for  floral visitor species f in gamma_FV_if, here f = species f
        gamma_FV_if[i,fv] = exp(gamma_FV[1] + gamma_FV_hat_if[fv]);
        
        // for two floral visitor species fv inbeta_iff, here f = floral visitor species f
    for(fv2 in 1:FV){
        beta_2FV_iff[fv,fv2] = exp(beta_2FV_generic[1] + beta_2FV_hat_iff[fv,fv2]);
}
matrix_beta_2FV_iff[i,fv] = sum(beta_2FV_iff[fv,] .* matrix_HOIs_iff[fv,]);

        // for one floral visitor species fv and one herbivore inbeta_ifh, here f = floral visitor species f and h= herbivore species h
    for(h in 1:H){
        beta_FvH_ifh[fv,h] = exp(beta_FvH_generic[1] + beta_FvH_hat_ifh[fv,h]);
}
matrix_beta_FvH_ifh[i,fv] = sum(beta_FvH_ifh[fv,] .* matrix_HOIs_ifh[fv,]);

}

     HOI_effects[i] = sum(matrix_beta_ijk[i,]);
     
     D_effects[i] = sum(alpha_ij[i,] .* SpMatrix[i,]) +  sum( gamma_H_ih[i,] .* SpMatrix_herbivore[i,]) + sum( gamma_FV_if[i,] .* SpMatrix_floralvis[i,]);
     
    
    F_hat[i] = lambda_i[i]/ D_effects[i]) + HOI_effects[i];
  }
  Fecundity ~ poisson(F_hat); // change to negative binomial
}
