// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

data{
  int<lower = 1> N; // Number of plots
  int<lower = 1> S; // Number of plant species
  int<lower = 1> H; // Number of herbivores species
  int<lower = 1> FV; // Number of floral visitors species
  
  int Fecundity[N];  // Fecundity of the focal species in each plot
  
  matrix[N,S] SpMatrix; // Matrix of abundances for each plant species 
  matrix[N,H] SpMatrix_H; // Matrix of abundances for each herbivores species 
  matrix[N,FV] SpMatrix_FV; // Matrix of abundances for each floral visitors species
  
  real matrix_HOIs_plant[N,S,S]; // Matrix of abundances for each plant species with each other plant species
  
  real matrix_HOIs_ijh[N,S,H]; // Matrix of abundances for each herbivores species and competitor
  real matrix_HOIs_ijf[N,S,FV]; // Matrix of abundances for each herbivores species and competitor
  
  real matrix_HOIs_iff[N,FV,FV]; // Matrix of abundances for each floral visitors species with each floral visitors species
  real matrix_HOIs_ihh[N,H,H]; // Matrix of abundances for each herbivores species with each herbivores species
  real matrix_HOIs_ifh[N,FV,H]; // Matrix of abundances for each floral visitors species with each herbivores species
  
  //vector[N] env;  // Environmental values for each plot
  int<lower = 0> Intra[S];  // Indicator boolean variable to identify the focal species (0 for non-focal and 1 for focal). Included for easier calculations
 
  int Inclusion_ij[S];  // Boolean indicator variables to identify the plant species
  int Inclusion_FV[FV] // Boolean indicator variables to identify the floral visitor species
  int Inclusion_H[H] // Boolean indicator variables to identify the herbivore species
  
  int beta_Inclusion_plant[S,S];  // Boolean indicator variables to identify the HOIs plant-plant
  int beta_Inclusion_FV[S,FV];  // Boolean indicator variables to identify the HOIs 2 plants-FV
  int beta_Inclusion_H[S,H];  // Boolean indicator variables to identify the HOIs 2 plants- H

  int beta_Inclusion_2FV[FV,FV];  // Boolean indicator variables to identify the HOIs 1 plant - 2 FV
  int beta_Inclusion_2H[H,H];  // Boolean indicator variables to identify the HOIs 1 plant - 2 H
  int beta_Inclusion_FvH[FV,H];  // Boolean indicator variables to identify the HOIs 1 plant - 2 H

}

parameters{
  vector[2] lambdas;
  
  vector[2] alpha_intra_tilde;
  vector[2] alpha_generic_tilde;
  vector[S] alpha_hat_ij;
  
  vector[2] gamma_H_generic_tilde; // direct interaction plants - herbivores ; generic
  vector[H] gamma_H_hat_ih; // direct interaction plants - herbivores ; species -specific term

  vector[2] gamma_FV_generic_tilde; // direct interaction plants - FV ; generic
  vector[FV] gamma_FV_hat_if; // direct interaction plants - FV ; species -specific term

  vector[2] beta_plant_generic_tilde; // HOIs plants
  matrix[S,S] beta_plant_hat_ijk; // HOIs plants
  
  vector[2] beta_H_generic_tilde; // HOIs herbivores
  matrix[S,H] beta_H_hat_ijh; // HOIs herbivores
  
  vector[2] beta_FV_generic_tilde; // HOIs floral visitors
  matrix[S,FV] beta_FV_hat_ijf; // HOIs floral visitors
  
  vector[2] beta_2FV_generic_tilde; // HOIs 2 floral visitors
  matrix[FV,FV] beta_2FV_hat_iff; // HOIs 2 floral visitors
  
  vector[2] beta_2H_generic_tilde; // HOIs 2 herbivores
  matrix[H,H] beta_2H_hat_ihh; // HOIs 2 herbivores 
  
  vector[2] beta_FvH_generic_tilde; // HOIs 1 floral visitor & 1 herbivore
  matrix[FV,H] beta_FvH_hat_ifh; // HOIs 1 floral visitor & 1 herbivore
}

transformed parameters{
  vector[2] alpha_generic;
  vector[2] alpha_intra;
  
  vector[2] gamma_FV_generic;
  vector[2] gamma_H_generic;
  
  vector[2] beta_plant_generic
  vector[2] beta_H_generic
  vector[2] beta_FV_generic
  vector[2] beta_2FV_generic
  vector[2] beta_2H_generic
  vector[2] beta_FvH_generic
  
  // scale the alpha values
  alpha_intra[1] = 3 * alpha_intra_tilde[1] - 6;
    alpha_intra[2] = 0.5 * alpha_intra_tilde[2];
    
  alpha_generic[1] = 3 * alpha_generic_tilde[1] - 6;
    alpha_generic[2] = 0.5 * alpha_generic_tilde[2];
    
  gamma_H_generic[1] = 3 * gamma_H_generic_tilde[1] - 6;
  
  gamma_FV_generic[1] = 3 * gamma_FV_generic_tilde[1] - 6;
  
  beta_plant_generic[1] = 3 * beta_plant_generic_tilde[1] - 6;
    beta_plant_generic[2] = 0.5 * beta_plant_generic_tilde[2];
  
  beta_H_generic[1] = 3 * beta_H_generic_tilde[1] - 6;
    beta_plant_generic[2] = 0.5 * beta_plant_generic_tilde[2];
    
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
  
  vector[N] lambda_i; //lambda of focal i for obs s
  
  matrix[S,S] matrix_HOIs;
  matrix[N,S] alpha_ij; //alpha of focal i
  matrix[N,H] gamma_H_ih; //gamma herbivores of focal i
  matrix[N,FV] gamma_FV_if; //gamma floral visitors of focal i
  
  // the matrix of interaction of plant-plant-plant HOIs was not created prior to the script
  matrix[S,S] beta_ijk; // HOIs plant 
    matrix[N,S] matrix_beta_ijk;
  matrix[S,FV] beta_F_ijf;
    matrix[N,S] matrix_beta_F_ijf;
  matrix[S,H] beta_H_ijh; 
    matrix[N,S] matrix_beta_H_ijh;
  matrix[FV,FV] beta_2F_iff;
    matrix[N,FV] matrix_beta_2FV_iff;
  matrix[H,H] beta_2H_ihh; 
    matrix[N,H] matrix_beta_2H_ihh;
  matrix[FV,H] beta_FvH_ifh; 
    matrix[N,FV] matrix_beta_FvH_ifh;

  // set regular priors
  alpha_generic_tilde ~ normal(0,1);
  alpha_intra_tilde ~ normal(0,1);
  
  gamma_H_generic_tilde ~ normal(0,1);
  gamma_FV_generic_tilde ~ normal(0,1);
  
  beta_plant_generic_tilde ~ normal(0,1);
  beta_H_generic_tilde ~ normal(0,1);
  beta_FV_generic_tilde ~ normal(0,1);
  beta_2FV_generic_tilde ~ normal(0,1);
  beta_2H_generic_tilde ~ normal(0,1);
  beta_FvH_generic_tilde ~ normal(0,1);
  lambdas ~ normal(0, 1);


// //---- Implement the biological model ----
  for(i in 1:N){ // for the observation N 
     lambda_i[i] = exp(lambdas[1]);
    for(s in 1:S){ // for one competing species j in alpha_ij, here s = species j
        alpha_ij[i,s] = exp((1-Intra[s]) * alpha_generic[1] + Intra[s] * alpha_intra[1] + Inclusion_ij[s] *(1-Intra[s]) * alpha_hat_ij[s]);
        
        for(k in 1:S){ // for all third competing species k in HOIs_ijk, here k = plant species  
        beta_ijk[s,k] = exp(beta_plant_generic[1] + beta_Inclusion_plant[s,k]*beta_plant_hat_ijk[s,k]) ;
        beta_ijk[s,k] = beta_ijk[s,k]*matrix_HOIs_plant[i,s,k];
        }
        
        matrix_beta_ijk[i,s] = sum(beta_ijk[s,]);
        
        // for one herbivore species h in beta_H_ijh, here h = herbivor species h and j = plant species
    for(h in 1:H){
        beta_H_ijh[s,h] = exp(beta_H_generic[1] + beta_Inclusion_H[s,h]*beta_H_hat_ijh[s,h]);
        beta_H_ijh[s,h] = beta_H_ijh[s,h]*matrix_HOIs_ifh[i,s,h];
        
}
         matrix_beta_H_ijh[i,s] = sum(beta_H_ijh[s,]);

        // for one floral visitor species h in beta_H_ijf, here f = floral visitor species and j = plant species
    for(fv in 1:FV){
        beta_F_ijf[s,fv] = exp(beta_FV_generic[1] + beta_Inclusion_FV[s,fv]*beta_FV_hat_ijf[s,fv]);
         beta_F_ijf[s,fv] =  beta_F_ijf[s,fv]*matrix_HOIs_ijf[i,s,fv];
        
}
        matrix_beta_F_ijf[i,s] = sum(beta_F_ijf[s,]);
}
    for(h in 1:H){ // for one herbivore species h in gamma_H_ih, here h = species h
        gamma_H_ih[i,h] = exp(gamma_H[1] + Inclusion_FV[h]*gamma_H_hat_ih[h]);
        
        // for two herbivores inbeta_ihh, here h = herbivore species h
    for(h2 in 1:H){
        beta_2H_ihh[h,h2] = exp(beta_2H_generic[1] + beta_Inclusion_2H[h,h2]*beta_2H_hat_ihh[h,h2]);
        beta_2H_ihh[h,h2] = beta_2H_ihh[h,h2]*matrix_HOIs_ihh[i,h,h2];
        
}

    matrix_beta_2H_ihh[i,h] = sum(beta_2H_ihh[h,]);

}
    for(fv in 1:FV){ // for  floral visitor species f in gamma_FV_if, here f = species f
        gamma_FV_if[i,fv] = exp(gamma_FV[1] + Inclusion_FV[fv]*gamma_FV_hat_if[fv]);
        
        // for two floral visitor species fv inbeta_iff, here f = floral visitor species f
    for(fv2 in 1:FV){
        beta_2F_iff[fv,fv2] = exp(beta_2FV_generic[1] + beta_Inclusion_2FV[fv,fv2]*beta_2FV_hat_iff[fv,fv2]);
        beta_2F_iff[fv,fv2] = beta_2F_iff[fv,fv2]*matrix_HOIs_iff[i,fv,fv2];
}
        matrix_beta_2FV_iff[i,fv] = sum(beta_2F_iff[fv,]);

        // for one floral visitor species fv and one herbivore inbeta_ifh, here f = floral visitor species f and h= herbivore species h
    for(h in 1:H){
        beta_FvH_ifh[fv,h] = exp(beta_FvH_generic[1] + beta_Inclusion_FvH[fv,h]*beta_FvH_hat_ifh[fv,h]);
        beta_FvH_ifh[fv,h] =  beta_FvH_ifh[fv,h]*matrix_HOIs_ifh[i,fv,h];
}
        matrix_beta_FvH_ifh[i,fv] = sum(beta_FvH_ifh[fv,]);
}

     HOI_effects[i] = sum(matrix_beta_ijk[i,]) + sum(matrix_beta_F_ijf[i,]) +  sum(matrix_beta_H_ijh[i,])  + sum(matrix_beta_2H_ihh[i,]) +  sum(matrix_beta_FvH_ifh[i,]) + sum(matrix_beta_2FV_iff[i,]);
     
     D_effects[i] = sum(alpha_ij[i,] .* SpMatrix[i,]) +  sum( gamma_H_ih[i,] .* SpMatrix_H[i,]) + sum( gamma_FV_if[i,] .* SpMatrix_FV[i,]);
     
    F_hat[i] = lambda_i[i]/ (1+ D_effects[i] + HOI_effects[i]);
  }
  Fecundity ~ poisson(F_hat); // change to negative binomial
}
