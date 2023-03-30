// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

data{
  int<lower = 1> N; // Number of plots
  int<lower = 1> S; // Number of plant species
  int<lower = 1> H; // Number of herbivores species
  int<lower = 1> FV; // Number of floral visitors species
  int<lower = 0> run_estimation;

  int Fecundity[N];  // Fecundity of the focal species in each plot
  matrix[N,S] SpMatrix;  // Matrix of abundances for each species (including abundances of non-focal individuals of the focal species)
  matrix[N,H] SpMatrix_H; // Matrix of abundances for each herbivores species 
  matrix[N,FV] SpMatrix_FV; // Matrix of abundances for each floral visitors species
  int<lower = 0> Intra[S];  // Indicator boolean variable to identify the focal species (0 for non-focal and 1 for focal). Included for easier calculations

   // vector telling which interactions to include
   
  int Inclusion_ij[1,S];  // Boolean indicator variables to identify the plant species
  int Inclusion_FV[1,FV]; // Boolean indicator variables to identify the floral visitor species
  int Inclusion_H[1,H]; // Boolean indicator variables to identify the herbivore species
  
}

parameters{
  real<lower=0> lambdas[1];
  real alpha_generic_tilde[1];
  real alpha_intra_tilde[1];
  vector[S] alpha_hat_ij;

  real gamma_H_generic_tilde[1]; // direct interaction plants - herbivores ; generic
  vector[H] gamma_H_hat_ih; // direct interaction plants - herbivores ; species -specific term

  real gamma_FV_generic_tilde[1]; // direct interaction plants - FV ; generic
  vector[FV] gamma_FV_hat_if; // direct interaction plants - FV ; species -specific term
vector[S] beta_plant_generic_tilde; // HOIs plants
  matrix[S,S] beta_plant_hat_ijk_tilde; // HOIs plants
  matrix<lower = 0>[S,S] beta_plant_local_shrinkage_ijk; // HOIs plants
  
  vector[S] beta_H_generic_tilde; // HOIs herbivores
  matrix[S,H] beta_H_hat_ijh_tilde; // HOIs herbivores 
  matrix<lower = 0>[S,H] beta_H_local_shrinkage_ijh; // HOIs herbivores
  
  vector[S] beta_FV_generic_tilde; // HOIs floral visitors
  matrix[S,FV] beta_FV_hat_ijf_tilde; // HOIs floral visitors
  matrix<lower = 0>[S,FV] beta_FV_local_shrinkage_ijf; // HOIs floral visitors
  
  vector[FV] beta_2FV_generic_tilde; // HOIs 2 floral visitors
  matrix[FV,FV] beta_2FV_hat_iff_tilde; // HOIs 2 floral visitors
  matrix<lower = 0>[FV,FV] beta_2FV_local_shrinkage_iff; // HOIs 2 floral visitors
  
  vector[H] beta_2H_generic_tilde; // HOIs 2 herbivores
  matrix[H,H] beta_2H_hat_ihh_tilde; // HOIs 2 herbivores 
  matrix<lower = 0>[H,H] beta_2H_local_shrinkage_ihh; // HOIs 2 herbivores
  
  vector[FV] beta_FvH_generic_tilde; // HOIs 1 floral visitor & 1 herbivore
  matrix[FV,H] beta_FvH_hat_ifh_tilde; // HOIs 1 floral visitor & 1 herbivore
  matrix<lower = 0>[FV,H] beta_FvH_local_shrinkage_ifh; // HOIs 1 floral visitor & 1 herbivore


  real<lower=0> disp_dev; // species-specific dispersion deviation parameter,
  // defined for the negative binomial distribution used to reflect seed production (perform)
  // disp_dev = 1/sqrt(phi)
  
}

transformed parameters{
 // Declare objects necessary for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //     a matrix of the species specific alpha values for each species and plot (interaction_effects), and a matrix
  //     of the the alpha*N values for each species.
  vector[N] F_hat;
  vector[N] interaction_effects;
  //vector[N] pollinator_effects;
  
  
  // loop parameters
  vector[N] lambda_ei;
  matrix[N,S] alpha_eij;
  matrix[N,H] gamma_H_eih; //gamma herbivores of focal i
  matrix[N,FV] gamma_FV_eif; //gamma floral visitors of focal i
  
    
  // scale  values
  vector[1] alpha_intra;
  vector[1] alpha_generic;
  vector[1] gamma_H;  // direct interaction plants - H
  vector[1] gamma_FV;  // direct interaction plants - FV

  alpha_generic[1] = alpha_generic_tilde[1]; 
    alpha_intra[1] = alpha_intra_tilde[1]; 
  gamma_H[1] = gamma_H_generic_tilde[1]; 
    gamma_FV[1] = gamma_FV_generic_tilde[1]; 

  
  
 // implement the biological model
  for(i in 1:N){
    lambda_ei[i] = lambdas[1];
    for(s in 1:S){
      alpha_eij[i,s] = (1-Intra[s]) * alpha_generic[1] + Intra[s] * alpha_intra[1] + (1-Intra[s]) *alpha_hat_ij[s]*Inclusion_ij[1,s];
      
        }
  
    for(h in 1:H){ // for one herbivore species h in gamma_H_ih, here h = species h
        gamma_H_eih[i,h] = gamma_H[1] + Inclusion_H[1,h]*gamma_H_hat_ih[h];
        
      }
      
    for(fv in 1:FV){ // for  floral visitor species f in gamma_FV_if, here f = species f
        gamma_FV_eif[i,fv] = gamma_FV[1] + Inclusion_FV[1,fv]*gamma_FV_hat_if[fv];
        
    }


    interaction_effects[i] = sum(alpha_eij[i,] .* SpMatrix[i,]) +  sum(gamma_H_eih[i,] .* SpMatrix_H[i,]) + sum( gamma_FV_eif[i,] .* SpMatrix_FV[i,]);
     
 
          F_hat[i] = exp(lambda_ei[i] + interaction_effects[i]);

  }

}

model{
  // set regular priors
  alpha_generic_tilde ~ normal(0,1);
  alpha_intra_tilde ~ normal(0,1);
  //gamma_generic_tilde ~ normal(0,1);
    
  lambdas ~ normal(0, 1);
  alpha_hat_ij ~ normal(0,1);
  gamma_FV_hat_if ~ normal(0,1);
  gamma_H_hat_ih ~ normal(0,1);
  disp_dev ~ cauchy(0, 1);  // safer to place prior on disp_dev than on phi
  

 for(i in 1:N){
  Fecundity[i] ~ neg_binomial_2(F_hat[i],(disp_dev^2)^(-1)); 
   }

}
generated quantities{
  vector[N] F_sim;
    if(run_estimation==1){
 for(i in 1:N){
    if(Fecundity[i] <= 0) break ;
    F_sim[i] = neg_binomial_2_rng(F_hat[i],(disp_dev^2)^(-1));
              }
    }
}
