// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

data{
  int<lower = 1> N; // Number of plots
  int<lower = 1> S; // Number of plant species
  int<lower = 1> H; // Number of herbivores species
  int<lower = 1> FV; // Number of floral visitors species

  int Fecundity[N];  // Fecundity of the focal species in each plot
  matrix[N,S] SpMatrix;  // Matrix of abundances for each species (including abundances of non-focal individuals of the focal species)
  matrix[N,H] SpMatrix_H; // Matrix of abundances for each herbivores species 
  matrix[N,FV] SpMatrix_FV; // Matrix of abundances for each floral visitors species
  int<lower = 0> Intra[S];  // Indicator boolean variable to identify the focal species (0 for non-focal and 1 for focal). Included for easier calculations
  
  matrix[S,S] matrix_HOIs_plant[N]; // Matrix of abundances for each plant species with each other plant species

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
  real<lower=0> lambdas[1];

  real<lower = -5, upper = 5> alpha_intra_tilde[1]; // direct interaction intra plants - plants;
    
  real<lower = -5, upper = 5> alpha_generic_tilde[1];
  vector<lower = -5, upper = 5>[S] alpha_hat_ij_tilde; // direct interaction inter plants - plants ; species -specific term 
  vector<lower = 0>[S] local_shrinkage_ij; // direct interaction inter plants - plants, shrinkage
   
  real<lower = -5, upper = 5> gamma_H_generic_tilde[1]; // direct interaction plants - herbivores ; generic
  vector<lower = -5, upper = 5>[H] gamma_H_hat_ih_tilde; // direct interaction plants - herbivores ; species -specific term
  vector<lower = 0>[H] gamma_H_local_shrinkage_ih; // direct interaction plants - herbivores ; shrinkage
  
  real<lower = -5, upper = 5> gamma_FV_generic_tilde[1]; // direct interaction plants - FV ; generic
  vector<lower = -5, upper = 5>[FV] gamma_FV_hat_if_tilde; // direct interaction plants - FV ; species -specific term
  vector<lower = 0>[FV] gamma_FV_local_shrinkage_if; // direct interaction plants - FV ; shrinkage
  
  real<lower = 0> c2_tilde;
  real<lower = 0> tau_tilde;
  
  real<lower=0> disp_dev; // species-specific dispersion deviation parameter,
  // defined for the negative binomial distribution used to reflect seed production (perform)
  // disp_dev = 1/sqrt(phi)
  
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

  // Declare objects necessary for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //a matrix of the species specific alpha values for each species and plot (interaction_effects), and a matrix
  //of the the alpha*N values for each species.
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

  alpha_generic[1] = 3*alpha_generic_tilde[1]-6; 
    alpha_intra[1] = 3*alpha_intra_tilde[1]-6;
  gamma_H[1] = 3*gamma_H_generic_tilde[1]-6;
    gamma_FV[1] = 3*gamma_FV_generic_tilde[1]-6; 

  tau = tau0*tau_tilde; 	// tau ~ cauchy(0, tau0)
  c2 = slab_scale2*c2_tilde;	// c2 ~ inv_gamma(half_slab_df, half_slab_df*slab_scale2)

  // Compute prior of sparsity approach - this calculation follows equation 2.8 in Piironen and Vehtari 2013

    for(s in 1:S){ // direct interaction INTER plants
      local_shrinkage_ij_tilde[s] = sqrt( c2 * square(local_shrinkage_ij[s]) / (c2 + square(tau) * square(local_shrinkage_ij[s])) );
      alpha_hat_ij[s] = tau * local_shrinkage_ij_tilde[s] * alpha_hat_ij_tilde[s];

       }
       
      for(h in 1:H){ // direct interaction plant - H
    gamma_H_local_shrinkage_ih_tilde[h] = sqrt( c2 * square(gamma_H_local_shrinkage_ih[h]) / (c2 + square(tau) * square(gamma_H_local_shrinkage_ih[h])) );
    gamma_H_hat_ih[h] = tau * gamma_H_local_shrinkage_ih_tilde[h] * gamma_H_hat_ih_tilde[h];

      }
    
    for(fv in 1:FV){ // direct interaction plant - FV
      gamma_FV_local_shrinkage_if_tilde[fv] = sqrt( c2 * square(gamma_FV_local_shrinkage_if[fv]) / (c2 + square(tau) * square(gamma_FV_local_shrinkage_if[fv])));
      gamma_FV_hat_if[fv] = tau * gamma_FV_local_shrinkage_if_tilde[fv] * gamma_FV_hat_if_tilde[fv];

      }
 // implement the biological model
  for(i in 1:N){
    lambda_ei[i] = lambdas[1];
    for(s in 1:S){
      alpha_eij[i,s] = (1-Intra[s]) * alpha_generic[1] + Intra[s] * alpha_intra[1] + alpha_hat_ij[s];
        

        }
  
    for(h in 1:H){ // for one herbivore species h in gamma_H_ih, here h = species h
        gamma_H_eih[i,h] = gamma_H[1] + gamma_H_hat_ih[h];
        

      }
      
    for(fv in 1:FV){ // for  floral visitor species f in gamma_FV_if, here f = species f
        gamma_FV_eif[i,fv] = gamma_FV[1] + gamma_FV_hat_if[fv];
        
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
  //alpha_hat_ij ~ normal(0,1);
  //gamma_hat_if ~ normal(0,1);
  disp_dev ~ cauchy(0, 1);  // safer to place prior on disp_dev than on phi
  
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

  
 for(i in 1:N){
  Fecundity[i] ~ neg_binomial_2(F_hat[i],(disp_dev^2)^(-1)); 
   }

}
