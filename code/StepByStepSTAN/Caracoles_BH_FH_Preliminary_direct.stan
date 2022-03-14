// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

data{
  int<lower = 1> N; // Number of plots
  int<lower = 1> S; // Number of plant species
  int<lower = 1> H; // Number of herbivores species
  int<lower = 1> FV; // Number of floral visitors species

  int Fecundity[N]; // Fecundity of the focal species in each plot
  //int plot[N];   // Indicator variable for plot
  matrix[N,S] SpMatrix; // Matrix of abundances for each plant species 
  matrix[N,H] SpMatrix_H; // Matrix of abundances for each herbivores species 
  matrix[N,FV] SpMatrix_FV; // Matrix of abundances for each floral visitors species
  
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
 
  
  vector[1] alpha_generic;  // direct interaction INTER plants
  vector[1] alpha_intra;  // direct interaction INTRA plants
  vector[1] gamma_H;  // direct interaction plants - H
  vector[1] gamma_FV;  // direct interaction plants - FV


  tau = tau0*tau_tilde; 	// tau ~ cauchy(0, tau0)
  c2 = slab_scale2*c2_tilde;	// c2 ~ inv_gamma(half_slab_df, half_slab_df*slab_scale2)

  // This calculation follows equation 2.8 in Piironen and Vehtari 2013

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


  // scale the lambdas,alphas,gammas and HOIs values
  alpha_generic[1] = 3 * alpha_generic_tilde[1] - 6;

  alpha_intra[1] = 3 * alpha_intra_tilde[1] - 6;

  gamma_H[1] = 3 * gamma_H_generic_tilde[1] - 6;

  gamma_FV[1] = 3 * gamma_FV_generic_tilde[1] - 6;


}

model{
  // Declare objects necessary for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //     a matrix of the species specific alpha values for each species and plot (interaction_effects), and a matrix
  //     of the the alpha*N values for each species.
  vector[N] F_hat;
  vector[N] D_effects;
  vector[N] HOI_effects;
  
  vector[N] lambda_i; //lambda of focal i for obs s
  
  matrix[N,S] alpha_ij; //alpha of focal i
  matrix[N,H] gamma_H_ih; //gamma herbivores of focal i
  matrix[N,FV] gamma_FV_if; //gamma floral visitors of focal i


  // set regular priors
  alpha_generic_tilde ~ normal(0,1);
  alpha_intra_tilde ~ normal(0,1);
  
  gamma_H_generic_tilde ~ normal(0,1);
  gamma_FV_generic_tilde ~ normal(0,1);

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



// //---- Implement the biological model ----
  for(i in 1:N){ // for the observation N 
     lambda_i[i] = exp(lambdas[1]);
    for(s in 1:S){ // for one competing species j in alpha_ij, here s = species j
        alpha_ij[i,s] = exp((1-Intra[s]) * alpha_generic[1] + Intra[s] * alpha_intra[1] + (1-Intra[s]) * alpha_hat_ij[s]);
  
}
    for(h in 1:H){ // for one herbivore species h in gamma_H_ih, here h = species h
        gamma_H_ih[i,h] = exp(gamma_H[1] + gamma_H_hat_ih[h]);
      
}
    for(fv in 1:FV){ // for  floral visitor species f in gamma_FV_if, here f = species f
        gamma_FV_if[i,fv] = exp(gamma_FV[1] + gamma_FV_hat_if[fv]);
        
}

     
     D_effects[i] = sum(alpha_ij[i,] .* SpMatrix[i,]) +  sum( gamma_H_ih[i,] .* SpMatrix_H[i,]) + sum( gamma_FV_if[i,] .* SpMatrix_FV[i,]);
     
    F_hat[i] = lambda_i[i]/ (D_effects[i]);
  }
  Fecundity ~ poisson(F_hat); // change to negative binomial
}
