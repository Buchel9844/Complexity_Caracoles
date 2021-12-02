// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

data{
  int<lower = 1> N; // Number of plots
  int<lower = 1> S; // Number of plant species
  int<lower = 1> H; // Number of herbivores species
  int<lower = 1> FV; // Number of floral visitors species
  int<lower = 1> H_comp; // Number of HOIs with herbivores species
  int<lower = 1> FV_comp; // Number of HOIs with floral visitors species
  int Fecundity[N]; // Fecundity of the focal species in each plot
  //int plot[N];   // Indicator variable for plot
  matrix[N,S] SpMatrix; // Matrix of abundances for each plant species 
  matrix[N,H] SpMatrix_herbivore; // Matrix of abundances for each herbivores species 
  matrix[N,FV] SpMatrix_floralvis; // Matrix of abundances for each floral visitors species

  matrix[N,H_comp] SpMatrix_herbivore_comp; // Matrix of abundances for each herbivores species and competitor
  matrix[N,FV_comp] SpMatrix_floralvis_comp; // Matrix of abundances for each herbivores species and competitor
  
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
  vector[2] lambdas;
  vector[2] alpha_generic_tilde; // direct interaction inter plants - plants ; generic
  vector[2] alpha_intra_tilde; // direct interaction intra plants - plants;
  vector[S] alpha_hat_ij_tilde; // direct interaction inter plants - plants ; species -specific term 
  vector<lower = 0>[S] local_shrinkage_ij; // direct interaction inter plants - plants, shrinkage
   
  vector[2] gamma_H_generic_tilde; // direct interaction plants - herbivores ; generic
  vector[H] gamma_H_hat_ih_tilde; // direct interaction plants - herbivores ; species -specific term
  vector<lower = 0>[H] gamma_H_local_shrinkage_ih; // direct interaction plants - herbivores ; shrinkage
  
  vector[2] gamma_FV_generic_tilde; // direct interaction plants - FV ; generic
  vector[FV] gamma_FV_hat_if_tilde; // direct interaction plants - FV ; species -specific term
  vector<lower = 0>[FV] gamma_FV_local_shrinkage_if; // direct interaction plants - FV ; shrinkage
  
  vector[2] beta_plant_generic_tilde; // HOIs plants
  matrix[S,S] beta_plant_hat_ijk_tilde; // HOIs plants
  matrix<lower = 0>[S,S] beta_plant_local_shrinkage_ijk; // HOIs plants
  
  vector[2] beta_H_generic_tilde; // HOIs herbivores
  matrix[S,H_comp] beta_H_hat_ijh_tilde; // HOIs herbivores 
  matrix<lower = 0>[S,H_comp] beta_H_local_shrinkage_ijh; // HOIs herbivores
  
  vector[2] beta_FV_generic_tilde; // HOIs floral visitors
  matrix[S,FV_comp] beta_FV_hat_ijf_tilde; // HOIs floral visitors
  matrix<lower = 0>[S,FV_comp] beta_FV_local_shrinkage_ijf; // HOIs floral visitors
  
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
  
   matrix[H_comp,H_comp] beta_H_hat_ijh; // HOIs herbivores
  matrix[H_comp,H_comp] beta_H_local_shrinkage_ijh_tilde;// HOIs herbivores
   matrix[FV_comp,FV_comp] beta_FV_hat_ijf; // HOIs floral visitors
  matrix[FV_comp,FV_comp] beta_FV_local_shrinkage_ijf_tilde; // HOIs floral visitors
  
  vector[2] alpha_generic;  // direct interaction INTER plants
  vector[2] alpha_intra;  // direct interaction INTRA plants
  vector[2] gamma_H;  // direct interaction plants - H
  vector[2] gamma_FV;  // direct interaction plants - FV
  vector[2] beta_plant_generic; //HOIs plant 
  vector[2] beta_H_generic; //HOIs herbivore
  vector[2] beta_FV_generic; //HOIs FV


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
       for(fv in 1:FV_comp){ //HOIs plant-plant-FV 
      beta_FV_local_shrinkage_ijf_tilde[s,fv] = sqrt( c2 * square(beta_FV_local_shrinkage_ijf[s,fv]) / (c2 + square(tau) * square(beta_FV_local_shrinkage_ijf[s,fv])) );
      beta_FV_hat_ijf[s,fv] = tau * beta_FV_local_shrinkage_ijf_tilde[s,fv] * beta_FV_hat_ijf_tilde[s,fv];

       }
       for(h in 1:H_comp){ //HOIs plant-plant-H
      beta_H_local_shrinkage_ijh_tilde[s,h] = sqrt( c2 * square(beta_H_local_shrinkage_ijh[s,h]) / (c2 + square(tau) * square(beta_H_local_shrinkage_ijh[s,h])) );
      beta_H_hat_ijh[s,h] = tau * beta_H_local_shrinkage_ijh_tilde[s,h] * beta_H_hat_ijh_tilde[s,h];
       }
    }
    
    for(h in 1:H){ # direct interaction plant - H
      gamma_H_local_shrinkage_ih_tilde[h] = sqrt( c2 * square(gamma_H_local_shrinkage_ih[h]) / (c2 + square(tau) * square(gamma_H_local_shrinkage_ih[h])) );
      gamma_H_hat_ih[h] = tau * gamma_H_local_shrinkage_ih_tilde[h] * gamma_H_hat_ih_tilde[h];

}

    for(fv in 1:FV){ # direct interaction plant - FV
      gamma_FV_local_shrinkage_if_tilde[fv] = sqrt( c2 * square(gamma_FV_local_shrinkage_if[fv]) / (c2 + square(tau) * square(gamma_FV_local_shrinkage_if[fv])));
      gamma_FV_hat_if[fv] = tau * gamma_FV_local_shrinkage_if_tilde[fv] * gamma_FV_hat_if_tilde[fv];

}


  // scale the lambdas,alphas,gammas and HOIs values
  alpha_generic[1] = 3 * alpha_generic_tilde[1] - 6;
    alpha_generic[2] = 0.5 * alpha_generic_tilde[2];
    
  alpha_intra[1] = 3 * alpha_intra_tilde[1] - 6;
  alpha_intra[2] = 0.5 * alpha_intra_tilde[2];
  
  
  gamma_H[1] = 3 * gamma_H_generic_tilde[1] - 6;
    gamma_H[2] = 0.5 * gamma_H_generic_tilde[2];
    
    gamma_FV[1] = 3 * gamma_FV_generic_tilde[1] - 6;
  gamma_FV[2] = 0.5 * gamma_FV_generic_tilde[2];
  
    beta_plant_generic[1] = 3 * beta_generic_tilde[1] - 6;
  beta_plant_generic[2] = 0.5 * beta_generic_tilde[2];
  
    beta_H_generic[1] = 3 * beta_generic_tilde[1] - 6;
  beta_H_generic[2] = 0.5 * beta_generic_tilde[2];
  
    beta_FV_generic[1] = 3 * beta_generic_tilde[1] - 6;
  beta_FV_generic[2] = 0.5 * beta_generic_tilde[2];
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
  
  
  matrix[S,S] beta_ijk;
  matrix[N,S] matrix_beta_ijk;
  matrix[1,S] Spmatrixi;
  matrix[1,S] Spmatrixi_HTL;


  // set regular priors
  alpha_generic_tilde ~ normal(0,1);
  alpha_intra_tilde ~ normal(0,1);
  
  gamma_H_generic_tilde ~ normal(0,1);
  gamma_FV_generic_tilde ~ normal(0,1);
  
  beta_generic_tilde ~ normal(0,1);
  lambdas ~ normal(0, 1);


  // set the hierarchical priors for the Finnish horseshoe (regularized horseshoe) (Piironen and Vehtari 2017)
  // Following the stan implementation from https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

    alpha_hat_ij_tilde ~ normal(0,1);
    local_shrinkage_ij ~ cauchy(0,1);
    
    gamma_H_hat_ih_tilde ~ normal(0,1);
    gamma_H_local_shrinkage_ih ~ cauchy(0,1);
    
    gamma_FV_hat_if_tilde ~ normal(0,1);
    gamma_FV_local_shrinkage_if ~ cauchy(0,1);

  
for (s in 1:S){
    beta_hat_ijk_tilde[s,] ~ normal(0,1); 
    beta_local_shrinkage_ijk[s,] ~ cauchy(0,1);

}
    

  tau_tilde ~ cauchy(0,1);
  c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);

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
        
    // creation of a matrix of  by herbivore of the interaction HOIs_ijh for herbivores on plants
      Spmatrixi_HTL_H[1,] = SpMatrix_herbivore[i,]; 
    for (h in 1:H_comp) {
        matrix_HOIs[n,m] = Spmatrixi_HTL_H[1,H_comp];
        }
        
    for(s in 1:S){ // for one competing species j in alpha_ij, here s = species j
        alpha_ij[i,s] = exp((1-Intra[s]) * alpha_generic[1] + Intra[s] * alpha_intra[1] + (1-Intra[s]) * alpha_hat_ij[s] + (1-Intra[s]) * alpha_generic[2]);
        
        for(k in 1:S){ // for all third competing species k in HOIs_ijk, here k = species k 
        beta_ijk[s,k] = exp(beta_generic[1] + beta_hat_ijk[s,k]) ;
        }
        
        matrix_beta_ijk[i,s] = sum(beta_ijk[s,] .* matrix_HOIs[s,]);
    }
    
    for(h in 1:H){ // for one herbivore species h in gamma_H_ih, here h = species h
        gamma_H_ih[i,h] = exp(gamma_H[1] + gamma_H_hat_ih[h] + gamma_H[2]);

}
    for(fv in 1:FV){ // for one floral visitor species f in gamma_FV_if, here f = species f
        gamma_FV_if[i,fv] = exp(gamma_FV[1] + gamma_FV_hat_if[fv] + gamma_FV[2]);

}
    for(h1 in 1:H_comp){ // for one herbivore species h in beta_H_ijh, here h = herbivor species h and j = floral visitor species f
    for(h2 in 1:H_comp){
        beta_H_ih[h1,h2] = exp(beta_H[1] + beta_H_hat_ih[h1,h2] + beta_H[2]);
}
}

    for(fv1 in 1:FV_comp){ // for one floral visitor species h in beta_H_ijh, here h = herbivor species h and j = floral visitor species f
    for(fv2 in 1:FV_comp){
        beta_FV_ih[fv1,fv2] = exp(beta_FV[1] + beta_FV_hat_ih[fv1,fv2] + beta_FV[2]);
}
}

     HOI_effects[i] = sum(matrix_beta_ijk[i,]);
     
     D_effects[i] = sum(alpha_ij[i,] .* SpMatrix[i,]) +  sum( gamma_H_ih[i,] .* SpMatrix_herbivore[i,]) + sum( gamma_FV_if[i,] .* SpMatrix_floralvis[i,]);
     
    
    F_hat[i] = lambda_i[i]/ D_effects[i]) + HOI_effects[i];
  }
  Fecundity ~ poisson(F_hat); // change to negative binomial
}
