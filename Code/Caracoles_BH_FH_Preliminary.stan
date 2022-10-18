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
  matrix[S,H] matrix_HOIs_ijh[N]; // Matrix of abundances for each herbivores species and competitor
  matrix[S,FV] matrix_HOIs_ijf[N]; // Matrix of abundances for each herbivores species and competitor
  
  matrix[FV,FV] matrix_HOIs_iff[N]; // Matrix of abundances for each floral visitors species with each floral visitors species
  matrix[H,H] matrix_HOIs_ihh[N]; // Matrix of abundances for each herbivores species with each herbivores species
  matrix[FV,H]  matrix_HOIs_ifh[N]; // Matrix of abundances for each floral visitors species with each herbivores species


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

  real<lower=-5,upper=5> alpha_intra_tilde[1]; // direct interaction intra plants - plants;
    
  real<lower=-5,upper=5> alpha_generic_tilde[1];
  vector[S] alpha_hat_ij_tilde; // direct interaction inter plants - plants ; species -specific term 
  vector<lower = 0>[S] local_shrinkage_ij; // direct interaction inter plants - plants, shrinkage
   
  real<lower=-5,upper=5> gamma_H_generic_tilde[1]; // direct interaction plants - herbivores ; generic
  vector[H] gamma_H_hat_ih_tilde; // direct interaction plants - herbivores ; species -specific term
  vector<lower = 0>[H] gamma_H_local_shrinkage_ih; // direct interaction plants - herbivores ; shrinkage
  
  real<lower=-5,upper=5> gamma_FV_generic_tilde[1]; // direct interaction plants - FV ; generic
  vector[FV] gamma_FV_hat_if_tilde; // direct interaction plants - FV ; species -specific term
  vector<lower = 0>[FV] gamma_FV_local_shrinkage_if; // direct interaction plants - FV ; shrinkage
  
  vector<lower = -5, upper = 5>[S] beta_plant_generic_tilde; // HOIs plants
  matrix<lower=-5,upper=5>[S,S] beta_plant_hat_ijk_tilde; // HOIs plants
  matrix<lower = 0>[S,S] beta_plant_local_shrinkage_ijk; // HOIs plants
  
  vector<lower = -5, upper = 5>[S] beta_H_generic_tilde; // HOIs herbivores
  matrix<lower=-5,upper=5>[S,H] beta_H_hat_ijh_tilde; // HOIs herbivores 
  matrix<lower = 0>[S,H] beta_H_local_shrinkage_ijh; // HOIs herbivores
  
  vector<lower = -5, upper = 5>[S] beta_FV_generic_tilde; // HOIs floral visitors
  matrix<lower=-5,upper=5>[S,FV] beta_FV_hat_ijf_tilde; // HOIs floral visitors
  matrix<lower = 0>[S,FV] beta_FV_local_shrinkage_ijf; // HOIs floral visitors
  
  vector<lower = -5, upper = 5>[FV] beta_2FV_generic_tilde; // HOIs 2 floral visitors
  matrix<lower=-5,upper=5>[FV,FV] beta_2FV_hat_iff_tilde; // HOIs 2 floral visitors
  matrix<lower = 0>[FV,FV] beta_2FV_local_shrinkage_iff; // HOIs 2 floral visitors
  
  vector<lower = -5, upper = 5>[H] beta_2H_generic_tilde; // HOIs 2 herbivores
  matrix<lower=-5,upper=5>[H,H] beta_2H_hat_ihh_tilde; // HOIs 2 herbivores 
  matrix<lower = 0>[H,H] beta_2H_local_shrinkage_ihh; // HOIs 2 herbivores
  
  vector<lower = -5, upper = 5>[FV] beta_FvH_generic_tilde; // HOIs 1 floral visitor & 1 herbivore
  matrix<lower=-5,upper=5>[FV,H] beta_FvH_hat_ifh_tilde; // HOIs 1 floral visitor & 1 herbivore
  matrix<lower = 0>[FV,H] beta_FvH_local_shrinkage_ifh; // HOIs 1 floral visitor & 1 herbivore
 
  
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
 
// Declare objects necessary for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //a matrix of the species specific alpha values for each species and plot (interaction_effects), and a matrix
  //of the the alpha*N values for each species.
  vector[N] F_hat;
  vector[N] interaction_effects;
  vector[N] HOI_effects;
  //vector[N] pollinator_effects;
  
  
  // loop parameters
  vector[N] lambda_ei;
  matrix[N,S] alpha_eij;
  matrix[N,H] gamma_H_eih; //gamma herbivores of focal i
  matrix[N,FV] gamma_FV_eif; //gamma floral visitors of focal i
  

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
  

  // scale  values
  vector[1] alpha_intra;
  vector[1] alpha_generic;
  vector[1] gamma_H;  // direct interaction plants - H
  vector[1] gamma_FV;  // direct interaction plants - FV
  vector[S] beta_plant_generic; //HOIs plant 
  vector[S] beta_H_generic; //HOIs herbivore
  vector[S] beta_FV_generic; //HOIs FV
  vector[FV] beta_2FV_generic; //HOIs 2 floral visitors
  vector[H] beta_2H_generic; //HOIs 2 herbivores
  vector[FV] beta_FvH_generic; //HOIs 1 floral visitor & 1 herbivore
  
  alpha_generic[1] = alpha_generic_tilde[1]; 
    alpha_intra[1] = alpha_intra_tilde[1]; 
  gamma_H[1] = gamma_H_generic_tilde[1]; 
    gamma_FV[1] = gamma_FV_generic_tilde[1]; 
  beta_plant_generic = beta_plant_generic_tilde; 
    beta_H_generic = beta_H_generic_tilde; 
      beta_FV_generic = beta_FV_generic_tilde; 
        beta_2FV_generic = beta_2FV_generic_tilde;
          beta_2H_generic = beta_2H_generic_tilde;
            beta_FvH_generic = beta_FvH_generic_tilde;
  
  
  tau = tau0*tau_tilde; 	// tau ~ cauchy(0, tau0)
  c2 = slab_scale2*c2_tilde;	// c2 ~ inv_gamma(half_slab_df, half_slab_df*slab_scale2)

  // Compute prior of sparsity approach - this calculation follows equation 2.8 in Piironen and Vehtari 2013

    for(s in 1:S){ // direct interaction INTER plants
      local_shrinkage_ij_tilde[s] = sqrt( c2 * square(local_shrinkage_ij[s]) / (c2 + square(tau) * square(local_shrinkage_ij[s])) );
      alpha_hat_ij[s] = tau * local_shrinkage_ij_tilde[s] * alpha_hat_ij_tilde[s];


      for(k in 1:S){ //HOIs plant-plant 
        beta_plant_local_shrinkage_ijk_tilde[s,k] = sqrt( c2 * square(beta_plant_local_shrinkage_ijk[s,k]) / (c2 + square(tau) * square(beta_plant_local_shrinkage_ijk[s,k])) );
        beta_plant_hat_ijk[s,k] = tau * beta_plant_local_shrinkage_ijk_tilde[s,k] * beta_plant_hat_ijk_tilde[s,k];
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
       
      for(h in 1:H){ // direct interaction plant - H
    gamma_H_local_shrinkage_ih_tilde[h] = sqrt( c2 * square(gamma_H_local_shrinkage_ih[h]) / (c2 + square(tau) * square(gamma_H_local_shrinkage_ih[h])) );
    gamma_H_hat_ih[h] = tau * gamma_H_local_shrinkage_ih_tilde[h] * gamma_H_hat_ih_tilde[h];
    
          for(h2 in 1:H){ //HOIs plant-H-H 
          beta_2H_local_shrinkage_ihh_tilde[h,h2] = sqrt( c2 * square(beta_2H_local_shrinkage_ihh[h,h2]) / (c2 + square(tau) * square(beta_2H_local_shrinkage_ihh[h,h2])) );
          beta_2H_hat_ihh[h,h2] = tau * beta_2H_local_shrinkage_ihh_tilde[h,h2] * beta_2H_hat_ihh_tilde[h,h2];
           }
      }
    
    for(fv in 1:FV){ // direct interaction plant - FV
      gamma_FV_local_shrinkage_if_tilde[fv] = sqrt( c2 * square(gamma_FV_local_shrinkage_if[fv]) / (c2 + square(tau) * square(gamma_FV_local_shrinkage_if[fv])));
      gamma_FV_hat_if[fv] = tau * gamma_FV_local_shrinkage_if_tilde[fv] * gamma_FV_hat_if_tilde[fv];
      
        for(fv2 in 1:FV){ //HOIs plant-FV-FV 
        beta_2FV_local_shrinkage_iff_tilde[fv,fv2] = sqrt( c2 * square(beta_2FV_local_shrinkage_iff[fv,fv2]) / (c2 + square(tau) * square(beta_2FV_local_shrinkage_iff[fv,fv2])) );
        beta_2FV_hat_iff[fv,fv2] = tau * beta_2FV_local_shrinkage_iff_tilde[fv,fv2] * beta_2FV_hat_iff_tilde[fv,fv2];
         }
       
        for(h in 1:H){ //HOIs plant-FV-H 
        beta_FvH_local_shrinkage_ifh_tilde[fv,h] = sqrt( c2 * square(beta_FvH_local_shrinkage_ifh[fv,h]) / (c2 + square(tau) * square(beta_FvH_local_shrinkage_ifh[fv,h])) );
        beta_FvH_hat_ifh[fv,h] = tau * beta_FvH_local_shrinkage_ifh_tilde[fv,h] * beta_FvH_hat_ifh_tilde[fv,h];
          }
      }
 // implement the biological model
  for(i in 1:N){
    lambda_ei[i] = lambdas[1];
    for(s in 1:S){
      alpha_eij[i,s] = (1-Intra[s]) * alpha_generic[1] + Intra[s] * alpha_intra[1] + alpha_hat_ij[s];
        
        for(k in 1:S){ // for all third competing species k in HOIs_ijk, here k = plant species  
        beta_ijk[s,k] = beta_plant_generic[s] + beta_plant_hat_ijk[s,k];
        }
      matrix_beta_ijk[i,s] = sum(beta_ijk[s,].* matrix_HOIs_plant[i,s]);
      
        for(h in 1:H){// for one herbivore species h in beta_H_ijh, here h = herbivor species h and j = plant species
        beta_H_ijh[s,h] = beta_H_generic[s] + beta_H_hat_ijh[s,h];
        }
      matrix_beta_H_ijh[i,s] = sum(beta_H_ijh[s,].* matrix_HOIs_ijh[i,s]);
    
        for(fv in 1:FV){
        beta_F_ijf[s,fv] = exp(beta_FV_generic[s] + beta_FV_hat_ijf[s,fv]);
        }
        matrix_beta_F_ijf[i,s] = sum(beta_F_ijf[s,].* matrix_HOIs_ijf[i,s]);
        }
  
    for(h in 1:H){ // for one herbivore species h in gamma_H_ih, here h = species h
        gamma_H_eih[i,h] = gamma_H[1] + gamma_H_hat_ih[h];
        
        // for two herbivores inbeta_ihh, here h = herbivore species h
        for(h2 in 1:H){
            beta_2H_ihh[h,h2] = beta_2H_generic[h] + beta_2H_hat_ihh[h,h2];
        }
        matrix_beta_2H_ihh[i,h] = sum(beta_2H_ihh[h,].* matrix_HOIs_ihh[i,h]);

      }
      
    for(fv in 1:FV){ // for  floral visitor species f in gamma_FV_if, here f = species f
        gamma_FV_eif[i,fv] = gamma_FV[1] + gamma_FV_hat_if[fv];
        
          // for two floral visitor species fv inbeta_iff, here f = floral visitor species f
        for(fv2 in 1:FV){
            beta_2F_iff[fv,fv2] = beta_2FV_generic[fv] + beta_2FV_hat_iff[fv,fv2];
          }
        matrix_beta_2FV_iff[i,fv] = sum(beta_2F_iff[fv,].* matrix_HOIs_iff[i,fv]);
    
              // for one floral visitor species fv and one herbivore inbeta_ifh, here f = floral visitor species f and h= herbivore species h
        for(h in 1:H){
            beta_FvH_ifh[fv,h] = beta_FvH_generic[fv] + beta_FvH_hat_ifh[fv,h];
        }
            matrix_beta_FvH_ifh[i,fv] = sum(beta_FvH_ifh[fv,].* matrix_HOIs_ifh[i,fv]);
    }


    HOI_effects[i] = sum(matrix_beta_ijk[i,]) + sum(matrix_beta_F_ijf[i,]) +  sum(matrix_beta_H_ijh[i,])  + sum(matrix_beta_2H_ihh[i,]) +  sum(matrix_beta_FvH_ifh[i,]) + sum(matrix_beta_2FV_iff[i,]);
     
    interaction_effects[i] = sum(alpha_eij[i,] .* SpMatrix[i,]) +  sum(gamma_H_eih[i,] .* SpMatrix_H[i,]) + sum( gamma_FV_eif[i,] .* SpMatrix_FV[i,]);
     
 
          F_hat[i] = exp(lambda_ei[i] + interaction_effects[i] + HOI_effects[i]);

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

  
    
  for(s in 1:S){
    beta_plant_generic_tilde[s] ~ normal(0,1);
      // beta_plant_hat_ijk[,s] ~ normal(0,1); // the ensemble of the species-specific beta related to one generic beta follow a normal distribution
    beta_plant_hat_ijk_tilde[s,] ~ normal(0,1);
    beta_plant_local_shrinkage_ijk[s,] ~ cauchy(0,1);
  
    beta_FV_generic_tilde[s] ~ normal(0,1);
    beta_FV_hat_ijf_tilde[s,] ~ normal(0,1);
    beta_FV_local_shrinkage_ijf[s,] ~ cauchy(0,1);
    
    beta_H_generic_tilde[s] ~ normal(0,1); 
    beta_H_hat_ijh_tilde[s,] ~ normal(0,1);
    beta_H_local_shrinkage_ijh[s,] ~ cauchy(0,1);
  }
  for (h in 1:H){
      beta_2H_generic_tilde[h] ~ normal(0,1);
      beta_2H_hat_ihh_tilde[h,] ~ normal(0,1);
      beta_2H_local_shrinkage_ihh[h,] ~ cauchy(0,1);
    }
    
  for (fv in 1:FV){
    beta_2FV_generic_tilde[fv] ~ normal(0,1);
    beta_2FV_hat_iff_tilde[fv,] ~ normal(0,1);
    beta_2FV_local_shrinkage_iff[fv,] ~ cauchy(0,1);
    
    beta_FvH_generic_tilde[fv] ~ normal(0,1);
    beta_FvH_hat_ifh_tilde[fv,] ~ normal(0,1);
    beta_FvH_local_shrinkage_ifh[fv,] ~ cauchy(0,1);
    }

 for(i in 1:N){
  Fecundity[i] ~ neg_binomial_2(F_hat[i],(disp_dev^2)^(-1)); 
   }

}
