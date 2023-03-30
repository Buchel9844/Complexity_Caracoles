// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

data{
  int<lower = 1> N; // Number of plots
  int<lower = 1> S; // Number of plant species
  int<lower = 1> H; // Number of herbivores species
  int<lower = 1> FV; // Number of floral visitors species
  int<lower = 1> FvH_H; // Number of herbivores species
  int<lower = 1> FvH_FV; // Number of floral visitors species
    int RemoveFvH; // Remove Higher trophic level 

  int<lower = 0> run_estimation;

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
  matrix[FvH_FV,FvH_H]  matrix_HOIs_ifh[N]; // Matrix of abundances for each floral visitors species with each herbivores species

   // vector telling which interactions to include
   
  int Inclusion_ij[1,S];  // Boolean indicator variables to identify the plant species
  int Inclusion_FV[1,FV]; // Boolean indicator variables to identify the floral visitor species
  int Inclusion_H[1,H]; // Boolean indicator variables to identify the herbivore species
  
  int beta_Inclusion_plant[S,S];  // Boolean indicator variables to identify the HOIs plant-plant
  int beta_Inclusion_FV[S,FV];  // Boolean indicator variables to identify the HOIs 2 plants-FV
  int beta_Inclusion_H[S,H];  // Boolean indicator variables to identify the HOIs 2 plants- H

  int beta_Inclusion_2FV[FV,FV];  // Boolean indicator variables to identify the HOIs 1 plant - 2 FV
  int beta_Inclusion_2H[H,H];  // Boolean indicator variables to identify the HOIs 1 plant - 2 H
  int beta_Inclusion_FvH[FV,H];  // Boolean indicator variables to identify the HOIs 1 plant - 2 H


}

parameters{
  real<lower=0> lambdas[1];
  real<lower=-5,upper=5> alpha_generic_tilde[1];
  real<lower=-5,upper=5> alpha_intra_tilde[1];
  vector<lower=-5,upper=5>[S] alpha_hat_ij;
  
  vector<lower=-5,upper=5>[S] beta_plant_generic_tilde; // HOIs plants
  matrix<lower=-5,upper=5>[S,S] beta_plant_hat_ijk; // HOIs plants

  real<lower=-5,upper=5> gamma_H_generic_tilde[1]; // direct interaction plants - herbivores ; generic
  vector<lower=-5,upper=5>[H] gamma_H_hat_ih; // direct interaction plants - herbivores ; species -specific term

  real<lower=-5,upper=5> gamma_FV_generic_tilde[1]; // direct interaction plants - FV ; generic
  vector<lower=-5,upper=5>[FV] gamma_FV_hat_if; // direct interaction plants - FV ; species -specific term

  vector<lower=-5,upper=5>[S] beta_H_generic_tilde; // HOIs herbivores
  matrix<lower=-5,upper=5>[S,H] beta_H_hat_ijh; // HOIs herbivores 

  vector[S] beta_FV_generic_tilde; // HOIs floral visitors
  matrix[S,FV] beta_FV_hat_ijf; // HOIs floral visitors

  vector<lower=-5,upper=5>[FV] beta_2FV_generic_tilde; // HOIs 2 floral visitors
  matrix<lower=-5,upper=5>[FV,FV] beta_2FV_hat_iff; // HOIs 2 floral visitors

  vector<lower=-5,upper=5>[H] beta_2H_generic_tilde; // HOIs 2 herbivores
  matrix<lower=-5,upper=5>[H,H] beta_2H_hat_ihh; // HOIs 2 herbivores 

  vector<lower=-5,upper=5>[FvH_FV] beta_FvH_generic_tilde; // HOIs 1 floral visitor & 1 herbivore
  matrix<lower=-5,upper=5>[FvH_FV,FvH_H] beta_FvH_hat_ifh; // HOIs 1 floral visitor & 1 herbivore


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
  matrix[FvH_FV,FvH_H] beta_FvH_ifh; 
    matrix[N,FvH_FV] matrix_beta_FvH_ifh;
    
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
  vector[FvH_FV] beta_FvH_generic; //HOIs 1 floral visitor & 1 herbivore
  
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
  
  
 // implement the biological model
  for(i in 1:N){
    lambda_ei[i] = lambdas[1];
    for(s in 1:S){
      alpha_eij[i,s] = (1-Intra[s]) * alpha_generic[1] + Intra[s] * alpha_intra[1] + (1-Intra[s]) *alpha_hat_ij[s]*Inclusion_ij[1,s];
        
        for(k in 1:S){ // for all third competing species k in HOIs_ijk, here k = plant species  
        beta_ijk[s,k] = beta_plant_generic[s] + beta_Inclusion_plant[s,k]* beta_plant_hat_ijk[s,k];
        }
      matrix_beta_ijk[i,s] = sum(beta_ijk[s,].* matrix_HOIs_plant[i,s]);
       if(RemoveFvH ==1){
        for(h in 1:H){// for one herbivore species h in beta_H_ijh, here h = herbivor species h and j = plant species
        beta_H_ijh[s,h] = beta_H_generic[s] + beta_Inclusion_H[s,h]*beta_H_hat_ijh[s,h];
        }
      matrix_beta_H_ijh[i,s] = sum(beta_H_ijh[s,].* matrix_HOIs_ijh[i,s]);
   
        for(fv in 1:FV){
        beta_F_ijf[s,fv] = exp(beta_FV_generic[s] + beta_Inclusion_FV[s,fv]*beta_FV_hat_ijf[s,fv]);
        }
        matrix_beta_F_ijf[i,s] = sum(beta_F_ijf[s,].* matrix_HOIs_ijf[i,s]);
       }
        }
  if(RemoveFvH ==1){
    for(h in 1:H){ // for one herbivore species h in gamma_H_ih, here h = species h
        gamma_H_eih[i,h] = gamma_H[1] + Inclusion_H[1,h]*gamma_H_hat_ih[h];
        
        // for two herbivores inbeta_ihh, here h = herbivore species h
        for(h2 in 1:H){
            beta_2H_ihh[h,h2] = beta_2H_generic[h] + beta_Inclusion_2H[h,h2].*beta_2H_hat_ihh[h,h2];
        }
        matrix_beta_2H_ihh[i,h] = sum(beta_2H_ihh[h,].* matrix_HOIs_ihh[i,h]);

      }
      
    for(fv in 1:FV){ // for  floral visitor species f in gamma_FV_if, here f = species f
        gamma_FV_eif[i,fv] = gamma_FV[1] + Inclusion_FV[1,fv]*gamma_FV_hat_if[fv];
        
          // for two floral visitor species fv inbeta_iff, here f = floral visitor species f
        for(fv2 in 1:FV){
            beta_2F_iff[fv,fv2] = beta_2FV_generic[fv] + beta_Inclusion_2FV[fv,fv2]*beta_2FV_hat_iff[fv,fv2];
          }
        matrix_beta_2FV_iff[i,fv] = sum(beta_2F_iff[fv,].* matrix_HOIs_iff[i,fv]);
    
    }
    
    // for one floral visitor species fv and one herbivore inbeta_ifh, here f = floral visitor species f and h= herbivore species h
         for(fv in 1:FvH_FV){
           for(h in 1:FvH_H){
            beta_FvH_ifh[fv,h] = beta_FvH_generic[fv] + beta_Inclusion_FvH[fv,h]*beta_FvH_hat_ifh[fv,h];
        }
            matrix_beta_FvH_ifh[i,fv] = sum(beta_FvH_ifh[fv,].* matrix_HOIs_ifh[i,fv]);
    }
  

    HOI_effects[i] = sum(matrix_beta_ijk[i,]);
     
    interaction_effects[i] = sum(alpha_eij[i,] .* SpMatrix[i,]);
     
    }else{
    HOI_effects[i] = sum(matrix_beta_ijk[i,]) + sum(matrix_beta_F_ijf[i,]) +  sum(matrix_beta_H_ijh[i,])  + sum(matrix_beta_2H_ihh[i,]) +  sum(matrix_beta_FvH_ifh[i,]) + sum(matrix_beta_2FV_iff[i,]);
     
    interaction_effects[i] = sum(alpha_eij[i,] .* SpMatrix[i,]) +  sum(gamma_H_eih[i,] .* SpMatrix_H[i,]) + sum( gamma_FV_eif[i,] .* SpMatrix_FV[i,]);
    }
 
          F_hat[i] = exp(lambda_ei[i] + interaction_effects[i] + HOI_effects[i]);

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

  for(s in 1:S){
      beta_plant_generic_tilde[s] ~ normal(0,1);
  beta_plant_hat_ijk[,s] ~ normal(0,1); // the ensemble of the species-specific beta related to one generic beta follow a normal distribution
 
     beta_FV_generic_tilde[s] ~ normal(0,1);
  beta_FV_hat_ijf[s,] ~ normal(0,1);
    
    beta_H_generic_tilde[s] ~ normal(0,1); 
    beta_H_hat_ijh[s,] ~ normal(0,1);
  }
  for (h in 1:H){
      beta_2H_generic_tilde[h] ~ normal(0,1);
      beta_2H_hat_ihh[h,] ~ normal(0,1);
    }
    
  for (fv in 1:FV){
    beta_2FV_generic_tilde[fv] ~ normal(0,1);
    beta_2FV_hat_iff[fv,] ~ normal(0,1);
    }
    
    for (fv in 1:FvH_FV){
    beta_FvH_generic_tilde[fv] ~ normal(0,1);
    beta_FvH_hat_ifh[fv,] ~ normal(0,1);
    }
    
 for(i in 1:N){
  Fecundity[i] ~ neg_binomial_2(F_hat[i],(disp_dev^2)^(-1)); 
   }

}
generated quantities{
  vector[N] F_sim;
    if(run_estimation==1){
 for(i in 1:N){
    if(Fecundity[i] <= 0){
      F_sim[i] = 0;
    }
    F_sim[i] = neg_binomial_2_rng(  Fecundity[i],(disp_dev^2)^(-1));
              }
    }
}
