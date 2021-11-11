// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

data{
  int<lower = 1> N;  // Number of observations
  int<lower = 1> S;  // Number of species
  int Fecundity[N];  // Fecundity of the focal species in each plot
  matrix[N,S] SpMatrix;  // Matrix of abundances for each species (including abundances of non-focal individuals of the focal species)
  //vector[N] env;  // Environmental values for each plot
  int<lower = 0> Intra[S];  // Indicator boolean variable to identify the focal species (0 for non-focal and 1 for focal). Included for easier calculations
  int Inclusion_ij[S];  // Boolean indicator variables to identify the species x reserve intercept parameters identified for inclusion in the final model
  int beta_Inclusion[S,S];  // Boolean indicator variables to identify the species x reserve intercept parameters identified for inclusion in the final model
}

parameters{
  vector[2] lambdas;
  vector[2] alpha_generic_tilde;
  vector[2] alpha_intra_tilde;
  vector[S] alpha_hat_ij;

   vector[2] beta_generic_tilde;
   matrix[S,S] beta_hat_ij;
}

transformed parameters{
  vector[2] alpha_generic;
  vector[2] alpha_intra;
  
  vector[2] beta_generic;
  // scale the alpha values
  alpha_generic[1] = 3 * alpha_generic_tilde[1] - 6;
  alpha_intra[1] = 3 * alpha_intra_tilde[1] - 6;
  alpha_generic[2] = 0.5 * alpha_generic_tilde[2];
  alpha_intra[2] = 0.5 * alpha_intra_tilde[2];
  
  beta_generic[1] = 3 * beta_generic_tilde[1] - 6;
    beta_generic[2] = 0.5 * beta_generic_tilde[2];

}

model{
  // Declare objects necessary for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //     a matrix of the species specific alpha values for each species and plot (interaction_effects), and a matrix
  //     of the the alpha*N values for each species.
  vector[N] F_hat;
  vector[N] interaction_effects;
  vector[N] HOI_effects;
  matrix[N,S] alpha_eij;
  matrix[N,S] beta_eij;
  vector[N] lambda_ei;

matrix[1,S] Spmatrixi; // to transform spmatrix to matrix S by S with all species interactions for each obs
matrix[N,S] matrix_beta_eij;
matrix[S,S] matrix_HOIs;

  // set regular priors
  alpha_generic_tilde ~ normal(0,1);
  alpha_intra_tilde ~ normal(0,1);
  beta_generic_tilde ~ normal(0,1);

    lambdas ~ normal(0, 1);
    alpha_hat_ij ~ normal(0,1);
    
    for (s in 1:S){
    beta_hat_ij[s,] ~ normal(0,1);
    }

  // implement the biological model
  for(i in 1:N){
    lambda_ei[i] = exp(lambdas[1]);
    Spmatrixi[1,] = SpMatrix[i,];
      for (n in 1:S) {
          for (m in 1:S) {
            if (m <= n){
              matrix_HOIs[n,m] = 0;
            }
        matrix_HOIs[n,m] = Spmatrixi[1,n]* Spmatrixi[1,m]; 
        }
        }
    for(s in 1:S){
      alpha_eij[i,s] = exp((1-Intra[s]) * alpha_generic[1] + Intra[s] * alpha_intra[1] + Inclusion_ij[s] * alpha_hat_ij[s]);
      
      for(k in 1:S){ // for all third competing species k in HOIs_ijk, here k = species k 

        beta_eij[s,k] = exp(beta_generic[1] + beta_Inclusion[s]*beta_hat_ij[s,k]) ;
        
        }
        
        
        matrix_beta_eij[i,s] = sum(beta_eij[s,] .* matrix_HOIs[s,]);
    }
    HOI_effects[i] = sum(matrix_beta_eij[i,]);
    interaction_effects[i] = sum(alpha_eij[i,] .* SpMatrix[i,]);
    F_hat[i] = lambda_ei[i] / (1 + interaction_effects[i] + HOI_effects[i]);
  }
  Fecundity ~ poisson(F_hat);
}
