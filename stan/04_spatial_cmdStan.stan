//fxn for GPR from https://peter-stewart.github.io/blog/gaussian-process-occupancy-tutorial/
functions{
    matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;
        for (i in 1:(N-1)) {
          K[i, i] = sq_alpha + delta;
          for (j in (i + 1):N) {
            K[i, j] = sq_alpha * exp(-sq_rho * x[i,j] );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = sq_alpha + delta;
        return K;
    }
}


//test with random effects for each source for growth rate, intercept, and initial resistance (TO DO EDIT FOR RE at country level)
data {
  int<lower=1> N;  // number of observations
  int<lower=1> C; //number of countries
  array[N] int<lower=1, upper=C> country;
  vector[N] t;     // time points
  int calc_likelihood; //calculate likelihood
  array[N] int<lower=0> y; //observed values as an int for the number of resistant cases 
  array[N] int<lower=1> cases; //the total number of cases as an integer as an input into the bionomial likelihood
  
  //add test with spatial correlation variables
  int<lower=0> N_edges;             // number of edges
  array[N_edges] int<lower=1, upper=C> node1;  // node1[i], node2[i] are neighbors
  array[N_edges]int<lower=1, upper=C> node2;
  
  // Distance matrix between countries
  matrix[C, C] D;  // Distance matrix between countries
  

}

parameters {
  
  //Initial Resistance - non hierarchical contrained btw 0 and 1 (works better 
  // if we start at previous year with no data ex 2000 start)
  real<lower=0, upper=1> I0_prob;
  
  
  //reparam with beta gamma I0
  real log_beta;
  real log_gamma;
  
  // Gaussian process parameters
  vector[C] z; // z-scores for intercept term (for non-centred parameterisation)
  //real<lower=0> etasq; // Maximum covariance between sites
  real<lower=1e-60> etasq;
  //real<lower=0> rhosq; // Rate of decline in covariance with distance
  real<lower=1e-6> rhosq;
  
  vector[C] z_I0; // z-scores for intercept term (for non-centred parameterisation)
  real<lower=0> etasq_I0; // Maximum covariance between sites
 
  vector[C] z_gamma; // z-scores for intercept term (for non-centred parameterisation)
  real<lower=0> etasq_gamma; // Maximum covariance between sites 

  
}

transformed parameters {

  matrix[C, C] L_SIGMA; // Cholesky-decomposed covariance matrix
  matrix[C, C] SIGMA; // Covariance matrix
  vector[C] k; // Intercept term for each site (perturbation from k_bar)
  
  vector[C] k_I0; // Intercept term for each site (perturbation from k_bar)
  vector[C] k_gamma; // Intercept term for each site (perturbation from k_bar)

  // Gaussian process - non-centred
  //SIGMA = cov_GPL2(D2, etasq, rhosq, 0.05); //ball modified
  //L_SIGMA = cholesky_decompose(SIGMA); //ball modified
  //k = L_SIGMA * z; //ball modified
 
  SIGMA = cov_GPL2(D, 1, rhosq, 0.05);
  L_SIGMA = cholesky_decompose(SIGMA);
  
  k = sqrt(etasq) * (L_SIGMA * z);
  k_I0 = sqrt(etasq_I0) * (L_SIGMA * z_I0);
  k_gamma = sqrt(etasq_gamma) * (L_SIGMA * z_gamma);
   
  
  //matrix[C, C] L_SIGMA_I0; // Cholesky-decomposed covariance matrix
  //matrix[C, C] SIGMA_I0; // Covariance matrix
  //vector[C] k_I0; // Intercept term for each site (perturbation from k_bar)

  // Gaussian process - non-centred
  //SIGMA_I0 = cov_GPL2(D2, etasq_I0, rhosq, 0.05); //ball modified
  //L_SIGMA_I0 = cholesky_decompose(SIGMA_I0); //ball modified
  //k_I0 = L_SIGMA_I0 * z_I0; //ball modified
  
  //
  //matrix[C, C] L_SIGMA_gamma; // Cholesky-decomposed covariance matrix
  //matrix[C, C] SIGMA_gamma; // Covariance matrix
  //vector[C] k_gamma; // Intercept term for each site (perturbation from k_bar)

  // Gaussian process - non-centred
  //SIGMA_gamma = cov_GPL2(D2, etasq_gamma, rhosq, 0.05); //ball modified
  //L_SIGMA_gamma = cholesky_decompose(SIGMA_gamma); //ball modified
  //k_gamma = L_SIGMA_gamma * z_gamma; //ball modified

}


model {
  
  //Initial Resistance
  I0_prob ~ beta(2,2);


  //reparam with beta gamma and I0
  log_beta ~ normal(0,1);
  log_gamma ~ normal(0,1);

  //rhosq ~ exponential(0.5); //determines the rate of decay between neighbours
  //etasq ~ exponential(1); // determines maximum covariance between two societies
  rhosq ~ lognormal(0,1); // pause for testing additional priors
  etasq ~ lognormal(0,1); // pause for testing additional priors
  
  
  z ~ normal(0, 1);
  
  //etasq_I0 ~ exponential(1);
  etasq_I0 ~ lognormal(0,1); // pause for testing additional priors
  z_I0 ~ normal(0, 1);
  
  //etasq_gamma ~ exponential(1);
  etasq_gamma ~ lognormal(0,1); // pause for testing additional priors
  z_gamma ~ normal(0, 1);


if (calc_likelihood == 1){
  
  vector[N] prob;
  
  int lastCountry = -1;
  real lastYear = -9999;
  real lastI = -9999;
  real K_c_country = -9999;
  
  // Likelihood
  for(i in 1:N){
   
    //can change this so if you have unordered data it should still work or signpost better that the data has to be ordered!!!!
    if(country[i] != lastCountry){
        
        
          //resistance rate at 2000 or whatever start year and use that to calculate the resistance rate at that first train year
          real prob_I0 = inv_logit(logit(I0_prob) + k_I0[country[i]]);  //use this for option 2 initial resistance 

          real current_beta = exp(log_beta + k[country[i]]);
          real current_gamma = exp(log_gamma + k_gamma[country[i]]);

          prob[i] = ((current_beta - current_gamma) / current_beta ) / (1 + (((current_beta - current_gamma) / current_beta )/prob_I0 - 1) * exp(-(current_beta - current_gamma)*(t[i])));

        
        
    }else{
     
      
      real current_beta = exp(log_beta + k[country[i]]);
      real current_gamma = exp(log_gamma + k_gamma[country[i]]);

      //handles repeated years 
      if(t[i] - lastYear < pow(10,-4) ){
        prob[i] = lastI;
      }else{
      //in terms of gamma, beta, and i0 ((β - γ) / β ) / (1 + (((β - γ) / β )/i₀ - 1) * e^(-(β - γ)t))
      prob[i] = ((current_beta - current_gamma) / current_beta ) / (1 + (((current_beta - current_gamma) / current_beta )/lastI - 1) * exp(-(current_beta - current_gamma)*(t[i] - lastYear)));
      }
      

      
      
    }
    
    lastI = prob[i];
    lastCountry = country[i];
    lastYear = t[i];
  
 
 y[i] ~ binomial(cases[i], fmax( fmin(prob[i], 0.99999 ),0.00001));
 
  }
}
}

//test with temp eqn - should just be similar structure to previous block 
generated quantities {
  vector[N] y_pred;
  vector[N] log_lik;

  
  int lastCountry = -1;
  real lastYear = -9999;
  real lastI = -9999;
  real K_c_country = -9999;
  for (i in 1:N) {
    
    if(country[i] != lastCountry){
      
       
         //calculate resistance rate at 2000 and use that to calculate the resistance rate at that first train year
          real prob_I0 = inv_logit(logit(I0_prob) + k_I0[country[i]]);
          real current_beta = exp(log_beta + k[country[i]]);
          real current_gamma = exp(log_gamma + k_gamma[country[i]]);

          y_pred[i] = ((current_beta - current_gamma) / current_beta ) / (1 + (((current_beta - current_gamma) / current_beta )/prob_I0 - 1) * exp(-(current_beta - current_gamma)*(t[i])));

         
       
    }else{

      
      
      real current_beta = exp(log_beta + k[country[i]]);
      real current_gamma = exp(log_gamma + k_gamma[country[i]]);
     
      if(t[i] - lastYear < pow(10,-4)){
        y_pred[i] = lastI;
      }else{
      y_pred[i] = ((current_beta - current_gamma) / current_beta ) / (1 + (((current_beta - current_gamma) / current_beta )/lastI - 1) * exp(-(current_beta - current_gamma)*(t[i] - lastYear)));
        }

    } 
    
    

    
    lastI = y_pred[i];
    lastCountry = country[i];
    lastYear = t[i];
    
    log_lik[i] = binomial_lpmf(y[i]|cases[i], fmax(fmin(y_pred[i], 0.999999 ),0.000001));
  }
}
