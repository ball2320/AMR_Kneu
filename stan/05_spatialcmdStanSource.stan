

//fxn for GPR from https://peter-stewart.github.io/blog/gaussian-process-occupancy-tutorial/
functions{
    matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;
        for (i in 1:(N-1)) {
          K[i, i] = sq_alpha + delta;
          for (j in (i + 1):N) {
            K[i, j] = sq_alpha * exp(-sq_rho * square(x[i,j]) );
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
  //int<lower=1, upper=C> country[N]; // country indicator for each data point
  vector[N] t;     // time points
  //vector[N] y;     // observed values
  int calc_likelihood; //calculate likelihood
  array[N] int<lower=0> y; //observed values as an int for the number of resistant cases 
  array[N] int<lower=1> cases; //the total number of cases as an integer as an input into the bionomial likelihood
  
  //lets test with spatial correlation variables
  int<lower=0> N_edges;             // number of edges
  array[N_edges] int<lower=1, upper=C> node1;  // node1[i], node2[i] are neighbors
  array[N_edges]int<lower=1, upper=C> node2;
  
  // Distance matrix between countries
  matrix[C, C] D;  // Distance matrix between countries
  
  //add in source level RE 
  int<lower=1> S; //number of sources
  array[N] int<lower=1, upper=S> source; // source indicator for each data point


}

parameters {

  //OPTION 1 for I0: Set hirerachial structure for I0 (works better if we start
  //at the first point in our data and works for most part for 2000 start)
  //NEW REPARAM
  //real I0_logit;
  //vector[C] I0_logit_c; // no bounds needed if on logit scale
  //real<lower=0> sigma_I0_logit;
  //real I0_logit_tilde[C];
  

  
  //OPTION 2 for I0: Set prior for I0  with no hirerachial structure (works better 
// if we start at previous year with no data ex 2000 start)
  real<lower=0, upper=1> I0_prob;
 


  
  //reparam with beta gamma I0 and R0
  //real<lower=0> beta;
  // vector<lower=0>[C] beta_c;
  // real<lower=0> sigma_beta;
  //real beta_tilde[C];

  //real<lower=0> gamma;
  // vector<lower=0>[C] gamma_c;
  // real<lower=0> sigma_gamma;
  //real gamma_tilde[C];

  //Spatial Correlation Parameters
  //vector[C] phi;                    // spatial effects
  //real<lower=0> tau_phi; // precision of spatial effects
  real log_beta;
  real log_gamma;
  
  // //TEST SPATIAL
  // real<lower=0> tau_phi_I0;
  //real<lower=0> tau_phi_beta;
  
  // real<lower=0> sigma_phi_I0;
  // real<lower=0> sigma_phi_beta;
  // real<lower=0> sigma_phi_gamma;
  //real eta_I0;
  // real eta_beta;
  // real eta_gamma;
  // vector[C] phi_I0;
  // vector[C] phi_beta;
  //vector[C] phi_gamma;
  //real<lower=0> tau_phi_gamma;
  //real<lower=0> tau_phi_I0;


  // Gaussian process parameters
  vector[C] z; // z-scores for intercept term (for non-centred parameterisation)
  real<lower=0> etasq; // Maximum covariance between sites
  real<lower=0> rhosq; // Rate of decline in covariance with distance
  
  vector[C] z_I0; // z-scores for intercept term (for non-centred parameterisation)
  real<lower=0> etasq_I0; // Maximum covariance between sites
  //real<lower=0> rhosq_I0; // Rate of decline in covariance with distance
 
  vector[C] z_gamma; // z-scores for intercept term (for non-centred parameterisation)
  real<lower=0> etasq_gamma; // Maximum covariance between sites
   //real<lower=0> rhosq_gamma;
   
  // vector[C] z_gamma; // z-scores for intercept term (for non-centred parameterisation)
  // real<lower=0> etasq_gamma; // Maximum covariance between sites
  // real<lower=0> rhosq_gamma; // Rate of decline in covariance with distance

 //source
 vector<lower=0, upper=1>[S] source_eta_S;
  
}


transformed parameters {

  //real<lower=0> sigma_phi = inv(sqrt(tau_phi)); // convert precision to sigma
  //TEST SPATIAL
  //real<lower=0> sigma_phi_I0 = inv(sqrt(tau_phi_I0));
  // real<lower=0> sigma_phi_beta = inv(sqrt(tau_phi_beta));
  //real<lower=0> sigma_phi_gamma = inv(sqrt(tau_phi_gamma));

  matrix[C, C] L_SIGMA; // Cholesky-decomposed covariance matrix
  matrix[C, C] SIGMA; // Covariance matrix
  vector[C] k; // Intercept term for each site (perturbation from k_bar)
  //vector[C] k_temp;

  // Gaussian process - non-centred
  SIGMA = cov_GPL2(D, etasq, rhosq, 0.05);
  L_SIGMA = cholesky_decompose(SIGMA);
  k = L_SIGMA * z;
  // k_temp = L_SIGMA * z;
  // k = k_temp - mean(k_temp);
 
  
  matrix[C, C] L_SIGMA_I0; // Cholesky-decomposed covariance matrix
  matrix[C, C] SIGMA_I0; // Covariance matrix
  vector[C] k_I0; // Intercept term for each site (perturbation from k_bar)
  //vector[C] k_I0_temp;

  // Gaussian process - non-centred
  //SIGMA_I0 = cov_GPL2(D, etasq_I0, rhosq_I0, 0.05);
  SIGMA_I0 = cov_GPL2(D, etasq_I0, rhosq, 0.05);
  L_SIGMA_I0 = cholesky_decompose(SIGMA_I0);
  k_I0 = L_SIGMA_I0 * z_I0;
  // k_I0_temp = L_SIGMA_I0 * z_I0;
  // k_I0 = k_I0_temp - mean(k_I0_temp);
  
  //
  matrix[C, C] L_SIGMA_gamma; // Cholesky-decomposed covariance matrix
  matrix[C, C] SIGMA_gamma; // Covariance matrix
  vector[C] k_gamma; // Intercept term for each site (perturbation from k_bar)
  //vector[C] k_I0_temp;

  // Gaussian process - non-centred
  //SIGMA_gamma = cov_GPL2(D, etasq_gamma, rhosq_gamma, 0.05);
  SIGMA_gamma = cov_GPL2(D, etasq_gamma, rhosq, 0.05);
  L_SIGMA_gamma = cholesky_decompose(SIGMA_gamma);
  k_gamma = L_SIGMA_gamma * z_gamma;

  // 
  // matrix[C, C] L_SIGMA_gamma; // Cholesky-decomposed covariance matrix
  // matrix[C, C] SIGMA_gamma; // Covariance matrix
  // vector[C] k_gamma; // Intercept term for each site (perturbation from k_bar)
  // 
  // // Gaussian process - non-centred
  // SIGMA_gamma = cov_GPL2(D, etasq_gamma, rhosq_gamma, 0.01);
  // L_SIGMA_gamma = cholesky_decompose(SIGMA_gamma);
  // k_gamma = L_SIGMA_gamma * z_gamma;

  

  
}


model {

    //OPTION 1 for I0:
  //I0_logit ~ normal(0,1);
  //sigma_I0_logit ~ exponential(0.5);
  //I0_logit_c ~ normal(I0_logit,sigma_I0_logit);
  //I0_logit_tilde ~ normal(0,1);
  
  //OPTION 2 for I0:
  I0_prob ~ beta(2,2);
  
  log_beta ~ normal(0,1);
  // sigma_beta ~ exponential(0.5);
  // log_beta_c ~ normal(log_beta,sigma_beta);
  
  log_gamma ~ normal(0,1);

  //gamma ~ exponential(3);
  // sigma_gamma ~ exponential(0.5);
  // gamma_c ~ normal(gamma, sigma_gamma);
  
  // target += -0.5*dot_self(phi[node1] - phi[node2]);
  // sum(phi) ~ normal(0, 0.001 * C);
  // tau_phi ~ gamma(1, 0.1);
  
  // //TEST FOR SPATIAL
  // tau_phi_I0 ~ gamma(1, 0.1);
  // tau_phi_beta ~ gamma(1, 0.1);
 
  // 
  // sigma_phi_I0 ~ normal(0,0.1);
  // sigma_phi_beta ~ normal(0,1);
  // sigma_phi_gamma ~ normal(0,1);
  //eta_I0 ~ normal(0,0.1);
  // eta_beta ~ normal(0,0.5);
  // eta_gamma ~ normal(0,0.5);
  
  // target += -0.5*dot_self(phi_I0[node1] - phi_I0[node2])*tau_phi_I0;
  // sum(phi_I0) ~ normal(0, 0.001 * C);
  // 
  // target += -0.5*dot_self(phi_beta[node1] - phi_beta[node2])*tau_phi_beta;
  // sum(phi_beta) ~ normal(0, 0.001 * C);
  // 
  // target += -0.5*dot_self(phi_gamma[node1] - phi_gamma[node2]);
  // sum(phi_gamma) ~ normal(0, 0.01 * C);
  //tau_phi_gamma ~ gamma(1, 0.1);
  //tau_phi_I0 ~ gamma(1, 0.1);
  
  
  rhosq ~ exponential(0.5); 
  etasq ~ exponential(1);
  z ~ normal(0, 1);
  //sum(z) ~ normal(0, 0.001 * sqrt(C));
  
  //rhosq_I0 ~ exponential(1);
  etasq_I0 ~ exponential(1); //3
  z_I0 ~ normal(0, 1);
  //sum(z_I0) ~ normal(0, 0.001 * sqrt(C));
  
  etasq_gamma ~ exponential(1); //3
  z_gamma ~ normal(0, 1);
  //rhosq_gamma ~ exponential(1);
  // etasq_gamma ~ exponential(3);
  // z_gamma ~ normal(0, 1);

  //source
  source_eta_S ~ beta(5,1);


if (calc_likelihood == 1){
  
  vector[N] prob;
  
  int lastCountry = -1;
  real lastYear = -9999;
  real lastI = -9999;
  real K_c_country = -9999;
  // Likelihood
  for(i in 1:N){
   
    //can change this so id you have unordered data it should still work or signpost better that the data has to be ordered!!!!
    if(country[i] != lastCountry){
        //prob[i] = inv_logit(I0_logit_c[country[i]]);
        
        //calculate resistance rate at 2000 and use that to calculate the resistance rate at that first train year

          real prob_I0 = inv_logit(logit(I0_prob) + k_I0[country[i]]);
          real current_beta = exp(log_beta + k[country[i]]);
          real current_gamma = exp(log_gamma + k_gamma[country[i]]);
          
          prob[i] = ((current_beta - current_gamma) / current_beta ) / (1 + (((current_beta - current_gamma) / current_beta )/prob_I0 - 1) * exp(-(current_beta - current_gamma)*(t[i])));
        
          
          //forward
         //prob[i] = inv_logit(I0_logit + k_I0[country[i]]);
         
        
    }else{
     

      
      //analytical beta and gamma, i(t) =((β - γ) / β ) / (1 + (((β - γ) / β )/i₀ - 1) * e^(-(β - γ)t))

      
      //TEST SPATIAL AGAIN
      real current_beta = exp(log_beta + k[country[i]]); 
      
      real current_gamma = exp(log_gamma + k_gamma[country[i]]);

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

    
  
  real source_prob = (prob[i])*source_eta_S[source[i]];

 
 //y[i] ~ binomial(cases[i], fmax( fmin(prob[i], 0.99999 ),0.00001));
 y[i] ~ binomial(cases[i], fmax( fmin(source_prob, 0.99999 ),0.00001));
 
  }
}
}

//test with temp eqn
generated quantities {
  vector[N] y_pred;
  //vector[C] r_c;
  vector[N] log_lik;
  //vector[N] y_pred_counts;
  //vector[N] y_pred_prop;


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
        
          
          //forward pred 
          //y_pred[i] = inv_logit(I0_logit + k_I0[country[i]]);

    }else{


      real current_beta = exp(log_beta + k[country[i]]);  //exp(log_beta + sigma_phi*phi[country[i]]); 
      //real current_gamma = exp(log_gamma + k[country[i]]); //gamma;//exp(log_gamma + sigma_phi_gamma*phi[country[i]]); 
      real current_gamma = exp(log_gamma + k_gamma[country[i]]); 
      if(t[i] - lastYear < pow(10,-4)){
        y_pred[i] = lastI;
      }else{
      y_pred[i] = ((current_beta - current_gamma) / current_beta ) / (1 + (((current_beta - current_gamma) / current_beta )/lastI - 1) * exp(-(current_beta - current_gamma)*(t[i] - lastYear)));
        }

    }

    // y_pred_counts[i] = binomial_rng(cases[i],fmax(fmin(y_pred[i], 0.999999 ),0.000001));
    // y_pred_prop[i] = y_pred_counts[i]/cases[i];

    lastI = y_pred[i];
    lastCountry = country[i];
    lastYear = t[i];
    
    real source_prob = (y_pred[i])*source_eta_S[source[i]];

    // //y_pred_counts[i] = binomial_rng(cases[i],fmax(fmin(y_pred[i], 0.999999 ),0.000001));
    // y_pred_counts[i] = binomial_rng(cases[i],fmax(fmin(source_prob, 0.999999 ),0.000001));
    // y_pred_prop[i] = y_pred_counts[i]/cases[i];
    
    log_lik[i] = binomial_lpmf(y[i]|cases[i], fmax(fmin(source_prob , 0.999999 ),0.000001));
  }
}
