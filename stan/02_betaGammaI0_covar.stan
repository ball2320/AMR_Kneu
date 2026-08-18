//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

//test with random effects for each source for growth rate, intercept, and initial resistance (TO DO EDIT FOR RE at country level)
data {
  int<lower=1> N;  // number of observations
  int<lower=1> C; //number of countries
  int<lower=1, upper=C> country[N]; // country indicator for each data point
  vector[N] t;     // time points
  //vector[N] y;     // observed values
  int calc_likelihood; //calculate likelihood
  int<lower=0> y[N]; //observed values as an int for the number of resistant cases 
  int<lower=1> cases[N]; //the total number of cases as an integer as an input into the bionomial likelihood
  
  
  int<lower=1> K; //number of parameters we are interested in getting so here is 2, intercept and slope 
  matrix[N, K] cov_mat; //covariate matrix 

}

parameters {

  
  //OPTION 1 for I0: Set hirerachial structure for I0 (works better if we start
  //at the first point in our data and works for most part for 2000 start)
  // real I0_logit;
  // vector[C] I0_logit_c; // no bounds needed if on logit scale
  // real<lower=0> sigma_I0_logit;
  // // // //real I0_logit_tilde[C];
  
  //OPTION 2 for I0: Set prior for I0  with no hirerachial structure (works better 
// if we start at previous year with no data ex 2000 start)
  vector<lower=0, upper=1>[C] I0_prob_c;
  

  // //reparam with beta gamma I0 and R0
  // real<lower=0> beta;
  // vector<lower=0>[C] beta_c;
  // real<lower=0> sigma_beta;
  // //real beta_tilde[C];
  
  real log_beta;
  vector[C] log_beta_c;
  real<lower=0> sigma_beta;
  //real beta_tilde[C];
  //vector[C] log_beta_c_raw;

  real<lower=0> gamma;
  vector<lower=0>[C] gamma_c;
  real<lower=0> sigma_gamma;
  //real gamma_tilde[C];
  
  //multivar mat for covariates
  //real<lower=0> sigma_eta;
  //vector<lower=0>[K] eta; // intercept and slope parameters 
  vector[K] eta;
  //TEsting with positive eta forced for eta forabx consumption 
  
  //testing with different covar for the two parameters
  //vector[K] eta_growth; // intercept and slope parameters
  //vector[K] eta_K; // intercept and slope parameters

  
}


model {

  //OPTION 1 for I0:
  // I0_logit ~ normal(0,1);
  // sigma_I0_logit ~ exponential(2); //0.5 //3
  // I0_logit_c ~ normal(I0_logit,sigma_I0_logit);
  // // // //I0_logit_tilde ~ normal(0,1);
  
  //OPTION 2 for I0:
   I0_prob_c ~ beta(2,2);
   

  // //reparam with beta gamma and I0
  // beta ~ exponential(3);
  // sigma_beta ~ exponential(0.5);
  // beta_c ~ normal(beta,sigma_beta);
  // //beta_tilde ~ normal(0,1);
  
  log_beta ~ normal(0,1);
  sigma_beta ~ exponential(0.5);
  log_beta_c ~ normal(log_beta,sigma_beta);
  //log_beta_c_raw ~ normal(0, 1);

  gamma ~ exponential(3);
  sigma_gamma ~ exponential(0.5);
  gamma_c ~ normal(gamma, sigma_gamma);
  
  //coeff covar
  eta ~ normal(0,0.10); 

  
  //eta_growth ~ normal(0,0.1);
  //eta_K ~ normal(0,0.1);


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
        
        
        //calculate resistance rate from set year and use that to calculate the resistance rate at that first train year
          //real prob_I0 = inv_logit(I0_logit_c[country[i]]);
          real prob_I0 = I0_prob_c[country[i]];

          real current_beta =  exp(log_beta_c[country[i]] + cov_mat[i, ] * eta);
          //real current_beta =  beta_c[country[i]];
          real current_gamma = gamma_c[country[i]];
          real growth_r = current_beta - current_gamma;

          prob[i] = ((current_beta - current_gamma) / current_beta) / (1 + (((current_beta - current_gamma) / current_beta )/prob_I0 - 1) * exp(-(growth_r)*(t[i])));

        
    }else{
     

      
      //analytical beta and gamma, i(t) =((β - γ) / β ) / (1 + (((β - γ) / β )/i₀ - 1) * e^(-(β - γ)t))
      
      real current_beta = exp(log_beta_c[country[i]]+ cov_mat[i, ] * eta);
      //real current_beta = beta_c[country[i]];
      real current_gamma = gamma_c[country[i]];
      real growth_r = current_beta- current_gamma;// + cov_mat[i, ] * eta;


      if(t[i] - lastYear < pow(10,-4)){
        prob[i] = lastI;
      }else{
      //in terms of gamma, beta, and i0 ((β - γ) / β ) / (1 + (((β - γ) / β )/i₀ - 1) * e^(-(β - γ)t))
      prob[i] = ((current_beta - current_gamma) / current_beta) / (1 + (((current_beta - current_gamma) / current_beta )/lastI - 1) * exp(-(growth_r)*(t[i] - lastYear)));
      }
      

      
      
    }
    
    lastI = prob[i];
    lastCountry = country[i];
    lastYear = t[i];

    
  
 

 
 y[i] ~ binomial(cases[i], fmax( fmin(prob[i], 0.99999 ),0.00001));
 
  }
}
}

//test with temp eqn
generated quantities {
  vector[N] y_pred;
  //vector[C] r_c;
  vector[N] log_lik;
  vector[N] y_pred_counts;
  vector[N] y_pred_prop;
  
  
  int lastCountry = -1;
  real lastYear = -9999;
  real lastI = -9999;
  real K_c_country = -9999;
  for (i in 1:N) {
    
    if(country[i] != lastCountry){
       
         //calculate resistance rate at 2000 and use that to calculate the resistance rate at that first train year
          //real prob_I0 = inv_logit(I0_logit_c[country[i]]);

          real prob_I0 = I0_prob_c[country[i]];

          real current_beta =  exp(log_beta_c[country[i]] + cov_mat[i, ] * eta);
          //real current_beta =  beta_c[country[i]];// + cov_mat[i, ] * eta;
          real current_gamma = gamma_c[country[i]];
          real growth_r = current_beta - current_gamma;// + cov_mat[i, ] * eta;

          y_pred[i] = ((current_beta - current_gamma) / current_beta ) / (1 + (((current_beta - current_gamma) / current_beta )/prob_I0 - 1) * exp(-(growth_r)*(t[i])));

       
    }else{

      real current_beta =  exp(log_beta_c[country[i]]+ cov_mat[i, ] * eta);
      //test with analytical solution gamma and beta
      real current_gamma = gamma_c[country[i]] ;
      //real current_beta = beta_c[country[i]];//+ cov_mat[i, ] * eta;;// + cleanWeight*cleanHosp[i] + hygieneWeight*hygieneHosp[i]+ saniHospWeight*sanitationHosp[i] + waterWeight*waterHosp[i] + wasteWeight*wasteHosp[i] + tracssWeight*tracss[i] + tempWeight*temp[i] + precipWeight*precip[i] + handwashWeight*WHOhandwash[i] + abxConsumptionWeight*ddd_per_1000[i]; //could move the r_c term to other side
      real growth_r = current_beta- current_gamma;// + cov_mat[i, ] * eta;
      //if(t[i] - lastYear < pow(10,-4) || current_beta == current_gamma ){
      if(t[i] - lastYear < pow(10,-4)){
        y_pred[i] = lastI;
      }else{
      y_pred[i] = ((current_beta - current_gamma) / current_beta ) / (1 + (((current_beta - current_gamma) / current_beta )/lastI - 1) * exp(-(growth_r)*(t[i] - lastYear)));
        }

      
  
    } 
    
    y_pred_counts[i] = binomial_rng(cases[i],fmax(fmin(y_pred[i], 0.999999 ),0.000001));
    y_pred_prop[i] = y_pred_counts[i]/cases[i];
    
    lastI = y_pred[i];
    lastCountry = country[i];
    lastYear = t[i];
    
    log_lik[i] = binomial_lpmf(y[i]|cases[i], fmax(fmin(y_pred[i], 0.999999 ),0.000001));
  }
}

// Edit by ball.

// The input data is a vector 'y' of length 'N'.
// data {
 // int<lower=0> N;
//  vector[N] y;
//} 

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
//parameters {
//  real mu;
//  real<lower=0> sigma;
//}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
//model {
//  y ~ normal(mu, sigma);
//}

