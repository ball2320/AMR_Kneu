
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


}

parameters {

  
  //OPTION 1: Initial Resistance - hierarchical in logit space (works better if we start
  //at the first point in our data and works for most part for 2000 start)
  // real I0_logit;
  // vector[C] I0_logit_c; // no bounds needed if on logit scale
  // real<lower=0> sigma_I0_logit;
  // //real I0_logit_tilde[C];
  
   //OPTION 2: Initial Resistance - non hierarchical contrained btw 0 and 1 (works better 
// if we start at previous year with no data ex 2000 start)
  vector<lower=0, upper=1>[C] I0_prob_c;
  
  
  //reparam with beta gamma I0 
  real<lower=0> beta;
  vector<lower=0>[C] beta_c;
  real<lower=0> sigma_beta;
  //real beta_tilde[C];

  real<lower=0> gamma;
  vector<lower=0>[C] gamma_c;
  real<lower=0> sigma_gamma;
  //real gamma_tilde[C];

  
}


model {

  //OPTION 1: Initial Resistance 
  // I0_logit ~ normal(0,1);
  // sigma_I0_logit ~ exponential(2); //0.5 //3
  // I0_logit_c ~ normal(I0_logit,sigma_I0_logit);
  // //I0_logit_tilde ~ std_normal();
  
  //OPTION 2: Initial Resistance
  I0_prob_c ~ beta(2,2);


  //reparam with beta gamma and I0
  beta ~ exponential(3);
  sigma_beta ~ exponential(0.5);
  beta_c ~ normal(beta,sigma_beta);
  //beta_tilde ~ normal(0,1);

  gamma ~ exponential(3);
  sigma_gamma ~ exponential(0.5);
  gamma_c ~ normal(gamma, sigma_gamma);


if (calc_likelihood == 1){
  
  vector[N] prob;
  
  int lastCountry = -1;
  real lastYear = -9999;
  real lastI = -9999;
  
  // Likelihood
  for(i in 1:N){
   
    //can change this so if you have unordered data it should still work or signpost better that the data has to be ordered!!!!
    if(country[i] != lastCountry){
        
        
          //resistance rate at 2000 or whatever start year and use that to calculate the resistance rate at that first train year
          //real prob_I0 = inv_logit(I0_logit_c[country[i]]); //use this for option 1 initial resistance
          real prob_I0 = I0_prob_c[country[i]];  //use this for option 2 initial resistance 

          real current_beta =  beta_c[country[i]];
          real current_gamma = gamma_c[country[i]];

          prob[i] = ((current_beta - current_gamma) / current_beta ) / (1 + (((current_beta - current_gamma) / current_beta )/prob_I0 - 1) * exp(-(current_beta - current_gamma)*(t[i])));

          
          //now just fitting with no backcast so now just estimate the initial resistance value
          //prob[i] = inv_logit(I0_logit_c[country[i]]);
        
        
    }else{
     

      
      real current_beta = beta_c[country[i]];
      real current_gamma = gamma_c[country[i]];

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
      
       
         // //calculate resistance rate at 2000 and use that to calculate the resistance rate at that first train year
          //real prob_I0 = inv_logit(I0_logit_c[country[i]]); 
          real prob_I0 = I0_prob_c[country[i]]; 
 
          real current_beta =  beta_c[country[i]];
          real current_gamma = gamma_c[country[i]];

          y_pred[i] = ((current_beta - current_gamma) / current_beta ) / (1 + (((current_beta - current_gamma) / current_beta )/prob_I0 - 1) * exp(-(current_beta - current_gamma)*(t[i])));

         
       
    }else{

      
      
      real current_gamma = gamma_c[country[i]];
      real current_beta = beta_c[country[i]];
     
      if(t[i] - lastYear < pow(10,-4)){
        y_pred[i] = lastI;
      }else{
      y_pred[i] = ((current_beta - current_gamma) / current_beta ) / (1 + (((current_beta - current_gamma) / current_beta )/lastI - 1) * exp(-(current_beta - current_gamma)*(t[i] - lastYear)));
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

