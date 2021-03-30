
functions {
  real log_RichLin(real lr, real lh, real pe, real ls, real ti) {
    real h = exp(lh);
    real s = exp(ls);
    real lout;
    
    real lhpt = log(1+exp(h * (pe-ti)));
    if (is_inf(lhpt))
    {
        lhpt = h * (pe-ti);
    }
    
    lout = lr + ls + lh + h*(pe-ti) - (s+1) * lhpt;

    return lout;
  }
}


data {
        // Data block
        int<lower=0> N; // Number of obs
        int Y[N]; // Vector of Obs

        int<lower=0> Nreg;  // Number of regions
        int<lower=0> Ntimes;  // Number of times
        real<lower=0> t[Ntimes]; // Times
        
        vector[N] lOff;   //Offset
        
        int rId[N]; // Associate obs to region
        int tId[N]; // Associate obs to time

    }
      
parameters {  
    // Parameters block
    
        // Baseline parameter
        real logbase[Nreg];
        // Richard's parameter
        real logr[Nreg];
        real logh[Nreg];
        real p[Nreg];
        real logs[Nreg];
        
        vector[Nreg] theta[Ntimes];
        real<lower=0> taut;
        real<lower=-1, upper=1> rho;
    }

transformed parameters {
        // Transformed parameters block
        vector[N] lm;
        real<lower=0> sigmat;
        
        // SD from Precision
        sigmat = 1/taut;
        
        // Mean function
        for (i in 1:N) {
            lm[i] = theta[tId[i]][rId[i]]*sigmat + 
                    log(exp(logbase[rId[i]]) +
                    exp(log_RichLin(logr[rId[i]], logh[rId[i]], p[rId[i]], logs[rId[i]], t[tId[i]])));
        }
        lm = lOff + lm;
    } 


    model {
        // Model block
        // Baseline priors
        logbase ~ normal(0, 1);
        // Richards Priors
        logr ~ normal(0, 1);
        logh ~ normal(0, 1);
        p ~ normal(t[Ntimes]/2, t[Ntimes]/(2*1.96));
        logs ~ normal(0, 1);
        
        // Random effects
        theta[1] ~ normal(0, 1);
        for (i in 2:Ntimes){
          theta[i] ~ normal(rho*theta[i-1], 1);
        }
        taut ~ gamma(2, 2);

        // likelihood
        Y ~ poisson_log(lm);
    }

    generated quantities{
        // Output block
        real Y_pred[N];
        vector[N] log_lik;
        
        Y_pred = poisson_log_rng(lm);  // Posterior predictive distribution
        for (i in 1:N){
              log_lik[i] = poisson_log_lpmf(Y[i] | lm[i]);
        }
    }

