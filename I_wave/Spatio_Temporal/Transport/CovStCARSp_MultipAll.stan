
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
  
  real sparse_car_lpdf(vector phi, vector phiOld, real rho, real tau, real alpha, 
    int[,] W_sparse, vector W_weight, vector D_sparse, vector lambda, int n, int W_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_W; // phi' * W
      vector[n] ldet_terms;
      vector[n] phiNew = phi-rho*phiOld;
      
      phit_D = (phiNew .* D_sparse)';
      phit_W = rep_row_vector(0, n);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + W_weight[i]*phiNew[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + W_weight[i]*phiNew[W_sparse[i, 1]];
      }
    
      for (i in 1:n) ldet_terms[i] = log1m(alpha * lambda[i]);
      return 0.5 * (n * log(tau)
                    + sum(ldet_terms)
                    - tau * (phit_D * phiNew - alpha * (phit_W * phiNew)));
  }

  // real icar_normal_lpdf(vector phi, int N, int[] node1, int[] node2){
  //   return -0.5 * dot_self(phi[node1] - phi[node2]) + normal_lpdf(sum(phi) | 0, 0.001 * N);
  // }
  
}


data {
    // Data block
        int<lower=0> N; // Number of obs
        int Y[N]; // Vector of Obs

        int<lower=0> Nreg;  // Number of regions
        int<lower=0> Ntimes;  // Number of times
        real<lower=0> t[Ntimes]; // Times
        
        matrix<lower = 0>[Nreg, Nreg] W; // Adjacency
        int W_n;                         // Number of adjacent region pairs
        
        vector[N] lOff;   //Offset
        
        int<lower = 1> k; //Covariates
        matrix[N, k] X;
        
        int rId[N]; // Associate obs to region
        int tId[N]; // Associate obs to time

    }

transformed data {
        // Vector of zeros
        vector[Nreg] zeros;
        // Data for sparse car
        int W_sparse[W_n, 2];   // adjacency pairs
        vector[W_n] W_weight;     // Connection weights
        vector[Nreg] D_sparse;     // diagonal of D (number of neigbors for each site)
        vector[Nreg] lambda;       // eigenvalues of invsqrtD * W * invsqrtD
        
        zeros = rep_vector(0, Nreg);
  
        { // generate sparse representation for W
          int counter;
          counter = 1;
          // loop over upper triangular part of W to identify neighbor pairs
          for (i in 1:(Nreg - 1)) {
            for (j in (i + 1):Nreg) {
              if (W[i, j] > 0) {
                W_sparse[counter, 1] = i;
                W_sparse[counter, 2] = j;
                W_weight[counter] = W[i, j];
                counter = counter + 1;
              }
            }
          }
        }
        
        for (i in 1:Nreg) D_sparse[i] = sum(W[i]);
        
        {
          vector[Nreg] invsqrtD;  
          for (i in 1:Nreg) {
            invsqrtD[i] = 1 / sqrt(D_sparse[i]);
          }
          lambda = eigenvalues_sym(quad_form(W, diag_matrix(invsqrtD)));
        }
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
        
        // Coefficients
        vector[k] beta;
        // CAR parameters
        vector[Nreg] phi[Ntimes];
        real<lower=0> tau;
        real<lower = 0, upper = 1> alpha;
        
        // AR parameters
        real<lower=-1, upper=1> rho;

    }

transformed parameters {
        // Transformed parameters block
        real<lower=0> sigma;
        vector[N] lm;
        
        // SD from Precision
        sigma = 1/tau;

        // Mean function
        for (i in 1:N) {
            lm[i] = phi[tId[i]][rId[i]] + 
                    log(exp(logbase[rId[i]]) +
                    exp(log_RichLin(logr[rId[i]], logh[rId[i]], p[rId[i]], logs[rId[i]], t[tId[i]])));
        }
        lm = lOff + X*beta + lm;
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
        
        // Coefficient prior
        beta ~ normal(0, 10);
        // CAR prior
        phi[1] ~ sparse_car(zeros, rho, tau, alpha, 
                            W_sparse, W_weight, D_sparse, lambda, Nreg, W_n);
        for (i in 2:Ntimes){
              phi[i] ~ sparse_car(phi[i-1], rho, tau, alpha, 
                                  W_sparse, W_weight, D_sparse, lambda, Nreg, W_n);
        }
        alpha ~ beta(.5, .5);
        tau ~ gamma(2, 2);

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

