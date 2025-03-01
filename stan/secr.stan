data {
    int<lower = 1> N;
    int<lower = 1> n_trap;
    int<lower = 0> yred[N,n_trap]; // detection sum
    int<lower = 0, upper = 1> everdetcted[N]; // ever detected
    vector<lower=0, upper=1>[N] sex; // sex of the individuals
    int<lower = 0> deployred[n_trap]; //total number of deployments
    matrix[n_trap, 2] X; //trap locations

    int<lower = 1> n_grid;
    matrix[n_grid, 2] grid_pts;
    
    matrix[n_grid, n_trap] distsqr; // distance between grid points and camera traps

    // get enviromental variables at grid points
    int<lower = 1> n_env; // without intercept!
    matrix[n_grid, n_env] envX;
}

parameters {
    // whether alive or not 
    vector[2] log_psi_mean; // distance decay of survival probability
    vector[2] log_psi_var; // sd of beta
    vector[N] log_psi; // survival probability

    
    
    // detection related 
    vector[2] log_sigma_mean; // distance decay of detection probability
    vector[2] log_sigma_var; // sd of beta
    vector[N] log_sigma; // sd of detection probability


    vector[2] log_p0_mean; // detection probability at trap
    vector[2] log_p0_var; // sd of beta
    vector[N] log_p0; // detection probability at trap

    // activitiy center realted
    vector[n_env] betaenv;
}

model {
    {
        vector[n_grid] log_softmax_Xbeta;
        vector[n_grid] logliklocal_occu; // log likelihood of this individual, condition on AC
        vector[n_grid] loglikdetprob; // log likelihood of this individual, condition on AC
        vector[n_grid] log_det; // log detection probability
        real loglik_ind_tmp;

        loglikdetprob = rep_vector(0, n_grid);
        logliklocal_occu = rep_vector(0, n_grid);

        log_softmax_Xbeta = log_softmax(envX * betaenv); // prior of activity center
        for(i in 1:N){
            for(j in 1:n_trap){
                log_det = log_p0[i] - exp(log_sigma[i]) * distsqr[,j];
                loglikdetprob += yred[i,j] * log_det + (deployred[j] - yred[i,j]) * log1m_exp(log_det);
            }

            logliklocal_occu = log_softmax_Xbeta + loglikdetprob + log_psi[i];
            loglik_ind_tmp = log_sum_exp(logliklocal_occu);
            if (everdetcted[i] == 1){
                target += loglik_ind_tmp;
            }
            else{
                target += log_sum_exp(loglik_ind_tmp, log1m_exp(log_psi[i]));
            }

            log_p0[i] ~ normal(log_p0_mean[1] * (1-sex[i]) + log_p0_mean[2] * sex[i],
                                log_p0_var[1] * (1-sex[i]) + log_p0_var[2] * sex[i]); 
            log_psi[i] ~ normal(log_psi_mean[1] * (1-sex[i]) + log_psi_mean[2] * sex[i],
                                log_psi_var[1] * (1-sex[i]) + log_psi_var[2] * sex[i]);

            log_sigma[i] ~ normal(log_sigma_mean[1] * (1-sex[i]) + log_sigma_mean[2] * sex[i],
                                log_sigma_var[1] * (1-sex[i]) + log_sigma_var[2] * sex[i]);
            
        }
        // priors
        log_psi_mean ~ normal(0, 1);
        log_psi_var ~ normal(0, .1);
        log_sigma_mean ~ normal(0, 1);
        log_sigma_var ~ normal(0, .1);
        log_p0_mean ~ normal(0, 1);
        log_p0_var ~ normal(0, .1);
        betaenv ~ normal(0, 100);
    }
}

generated quantities {
   vector[N] z; // alive or not
   vector[N] s; // activity center
   {
        vector[n_grid] log_softmax_Xbeta;
        vector[n_grid] logliklocal_occu; // log likelihood of this individual, condition on AC
        vector[n_grid] loglikdetprob; // log likelihood of this individual, condition on AC
        vector[n_grid] log_det; // log detection probability

        real loglik_ind_tmp;

        loglikdetprob = rep_vector(0, n_grid);
        logliklocal_occu = rep_vector(0, n_grid);

        log_softmax_Xbeta = log_softmax(envX * betaenv); // prior of activity center
        for (i in 1:N){
            for (j in 1:n_trap){
                log_det = log_p0[i] - exp(log_sigma[i]) * distsqr[,j];
                loglikdetprob += yred[i,j] * log_det + (deployred[j] - yred[i,j]) * log1m_exp(log_det);
            }

            logliklocal_occu = log_softmax_Xbeta + loglikdetprob + log_psi[i];
            loglik_ind_tmp = log_sum_exp(logliklocal_occu);
            if (everdetcted[i] == 1){
                z[i] = 1;
                s[i] = categorical_rng(exp(logliklocal_occu));
            }
            else{
                z[i] = bernoulli_rng(exp(log_psi[i]));
                s[i] = categorical_rng(exp(logliklocal_occu));
            }
            
        }
    }
}

