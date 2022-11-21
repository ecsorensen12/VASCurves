#' Fitting curve model and drawing samples
#' 
#' Fit the appropriate VAS curve model to the formatted data.
#' @param stan.data An object of class VASCurves_Raw. See setup.data function.
#' @param n.chain The number of Markov chains to be run.
#' @param iter The number of iterations for sampling.
#' @param warmup The number of warm-up iterations to perform.
#' @param cores The number of cores available for running model.
#' @param initial.values Optional list of starting values for sampler.
#' @param gen.inits Logical value for whether initial values should be generated automatically using input VASCurves_Raw object. Only applies if initial.values is also NULL. See details.
#' @param control.param List of control settings for performing sampling.
#' @param verbose Would you like printed messages from sampler?
#' @param seed Seed for random number generation.
#' @param ... Other parameters passed to rstan::sampling function.
#' @details 
#' gen.inits generates initial values for each chain by initializing all individual and continuum effects to zero 
#' while simultaneously using a random multiplier between 0.9 and 1.1 applied to the prior specification for betameans and Bmeans.
#' Diagnostics in the VASCurves_Res object are calculated using the gelman.diag function in the coda package if n.chain > 1. Otherwise, the Rhat statistic from the rstan package is given for each parameter.
#' @return An object of class VASCurves_Res that includes the resultant stan object from the model fit, number of chains, VASCurves_Raw input object, diagnostics of the parameters traced, and the summary of the model fit. See details.
#' @import rstan
#' @import coda
#' @export
fit.model <- function(stan.data,
                      n.chain = 1,
                      iter = 20000,
                      warmup = 10000,
                      cores = 1,
                      initial.values = NULL,
                      gen.inits = T,
                      control.param = list(adapt_delta = 0.995,
                                           max_treedepth = 12,
                                           stepsize = 0.01),
                      verbose = T,
                      seed = NULL,
                      ...){
  
  if(class(stan.data) != "VASCurves_Raw"){
    stop("Class not recognized! Expecting VASCurves_Raw object returned from setup.data.")
  }
  
  ### Check Initial Values
  
  if(is.null(initial.values) & !gen.inits){
    warning("No initial values are provided and gen.inits is FALSE. Initial values will be chosen by STAN.")
  }else if(!is.null(initial.values)){
    if(!is.list(initial.values)){
      stop("Initial values not list!")
    }else if(n.chain != length(initial.values)){
      stop("Initial values not specified for all chains!")
    }else if(T %in% unlist(lapply(initial.values,function(x){!is.list(x)}))){
      stop("Initial values for at least one chain not a list!")
    }
  }
  
  ### Generate Initial Values
  if(is.null(initial.values) & gen.inits){
    if(stan.data$type == "No Continuum"){
      multipliers <- c(1, runif(n.chain-1, 0.9, 1.1))
      initial.values <- list()
      for(i in 1:n.chain){
        initial.values[[i]] <- list(usig_i = rep(0, stan.data$N_subj),
                                    B = multipliers[i]*stan.data$Bmeans,
                                    beta = multipliers[i]*stan.data$betameans,
                                    U = lapply(1:stan.data$N_subj, 
                                               function(x){c(0,0,0,0)}))
      }
    }else if(stan.data$type == "Continuum"){
      multipliers <- c(1, runif(n.chain-1, 0.9, 1.1))
      initial.values <- list()
      for(i in 1:n.chain){
        initial.values[[i]] <- list(beta = stan.data$betameans*multipliers[i],
                                    B = stan.data$Bmeans*multipliers[i],
                                    C = lapply(1:stan.data$N_cont, function(x){c(0,0,0,0)}),
                                    U = lapply(1:stan.data$N_subj, function(x){c(0,0,0,0)}),
                                    UC = lapply(1:stan.data$N_subj_cont, function(x){c(0,0,0,0)}),
                                    Csig = rep(0,stan.data$N_cont),
                                    Usig = rep(0,stan.data$N_subj),
                                    UCsig = rep(0,stan.data$N_subj_cont))
      }
    }else{
      warning("Model type not recognized! Expecting Continuum or No Continuum.")
    }
  }
  
  ### Writing model and setting parameters to trace
  
  smod <- rstan::stan_model(model_code = choose.model(stan.data$type), 
                            auto_write = TRUE)
  
  if(stan.data$type == "No Continuum"){
    pars_follow <- c("beta", "B", "U", "usig_i")
  }else if(stan.data$type == "Continuum"){
    pars_follow <- c("beta", "B", "C", "U", "UC", "Csig", "Usig", "UCsig")
  }else{
    stop("Specified model name not recognized!")
  }
  
  ### Fitting Model
  
  result <- rstan::sampling(smod, 
                            data = stan.data,
                            chain = n.chain,
                            iter = iter,
                            warmup = warmup,
                            cores = cores,
                            init = initial.values,
                            control = control.param,
                            pars = pars_follow, 
                            verbose = verbose,
                            ...)
  result.summary <- summary(result)$summary
  
  ### Diagnostics
  
  if(n.chain > 1){
    param_names <- names(result)
    param_names <- param_names[which(param_names != "lp__")] 
    samples <- extract(result, permuted = F)
    results <- matrix(nrow = length(param_names),
                      ncol = 2)
    colnames(results) <- c("Point Est.", "Upper C.I.")
    for(i in 1:length(param_names)){
      temp_sample <- lapply(as.list(1:n.chain),
                            function(x){
                              as.mcmc(samples[,x,param_names[i]])
                            })
      param_mcmc <- as.mcmc.list(temp_sample)
      geldiag <- gelman.diag(param_mcmc)
      results[i,] <- geldiag$psrf
    }
    rownames(results) <- param_names
  }else{
    results <- result.summary[,"Rhat"]
  }
  
  ### Return Object
  
  res_VASCurves <- list(data = stan.data,
                        result = result,
                        n.chain = n.chain,
                        diagnostics = results,
                        result.summary = result.summary)
  class(res_VASCurves) <- "VASCurves_Res"
  
  return(res_VASCurves)
  
}

choose.model <- function(type){
  
  if(type == "No Continuum"){
    
    model.return <-  "
    data{
    
      int<lower=0> N_subj;
      int<lower=0> N_obs;
      int<lower=0> N_timecol;
      int<lower=0> N_time;
      real Y[N_obs];
      vector[N_time] time;
      vector[4] betameans;
      vector[4] umeans;
      vector[2] Bmeans;
      real usigmean;
      matrix[4,4] betavar;
      matrix[4,4] uvar;
      matrix[2,2] Bvar;
      real usigsd;
      int<lower=0> timeindex[N_subj,N_timecol];
      int<lower=0> timeindex2[N_subj,2];
      int<lower=0> muindex[N_subj,2];
      int<lower=0> subj_idx[N_obs];
    }
    parameters{
    
      real usig_i[N_subj];
      vector[2] B;
      vector[4] beta;
      vector[4] U[N_subj];
      
    }
    transformed parameters{
    
      vector[N_obs] mu;
      vector[N_time] temp;
      vector<lower=0>[N_subj] sigma;
      vector[N_subj] sigma_sqrt;
      vector[4] beta_i[N_subj];
      
      for(i in 1:N_subj){
        beta_i[i] = beta + U[i];     
        temp = rep_vector((beta_i[i,2] - beta_i[i,1]), N_time);
        temp ./= (1+exp((4*beta_i[i,3]/(beta_i[i,2] - beta_i[i,1]))*(beta_i[i,4] - time)));
        temp += beta_i[i,1];
        mu[muindex[i,1]:muindex[i,2]] = temp[timeindex[i,timeindex2[i,1]:timeindex2[i,2]]];
        sigma[i] = exp(B[1] + (B[2]*beta_i[i,3]) + usig_i[i]);
      }
      
      sigma_sqrt = sqrt(sigma);
    
    }
    model{
    
      
      beta ~ multi_normal(betameans, betavar);
      B ~ multi_normal(Bmeans, Bvar);
      Y ~ normal(mu, sigma_sqrt[subj_idx]);
      U ~ multi_normal(umeans, uvar);
      usig_i ~ normal(usigmean, usigsd);
      
    }
    "
    
  }else if(type == "Continuum"){
    
    model.return <- "
    data{
    
      int<lower=0> N_subj_cont;
      int<lower=0> N_cont;
      int<lower=0> N_subj;
      int<lower=0> N_obs;
      int<lower=0> N_timecol;
      int<lower=0> N_time;
      real Y[N_obs];
      vector[N_time] time;
      vector[4] betameans;
      vector[4] cmeans;
      vector[4] umeans;
      vector[2] Bmeans;
      real usigmean;
      real csigmean;
      matrix[4,4] betavar;
      matrix[4,4] cvar;
      matrix[4,4] ucvar;
      matrix[4,4] uvar;
      matrix[2,2] Bvar;
      real usigsd;
      real csigsd;
      real ucsigsd;
      int<lower=0> timeindex[N_subj_cont,N_timecol];
      int<lower=0> timeindex2[N_subj_cont,2];
      int<lower=0> muindex[N_subj_cont,2];
      int<lower=0> subjcont_idx[N_obs];
      int<lower=0> subj_idx[N_subj_cont];
      int<lower=0> cont_idx[N_subj_cont];
    }
    parameters{
    
      vector[2] B;
      vector[4] beta;
      vector[4] UC[N_subj_cont];
      vector[4] C[N_cont];
      vector[4] U[N_subj];
      real Csig[N_cont];
      real UCsig[N_subj_cont];
      real Usig[N_subj];
      
    }
    transformed parameters{
    
      vector[N_obs] mu;
      vector[N_time] temp;
      vector<lower=0>[N_subj_cont] sigma;
      vector[N_subj_cont] sigma_sqrt;
      vector[4] beta_li[N_subj_cont];
      
      for(i in 1:N_subj_cont){
        beta_li[i] = beta + C[cont_idx[i]] + UC[i];     
        temp = rep_vector((beta_li[i,2] - beta_li[i,1]), N_time);
        temp ./= (1+exp((4*beta_li[i,3]/(beta_li[i,2] - beta_li[i,1]))*(beta_li[i,4] - time)));
        temp += beta_li[i,1];
        mu[muindex[i,1]:muindex[i,2]] = temp[timeindex[i,timeindex2[i,1]:timeindex2[i,2]]];
        sigma[i] = exp(B[1] + (B[2]*beta_li[i,3]) + Csig[cont_idx[i]] + UCsig[i]);
      }
      
      sigma_sqrt = sqrt(sigma);
    
    }
    model{
      
      beta ~ multi_normal(betameans, betavar);
      C ~ multi_normal(cmeans, cvar);
      UC ~ multi_normal(U[subj_idx], ucvar);
      U ~ multi_normal(umeans, uvar);
      B ~ multi_normal(Bmeans, Bvar);
      Csig ~ normal(csigmean, csigsd);
      UCsig ~ normal(Usig[subj_idx], ucsigsd);
      Usig ~ normal(usigmean, usigsd);
      Y ~ normal(mu, sigma_sqrt[subjcont_idx]);
      
    }
    "
    
  }else{
    
    stop("Specified model type is not recognized!")
    
  }
  
  return(model.return)
}