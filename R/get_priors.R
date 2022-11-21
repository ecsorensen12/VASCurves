#' Generate prior values for model
#' 
#' Helper function used for calculating priors in setup.data function.
#' @param stan.data A list of data created in the process of running setup.data.
#' @param seed Seed for random number generation.
#' @return A named list of hyperparameter values for model.
#' @import saemix
#' @import dplyr
#' @export
get.priors <- function(stan.data, seed = 123){
  
  vas.f <- function(min,max,s,x0,time){(max-min)/(1+exp((4*s/(max-min))*(x0-time)))+min}
  
  min_pop <- c()
  max_pop <- c()
  slope_pop <- c()
  x0_pop <- c()
  sigma_overall <- c()
  
  if(stan.data$type == "No Continuum"){
    subj_idx <- stan.data$subj_idx
  }else{
    subj_idx <- stan.data$subjcont_idx
  }
  
  time_full <- NULL
  for(i in 1:nrow(stan.data$timeindex)){
    time_full <- c(time_full,
                   stan.data$time[stan.data$timeindex[i,stan.data$timeindex2[i,1]:stan.data$timeindex2[i,2]]])
  }
  
  ### Subset data to persons in group
  
  temp_data <- data.frame(Y = stan.data$Y,
                          time = time_full,
                          subj_idx = subj_idx)
  
  temp_data <- temp_data %>%
    group_by(subj_idx) %>%
    mutate(id = dplyr::cur_group_id())
  
  ### Define saemix data
  
  saemix_temp_data <- saemix::saemixData(name.data = temp_data,
                                         name.group = "id",
                                         name.predictors = "time",
                                         name.response = "Y",
                                         verbose = F)
  
  ### Define model
  
  model <- function(psi,id,x) {
    t <- x[,1]
    min <- psi[id,1]
    max <- psi[id,2]
    slope <- psi[id,3]
    x0 <- psi[id,4]
    fpred <- vas.f(min, max, slope, x0, t)
    return(fpred)
  }
  
  ### Define saemix model
  
  saemix_temp_model <- saemix::saemixModel(model = model,
                                           psi0 = c(min=min(temp_data$Y),
                                                    max=max(temp_data$Y),
                                                    slope=(max(temp_data$Y) - min(temp_data$Y))/
                                                      (max(temp_data$time) - min(temp_data$time)),
                                                    x0=median(temp_data$time)),
                                           verbose = F)
  
  ### Define options to remove printing
  
  saemix_temp_options <- list(map=TRUE,
                              fim=F,
                              ll.is=FALSE,
                              displayProgress=F,
                              seed=seed,
                              print = F,
                              save = F,
                              save.graphs = F)
  
  ### Fit model
  
  saemix_temp_fit <- saemix::saemix(saemix_temp_model,
                                    saemix_temp_data,
                                    saemix_temp_options)
  
  ### Get empirical priors
  
  min_pop <- saemix_temp_fit@results@fixed.effects[1]
  max_pop <- saemix_temp_fit@results@fixed.effects[2]
  slope_pop <- saemix_temp_fit@results@fixed.effects[3]
  x0_pop <- saemix_temp_fit@results@fixed.effects[4]
  B0 <- log(saemix_temp_fit@results@respar[1])
  Var <- as.matrix(saemix_temp_fit@results@omega)
  
  if(stan.data$type == "No Continuum"){
    priors <- list(betameans = c(min_pop, max_pop, slope_pop, x0_pop),
                   umeans = c(0,0,0,0),
                   Bmeans = c(B0,0),
                   usigmean = 0,
                   betavar = Var,
                   uvar = Var,
                   Bvar = diag(c(B0,0.5)),
                   usigsd = 1)
  }else{
    priors <- list(betameans = c(min_pop, max_pop, slope_pop, x0_pop),
                   cmeans = c(0,0,0,0),
                   umeans = c(0,0,0,0),
                   Bmeans = c(B0,0),
                   usigmean = 0,
                   csigmean = 0,
                   betavar = Var,
                   cvar = Var,
                   ucvar = Var,
                   uvar = Var,
                   Bvar = diag(c(B0,0.5)),
                   usigsd = 0.5,
                   csigsd = 0.5,
                   ucsigsd = 0.5)
  }
  
  return(priors)
  
}