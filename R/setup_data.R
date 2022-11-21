#' Format given data for model fitting
#'
#' For a given set of data, format it appropriately into a list with parameters expected from model fitting functions.
#' @param input.data A data.frame containing data. Expecting columns Response (observed values), Time, and Subject. If model.type is Continuum, also expecting Continuum column.
#' @param model.name Character string detailing model to be fit. Current options are Continuum and No Continuum.
#' @param priors A list of hyperparameter specifications needed for model fitting. Default option of NULL uses get.priors function to generate reasonable priors based on data. See details for more information.
#' @param seed Seed for random number generation.
#' @details
#' For model.type = "No Continuum", priors is expecting a list with elements:
#' betameans (4-element vector for population curve parameters),
#' umeans (4-element vector for individual effects),
#' Bmeans (two element vector for log variance model),
#' usigmean (real value for individual effect for log variance),
#' betavar (4x4 covariance matrix for population parameters),
#' uvar (4x4 covariance matrix for individual curve parameters),
#' Bvar (2x2 covariance matrix for parameters of log variance model),
#' usigsd (positive value for variance of individual effect in log variance model).
#'
#' For model.type = "Continuum", priors is expecting a list with elements:
#' betameans (4-element vector for population curve parameters),
#' cmeans (4-element vector for continuum effects),
#' umeans (4-element vector for individual effects),
#' Bmeans (two element vector for log variance model),
#' usigmean (real value for individual effect for log variance),
#' csigmean (real value for individual effect for log variance),
#' betavar (4x4 covariance matrix for population parameters),
#' cvar (4x4 covariance matrix for continuum parameters),
#' ucvar (4x4 covariance matrix for individual-continuum parameters),
#' uvar (4x4 covariance matrix for individual parameters),
#' Bvar (2x2 covariance matrix for parameters of log variance model),
#' usigsd (positive value for variance of individual effect in log variance model),
#' csigsd (positive value for variance of continuum effect in log variance model),
#' ucsigsd (positive value for variance of individual-continuum effect in log variance model)
#' @return An object of class VASCurves_Raw.
#' @examples
#' \dontrun{
#' ### No Continuum Curve Model
#' raw_data <- setup.data(input.data, "No Continuum")
#' ### Continuum Model
#' raw_data <- setup.data(input.data, "Continuum")
#' }
#' @import dplyr
#' @import saemix
#' @export
setup.data <- function(input.data, model.type, priors = NULL, seed = 123){

  if(model.type == "Continuum"){

    if(!all(c("Subject", "Continuum", "Time", "Response") %in% colnames(input.data))){
      stop("Input data doesn't have correct columns; expecting Subject, Time, Response")
    }

    warning("Data will be subset to complete cases!")

    input.data <- input.data %>%
      select(Subject, Continuum, Time, Response) %>%
      filter(complete.cases(.)) %>%
      mutate(SubjCont = paste0(Subject, "_", Continuum)) %>%
      arrange(SubjCont, Time)

    counts <- c(table(input.data$SubjCont))
    input.data.sub <- input.data %>%
      select(Subject, Continuum, SubjCont) %>%
      distinct()

    subject_continuum_map <- data.frame(Subject = input.data.sub$Subject,
                                        NewSubj = as.numeric(factor(input.data.sub$Subject,
                                                                    levels = unique(input.data.sub$Subject))),
                                        Continuum = input.data.sub$Continuum,
                                        NewCont = as.numeric(factor(input.data.sub$Continuum,
                                                                    levels = unique(input.data.sub$Continuum))),
                                        SubjCont = names(counts),
                                        NewSubjCont = 1:length(counts))
    names(counts) <- 1:length(counts)
    input.data$SubjCont <- as.numeric(unlist(lapply(names(counts),
                                                    function(x){rep(x,counts[x])})))

    muindex <- input.data %>%
      mutate(rowindexmin = 1:n(),
             rowindexmax = 1:n()) %>%
      group_by(SubjCont) %>%
      filter(row_number()==1 | row_number()==n()) %>%
      mutate(rowindexmin = min(rowindexmin),
             rowindexmax = max(rowindexmax)) %>%
      ungroup() %>%
      select(SubjCont, rowindexmin, rowindexmax) %>%
      distinct() %>%
      select(rowindexmin, rowindexmax) %>%
      data.frame(.)

    all_times <- unique(input.data$Time)
    time_positions <- unlist(lapply(input.data$Time, function(x){which(x == all_times)}))
    timeindex <- matrix(nrow = nrow(muindex),
                        ncol = max(table(input.data$SubjCont)))
    for(i in 1:nrow(timeindex)){
      timeindex[i,1:table(input.data$SubjCont)[i]] <- time_positions[muindex[i,1]:muindex[i,2]]
    }
    timeindex[is.na(timeindex)] <- max(time_positions) + 1

    timeindex2 <- matrix(nrow = nrow(muindex),
                         ncol = 2)
    timeindex2[,1] <- 1
    timeindex2[,2] <- table(input.data$SubjCont)

    stan.data.return <- list(N_subj_cont = nrow(subject_continuum_map),
                             N_subj = max(subject_continuum_map$NewSubj),
                             N_cont = max(subject_continuum_map$NewCont),
                             N_obs = nrow(input.data),
                             Y = input.data$Response,
                             time = all_times,
                             N_time = length(all_times),
                             timeindex = timeindex,
                             timeindex2 = timeindex2,
                             N_timecol = ncol(timeindex),
                             muindex = as.matrix(muindex),
                             subj_idx = subject_continuum_map$NewSubj,
                             cont_idx = subject_continuum_map$NewCont,
                             subjcont_idx = input.data$SubjCont,
                             Map = subject_continuum_map,
                             type = "Continuum")

    if(is.null(priors)){
      priors <- get.priors(stan.data.return, seed)
    }else{
      warning("Using user input priors. Please see details of setup.data to make sure it is the correct structure!")
    }

  }else if(model.type == "No Continuum"){

    if(!all(c("Subject", "Time", "Response") %in% colnames(input.data))){
      stop("Input data doesn't have correct columns; expecting Subject, Time, Response")
    }

    warning("Data will be subset to complete cases!")

    input.data <- input.data %>%
      select(Subject, Time, Response) %>%
      arrange(Subject, Time) %>%
      filter(complete.cases(.))

    counts <- c(table(input.data$Subject))
    subject_map <- data.frame(Original = names(counts),
                              New = 1:length(counts))
    names(counts) <- 1:length(counts)
    input.data$Subject <- as.numeric(unlist(lapply(names(counts),
                                                   function(x){rep(x,counts[x])})))

    muindex <- input.data %>%
      mutate(rowindexmin = 1:n(),
             rowindexmax = 1:n()) %>%
      group_by(Subject) %>%
      filter(row_number()==1 | row_number()==n()) %>%
      mutate(rowindexmin = min(rowindexmin),
             rowindexmax = max(rowindexmax)) %>%
      ungroup() %>%
      select(Subject, rowindexmin, rowindexmax) %>%
      distinct() %>%
      select(rowindexmin, rowindexmax) %>%
      data.frame(.)

    all_times <- unique(input.data$Time)
    time_positions <- unlist(lapply(input.data$Time, function(x){which(x == all_times)}))
    timeindex <- matrix(nrow = nrow(muindex),
                        ncol = max(table(input.data$Subject)))
    for(i in 1:nrow(timeindex)){
      timeindex[i,1:table(input.data$Subject)[i]] <- time_positions[muindex[i,1]:muindex[i,2]]
    }
    timeindex[is.na(timeindex)] <- max(time_positions) + 1

    timeindex2 <- matrix(nrow = nrow(muindex),
                         ncol = 2)
    timeindex2[,1] <- 1
    timeindex2[,2] <- table(input.data$Subject)

    stan.data.return <- list(N_subj = nrow(subject_map),
                             N_obs = nrow(input.data),
                             Y = input.data$Response,
                             time = all_times,
                             N_time = length(all_times),
                             timeindex = timeindex,
                             timeindex2 = timeindex2,
                             N_timecol = ncol(timeindex),
                             muindex = as.matrix(muindex),
                             subj_idx = input.data$Subject,
                             Map = subject_map,
                             type = "No Continuum")

    if(is.null(priors)){
      priors <- get.priors(stan.data.return, seed)
    }else{
      warning("Using user input priors. Please see details of setup.data to make sure it is the correct structure!")
    }

  }else{
    stop("model.type not recognized! Expecting Continuum or No Continuum")
  }

  ### Define class of return object
  stan.data.return <- c(stan.data.return, priors)
  class(stan.data.return) <- "VASCurves_Raw"

  return(stan.data.return)

}
