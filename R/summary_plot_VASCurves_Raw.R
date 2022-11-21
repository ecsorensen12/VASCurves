VASCurves_Raw <- function(object, ...){
  UseMethod("VASCurves_Raw")
}

#' Method for summary function in VASCurves
#' 
#' Summary method for an object of class VASCurves_Raw.
#' @param object An object of class VASCurves_Raw.
#' @param digits Number of digits to print for numeric variables.
#' @param ... Other parameters.
#' @return An object of summary statistics of the data.
#' @method summary VASCurves_Raw
#' @export
summary.VASCurves_Raw <- function(object, digits = 3, ...){
  
  if(object$type == "No Continuum"){
    return_list <- list(object$type,
                        object$N_obs,
                        object$N_subj,
                        length(object$time),
                        c("mean" = round(mean(object$time), digits),
                          "sd" = round(sd(object$time), digits),
                          "min" = round(min(object$time), digits),
                          "max" = round(max(object$time), digits)),
                        c("mean" = round(mean(object$Y), digits),
                          "sd" = round(sd(object$Y), digits),
                          "min" = round(min(object$Y), digits),
                          "max" = round(max(object$Y), digits)))
    names(return_list) <- c("Type", "N_obs", "N_subj", "N_time", "Time_Stats", "Overall_Stats")
  }else{
    return_list <- list(object$type,
                        object$N_obs,
                        object$N_subj,
                        object$N_cont,
                        length(object$time),
                        c("mean" = round(mean(object$time), digits),
                          "sd" = round(sd(object$time), digits),
                          "min" = round(min(object$time), digits),
                          "max" = round(max(object$time), digits)),
                        c("mean" = round(mean(object$Y), digits),
                          "sd" = round(sd(object$Y), digits),
                          "min" = round(min(object$Y), digits),
                          "max" = round(max(object$Y), digits)))
    names(return_list) <- c("Type", "N_obs", "N_subj", "N_cont", "N_time", "Time_Stats", "Overall_Stats")
  }
  
  class(return_list) <- "summary.VASCurves_Raw"
  
  return(return_list)
}

#' Print method for summary function in VASCurves
#' 
#' Method to print summary object returned from summary.VASCurves_Raw.
#' @param object An object of class VASCurves_Raw.
#' @param ... Other parameters.
#' @return An object of summary statistics of the data.
#' @method print summary.VASCurves_Raw
#' @export  
print.summary.VASCurves_Raw <- function(object, ...){
  
  if(object$Type == "No Continuum"){
    
    obs_num_string <- paste0("Number of Observations: ", object$N_obs)
    sub_num_string <- paste0("Number of Subjects: ", object$N_subj)
    time_num_string <- paste0("Number of Times: ", object$N_time)
    time_string_overall <- paste0("Average time: ", object$Time_Stats[1], 
                                  " (SD: ",object$Time_Stats[2], 
                                  ", Min: ", object$Time_Stats[3],
                                  ", Max: ", object$Time_Stats[4],
                                  ")")
    response_string_overall <- paste0("Overall mean response: ", object$Overall_Stats[1], 
                                      " (SD: ",object$Overall_Stats[2], 
                                      ", Min: ", object$Overall_Stats[3],
                                      ", Max: ", object$Overall_Stats[4],
                                      ")")
    
    cat(obs_num_string,
        "\n",
        sub_num_string,
        "\n",
        time_num_string,
        "\n",
        time_string_overall,
        "\n",
        response_string_overall,
        sep = "")
    
  }else{
    
    obs_num_string <- paste0("Number of Observations: ", object$N_obs)
    sub_num_string <- paste0("Number of Subjects: ", object$N_subj)
    cont_num_string <- paste0("Number of Continua: ", object$N_cont)
    time_num_string <- paste0("Number of Times: ", object$N_time)
    time_string_overall <- paste0("Average time: ", object$Time_Stats[1], 
                                  " (SD: ",object$Time_Stats[2], 
                                  ", Min: ", object$Time_Stats[3],
                                  ", Max: ", object$Time_Stats[4],
                                  ")")
    response_string_overall <- paste0("Overall mean response: ", object$Overall_Stats[1], 
                                      " (SD: ",object$Overall_Stats[2], 
                                      ", Min: ", object$Overall_Stats[3],
                                      ", Max: ", object$Overall_Stats[4],
                                      ")")
    
    cat(obs_num_string,
        "\n",
        sub_num_string,
        "\n",
        cont_num_string,
        "\n",
        time_num_string,
        "\n",
        time_string_overall,
        "\n",
        response_string_overall,
        sep = "")
  }
}

#' Method for plot function in VASCurves
#' 
#' Plotting method for an object of class VASCurves_Raw. Plots the raw data for each individual in the dataset.
#' @param object An object of class VASCurves_Raw.
#' @param subj_ids Vector of string values detailing which subject ids to plot respective data. Defaults to "all" which gives a plot for all subjects.
#' @param ... Other graphical parameters.
#' @return Plots of the raw data for specified individuals in the input VASCurves_Raw object.
#' @method plot VASCurves_Raw
#' @import ggplot2
#' @export
plot.VASCurves_Raw <- function(object, subj_ids = "all", ...){
  
  if(object$type == "No Continuum"){
    
    if("all" %in% subj_ids){
      subj_idx <- 1:object$N_subj
    }else{
      subj_ids <- as.character(subj_ids)
      if(all(subj_ids %in% object$Map$Original)){
        subj_idx <- object$Map$New[object$Map$Original %in% subj_ids]
      }else{
        stop("Input IDs not found in VASCurves_Raw input object!")
      }
    }
    
    ### Print plot of raw data for subjects
    for(i in subj_idx){
      
      temp_data <- list()
      
      temp_data$Y <- object$Y[object$subj_idx == i]
      temp_data$time <- object$time[object$timeindex[i,object$timeindex2[i,1]:object$timeindex2[i,2]]]
      
      print(ggplot(data = data.frame(Value = temp_data$Y,
                                     Time = temp_data$time),
                   aes(x = Time, y = Value)) +
              geom_point() +
              xlim(min(object$time, na.rm = T), 
                   max(object$time, na.rm = T)) +
              ylim(min(object$Y, na.rm = T), 
                   max(object$Y, na.rm = T)) +
              ggtitle(paste0(object$Map$Original[i],
                             " Raw Data")))
      
    }
    
  }else{
    
    if("all" %in% subj_ids){
      subj_idx <- 1:object$N_subj
    }else{
      subj_ids <- as.character(subj_ids)
      if(all(subj_ids %in% object$Map$Subject)){
        subj_idx <- unique(object$Map$NewSubj[object$Map$Subject %in% subj_ids])
      }else{
        stop("Input IDs not found in VASCurves_Raw input object!")
      }
    }
    
    ### Print plot of raw data for subjects across continua
    for(i in subj_idx){
      
      temp_data <- list()
      
      subjcont_idxs <- object$Map$NewSubjCont[object$Map$NewSubj == i]
      continuums <- object$Map$Continuum[object$Map$NewSubj == i]
      temp_data$Y <- object$Y[object$subjcont_idx %in% subjcont_idxs]
      
      
      times <- c()
      continuums_full <- c()
      for(j in 1:length(subjcont_idxs)){
        times <- c(times,
                   object$time[object$timeindex[subjcont_idxs[j],
                                                object$timeindex2[subjcont_idxs[j],1]:object$timeindex2[subjcont_idxs[j],2]]])
        continuums_full <- c(continuums_full,
                             rep(continuums[j], 
                                 length(object$timeindex2[subjcont_idxs[j],1]:object$timeindex2[subjcont_idxs[j],2])))
      }
      temp_data$time <- times
      temp_data$continuum <- continuums_full
      
      print(ggplot(data = data.frame(Value = temp_data$Y,
                                     Time = temp_data$time,
                                     Continuum = temp_data$continuum),
                   aes(x = Time, y = Value, col = Continuum)) +
              geom_point() +
              xlim(min(object$time, na.rm = T), 
                   max(object$time, na.rm = T)) +
              ylim(min(object$Y, na.rm = T), 
                   max(object$Y, na.rm = T)) +
              ggtitle(paste0(unique(object$Map$Subject[object$Map$NewSubj == i]),
                             " Raw Data")))
      
    }
    
  }
  
}