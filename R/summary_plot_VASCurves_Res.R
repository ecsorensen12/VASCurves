VASCurves_Res <- function(object, ...){
  UseMethod("VASCurves_Res")
}

#' Method for summary function in VASCurves
#' 
#' Summary method for an object of class VASCurves_Res.
#' @param object An object of class VASCurves_Res.
#' @param digits Number of digits to print for numeric variables.
#' @param ... Other parameters.
#' @return An object of summary statistics of the data.
#' @method summary VASCurves_Res
#' @export
summary.VASCurves_Res <- function(object, digits = 3, ...){
  
  res <- apply(object$result.summary, 1, function(x){
    paste0(round(x["mean"], digits), " (", round(x["2.5%"], digits), ", ", round(x["97.5%"], digits), ")")
  })
  full_params <- object$result.summary
  beta <- full_params[grepl("beta", rownames(full_params)),1]
  B0 <- full_params["B[1]",1]
  B1 <- full_params["B[2]",1]
  
  if(object$data$type == "No Continuum"){
    
    indiv_params_full <- data.frame(matrix(nrow = object$data$N_subj, ncol = 6))
    for(i in 1:object$data$N_subj){
      tryCatch(expr = {
        U_i <- full_params[paste0("U[", i, ",", 1:4, "]"),1]
        Usig_i <- full_params[paste0("usig_i[", i, "]"),1]
        indiv_params <- beta + U_i
        indiv_var <- exp(B0 + B1*indiv_params[3] + Usig_i)
        indiv_params_full[i,] <- c(object$data$Map$Original[object$data$Map$New == i],
                                   round(indiv_params, digits),
                                   round(indiv_var, digits))
      },
      error = function(e){NA})
    }
    colnames(indiv_params_full) <- c("ID", "Min", "Max", "Slope", "X0", "Variance")
    
    return_list <- list(unname(res["beta[1]"]),
                        unname(res["beta[2]"]),
                        unname(res["beta[3]"]),
                        unname(res["beta[4]"]),
                        unname(res["B[1]"]),
                        unname(res["B[2]"]),
                        indiv_params_full)
    names(return_list) <- c("Min", "Max", "Slope", "X0", "B0", "B1", "Indiv_Params")
    
  }else{
    
    cont_params_full <- data.frame(matrix(nrow = object$data$N_cont, ncol = 5))
    for(j in 1:nrow(cont_params_full)){
      Cj <- full_params[paste0("C[", j, ",", 1:4, "]"),1]
      cont_params_full[j,] <- c(unique(object$data$Map$Continuum[object$data$Map$NewCont == j]),
                                round(beta + Cj, digits))
    }
    colnames(cont_params_full) <- c("Continuum", "Min", "Max", "Slope", "X0")
    
    indiv_cont_params_full <- data.frame(matrix(nrow = object$data$N_subj_cont, ncol = 7))
    for(i in 1:object$data$N_subj_cont){
      tryCatch(expr = {
        Cj <- full_params[paste0("C[", object$data$cont_idx[i], ",", 1:4, "]"),1]
        UC_ij <- full_params[paste0("UC[", i, ",", 1:4, "]"),1]
        Csigj <- full_params[paste0("Csig[", object$data$cont_idx[i], "]"),1]
        Usig_ij <- full_params[paste0("UCsig[", i, "]"),1]
        indiv_params <- beta + Cj + UC_ij
        indiv_var <- exp(B0 + B1*indiv_params[3] + Csigj + Usig_ij)
        indiv_cont_params_full[i,] <- c(object$data$Map$Subject[object$data$Map$NewSubjCont == i],
                                        object$data$Map$Continuum[object$data$Map$NewSubjCont == i],
                                        round(indiv_params, digits),
                                        round(indiv_var, digits))
      },
      error = function(e){NA})
    }
    colnames(indiv_cont_params_full) <- c("ID", "Continuum", "Min", "Max", "Slope", "X0", "Variance")
    
    return_list <- list(unname(res["beta[1]"]),
                        unname(res["beta[2]"]),
                        unname(res["beta[3]"]),
                        unname(res["beta[4]"]),
                        unname(res["B[1]"]),
                        unname(res["B[2]"]),
                        cont_params_full,
                        indiv_cont_params_full)
    names(return_list) <- c("Min", "Max", "Slope", "X0", "B0", "B1", "Cont_Params", "Indiv_Params")
  }
  
  class(return_list) <- "summary.VASCurves_Res"
  
  return(return_list)
}

#' Print method for summary function in VASCurves
#' 
#' Method to print summary object returned from summary.VASCurves_Res.
#' @param object An object of class VASCurves_Res.
#' @param ... Other parameters.
#' @return An object of summary statistics of the data.
#' @method print summary.VASCurves_Res
#' @export  
print.summary.VASCurves_Res <- function(object, ...){
  
  min_string <- paste0("Min: ", object$Min)
  max_string <- paste0("Max: ", object$Max)
  slope_string <- paste0("Slope: ", object$Slope)
  x0_string <- paste0("X0: ", object$X0)
  B0_string <- paste0("B0 (Log Var. Intercept): ", object$B0)
  B1_string <- paste0("B1 (Log Var. Slope Effect): ", object$B1)
  
  cat(min_string,
      "\n",
      max_string,
      "\n",
      slope_string,
      "\n",
      x0_string,
      "\n",
      B0_string,
      "\n",
      B1_string,
      sep = "")
}

#' Method for plot function in VASCurves
#' 
#' Plotting method for an object of class VASCurves_Res. Plots the raw data for each individual in the dataset along with fitted curve(s).
#' @param object An object of class VASCurves_Res.
#' @param subj_ids Vector of string values detailing which subject ids to plot respective data. Defaults to "all" which gives a plot for all subjects.
#' @param ... Other graphical parameters.
#' @return Plots of the raw data for specified individuals in the input VASCurves_Res object as well as fitted curves.
#' @method plot VASCurves_Res
#' @import ggplot2
#' @export
plot.VASCurves_Res <- function(object, subj_ids = "all", ...){
  
  full_params <- object$result.summary
  beta <- full_params[grepl("beta", rownames(full_params)),1]
  B0 <- full_params["B[1]",1]
  B1 <- full_params["B[2]",1]
  
  vas.f <- function(min,max,s,x0,time){(max-min)/(1+exp((4*s/(max-min))*(x0-time)))+min}
  vas_wrapper <- function(beta, time){
    return(vas.f(beta[1],
                 beta[2],
                 beta[3],
                 beta[4],
                 time))
  }
  
  if(object$data$type == "No Continuum"){
    
    if("all" %in% subj_ids){
      subj_idx <- 1:object$data$N_subj
    }else{
      subj_ids <- as.character(subj_ids)
      if(all(subj_ids %in% object$data$Map$Original)){
        subj_idx <- object$data$Map$New[object$data$Map$Original %in% subj_ids]
      }else{
        stop("Input IDs not found in VASCurves_Res input object!")
      }
    }
    
    ### Print plot of raw data for subjects
    for(i in subj_idx){
      
      temp_data <- list()
      
      temp_data$Y <- object$data$Y[object$data$subj_idx == i]
      temp_data$time <- object$data$time[object$data$timeindex[i,object$data$timeindex2[i,1]:object$data$timeindex2[i,2]]]
      U_i <- full_params[paste0("U[", i, ",", 1:4, "]"),1]
      indiv_params <- beta + U_i
      temp_data$Pred_Y <- vas_wrapper(indiv_params, seq(min(temp_data$time), 
                                                        max(temp_data$time),
                                                        (max(temp_data$time) - min(temp_data$time))/1000))
      temp_data$Pred_Time <- seq(min(temp_data$time), 
                                 max(temp_data$time),
                                 (max(temp_data$time) - min(temp_data$time))/1000)
      
      print(ggplot() +
              geom_line(data = data.frame(Time = temp_data$Pred_Time,
                                          Value = temp_data$Pred_Y),
                        aes(x = Time, y = Value), col = "red") +
              geom_point(aes(x = temp_data$time, y = temp_data$Y), col = "black") +
              xlim(min(temp_data$time, na.rm = T), 
                   max(temp_data$time, na.rm = T)) +
              ylim(min(c(temp_data$Y, temp_data$Pred_Y), na.rm = T), 
                   max(c(temp_data$Y, temp_data$Pred_Y), na.rm = T)) +
              ggtitle(paste0(object$data$Map$Original[i],
                             " Raw Data and Estimated Curve")))
      
    }
    
  }else{
    
    if("all" %in% subj_ids){
      subj_idx <- 1:object$data$N_subj
    }else{
      subj_ids <- as.character(subj_ids)
      if(all(subj_ids %in% object$data$Map$Subject)){
        subj_idx <- unique(object$data$Map$NewSubj[object$data$Map$Subject %in% subj_ids])
      }else{
        stop("Input IDs not found in VASCurves_Res input object!")
      }
    }
    
    ### Print plot of raw data for subjects across continua
    for(i in subj_idx){
      
      temp_data <- list()
      
      subjcont_idxs <- object$data$Map$NewSubjCont[object$data$Map$NewSubj == i]
      continuums <- object$data$Map$Continuum[object$data$Map$NewSubj == i]
      temp_data$Y <- object$data$Y[object$data$subjcont_idx %in% subjcont_idxs]
      
      times <- c()
      continuums_full <- c()
      Pred_Y <- c()
      Pred_Time <- c()
      Pred_Continuum <- c()
      for(j in 1:length(subjcont_idxs)){
        times <- c(times,
                   object$data$time[object$data$timeindex[subjcont_idxs[j],
                                                          object$data$timeindex2[subjcont_idxs[j],1]:object$data$timeindex2[subjcont_idxs[j],2]]])
        continuums_full <- c(continuums_full,
                             rep(continuums[j], 
                                 length(object$data$timeindex2[subjcont_idxs[j],1]:object$data$timeindex2[subjcont_idxs[j],2])))
        
        Cj <- full_params[paste0("C[", object$data$cont_idx[subjcont_idxs[j]], ",", 1:4, "]"),1]
        UC_ij <- full_params[paste0("UC[", subjcont_idxs[j], ",", 1:4, "]"),1]
        indiv_params <- beta + Cj + UC_ij
        time_seq <- seq(min(times), 
                        max(times),
                        (max(times) - min(times))/1000)
        Pred_Y <- c(Pred_Y, vas_wrapper(indiv_params, time_seq))
        Pred_Time <- c(Pred_Time, time_seq)
        Pred_Continuum <- c(Pred_Continuum,
                            rep(continuums[j],
                                length(time_seq)))
      }
      temp_data$time <- times
      temp_data$continuum <- continuums_full
      temp_data$Pred_Y <- Pred_Y
      temp_data$Pred_Time <- Pred_Time
      temp_data$Pred_Continuum <- Pred_Continuum
      
      print(ggplot() +
              geom_line(data = data.frame(Value = temp_data$Pred_Y,
                                          Time = temp_data$Pred_Time,
                                          Continuum = temp_data$Pred_Continuum),
                        aes(x = Time, y = Value, col = Continuum)) +
              geom_point(data = data.frame(Value = temp_data$Y,
                                           Time = temp_data$time,
                                           Continuum = temp_data$continuum),
                         aes(x = Time, y = Value, col = Continuum)) +
              xlim(min(temp_data$time, na.rm = T), 
                   max(temp_data$time, na.rm = T)) +
              ylim(min(c(temp_data$Y, temp_data$Pred_Y), na.rm = T), 
                   max(c(temp_data$Y, temp_data$Pred_Y), na.rm = T)) +
              ggtitle(paste0(unique(object$data$Map$Subject[object$data$Map$NewSubj == i]),
                             " Raw Data and Estimated Curves")))
      
    }
    
  }
  
}