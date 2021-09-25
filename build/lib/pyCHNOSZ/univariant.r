# Univariant Curve option in R --------------------------------------------------------------------
# Summary: functions to calculate temperatures at given logK and pressures, calculate
#          pressures at given logK and temperatures, and to export and plot results.
# Author: Grayson Boyer
# Last updated: September 20, 2021

########### UNIVARIANT CURVE FUNCTIONS ###########

uc_solveT <- function(logK, species, state, coeff,
                      pressures = 1, IS=0, minT = 0.1, maxT = 100, tol=0.00001){

  user_minT <- minT
  user_maxT <- maxT

  # create a data frame to hold results
  PlH <- rep(NA, length(pressures)) # placeholder

  df <- data.frame(P = sprintf("%.3f", round(pressures, 3)),
                   T = PlH,
                   rho = PlH,
                   logK = PlH,
                   G = PlH,
                   H = PlH,
                   S = PlH,
                   V = PlH,
                   Cp = PlH,
                   stringsAsFactors=FALSE
                 )

  for (pressure in pressures){

    minT <- user_minT
    maxT <- user_maxT
    completed_iter <- 0
    converged <- FALSE

    while(!converged){
      if(completed_iter > 100){
        warning("Too many iterations (>100). Terminating calculation.")
        break
      }

      # perform a subcrt calculation using a guess
      found_guesslogK <- FALSE
      iter_guesslogK <- 0
      while(!found_guesslogK){
        if(iter_guesslogK == 0){
          guessT <- mean(c(minT, maxT))
        }else if(iter_guesslogK <= 20){
          guessT <- seq(minT, maxT, length.out=20)[iter_guesslogK]
        }else{
          warning("Could not find a viable temperature guess in this range after 20 tries. Terminating calculation.")
          break
        }
            
        guess_calc <- suppressMessages(subcrt(species, state, coeff, T = guessT, P = pressure, IS = IS, exceed.Ttr = T))
        guesslogK <- guess_calc$out$logK
          
        if(!is.na(guesslogK) && is.finite(guesslogK)){
          found_guesslogK <- TRUE
        }else{
          iter_guesslogK <- iter_guesslogK + 1
        }
      }
    
      if(is.na(guesslogK) | !is.finite(guesslogK)){
        break
      }

      # check if initial guess is close enough to actual logK within tolerance
      if(abs(logK - guesslogK) < tol){
        this_P <- sprintf("%.3f", round(pressure, 3))
        df[df[, "P"] == this_P, "T"] <- sprintf("%.3f", round(guessT, 3))
        
        if("rho" %in% names(guess_calc$out)){
          df[df[, "P"] == this_P, "rho"] <- sprintf("%.3f", round(guess_calc$out$rho, 3))
        }
        df[df[, "P"] == this_P, "logK"] <- sprintf("%.3f", round(guesslogK, 3))
        df[df[, "P"] == this_P, "G"] <- sprintf("%.0f", round(guess_calc$out$G, 0))
        df[df[, "P"] == this_P, "H"] <- sprintf("%.0f", round(guess_calc$out$H, 0))
        df[df[, "P"] == this_P, "S"] <- sprintf("%.1f", round(guess_calc$out$S, 1))
        df[df[, "P"] == this_P, "V"] <- sprintf("%.1f", round(guess_calc$out$V, 1))
        df[df[, "P"] == this_P, "Cp"] <- sprintf("%.1f", round(guess_calc$out$Cp, 1))

        result <- list(reaction = guess_calc$reaction, out = df)

        converged <- TRUE
        break
      }

      # perform an initial calculation across a range of temperatures bounded by current minT and maxT
      init_calc <- suppressMessages(subcrt(species, state, coeff, T = seq(minT, maxT, length.out = 200), P = pressure, IS=IS, exceed.Ttr=T)$out)

      # Check if logK falls between any of the temperature iterations in the initial calculation
      logK_check_complete <- FALSE
      for (i in 1:(length(init_calc$logK)-1)){

        logKmin <- init_calc$logK[i]
        logKmax <- init_calc$logK[i+1]
        
        
        if ((!is.finite(logKmin)) | (!is.finite(logKmax))){
          if(!("Warning" %in% colnames(df))){
            df$Warning <-PlH # create 'Warning' column
          }

          this_P <- sprintf("%.3f", round(pressure, 3))
          df[df[, "P"] == this_P, "Warning"] <- paste("Could not converge on T for this P within", user_minT, "and", user_maxT, "degC")
            
          result <- list(reaction = guess_calc$reaction, out = df)
          converged <- TRUE
        } else if ((logK <= 0) & (logK <= logKmin) & (logK >= logKmax)){
          minT <- init_calc$T[i]
          maxT <- init_calc$T[i+1]
          logK_check_complete <- TRUE
          break  
        } else if ((logK >= 0) & (logK <= logKmin) & (logK >= logKmax)){
          minT <- init_calc$T[i]
          maxT <- init_calc$T[i+1]
          logK_check_complete <- TRUE
          break
        } else if ((logK < 0) & (logK >= logKmin) & (logK <= logKmax)){
          minT <- init_calc$T[i]
          maxT <- init_calc$T[i+1]
          logK_check_complete <- TRUE
          break
        } else if ((i == (length(init_calc$logK)-1)) & (!logK_check_complete)){
          if(!("Warning" %in% colnames(df))){
            df$Warning <-PlH # create 'Warning' column
          }

          this_P <- sprintf("%.3f", round(pressure, 3))
          df[df[, "P"] == this_P, "Warning"] <- paste("Could not converge on T for this P within", user_minT, "and", user_maxT, "degC")
            
          result <- list(reaction = guess_calc$reaction, out = df)
          converged <- TRUE
        }


      } # end logK for loop

      completed_iter <- completed_iter + 1 # increase iteration counter

    } # end !converged while loop

  } # end pressure for loop
    
  # convert columns to numeric
  df <- result[["out"]]
  cols.num <- c("P", "T", "logK", "G", "H", "S", "V", "Cp")

  if(!("rho" %in% names(guess_calc$out))){
    df <- df[ , -which(names(df) %in% c("rho"))] # drop the empty rho column
  }else{
    cols.num <- c(cols.num, "rho")
  }

  df[cols.num] <- sapply(df[cols.num], as.numeric)
  result[["out"]] <- df
    
  return(result)

} # end uc_solveT function




uc_solveP <- function(logK, species, state, coeff,
                      temperatures = 1, IS=0,
                      minP = 1, maxP = 500, tol=0.00001){
    
  user_minP <- minP
  user_maxP <- maxP

  # create a data frame to hold results
  PlH <- rep(NA, length(temperatures)) # placeholder

  df <- data.frame(P = PlH,
                   T = sprintf("%.3f", round(temperatures, 3)),
                   rho = PlH,
                   logK = PlH,
                   G = PlH,
                   H = PlH,
                   S = PlH,
                   V = PlH,
                   Cp = PlH,
                   stringsAsFactors=FALSE
                 )

  for (temperature in temperatures){

    minP <- user_minP
    maxP <- user_maxP
    completed_iter <- 0
    converged <- FALSE

    while(!converged){
      if(completed_iter > 100){
        print("Too many iterations (>100). Terminating calculation.")
        break
      }

      # perform a subcrt calculation using a guess
      guessP <- mean(c(minP, maxP))
      guess_calc <- suppressMessages(subcrt(species, state, coeff, P = guessP, T = temperature, IS=IS))
      guesslogK <- guess_calc$out$logK

      # check if initial guess is close enough to actual logK within tolerance
      if(abs(logK - guesslogK) < tol){
        this_T <- sprintf("%.3f", round(temperature, 3))
        df[df[, "T"] == this_T, "P"] <- sprintf("%.3f", round(guessP, 3))
        df[df[, "T"] == this_T, "rho"] <- sprintf("%.3f", round(guess_calc$out$rho, 3))
        df[df[, "T"] == this_T, "logK"] <- sprintf("%.3f", round(guesslogK, 3))
        df[df[, "T"] == this_T, "G"] <- sprintf("%.0f", round(guess_calc$out$G, 0))
        df[df[, "T"] == this_T, "H"] <- sprintf("%.0f", round(guess_calc$out$H, 0))
        df[df[, "T"] == this_T, "S"] <- sprintf("%.1f", round(guess_calc$out$S, 1))
        df[df[, "T"] == this_T, "V"] <- sprintf("%.1f", round(guess_calc$out$V, 1))
        df[df[, "T"] == this_T, "Cp"] <- sprintf("%.1f", round(guess_calc$out$Cp, 1))

        result <- list(reaction = guess_calc$reaction, out = df)

        converged <- TRUE
        break
      }

      # perform an initial calculation across a range of temperatures bounded by current minT and maxT
      init_calc <- suppressMessages(subcrt(species, state, coeff, P = seq(minP, maxP, length.out = 10), T = temperature, IS=IS)$out)


      # Check if logK falls between any of the temperature iterations in the initial calculation
      logK_check_complete <- FALSE
      for (i in 1:(length(init_calc$logK)-1)){

        logKmin <- init_calc$logK[i]
        logKmax <- init_calc$logK[i+1]
        if ((!is.finite(logKmin)) | (!is.finite(logKmax))){
          if(!("Warning" %in% colnames(df))){
            df$Warning <-PlH # create 'Warning' column
          }

          this_T <- sprintf("%.3f", round(pressure, 3))
          df[df[, "T"] == this_T, "Warning"] <- paste("Could not converge on P for this T within", user_minP, "and", user_maxP, "bars")
            
          result <- list(reaction = guess_calc$reaction, out = df)
          converged <- TRUE
        }else if ((logK >= 0) & (logK <= logKmin) & (logK >= logKmax)){
          minP <- init_calc$P[i]
          maxP <- init_calc$P[i+1]
          logK_check_complete <- TRUE
          break
        } else if ((logK <= 0) & (logK <= logKmin) & (logK >= logKmax)){
          minP <- init_calc$P[i]
          maxP <- init_calc$P[i+1]
          logK_check_complete <- TRUE
          break
        } else if ((logK > 0) & (logK >= logKmin) & (logK <= logKmax)){
          minP <- init_calc$P[i]
          maxP <- init_calc$P[i+1]
          logK_check_complete <- TRUE
          break
        } else if ((i == (length(init_calc$logK)-1)) & (!logK_check_complete)){
          if(!("Warning" %in% colnames(df))){
            df$Warning <-PlH # create 'Warning' column
          }

          this_T <- sprintf("%.3f", round(temperature, 3))
          df[df[, "T"] == this_T, "Warning"] <- paste("Could not converge on P for this T within", user_minP, "and", user_maxP, "bars")

          result <- list(reaction = guess_calc$reaction, out = df)
          converged <- TRUE
        }


      } # end logK for loop

      completed_iter <- completed_iter + 1 # increase iteration counter

    } # end !converged while loop

  } # end pressure for loop
    
  # convert columns to numeric
  df <- result[["out"]]
  cols.num <- c("P", "T", "logK", "G", "H", "S", "V", "Cp")
    
  if(!("rho" %in% names(guess_calc$out))){
    df <- df[ , -which(names(df) %in% c("rho"))] # drop the empty rho column
  }else{
    cols.num <- c(cols.num, "rho")
  }
    
  df[cols.num] <- sapply(df[cols.num], as.numeric)
  result[["out"]] <- df
    
  return(result)

} # end uc_solveP function

write_csv_output <- function(result, create_output_csv=F){
    # write csv table
    if(create_output_csv){
      write.csv(result$out, csv_filename, row.names=F)
    }
}

### Create plot temperature-based
create_output_plot_T <- function(logK, species, state, coeff, pressures, minT, maxT, res=300){
    # print(as.numeric(as.character(result$out$T)))
    calc <- subcrt(species, state, coeff, T=seq(minT, maxT, length.out=res), P=pressures[1])$out$logK

    if(min(calc, na.rm=T) != Inf){
        if(logK < min(calc, na.rm=T)){
            this_ylim <- c(logK-1, max(calc, na.rm=T))
        } else if(logK > max(calc, na.rm=T)){
            this_ylim <- c(min(calc, na.rm=T), logK+1)
        } else {
            this_ylim <- c(min(calc, na.rm=T), max(calc, na.rm=T))
        }
    }

    plot(x=NA,
         y=NA,
         xlim=c(minT-10, maxT+10),
         ylim=this_ylim,
         ylab="logK",
         xlab=expression("Temperature " ( degree*C)),
         type="l")
    grid (NULL, NULL, lty = 1, col = "lightgray")
    lines(x=seq(minT, maxT, length.out=res), y=rep(logK, res), col="red", lwd=3)
    for(pressure in pressures){
      lines(x=seq(minT, maxT, length.out=res),
            y=subcrt(species, state, coeff, T=seq(minT, maxT, length.out=res), P=pressure)$out$logK)
    }
    
}

### Create plot pressure-based
create_output_plot_P <- function(logK, species, state, coeff, temperatures, minP, maxP, res=300){
    #print(as.numeric(as.character(result$out$T)))
    calc <- subcrt(species, state, coeff, P=seq(minP, maxP, length.out=res), T=temperatures[1])$out$logK

    if(min(calc, na.rm=T) != Inf){
        if(logK < min(calc, na.rm=T)){
            this_ylim <- c(logK-1, max(calc, na.rm=T))
        } else if(logK > max(calc, na.rm=T)){
            this_ylim <- c(min(calc, na.rm=T), logK+1)
        } else {
            this_ylim <- c(min(calc, na.rm=T), max(calc, na.rm=T))
        }
    }

    plot(x=NA,
         y=NA,
         xlim=c(minP-1, maxP+1),
         ylim=this_ylim,
         ylab="logK",
         xlab="Pressure (bars)",
         type="l")
    grid (NULL, NULL, lty = 1, col = "lightgray")
    lines(x=seq(minP, maxP, length.out=res), y=rep(logK, res), col="red", lwd=3)
    for(temperature in temperatures){
      lines(x=seq(minP, maxP, length.out=res),
            y=subcrt(species, state, coeff, P=seq(minP, maxP, length.out=res), T=temperature)$out$logK)
    }
    
}


thermoinfo <- function(species){
    return(info(info(species)))
}