# Utilities for tauFisher

########################################
# Convert time to 24-hours
########################################

#' Converts time to [0,24)
#'
#' @param time vector of time
#' @return vector of time values within 24 hours
#' @export
#'
convert_to_24hr <- function(time){
  time = as.numeric(time)
  time24 = time - 24*floor(time/24)

  return(time24)
}


########################################
# Convert (x,y) to angular coordinates
########################################

#' Get angles of Cartesian coordinates
#'
#' @param x x values
#' @param y y values
#' @return angles of (x,y)
ang <- function(x,y) {
  z <- x + 1i * y
  res <- 90 - Arg(z) / pi * 180
  return(res %% 360)
}

#' Converts Cartesian coordinates to polar
#'
#' @param x x values
#' @param y y values
#' @return data frame with columns angles, angular x and y coordinates
#' @export
convert_to_polar <- function(x, y){
  angles = ang(x, y)
  ang_x = cos(pi*angles/180)
  ang_y = sin(pi*angles/180)

  ang_df = as.data.frame(cbind(angles, ang_x, ang_y))
  return(ang_df)
}

########################################
# Calculate error
########################################

#' Calculate error of time prediction
#'
#' @param truth vector of true circadian times
#' @param pred vector of predicted circadian times
#' @return List of
#'         (1) data frame with the true time, predicted time, and differences
#'         (2) root mean squared error
#'         (3) accuracy, with correct defined as predicted within 2hrs of truth
#' @export
calc_error <- function(truth, pred){

  truth_24 = truth - 24*floor(truth/24)

  # convert to possible overlap times and subtract
  diff = truth_24 - pred
  for (i in 1:length(diff)){
    if (diff[i] >12){
      diff[i] = diff[i]-24
    } else{
      if (diff[i] < -12)
        diff[i] = diff[i]+24
    }
  }

  # calc RMSE
  rmse = sqrt(mean(diff^2))

  # calc accuracy (correct defined as <=2 hours)
  acceptable = ifelse(abs(diff) <=2, 1, 0)
  acceptable_pct = sum(acceptable)/length(acceptable)

  d_df = data.frame("truth" = truth, "truth_24" = truth_24,
                    "predicted" = pred, "differences" = diff)

  return(list("differences" = d_df, "RMSE" = rmse, "Accuracy" = acceptable_pct))
}
