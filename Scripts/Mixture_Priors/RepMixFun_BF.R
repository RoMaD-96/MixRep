#   ____________________________________________________________________________
#   Libraries                                                               ####

library(ggplot2)
library(ggpubr)
library(colorspace)
library(spatstat)
library(repmix)

#   ____________________________________________________________________________
#   Median Function                                                         ####

median_fun <- function(theta, tr, sr, to, so, m, v, w) {
  # Ensure the data is sorted by theta
  density <- thetaposteriormix(theta = thetaseq, tr = tr, sr = sr, to = to, so = so,
                               m = mu_UIP, v = tau_UIP, w = w)
  data <- data.frame(theta, density)
  data <- data[order(data$theta),]
  w_median <- weighted.median(data[,1], data[,2], na.rm = TRUE)
  
  # Return the median value of theta
  return(w_median)
}


#   ____________________________________________________________________________
#   Bayes Factor                                                            ####


##  ............................................................................
##  Format Bayes Factor                                                     ####

formatBF <- function(BF, digits = "default") {
  ## check inputs
  stopifnot(
    length(BF) == 1,
    is.numeric(BF),
    (is.finite(BF) && 0 < BF) || is.na(BF),
    
    length(digits) == 1,
    (is.character(digits) && digits == "default") ||
      (is.numeric(digits) && 0 <= digits)
  )
  ## return NA if input NA/NaN
  if (is.na(BF) || is.nan(BF))
    result <- NA
  else {
    ## format BF
    if (digits == "default") {
      if (BF < 1/1000)
        result <- "< 1/1000"
      if ((BF >= 1/1000) & (BF <= 1/10))
        result <- paste0("1/", as.character(round(1/BF)))
      if ((BF > 1/10) & (BF < 1))
        result <- paste0("1/", as.character(round(1/BF, digits = 1)))
      if ((BF < 10) & (BF >= 1))
        result <- as.character(round(BF, digits = 1))
      if ((BF >= 10) & (BF <= 1000))
        result <- as.character(round(BF))
      if (BF > 1000)
        result <- "> 1000"
    } else {
      if (BF < 1)
        result <- paste0("1/", as.character(round(1/BF, digits = digits)))
      else
        result <- as.character(round(BF, digits = digits))
    }
    ## when 1/1 return 1
    if (result == "1/1") result <- "1"
  }
  return(result)
}
formatBF_vec <- Vectorize(FUN = formatBF)


##  ............................................................................
##  BF Theta                                                                ####

bf_theta_mix <- function(tr, sr, to, so, x = 1, y = 1, m = m,
                         v = v, w = NA) {
  ## marginal density under H0
  null_H_0 <- dnorm(x = tr, mean = 0, sd = sr)

  ## marginal density under H1
  if (!is.na(w)) {
    alt_H_1 <- w * dnorm(x = tr, mean = to, sd = sqrt(sr^2 + so^2)) + 
              (1 - w) * dnorm(x = tr, mean = m, sd = sqrt(sr^2 + v))
  } else {
    alt_H_1 <- marglik(
      tr = tr, to = to, sr = sr, so = so, m = m,
      v = v, x = x, y = y)
  }

  ## observed BF
  obs_bf <- null_H_0 / alt_H_1
  return(obs_bf)
}


##  ............................................................................
##  BF Omega                                                                ####


bf_omega_mix <- function(tr, sr, to, so, x = x, y = y, m = m,
                         v = v, w_null = NA, w_alt = NA){

  
  if (!is.na(w_null) && !is.na(w_alt)) {
    
    ## marginal density under H1: omega = 1
    alt_H_1 <- marglik(tr = tr, to = to, sr = sr, so = so, m = m,
                       v = v, w = w_alt)
    
    ## marginal density under H0: omega = 0
    null_H_0 <- marglik(tr = tr, to = to, sr = sr, so = so, m = m,
                        v = v, w = w_null)
    
  } else if (x == 1 && y > 1) {
    ## marginal density under H1: omega = 1
    alt_H_1 <- dnorm(x = tr, mean = to, sd = sqrt(sr^2 + so^2))
    
    ## marginal density under H0
    ## integrating likelihood over alpha|H1 ~ Beta(1, y) prior
    null_H_0 <- marglik(tr = tr, to = to, sr = sr, so = so, m = m,
                        v = v, x = 1, y = y)
  } else if (x > 1 && y == 1) {
    ## marginal density under H1: omega = 0
    alt_H_1 <- dnorm(x = tr, mean = m, sd = sqrt(sr^2 + v))
    
    ## marginal density under H0
    ## integrating likelihood over alpha|H1 ~ Beta(x, 1) prior
    null_H_0 <- marglik(tr = tr, to = to, sr = sr, so = so, m = m,
                        v = v, x = x, y = 1)
  } else {
    alt_H_1 <- NULL
    null_H_0 <- NULL
  }
  
  
  ## compute BF
  obs_bf <- null_H_0/alt_H_1
  return(obs_bf)
}












