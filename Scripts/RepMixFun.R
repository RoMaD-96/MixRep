#   ____________________________________________________________________________
#   Libraries                                                               ####

library(ggplot2)
library(ggpubr)
library(colorspace)
library(tipmap)
library(spatstat)

#   ____________________________________________________________________________
#   Median Function                                                         ####

median_fun <- function(theta, tr, sr, to, so, null, priorsd, w) {
  # Ensure the data is sorted by theta
  density <-  rmapPostFix(theta = theta, tr = tr, sr = sr, to = to, so = so,
                          null = null, priorsd = priorsd, w = w)
  data <- data.frame(theta, density)
  data <- data[order(data$theta),]
  w_median <- weighted.median(data[,1], data[,2], na.rm = TRUE)
  
  # Return the median value of theta
  return(w_median)
}


#   ____________________________________________________________________________
#   Normalizing Constant Function (No Integration)                          ####

normConst <- function(tr, sr, to, so, null, priorsd, x, y) {
  mean_beta <- x / (x + y)
  marginal_lik <- (mean_beta * (dnorm(x = tr,
                                      mean = to,
                                      sd = sqrt(sr ^ 2 + so ^ 2)
  ) 
  - dnorm(x = tr,
          mean = null,
          sd = sqrt(sr ^ 2 + priorsd ^ 2)
  )) 
  + dnorm(x = tr,
          mean = null,
          sd = sqrt(sr ^ 2 + priorsd ^ 2))) 
  
  return(marginal_lik)
}



#   ____________________________________________________________________________
#   Posterior Density Functions                                             ####


##  ............................................................................
##  Posterior Distribution (Fixed Weights)                                  ####

rmapPostFix <- function(theta, w,  tr, sr, to, so, null, priorsd, x, y) {
  v1 <- 1/(1/so^2 + 1/sr^2)
  m1 <- (to/so^2 + tr/sr^2)*v1
  v2 <- 1/(1/priorsd^2 + 1/sr^2)
  m2 <- (null/priorsd^2 + tr/sr^2)*v2
  bf <- dnorm(x = tr, mean = to, sd = sqrt(sr^2 + so^2))/
    dnorm(x = tr, mean = null, sd = sqrt(sr^2 + priorsd^2))
  wUpdate <- 1/(1 + (1/w - 1)/bf)                     # Update based on the BF ?
  wUpdate*dnorm(x = theta, mean = m1, sd = sqrt(v1)) +
    (1 - wUpdate)*dnorm(x = theta, mean = m2, sd = sqrt(v2))
}


rmapPostFix_alt <- function(theta, w,  tr, sr, to, so, null, priorsd, x, y) {
  
  num_post <- dnorm( x = tr, mean = theta, sd =  sr)*(w * dnorm( x = theta, mean = to, sd = so) +
                                                        (1 - w) * dnorm(x = theta, mean = null, sd = priorsd))
  den_post <- w*dnorm(x = tr, mean = to, sd = sqrt(sr^2 + so^2)) + (1-w)*dnorm(x = tr, mean = null, sd = sqrt(sr^2 + priorsd^2))
  post <- num_post/den_post
}


##  ............................................................................
##  Marginal Likelihood                                                     ####

marginal_lik_fix_w <- function(theta, w,  tr, sr, to, so, null, priorsd, x, y) {
  
  den_post <- w*dnorm(x = tr, mean = to, sd = sqrt(sr^2 + so^2)) + (1-w)*dnorm(x = tr, mean = null, sd = sqrt(sr^2 + priorsd^2))
}

##  ............................................................................
##  Joint Posterior Distribution                                            ####

rmapPost <- function(theta, w,  tr, sr, to, so, null, priorsd, x, y) {
  num_post <- dnorm( x = tr, mean = theta, sd =  sr)*(w * dnorm( x = theta, mean = to, sd = so) +
              (1 - w) * dnorm(x = theta, mean = null, sd = priorsd)) * dbeta(x = w, shape1 = x, shape2 = y)
  den_post <- normConst(tr, sr, to, so, null, priorsd, x, y)
  joint_post <- num_post/den_post
  return(joint_post)
}


##  ............................................................................
##  Mixture Weights Marginal Posterior                                      ####

m_post_weights <- function(w, tr, sr, to, so, null, priorsd, x, y) {
  num_post_weights <-  (w*dnorm(x = tr, mean = to, sd = sqrt(so^2 + sr^2)) +
                        (1 - w)*dnorm(x = tr, mean = null, sd = sqrt(priorsd^2 + sr^2)))*dbeta(x = w, shape1 = x, shape2 = y)
  den_post_weights <- normConst(tr, sr, to, so, null, priorsd, x, y)
  m_post <- num_post_weights/den_post_weights
  return(m_post)
}



##  ............................................................................
##  Effect Size Marginal Posterior                                          ####

m_post_theta <- function(theta, tr, sr, to, so, null, priorsd, x, y) {
  mean_beta <- x / (x + y)
  num_post_theta <- (dnorm(x = tr, mean = theta, sd = sr) *
                      (mean_beta *(dnorm(x = theta, mean = to, sd = so)
                      - dnorm(x = theta, mean = null, sd = priorsd))
                      + dnorm(x = theta, mean = null, sd = priorsd)))
  den_post_theta <- normConst(tr, sr, to, so, null, priorsd, x, y)
  m_post <- num_post_theta / den_post_theta
  return(m_post)
}


#   ____________________________________________________________________________
#   HPDI Computation                                                        ####



##  ............................................................................
##  HPDI using Fixed Weights                                                ####

HPDI_post_m_theta_fix <- function(level, tr, sr, to, so, w, null, priorsd,
                              thetaRange = tr + c(-1, 1)*qnorm(p = (1 + level)/2)*sr*3,
                              quantileRange = c((1 - level)*0.2, (1 - level)*0.8)) {
  # Posterior Quantile Function
  quantile_fun <- function(q) {
    if (q == 0) {
      res <- -Inf
    } else if (q == 1) {
      res <- Inf
    } else {
      theta_dens_fun <- function(theta) {
        rmapPostFix(theta = theta, tr = tr, sr = sr, to = to, so = so,
                    null = null, priorsd = priorsd, w = w)
      }
      root_fun <- function(x) {
        integrate(f = theta_dens_fun, lower = -Inf, upper = x)$value - q
      }
      res <- uniroot(f = root_fun, interval = thetaRange)$root
    }
    return(res)
  }
  
  # Smallest HPDI
  opt_fun_scalar <- function(lower_q) {
    width <- quantile_fun(q = lower_q + level) - quantile_fun(q = lower_q)
    return(width)
  }
  opt_fun_vec <- Vectorize(FUN = opt_fun_scalar)
  opt_min_lower <- try(optim(par = (1 - level)/2, fn = opt_fun_vec,
                             method = "L-BFGS-B", lower = quantileRange[1],
                             upper = quantileRange[2])$par)
  if (inherits(opt_min_lower, "try-error")) {
    cr_int <- c("lower" = NaN, "upper" = NaN)
  } else {
    cr_int <- c("lower" = quantile_fun(q = opt_min_lower),
                "upper" = quantile_fun(q = opt_min_lower + level))
  }
  return(cr_int)
}

##  ............................................................................
##  Weights HPDI                                                            ####

HPDI_post_m_weights <- function(tr, sr, to, so, x, y, null, priorsd, level = 0.95) {
  # Posterior Quantile Function
  quantile_fun <- function(q) {
    if (q == 0) {
      res <- 0
    } else if (q == 1) {
      res <- 1
    } else {
      m_dens_fun <- function(w) {
        m_post_weights(w = w, tr = tr, sr = sr, to = to, so = so,
                       x = x, y = y, null = null, priorsd = priorsd)
      }
      root_fun <- function(x) {
        integrate(f = m_dens_fun, lower = 0, upper = x)$value - q
      }
      res <- uniroot(f = root_fun, interval = c(0, 1))$root
    }
    return(res)
  }
  
  # Smallest HPDI
  opt_fun_scalar <- function(lower_q) {
    width <- quantile_fun(q = lower_q + level) - quantile_fun(q = lower_q)
    return(width)
  }
  opt_fun_vec <- Vectorize(FUN = opt_fun_scalar)
  opt_min_lower <- try(optim(par = (1 - level)/2, fn = opt_fun_vec,
                             method = "L-BFGS-B", lower = 0,
                             upper = 1 - level)$par)
  if (inherits(opt_min_lower, "try-error")) {
    cr_int <- c("lower" = NaN, "upper" = NaN)
  } else {
    cr_int <- c("lower" = quantile_fun(q = opt_min_lower),
            "upper" = quantile_fun(q = opt_min_lower + level))
  }
  return(cr_int)
}


##  ............................................................................
##  Effect Size HPDI                                                        ####

HPDI_post_m_theta <- function(level, tr, sr, to, so, x = 1, y = 1, null, priorsd,
                           thetaRange = tr + c(-1, 1)*qnorm(p = (1 + level)/2)*sr*3,
                           quantileRange = c((1 - level)*0.2, (1 - level)*0.8)) {
  # Posterior Quantile Function
  quantile_fun <- function(q) {
    if (q == 0) {
      res <- -Inf
    } else if (q == 1) {
      res <- Inf
    } else {
      m_dens_fun <- function(theta) {
        m_post_theta(theta = theta, tr = tr, sr = sr, to = to, so = so,
                     x = x, y = y,  null = null, priorsd = priorsd)
      }
      root_fun <- function(x) {
        integrate(f = m_dens_fun, lower = -Inf, upper = x)$value - q
      }
      res <- uniroot(f = root_fun, interval = thetaRange)$root
    }
    return(res)
  }
  
  # Smallest HPDI
  opt_fun_scalar <- function(lower_q) {
    width <- quantile_fun(q = lower_q + level) - quantile_fun(q = lower_q)
    return(width)
  }
  opt_fun_vec <- Vectorize(FUN = opt_fun_scalar)
  opt_min_lower <- try(optim(par = (1 - level)/2, fn = opt_fun_vec,
                               method = "L-BFGS-B", lower = quantileRange[1],
                               upper = quantileRange[2])$par)
  if (inherits(opt_min_lower, "try-error")) {
    cr_int <- c("lower" = NaN, "upper" = NaN)
  } else {
    cr_int <- c("lower" = quantile_fun(q = opt_min_lower),
            "upper" = quantile_fun(q = opt_min_lower + level))
  }
  return(cr_int)
}
