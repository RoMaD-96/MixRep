#   ____________________________________________________________________________
#   Libraries                                                               ####

library(dplyr)
library(knitr)
library(xtable)
library(repmix)

source("Scripts/Mixture_Priors/RepMixFun_BF.R")

#   ____________________________________________________________________________
#   Random Weights                                                          ####

load("credentials_data.RData")


##  ............................................................................
##  Parameter Setting                                                       ####

# Original and Replicated Studies
to <- data %>%
  dplyr::filter(type == "original") %>%
  dplyr::pull(fis) %>%
  as.numeric()
so <- data %>%
  dplyr::filter(site == "original") %>%
  dplyr::pull(se_fis) %>%
  as.numeric()

trep <- data %>%
  dplyr::filter(site %in% c("University of Toronto", "Montana State University", "Ashland University")) %>%
  dplyr::pull(fis) %>%
  as.numeric()

srep <- data %>%
  dplyr::filter(site %in% c("University of Toronto", "Montana State University", "Ashland University")) %>%
  dplyr::pull(se_fis) %>%
  as.numeric()

tp <- round(sum(trep/srep^2)/sum(1/srep^2),2)
sp <- round(sqrt(1/sum(1/srep^2)),2)

tr <- c(trep,tp)
sr <- c(srep,sp)

# Mean and Variance Unit Informative Prior
mu_UIP <- 0
tau_UIP <- 2

# Parameter Grid
n_weights <- 400
n_theta <- 400
wseq <- seq(0, 1, length.out = n_weights)
thetaseq <- seq(-0.9, 0.9, length.out = n_theta)
par_grid <- expand.grid(omega = wseq, theta = thetaseq)

# Uniform Prior for the Weight
eta <- 1
nu <- 1

# Replication Number
rep_number <- c(1,2,3,4)


##  ............................................................................
##  Joint Posterior                                                         ####

postdens <- posteriormix(theta = par_grid$theta, w = par_grid$omega, tr = tr[1], sr = sr[1], 
                         to = to, so = so, x = eta, y = nu, m = mu_UIP, v = tau_UIP)

postdens_wrapper <- do.call("rbind", lapply(X = seq(1:length(tr)), FUN = function(index) {
  post_dens <- posteriormix(
    theta = par_grid$theta,
    w = par_grid$omega,
    tr = tr[index],
    sr = sr[index],
    to = to,
    so = so,
    x = eta,
    y = nu,
    m = mu_UIP,
    v = tau_UIP
  )
  par_grid$tr <- tr[index]
  par_grid$sr <- sr[index]
  par_grid$rnumber <- index
  par_grid$density <- post_dens
  return(par_grid)
}))



##  ............................................................................
##  Marginal Posterior for w                                                ####

weights_m_post <- do.call("rbind", lapply(X = seq(1:length(tr)), FUN = function(index) {
  marg_p_dens <- wposteriormix(w = wseq,
                                tr = tr[index],
                                sr = sr[index],
                                to = to,
                                so = so,
                                x = eta,
                                y = nu,
                                m = mu_UIP,
                                v = tau_UIP
  )
  out <- data.frame(x = wseq, density = marg_p_dens, rnumber = rep_number[index],
                    parameter = "'Weight parameter' ~ omega", tr = tr[index], sr = sr[index])
  return(out)
}))


##  ............................................................................
##  Weights Marginal Posterior Plot                                         ####



### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Credibility Intervals                                                   ####


HPDI_weights <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  hpdi <- wHPD(level = 0.95, tr = tr[i], sr = sr[i], to = to,
                              so = so, m = mu_UIP, v = tau_UIP, x = eta, y = nu)
  out <- data.frame(y = max(weights_m_post$density)*(1 + 0.05*i),
                    lower = hpdi[1], upper = hpdi[3], rnumber = rep_number[i],
                    parameter = "'Weight parameter' ~ omega", tr = tr[i],
                    sr = sr[i], height = 0.2)
  return(out)
}))



##  ............................................................................
##  Theta Marginal Posterior                                                ####


theta_m_post <- do.call("rbind", lapply(X = seq(1:length(tr)), FUN = function(index) {
  marg_p_dens <- thetaposteriormix(  theta = thetaseq,
                                     tr = tr[index],
                                     sr = sr[index],
                                     to = to,
                                     so = so,
                                     x = eta,
                                     y = nu,
                                     m = mu_UIP,
                                     v = tau_UIP
  )
  out <- data.frame(x = thetaseq, density = marg_p_dens, rnumber = rep_number[index],
                    parameter = "'Effect size' ~ theta", tr = tr[index], sr = sr[index])
  return(out)
}))


## Posterior of effect size without using original data
theta_m_post_2 <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  pDens <- dnorm(x = thetaseq, mean = tr[i], sd = sr[i])
  out <- data.frame(x = thetaseq, density = pDens, rnumber = rep_number[i],
                    parameter = "'Effect size' ~ theta", tr = tr[i], sr = sr[i])
  return(out)
}))


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Credibility Intervals                                                   ####

HPDI_theta <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  hpd <- thetaHPD(level = 0.95, tr = tr[i], sr = sr[i], to = to,
                        so = so, m = mu_UIP, v = tau_UIP, x = eta, y = nu)
  out <- data.frame(y = max(c(theta_m_post_2$density, theta_m_post$density))*(1 + 0.06*i),
                    lower = hpd[1], upper = hpd[3], rnumber = rep_number[i],
                    parameter = "'Effect size' ~ theta", tr = tr[i],
                    sr = sr[i], height = 0.6)
  return(out)
}))

HPDI_theta_2 <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  hpd <- tr[i] + c(-1, 1)*qnorm(p = 0.975)*sr[i]
  out <- data.frame(y = max(c(theta_m_post_2$density, theta_m_post$density))*(1 + 0.05*i),
                    lower = hpd[1], upper = hpd[2], rnumber = rep_number[i],
                    parameter = "'Effect size' ~ theta", tr = tr[i],
                    sr = sr[i], height = 0.6)
  return(out)
}))
HPDI_theta_2$trFormat <- paste0("{hat(theta)[italic('r')*", HPDI_theta_2$rnumber, "] == ",
                             round(HPDI_theta_2$tr, 2), "}*',' ~ sigma[italic('r')*",
                             HPDI_theta_2$rnumber, "] == ", round(HPDI_theta_2$sr, 2))

  
#   ____________________________________________________________________________
#   Bayes Factor                                                            ####

format_bf <- function(BF, digits = "default") {
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
format_bf_vec <- Vectorize(FUN = format_bf)



rnumber <- c(1, 2, 3, 4)

bf_df <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  bf_theta_random <- bf_theta_mix(tr = tr[i], sr = sr[i], to = to, so = so,
                    x = 1, y = 1, m = mu_UIP, v = tau_UIP)
  bf_theta <- bf_theta_mix(tr = tr[i], sr = sr[i], to = to, so = so, x = 1, y = 1, 
                           m = mu_UIP, v = tau_UIP, w = 1)
  bf_omega <- bf_omega_mix(tr = tr[i], sr = sr[i], to = to, so = so,
                       x = 1, y = 1, m = mu_UIP, v = tau_UIP, w_null = 0, w_alt = 1)
  bf_random_omega_1 <- bf_omega_mix(tr = tr[i], sr = sr[i], to = to, so = so,
                             x = 1, y = 2, m = mu_UIP, v = tau_UIP, w_null = NA, w_alt = NA)
  bf_random_omega_2 <- bf_omega_mix(tr = tr[i], sr = sr[i], to = to, so = so,
                               x = 2, y = 1, m = mu_UIP, v = tau_UIP, w_null = NA, w_alt = NA)
  out <- data.frame(number = rnumber[i], tr = tr[i], sr = sr[i], bf_theta = bf_theta,
                    bf_theta_random = bf_theta_random, bf_omega = bf_omega, 
                    bf_random_omega_1 = bf_random_omega_1, bf_random_omega_2 = bf_random_omega_2)
  return(out)
}))



# Create LaTeX table for theta
dfTab_theta <- bf_df[,1:5] %>%
  mutate(bf_theta = format_bf_vec(bf_theta),
         bf_theta_random = format_bf_vec(bf_theta_random),
         tr = round(tr, 2),
         sr = round(sr, 2),
         number = as.integer(number)) %>%
  arrange(number)
xtab_theta <- xtable(dfTab_theta)
colnames(xtab_theta) <- c("",
                    "$\\hat{\\theta}_r$",
                    "$\\sigma_r$",
                    paste0("$\\mathrm{BF}_{01}\\{\\hat{\\theta}_r \\mid \\mathcal{H}_{1} : \\omega \\sim \\mathrm{Beta}(",
                           eta, ", ", nu, ")\\}$"),
                    "$\\mathrm{BF}_{01}(\\hat{\\theta}_r \\mid \\mathcal{H}_{1} : \\omega = 1)$"
                    )
align(xtab_theta) <- rep("c", length(colnames(xtab_theta)) + 1)

# Add multicolumns for effet size test and power parameter test
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\toprule'

print(xtab_theta, floating = FALSE, include.rownames = FALSE, add.to.row = addtorow,
      sanitize.text.function = function(x){x}, booktabs = TRUE, hline.after = c(0, nrow(xtab_theta)))



# Create LaTeX table for omega
dfTab_omega <- bf_df[,c(1:3,6:8)] %>%
  mutate(bf_omega = format_bf_vec(bf_omega),
         bf_random_omega_1 = format_bf_vec(bf_random_omega_1),
         bf_random_omega_2 = format_bf_vec(bf_random_omega_2),
         tr = round(tr, 2),
         sr = round(sr, 2),
         number = as.integer(number)) %>%
  arrange(number)
xtab_omega <- xtable(dfTab_omega)
colnames(xtab_omega) <- c("",
                    "$\\hat{\\theta}_r$",
                    "$\\sigma_r$",
                    paste0("$\\mathrm{BF}_{\\text{dc}}(\\hat{\\theta}_r \\mid \\mathcal{H}_d \\: \\omega = ",
                           0, ")$"),
                    paste0("$\\mathrm{BF}_{\\text{dc}}\\{\\hat{\\theta}_r \\mid \\mathcal{H}_d \\: \\omega \\sim \\mathrm{Beta}(",
                           1, ", ", 2, ")\\}$"),
                    paste0("$\\mathrm{BF}_{\\text{dc}}\\{\\hat{\\theta}_r \\mid \\mathcal{H}_d \\: \\omega \\sim \\mathrm{Beta}(",
                           2, ", ", 1, ")\\}$")
)
align(xtab_omega) <- rep("c", length(colnames(xtab_omega)) + 1)

# Add multicolumns for effet size test and power parameter test
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\toprule'

print(xtab_omega, floating = FALSE, include.rownames = FALSE, add.to.row = addtorow,
      sanitize.text.function = function(x){x}, booktabs = TRUE, hline.after = c(0, nrow(xtab_omega)))

