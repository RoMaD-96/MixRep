#   ____________________________________________________________________________
#   Libraries                                                               ####

library(ggplot2)
library(ggthemes)
library(ggpubr)
library(colorspace)

source("Scripts/RepMixFun.R")

#   ____________________________________________________________________________
#   Random Weights                                                          ####


##  ............................................................................
##  Parameter Setting                                                       ####

# Original and Replicated Studies
to <- 0.21
so <- 0.05
tr <- c(0.09, 0.21, 0.44)
sr <- c(0.05, 0.06, 0.04)
null <- 0
priorsd <- 2

# Parameter Grid
n_weights <- 300
n_theta <- 300
wseq <- seq(0, 1, length.out = n_weights)
thetaseq <- seq(-0.2, 0.6, length.out = n_theta)
par_grid <- expand.grid(omega = wseq, theta = thetaseq)

# Uniform Prior 
alpha <- 1
beta <- 1

# Replication Number
rep_number <- c(1,2,3)


##  ............................................................................
##  Joint Posterior                                                         ####

postdens <- rmapPost(theta = par_grid$theta, w = par_grid$omega, tr = tr[1], sr = sr[1], to = to, so = so,
                     null = null, priorsd = priorsd, x = alpha, y = beta)

postdens_wrapper <- do.call("rbind", lapply(X = seq(1:length(tr)), FUN = function(index) {
  post_dens <- rmapPost(
    theta = par_grid$theta,
    w = par_grid$omega,
    tr = tr[index],
    sr = sr[index],
    to = to,
    so = so,
    null = null,
    priorsd = priorsd,
    x = alpha,
    y = beta
  )
  par_grid$tr <- tr[index]
  par_grid$sr <- sr[index]
  par_grid$rep_number <- index
  par_grid$density <- post_dens
  return(par_grid)
}))



##  ............................................................................
##  Marginal Posterior for w                                                ####

weights_m_post <- do.call("rbind", lapply(X = seq(1:length(tr)), FUN = function(index) {
  marg_p_dens <- m_post_weights(w = wseq,
                                tr = tr[index],
                                sr = sr[index],
                                to = to,
                                so = so,
                                null = null,
                                priorsd = priorsd,
                                x = alpha,
                                y = beta
  )
  out <- data.frame(x = wseq, density = marg_p_dens, rep_number = rep_number[index],
                    parameter = "'Weight parameter' ~ omega", tr = tr[index], sr = sr[index])
  return(out)
}))


##  ............................................................................
##  Weights Marginal Posterior Plot                                         ####



### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Credibility Intervals                                                   ####

HPDI_weights <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  hpdi <- HPDI_post_m_weights(level = 0.95, tr = tr[i], sr = sr[i], to = to,
                              so = so, x = alpha, y = beta, null = null, priorsd = priorsd)
  out <- data.frame(y = max(weights_m_post$density)*(1 + 0.05*i),
                    lower = hpdi[1], upper = hpdi[2], rep_number = rep_number[i],
                    parameter = "'Weight parameter' ~ omega", tr = tr[i],
                    sr = sr[i], height = 0.2)
  return(out)
}))


# plot_m_weight <- ggplot(data=weights_m_post, aes(x=x, y=density, group=rep_number, color=factor(rep_number))) +
#   geom_line(size = 1) +
#   scale_color_manual(
#     values = c("1" = "#E69F00", "2" = "#009E20", "3" = "#0072B2"),
#     labels = c(expression(" "~hat(theta)[r*1] == 0.09 ~ ", " ~ sigma[r*1] == 0.04),
#                expression(" "~hat(theta)[r * 2] == 0.21 ~ ", " ~ sigma[r*2] == 0.06),
#                expression(" "~hat(theta)[r * 3] == 0.44 ~ ", " ~ sigma[r*3] == 0.04)))+
#   labs(
#     subtitle = "Marginal Posterior Densities of Weight Parameter",
#     x = expression(omega~" Values"),
#     y = "Density"
#   ) +
#   theme_minimal() +
#   guides(color=guide_legend(title="Replicated Experiment")) 
# 
# print(plot_m_weight)





##  ............................................................................
##  Theta Marginal Posterior                                                ####

theta_m_post <- do.call("rbind", lapply(X = seq(1:length(tr)), FUN = function(index) {
  marg_p_dens <- m_post_theta(  theta = thetaseq,
                                tr = tr[index],
                                sr = sr[index],
                                to = to,
                                so = so,
                                null = null,
                                priorsd = priorsd,
                                x = alpha,
                                y = beta
  )
  out <- data.frame(x = thetaseq, density = marg_p_dens, rep_number = rep_number[index],
                    parameter = "'Effect size' ~ theta", tr = tr[index], sr = sr[index])
  return(out)
}))


## Posterior of effect size without using original data
theta_m_post_2 <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  pDens <- dnorm(x = thetaseq, mean = tr[i], sd = sr[i])
  out <- data.frame(x = thetaseq, density = pDens, rep_number = rep_number[i],
                    parameter = "'Effect size' ~ theta", tr = tr[i], sr = sr[i])
  return(out)
}))


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Credibility Intervals                                                   ####

HPDI_theta <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  hpd <- HPDI_post_m_theta(level = 0.95, tr = tr[i], sr = sr[i], to = to,
                        so = so, x = alpha, y = beta, null = null, priorsd = priorsd)
  out <- data.frame(y = max(c(theta_m_post_2$density, theta_m_post$density))*(1 + 0.06*i),
                    lower = hpd[1], upper = hpd[2], rep_number = rep_number[i],
                    parameter = "'Effect size' ~ theta", tr = tr[i],
                    sr = sr[i], height = 0.6)
  return(out)
}))

HPDI_theta_2 <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  hpd <- tr[i] + c(-1, 1)*qnorm(p = 0.975)*sr[i]
  out <- data.frame(y = max(c(theta_m_post_2$density, theta_m_post$density))*(1 + 0.05*i),
                    lower = hpd[1], upper = hpd[2], rep_number = rep_number[i],
                    parameter = "'Effect size' ~ theta", tr = tr[i],
                    sr = sr[i], height = 0.6)
  return(out)
}))
HPDI_theta_2$trFormat <- paste0("{hat(theta)[italic('r')*", HPDI_theta_2$rnumber, "] == ",
                             round(HPDI_theta_2$tr, 2), "}*',' ~ sigma[italic('r')*",
                             HPDI_theta_2$rnumber, "] == ", round(HPDI_theta_2$sr, 2))

  
#   ____________________________________________________________________________
#   Bayes Factor                                                            ####

rnumber <- c(1, 2, 3)

bf_df <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  bf_theta <- bf_theta_mix(tr = tr[i], sr = sr[i], to = to, so = so,
                    x = 1, y = 1, null = null, priorsd = priorsd)
  bf_theta_random <- bf_theta_mix(tr = tr[i], sr = sr[i], to = to, so = so, x = 1, y = 1, 
                   null = null,priorsd = priorsd, w = 1)
  bf_omega <- bf_omega_mix(tr = tr[i], sr = sr[i], to = to, so = so,
                       x = 1, y = 1, null = null, priorsd = priorsd, w_null = 0, w_alt = 1)
  bf_random_omega_1 <- bf_omega_mix(tr = tr[i], sr = sr[i], to = to, so = so,
                             x = 1, y = 2, null = null, priorsd = priorsd, w_null = NA, w_alt = NA)
  bf_random_omega_2 <- bf_omega_mix(tr = tr[i], sr = sr[i], to = to, so = so,
                               x = 2, y = 1, null = null, priorsd = priorsd, w_null = NA, w_alt = NA)
  out <- data.frame(number = rnumber[i], tr = tr[i], sr = sr[i], bf_theta = bf_theta,
                    bf_theta_random = bf_theta_random, bf_omega = bf_omega, 
                    bf_random_omega_1 = bf_random_omega_1, bf_random_omega_2 = bf_random_omega_2)
  return(out)
}))
