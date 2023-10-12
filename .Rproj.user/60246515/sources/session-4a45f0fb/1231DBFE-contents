#   ____________________________________________________________________________
#   Libraries                                                               ####

library(ggplot2)
library(ggpubr)
library(colorspace)

source("RepMixFun.R")

#   ____________________________________________________________________________
#   Random Weights                                                          ####


##  ............................................................................
##  Parameter Setting                                                       ####

# Original and Replicated Studies
to <- 0.21
so <- 0.05
tr <- c(0.09, 0.21, 0.44)
sr <- c(0.045, 0.06, 0.04)
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
##  Contour plot of joint posterior                                         ####

postdens_wrapper$rep_setting <-paste0( "{hat(theta)[italic('r')*",
                                       postdens_wrapper$rep_number,
                                       "] == ",
                                       round(postdens_wrapper$tr, 2),
                                       "}*',' ~ sigma[italic('r')*",
                                       postdens_wrapper$rep_number,
                                       "] == ",
                                       round(postdens_wrapper$sr, 2)
)

plot_joint <- ggplot(data = postdens_wrapper, aes(x = theta, y = omega, fill = density)) +
  facet_wrap(~ rep_setting, labeller = label_parsed) +
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(z = density), breaks = seq(0, 24, 3),  alpha = 0.35, size = 0.5) +
  scale_fill_continuous_sequential(palette = "Blues 3", rev = TRUE) +
  labs(
    subtitle = "Contour Plot considering three replications scenario",
    x = bquote("Effect Size" ~ theta),
    y = bquote("Weight Parameter" ~ omega),
    fill = "Posterior \n Density"
  ) +
  guides(fill = guide_colorbar(barheight = 12, barwidth = 0.9, title.position = "top")) +
  theme_minimal()

print(plot_joint)



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
  out <- data.frame(x = wseq, density = marg_p_dens, rep_exp = rep_number[index],
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
                    parameter = "'Weight Parameter' ~ omega", tr = tr[i],
                    sr = sr[i], height = 0.2)
  return(out)
}))


# plot_m_weight <- ggplot(data=weights_m_post, aes(x=x, y=density, group=rep_exp, color=factor(rep_exp))) +
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



plot_weights_m_hpd <- ggplot() +
  geom_errorbarh(data = HPDI_weights,
                 aes(xmin = lower, xmax = upper, y = y*1.05, color = factor(rep_number),
                     height = height), alpha = 0.8, size = 1.2) +
  geom_line(data=weights_m_post, aes(x=x, y=density, group=rep_exp, color=factor(rep_exp)),
            lty = 1, alpha = 0.9, size = 1.2) +
  scale_color_manual(
    values = c("1" = "#E69F00", "2" = "#009E20", "3" = "#0072B2"),
    labels = c(expression(" "~hat(theta)[r*1] == 0.09 ~ ", " ~ sigma[r*1] == 0.04),
               expression(" "~hat(theta)[r * 2] == 0.21 ~ ", " ~ sigma[r*2] == 0.06),
               expression(" "~hat(theta)[r * 3] == 0.44 ~ ", " ~ sigma[r*3] == 0.04))) +
  labs(
    subtitle = "Marginal Posterior Densities of Weight Parameter",
    x = expression(omega~" Values"),
    y = "Density"
  ) +
  theme_minimal() +
  guides(color=guide_legend(title="Replicated Experiment")) 

print(plot_weights_m_hpd)


#Flat prior and linear combination of omegas

##  ............................................................................
##  Theta Marginal Posterior Plot                                                 ####

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
                    parameter = "'Effect parameter' ~ theta", tr = tr[index], sr = sr[index])
  return(out)
}))


## Posterior of effect size without using original data
thetaplotDF2 <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  pDens <- dnorm(x = thetaseq, mean = tr[i], sd = sr[i])
  out <- data.frame(x = thetaseq, density = pDens, rep_number = rep_number[i],
                    parameter = "'Effect size' ~ theta", tr = tr[i], sr = sr[i])
  return(out)
}))


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Credibility Intervals                                                   ####

thetaHPD <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  hpd <- HPDI_post_m_theta(level = 0.95, tr = tr[i], sr = sr[i], to = to,
                        so = so, x = alpha, y = beta, null = null, priorsd = priorsd)
  out <- data.frame(y = max(c(thetaplotDF2$density, theta_m_post$density))*(1 + 0.06*i),
                    lower = hpd[1], upper = hpd[2], rep_number = rep_number[i],
                    parameter = "'Effect size' ~ theta", tr = tr[i],
                    sr = sr[i], height = 0.6)
  return(out)
}))

theta2HPD <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  hpd <- tr[i] + c(-1, 1)*qnorm(p = 0.975)*sr[i]
  out <- data.frame(y = max(c(thetaplotDF2$density, theta_m_post$density))*(1 + 0.05*i),
                    lower = hpd[1], upper = hpd[2], rep_number = rep_number[i],
                    parameter = "'Effect size' ~ theta", tr = tr[i],
                    sr = sr[i], height = 0.6)
  return(out)
}))
theta2HPD$trFormat <- paste0("{hat(theta)[italic('r')*", theta2HPD$rnumber, "] == ",
                             round(theta2HPD$tr, 2), "}*',' ~ sigma[italic('r')*",
                             theta2HPD$rnumber, "] == ", round(theta2HPD$sr, 2))




plot_m_theta <- ggplot() +
  geom_errorbarh(data = thetaHPD,
                 aes(xmin = lower, xmax = upper, y = y*1.05, color = factor(rep_number),
                     height = height), alpha = 0.8, size = 1.0) +
  geom_errorbarh(data = theta2HPD,
                 aes(xmin = lower, xmax = upper, y = y*1.05, color = factor(rep_number),
                     height = height), alpha = 0.7, linetype = "22", size = 1.0) +
  geom_line(data = thetaplotDF2, aes(x = x, y = density, color = factor(rep_number)),
            lty = 2, alpha = 0.5, size = 1.0) +
  geom_line(data = theta_m_post, aes(x = x, y = density, color = factor(rep_number)),
            alpha = 0.9, size = 1.0) +
  scale_color_manual(
    values = c("1" = "#E69F00", "2" = "#009E20", "3" = "#0072B2"),
    labels = c(expression(" "~hat(theta)[r*1] == 0.09 ~ ", " ~ sigma[r*1] == 0.04),
               expression(" "~hat(theta)[r * 2] == 0.21 ~ ", " ~ sigma[r*2] == 0.06),
               expression(" "~hat(theta)[r * 3] == 0.44 ~ ", " ~ sigma[r*3] == 0.04))) +
  labs(
    subtitle = "Marginal Posterior Densities of Effect Parameter",
    x = expression(theta~" Values"),
    y = "Density"
  ) +
  theme_minimal() +
  guides(color=guide_legend(title="Replicated Experiment")) 


print(plot_m_theta)

