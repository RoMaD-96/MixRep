#   ____________________________________________________________________________
#   Libraries                                                               ####

library(dplyr)
library(knitr)
library(xtable)
library(repmix)
library(hypergeo)
library(ppRep)
library(bayesmeta)
library(ggrepel)

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

tp <- round(sum(trep / srep^2) / sum(1 / srep^2), 2)
sp <- round(sqrt(1 / sum(1 / srep^2)), 2)

tr <- c(trep, tp)
sr <- c(srep, sp)


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Mixture Prior                                                           ####

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


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Power Prior                                                             ####

# Uniform prior for alpha
x <- 1
y <- 1


# Parameter grid to compute posterior density
n_alpha <- 400
n_theta <- 400
alphaseq <- seq(0, 1, length.out = n_alpha)
thetaseq <- seq(-0.9, 0.9, length.out = n_theta)
par_grid <- expand.grid(alpha = alphaseq, theta = thetaseq)
m <- 0
v <- Inf


# Replication Number
rnumber <- c(1, 2, 3, 4)



#   ____________________________________________________________________________
#   Mixture Prior Analysis                                                  ####

##  ............................................................................
##  Theta Marginal Posterior                                                ####


theta_m_post <- do.call("rbind", lapply(X = seq(1:length(tr)), FUN = function(index) {
  marg_p_dens <- thetaposteriormix(
    theta = thetaseq,
    tr = tr[index],
    sr = sr[index],
    to = to,
    so = so,
    x = eta,
    y = nu,
    m = mu_UIP,
    v = tau_UIP
  )
  out <- data.frame(
    x = thetaseq, density = marg_p_dens, rnumber = rnumber[index],
    parameter = "'Effect size' ~ theta", tr = tr[index], sr = sr[index]
  )
  return(out)
}))


## Posterior of effect size without using original data
theta_unif <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(index) {
  pDens <- dnorm(x = thetaseq, mean = tr[index], sd = sr[index])
  out <- data.frame(
    x = thetaseq, density = pDens, rnumber = rnumber[index],
    parameter = "'Effect size' ~ theta", tr = tr[index], sr = sr[index]
  )
  return(out)
}))


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Credibility Intervals                                                   ####

HPDI_theta_mix <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(index) {
  hpd <- thetaHPD(
    level = 0.95, tr = tr[index], sr = sr[index], to = to,
    so = so, m = mu_UIP, v = tau_UIP, x = eta, y = nu
  )
  out <- data.frame(
    y = 11.1,
    lower = hpd[1], upper = hpd[3], rnumber = rnumber[index],
    parameter = "'Effect size' ~ theta", tr = tr[index],
    sr = sr[index], height = 0.6
  )
  return(out)
}))

HPDI_unif <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(index) {
  hpd <- tr[index] + c(-1, 1) * qnorm(p = 0.975) * sr[index]
  out <- data.frame(
    y = 10.4,
    lower = hpd[1], upper = hpd[2], rnumber = rnumber[index],
    parameter = "'Effect size' ~ theta", tr = tr[index],
    sr = sr[index], height = 0.6
  )
  return(out)
}))



#   ____________________________________________________________________________
#   Power Prior Analysis                                                    ####

##  ............................................................................
##  Theta Marginal Posterior                                                ####

marg_theta_dens <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(index) {
  pp_joint_post <- postPPtheta(
    theta = thetaseq, tr = tr[index], sr = sr[index], to = to,
    so = so, x = x, y = y, m = m, v = v
  )
  out <- data.frame(
    x = thetaseq, density = pp_joint_post, rnumber = rnumber[index],
    parameter = "'Effect size' ~ theta", tr = tr[index], sr = sr[index]
  )
  return(out)
}))


# # Posterior of effect size without using original data
# marg_theta_dens2 <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
#   pp_joint_post <- dnorm(x = thetaseq, mean = tr[i], sd = sr[i])
#   out <- data.frame(
#     x = thetaseq, density = pp_joint_post, rnumber = rnumber[i],
#     parameter = "'Effect size' ~ theta", tr = tr[i], sr = sr[i]
#   )
#   return(out)
# }))



### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Credibility Intervals                                                   ####

HPDI_theta_pp <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(index) {
  hpd <- postPPthetaHPD(
    level = 0.95, tr = tr[index], sr = sr[index], to = to,
    so = so, x = x, y = y, m = m, v = v
  )
  out <- data.frame(
    y = 11.8,
    lower = hpd[1], upper = hpd[2], rnumber = rnumber[index],
    parameter = "'Effect size' ~ theta", tr = tr[index],
    sr = sr[index], height = 0.6
  )
  return(out)
}))





#   ____________________________________________________________________________
#   Hierarchical models                                                     ####


# Define labels for the replication sites
rep_names <- c(
  "University of Toronto", "Montana State University",
  "Ashland University", "Pooled"
)

##  ............................................................................
##  Theta Marginal Posterior                                                ####

theta_m_post_hm <- do.call("rbind", lapply(X = seq(1:length(tr)), FUN = function(index) {
  y_pair <- c(to, tr[index])
  sigma_pair <- c(so, sr[index])

  # Bayesian two-study "meta-analysis" (replicated and original study)
  result_pair <- bayesmeta(
    y = y_pair,
    sigma = sigma_pair,
    labels = c("Original", rep_names[index]),
    tau.prior = function(t) dhalfnormal(t, scale = 0.1) # Half-normal prior for heterogeneity
  )
  theta_post <- result_pair$dposterior(mu = thetaseq)
  out <- data.frame(
    x = thetaseq, density = theta_post, rnumber = rnumber[index],
    parameter = "'Effect size' ~ theta", tr = tr[index], sr = sr[index]
  )
  return(out)
}))


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Credibility Intervals                                                   ####

HPDI_theta_hm <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(index) {
  y_pair <- c(to, tr[index])
  sigma_pair <- c(so, sr[index])

  result_pair <- bayesmeta(
    y = y_pair,
    sigma = sigma_pair,
    labels = c("Original", rep_names[index]),
    tau.prior = function(t) dhalfnormal(t, scale = 0.1) # Half-normal prior for heterogeneity
  )

  hpd <- result_pair$post.interval(mu.level = 0.95)

  out <- data.frame(
    y = 12.5,
    lower = hpd[[1]], upper = hpd[[2]], rnumber = rnumber[index],
    parameter = "'Effect size' ~ theta", tr = tr[index],
    sr = sr[index], height = 0.6
  )
  return(out)
}))



#   ____________________________________________________________________________
#   Plots                                                                   ####

# Densities
mix_theta <- theta_m_post %>% mutate(prior = "Mixture")
hm_theta <- theta_m_post_hm %>% mutate(prior = "Hierarchical")
unif_theta <- theta_unif %>% mutate(prior = "Uniform")
pp_theta <- marg_theta_dens %>% mutate(prior = "Power")

dens_orig <- bind_rows(mix_theta, pp_theta, hm_theta)

dens_orig$rnumber <- factor(dens_orig$rnumber,
  levels = c("1", "2", "3", "4"),
  labels = c(
    expression(" " ~ hat(theta)[r * 1] == 0.29 ~ ", " ~ sigma[r * 1] == 0.11),
    expression(" " ~ hat(theta)[r * 2] == 0.25 ~ ", " ~ sigma[r * 2] == 0.09),
    expression(" " ~ hat(theta)[r * 3] == -0.18 ~ ", " ~ sigma[r * 3] == 0.11),
    expression(" " ~ hat(theta)[r * p] == 0.14 ~ ", " ~ sigma[r * p] == 0.06)
  )
)

unif_theta$rnumber <- factor(unif_theta$rnumber,
  levels = c("1", "2", "3", "4"),
  labels = c(
    expression(" " ~ hat(theta)[r * 1] == 0.29 ~ ", " ~ sigma[r * 1] == 0.11),
    expression(" " ~ hat(theta)[r * 2] == 0.25 ~ ", " ~ sigma[r * 2] == 0.09),
    expression(" " ~ hat(theta)[r * 3] == -0.18 ~ ", " ~ sigma[r * 3] == 0.11),
    expression(" " ~ hat(theta)[r * p] == 0.14 ~ ", " ~ sigma[r * p] == 0.06)
  )
)

# HPDIs
HPDI_mix <- HPDI_theta_mix %>% mutate(prior = "Mixture")
HPDI_hm <- HPDI_theta_hm %>% mutate(prior = "Hierarchical")
HPDI_unif <- HPDI_unif %>% mutate(prior = "Uniform")
HPDI_pp <- HPDI_theta_pp %>% mutate(prior = "Power")


hpdi_orig <- bind_rows(HPDI_mix, HPDI_pp, HPDI_hm)

hpdi_orig$rnumber <- factor(hpdi_orig$rnumber,
  levels = c("1", "2", "3", "4"),
  labels = c(
    expression(" " ~ hat(theta)[r * 1] == 0.29 ~ ", " ~ sigma[r * 1] == 0.11),
    expression(" " ~ hat(theta)[r * 2] == 0.25 ~ ", " ~ sigma[r * 2] == 0.09),
    expression(" " ~ hat(theta)[r * 3] == -0.18 ~ ", " ~ sigma[r * 3] == 0.11),
    expression(" " ~ hat(theta)[r * p] == 0.14 ~ ", " ~ sigma[r * p] == 0.06)
  )
)

HPDI_unif$rnumber <- factor(HPDI_unif$rnumber,
  levels = c("1", "2", "3", "4"),
  labels = c(
    expression(" " ~ hat(theta)[r * 1] == 0.29 ~ ", " ~ sigma[r * 1] == 0.11),
    expression(" " ~ hat(theta)[r * 2] == 0.25 ~ ", " ~ sigma[r * 2] == 0.09),
    expression(" " ~ hat(theta)[r * 3] == -0.18 ~ ", " ~ sigma[r * 3] == 0.11),
    expression(" " ~ hat(theta)[r * p] == 0.14 ~ ", " ~ sigma[r * p] == 0.06)
  )
)


hpdi_all <- bind_rows(
  HPDI_unif,
  hpdi_orig
)

# Create unif_len dataframe with uniform prior interval lengths
unif_len <- HPDI_unif %>%
  mutate(
    len_unif = upper - lower,
    rnumber = as.character(rnumber)
  ) %>%
  select(rnumber, len_unif)

# Each method's length and percentage reduction
reduct_df <- hpdi_all %>%
  filter(prior != "Uniform") %>%
  mutate(
    len = upper - lower
  ) %>%
  left_join(unif_len, by = "rnumber") %>%
  mutate(
    pct_red = (len - len_unif) / len_unif * 100,
    x_pos   = upper + 0.1
  )

comp_plot <- ggplot() +
  geom_errorbarh(
    data = hpdi_orig,
    aes(
      xmin = lower, xmax = upper, y = y, color = factor(prior),
      height = height
    ), alpha = 0.8, linewidth = 1.0
  ) +
  geom_errorbarh(
    data = HPDI_unif,
    aes(
      xmin = lower, xmax = upper, y = y, color = factor(prior),
      height = height
    ), alpha = 0.7, linewidth = 1.0, linetype = "22"
  ) +
  geom_line(
    data = unif_theta, aes(x = x, y = density, color = factor(prior)),
    lty = 22, alpha = 0.9, linewidth = 1.2
  ) +
  geom_line(
    data = dens_orig, aes(x = x, y = density, color = factor(prior)),
    lty = 1, alpha = 0.9, linewidth = 1.2
  ) +
  geom_text(
    data = reduct_df,
    aes(
      x     = x_pos,
      y     = y,
      label = paste0(round(pct_red, 1), "%"),
      color = prior # if you want them colored by method
    ),
    hjust = 0, # left-align at x_pos
    vjust = 0.5, # center-vertically on the y value
    size = 5,
    show.legend = FALSE
  ) +
  facet_wrap(~rnumber,
    ncol = 4,
    labeller = label_parsed
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      Hierarchical = "#332288",
      Mixture      = "#44AA99",
      Power        = "#E69F00",
      Uniform      = "black"
    ),
    breaks = c("Hierarchical", "Mixture", "Power", "Uniform"),
    labels = c(
      "Hierarchical Model",
      "Mixture Prior",
      "Power Prior",
      "Uniform Prior"
    )
  ) +
  labs(
    x = expression("Effect Size" ~ theta),
    y = "Density"
  ) +
  theme_bw() +
  theme(
    strip.text.x = element_text(size = 18),
    legend.position = "top",
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 22),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16)
  )

print(comp_plot)

ggsave(
  filename = "comp_plot.pdf", path = "Plots",
  plot = comp_plot,
  width = 17, height = 7.5, device = "pdf", dpi = 500, useDingbats = FALSE
)
