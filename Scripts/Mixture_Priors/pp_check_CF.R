#   ____________________________________________________________________________
#   Libraries                                                               ####

library(repmix)
library(ggplot2)
library(dplyr)

set.seed(2333)

#   ____________________________________________________________________________
#   Posterior predictive functions                                          ####

# Posterior parameters
normal2mixposterior <- function(y, se, m1, v1, m2, v2, w) {
  ## updated weight
  wpost <- 1 /
    (1 +
      ## priors odds
      (1 - w) / w *
        ## Bayes factor
        stats::dnorm(x = y, mean = m2, sd = sqrt(se^2 + v2)) /
        stats::dnorm(x = y, mean = m1, sd = sqrt(se^2 + v1))

    )
  ## updated component means and variances
  v1post <- 1 / (1 / v1 + 1 / se^2)
  m1post <- (m1 / v1 + y / se^2) * v1post
  v2post <- 1 / (1 / v2 + 1 / se^2)
  m2post <- (m2 / v2 + y / se^2) * v2post
  res <- list(
    "w" = wpost, "m1" = m1post, "v1" = v1post, "m2" = m2post,
    "v2" = v2post
  )
  return(res)
}

# Closed-form posterior predictive distribution
posterior_predictive_rep <- function(theta, tr, sr, to, so,
                                     w = 0.5, m = 0, v = 1) {
  # Posterior parameters
  post_mix_pars <- normal2mixposterior(y = tr, sr, to, so, mu_UIP, tau_UIP, w)


  w_post <- post_mix_pars$w
  m1_post <- post_mix_pars$m1
  v1_post <- post_mix_pars$v1
  m2_post <- post_mix_pars$m2
  v2_post <- post_mix_pars$v2

  post_comp_1 <- dnorm(theta, m1_post, sqrt(v1_post + sr^2))
  post_comp_2 <- dnorm(theta, m2_post, sqrt(v2_post + sr^2))

  # Posterior predictive replicated data
  tr_rep <- w_post * post_comp_1 + (1 - w_post) * post_comp_2

  return(tr_rep)
}


# Posterior predictive mean
posterior_predictive_mean <- function(tr, sr, to, so, w = 0.5) {
  post_mix_pars <- normal2mixposterior(
    y = tr, se = sr, m1 = to, v1 = so,
    m2 = mu_UIP, v2 = tau_UIP, w = w
  )

  w_post <- post_mix_pars$w
  m1_post <- post_mix_pars$m1
  m2_post <- post_mix_pars$m2

  # Mean posterior predictive replicated data
  mean_post_rep <- w_post * m1_post + (1 - w_post) * m2_post

  return(mean_post_rep)
}


# Posterior predictive CDF
posterior_predictive_cdf <- function(x, tr, sr, to, so, w = 0.5, m = 0, v = 1) {
  post <- normal2mixposterior(y = tr, se = sr, m1 = to, v1 = so, m2 = m, v2 = v, w = w)
  w_post <- post$w
  m1_post <- post$m1
  v1_post <- post$v1
  m2_post <- post$m2
  v2_post <- post$v2
  sd1_post <- sqrt(v1_post + sr^2)
  sd2_post <- sqrt(v2_post + sr^2)
  w_post * pnorm(x, mean = m1_post, sd = sd1_post) +
    (1 - w_post) * pnorm(x, mean = m2_post, sd = sd2_post)
}

#   ____________________________________________________________________________
#   Data                                                                    ####

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

# Mean and Variance Unit Informative Prior
mu_UIP <- 0
tau_UIP <- 2


w <- 0.5
thetaseq <- seq(-0.9, 0.9, length.out = 2500)

# Replication Number
rep_number <- c(1, 2, 3, 4)


#   ____________________________________________________________________________
#   Plot posterior predictive                                               ####

densities_ppd <- list()
for (i in 1:length(tr)) {
  densities_ppd[[i]] <- posterior_predictive_rep(
    theta = thetaseq, tr = tr[i], sr = sr[i], to = to, so = so,
    m = mu_UIP, v = tau_UIP, w = w
  )
}

df <- list()

for (i in 1:length(tr)) {
  # Turn the results into a data frame
  df[[i]] <- data.frame(
    theta = thetaseq,
    density = as.vector(densities_ppd[[i]])
  )
  df[[i]]$replication <- paste0("Replication ", i)
}

df_densities <- do.call(rbind, df)

df_densities$replication <- factor(df_densities$replication,
  levels = c("Replication 1", "Replication 2", "Replication 3", "Replication 4"),
  labels = c(
    expression(" " ~ hat(theta)[r * 1] == 0.29 ~ ", " ~ sigma[r * 1] == 0.11),
    expression(" " ~ hat(theta)[r * 2] == 0.25 ~ ", " ~ sigma[r * 2] == 0.09),
    expression(" " ~ hat(theta)[r * 3] == -0.18 ~ ", " ~ sigma[r * 3] == 0.11),
    expression(" " ~ hat(theta)[r * p] == 0.14 ~ ", " ~ sigma[r * p] == 0.06)
  )
)


# Posterior predictive mean for each replication
pp_means <- vapply(seq_along(tr), function(i) posterior_predictive_mean(tr[i], sr[i], to, so, w), numeric(1))

# Density height
dens_at_x <- function(theta_grid, dens_vals, x) {
  approx(x = theta_grid, y = dens_vals, xout = x, rule = 2)$y
}

lines_list <- vector("list", length(tr))
for (i in seq_along(tr)) {
  dens_i <- densities_ppd[[i]]
  y_tr <- dens_at_x(thetaseq, dens_i, tr[i])
  y_mean <- dens_at_x(thetaseq, dens_i, pp_means[i])

  lines_list[[i]] <- data.frame(
    replication = paste0("Replication ", i),
    x = c(tr[i], pp_means[i]),
    yend = c(y_tr, y_mean),
    type = c("Observed", "Predictive mean")
  )
}

lines_df <- do.call(rbind, lines_list)

# Facet labels
lines_df$replication <- factor(
  lines_df$replication,
  levels = c("Replication 1", "Replication 2", "Replication 3", "Replication 4"),
  labels = c(
    expression(" " ~ hat(theta)[r * 1] == 0.29 ~ ", " ~ sigma[r * 1] == 0.11),
    expression(" " ~ hat(theta)[r * 2] == 0.25 ~ ", " ~ sigma[r * 2] == 0.09),
    expression(" " ~ hat(theta)[r * 3] == -0.18 ~ ", " ~ sigma[r * 3] == 0.11),
    expression(" " ~ hat(theta)[r * p] == 0.14 ~ ", " ~ sigma[r * p] == 0.06)
  )
)



#   ____________________________________________________________________________
#   Posterior predictive p-values                                           ####

ppp_1 <- 1 - posterior_predictive_cdf(tr[1], tr[1], sr[1], to, so, w)
ppp_2 <- 1 - posterior_predictive_cdf(tr[2], tr[2], sr[2], to, so, w)
ppp_3 <- 1 - posterior_predictive_cdf(tr[3], tr[3], sr[3], to, so, w)
ppp_4 <- 1 - posterior_predictive_cdf(tr[4], tr[4], sr[4], to, so, w)

round(c(ppp_1, ppp_2, ppp_3, ppp_4), 3)


#   ____________________________________________________________________________
#   Plot with colored density                                               ####

ppp_vals <- round(c(ppp_1, ppp_2, ppp_3, ppp_4), 3)

ppp_facets <- data.frame(
  replication = factor(
    c("Replication 1", "Replication 2", "Replication 3", "Replication 4"),
    levels = c("Replication 1", "Replication 2", "Replication 3", "Replication 4"),
    labels = c(
      expression(" " ~ hat(theta)[r * 1] == 0.29 ~ ", " ~ sigma[r * 1] == 0.11),
      expression(" " ~ hat(theta)[r * 2] == 0.25 ~ ", " ~ sigma[r * 2] == 0.09),
      expression(" " ~ hat(theta)[r * 3] == -0.18 ~ ", " ~ sigma[r * 3] == 0.11),
      expression(" " ~ hat(theta)[r * p] == 0.14 ~ ", " ~ sigma[r * p] == 0.06)
    )
  ),
  x = -Inf, # top-left corner
  y = 5,
  label = paste0("PPP = ", sprintf("%.3f", ppp_vals))
)

colored_density <- list()
for (i in 1:length(tr)) {
  dens_vals <- densities_ppd[[i]]

  # Color density where trep > tr
  upper_idx <- which(thetaseq >= tr[i])
  if (length(upper_idx) > 0) {
    colored_density[[i]] <- data.frame(
      theta = c(tr[i], thetaseq[upper_idx], max(thetaseq[upper_idx])),
      density = c(0, dens_vals[upper_idx], 0),
      replication = paste0("Replication ", i)
    )
  } else {
    colored_density[[i]] <- data.frame()
  }
}

colored_density_df <- do.call(rbind, colored_density)

# Factor levels
colored_density_df$replication <- factor(
  colored_density_df$replication,
  levels = c("Replication 1", "Replication 2", "Replication 3", "Replication 4"),
  labels = c(
    expression(" " ~ hat(theta)[r * 1] == 0.29 ~ ", " ~ sigma[r * 1] == 0.11),
    expression(" " ~ hat(theta)[r * 2] == 0.25 ~ ", " ~ sigma[r * 2] == 0.09),
    expression(" " ~ hat(theta)[r * 3] == -0.18 ~ ", " ~ sigma[r * 3] == 0.11),
    expression(" " ~ hat(theta)[r * p] == 0.14 ~ ", " ~ sigma[r * p] == 0.06)
  )
)

# Labels
ppp_labels <- data.frame(
  replication = factor(
    c("Replication 1", "Replication 2", "Replication 3", "Replication 4"),
    levels = c("Replication 1", "Replication 2", "Replication 3", "Replication 4"),
    labels = c(
      expression(" " ~ hat(theta)[r * 1] == 0.29 ~ ", " ~ sigma[r * 1] == 0.11),
      expression(" " ~ hat(theta)[r * 2] == 0.25 ~ ", " ~ sigma[r * 2] == 0.09),
      expression(" " ~ hat(theta)[r * 3] == -0.18 ~ ", " ~ sigma[r * 3] == 0.11),
      expression(" " ~ hat(theta)[r * p] == 0.14 ~ ", " ~ sigma[r * p] == 0.06)
    )
  ),
  x = 0.7,
  y = max(sapply(densities_ppd, max)) * 0.9,
  label = paste0("p = ", round(c(ppp_1, ppp_2, ppp_3, ppp_4), 3))
)


plot_ppd_color <- ggplot() +
  geom_polygon(
    data = colored_density_df,
    aes(x = theta, y = density),
    fill = "grey40",
    alpha = 0.3
  ) +
  geom_text(
    data = ppp_facets,
    aes(x = x, y = y, label = label),
    hjust = -0.1, vjust = 1.1,
    size = 5.5,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = df_densities,
    aes(x = theta, y = density),
    linewidth = 1
  ) +
  geom_segment(
    data = lines_df,
    aes(x = x, xend = x, y = 0, yend = yend, color = type, linetype = type),
    linewidth = 1
  ) +
  scale_color_manual(
    values = c("Observed" = "firebrick2", "Predictive mean" = "black"),
    labels = c(
      "Observed" = "Observed effect size estimate",
      "Predictive mean" = "Predictive mean"
    )
  ) +
  scale_linetype_manual(
    values = c("Observed" = "twodash", "Predictive mean" = "dotted"),
    labels = c(
      "Observed" = "Observed effect size estimate",
      "Predictive mean" = "Predictive mean"
    )
  ) +
  facet_wrap(~replication, labeller = label_parsed, ncol = 4) +
  labs(
    x = expression("Effect size estimate" ~ hat(theta)["r,"]["new"]),
    y = "Density",
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text.x = element_text(size = 18),
    legend.position = "top",
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 22),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    plot.caption = element_text(size = 10, hjust = 0.5)
  )

print(plot_ppd_color)

ggsave(
  filename = "plot_ppd.pdf", path = "Plots/Mixture_Prior",
  plot = plot_ppd_color,
  width = 17, height = 7.5, device = "pdf", dpi = 500, useDingbats = FALSE
)
