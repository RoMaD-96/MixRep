#   ____________________________________________________________________________
#   Libraries                                                               ####

library(ggplot2)
library(ggpubr)
library(colorspace)
library(spatstat)
library(repmix)
library(dplyr)
library(RColorBrewer)

source("Scripts/Mixture_Priors/RepMixFun_BF.R")


#   ____________________________________________________________________________
#   Fixed Weights                                                           ####


##  ............................................................................
##  Parameter Setting                                                       ####

load("credentials_data.RData")


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

# Parameter Grid
n_weights <- 300
n_theta <- 300
wseq <- seq(0, 1, by = 0.1)
thetaseq <- seq(-0.9, 0.9, length.out = 2500)
par_grid <- expand.grid(omega = wseq, theta = thetaseq)




#   ____________________________________________________________________________
#   Posterior with varying weights using ggplot 2                           ####


cols <- hcl.colors(n = length(wseq), palette = "viridis", alpha = 0.9, rev = TRUE)


# Create a data frame for ggplot
densities <- list()
for (i in 1:length(tr)) {
  densities[[i]] <- sapply(X = wseq, FUN = function(w) {
    thetaposteriormix(
      theta = thetaseq, tr = tr[i], sr = sr[i], to = to, so = so,
      m = mu_UIP, v = tau_UIP, w = w
    )
  })
}


df <- list()
additional_lines <- list()

for (i in 1:length(tr)) {
  # Turn the results into a data frame
  df[[i]] <- data.frame(
    theta = rep(thetaseq, times = length(wseq)),
    density = as.vector(densities[[i]]),
    w = factor(rep(wseq, each = length(thetaseq)))
  )

  # Add the additional lines
  df[[i]]$likelihood <- dnorm(x = df[[i]]$theta, mean = tr[i], sd = sr[i])
  df[[i]]$prior_original <- dnorm(x = df[[i]]$theta, mean = to, sd = so)
  df[[i]]$prior_robust <- dnorm(x = df[[i]]$theta, mean = mu_UIP, sd = sqrt(tau_UIP))
  df[[i]]$replication <- paste0("Replication ", i)

  # Create a separate data frame for the additional lines to help in creating the legend
  additional_lines[[i]] <- data.frame(
    theta = rep(thetaseq, 3),
    value = c(
      dnorm(x = thetaseq, mean = tr[i], sd = sr[i]),
      dnorm(x = thetaseq, mean = to, sd = so),
      dnorm(x = thetaseq, mean = mu_UIP, sd = sqrt(tau_UIP))
    ),
    replication = paste0("Replication ", i),
    linetype = factor(rep(c("Likelihood", "Non-Informative Prior", "Prior Original"), each = length(thetaseq)))
  )
}


df_densities <- do.call(rbind, df)
df_additional_lines <- do.call(rbind, additional_lines)


df_densities$replication <- factor(df_densities$replication,
  levels = c("Replication 1", "Replication 2", "Replication 3", "Replication 4"),
  labels = c(
    expression(" " ~ hat(theta)[r * 1] == 0.29 ~ ", " ~ sigma[r * 1] == 0.11),
    expression(" " ~ hat(theta)[r * 2] == 0.25 ~ ", " ~ sigma[r * 2] == 0.09),
    expression(" " ~ hat(theta)[r * 3] == -0.18 ~ ", " ~ sigma[r * 3] == 0.11),
    expression(" " ~ hat(theta)[r * p] == 0.14 ~ ", " ~ sigma[r * p] == 0.06)
  )
)

df_additional_lines$replication <- factor(df_additional_lines$replication,
  levels = c("Replication 1", "Replication 2", "Replication 3", "Replication 4"),
  labels = c(
    expression(" " ~ hat(theta)[r * 1] == 0.29 ~ ", " ~ sigma[r * 1] == 0.11),
    expression(" " ~ hat(theta)[r * 2] == 0.25 ~ ", " ~ sigma[r * 2] == 0.09),
    expression(" " ~ hat(theta)[r * 3] == -0.18 ~ ", " ~ sigma[r * 3] == 0.11),
    expression(" " ~ hat(theta)[r * p] == 0.14 ~ ", " ~ sigma[r * p] == 0.06)
  )
)
# The ggplot
plot_post_fix <- ggplot() +
  geom_line(
    data = df_densities, aes(x = theta, y = density, color = w),
    linewidth = 1
  ) +
  geom_line(
    data = df_additional_lines, aes(x = theta, y = value, linetype = linetype),
    linewidth = 0.8
  ) +
  scale_color_manual(values = cols) +
  scale_linetype_manual(
    values = c("dashed", "dotted", "dotdash"),
    labels = c("Likelihood (Replication)", "Prior (Original component)", "Prior (Non-Informative component)")
  ) +
  labs(x = expression("Effect Size" ~ theta), y = "Density") +
  theme_bw() +
  guides(
    linetype = guide_legend(title = "Density: ", position = "top"),
    color = guide_legend(title = "Weight: ", position = "left")
  ) +
  facet_wrap(~replication, labeller = label_parsed, ncol = 4) +
  scale_x_continuous(limits = c(-0.7, 0.7)) +
  theme(
    strip.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 18),
    axis.title.x = element_text(size = 22),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 19)
  )


print(plot_post_fix)

ggsave(
  filename = "plot_post_fix.pdf", path = "Plots/Mixture_Prior", plot = plot_post_fix,
  width = 17, height = 7.5, device = "pdf", dpi = 500, useDingbats = FALSE
)


#   ____________________________________________________________________________
#   Effect Size HPDI                                                        ####

HPDI_theta_w_median <- do.call("rbind", lapply(seq(1, length(wseq)), function(i) {
  # Nested lapply for each element of tr and sr
  results <- lapply(seq_along(tr), function(j) {
    hpd <- thetaHPD(
      level = 0.95, tr = tr[j], sr = sr[j], to = to,
      so = so, w = wseq[i], m = mu_UIP, v = tau_UIP
    )
    median_vect <- median_fun(thetaseq,
      tr = tr[j], sr = sr[j], to = to,
      so = so, w = wseq[i], m = mu_UIP, v = tau_UIP
    )
    out <- data.frame(
      lower = hpd[1], upper = hpd[3], median = median_vect, weight = wseq[i],
      tr_val = tr[j], sr_val = sr[j], # Add tr and sr values to the output
      width_int = (hpd[3] - hpd[1]), rnumber = j
    )
    return(out)
  })
  do.call("rbind", results) # Combine results for each tr and sr pair
}))
HPDI_theta_w_median <- HPDI_theta_w_median[order(HPDI_theta_w_median$rnumber), ]
HPDI_theta_w_median$rep_order <- ifelse(HPDI_theta_w_median$rnumber == 4, "p",
  as.character(HPDI_theta_w_median$rnumber)
)

# Create rep_setting labels for the HPDI object (which will be parsed as mathematical expressions)
HPDI_theta_w_median$rep_setting <- paste0(
  "{hat(theta)[r*",
  ifelse(HPDI_theta_w_median$rnumber == 4, "p", HPDI_theta_w_median$rnumber),
  "] == ",
  round(HPDI_theta_w_median$tr, 2),
  "}*',' ~ sigma[r*",
  ifelse(HPDI_theta_w_median$rnumber == 4, "p", HPDI_theta_w_median$rnumber),
  "] == ",
  round(HPDI_theta_w_median$sr, 2)
)
#   ____________________________________________________________________________
#   Tipping Point Analysis                                                  ####


# Confidence Interval Replication and Original Study


# Compute the z-value for 95% confidence intervals
z_value <- qnorm(0.975)

# Confidence interval for the original study
ci_o <- c(to - z_value * so, to + z_value * so)

# Confidence intervals for replication studies (and pooled)
ci_r1 <- c(tr[1] - z_value * sr[1], tr[1] + z_value * sr[1])
ci_r2 <- c(tr[2] - z_value * sr[2], tr[2] + z_value * sr[2])
ci_r3 <- c(tr[3] - z_value * sr[3], tr[3] + z_value * sr[3])
ci_rp <- c(tr[4] - z_value * sr[4], tr[4] + z_value * sr[4])

# Create a data frame for the replication study CIs
data_ci_rep <- data.frame(
  median = c(tr[1], tr[2], tr[3], tr[4]),
  sr_val = c(sr[1], sr[2], sr[3], sr[4]),
  ymin = c(ci_r1[1], ci_r2[1], ci_r3[1], ci_rp[1]),
  ymax = c(ci_r1[2], ci_r2[2], ci_r3[2], ci_rp[2]),
  weight = rep(-0.15, 4)
)

# Construct rep_setting labels for replication CIs (using dynamic values)
data_ci_rep$rep_setting <- paste0(
  "{hat(theta)[r*",
  ifelse(1:4 == 4, "p", 1:4),
  "] == ",
  round(data_ci_rep$median, 2),
  "}*',' ~ sigma[r*",
  ifelse(1:4 == 4, "p", 1:4),
  "] == ",
  round(data_ci_rep$sr_val, 2)
)

# Create a data frame for the original study CIs (repeating the same original value for each facet)
data_ci_orig <- data.frame(
  median = rep(to, 4),
  sr_val = rep(so, 4),
  ymin = rep(ci_o[1], 4),
  ymax = rep(ci_o[2], 4),
  weight = rep(1.15, 4)
)

# Construct rep_setting labels for the original study
data_ci_orig$rep_setting <- paste0(
  "{hat(theta)[r*",
  ifelse(1:4 == 4, "p", 1:4),
  "] == ",
  round(data_ci_rep$median, 2),
  "}*',' ~ sigma[r*",
  ifelse(1:4 == 4, "p", 1:4),
  "] == ",
  round(data_ci_rep$sr_val, 2)
)


palette_colors <- rep(cols, 4)


# Plotting the error bars and points
plot_HPDI_median_rep <- ggplot(HPDI_theta_w_median, aes(x = weight, y = median)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.05, size = 1.4, color = palette_colors) +
  geom_errorbar(
    data = data_ci_rep, aes(ymin = ymin, ymax = ymax),
    width = 0.05, size = 1.4, linetype = 5, color = "black"
  ) +
  geom_errorbar(
    data = data_ci_orig, aes(ymin = ymin, ymax = ymax),
    width = 0.05, size = 1.4, linetype = 5, color = "firebrick4"
  ) +
  geom_point(data = data_ci_rep, aes(x = weight, y = median), shape = 16, size = 4, color = "black") +
  geom_point(data = data_ci_orig, aes(x = weight, y = median), shape = 16, size = 4, color = "firebrick4") +
  geom_point(shape = 16, size = 4, color = palette_colors) +
  geom_hline(yintercept = 0, linetype = "dotdash", color = "black") +
  labs(x = "Prior Weight", y = "Effect Size Posterior") +
  theme_bw() +
  theme(
    strip.text.x = element_text(size = 18),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 22),
    legend.text = element_text(size = 18)
  ) +
  facet_wrap(~rep_setting, labeller = label_parsed, ncol = 4) +
  scale_x_continuous(
    breaks = c(-0.15, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.15),
    labels = c("Rep.", 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, "Orig.")
  )

# Display the plot
print(plot_HPDI_median_rep)


ggsave(
  filename = "plot_HPDI_median_rep.pdf", path = "Plots/Mixture_Prior", plot = plot_HPDI_median_rep,
  width = 17, height = 7.5, device = "pdf", dpi = 500, useDingbats = FALSE
)
