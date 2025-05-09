#   ____________________________________________________________________________
#   Libraries                                                               ####

library(ggplot2)
library(ggthemes)
library(ggpubr)
library(colorspace)


#   ____________________________________________________________________________
#   Sources                                                                 ####
source("Scripts/Mixture_Priors/Random_Weights.R")

#   ____________________________________________________________________________
#   Plots                                                                   ####


##  ............................................................................
##  Plot Effect Estimate and 95% CI                                         ####

# z-value for 95% CI
z_value <- qnorm(0.975)

# Calculate 95% CI for each group
ci_o <- c(to - z_value * so, to + z_value * so)
ci_r1 <- c(tr[1] - z_value * sr[1], tr[1] + z_value * sr[1])
ci_r2 <- c(tr[2] - z_value * sr[2], tr[2] + z_value * sr[2])
ci_r3 <- c(tr[3] - z_value * sr[3], tr[3] + z_value * sr[3])
ci_r4 <- c(tr[4] - z_value * sr[4], tr[4] + z_value * sr[4])


# Create a data frame
data <- data.frame(
  group = c("Original", "Replication 1", "Replication 2", "Replication 3", "Replication Pooled"),
  estimate = c(to, tr[1], tr[2], tr[3], tr[4]),
  ymin = c(ci_o[1], ci_r1[1], ci_r2[1], ci_r3[1], ci_r4[1]),
  ymax = c(ci_o[2], ci_r1[2], ci_r2[2], ci_r3[2], ci_r4[2]),
  color = c("#8A0404", "#E69F00", "#009E20", "#0072B2", "#AA4499")
)

# Plot

plot_theta <- ggplot(data) +
  geom_point(aes(x = group, y = estimate, color = group), shape = 16, fill = "white", size = 4.5) +
  geom_errorbar(aes(x = group, ymin = ymin, ymax = ymax, color = group), width = 0.15, size = 1.8) +
  scale_color_manual(values = data$color) +
  labs(
    y = "Effect Size Estimate",
    x = "Study"
  ) +
  scale_color_manual(
    values = c(
      "Original" = "#8A0404", "Replication 1" = "#E69F00", "Replication 2" = "#009E20", "Replication 3" = "#0072B2",
      "Replication Pooled" = "#AA4499"
    ),
    labels = c(
      expression(" " ~ hat(theta)[o] == 0.21 ~ ", " ~ sigma[o] == 0.06),
      expression(" " ~ hat(theta)[r * 1] == 0.29 ~ ", " ~ sigma[r * 1] == 0.11),
      expression(" " ~ hat(theta)[r * 2] == 0.25 ~ ", " ~ sigma[r * 2] == 0.09),
      expression(" " ~ hat(theta)[r * 3] == -0.18 ~ ", " ~ sigma[r * 3] == 0.11),
      expression(" " ~ hat(theta)[r * p] == 0.14 ~ ", " ~ sigma[r * p] == 0.06)
    )
  ) +
  theme_bw() +
  theme(
    # axis.title.x = element_blank(),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 18),
    axis.title.x = element_text(size = 22),
    legend.text = element_text(size = 18)
  )


# pdf("ggplot.pdf", pointsize=100, width=20, height=10)
print(plot_theta)

# dev.off()

ggsave(
  filename = "plot_theta.pdf", path = "Plots", plot = plot_theta,
  width = 17, height = 7.5, device = "pdf", dpi = 500, useDingbats = FALSE
)

##  ............................................................................
##  Contour plot of joint posterior                                         ####

# Create a new column for ordering
postdens_wrapper$rep_order <- ifelse(postdens_wrapper$rnumber == 4, "p", as.character(postdens_wrapper$rnumber))

# Create the original rep_setting column with proper labels
postdens_wrapper$rep_setting <- paste0(
  "{hat(theta)[r*",
  ifelse(postdens_wrapper$rnumber == 4, "p", postdens_wrapper$rnumber),
  "] == ",
  round(postdens_wrapper$tr, 2),
  "}*',' ~ sigma[r*",
  ifelse(postdens_wrapper$rnumber == 4, "p", postdens_wrapper$rnumber),
  "] == ",
  round(postdens_wrapper$sr, 2)
)

# Convert the new rep_order column to a factor to specify the order
postdens_wrapper$rep_order <- factor(postdens_wrapper$rep_order, levels = c("1", "2", "3", "p"))


plot_joint <- ggplot(data = postdens_wrapper, aes(x = theta, y = omega, fill = density)) +
  facet_wrap(~rep_setting, labeller = label_parsed, ncol = 4) +
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(z = density), breaks = seq(0, 24, 3), alpha = 0.35, size = 0.5) +
  scale_fill_continuous_sequential(palette = "Blues 3", rev = TRUE) +
  labs(
    x = bquote("Effect Size" ~ theta),
    y = bquote("Weight Parameter" ~ omega),
    fill = "Posterior \n Density"
  ) +
  guides(fill = guide_colorbar(barheight = 20, barwidth = 1.3, title.position = "top")) +
  scale_x_continuous(limits = c(-0.5, 0.6), expand = c(0.04, 0.005)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0.005, 0.005)) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.text.x = element_text(size = 18),
    legend.position = "right",
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 22),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.key.size = unit(10, "cm")
  )

print(plot_joint)

ggsave(
  filename = "plot_joint.pdf", path = "Plots/Mixture_Prior", plot = plot_joint,
  width = 17, height = 7.5, device = "pdf", dpi = 500, useDingbats = FALSE
)


##  ............................................................................
##  Plot Marginal Posterior weights                                         ####

plot_weights_m_hpd <- ggplot() +
  geom_errorbarh(
    data = HPDI_weights,
    aes(
      xmin = lower, xmax = upper, y = y * 1.05, color = factor(rnumber),
      height = height
    ), alpha = 0.8, size = 1.2
  ) +
  geom_line(
    data = weights_m_post, aes(x = x, y = density, group = rnumber, color = factor(rnumber)),
    lty = 1, alpha = 0.9, size = 1.2
  ) +
  scale_color_manual(
    values = c("1" = "#E69F00", "2" = "#009E20", "3" = "#0072B2", "4" = "#AA4499"),
    labels = c(
      expression(" " ~ hat(theta)[r * 1] == 0.29 ~ ", " ~ sigma[r * 1] == 0.11),
      expression(" " ~ hat(theta)[r * 2] == 0.25 ~ ", " ~ sigma[r * 2] == 0.09),
      expression(" " ~ hat(theta)[r * 3] == -0.18 ~ ", " ~ sigma[r * 3] == 0.11),
      expression(" " ~ hat(theta)[r * p] == 0.14 ~ ", " ~ sigma[r * p] == 0.06)
    )
  ) +
  labs(
    x = expression(omega ~ " Values"),
    y = "Density"
  ) +
  theme_bw() +
  theme(
    strip.placement = "outside", # format to look like title
    strip.background = element_blank(),
    strip.text.x = element_text(size = 18),
    legend.position = "top",
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 22),
    legend.title = element_blank(),
    legend.text = element_text(size = 16)
  ) +
  guides(color = guide_legend(title = "Replicated Experiment"))

print(plot_weights_m_hpd)


# Flat prior and linear combination of omegas



##  ............................................................................
##  Marginal Posterior Theta                                                ####

plot_m_theta <- ggplot() +
  geom_errorbarh(
    data = HPDI_theta,
    aes(
      xmin = lower, xmax = upper, y = y * 1.05, color = factor(rnumber),
      height = height
    ), alpha = 0.8, size = 1.0
  ) +
  geom_errorbarh(
    data = HPDI_theta_2,
    aes(
      xmin = lower, xmax = upper, y = y * 1.05, color = factor(rnumber),
      height = height
    ), alpha = 0.7, linetype = "22", size = 1.0
  ) +
  geom_line(
    data = theta_m_post_2, aes(x = x, y = density, color = factor(rnumber)),
    lty = 22, alpha = 0.5, size = 1.0
  ) +
  geom_line(
    data = theta_m_post, aes(x = x, y = density, color = factor(rnumber)),
    alpha = 0.9, size = 1.0
  ) +
  scale_color_manual(
    values = c("1" = "#E69F00", "2" = "#009E20", "3" = "#0072B2", "4" = "#AA4499"),
    labels = c(
      expression(" " ~ hat(theta)[r * 1] == 0.29 ~ ", " ~ sigma[r * 1] == 0.11),
      expression(" " ~ hat(theta)[r * 2] == 0.25 ~ ", " ~ sigma[r * 2] == 0.09),
      expression(" " ~ hat(theta)[r * 3] == -0.18 ~ ", " ~ sigma[r * 3] == 0.11),
      expression(" " ~ hat(theta)[r * p] == 0.14 ~ ", " ~ sigma[r * p] == 0.06)
    )
  ) +
  labs(
    x = expression(theta ~ " Values"),
    y = "Density"
  ) +
  xlim(-0.6, 0.6) +
  theme_bw() +
  theme(
    strip.placement = "outside", # format to look like title
    strip.background = element_blank(),
    strip.text.x = element_text(size = 18),
    legend.position = "top",
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 22),
    legend.title = element_blank(),
    legend.text = element_text(size = 16)
  ) +
  guides(color = guide_legend(title = "Replicated Experiment"))


print(plot_m_theta)


#   ____________________________________________________________________________
#   Arrange of Fig. 6 (Marginal Posterior Densities)                        ####


plot_marg_post_joint <- ggplot() +
  geom_errorbarh(
    data = HPDI_theta,
    aes(
      xmin = lower, xmax = upper, y = y * 1.05, color = factor(rnumber),
      height = height
    ), alpha = 0.8, size = 1.0
  ) +
  geom_errorbarh(
    data = HPDI_theta_2,
    aes(
      xmin = lower, xmax = upper, y = y * 1.05, color = factor(rnumber),
      height = height
    ), alpha = 0.7, linetype = "22", size = 1.0
  ) +
  geom_line(
    data = theta_m_post_2, aes(x = x, y = density, color = factor(rnumber)),
    lty = 22, alpha = 0.5, size = 1.0
  ) +
  geom_line(
    data = theta_m_post, aes(x = x, y = density, color = factor(rnumber)),
    alpha = 0.9, size = 1.0
  ) +
  geom_errorbarh(
    data = HPDI_weights,
    aes(
      xmin = lower, xmax = upper, y = y * 1.05, color = factor(rnumber),
      height = height
    ), alpha = 0.8, size = 1.2
  ) +
  geom_line(
    data = weights_m_post, aes(x = x, y = density, group = rnumber, color = factor(rnumber)),
    lty = 1, alpha = 0.9, size = 1.2
  ) +
  labs(
    x = NULL,
    y = "Marginal Posterior Density",
    color = ""
  ) +
  scale_color_manual(
    values = c("1" = "#E69F00", "2" = "#009E20", "3" = "#0072B2", "4" = "#AA4499"),
    labels = c(
      expression(" " ~ hat(theta)[r * 1] == 0.29 ~ ", " ~ sigma[r * 1] == 0.11),
      expression(" " ~ hat(theta)[r * 2] == 0.25 ~ ", " ~ sigma[r * 2] == 0.09),
      expression(" " ~ hat(theta)[r * 3] == -0.18 ~ ", " ~ sigma[r * 3] == 0.11),
      expression(" " ~ hat(theta)[r * p] == 0.14 ~ ", " ~ sigma[r * p] == 0.06)
    )
  ) +
  facet_wrap(~parameter,
    scales = "free", labeller = label_parsed,
    strip.position = "bottom"
  ) +
  theme_bw() +
  theme(
    strip.placement = "outside", # format to look like title
    strip.background = element_blank(),
    strip.text.x = element_text(size = 22),
    legend.position = "top",
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 18),
    axis.title.x = element_text(size = 22),
    legend.title = element_blank(),
    legend.text = element_text(size = 18)
  ) +
  guides(color = guide_legend(title = "Replicated Experiment"))

print(plot_marg_post_joint)

ggsave(
  filename = "plot_marg_post_joint.pdf", path = "Plots/Mixture_Prior",
  plot = plot_marg_post_joint,
  width = 17, height = 7.5, device = "pdf", dpi = 500, useDingbats = FALSE
)

