#   ____________________________________________________________________________
#   Libraries                                                               ####

library(ppRep)
library(ggplot2)
library(colorspace)
library(xtable)
library(dplyr)
library(hypergeo)
library(ggpubr)


#   ____________________________________________________________________________
#   Sources                                                                 ####

source("Scripts/Power_Priors/Random_Power_Prior.R")


#   ____________________________________________________________________________
#   Plots                                                                   ####


##  ............................................................................
##  Contour plot of joint posterior                                         ####

pp_plot_joint <- ggplot(data = joint_plot_df, aes(x = theta, y = alpha, fill = density)) +
  facet_wrap(~rep_setting,
    labeller = label_parsed, ncol = 4
  ) +
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(z = density),
    breaks = seq(0, 30, 5), alpha = 0.25, col = 1,
    linewidth = 0.3
  ) +
  scale_fill_continuous_sequential(palette = "Blues 3", rev = TRUE) +
  labs(
    x = bquote("Effect size" ~ theta),
    y = bquote("Power parameter" ~ alpha),
    fill = "Posterior \ndensity"
  ) +
  guides(fill = guide_colorbar(barheight = 10, barwidth = 0.5)) +
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

ggsave(
  filename = "pp_plot_joint.pdf", path = "Plots/Power_Prior", plot = pp_plot_joint,
  width = 17, height = 7.5, device = "pdf", dpi = 500, useDingbats = FALSE
)

##  ............................................................................
##  Plot marginal posteriors                                                ####

ncat <- length(unique(marg_plot_df$rep_setting))
colblind <- c("#E69F00", "#009E20", "#0072B2", "#AA4499")
names(colblind) <- unique(marg_plot_df$rep_setting)
colblind <- colblind[order(names(colblind))]


plot_marg_post_joint <- ggplot() +
  geom_errorbarh(
    data = HPDI_df,
    aes(
      xmin = lower, xmax = upper, y = y * 1.05, color = rep_setting,
      height = height
    ), alpha = 0.8, size = 1.0
  ) +
  geom_errorbarh(
    data = HPDI_theta_2,
    aes(
      xmin = lower, xmax = upper, y = y * 1.05, color = rep_setting,
      height = height
    ), alpha = 0.7, linetype = "22", size = 1.0
  ) +
  geom_line(
    data = marg_theta_dens2, aes(x = x, y = density, color = rep_setting),
    lty = 2, alpha = 0.5, linewidth = 1.0
  ) +
  geom_line(
    data = marg_plot_df, aes(x = x, y = density, color = rep_setting),
    alpha = 0.9, linewidth = 1.0
  ) +
  # geom_line(data = alpha_limit_df, aes(x = x, y = density), col = 1, lty = 3, alpha = 0.9, size = 1.0) +
  facet_wrap(~parameter,
    scales = "free", labeller = label_parsed,
    strip.position = "bottom"
  ) +
  theme_bw() +
  labs(x = NULL, y = "Marginal posterior density", color = "") +
  scale_color_manual(values = colblind, labels = scales::parse_format()) +
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


ggsave(
  filename = "plot_marg_post_joint.pdf", path = "Plots/Power_Prior", plot = plot_marg_post_joint,
  width = 17, height = 7.5, device = "pdf", dpi = 500, useDingbats = FALSE
)
## Combine all plots
ggpubr::ggarrange(pp_plot_joint, plot_marg_post_joint, ncol = 1)
