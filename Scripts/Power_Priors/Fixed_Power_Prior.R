#   ____________________________________________________________________________
#   Libraries                                                               ####

library(ppRep)
library(ggplot2)
library(colorspace)
library(xtable)
library(dplyr)
library(hypergeo)
library(ggpubr)
library(spatstat)
library(RColorBrewer)

#   ____________________________________________________________________________
#   Parameter Setting                                                       ####

# Load data and credentials
load("credentials_data.RData")

# Extract effect sizes and standard errors for original study
to <- data %>%
  dplyr::filter(type == "original") %>%
  dplyr::pull(fis) %>%
  as.numeric()
so <- data %>%
  dplyr::filter(site == "original") %>%
  dplyr::pull(se_fis) %>%
  as.numeric()

# Extract effect sizes and standard errors for replication studies
trep <- data %>%
  dplyr::filter(site %in% c("University of Toronto", "Montana State University", "Ashland University")) %>%
  dplyr::pull(fis) %>%
  as.numeric()
srep <- data %>%
  dplyr::filter(site %in% c("University of Toronto", "Montana State University", "Ashland University")) %>%
  dplyr::pull(se_fis) %>%
  as.numeric()

# Compute a pooled replication effect and standard error
tp <- round(sum(trep / srep^2) / sum(1 / srep^2), 2)
sp <- round(sqrt(1 / sum(1 / srep^2)), 2)

# Combine replication and pooled values
tr <- c(trep, tp)
sr <- c(srep, sp)
rnumber <- c(1, 2, 3, 4)

# Uniform prior for alpha
x <- 1
y <- 1

#   ____________________________________________________________________________
#   Parameter Grid and Other Settings                                       ####

n_alpha <- 400
n_theta <- 400
alphaseq <- seq(0, 1, by = 0.1)
thetaseq <- seq(-0.9, 0.9, length.out = n_theta)
par_grid <- expand.grid(alpha = alphaseq, theta = thetaseq)
m <- 0
v <- Inf

# Define a color palette; ensure that 'wseq' is defined in your workspace.
cols <- hcl.colors(n = length(alphaseq), palette = "viridis", alpha = 0.9, rev = TRUE)

#   ____________________________________________________________________________
#   Median Function                                                         ####

median_fun_pp <- function(theta, tr, sr, to, so, m = 0, v = Inf, alpha) {
  # Calculate posterior density for theta
  density <- postPPtheta(
    theta = thetaseq, tr = tr, sr = sr, to = to,
    so = so, x = x, y = y, m = m, v = v, alpha = alpha
  )
  dat <- data.frame(theta, density)
  dat <- dat[order(dat$theta), ]
  # Compute weighted median
  w_median <- weighted.median(dat[, 1], dat[, 2], na.rm = TRUE)
  return(w_median)
}

#   ____________________________________________________________________________
#   Effect Size HPDI                                                        ####

HPDI_theta_pp_median <- do.call("rbind", lapply(seq_along(alphaseq), function(i) {
  # For each value of alpha, compute HPDI and median for each replication (and pooled) study
  results <- lapply(seq_along(tr), function(j) {
    hpd <- postPPthetaHPD(level = 0.95, tr = tr[j], sr = sr[j], to = to, so = so,
                          alpha = alphaseq[i])
    median_vect <- median_fun_pp(thetaseq, tr = tr[j], sr = sr[j], to = to,
                                 so = so, alpha = alphaseq[i])
    out <- data.frame(lower = hpd[1],
                      upper = hpd[2],
                      median = median_vect,
                      alpha = alphaseq[i],
                      tr = tr[j],
                      sr = sr[j],
                      width_int = (hpd[2] - hpd[1]),
                      rnumber = j)
    return(out)
  })
  do.call("rbind", results)
}))

# Order the results and create a label for the replication order (with 4 labeled as "p")
HPDI_theta_pp_median <- HPDI_theta_pp_median[order(HPDI_theta_pp_median$rnumber), ]
HPDI_theta_pp_median$rep_order <- ifelse(HPDI_theta_pp_median$rnumber == 4, "p",
                                         as.character(HPDI_theta_pp_median$rnumber))

# Create rep_setting labels for the HPDI object
HPDI_theta_pp_median$rep_setting <- paste0(
  "{hat(theta)[r*",
  ifelse(HPDI_theta_pp_median$rnumber == 4, "p", HPDI_theta_pp_median$rnumber),
  "] == ",
  round(HPDI_theta_pp_median$tr, 2),
  "}*',' ~ sigma[r*",
  ifelse(HPDI_theta_pp_median$rnumber == 4, "p", HPDI_theta_pp_median$rnumber),
  "] == ",
  round(HPDI_theta_pp_median$sr, 2)
)

#   ____________________________________________________________________________
#   Tipping Point Analysis: Confidence Intervals                          ####

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
  alpha = rep(-0.15, 4)
)

# Construct rep_setting labels for replication CIs 
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

# Create a data frame for the original study CIs
data_ci_orig <- data.frame(
  median = rep(to, 4),
  sr_val = rep(so, 4),
  ymin = rep(ci_o[1], 4),
  ymax = rep(ci_o[2], 4),
  alpha = rep(1.15, 4)
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

#   ____________________________________________________________________________
#   Plotting                                                                ####

palette_colors <- rep(cols, 4)

plot_HPDI_median_rep_pp <- ggplot(HPDI_theta_pp_median, aes(x = alpha, y = median)) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.05, size = 1.4, color = palette_colors) +
  geom_errorbar(data = data_ci_rep, aes(ymin = ymin, ymax = ymax),
                width = 0.05, size = 1.4, linetype = 5, color = "black") +
  geom_errorbar(data = data_ci_orig, aes(ymin = ymin, ymax = ymax),
                width = 0.05, size = 1.4, linetype = 5, color = "firebrick4") +
  geom_point(data = data_ci_rep, aes(x = alpha, y = median),
             shape = 16, size = 4, color = "black") +
  geom_point(data = data_ci_orig, aes(x = alpha, y = median),
             shape = 16, size = 4, color = "firebrick4") +
  geom_point(shape = 16, size = 4, color = palette_colors) +
  geom_hline(yintercept = 0, linetype = "dotdash", color = "black") +  
  labs(x = "Prior alpha", y = "Effect Size Posterior") +
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
  # Facet on the rep_setting label; label_parsed will render the math expressions
  facet_wrap(~ rep_setting, labeller = label_parsed, ncol = 4) +
  scale_x_continuous(breaks = c(-0.15, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.15),
                     labels = c("Rep.", 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, "Orig."))

# Display the plot
print(plot_HPDI_median_rep_pp)


ggsave(filename = "plot_HPDI_median_rep_pp.pdf",path = "Plots/Power_Prior", plot = plot_HPDI_median_rep_pp,
       width = 17, height = 7.5, device='pdf', dpi=500, useDingbats = FALSE)

