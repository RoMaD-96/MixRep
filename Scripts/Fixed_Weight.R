#   ____________________________________________________________________________
#   Libraries                                                               ####

library(ggplot2)
library(ggpubr)
library(colorspace)
library(tipmap)
library(spatstat)

source("Scripts/RepMixFun.R")


#   ____________________________________________________________________________
#   Fixed Weights                                                           ####


##  ............................................................................
##  Parameter Setting                                                       ####

# Original and Replicated Studies
to <- 0.21
so <- 0.05
tr <- c(0.09, 0.21, 0.44)
sr <- c(0.045, 0.06, 0.04)
null <- 0
priorsd <- 2
wseq <- seq(0, 1, 0.1)
thetaseq <- seq(-6, 6, length.out = 2500)




#   ____________________________________________________________________________
#   Posterior with varying weights using ggplot 2                           ####


cols <- hcl.colors(n = length(wseq), palette = "viridis", alpha = 0.9, rev = TRUE)


# Now let's create a data frame for ggplot
densities <- list()
for (i in 1:length(tr)) {
densities[[i]] <- sapply(X = wseq, FUN = function(w) {
  rmapPostFix_alt(theta = thetaseq, tr = tr[i], sr = sr[i], to = to, so = so,
           null = null, priorsd = priorsd, w = w)
  })
}


df <- list()
additional_lines <- list()

for (i in 1:length(tr)) {
  # Turn the results into a data frame
  df[[i]] <- data.frame(theta = rep(thetaseq, times = length(wseq)),
                   density = as.vector(densities[[i]]),
                   w = factor(rep(wseq, each = length(thetaseq))))
  
  # Add the additional lines
  df[[i]]$likelihood <- dnorm(x = df[[i]]$theta, mean = tr[i], sd = sr[i])
  df[[i]]$prior_original <- dnorm(x = df[[i]]$theta, mean = to, sd = so)
  df[[i]]$prior_robust <- dnorm(x = df[[i]]$theta, mean = null, sd = priorsd)
  df[[i]]$replication <- paste0( "Replication ", i)
  
  # Create a separate data frame for the additional lines to help in creating the legend
  additional_lines[[i]] <- data.frame(
    theta = rep(thetaseq, 3),
    value = c(dnorm(x = thetaseq, mean = tr[i], sd = sr[i]),
              dnorm(x = thetaseq, mean = to, sd = so),
              dnorm(x = thetaseq, mean = null, sd = priorsd)),
    replication = paste0( "Replication ", i),
    linetype = factor(rep(c("Likelihood", "Prior Original", "Prior Robust"), each = length(thetaseq)))
  )
  
}


df_densities <- do.call(rbind, df)
df_additional_lines <- do.call(rbind, additional_lines)


df_densities$replication <- factor(df_densities$replication,
                                     levels = c("Replication 1", "Replication 2", "Replication 3"), 
                                     labels = c(expression(" "~hat(theta)[r*1] == 0.09 ~ ", " ~ sigma[r*1] == 0.045),
                                                expression(" "~hat(theta)[r * 2] == 0.21 ~ ", " ~ sigma[r*2] == 0.06),
                                                expression(" "~hat(theta)[r * 3] == 0.44 ~ ", " ~ sigma[r*3] == 0.04)))

df_additional_lines$replication <- factor(df_additional_lines$replication,
                                       levels = c("Replication 1", "Replication 2", "Replication 3"), 
                                       labels = c(expression(" "~hat(theta)[r*1] == 0.09 ~ ", " ~ sigma[r*1] == 0.045),
                                                  expression(" "~hat(theta)[r * 2] == 0.21 ~ ", " ~ sigma[r*2] == 0.06),
                                                  expression(" "~hat(theta)[r * 3] == 0.44 ~ ", " ~ sigma[r*3] == 0.04)))
# The ggplot
plot_post_fix <- ggplot() + 
  geom_line(data = df_densities, aes(x = theta, y = density, color = w), size = 1) + 
  geom_line(data = df_additional_lines, aes(x = theta, y = value, linetype = linetype), size = 0.8) +
  scale_color_manual(values = cols) +
  scale_linetype_manual(values = c("dashed", "dotted", "dotdash"),
                        labels = c("Likelihood (Replication)", "Prior (Original component)", "Prior (Robust component)")) +
  labs(x = expression("Effect Size" ~ theta), y = "Density") +
  theme_bw() +
  guides(linetype = guide_legend(title = "Density"),
         color = guide_legend(title = "Weight")
  ) +
  facet_wrap(~ replication, labeller = label_parsed)

# To make sure that our additional lines are represented in the legend, we need to add them to the plot
plot_post_fix <- plot_post_fix + geom_line(aes(linetype = "Likelihood (Replication)"), linetype = "dashed", color = "black") +
  geom_line(aes(linetype = "Prior (Original component)"), linetype = "dotted", color = "black") +
  geom_line(aes(linetype = "Prior (Robust component)"), linetype = "dotdash", color = "black") +
  scale_x_continuous(limits=c(-0.10, 0.6)) +
  theme(legend.position = "right",
        strip.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_text(size = 22),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 19))

ggsave(filename = "plot_post_fix.pdf",path = "Plots", plot = plot_post_fix,
       width = 16.5, height = 8.5, device='pdf', dpi=500, useDingbats = FALSE)


#   ____________________________________________________________________________
#   Effect Size HPDI                                                        ####

HPDI_theta_w_median <- do.call("rbind", lapply(seq(1, length(wseq)), function(i) {
  # Nested lapply for each element of tr and sr
  results <- lapply(seq_along(tr), function(j) {
    hpd <- HPDI_post_m_theta_fix(level = 0.95, tr = tr[j], sr = sr[j], to = to,
                                 so = so, w = wseq[i], null = null, priorsd = priorsd)
    median_vect <- median_fun(thetaseq, tr = tr[j], sr = sr[j], to = to,
                              so = so, w = wseq[i], null = null, priorsd = priorsd)
    out <- data.frame(lower = hpd[1], upper = hpd[2], median = median_vect, weight = wseq[i],
                      tr_val = tr[j], sr_val = sr[j], # Add tr and sr values to the output
                      width_int = (hpd[2]-hpd[1]), replication = paste0( "Replication ", j))
    return(out)
  })
  do.call("rbind", results) # Combine results for each tr and sr pair
}))
HPDI_theta_w_median <- HPDI_theta_w_median[order(HPDI_theta_w_median$replication), ]

HPDI_theta_w_median$replication  <-  factor(HPDI_theta_w_median$replication,
                                            levels = c("Replication 1", "Replication 2", "Replication 3"), 
                                            labels = c(expression(" "~hat(theta)[r*1] == 0.09 ~ ", " ~ sigma[r*1] == 0.045),
                                                       expression(" "~hat(theta)[r * 2] == 0.21 ~ ", " ~ sigma[r*2] == 0.060),
                                                       expression(" "~hat(theta)[r * 3] == 0.44 ~ ", " ~ sigma[r*3] == 0.040)))



#   ____________________________________________________________________________
#   Tipping Point Analysis                                                  ####


# Confidence Interval Replication and Original Study

# z-value for 95% CI
z_value <- qnorm(0.975)

# Calculate 95% CI for each group
ci_o <- c(to - z_value * so, to + z_value * so)
ci_r1 <- c(tr[1] - z_value * sr[1], tr[1] + z_value * sr[1])
ci_r2 <- c(tr[2] - z_value * sr[2], tr[2] + z_value * sr[2])
ci_r3 <- c(tr[3] - z_value * sr[3], tr[3] + z_value * sr[3])

# Create a data frame for the replication study
data_ci_rep <- data.frame(
  median = c(tr[1], tr[2], tr[3]),
  sr_val = c(sr[1], sr[2], sr[3]),
  ymin = c(ci_r1[1], ci_r2[1], ci_r3[1]),
  ymax = c(ci_r1[2], ci_r2[2], ci_r3[2]),
  weight = c(rep(-0.15,3)),
  replication = factor(c("Replication 1", "Replication 2", "Replication 3"),
                       levels = c("Replication 1", "Replication 2", "Replication 3"), 
                       labels = c(expression(" "~hat(theta)[r*1] == 0.09 ~ ", " ~ sigma[r*1] == 0.045),
                                  expression(" "~hat(theta)[r * 2] == 0.21 ~ ", " ~ sigma[r*2] == 0.060),
                                  expression(" "~hat(theta)[r * 3] == 0.44 ~ ", " ~ sigma[r*3] == 0.040)))
)

# Create a data frame for the original study
data_ci_orig <- data.frame(
  median = c(rep(to,3)),
  sr_val = c(rep(so,3)),
  ymin = c(rep(ci_o[1], 3)),
  ymax = c(rep(ci_o[2], 3)),
  weight = c(rep(1.15,3)),
  replication = factor(c("Replication 1", "Replication 2", "Replication 3"),
                       levels = c("Replication 1", "Replication 2", "Replication 3"), 
                       labels = c(expression(" "~hat(theta)[r*1] == 0.09 ~ ", " ~ sigma[r*1] == 0.045),
                                  expression(" "~hat(theta)[r * 2] == 0.21 ~ ", " ~ sigma[r*2] == 0.060),
                                  expression(" "~hat(theta)[r * 3] == 0.44 ~ ", " ~ sigma[r*3] == 0.040)))
)



# Plotting the error bars and points
plot_HPDI_median_rep <- ggplot(HPDI_theta_w_median, aes(x = weight, y = median)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = replication), width = 0.05, size = 1.4) +
  geom_errorbar(data = data_ci_rep, aes(ymin = ymin, ymax = ymax, color = "black"), 
                width = 0.05, size = 1.4, linetype = 1) +
  geom_errorbar(data = data_ci_orig, aes(ymin = ymin, ymax = ymax, color = "firebrick4"), 
                width = 0.05, size = 1.4, linetype = 1) +
  geom_point(data = data_ci_rep, aes(x = weight, y = median), shape = 16, size = 4, color = "black") +
  geom_point(data = data_ci_orig, aes(x = weight, y = median), shape = 16, size = 4, color = "firebrick4") +
  geom_point(aes(color = replication), shape = 16, size = 4) +
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
  facet_wrap(~ replication, labeller = label_parsed) +
  scale_color_manual(values = c("#E69F00", "#009E20", "#0072B2", "black", "firebrick4")) +
  scale_x_continuous(breaks = c(-0.15, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 ,1.15), 
                     labels = c("Rep.", 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, "Orig."))
# Display the plot
print(plot_HPDI_median_rep)

ggsave(filename = "plot_HPDI_median_rep.pdf",path = "Plots", plot = plot_HPDI_median_rep,
       width = 15, height = 8, device='pdf', dpi=500, useDingbats = FALSE)



#   ____________________________________________________________________________
#   Empirical Bayes                                                         ####


empirical_bayes <- list()
for (i in 1:length(tr)) {
  empirical_bayes[[i]] <- sapply(X = wseq, FUN = function(w) {
    rmapPostFix_alt(theta = thetaseq, tr = tr[i], sr = sr[i], to = to, so = so,
                    null = null, priorsd = priorsd, w = w)
  })
}



bfPPtheta <- function(tr, sr, to, so, x = 1, y = 1, alpha = NA, ...) {
  ## input checks
  stopifnot(
    length(tr) == 1,
    is.numeric(tr),
    is.finite(tr),
    
    length(to) == 1,
    is.numeric(to),
    is.finite(to),
    
    length(sr) == 1,
    is.numeric(sr),
    is.finite(sr),
    0 < sr,
    
    length(so) == 1,
    is.numeric(so),
    is.finite(so),
    0 < so,
    
    length(x) == 1,
    is.numeric(x),
    is.finite(x),
    0 <= x,
    
    length(y) == 1,
    is.numeric(y),
    is.finite(y),
    0 <= y,
    
    length(alpha) == 1,
    (is.na(alpha) |
       ((is.numeric(alpha)) &
          (is.finite(alpha)) &
          (0 <= alpha) &
          (alpha <= 1)))
  )
  
  ## marginal density under H0
  fH0 <- stats::dnorm(x = tr, mean = 0, sd = sr)
  
  ## marginal density under H1
  if (!is.na(alpha)) {
    fH1 <- stats::dnorm(x = tr, mean = to, sd = sqrt(sr^2 + so^2/alpha))
  } else {
    fH1 <- margLik(tr = tr, to = to, sr = sr, so = so, x = x, y = y,
                   ... = ...)
  }
  
  ## compute BF
  bf01 <- fH0/fH1
  return(bf01)
}
