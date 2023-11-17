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
thetaseq <- seq(-0.2, 0.6, length.out = 800)




#   ____________________________________________________________________________
#   Posterior with varying weights using ggplot 2                           ####


cols <- hcl.colors(n = length(wseq), palette = "viridis", alpha = 0.9, rev = TRUE)


# Now let's create a data frame for ggplot
densities <- list()
for (i in 1:length(tr)) {
densities[[i]] <- sapply(X = wseq, FUN = function(w) {
  rmapPostFix(theta = thetaseq, tr = tr[i], sr = sr[i], to = to, so = so,
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






# Plotting the error bars and points
plot_HPDI_median_rep <- ggplot(HPDI_theta_w_median, aes(x = weight, y = median)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = replication), width = 0.05, size = 1.4) +
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
  scale_color_manual(values = c("#E69F00", "#009E20", "#0072B2")) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1))

# Display the plot
print(plot_HPDI_median_rep)

ggsave(filename = "plot_HPDI_median_rep.pdf",path = "Plots", plot = plot_HPDI_median_rep,
       width = 15, height = 8, device='pdf', dpi=500, useDingbats = FALSE)



#   ____________________________________________________________________________
#   Tipping Point Analysis                                                  ####
