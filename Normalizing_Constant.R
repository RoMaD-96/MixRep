#   ____________________________________________________________________________
#   Libraries                                                               ####

library(ggplot2)
library(ggpubr)
library(colorspace)


#   ____________________________________________________________________________
#   Fixed Weights                                                           ####

normConst <- function(tr, sr, to, so, null, priorsd, x, y) {
  intFun <- function(w) {
    (w*dnorm(x = tr, mean = to, sd = sqrt(so^2 + sr^2)) +
       (1 - w)*dnorm(x = tr, mean = null, sd = sqrt(priorsd^2 + sr^2)))*dbeta(x = w, shape1 = x, shape2 = y)
  }
  res <- try(integrate(f = intFun, lower = 0, upper = 1)$value)
  if (inherits(res, "try-error")) {
    out <- NaN
  } else {
    out <- res
  }
  return(out)
}


## posterior of theta conditional on w
rmapPost <- function(theta, w,  tr, sr, to, so, null, priorsd, x, y) {
  v1 <- 1/(1/so^2 + 1/sr^2)
  m1 <- (to/so^2 + tr/sr^2)*v1
  v2 <- 1/(1/priorsd^2 + 1/sr^2)
  m2 <- (null/priorsd^2 + tr/sr^2)*v2
  bf <- dnorm(x = tr, mean = to, sd = sqrt(sr^2 + so^2))/
    dnorm(x = tr, mean = null, sd = sqrt(sr^2 + priorsd^2))
  wUpdate <- 1/(1 + (1/w - 1)/bf)                     # Update based on the BF ?
  wUpdate*dnorm(x = theta, mean = m1, sd = sqrt(v1)) +
    (1 - wUpdate)*dnorm(x = theta, mean = m2, sd = sqrt(v2))
}

to <- 0.21
so <- 0.05
tr <- 0.09
sr <- 0.045
null <- 0
priorsd <- 2
wseq <- seq(0, 1, 0.1)
thetaseq <- seq(-0.2, 0.4, length.out = 500)
postdens <- sapply(X = wseq, FUN = function(w) {
  rmapPost(theta = thetaseq, tr = tr, sr = sr, to = to, so = so,
           null = null, priorsd = priorsd, w = w)})


# intFun <- function(x) {
#     rmapPost(theta = x, tr = tr, sr = sr, to = to, so = so,
#              null = null, priorsd = priorsd, w = w)
# }
# integrate(f = intFun, lower = -Inf, upper = Inf)

cols <- hcl.colors(n = length(wseq), palette = "Viridis", alpha = 0.7, rev = TRUE)
matplot(thetaseq, postdens, type = "l", xlab = bquote(theta), ylab = "Density",
        las = 1, lty = 1, lwd = 1.5, col = cols)
lines(thetaseq, dnorm(x = tr, mean = thetaseq, sd = sr), lty = 2, lwd = 1.5)
lines(thetaseq, dnorm(x = thetaseq, mean = to, sd = so), lty = 3, lwd = 1.5)
lines(thetaseq, dnorm(x = thetaseq, mean = null, sd = priorsd), lty = 4, lwd = 1.5)
legend("topleft", c("Posterior", "Likelihood (Replication)",
                    "Prior (Original component)", "Prior (Robust component)"),
       lty = c(1, 2, 3, 4), lwd = 1.5)
legend("topright", legend = rev(wseq), lty = 1, col = rev(cols), title="weight")


##  ............................................................................
##  Using ggplot 2                                                          ####


cols <- hcl.colors(n = length(wseq), palette = "viridis", alpha = 0.9, rev = TRUE)


# Now let's create a data frame for ggplot
densities <- sapply(X = wseq, FUN = function(w) {
  rmapPost(theta = thetaseq, tr = tr, sr = sr, to = to, so = so,
           null = null, priorsd = priorsd, w = w)
})

# Turn the results into a data frame
df <- data.frame(theta = rep(thetaseq, times = length(wseq)),
                 density = as.vector(densities),
                 w = factor(rep(wseq, each = length(thetaseq))))

# Add the additional lines
df$likelihood <- dnorm(x = df$theta, mean = tr, sd = sr)
df$prior_original <- dnorm(x = df$theta, mean = to, sd = so)
df$prior_robust <- dnorm(x = df$theta, mean = null, sd = priorsd)


# Create a separate data frame for the additional lines to help in creating the legend
additional_lines <- data.frame(
  theta = rep(thetaseq, 3),
  value = c(dnorm(x = thetaseq, mean = tr, sd = sr),
            dnorm(x = thetaseq, mean = to, sd = so),
            dnorm(x = thetaseq, mean = null, sd = priorsd)),
  linetype = factor(rep(c("Likelihood", "Prior Original", "Prior Robust"), each = length(thetaseq)))
)

# The ggplot
plot_post_fix <- ggplot() + 
  geom_line(data = df, aes(x = theta, y = density, color = w), size = 1) + 
  geom_line(data = additional_lines, aes(x = theta, y = value, linetype = linetype), size = 0.8) +
  scale_color_manual(values = cols) +
  scale_linetype_manual(values = c("dashed", "dotted", "dotdash"),
                        labels = c("Likelihood (Replication)", "Prior (Original component)", "Prior (Robust component)")) +
  labs(x = expression("Effect Size" ~ theta), y = "Density") +
  theme_bw() +
  guides(linetype = guide_legend(title = "Density"),
         color = guide_legend(title = "Weight")
         )

# To make sure that our additional lines are represented in the legend, we need to add them to the plot
plot_post_fix <- plot_post_fix + geom_line(aes(linetype = "Likelihood (Replication)"), linetype = "dashed", color = "black") +
  geom_line(aes(linetype = "Prior (Original component)"), linetype = "dotted", color = "black") +
  geom_line(aes(linetype = "Prior (Robust component)"), linetype = "dotdash", color = "black") +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_text(size = 22),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 19))

ggsave(filename = "plot_post_fix.pdf",path = "Plots", plot = plot_post_fix,
       width = 16, height = 8.5, device='pdf', dpi=500, useDingbats = FALSE)
#   ____________________________________________________________________________
#   Random Weights                                                          ####


##  ............................................................................
##  Normalizing Constant Function                                           ####

normConst <- function(tr, sr, to, so, null, priorsd, x, y) {
  intFun <- function(w) {
    (w*dnorm(x = tr, mean = to, sd = sqrt(so^2 + sr^2)) +
       (1 - w)*dnorm(x = tr, mean = null, sd = sqrt(priorsd^2 + sr^2)))*dbeta(x = w, shape1 = x, shape2 = y)
  }
  res <- try(integrate(f = intFun, lower = 0, upper = 1)$value)
  if (inherits(res, "try-error")) {
    out <- NaN
  } else {
    out <- res
  }
  return(out)
}


##  ............................................................................
##  Posterior Density                                                       ####

rmapPost <- function(theta, w,  tr, sr, to, so, null, priorsd, x, y) {
    num_post <- dnorm( x = tr, mean = theta, sd =  sr)*(w * dnorm( x = theta, mean = to, sd = so) +
      (1 - w) * dnorm(x = theta, mean = null, sd = priorsd)) * dbeta(x = w, shape1 = x, shape2 = y)
    den_post <- normConst(tr, sr, to, so, null, priorsd, x, y)
    joint_post <- num_post/den_post
    return(joint_post)
}




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

m_post_weights <- function(w, tr, sr, to, so, null, priorsd, x, y) {
  num_post_weights <-  (w*dnorm(x = tr, mean = to, sd = sqrt(so^2 + sr^2)) +
                       (1 - w)*dnorm(x = tr, mean = null, sd = sqrt(priorsd^2 + sr^2)))*dbeta(x = w, shape1 = x, shape2 = y)
  den_post_weights <- normConst(tr, sr, to, so, null, priorsd, x, y)
  m_post <- num_post_weights/den_post_weights
  return(m_post)
}

# m_postdens_w_wrapper <- do.call("rbind", lapply(X = seq(1:length(tr)), FUN = function(index) {
#   m_post_dens <- m_post_weights(
#     w = par_grid$omega,
#     tr = tr[index],
#     sr = sr[index],
#     to = to,
#     so = so,
#     null = null,
#     priorsd = priorsd,
#     x = alpha,
#     y = beta
#   )
#   par_grid$density_weights <- m_post_dens
#   par_grid$rep_number <- index
#   return(par_grid)
# }))


rep_number <- c(1,2,3)

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
# 
# weights_m_post$rep_setting <-paste0( "{hat(theta)[italic('r')*",
#                                        postdens_wrapper$rep_number,
#                                        "] == ",
#                                        round(postdens_wrapper$tr, 2),
#                                        "}*',' ~ sigma[italic('r')*",
#                                        postdens_wrapper$rep_number,
#                                        "] == ",
#                                        round(postdens_wrapper$sr, 2)
# )


ggplot(data=weights_m_post, aes(x=x, y=density, group=rep_exp, color=factor(rep_exp))) +
  geom_line() +
  scale_color_manual(
    values = c("1" = "#E69F00", "2" = "#009E20", "3" = "#0072B2"),
    labels = c(expression(" "~hat(theta)[r*1] == 0.09 ~ ", " ~ sigma[r*1] == 0.04),
               expression(" "~hat(theta)[r * 2] == 0.21 ~ ", " ~ sigma[r*2]^2 == 0.06),
               expression(" "~hat(theta)[r * 3] == 0.44 ~ ", " ~ sigma[r*3]^2 == 0.04)))+
  labs(
    subtitle = "Marginal Posterior Densities of Weight Parameter",
    x = expression(omega~" Values"),
    y = "Density"
  ) +
  theme_minimal() +
  guides(color=guide_legend(title="Replicated Experiment")) 




##  ............................................................................
##  Marginal Posterior for theta                                            ####


m_post_theta <- function(theta, tr, sr, to, so, null, priorsd, x, y) {
    mean_beta <- x / (x + y)
    num_post_weights <- dnorm(x = tr, mean = theta, sd = sr) * (mean_beta * 
        (dnorm(x = theta, mean = to, sd = so) - dnorm(x = theta, mean = null, sd = priorsd))
       + dnorm(x = theta, mean = null, sd = priorsd))
    den_post_weights <- normConst(tr, sr, to, so, null, priorsd, x, y)
    m_post <- num_post_weights / den_post_weights
    return(m_post)
  }




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
  out <- data.frame(x = thetaseq, density = marg_p_dens, rep_exp = rep_number[index],
                    parameter = "'Weight parameter' ~ omega", tr = tr[index], sr = sr[index])
  return(out)
}))


ggplot(data=theta_m_post, aes(x=x, y=density, group=rep_exp, color=factor(rep_exp))) +
  geom_line(size = 1.0) +
  scale_color_manual(
    values = c("1" = "#E69F00", "2" = "#009E20", "3" = "#0072B2"),
    labels = c(expression(" "~hat(theta)[r*1] == 0.09 ~ ", " ~ sigma[r*1] == 0.04),
               expression(" "~hat(theta)[r * 2] == 0.21 ~ ", " ~ sigma[r*2]^2 == 0.06),
               expression(" "~hat(theta)[r * 3] == 0.44 ~ ", " ~ sigma[r*3]^2 == 0.04)))+
  labs(
    subtitle = "Marginal Posterior Densities of Weight Parameter",
    x = expression("Effect Size"~theta),
    y = "Marginal Posterior Density"
  ) +
  theme_minimal() +
  guides(color=guide_legend(title="Replicated Experiment")) 

