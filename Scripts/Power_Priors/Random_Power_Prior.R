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
#   Parameter Setting                                                       ####

# Original and Replicated Studies

load("credentials_data.RData")

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
rnumber <- c(1, 2, 3, 4)

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


#   ____________________________________________________________________________
#   Joint Posterior                                                         ####

joint_plot_df <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  pp_joint_post <- postPP(
    theta = par_grid$theta, alpha = par_grid$alpha, tr = tr[i],
    sr = sr[i], to = to, so = so, x = x, y = y, m = m, v = v
  )
  par_grid$density <- pp_joint_post
  par_grid$tr <- tr[i]
  par_grid$sr <- sr[i]
  par_grid$rnumber <- rnumber[i]
  return(par_grid)
}))

# Create a new column for ordering
joint_plot_df$rep_order <- ifelse(joint_plot_df$rnumber == 4, "p", as.character(joint_plot_df$rnumber))

# Create the original rep_setting column with proper labels
joint_plot_df$rep_setting <- paste0(
  "{hat(theta)[r*",
  ifelse(joint_plot_df$rnumber == 4, "p", joint_plot_df$rnumber),
  "] == ",
  round(joint_plot_df$tr, 2),
  "}*',' ~ sigma[r*",
  ifelse(joint_plot_df$rnumber == 4, "p", joint_plot_df$rnumber),
  "] == ",
  round(joint_plot_df$sr, 2)
)

# Convert the new rep_order column to a factor to specify the order
joint_plot_df$rep_order <- factor(joint_plot_df$rep_order, levels = c("1", "2", "3", "p"))


#   ____________________________________________________________________________
#   Marginal Posteriors                                                     ####

# Marginal posterior of alpha
marg_alpha_dens <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  pp_joint_post <- postPPalpha(
    alpha = alphaseq, tr = tr[i], sr = sr[i], to = to,
    so = so, x = x, y = y, m = m, v = v
  )
  out <- data.frame(
    x = alphaseq, density = pp_joint_post, rnumber = rnumber[i],
    parameter = "'Power parameter' ~ alpha", tr = tr[i], sr = sr[i]
  )
  return(out)
}))

# Marginal posterior of theta
marg_theta_dens <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  pp_joint_post <- postPPtheta(
    theta = thetaseq, tr = tr[i], sr = sr[i], to = to,
    so = so, x = x, y = y, m = m, v = v
  )
  out <- data.frame(
    x = thetaseq, density = pp_joint_post, rnumber = rnumber[i],
    parameter = "'Effect size' ~ theta", tr = tr[i], sr = sr[i]
  )
  return(out)
}))

# Both marginal distributions
marg_plot_df <- rbind(marg_alpha_dens, marg_theta_dens)

# Create a new column for ordering
marg_plot_df$rep_order <- ifelse(marg_plot_df$rnumber == 4, "p", as.character(marg_plot_df$rnumber))

# Create the original rep_setting column with proper labels
marg_plot_df$rep_setting <- paste0(
  "{hat(theta)[r*",
  ifelse(marg_plot_df$rnumber == 4, "p", marg_plot_df$rnumber),
  "] == ",
  round(marg_plot_df$tr, 2),
  "}*',' ~ sigma[r*",
  ifelse(marg_plot_df$rnumber == 4, "p", marg_plot_df$rnumber),
  "] == ",
  round(marg_plot_df$sr, 2)
)

# Posterior of effect size without using original data
marg_theta_dens2 <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  pp_joint_post <- dnorm(x = thetaseq, mean = tr[i], sd = sr[i])
  out <- data.frame(
    x = thetaseq, density = pp_joint_post, rnumber = rnumber[i],
    parameter = "'Effect size' ~ theta", tr = tr[i], sr = sr[i]
  )
  return(out)
}))

# Create a new column for ordering
marg_theta_dens2$rep_order <- ifelse(marg_theta_dens2$rnumber == 4, "p", as.character(marg_theta_dens2$rnumber))

# Create the original rep_setting column with proper labels
marg_theta_dens2$rep_setting <- paste0(
  "{hat(theta)[r*",
  ifelse(marg_theta_dens2$rnumber == 4, "p", marg_theta_dens2$rnumber),
  "] == ",
  round(marg_theta_dens2$tr, 2),
  "}*',' ~ sigma[r*",
  ifelse(marg_theta_dens2$rnumber == 4, "p", marg_theta_dens2$rnumber),
  "] == ",
  round(marg_theta_dens2$sr, 2)
)

# 95% HPDI
HPDI_alpha <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  hpd <- postPPalphaHPD(
    level = 0.95, tr = tr[i], sr = sr[i], to = to,
    so = so, x = x, y = y, m = m, v = v
  )
  out <- data.frame(
    y = max(marg_alpha_dens$density) * (1 + 0.05 * i),
    lower = hpd[1], upper = hpd[2], rnumber = rnumber[i],
    parameter = "'Power parameter' ~ alpha", tr = tr[i],
    sr = sr[i], height = 0.2
  )
  return(out)
}))
HPDI_theta <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  hpd <- postPPthetaHPD(
    level = 0.95, tr = tr[i], sr = sr[i], to = to,
    so = so, x = x, y = y, m = m, v = v
  )
  out <- data.frame(
    y = max(c(marg_theta_dens2$density, marg_theta_dens$density)) * (1 + 0.06 * i),
    lower = hpd[1], upper = hpd[2], rnumber = rnumber[i],
    parameter = "'Effect size' ~ theta", tr = tr[i],
    sr = sr[i], height = 0.6
  )
  return(out)
}))
HPDI_df <- rbind(HPDI_alpha, HPDI_theta)

# Create a new column for ordering
HPDI_df$rep_order <- ifelse(HPDI_df$rnumber == 4, "p", as.character(HPDI_df$rnumber))

# Create the original rep_setting column with proper labels
HPDI_df$rep_setting <- paste0(
  "{hat(theta)[r*",
  ifelse(HPDI_df$rnumber == 4, "p", HPDI_df$rnumber),
  "] == ",
  round(HPDI_df$tr, 2),
  "}*',' ~ sigma[r*",
  ifelse(HPDI_df$rnumber == 4, "p", HPDI_df$rnumber),
  "] == ",
  round(HPDI_df$sr, 2)
)

HPDI_theta_2 <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  hpd <- tr[i] + c(-1, 1) * qnorm(p = 0.975) * sr[i]
  out <- data.frame(
    y = max(c(marg_theta_dens2$density, marg_theta_dens$density)) * (1 + 0.05 * i),
    lower = hpd[1], upper = hpd[2], rnumber = rnumber[i],
    parameter = "'Effect size' ~ theta", tr = tr[i],
    sr = sr[i], height = 0.6
  )
  return(out)
}))

# Create a new column for ordering
HPDI_theta_2$rep_order <- ifelse(HPDI_theta_2$rnumber == 4, "p", as.character(HPDI_theta_2$rnumber))

# Create the original rep_setting column with proper labels
HPDI_theta_2$rep_setting <- paste0(
  "{hat(theta)[r*",
  ifelse(HPDI_theta_2$rnumber == 4, "p", HPDI_theta_2$rnumber),
  "] == ",
  round(HPDI_theta_2$tr, 2),
  "}*',' ~ sigma[r*",
  ifelse(HPDI_theta_2$rnumber == 4, "p", HPDI_theta_2$rnumber),
  "] == ",
  round(HPDI_theta_2$sr, 2)
)

## Limitting density for perfectly agreeing effect estimates with c = so^2/sr^2 -> infty
alpha_limit_df <- data.frame(
  x = alphaseq,
  density = dbeta(x = alphaseq, x + 0.5, y),
  parameter = "'Power parameter' ~ alpha"
)


#   ____________________________________________________________________________
#   Bayes Factor                                                            ####

# Parameters for Bayes factors
k <- sqrt(2) # unit-information standard deviation
x <- 1 # uniform prior for effect size BF
y <- 1 # uniform prior for effect size BF
yd <- 2 # monotonically decreasing prior for power parameter BF


# Function to format Bayes factors

format_bf <- function(BF, digits = "default") {
  ## check inputs
  stopifnot(
    length(BF) == 1,
    is.numeric(BF),
    (is.finite(BF) && 0 < BF) || is.na(BF),
    length(digits) == 1,
    (is.character(digits) && digits == "default") ||
      (is.numeric(digits) && 0 <= digits)
  )
  ## return NA if input NA/NaN
  if (is.na(BF) || is.nan(BF)) {
    result <- NA
  } else {
    ## format BF
    if (digits == "default") {
      if (BF < 1 / 1000) {
        result <- "< 1/1000"
      }
      if ((BF >= 1 / 1000) & (BF <= 1 / 10)) {
        result <- paste0("1/", as.character(round(1 / BF)))
      }
      if ((BF > 1 / 10) & (BF < 1)) {
        result <- paste0("1/", as.character(round(1 / BF, digits = 1)))
      }
      if ((BF < 10) & (BF >= 1)) {
        result <- as.character(round(BF, digits = 1))
      }
      if ((BF >= 10) & (BF <= 1000)) {
        result <- as.character(round(BF))
      }
      if (BF > 1000) {
        result <- "> 1000"
      }
    } else {
      if (BF < 1) {
        result <- paste0("1/", as.character(round(1 / BF, digits = digits)))
      } else {
        result <- as.character(round(BF, digits = digits))
      }
    }
    ## when 1/1 return 1
    if (result == "1/1") result <- "1"
  }
  return(result)
}
format_bf_vec <- Vectorize(FUN = format_bf)

rnumber <- c(1, 2, 3, 4)

# Compute BFs for effect sizes and power parameter
bf_df <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  bf_theta_random <- bfPPtheta(tr = tr[i], sr = sr[i], to = to, so = so, x = x, y = y)
  bf_theta <- bfPPtheta(tr = tr[i], sr = sr[i], to = to, so = so, alpha = 1)
  bf_alpha <- bfPPalpha(tr = tr[i], sr = sr[i], to = to, so = so, uv = k^2)
  bf_random_alpha <- bfPPalpha(tr = tr[i], sr = sr[i], to = to, so = so, x = 1, y = yd)
  out <- data.frame(
    number = rnumber[i], tr = tr[i], sr = sr[i],
    bf_theta_random = bf_theta_random,
    bf_theta = bf_theta, bf_alpha = bf_alpha, bf_random_alpha = bf_random_alpha
  )
  return(out)
}))

## Create LaTeX table for theta
dfTab_theta <- bf_df[, 1:5] %>%
  mutate(
    bf_theta = format_bf_vec(bf_theta),
    bf_theta_random = format_bf_vec(bf_theta_random),
    tr = round(tr, 2),
    sr = round(sr, 2),
    number = as.integer(number)
  ) %>%
  arrange(number)
xtab_theta <- xtable(dfTab_theta)
colnames(xtab_theta) <- c(
  "",
  "$\\hat{\\theta}_r$",
  "$\\sigma_r$",
  paste0(
    "$\\mathrm{BF}_{01}\\{\\hat{\\theta}_r \\mid \\mathcal{H}_{1} : \\alpha \\sim \\mathrm{Beta}(",
    x, ", ", y, ")\\}$"
  ),
  "$\\mathrm{BF}_{01}(\\hat{\\theta}_r \\mid \\mathcal{H}_{1} : \\alpha = 1)$"
)
align(xtab_theta) <- rep("c", length(colnames(xtab_theta)) + 1)
# Add multicolumns for effet size test and power parameter test
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- "\\toprule"

print(xtab_theta,
  floating = FALSE, include.rownames = FALSE, add.to.row = addtorow,
  sanitize.text.function = function(x) {
    x
  }, booktabs = TRUE, hline.after = c(0, nrow(xtab_theta))
)



# Create LaTeX table for alpha
dfTab_alpha <- bf_df[, c(1:3, 6:7)] %>%
  mutate(
    bf_alpha = format_bf_vec(bf_alpha),
    bf_random_alpha = format_bf_vec(bf_random_alpha),
    tr = round(tr, 2),
    sr = round(sr, 2),
    number = as.integer(number)
  ) %>%
  arrange(number)
xtab_alpha <- xtable(dfTab_alpha)
colnames(xtab_alpha) <- c(
  "",
  "$\\hat{\\theta}_r$",
  "$\\sigma_r$",
  paste0(
    "$\\mathrm{BF}_{\\text{dc}}(\\hat{\\theta}_r \\mid \\mathcal{H}_d : \\alpha = ",
    0, ")$"
  ),
  paste0(
    "$\\mathrm{BF}_{\\text{dc}}\\{\\hat{\\theta}_r \\mid \\mathcal{H}_d : \\alpha \\sim \\mathrm{Beta}(",
    1, ", ", 2, ")\\}$"
  )
)
align(xtab_alpha) <- rep("c", length(colnames(xtab_alpha)) + 1)

# Add multicolumns for effet size test and power parameter test
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- "\\toprule"

print(xtab_alpha,
  floating = FALSE, include.rownames = FALSE, add.to.row = addtorow,
  sanitize.text.function = function(x) {
    x
  }, booktabs = TRUE, hline.after = c(0, nrow(xtab_alpha))
)
