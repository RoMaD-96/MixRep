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



## ----"BF-parameters"----------------------------------------------------------
## Parameters for Bayes factors
k <- sqrt(2) # unit-information standard deviation
x <- 1 # uniform prior for effect size BF
y <- 1 # uniform prior for effect size BF
yd <- 2 # monotonically decreasing prior for power parameter BF


## ----"table-Bayes-factors", results = "asis"----------------------------------
## Function to nicely format Bayes factors
.formatBF_ <- function(BF, digits = "default") {
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
  if (is.na(BF) || is.nan(BF))
    result <- NA
  else {
    ## format BF
    if (digits == "default") {
      if (BF < 1/1000)
        result <- "< 1/1000"
      if ((BF >= 1/1000) & (BF <= 1/10))
        result <- paste0("1/", as.character(round(1/BF)))
      if ((BF > 1/10) & (BF < 1))
        result <- paste0("1/", as.character(round(1/BF, digits = 1)))
      if ((BF < 10) & (BF >= 1))
        result <- as.character(round(BF, digits = 1))
      if ((BF >= 10) & (BF <= 1000))
        result <- as.character(round(BF))
      if (BF > 1000)
        result <- "> 1000"
    } else {
      if (BF < 1)
        result <- paste0("1/", as.character(round(1/BF, digits = digits)))
      else
        result <- as.character(round(BF, digits = digits))
    }
    ## when 1/1 return 1
    if (result == "1/1") result <- "1"
  }
  return(result)
}
formatBF <- Vectorize(FUN = .formatBF_)

## Compute BFs for effect sizes and power parameter
bfDF <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  bf01 <- bfPPtheta(tr = tr[i], sr = sr[i], to = to, so = so, x = x, y = y)
  bfr <- bfPPtheta(tr = tr[i], sr = sr[i], to = to, so = so, alpha = 1)
  bfdc <- bfPPalpha(tr = tr[i], sr = sr[i], to = to, so = so, uv = k^2)
  bfdcRandom <- bfPPalpha(tr = tr[i], sr = sr[i], to = to, so = so, y = yd)
  out <- data.frame(number = rnumber[i], tr = tr[i], sr = sr[i], bf = bf01,
                    bfr = bfr, bfdc = bfdc, bfdcRandom = bfdcRandom)
  return(out)
}))

## Create LaTeX table
dfTab <- bfDF %>%
  mutate(bf = formatBF(bf),
         bfr = formatBF(bfr),
         bfdc = formatBF(bfdc),
         bfdcRandom = formatBF(bfdcRandom),
         tr = round(tr, 2),
         sr = round(sr, 2),
         number = as.integer(number)) %>%
  arrange(tr)
xtab <- xtable(dfTab)
colnames(xtab) <- c("",
                    "$\\hat{\\theta}_r$",
                    "$\\sigma_r$",
                    paste0("$\\BF_{01}\\{\\hat{\\theta}_r \\given \\h{1} \\colon \\alpha \\sim \\Be(",
                           x, ", ", y, ")\\}$"),
                    "$\\BF_{01}(\\hat{\\theta}_r \\given \\h{1} \\colon \\alpha = 1)$",
                    ## paste0("$\\BF_{\\text{dc}}(\\hat{\\theta}_r \\given \\kappa^2 =",
                    ##        k^2, ")$"),
                    paste0("$\\BF_{\\text{dc}}(\\hat{\\theta}_r \\given \\h{\\text{d}} \\colon \\alpha = 0)$"),
                    paste0("$\\BF_{\\text{dc}}\\{\\hat{\\theta}_r \\given \\h{\\text{d}} \\colon \\alpha \\sim \\Be(1, ",
                           yd, ")\\}$"))
align(xtab) <- rep("c", length(colnames(xtab)) + 1)

## add multicolumns for effet size test and power parameter test
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\toprule & & & \\multicolumn{2}{c}{Tests about the effect size $\\theta$} & \\multicolumn{2}{c}{Tests about the power parameter $\\alpha$} \\\\ \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}'


print(xtab, floating = FALSE, include.rownames = FALSE, add.to.row = addtorow,
      sanitize.text.function = function(x){x}, booktabs = TRUE, hline.after = c(0, nrow(xtab)))


## ----"Fig2-BF-sensitivity", fig.height = 3.25---------------------------------
## Compute BFs for a grid of prior parameters to see how robust results are
yseq <- xseq <- exp(seq(log(1.001), log(100), 0.1))
bfDF2 <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  ## simple vs. composite BF
  bfdcRandom <- sapply(X = yseq, FUN = function(y) {
    bfPPalpha(tr = tr[i], sr = sr[i], to = to, so = so, y = y)
  })
  ## composite vs. composite BF
  bfdcRandomFull <- sapply(X = xseq, FUN = function(x) {
    fd <- margLik(tr = tr[i], sr = sr[i], to = to, so = so, x = 1, y = 2)
    fc <- margLik(tr = tr[i], sr = sr[i], to = to, so = so, x = x, y = 1)
    fd/fc
  })
  out <- data.frame(number = rnumber[i], tr = tr[i], sr = sr[i],
                    bfdcRandom = bfdcRandom, bfdcRandomFull = bfdcRandomFull,
                    y = yseq)
  return(out)
}))
bfDF2$rep_setting <- paste0("{hat(theta)[italic('r')*", bfDF2$number, "] == ",
                         round(bfDF2$tr, 2), "}*',' ~ sigma[italic('r')*",
                         bfDF2$number, "] == ", round(bfDF2$sr, 2))

## plot BFs as a function of prior parameters
bfbks <- c(1/100, 1/30, 1/10, 1/3, 1, 3, 10, 30, 100)
bflabs <- formatBF(BF = bfbks)
ploty <- ggplot(data = bfDF2, aes(x = y, y = bfdcRandom, color = rep_setting)) +
  annotate(geom = "segment", x = 0.95, xend = 0.95, y = 1.1, yend = 11,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")), alpha = 0.95,
           color = "darkgrey") +
  annotate(geom = "segment", x = 0.95, xend = 0.95, y = 1/1.1, yend = 1/11,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")), alpha = 0.95,
           color = "darkgrey") +
  annotate(geom = "text", x = 0.88, y = 3.5, label = "incompatible", angle = 90,
           alpha = 0.7, size = 2.5) +
  annotate(geom = "text", x = 0.88, y = 1/3.5, label = "compatible", angle = 90,
           alpha = 0.7, size = 2.5) +
  geom_hline(yintercept = 1, lty = 2, alpha = 0.2) +
  geom_line(linewidth = 0.7) +
  scale_y_log10(breaks = bfbks, labels = bflabs) +
  scale_x_log10(breaks = bfbks) +
  coord_cartesian(xlim = c(min(yseq), max(yseq)), ylim = c(1/45, 45)) +
  scale_color_manual(values = colblind, labels = scales::parse_format()) +
  labs(x = bquote("Prior parameter" ~ italic(y)),
       y = bquote("BF"["dc"] * "{"* hat(theta)[italic("ri")] ~ "|" * italic("H")["d"] ~
                    ":" ~ alpha ~
                    "~ Be(1," * italic(y) * ")}"),
       color = "") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top")
plotx <- ggplot(data = bfDF2, aes(x = y, y = bfdcRandomFull, color = rep_setting)) +
  annotate(geom = "segment", x = 0.95, xend = 0.95, y = 1.1, yend = 11,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")), alpha = 0.95,
           color = "darkgrey") +
  annotate(geom = "segment", x = 0.95, xend = 0.95, y = 1/1.1, yend = 1/11,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")), alpha = 0.95,
           color = "darkgrey") +
  annotate(geom = "text", x = 0.88, y = 3.5, label = "incompatible", angle = 90,
           alpha = 0.7, size = 2.5) +
  annotate(geom = "text", x = 0.88, y = 1/3.5, label = "compatible", angle = 90,
           alpha = 0.7, size = 2.5) +
  geom_hline(yintercept = 1, lty = 2, alpha = 0.2) +
  geom_line(linewidth = 0.7) +
  scale_y_log10(breaks = bfbks, labels = bflabs) +
  scale_x_log10(breaks = bfbks) +
  coord_cartesian(xlim = c(min(yseq), max(yseq)), ylim = c(1/45, 45)) +
  scale_color_manual(values = colblind, labels = scales::parse_format()) +
  labs(x = bquote("Prior parameter" ~ italic(x)),
       y = bquote("BF"["dc"] * "("* hat(theta)[italic("ri")] ~ "|" ~ italic(x) * ")"),
       color = "") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top")
ploty
## ggpubr::ggarrange(ploty, plotx, common.legend = TRUE)


## ----"asymptotic-Bayes-factor"------------------------------------------------
## Compute asymptotic bounds for the Bayes factors
k <- sqrt(2)
g <- k^2/so^2
s <- g/(1 + g)
trueESmin <- to ## minimized for true effect size = original effect estimate
bfBound <- dnorm(x = trueESmin, mean = 0, sd = k)/
  dnorm(x = trueESmin, mean = s*to, sd = sqrt(s)*so)
bfBound2 <- beta(3/2, 2)/beta(1, 2)


## ----"Fig3-Bayes-factor-design", fig.height = 4.5-----------------------------
## Function to compute probability of replication success
powerFun <- function(sr, to, so, k, level, mi, vi) {
  s <- k^2 / (so^2 + k^2)
  X <- (sr^2 + k^2) * (sr^2 + s * so^2) / (k^2 - s * so^2) *
    (2 * log(level) - log((sr^2 + s * so^2) / (sr^2 + k^2)) -
       s^2*to^2/(s*so^2 - k^2))
  ## lambdai <- (mi - s * to * (sr^2 + k^2) / (k^2 - s * so^2))^2/vi
  lambdai <- (mi - to * (sr^2 + k^2) / k^2)^2/vi
  p <- stats::pchisq(q = X/vi, df = 1, ncp = lambdai, lower.tail = TRUE)
  return(p)
}

## Compute proability of replication success for a grid of relative
## variances (c = so^2/sr^2 =~ nr/no)
cSeq <- exp(seq(log(1/12), log(12), length.out = 1000))
cbks <- c(1/10, 1/3, 1, 3, 10)
level <- 1/10
plotDF <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  to <- tr[i]
  so <- sr[i]
  srSeq <- so / sqrt(cSeq)
  ## mean and variance under Hd and Hc
  md <- 0
  vd <- srSeq^2 + k^2
  s <- k^2/so^2 / (1 + k^2/so^2)
  mc <- s*to
  vc <- srSeq^2 + s*so^2
  ## power to achieve BF < level
  p1Hd <- powerFun(sr = srSeq, to = to, so = so, k = k, level = level,
                   mi = md, vi = vd)
  p1Hc <- powerFun(sr = srSeq, to = to, so = so, k = k, level = level,
                   mi = mc, vi = vc)
  outDF1 <- rbind(data.frame(p = p1Hd, level = level, type = "italic(H['d'])",
                             sr = srSeq, c = cSeq, to = to, so = so,
                             direction = "<="),
                  data.frame(p = p1Hc, level = level, type = "italic(H['c'])",
                             sr = srSeq, c = cSeq, to = to, so = so,
                             direction = "<="))
  ## power to achieve BF > 1/level
  p2Hd <- 1 - powerFun(sr = srSeq, to = to, so = so, k = k, level = 1/level,
                       mi = md, vi = vd)
  p2Hc <- 1 - powerFun(sr = srSeq, to = to, so = so, k = k, level = 1/level,
                       mi = mc, vi = vc)
  outDF2 <- rbind(data.frame(p = p2Hd, level = 1/level, type = "italic(H['d'])",
                             sr = srSeq, c = cSeq, to = to, so = so,
                             direction = ">="),
                  data.frame(p = p2Hc, level = 1/level, type = "italic(H['c'])",
                             sr = srSeq, c = cSeq, to = to, so = so,
                             direction = ">="))
  outDF <- rbind(outDF1, outDF2)
  outDF$yFacetLab <- paste0("'Pr(BF'['dc']", outDF$direction,
                            formatBF(outDF$level),
                            "~ '|' ~ italic(H['i']) *", "')'")
  outDF$xFacetLab <- paste0("{hat(theta)[italic(o)] ==",
                            round(outDF$to, 2),
                            "}*', '*sigma[italic(o)] == ",
                            round(outDF$so, 2))
  return(outDF)
}))

## Determine replication standard error such that certain power is achieved
pow <- 0.8
ssDF <- do.call("rbind", lapply(X = seq(1, length(tr)), FUN = function(i) {
  to <- tr[i]
  so <- sr[i]
  ## standard error to achieve P(BF < level | Hc, sr) = 0.8
  rootFunHc <- function(c) {
    sr <- so/sqrt(c)
    mc <- s*to
    vc <- sr^2 + s*so^2
    powerFun(sr = sr, to = to, so = so, k = k, level = level,
             mi = mc, vi = vc) - pow
  }
  resHd <- try(uniroot(f = rootFunHc, interval = c(1/100, 100))$root)
  if (class(resHd) == "try-error") {
    srHc <- NaN
  } else {
    srHc <- so/sqrt(resHd)
  }
  ## standard error to achieve P(BF > 1/level | Hd, sr) = 0.8
  rootFunHd <- function(c) {
    sr <- so/sqrt(c)
    md <- 0
    vd <- sr^2 + k^2
    (1 - powerFun(sr = sr, to = to, so = so, k = k, level = 1/level,
                  mi = md, vi = vd)) - pow
  }
  resHd <- try(uniroot(f = rootFunHd, interval = c(1/100, 100))$root)
  if (class(resHd) == "try-error") {
    srHd <- NaN
  } else {
    srHd <- so/sqrt(resHd)
  }
  outDF <- rbind(data.frame(level = 1/level, type = "italic(H['d'])",
                            sr = srHd, c = so^2/srHd^2, to = to, so = so,
                            direction = ">=", power = pow),
                 data.frame(level = level, type = "italic(H['c'])",
                            sr = srHc, c = so^2/srHc^2, to = to, so = so,
                            direction = "<=", power = pow))
  outDF$yFacetLab <- paste0("'Pr(BF'['dc']", outDF$direction,
                            formatBF(outDF$level),
                            "~ '|' ~ italic(H['i']) *", "')'")
  outDF$xFacetLab <- paste0("{hat(theta)[italic(o)] ==",
                            round(outDF$to, 2),
                            "}*', '*sigma[italic(o)] == ",
                            round(outDF$so, 2))
  return(outDF)
}))

## colorblind friendly scale
ncat <- length(unique(plotDF$type))
colblind <- palette.colors(n = ncat + 1, palette = "Okabe-Ito")[2:(ncat + 1)]
names(colblind) <- unique(plotDF$type)

## Plot power curves and required standard errors
powbks <- seq(0, 1, 0.2)
ggplot(data = plotDF, aes(x = c, y = p, color = type)) +
  facet_grid(yFacetLab ~ xFacetLab, labeller = label_parsed,
             switch = "y") +
  geom_hline(yintercept = pow, lty = 2, alpha = 0.1) +
  geom_line(alpha = 0.9, size = 0.7) +
  geom_point(data = ssDF, aes(x = c, y = power), size = 0.8,
             show.legend = FALSE) +
  geom_segment(data = ssDF, aes(x = c, xend = c, y = power, yend = 0),
               alpha = 0.3, arrow = arrow(length = unit(0.15, "cm")),
               show.legend = FALSE, size = 0.5) +
  labs(x = bquote("Relative variance" ~ sigma[italic(o)]^2/sigma[italic(r)]^2  %~~%
                    italic(n[r]) / italic(n[o])),
       y = NULL, color = NULL) +
  scale_y_continuous(breaks = powbks, labels = scales::percent,
                     limits = c(0, 1)) +
  ## scale_color_brewer(palette = "Dark2", labels = scales::parse_format()) +
  scale_color_manual(values = colblind, labels = scales::parse_format()) +
  scale_x_log10(breaks = cbks, labels = formatBF(cbks)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        strip.background.y = element_blank(), strip.placement = "outside")


## ----"Fig4-I2-alpha", fig.height = 3------------------------------------------
## Compute alpha from I^2 and I
I2seq <- seq(from = 0, to = 1, length.out = 1000)
alphaseq <- (1 - I2seq)/(1 + I2seq)
I2DF <- data.frame(I2 = I2seq, alpha = alphaseq, tau2 = so^2/2*(1/alphaseq - 1))
I2DF$x <- I2DF$I2
I2DF$type <- "italic('I')^2"
IDF <- I2DF
IDF$x <- sqrt(I2DF$x)
IDF$type <- "italic('I')"

## compute alpha from tau^2
tau2seq <- seq(0, 0.1, length.out = 1000)
tauDF <- IDF
tauDF$x <- tau2seq
tauDF$alpha <-  so^2/(2*tau2seq + so^2)
tauDF$type <- "tau^2"

## Plot alpha as a function of tau, I, and I^2
plotDF <- rbind(I2DF, IDF, tauDF)
plotDF$type <- factor(plotDF$type, levels = c("tau^2", "italic('I')^2", "italic('I')"))
ggplot(data = filter(plotDF, type != "italic('I')"),
       aes(x = x, y = alpha)) +
  ## geom_abline(intercept = 1, slope = -1, alpha = 0.1) +
  geom_line() +
  facet_wrap(~ type, labeller = label_parsed,
             strip.position = "bottom", scales = "free_x") +
  theme_bw() +
  labs(x = NULL, y = bquote(alpha)) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        strip.background.x = element_blank(), strip.placement = "outside",
        strip.text.x = element_text(size = 11))



## ----"Fig5-corresponding-priors", fig.height = 5------------------------------
## Function for computing prior density of tau^2 based on the corresponding
## Beta(a, b) prior on power parameter alpha
ftau2 <- function(tau2, a, b, so) {
  dbeta(x = so^2/(2*tau2 + so^2), shape1 = a, shape2 = b) *
    2*so^2/(2*tau2 + so^2)^2
}

## check that integrates to one
## integrate(f = ftau2, lower = 0, upper = Inf, a = 0.5, b = 2, so = 0.05)

## Function for computing prior density of I^2 based on the corresponding
## Beta(a, b) prior on power parameter alpha
fi2 <- function(i2, a, b) {
  dbeta(x = (1 - i2)/(1 + i2), shape1 = a, shape2 = b) *
    2 / (1 + i2)^2
  ## VGAM::dlino(x = i2, shape1 = b, shape2 = a, lambda = 2)
}

## check that integrates to one
## integrate(f = fi2, lower = 0, upper = 1, a = 0.5, b = 2)

## Compute prior density of tau^2 for different Beta priors on alpha
paramsGrid <- data.frame(a = c(1, 2, 1), b = c(1, 1, 2))
aseq <- seq(0, 1, length.out = 500)
tau2seq <- seq(0, 0.09, length.out = 500)^2
alphaDF <- do.call("rbind", lapply(X = seq(1, nrow(paramsGrid)), FUN = function(i) {
  a <- paramsGrid$a[i]
  b <- paramsGrid$b[i]
  dens <- dbeta(x = aseq, shape1 = a, shape = b)
  out <- data.frame(x = aseq, xlab = "alpha", density = dens, a = a, b = b,
                    ylab = paste0("{italic(x) ==", a, "}*','~ italic(y) ==", b))
  return(out)
}))
tau2DF <- do.call("rbind", lapply(X = seq(1, nrow(paramsGrid)), FUN = function(i) {
  a <- paramsGrid$a[i]
  b <- paramsGrid$b[i]
  dens <- ftau2(tau2 = tau2seq, a = a, b = b, so = so)
  out <- data.frame(x = tau2seq, xlab = "tau^2", density = dens, a = a, b = b,
                    ylab = paste0("{italic(x) ==", a, "}*','~ italic(y) ==", b))
  return(out)
}))
i2DF <- do.call("rbind", lapply(X = seq(1, nrow(paramsGrid)), FUN = function(i) {
  a <- paramsGrid$a[i]
  b <- paramsGrid$b[i]
  dens <- fi2(i2 = aseq, a = a, b = b)
  out <- data.frame(x = aseq, xlab = "italic('I')^2",
                    density = dens, a = a, b = b,
                    ylab = paste0("{italic(x) ==", a, "}*','~ italic(y) ==", b))
  return(out)
}))
plotDF <- rbind(alphaDF, tau2DF, i2DF)
alphaDF$ylabFac <- factor(x = alphaDF$ylab, levels = unique(plotDF$ylab))
tau2DF$ylabFac <- factor(x = tau2DF$ylab, levels = unique(plotDF$ylab))
i2DF$ylabFac <- factor(x = i2DF$ylab, levels = unique(plotDF$ylab))

## Plot beta priors for alpha and corresponding priors for tau^2
plotalpha <- ggplot(data = alphaDF, aes(x = x, y = density)) +
  facet_grid(ylabFac ~ ., labeller = label_parsed, scales = "free",
             switch = "x") +
  geom_line() +
  labs(x = bquote(alpha), y = NULL) +
  expand_limits(y = c(0, 3)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        strip.background.x = element_blank(), strip.placement = "outside")

ploti2 <- ggplot(data = i2DF, aes(x = x, y = density)) +
  facet_grid(ylabFac ~ ., labeller = label_parsed,
             switch = "x") +
  geom_line() +
  labs(x = bquote(italic("I")^2), y = NULL) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), strip.text.y = element_blank(),
        strip.background.x = element_blank(), strip.placement = "outside")

plottau <- ggplot(data = tau2DF, aes(x = x, y = density)) +
  facet_grid(ylabFac ~ ., labeller = label_parsed, scales = "free",
             switch = "x") +
  geom_line() +
  labs(x = bquote(tau^2), y = "Density") +
  ## scale_x_continuous(sec.axis = sec_axis(trans = ~ ./so^2,
  ##                                        name = bquote(tau^2/sigma[italic("o")]^2))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), strip.text.y = element_blank(),
        strip.background.x = element_blank(), strip.placement = "outside")


## Combine plots
ggpubr::ggarrange(plottau, ploti2, plotalpha, ncol = 3, align = "h",
                  widths = c(1.25, 1.05, 1.1))


## ----"Fig6-consistent-prior", fig.height = 5----------------------------------
## prior for alpha that leads to consistent test
falpha <- function(alpha, q, r, so) {
  r^q/gamma(x = q)*alpha^(q - 1)/(1 - alpha)^(q + 1)*(2/so^2)^q*
    exp(-2*r*alpha/so^2/(1 - alpha))
}
applyGrid <- expand.grid(q = c(0.5, 1, 2), r = c(0.5, 1, 2),
                         so = c(0.05, 0.5, 1))
aseq <- c(seq(0, 0.1, 0.00001), seq(0.101, 1, 0.001))
plotDF <- lapply(X = seq(1, nrow(applyGrid)), FUN = function(i) {
  r <- applyGrid$r[i]
  q <- applyGrid$q[i]
  so <- applyGrid$so[i]
  ## int <- integrate(f = falpha, lower = 0, upper = 1, q = q, r = r, so = so)$value
  data.frame(alpha = aseq,
             density = falpha(alpha = aseq, q = q, r = r, so = so),
             q = q, r = r, so = round(so, 2))
}) %>%
  bind_rows()

ggplot(data = plotDF, aes(x = alpha, y = density, color = factor(so))) +
  facet_grid(q ~ r, scales = "free",
             labeller = label_bquote(cols = italic(r) == .(r),
                                     rows = italic(q) == .(q))) +
  geom_line(alpha = 0.9, linewidth = 0.7) +
  coord_cartesian(ylim = c(0, 9)) +
  scale_color_viridis_d() +
  theme_bw() +
  labs(x = bquote("Power parameter" ~ alpha), y = "Density",
       color = bquote("Original standard error" ~ sigma[o])) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top")