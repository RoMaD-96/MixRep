#   ____________________________________________________________________________
#   Libraries                                                               ####

library(bayesmeta)
library(dplyr)
library(knitr)
library(xtable)
library(repmix)

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

# Define labels for the replication sites
rep_names <- c("University of Toronto", "Montana State University", "Ashland University", "Pooled")
bf_list <- list()
models <- list()
for (i in seq_along(tr)) {
  y_pair <- c(to, tr[i])
  sigma_pair <- c(so, sr[i])

  # Bayesian two-study "meta-analysis" (replicated and original study)
  result_pair <- bayesmeta(
    y = y_pair,
    sigma = sigma_pair,
    mu_prior = c("mean" = 0, "sd" = 0.1),
    labels = c("Original", rep_names[i]),
    tau.prior = function(t) dhalfnormal(t, scale = 0.5) # Half-normal prior for heterogeneity
    # tau.prior = function(x){return(dhalfcauchy(x, scale=0.5))}
    # tau.prior = function(t) dinvgamma(t, shape =  0.5, rate =  0.5, log = FALSE)  # Inverse gamma prior for heterogeneity
    # tau.prior = "Jeffreys"
  )

  print(result_pair)
  models[[i]] <- result_pair
  bf_list[[i]] <- result_pair$bayesfactor
  # plot(result_pair)
}

# plot(thetaseq, models[[4]]$dposterior(mu = thetaseq))
#
# plot(models[[4]], which = 3)
#
# plot(range(thetaseq), c(0,10), type="n", xlab=expression(thetaseq[i]), ylab="")
# for (i in 1:models[[2]]$k) {
#   # draw effect's posterior distribution:
#   lines(thetaseq, models[[2]]$dposterior(theta=thetaseq, indiv=2), col=1, lty="solid")
#   lines(thetaseq, models[[1]]$dposterior(theta=thetaseq, indiv=2), col=2, lty="solid")
#   lines(thetaseq, models[[3]]$dposterior(theta=thetaseq, indiv=2), col=3, lty="solid")
#   lines(thetaseq, models[[4]]$dposterior(theta=thetaseq, indiv=2), col=4, lty="solid")
# }
# abline(h=0)
# legend(max(thetaseq), 10, legend=models[[1]]$label, col=(1:models[[1]]$k)+1, pch=15, xjust=1, yjust=1)

