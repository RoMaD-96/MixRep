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

tp <- round(sum(trep/srep^2)/sum(1/srep^2),2)
sp <- round(sqrt(1/sum(1/srep^2)),2)

tr <- c(trep,tp)
sr <- c(srep,sp)

# Define labels for the replication sites
rep_names <- c("University of Toronto", "Montana State University", "Ashland University")
bf_list <- list()
for(i in seq_along(trep)){

  y_pair <- c(to, trep[i])
  sigma_pair <- c(so, srep[i])
  
  # Bayesian two-study "meta-analysis" (replicated and original study)
  result_pair <- bayesmeta(
    y = y_pair,
    sigma = sigma_pair,
    labels = c("Original", rep_names[i]),
    # mu.prior.mean = 0,      # Vague prior centered at 0
    # mu.prior.sd = 2,        # unit information
    # tau.prior = function(t) dhalfnormal(t, scale = 0.5)  # Half-normal prior for heterogeneity
    # tau.prior = function(t) dinvgamma(t, shape =  0.5, rate =  0.5, log = FALSE)  # Inverse gamma prior for heterogeneity
    tau.prior = "Jeffreys"
  )

  print(result_pair)
  bf_list[[i]] <- result_pair$bayesfactor 
  #plot(result_pair)
}

