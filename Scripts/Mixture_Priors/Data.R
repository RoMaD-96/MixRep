
#   ____________________________________________________________________________
#   Data                                                                    ####

# Replication data
rdat <- read.csv(file = "https://osf.io/download/36ed5/")

# Original data (see https://osf.io/g37tm)
 eta2o <- 0.043 # eta^2 effect size
 fiso <- asinh(sqrt(eta2o)) # convert to Fisher's z
 eta2oCIupper <- 0.103  # CI upper bound on eta^2 scale
 so <- (asinh(sqrt(eta2oCIupper)) - asinh(sqrt(eta2o)))/stats::qnorm(p = 0.975)

# Combine in a data frame
 data <- data.frame(type = c("original", rep("replication", times = nrow(rdat))),
                           site = c("original", rdat$Site),
                           r = c(sqrt(eta2o), rdat$rMain),
                           n = c(1/so^2 + 3, # only an implied sample size
                           rdat$NT), # actual sample size
                           fis = c(fiso, asinh(rdat$rMain)),
                           se_fis = c(so, 1/sqrt(rdat$NT - 3)))

 save(data, file="credentials_data.RData") 
 