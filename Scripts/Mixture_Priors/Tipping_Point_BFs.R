#   ____________________________________________________________________________
#   Libraries                                                               ####

library(dplyr)
library(knitr)
library(xtable)
library(repmix)

source("Scripts/Mixture_Priors/RepMixFun_BF.R")

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

# Mean and Variance Unit Informative Prior
mu_UIP <- 0
tau_UIP <- 2


# Sequence of omega values
omega_seq <- seq(0, 1, by = 0.1)

#   ____________________________________________________________________________
#   Tipping Point BF for theta                                              ####


bf_df_list <- list()

for (lab in seq_along(tr)) {
  
  # Compute BF_{0:1}(ω) over the grid
  bf_theta_vals <- sapply(omega_seq, function(w) {
    bf_theta_mix(
      tr        = tr[lab],
      sr        = sr[lab],
      to = to,
      so    = so,
      m      = mu_UIP,
      v = tau_UIP,
      w = w)
  })
  
  bf_omega_vals <- sapply(omega_seq, function(w) {
    bf_omega_mix(
      tr        = tr[lab],
      sr        = sr[lab],
      to = to,
      so    = so,
      m      = mu_UIP,
      v = tau_UIP,
      w_null = 1,
      w_alt = w)
  })
  
  # Store in a data frame
  bf_df <- data.frame(
    omega = omega_seq,
    BF_theta    = bf_theta_vals,
    BF_omega    = bf_omega_vals,
    tr        = tr[lab],
    sr        = sr[lab],
    rnumber   = as.factor(lab)
  )
  bf_df_list[[lab]] <- bf_df
}

# Combine into one data frame
bf_all <- do.call(rbind, bf_df_list)

# Create a new column for ordering
bf_all$rep_order <- ifelse(bf_all$rnumber == 4, "p", as.character(bf_all$rnumber))

# Create the original rep_setting column with proper labels
bf_all$rep_setting <- paste0(
  "{hat(theta)[r*",
  ifelse(bf_all$rnumber == 4, "p", bf_all$rnumber),
  "] == ",
  round(bf_all$tr, 2),
  "}*',' ~ sigma[r*",
  ifelse(bf_all$rnumber == 4, "p", bf_all$rnumber),
  "] == ",
  round(bf_all$sr, 2)
)

bf_all$rep_order <- factor(bf_all$rep_order, levels = c("1", "2", "3", "p"))


##  ............................................................................
##  Plot                                                                    ####

tipp_BF_theta <- ggplot(bf_all, aes(x = omega, y = BF_theta,  color = factor(rnumber))) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_y_continuous(trans = "log10") +            # log‐scale for wide range
  scale_x_continuous(
    breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
    labels = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) 
  ) +
  labs(
    x =expression("Prior Weight" ~ omega [H][1]),
    y = expression(
      BF[0][1] *                 
        "(" *                    
        hat(theta)[r] *          
        "|" *                    
        H[1] *                   
        ": " *                   
        omega == omega[H[1]] *   
        ")"                      
    )
  ) +
  scale_color_manual(
    values = c("1" = "#E69F00", "2" = "#009E20", "3" = "#0072B2", "4" = "#AA4499"),
    labels = c(
      expression(" " ~ hat(theta)[r * 1] == 0.29 ~ ", " ~ sigma[r * 1] == 0.11),
      expression(" " ~ hat(theta)[r * 2] == 0.25 ~ ", " ~ sigma[r * 2] == 0.09),
      expression(" " ~ hat(theta)[r * 3] == -0.18 ~ ", " ~ sigma[r * 3] == 0.11),
      expression(" " ~ hat(theta)[r * p] == 0.14 ~ ", " ~ sigma[r * p] == 0.06)
    )
  ) + 
  theme_bw() +
  theme(
    strip.placement = "outside", 
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



print(tipp_BF_theta)


ggsave(
  filename = "tipp_BF_theta.pdf", path = "Plots/Mixture_Prior", plot = tipp_BF_theta,
  width = 17, height = 7.5, device = "pdf", dpi = 500, useDingbats = FALSE
)


tipp_BF_omega <- ggplot(bf_all, aes(x = omega, y = BF_omega,  color = factor(rnumber))) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_y_continuous(trans = "log10") +  # log‐scale for wide range
  scale_x_continuous(
    breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
    labels = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) 
  ) +
  labs(
    x =expression("Prior Weight" ~ omega [H][d]),
    y = expression(
      BF[d][c] *                 
        "(" *                    
        hat(theta)[r] *          
        "|" *                    
        H[d] *                   
        ": " *                   
        omega == omega[H[d]] *   
        ")"                      
    )
  ) +
  scale_color_manual(
    values = c("1" = "#E69F00", "2" = "#009E20", "3" = "#0072B2", "4" = "#AA4499"),
    labels = c(
      expression(" " ~ hat(theta)[r * 1] == 0.29 ~ ", " ~ sigma[r * 1] == 0.11),
      expression(" " ~ hat(theta)[r * 2] == 0.25 ~ ", " ~ sigma[r * 2] == 0.09),
      expression(" " ~ hat(theta)[r * 3] == -0.18 ~ ", " ~ sigma[r * 3] == 0.11),
      expression(" " ~ hat(theta)[r * p] == 0.14 ~ ", " ~ sigma[r * p] == 0.06)
    )
  ) + 
  theme_bw() +
  theme(
    strip.placement = "outside", 
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


print(tipp_BF_omega)


ggsave(
  filename = "tipp_BF_omega.pdf", path = "Plots/Mixture_Prior", plot = tipp_BF_omega,
  width = 17, height = 7.5, device = "pdf", dpi = 500, useDingbats = FALSE
)
