
#   ____________________________________________________________________________
#   Libraries                                                               ####
library(ggplot2)
library(dplyr)
library(tidyr)
library(repmix)

source("Scripts/Mixture_Priors/RepMixFun_BF.R")

#   ____________________________________________________________________________
#   Parameters Setting                                                      ####

to       <- 0.2
m        <- 0.0
v        <- 2.0
tr_equal <- to
tr_diff  <- to + 0.15

sigma_vals <- 10^seq(-1.47,0, length.out = 200)


#   ____________________________________________________________________________
#   Bayes Factors for omega                                                 ####


bf_equal_equal  <- sapply(sigma_vals, function(s) {
  bf_omega_mix(tr     = tr_equal,
               sr     = s,
               to     = to,
               so     = s,
               m      = m,
               v      = v,
               w_null = 0, w_alt = 1)
})
bf_equal_diff   <- sapply(sigma_vals, function(s) {
  bf_omega_mix(tr     = tr_diff,
               sr     = s,
               to     = to,
               so     = s,
               m      = m,
               v      = v,
               w_null = 0, w_alt = 1)
})



df <- tibble(
  sigma     = sigma_vals,
  `tr = to`       = bf_equal_equal,
  `tr != to` = bf_equal_diff
) %>%
  pivot_longer(
    cols      = -sigma,
    names_to  = "scenario",
    values_to = "BF"
  )




bf_asymp <- ggplot(df, aes(x = sigma, y = BF,
                           colour = scenario, linetype = scenario)) +
  geom_line(size = 1) +
  geom_line(
    size        = 1,
    arrow       = arrow(length = unit(0.35, "cm"),
                        ends   = "first",
                        type   = "closed"),
    show.legend = FALSE
  ) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::label_number(accuracy = 0.1)  ) +
  scale_colour_manual(
    values = c("grey40","orange2"),
    labels = c(
      expression(hat(theta)[o] != hat(theta)[r]),
      expression(hat(theta)[o] == hat(theta)[r])
    )
  ) +
  scale_linetype_manual(
    values = c("dashed","solid"),
    labels = c(
      expression(hat(theta)[o] != hat(theta)[r]),
      expression(hat(theta)[o] == hat(theta)[r])
    )
  ) +
  labs(
    x      = expression(paste(sigma[r], ", ", sigma[o])),
    y      = expression(BF[dc] ~~ (hat(theta)[r])),
    colour = NULL, linetype = NULL
  )  +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position  = "top",
    strip.placement  = "outside",
    strip.background = element_blank(),
    strip.text.x     = element_text(size = 22),
    axis.text        = element_text(size = 18),
    axis.title       = element_text(size = 22),
    legend.text      = element_text(size = 18)
  )

ggsave(
  filename   = "bf_asymp.pdf",
  path       = "Plots/Mixture_Prior",
  plot       = bf_asymp,
  width      = 17, height = 7.5,
  device     = "pdf", dpi = 500,
  useDingbats = FALSE
)



