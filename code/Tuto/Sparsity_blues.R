# Sparsity Blues
# Michael Betancourt
# May 2021
# https://betanalpha.github.io/assets/case_studies/modeling_sparsity.html
#########################################################################
# 2. Sparse Population Models
#########################################################################

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

util <- new.env()
source('Code/stan_utility.R', local=util)

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

c_light_trans <- c("#DCBCBC80")
c_dark_trans <- c("#8F272780")
c_green_trans <- c("#00FF0080")

par(family="CMU Serif", las=1, bty="l", cex.axis=1, cex.lab=1, cex.main=1,
    xaxs="i", yaxs="i", mar = c(5, 5, 3, 1))

set.seed(58533858)

#  simple simulation - six cluster below an inner scale of 0.1 and three venture out past 10.
K <- 9

theta_true <- c(-0.09604655,  0.02063505,  0.01246716,
                -0.04487767, -0.13519031,  0.09421304,
                -29.449233,  32.997172,  18.517443)

#  inform the location of a normal observational model
N <- K
sigma <- 0.5
context_idx <- rep(1:K, 1)

y <- rnorm(K, theta_true[context_idx], sigma) # mean is 0 and sd is 0.5

data <- list("N" = N, "K" = K, "context_idx" = context_idx,
             "y" = y, "sigma" = sigma)
# normal distribution can not accomodate both the inner and outer space of parameters

writeLines(readLines("Code/normal_narrow.stan"))

fit <- stan(file='Code/normal_narrow.stan', data=data,
            seed=4938483, refresh=1000)

util$check_all_diagnostics(fit)
