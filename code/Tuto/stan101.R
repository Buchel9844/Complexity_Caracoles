library(rstan)
rstan_options(auto_write = TRUE)            # Cache compiled Stan programs
options(mc.cores = parallel::detectCores()) # Parallelize chains
util <- new.env()
source('Code/stan_utility.R', local=util)

#input_data <- read_rdump("data/intro.data.R") # Formats the observed data into a list which will then get consumed by RStan.
input_data <- list( N = 5, 
                    y = c(1.22819583757764, 1.03747386101563,
                          1.29643993528305, 1.26174973066331, 0.587372610448683))

# Call the main RStan command which conditions our joint model on the observed data and 
# runs Hamiltonian Monte Carlo to estimate expect values. 
# To ensure reproducible behavior I make sure to set an explicit seed for Stan’s internal
# pseudo random number generator.

fit <- stan(file='Code/stan101.stan', data=input_data, seed=4938483)

util$check_all_diagnostics(fit) #check the available diagnostics for any indication of bias
# no problem 
# Use the samples output by Hamiltonian Monte Carlo to construct Markov chain Monte Carlo estimators. 
print(fit)
# Note: lp__ is the log target (unnormalized) density function evaluated at each sample
params <- extract(fit) # more bespoke estimators -> access the samples
names(params)
dim(params)
length(params$theta)

# visualize the marginal posterior distributions for each component of the model configuration space,
# in this case just the one theta
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

hist(params$theta, breaks=seq(-4, 4, 0.25),
     col=c_dark, main="", border=c_dark_highlight,
     xlim=c(-4, 4), xlab="theta", yaxt='n', ylab="")

hist(params$y[,1], breaks=seq(-7, 7, 0.5), prob=T,
     col=c_dark, border=c_dark_highlight, main="",
     xlim=c(-7, 7), xlab="y[1]", yaxt='n', ylab="")
