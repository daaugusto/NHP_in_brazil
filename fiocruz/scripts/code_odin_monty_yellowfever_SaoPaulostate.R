# This is a demonstration of key features of the updated modelling packages from Imperial College London:
# monty - A new package for parameter estimation, an update of mcstate - https://mrc-ide.github.io/monty/
# odin2 - Update of odin; domain specific language (DSL) for differential equation implementation - https://mrc-ide.github.io/odin2/
# dust2 - Update of dust; engine for running dynamic systems - https://mrc-ide.github.io/dust2/
# Note that eventually odin2 and dust2 will be folded into the existing odin and dust packages
# 
# For tutorials, see https://mrc-ide.github.io/odin-monty/odin.html

# INSTALLATION------------------------------------------------------------------
# Requires devtools package (run install.packages("devtools") if not installed)

#update.packages(ask = FALSE)
#devtools::install_github("mrc-ide/monty")
#devtools::install_github("mrc-ide/dust2")
#devtools::install_github("mrc-ide/odin2")

# Note the order - odin2 requires dust2 and monty to be installed

############################################################################
# clean global environment
#rm(list = ls())

# clean specific variables, for example:
#rm(beta, gamma, temp_scale, posterior_samples, posterior_summary_last)
##############################################################################
#set.seed(531)
# MODEL CODE--------------------------------------------------------------------
# This code is a simple SIR model with incidence as an output
# Differences from implementation of the same model in odin are noted in comments
# Note that not much has changed in the implementation of this model!
# This code is reproduced in the file sir_odin2.R in this folder
# To load from the separate file, run sir <- odin2::odin("new_packages_demo/sir_odin2.R")

library(odin2)
library(dust2)
library(monty)

# Load common functions
source("../R/common.R")

sir <- odin2::odin({
  N <- parameter(1000) # In odin2, parameter() is used for inputs where user() was used in odin
  I0 <- parameter(10)
  beta <- parameter(0.4)
  gamma <- parameter(0.1)
  
  p_SI <- 1 - exp(-beta * I / N * dt) # Note parameter dt - this is a time interval set when the model is set up via dust_system_create
  p_IR <- 1 - exp(-gamma * dt)
  n_SI <- Binomial(S, p_SI) # In odin2, Binomial() is used instead of rbinom()
  # See https://mrc-ide.github.io/odin2/articles/functions.html#distribution-functions for available distributions
  n_IR <- Binomial(I, p_IR)
  
  update(S) <- S - n_SI
  update(I) <- I + n_SI - n_IR
  update(R) <- R + n_IR
  #update(incidence) <- n_SI
  update(incidence) <- max(n_SI, 0.1)  
  
  initial(S) <- N - I0
  initial(I) <- I0
  initial(R) <- 0
  initial(incidence, zero_every = 1) <- 0 #ensures that the incidence will always start from zero, without "accumulating" previous cases


  #cases <- data()
  #cases ~ Poisson(incidence)
  
 
  # Negative Binomial
  rho <- parameter(0.1)
  eta <- parameter(10)
  cases <- data()
  cases ~ NegativeBinomial(size = eta, mu = rho * incidence)
  
})

###########################################################################################################


# load example data
data<- read.csv("../data/cases_per_state_pnh_sp_mg_pr.csv", header= TRUE, sep = ",") # this is just example data from the package
#data <- read.csv(system.file("sir_incidence.csv", package="mcstate")) # this is just example data from the package 
head(data)

data_week_selected <- data[, c(1,2)] #to get data only from SP for now
head(data_week_selected)


plot(data_week_selected[,1], data_week_selected[,2], 
     xlab = "day",  
     ylab = "Number of Cases SP",  
     pch = 16, ylim = c(0,30), cex =0.7)                      

data_week_selected$week <- ceiling(data_week_selected$day/30) #grouping monthly
head(data_week_selected)

data<- aggregate(cases_A ~ week, data = data_week_selected, sum)
head(data)

library(dplyr)

data <- data %>%
  rename(cases = cases_A)  #just in case if your columns have different names

#data <- data %>%
 # rename(time = day)

data <- data %>%
  rename(time = week)

##################################################################################################


######################## FUNCTION  FIT DATA  ###########################

results <- adjust_distributions(data$cases, titulo = "Casos SP - Likelihood Fit")



pars <- list(N = 2000, I0 =1, beta=0.4, gamma = 0.1, rho = 0.1, eta = 10)

# Setting up the model
# Note the parameter dt - this defines the time interval of each calculation step
# dt can be in any units but must have value 1 or less and be the inverse of an integer; here dt = 1 day
# For more flexibility (e.g. calculation steps 5 days apart), dt can have alternative units with a separate time interval defined
# Model can be set to deterministic mode by setting deterministic = TRUE
# Parameter seed can be used to produce reproducible stochastic output; set to NULL to randomize
# Parameter preserve_particle_dimension is set to TRUE so that output has same dimensionality even if n_particles = 1

sys <- dust2::dust_system_create(sir, pars, time = 0, dt = 1, deterministic = FALSE, n_particles = 1, n_threads = 1, seed = 1,
                                 preserve_particle_dimension = TRUE)


# The output time sequence has the same units as the calculation time interval dt
# Note that the separation of the output time points needs to be compatible with dt
# For example, if dt = 1, the output time sequence cannot have fractional points
#t <- seq(1, 100, by = 1)  #mudar de acordo com nossos dados (numero de meses ou semanas)

t <- c(0,data$time) #mudar de acordo com nossos dados (numero de meses ou semanas)
print(t)
#length(data)
# Initialize model to starting conditions
dust2::dust_system_set_state_initial(sys)

# Run model over time sequence
# Note that if you want to re-run with a different time sequence, you must go back and re-run from line 63
# Output is an array with dimensions (number of output values, number of particles, number of time points)
y <- dust2::dust_system_simulate(sys, t)
dim(y)

# Convert model output to labelled values
results <- dust2::dust_unpack_state(sys, y)

#Plot SIR on graph
matplot(x = t, y = t(results$S), type="p", pch = 16, col=2, xlab = "Day", ylab = "", ylim=c(0, pars$N))
matplot(x = t, y = t(results$I), type="p", pch = 16, col=3, add=TRUE)
matplot(x = t, y = t(results$R), type="p", pch = 16, col=4, add=TRUE)
legend("topright",c("S","I","R"),pch=c(16,16,16),col=c(2,3,4))

#Plot incidence on graph
matplot(x = t, y = t(results$incidence), type="p", pch = 16, col=1, xlab = "Day", ylab= "Incidence")

plot(t, dust2::dust_unpack_state(sys, y)$incidence,
     type = "l",
     xlab = "Time",
     ylab = "Incidence of cases",
     col = "blue",
     ylim = c(0, max(c(dust2::dust_unpack_state(sys, y)$incidence, data$cases)) * 1.2))
points(data$time, data$cases,
       pch = 19,
       col = "red")

# Obtaining individual outputs
# Index of all outputs can be obtained; this may be useful if number of outputs is large and only certain values are needed
index = dust2::dust_unpack_index(sys)
index
S = y[index$S,,]
I = y[index$I,,]
R = y[index$R,,]
incidence = y[index$incidence,,]

plot(t, dust2::dust_unpack_state(sys, y)$incidence, type = "l", col = "#000088ff")


################################################################################################

# Create filter
filter <- dust2::dust_filter_create(sir, data = data, time_start = 0, n_particles = 20)
dust2::dust_likelihood_run(filter, pars)

# Create packer - divide input parameters into estimated (beta + gamma) and fixed (N + I0)
# Significantly simplified from mcstate!
#packer <- monty::monty_packer(c("beta","gamma"),fixed=list(N = 1000, I0=10))  

packer <- monty::monty_packer(c("beta","gamma","rho","eta"),fixed=list(N = 1000, I0=10))  

#likelihood <- dust2::dust_likelihood_monty(filter, packer)

likelihood <- dust2::dust_likelihood_monty(filter, packer, save_trajectories = TRUE) #save_trajectories = TRUE is optional and requires much memory

# Set prior likelihood distributions for estimated parameters
# Here a simple uniform distribution is used for beta and gamma, with permitted minimum/maximum values
# See https://mrc-ide.github.io/odin2/articles/functions.html#distribution-functions for available distributions
prior <- monty::monty_dsl({
 
  beta ~ Uniform(min=0.2,max=0.6)
  gamma ~ Uniform(min=0.01,max=0.20)
  
  rho ~ Uniform(0, 1)
  eta ~ Uniform(10, 1000)
  #eta ~ Exponential(mean = 100)
})

prior$domain #Check limits (can adjust manually by adjusting values in prior$domain)

###################################################################################################


n_streams <- 100
r <- monty::monty_rng_create(n_streams = n_streams)
prior_samples <- matrix(monty::monty_model_direct_sample(prior, r), nrow = n_streams)
colnames(prior_samples) <- prior$parameters

bayesplot::mcmc_pairs(prior_samples)

####################################################################################################


#######################################################################################################
posterior <- likelihood + prior

pars <- list(gamma=0.1, beta = 0.45,  rho = 0.1, eta = 10)
no_param <- length(pars)


vcv <- matrix(0.005, 4, 4)  # Fill everything with the off-diagonal value

diag(vcv) <- 0.01           # Changes the elements of the main diagonal

# Random walk sampler (other samplers are available)
sampler <- monty::monty_sampler_random_walk(vcv)

# Can't use HMC with stochastic models because require a deterministic density
#sampler <- monty::monty_sampler_hmc(epsilon = 0.1, n_integration_steps = 10)

n_chains=5
n_iterations=500

#samples <- monty::monty_sample(posterior,sampler,n_iterations,initial=array(rep(c(0.05,0.25),n_chains),
#                                                                    dim=c(2,n_chains)),n_chains=n_chains)

samples <- monty::monty_sample(posterior,sampler,n_iterations,n_chains=n_chains)
#if we leave it like this, without the initial conditions, it starts with random values within the range of the prior

samples

dim(samples$pars)
dim(samples$observations$trajectories)
object.size(samples$observations$trajectories)

library(ggplot2)
library(dplyr)
library(tidyr)

str(samples)  # Estrutura do objeto
names(samples)  # Lista dos elementos dentro de samples


########################################################################################################


######################### PLOT THE LAST TRAJECTORY OF EACH CHAI N############################

library(dplyr)
library(tidyr)
library(ggplot2)

# Total number of chains and iterations
n_chains <- dim(samples$observations$trajectories)[4]
n_iterations <- dim(samples$observations$trajectories)[3]
n_time <- dim(samples$observations$trajectories)[2]

# Getting the trajectories from the last iteration for all chains

trajectories <- samples$observations$trajectories[4, , n_iterations, ]  # (81 x n_chains)

# dataframe
df_plot <- as.data.frame(trajectories)
df_plot$time <- 1:nrow(df_plot)  # Adding the time column correctly

# Transforming to long format
df_long <- df_plot %>%
  pivot_longer(cols = -time, names_to = "chain", values_to = "cases") %>%
  mutate(chain = as.factor(chain))  # Convert to factor to differentiate colors


ggplot() +
  geom_line(data = df_long, aes(x = time, y = cases, group = chain, color = chain), alpha = 0.8) +
  geom_point(data = data, aes(x = time, y = cases), color = "red", size = 2) +
  theme_minimal() +
  labs(title = "Last Trajectory: Simulations vs. Real Data",
       x = "Time (months)", 
       y = "Cases (Incidence)",
       color = "Chain") +
  theme(legend.position = "right")


#################################Log posterior probability density################################


matplot(samples$density, type = "l", lty = 1,
        xlab = "Sample", ylab = "Log posterior probability density")

matplot(
  samples$density, 
  type = "l", 
  lty = 1, 
  lwd = 2,             
  xlab = "Sample", 
  ylab = "Log posterior probability density", 
  cex.lab = 1.5,       
  cex.axis = 1.5,      
  col = rainbow(ncol(samples$density)) # Different colors for each chain
)

box(lwd = 1)  
legend("bottomright", 
       legend = paste("Chain", 1:ncol(samples$density)), 
       col = rainbow(ncol(samples$density)), 
       lty = 1, 
       lwd = 3, 
       cex = 1.2) 



##############################################################################################################################


samples_df <- posterior::as_draws_df(samples)
samples_df
posterior::summarise_draws(samples_df)



###########################################################################################
library(ggplot2)


n_iter <- dim(samples$pars)[2]  # Number of iterations
n_chain <- dim(samples$pars)[3]  # Number of chains

# data.frame
df_samples <- data.frame(
  iter = rep(1:n_iter, times = n_chain),  
  beta = c(samples$pars["beta", , ]),  
  gamma = c(samples$pars["gamma", , ]),  
  eta = c(samples$pars["eta", , ]),  
  rho = c(samples$pars["rho", , ]),  
  chain = rep(1:n_chain, each = n_iter)  # Chain index
)

# Trace plot for beta

ggplot(df_samples, aes(x = iter, y = beta, color = factor(chain))) +
  geom_line() +
  labs(title = "Evolution of the beta parameter by chain", color = "Chain") +
  labs(x = "Iteration", y = "Beta", color = "Chain") +
  theme_minimal()


# Trace plot for gamma
ggplot(df_samples, aes(x = iter, y = gamma, color = factor(chain))) +
  geom_line() +
  labs(title = "Evolution of the gamma parameter by chain", color = "Chain") +
  labs(x = "Iteration", y = "Gamma", color = "Chain") +
  theme_minimal()


# Trace plot para eta
ggplot(df_samples, aes(x = iter, y = eta, color = factor(chain))) +
  geom_line() +
  labs(title = "Evolution of the eta parameter by chain", color = "Chain") +
  labs(x = "Iteration", y = "Eta", color = "Chain") +
  theme_minimal()

# Trace plot para rho
ggplot(df_samples, aes(x = iter, y = rho, color = factor(chain))) +
  geom_line() +
  labs(title = "Evolution of the rho parameter by chain", color = "Chain") +
  labs(x = "Iteration", y = "rho", color = "Chain") +
  theme_minimal()

pars2 <- list(N = 1000, I0 = 10,  beta= as.numeric(samples$pars[1,n_iterations,1]), 
              gamma = as.numeric(samples$pars[2,n_iterations,1]))

likelihoods <- array(NA, dim = c(n_chains, nrow(data), 10))  

# Initial line with main chart
matplot(
  x = data$time, 
  y = data$cases, 
  #type = "l", 
  type = "p",                  
  col = 1, 
  pch = 16,  
  lwd = 1,  
  cex=1.5,
  xlab = "Month", 
  ylab = "Cases SP", 
  cex.axis = 1.5,         
  cex.lab = 1.5,          
  ylim = c(0, 150) 
)

colors <- rainbow(n_chains)
# Now, iterate over the chains and plot the last trajectory of each one
for (chain in 1:n_chains) {
  print(chain)
  pars2 <- list(
    N = 1000, 
    I0 = 10, 
    beta = as.numeric(samples$pars[1, n_iterations, chain]), 
    gamma = as.numeric(samples$pars[2, n_iterations, chain]),
    rho = as.numeric(samples$pars[3, n_iterations, chain]), 
    eta = as.numeric(samples$pars[4, n_iterations, chain])
  )
  
  sys2 <- dust2::dust_system_create(sir, pars2, time = 0, dt = 1, deterministic = FALSE, 
                                    n_particles = 10, n_threads = 1, seed = 3, 
                                    preserve_particle_dimension = TRUE)
  
  index2 <- dust2::dust_unpack_index(sys2)
  dust2::dust_system_set_state_initial(sys2)
  y2 <- dust2::dust_system_simulate(sys2, t = data$time)
  
  for (i in 1:nrow(data)) {
    likelihoods[chain, i, ] <- dust2::dust_system_compare_data(sys2, data[i, ])
  }
    
  # Access the last trajectory (last iteration) for the "incidence" variable  
  #last_trajectory <- y2[index2$incidence, n_iterations, chain]
  last_trajectory <- y2[index2$incidence, 10, ]
  
  # Plot the last trajectory of the chain  
  lines(data$time, last_trajectory, col = colors[chain], lwd = 2)
}

legend("topright", 
       legend = c("Real Data", paste("Chain", 1:n_chains)),   
       col = c(1, colors),                                      
       pch = c(16, rep(NA, n_chains)),                          
       lwd = c(NA, rep(2, n_chains)),                          
      cex = 1.2) 


########################################################




