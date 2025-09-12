

# INSTALLATION------------------------------------------------------------------
# Requires devtools package (run install.packages("devtools") if not installed)
devtools::install_github("mrc-ide/monty")
devtools::install_github("mrc-ide/dust2")
devtools::install_github("mrc-ide/odin2")
devtools::install_github("mrc-ide/dust")
devtools::install_github("mrc-ide/mcstate")

# Note the order - odin2 requires dust2 and monty to be installed
library(odin2)
library(dust2)
library(monty)
library(dust)
library(mcstate)

#packageVersion("odin2")
#packageVersion("dust2")
#packageVersion("monty")

#This code (SIR Odin) was based on an example from the Odin and Monty book ("4.7 Example: an age-structured SIR model with vaccination").

#Is this the correct way to implement a nested SIR model in Odin?
sir <- odin2::odin({
  update(S[]) <- S[i] - n_SI[i]
  update(I[]) <- I[i] + n_SI[i] - n_IR[i]
  update(R[]) <- R[i] + n_IR[i]
  update(incidence) <- incidence + sum(n_SI)
  
  # Individual probabilities of transition:
  p_SI[] <- 1 - exp(-lambda[i] * dt) # S to I
  p_IR <- 1 - exp(-gamma * dt) # I to R
  
  m <- parameter()
  
  # here s_ij[i, j] gives the mean number of contacts and individual in group
  # i will have with the currently infectious individuals of group j
  s_ij[, ] <- m[i, j] * I[j]
  
  # lambda[i] is the total force of infection on an individual in group i 
  lambda[] <- beta * sum(s_ij[i, ])
  
  # Draws from binomial distributions for numbers changing between
  # compartments:
  n_SI[] <- Binomial(S[i], p_SI[i])
  n_IR[] <- Binomial(I[i], p_IR)
  
  initial(S[]) <- S0[i]
  initial(I[]) <- I0[i]
  initial(R[]) <- 0
  initial(incidence, zero_every = 1) <- 0
  
  # User defined parameters - default in parentheses:
  S0 <- parameter()
  I0 <- parameter()
  beta <- parameter(0.6)
  gamma <- parameter(0.3)
  
  # Dimensions of arrays
  n_patch <- parameter(2)
  dim(S, S0, n_SI, p_SI) <- n_patch
  dim(I, I0, n_IR) <- n_patch
  dim(R) <- n_patch
  dim(m, s_ij) <- c(n_patch, n_patch)
  dim(lambda) <- n_patch

 # How should we include the likelihood here? How should we declare the "cases"?
   cases <- data()
   cases ~ Poisson(incidence)

})

###########################################################################################################

#What should the structure of the data look like?
# load example data
df_tidy <- read.csv(system.file("nested_sir_incidence.csv", package="mcstate"), 
                    stringsAsFactors = TRUE) # this is just example data from the package

library(dplyr)

data <- df_tidy %>%
  rename(time = day,      # já estava renomeando
         group = population) %>%  # renomeia population -> group
  mutate(group = as.character(group)) # garante que seja character

str(data)

data$group <- data$group

#How should we declare the set of parameters?
pars <- list(
#  n_patch = 2,
  group = c("A","B"),
  I0    = c(10, 10),
  S0    = c(1000, 1000),
  gamma = 0.3,
  beta  = 0.7,
  m  = matrix(c(5, 2, 2, 5), nrow = 2, ncol = 2)
)


#In this example:

#m[1, 1] = 5: Individuals in group 1 have 5 contacts with others in group 1.

#m[1, 2] = 2: Individuals in group 1 have 2 contacts with individuals in group 2.

#m[2, 1] = 2: Individuals in group 2 have 2 contacts with individuals in group 1.

#m[2, 2] = 5: Individuals in group 2 have 5 contacts with others in group 2.



 # pars <- list(
 #   list(
 #     I0 = 10,
 #     S0 = 1000,
 #     gamma = 0.1,
 #     beta  = 0.4,
 #     m     = matrix(4, nrow = 2, ncol = 2) # pode repetir em cada grupo
  #  ),
  #  list(
   #   I0 = 10,
  #    S0 = 1000,
  #    gamma = 0.1,
  #    beta  = 0.4,
  #    m     = matrix(4, nrow = 2, ncol = 2)
  #  )
 # )



library(dplyr)
library(tidyr)


print(data)

t <- c(0,data$time) 
print(t)

sys <- dust2::dust_system_create(generator = sir(), pars = pars,  time = 0, dt =1, 
                                 deterministic = FALSE, n_particles = 200, n_threads = 1, seed = 1, preserve_particle_dimension = TRUE)


t <- 0:max(data$time)
dust2::dust_system_set_state_initial(sys)
y <- dust2::dust_system_simulate(sys, t)

dim(y)

# Convert model output to labelled values
results <- dust2::dust_unpack_state(sys, y)

################################################################################################

#In dust_filter_create, should we include information about groups/populations?
# Create filter
filter <- dust2::dust_filter_create(sir, data = data, n_groups = 2, time_start = 0, n_threads = 1, n_particles = 200)
#filter <- dust2::dust_filter_create(sir, data = data, time_start = 0, n_particles = 200)



dust2::dust_likelihood_run(filter, pars)

library(ggplot2)

#Should we use monty::monty_packer or monty::monty_packer_grouped?
#packer <- monty::monty_packer(c("beta","gamma"),fixed=list(S0 = 1000,I0=10))  

packer <- monty::monty_packer_grouped(group = c("A", "B"),
                     scalar = c("beta", "gamma"),
                     shared = c("gamma","beta"),
                     fixed = list(S0 = c(1000, 1000), I0 = c(10, 10), m  = matrix(c(5, 2, 2, 5), nrow = 2, ncol = 2)))

likelihood <- dust2::dust_likelihood_monty(filter, packer, save_trajectories = TRUE) #save_trajectories = TRUE is optional and requires much memory


#We also have doubts about the priors, vcv, the sampler, monty_sample...
prior <- monty::monty_dsl({
 
  beta ~ Uniform(0.5, 0.9)
  gamma ~ Uniform(0.1, 0.9)
  
})


#######################################################################################################
posterior <- likelihood + prior


posterior

initial_chain_1 <- list(beta = 0.6, gamma = 0.2)
initial_chain_2 <- list(beta = 0.6, gamma = 0.2)

# Create the 'initial' list for the 'monty_sample'
pars_iniciais <- list(initial_chain_1, initial_chain_2)


vcv <- matrix(0.005, nrow = 2, ncol = 2)

diag(vcv) <- 0.01

print(vcv)

str(vcv)


sampler <- monty::monty_sampler_random_walk(vcv)

n_chains=2
n_iterations=1000

samples <- monty::monty_sample(posterior,sampler,n_iterations,n_chains=n_chains,initial = pars_iniciais)


######################### PLOT THE LAST TRAJECTORY OF EACH CHAIN############################

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

df_long$chain <- factor(df_long$chain, labels = paste("Chain", 1:length(unique(df_long$chain))))
ggplot() +
  geom_line(data = df_long, aes(x = time, y = cases, group = chain, color = chain), alpha = 0.8, lwd = 1.5, 
) +
  geom_point(data = data, aes(x = time, y = cases), color = "red", size = 2) +
  theme_minimal() +
  labs(title = "Last Trajectory: Simulations vs. Real Data",
       x = "Time (months)", 
       y = "Cases (Incidence)",
       color = "Chain") +
  theme(
    legend.position = c(0.95, 0.05),  # legenda no canto inferior direito (valores em proporção da área do gráfico)
    legend.justification = c("right", "bottom"),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  ) 



#################################Log posterior probability density################################
par(mfrow = c(1, 1)) # Reseta para 1 plot por vez

matplot(samples$density, type = "l", lty = 1,
        xlab = "Sample", ylab = "Log posterior probability density")

print(samples$density)


#################################################################

samples_df <- posterior::as_draws_df(samples)
samples_df
posterior::summarise_draws(samples_df)


unique(samples_df$.chain)
bayesplot::mcmc_pairs(samples_df)


###########################################################################################
library(ggplot2)


n_iter <- dim(samples$pars)[2]  # Number of iterations
n_chain <- dim(samples$pars)[3]  # Number of chains

# data.frame
df_samples <- data.frame(
  iter = rep(1:n_iter, times = n_chain),  
  temp_scale = c(samples$pars["temp_scale", , ]),  
  gamma = c(samples$pars["gamma", , ]), 
  I0 = c(samples$pars["I0", , ]),  
  eta = c(samples$pars["eta", , ]),  
  #rho = c(samples$pars["rho", , ]),  
  chain = rep(1:n_chain, each = n_iter)  # Chain index
)


########################################
# Trace plot for beta

ggplot(df_samples, aes(x = iter, y = beta, color = factor(chain))) +
  geom_line() +
  labs(title = "Evolution of the beta parameter by chain", color = "Chain") +
  labs(x = "Iteration", y = "Beta", color = "Chain") +
  theme_minimal() +
  scale_y_continuous(limits = c(0.0, 0.03))  


# Trace plot for gamma
ggplot(df_samples, aes(x = iter, y = gamma, color = factor(chain))) +
  geom_line() +
  labs(title = "Evolution of the gamma parameter by chain", color = "Chain") +
  labs(x = "Iteration", y = "Gamma", color = "Chain") +
  theme_minimal() +
  scale_y_continuous(limits = c(0.0, 10))  

#################################################################################
par(mfrow = c(1, 1)) # Reseta para 1 plot por vez
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
  xlab = "Time", 
  ylab = "Cases", 
  cex.axis = 1.5,         
  cex.lab = 1.5,          
)


colors <- rainbow(n_chains)
# Now, iterate over the chains and plot the last trajectory of each one
for (chain in 1:n_chains) {
  print(chain)
  pars2 <- list(
    N = 1000, 
    beta = as.numeric(samples$pars["beta", n_iterations, chain]), 
    gamma = as.numeric(samples$pars["gamma", n_iterations, chain]),
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



