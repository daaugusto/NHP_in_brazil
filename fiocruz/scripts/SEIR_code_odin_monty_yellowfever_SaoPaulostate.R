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

# INSTALLATION------------------------------------------------------------------
# Requires devtools package (run install.packages("devtools") if not installed)
devtools::install_github("mrc-ide/monty")
devtools::install_github("mrc-ide/dust2")
devtools::install_github("mrc-ide/odin2")
# Note the order - odin2 requires dust2 and monty to be installed


#install.packages("gamlss.dist", dependencies = TRUE)
#install.packages("arm", dependencies = TRUE)

#library("rstan")
#library("gamlss.dist")
#library("arm")
#rstan_options(auto_write = TRUE)
#options(warn = -1)
# Modelo SEIR Estocástico em Odin2

seir <- odin2::odin({
  # Parâmetros do Modelo
  N <- parameter(1000) # População total
  I0 <- parameter(10)  # Indivíduos infecciosos iniciais
  E0 <- parameter(10)   # Indivíduos expostos iniciais (geralmente 0 ou um pequeno número)
  
  beta <- parameter(0.4) # Taxa de infecção (S -> E)
  gamma <- parameter(0.1) # Taxa de recuperação (I -> R)
  sigma <- parameter(0.2) # Taxa de progressão de exposição para infecção (E -> I)

  p_SE <- max(0, min(1, 1 - exp(-beta * I / N * dt)))
  p_EI <- max(0, min(1, 1 - exp(-sigma * dt)))
  p_IR <- max(0, min(1, 1 - exp(-gamma * dt)))
  
  
  # Número de Indivíduos que se Movem Entre Classes (estocástico via Binomial)
  n_SE <- Binomial(S, p_SE) # Novas infecções (S -> E)
  n_EI <- Binomial(E, p_EI) # Novas infecciosidades (E -> I)
  n_IR <- Binomial(I, p_IR) # Novas recuperações (I -> R)
  
  # Equações de Atualização para o Próximo Passo de Tempo
  update(S) <- S - n_SE
  update(E) <- E + n_SE - n_EI
  update(I) <- I + n_EI - n_IR
  update(R) <- R + n_IR
  
  # Incidência: Novas Infecções (n_SE) ou Novas Infecciosidades (n_EI)
  # A "incidência de casos" pode ser interpretada como o número de pessoas
  # que se tornam infecciosas (n_EI). Se seus dados são de "novos infectados",
  # n_SE seria mais apropriado. Mas geralmente, em epidemiologia, "casos"
  # se referem a indivíduos que se tornaram sintomáticos ou infecciosos.
  # Por isso, n_EI é uma escolha comum para modelar "novos casos diários".
  update(incidence) <- max(n_EI, 0.1) # Garante que a lambda da Poisson seja > 0
  # É importante que seja o que você mede!
  
  # Estados Iniciais
  initial(S) <- N - I0 - E0 # S inicial ajustado para incluir E0
  initial(E) <- E0
  initial(I) <- I0
  initial(R) <- 0
  initial(incidence, zero_every = 1) <- 0 # Zera a incidência a cada passo para contagem incremental
  
  # Likelihood de Poisson para os casos observados
  #cases <- data()
  #cases ~ Poisson(incidence)
  
  # Negative Binomial
  #rho <- parameter(0.5)
  eta <- parameter(10)
  cases <- data()
  #cases ~ NegativeBinomial(size = eta, mu = rho*incidence)
  cases ~ NegativeBinomial(size = eta, mu = incidence) #odin and monty book's suggestion
  
})
###########################################################################################################

# load example data
data<- read.csv("fiocruz/data/cases_per_state_pnh_sp_mg_pr.csv", header= TRUE, sep = ",") # this is just example data from the package
#data<- read.csv("C:/Users/angel/OneDrive/Desktop/Angelica_R/sir_temp_example_epizotias/cases_per_state_pnh_sp_mg_pr.csv", header= TRUE, sep = ",") # this is just example data from the package
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
#  rename(time = day)

data <- data %>%
  rename(time = week)

#se quiser cortar os dados
#data <- data[data$time >= 20 & data$time <= 45, ]

#data$time <- data$time - min(data$time) +1

##################################################################################################


######################## FUNCTION  FIT DATA  ###########################

#source("fit_data_likelihood_poisson_nb_normal_binomial.R")  # Function script path

results <- adjust_distributions(data$cases, titulo = "Casos SP - Likelihood Fit")

pars <- list(N = 1000, I0 = 10, E0=10, beta=0.4, gamma = 0.1, sigma =0.2, eta=10)
#pars <- list(N = 1000, I0 =10, beta=0.4, gamma = 0.1, rho = 0.1, eta = 10)

# Setting up the model
# Note the parameter dt - this defines the time interval of each calculation step
# dt can be in any units but must have value 1 or less and be the inverse of an integer; here dt = 1 day
# For more flexibility (e.g. calculation steps 5 days apart), dt can have alternative units with a separate time interval defined
# Model can be set to deterministic mode by setting deterministic = TRUE
# Parameter seed can be used to produce reproducible stochastic output; set to NULL to randomize
# Parameter preserve_particle_dimension is set to TRUE so that output has same dimensionality even if n_particles = 1

#sys <- dust2::dust_system_create(sir, pars, time = 0, dt = 1, deterministic = FALSE, n_particles = 1, n_threads = 1, seed = 1,
#                                preserve_particle_dimension = TRUE)

sys <- dust2::dust_system_create(generator = seir, pars = pars, time = 0, dt =1, 
                                 deterministic = FALSE, n_particles = 200, n_threads = 1, seed = 1, preserve_particle_dimension = TRUE)


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
#matplot(x = t, y = t(results$S), type="p", pch = 16, col=2, xlab = "Day", ylab = "", ylim=c(0, pars$N))
#matplot(x = t, y = t(results$I), type="p", pch = 16, col=3, add=TRUE)
#matplot(x = t, y = t(results$R), type="p", pch = 16, col=4, add=TRUE)
#legend("topright",c("S","I","R"),pch=c(16,16,16),col=c(2,3,4))

#Plot incidence on graph
#matplot(x = t, y = t(results$incidence), type="p", pch = 16, col=1, xlab = "Day", ylab= "Incidence")


# Obtaining individual outputs
# Index of all outputs can be obtained; this may be useful if number of outputs is large and only certain values are needed
index = dust2::dust_unpack_index(sys)
index
S = y[index$S,,]
I = y[index$I,,]
R = y[index$R,,]
incidence = y[index$incidence,,]

#plot(t, dust2::dust_unpack_state(sys, y)$incidence, type = "l", col = "#000088ff")


################################################################################################

# Create filter
filter <- dust2::dust_filter_create(seir, data = data, time_start = 0, n_particles = 200)
dust2::dust_likelihood_run(filter, pars)


# Create packer - divide input parameters into estimated (beta + gamma) and fixed (N + I0)
# Significantly simplified from mcstate!
#packer <- monty::monty_packer(c("beta","gamma"),fixed=list(N = 1000, I0=10))  

packer <- monty::monty_packer(c("beta","gamma","sigma","eta", "E0", "I0"),fixed=list(N = 1000))  


likelihood <- dust2::dust_likelihood_monty(filter, packer, save_trajectories = TRUE) #save_trajectories = TRUE is optional and requires much memory

# Set prior likelihood distributions for estimated parameters
# Here a simple uniform distribution is used for beta and gamma, with permitted minimum/maximum values
# See https://mrc-ide.github.io/odin2/articles/functions.html#distribution-functions for available distributions
prior <- monty::monty_dsl({
  

  #beta ~ TruncatedNormal(0.5, 0.5, min = 0, max = 1)  
  #gamma ~ TruncatedNormal(0.5, 0.5, min = 0, max = 1)  
  #sigma ~ TruncatedNormal(0.5, 0.5, min = 0, max = 1)  
  #I0 ~ Uniform(min=5,max=10)
  
  beta ~ Uniform(0.1, 10.0)
  gamma ~ Uniform(0.01, 10.0)
  sigma ~ Uniform(0.1, 10.0)
  #gamma ~ Gamma(shape=1,scale=0.5) 
  #rho ~ Uniform(0.1, 1.0)
  eta ~ Exponential(mean = 1000)
  
  I0 ~ Uniform(2, 10)
  E0 ~ Uniform(2, 10)  
  
})


#######################################################################################################
posterior <- likelihood + prior

pars <- list(gamma=0.1, beta = 0.45, sigma=0.2, eta=10, E0=5, I0=5)
no_param <- length(pars)


vcv <- matrix(0.005, nrow = 6, ncol = 6)

diag(vcv) <- 0.01

#print(vcv)

# Random walk sampler (other samplers are available)

sampler <- monty::monty_sampler_random_walk(vcv)

# Can't use HMC with stochastic models because require a deterministic density
#sampler <- monty::monty_sampler_hmc(epsilon = 0.1, n_integration_steps = 10)
#sampler <- monty::monty_sampler_adaptive(vcv)

n_chains=10
n_iterations=10000

#samples <- monty::monty_sample(posterior,sampler,n_iterations,initial=array(rep(c(0.05,0.25),n_chains),
#                                                                  dim=c(2,n_chains)),n_chains=n_chains)

samples <- monty::monty_sample(posterior,sampler,n_iterations,n_chains=n_chains)
#if we leave it like this, without the initial conditions, it starts with random values within the range of the prior

samples
names(samples)
str(samples)

#?monty::monty_sample

dim(samples$pars)
dim(samples$observations$trajectories)
object.size(samples$observations$trajectories)

library(ggplot2)
library(dplyr)
library(tidyr)

str(samples)  # Estrutura do objeto
names(samples)  # Lista dos elementos dentro de samples

########################################################################################################


######################### PLOT THE LAST TRAJECTORY OF EACH CHAIN ############################

library(dplyr)
library(tidyr)
library(ggplot2)

# Total number of chains and iterations
n_chains <- dim(samples$observations$trajectories)[4]
print(n_chains)
n_iterations <- dim(samples$observations$trajectories)[3]
print(n_iterations)
n_time <- dim(samples$observations$trajectories)[2]

# Getting the trajectories from the last iteration for all chains

trajectories <- samples$observations$trajectories[5, , n_iterations, ]  # (81 x n_chains)

print(trajectories)
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
       ylim(0, 150) +
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

start_iteration <- 5000
end_iteration <- 10000

zoomed_density <- samples$density[start_iteration:end_iteration]

zoomed_iterations <- start_iteration:end_iteration

matplot(zoomed_iterations, zoomed_density, type = "l", lty = 1,
        xlab = "Iteração", ylab = "Densidade de Probabilidade Log-Posterior",
        main = paste("Zoom na Densidade da Log-Posterior (Iterações", start_iteration, "a", end_iteration, ")"))



##############################################################################################################################

samples_df <- posterior::as_draws_df(samples)
samples_df
posterior::summarise_draws(samples_df)
#bayesplot::mcmc_pairs(samples_df)


###########################################################################################
library(ggplot2)


n_iter <- dim(samples$pars)[2]  # Number of iterations
n_chain <- dim(samples$pars)[3]  # Number of chains

# data.frame
df_samples <- data.frame(
  iter = rep(1:n_iter, times = n_chain),  
  beta = c(samples$pars["beta", , ]),  
  gamma = c(samples$pars["gamma", , ]),
  sigma = c(samples$pars["sigma", , ]), 
  eta = c(samples$pars["eta", , ]),  
 # rho = c(samples$pars["rho", , ]), 
  I0 = c(samples$pars["I0", , ]), 
  E0 = c(samples$pars["E0", , ]), 
  chain = rep(1:n_chain, each = n_iter)  # Chain index
)


########################################
# Trace plot for beta

ggplot(df_samples, aes(x = iter, y = beta, color = factor(chain))) +
  geom_line() +
  labs(title = "Evolution of the beta parameter by chain", color = "Chain") +
  labs(x = "Iteration", y = "Beta", color = "Chain") +
  theme_minimal()

start_iteration <- 5000
end_iteration <- 10000

df_samples_zoomed <- df_samples[df_samples$iter >= start_iteration & df_samples$iter <= end_iteration, ]

ggplot(df_samples_zoomed, aes(x = iter, y = beta, color = factor(chain))) +
  geom_line() +
  labs(title = paste("Evolução do parâmetro Beta por cadeia (Iterações", start_iteration, "a", end_iteration, ")"),
       x = "Iteração",
       y = "Beta",
       color = "Cadeia") +
  theme_minimal()



# Trace plot for gamma
ggplot(df_samples, aes(x = iter, y = gamma, color = factor(chain))) +
  geom_line() +
  labs(title = "Evolution of the gamma parameter by chain", color = "Chain") +
  labs(x = "Iteration", y = "Gamma", color = "Chain") +
  theme_minimal()


# Trace plot para sigma
ggplot(df_samples, aes(x = iter, y = sigma, color = factor(chain))) +
  geom_line() +
  labs(title = "Evolution of the sigma parameter by chain", color = "Chain") +
  labs(x = "Iteration", y = "Sigma", color = "Chain") +
  theme_minimal()

# Trace plot para I0
ggplot(df_samples, aes(x = iter, y = I0, color = factor(chain))) +
  geom_line() +
  labs(title = "Evolution of the I0 parameter by chain", color = "Chain") +
  labs(x = "Iteration", y = "I0", color = "Chain") +
  theme_minimal()

# Trace plot para E0
ggplot(df_samples, aes(x = iter, y = E0, color = factor(chain))) +
  geom_line() +
  labs(title = "Evolution of the E0 parameter by chain", color = "Chain") +
  labs(x = "Iteration", y = "E0", color = "Chain") +
  theme_minimal()


# Trace plot para eta
#ggplot(df_samples, aes(x = iter, y = eta, color = factor(chain))) +
#  geom_line() +
#  labs(title = "Evolution of the eta parameter by chain", color = "Chain") +
#  labs(x = "Iteration", y = "Eta", color = "Chain") +
#  theme_minimal()

# Trace plot para rho
ggplot(df_samples, aes(x = iter, y = rho, color = factor(chain))) +
  geom_line() +
  labs(title = "Evolution of the rho parameter by chain", color = "Chain") +
  labs(x = "Iteration", y = "rho", color = "Chain") +
  theme_minimal()


acf(df_samples$beta[df_samples$chain == 1],
    lag.max = 10000,
    main = "Autocorrelation of beta (chain 1)")

acf(df_samples$beta[df_samples$chain == 2],
    lag.max = 10000,
    main = "Autocorrelation of beta (chain 2)")

acf_obj <- acf(df_samples$beta[df_samples$chain == 1],
               lag.max = 10000, plot = FALSE)

plot(acf_obj, main = "Autocorrelation of beta (chain 1)", ylim = c(-0.1, 0.1))

acf(df_samples$gamma[df_samples$chain == 1],
    lag.max = 10000,
    main = "Autocorrelation of gamma (chain 1)")

acf(df_samples$gamma[df_samples$chain == 2],
    lag.max = 10000,
    main = "Autocorrelation of gamma (chain 2)")

acf_obj <- acf(df_samples$gamma[df_samples$chain == 1],
               lag.max = 10000, plot = FALSE)

plot(acf_obj, main = "Autocorrelation of gamma (chain 1)", ylim = c(-0.1, 0.1))


################ DENSITY ####################

library(ggplot2)

# Select parameters of interest
samples_subset <- samples_df[, c("beta", "gamma","sigma","rho","eta","I0","E0")]

# Convert to long format
library(tidyr)
samples_long <- pivot_longer(samples_subset, cols = everything(),
                             names_to = "parameter", values_to = "value")

# Plot the density
ggplot(samples_long, aes(x = value, fill = parameter, color = parameter)) +
  geom_density(alpha = 0.4) +
  labs(title = "Distribuição a posteriori dos parâmetros",
       x = "Valor", y = "Densidade") +
  theme_minimal()


# Select parameters of interest
samples_subset <- samples_df[, c("beta")]

# Convert to long format
library(tidyr)
samples_long <- pivot_longer(samples_subset, cols = everything(),
                             names_to = "parameter", values_to = "value")

# Plot the density
ggplot(samples_long, aes(x = value, fill = parameter, color = parameter)) +
  geom_density(alpha = 0.4) +
  labs(title = "Distribuição a posteriori de beta",
       x = "Valor", y = "Densidade") +
  theme_minimal()


################# Plot for each chain ##############

library(ggplot2)
library(tidyr)
library(dplyr)

# Mantenha a coluna .chain
samples_subset <- samples_df[, c(".chain", "beta")]

# Formato longo
samples_long <- pivot_longer(samples_subset, cols = c("beta"),
                             names_to = "parameter", values_to = "value")

# Densidade por cadeia
ggplot(samples_long, aes(x = value, color = factor(.chain), fill = factor(.chain))) +
  geom_density(alpha = 0.4) +
  facet_wrap(~parameter, scales = "free") +
  labs(title = "Distribuição dos parâmetros por cadeia",
       x = "Valor", y = "Densidade", color = "Cadeia", fill = "Cadeia") +
  theme_minimal()




#par(mfrow = c(1, 1)) # Reseta para 1 plot por vez

########################################################

# Gráfico base com os dados reais
matplot(
  x = data$time, 
  y = data$cases, 
  type = "p",                  
  col = 1, 
  pch = 16,  
  lwd = 1,  
  cex = 1.5,
  xlab = "Month", 
  ylab = "Cases SP", 
  cex.axis = 1.5,         
  cex.lab = 1.5,          
  ylim = c(0, 150) 
)

colors <- rainbow(n_chains)

# Loop sobre as cadeias
for (chain in 1:n_chains) {
  print(paste("Rodando cadeia", chain))
  
  # Parâmetros da cadeia atual
  pars2 <- list(
    N = 1000, 
    beta = as.numeric(samples$pars["beta", n_iterations, chain]), 
    gamma = as.numeric(samples$pars["gamma", n_iterations, chain]),
    sigma = as.numeric(samples$pars["sigma", n_iterations, chain]), 
    I0 = as.numeric(samples$pars["I0", n_iterations, chain]), 
    E0 = as.numeric(samples$pars["E0", n_iterations, chain]), 
    eta = as.numeric(samples$pars["eta", n_iterations, chain])
    #rho = as.numeric(samples$pars["rho", n_iterations, chain])
  )
  
  # Criar sistema
  sys2 <- dust2::dust_system_create(
    seir, pars2, time = 0, dt = 1, deterministic = FALSE, 
    n_particles = 10, n_threads = 1, seed = 3, 
    preserve_particle_dimension = TRUE
  )
  
  # Indexar e simular
  index2 <- dust2::dust_unpack_index(sys2)
  dust2::dust_system_set_state_initial(sys2)
  y2 <- dust2::dust_system_simulate(sys2, t = data$time)
  
  # Calcular log-verossimilhança
  loglik_chain <- numeric(nrow(data))
  for (i in 1:nrow(data)) {
    ll <- dust2::dust_system_compare_data(sys2, data[i, ])
    loglik_chain[i] <- sum(safe_log(ll))
  }
  logposterior[chain] <- sum(loglik_chain)
  
  # Média das trajetórias simuladas (média sobre partículas para cada tempo)
  trajectory_mean <- apply(y2[index2$incidence, , ], 2, mean)
  
  
  # Diferença absoluta ponto a ponto entre simulação e dados
  abs_diff_total[chain] <- sum(abs(trajectory_mean - data$cases))
  
  # Plotar a trajetória média da simulação
  lines(data$time, trajectory_mean, col = colors[chain], lwd = 2)
}

# Adicionar legenda
legend("topright", 
       legend = c("Real Data", paste("Chain", 1:n_chains)),   
       col = c(1, colors),                                      
       pch = c(16, rep(NA, n_chains)),                          
       lwd = c(NA, rep(2, n_chains)),                          
       cex = 1.2)











