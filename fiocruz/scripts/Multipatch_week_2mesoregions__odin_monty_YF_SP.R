

# INSTALLATION------------------------------------------------------------------
# Requires devtools package (run install.packages("devtools") if not installed)
devtools::install_github("mrc-ide/monty")
devtools::install_github("mrc-ide/dust2")
devtools::install_github("mrc-ide/odin2")

remotes::install_github("mrc-ide/odin2", upgrade = FALSE)
remotes::install_github("mrc-ide/dust2", upgrade = FALSE)
remotes::install_github("mrc-ide/odin2", upgrade = FALSE)

#devtools::install_github("mrc-ide/dust")
#devtools::install_github("mrc-ide/mcstate")

# Note the order - odin2 requires dust2 and monty to be installed
library(odin2)
library(dust2)
library(monty)
#library(dust)
#library(mcstate)

packageVersion("odin2")
packageVersion("dust2")
packageVersion("monty")

#This code (SIR Odin) was based on an example from the Odin and Monty book ("4.7 Example: an age-structured SIR model with vaccination").

#Is this the correct way to implement a nested SIR model in Odin?
sir <- odin2::odin({
  
  # Dimensions of arrays
  n_patch <- parameter(2)
  S0 <- parameter()
  I0 <- parameter()
  beta <- parameter(0.0002)
  gamma <- parameter(0.0001)
  m <- parameter()
  
  
  dim(S, S0, n_SI, p_SI) <- n_patch
  dim(I, I0, n_IR) <- n_patch
  dim(R) <- n_patch
  dim(m, s_ij) <- c(n_patch, n_patch)
  dim(lambda) <- n_patch
  dim(incidence) <- n_patch
  
  
  initial(S[]) <- S0[i]
  initial(I[]) <- I0[i]
  initial(R[]) <- 0
  initial(incidence[], zero_every = 1) <- 0
  
  
  update(S[]) <- S[i] - n_SI[i]
  update(I[]) <- I[i] + n_SI[i] - n_IR[i]
  update(R[]) <- R[i] + n_IR[i]
  #update(incidence) <- incidence + sum(n_SI)
  update(incidence[]) <- max((incidence[i] + n_SI[i]),0.1)
  
  
  # Individual probabilities of transition:
  p_SI[] <- 1 - exp(-lambda[i] * dt) # S to I
  p_IR <- 1 - exp(-gamma * dt) # I to R
  
  
 #scaling factor
  s_ij[, ] <- m[i, j] * I[j]
  lambda[] <- beta * sum(s_ij[i, ])
  

  n_SI[] <- Binomial(S[i], p_SI[i])
  n_IR[] <- Binomial(I[i], p_IR)
  
  
  #######keith's suggestion##########
  cases <- data()
  dim(cases) <- n_patch
  
  cases[] ~ Poisson(incidence[i])  
  ##########################################
  
  
})

###########################################################################################################
library(dplyr)
library(tidyr)

# carregar os dados
df <- read.csv("C:/Users/angel/OneDrive/Desktop/Angelica_R/multipatch/cases_por_populacao_regioes_SP_semana.csv",
               stringsAsFactors = FALSE)

# garantir ordenação
df <- df %>% arrange(week, population)

# reorganizar os dados em formato wide
df_wide <- df %>%
  select(week, population, cases) %>%
  pivot_wider(names_from = population, values_from = cases)

# transformar cada linha em vetor (lista de vetores)
incidence_list <- split(as.matrix(df_wide[,-1]), seq(nrow(df_wide)))

# montar dataframe final
data <- data.frame(
  time = df_wide$week,
  cases = I(incidence_list)  # lista de vetores
)

# verificar estrutura
str(data)

print(data)

#######################################


#library(dplyr)

#data <- df_tidy %>%
#  rename(time = day,      # já estava renomeando
#         group = population) %>%  # renomeia population -> group
#  mutate(group = as.character(group)) # garante que seja character

#str(data)

#data$group <- data$group

#How should we declare the set of parameters?
pars <- list(
  n_patch = 2,
  I0    = c(2, 2),
  S0    = c(500, 500),
  gamma = 0.0002,
  beta  = 0.0001,
  m = matrix(c(
    1.0, 0.1,
    0.1, 1.0
  ), nrow = 2, ncol = 2, byrow = TRUE)
)


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

sys <- dust2::dust_system_create(generator = sir, pars = pars,  time = 0, dt =1, 
                                 deterministic = FALSE, n_particles = 200, n_threads = 1, seed = 1, preserve_particle_dimension = TRUE)


t <- 0:max(data$time)
dust2::dust_system_set_state_initial(sys)
y <- dust2::dust_system_simulate(sys, t)

dim(y) #An array with 3 dimensions (state x particle x time) - no nosso caso, 8 estados, S,I,R e incidence para cada patch
y[, 1, ] #tentar identificar qual estado é S,I,R e incidence de cada patch
head(y)                    
# Convert model output to labelled values
results <- dust2::dust_unpack_state(sys, y)

particle <- 200
times <- 0:max(data$time)

# Patch 1
states_patch1 <- y[, particle, ]  # todas as variáveis de estado
# supondo que S, I, R estão nas primeiras 3 linhas
S1 <- states_patch1[1, ]
I1 <- states_patch1[3, ]
R1 <- states_patch1[5, ]

# Patch 2
S2 <- states_patch1[2, ]  # ajuste se suas linhas são diferentes
I2 <- states_patch1[4, ]
R2 <- states_patch1[6, ]

# Patch 1
matplot(times, cbind(S1, I1, R1), type="l", lty=1, col=c("blue","red","green"),
        xlab="Time", ylab="Individuals", main="Patch 1")
legend("topright", legend=c("S","I","R"), col=c("blue","red","green"), lty=1)

# Patch 2
matplot(times, cbind(S2, I2, R2), type="l", lty=1, col=c("blue","red","green"),
        xlab="Time", ylab="Individuals", main="Patch 2")
legend("topright", legend=c("S","I","R"), col=c("blue","red","green"), lty=1)

particle <- 200
times <- 0:max(data$time)

inc_patch1 <- y[7, particle, ]  
inc_patch2 <- y[8, particle, ]

library(ggplot2)
library(reshape2)

# --- Número de patches e partículas ---
particle <- 200
times <- 0:(length(inc_patch1)-1)

# --- Criar df_sim_long com coluna 'type' ---
df_sim <- data.frame(
  time = times,
  patch1 = inc_patch1,
  patch2 = inc_patch2
)
df_sim_long <- melt(df_sim, id.vars = "time", variable.name = "patch", value.name = "cases")
df_sim_long$type <- "sim"

# --- Criar df_real_long com coluna 'type' ---
n_patch <- length(data$cases[[1]])  # número de patches
df_real <- data.frame(time = data$time)
for (p in 1:n_patch) {
  df_real[[paste0("patch", p)]] <- sapply(data$cases, function(x) x[p])
}
df_real_long <- melt(df_real, id.vars = "time", variable.name = "patch", value.name = "cases")
df_real_long$type <- "real"

# --- Combinar ---
df_plot <- rbind(df_real_long, df_sim_long)

# --- Plot com linetype mapeado ---
ggplot(df_plot, aes(x = time, y = cases, color = patch, linetype = type)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("patch1" = "red", "patch2" = "blue")) +
  scale_linetype_manual(values = c("real" = "solid", "sim" = "dashed")) +
  theme_minimal() +
  labs(x = "Time", y = "Incidence", title = "Simulação vs Dados Reais") +
  theme(legend.title = element_blank())

################################################################################################

#In dust_filter_create, should we include information about groups/populations?
# Create filter
filter <- dust2::dust_filter_create(generator=sir, data = data,  time_start = 0, n_threads = 1, n_particles = 500)
#filter <- dust2::dust_filter_create(sir, data = data, time_start = 0, n_particles = 200)

dust2::dust_likelihood_run(filter, pars)

library(ggplot2)

#Should we use monty::monty_packer or monty::monty_packer_grouped?
#packer <- monty::monty_packer(c("beta","gamma"),fixed=list(S0 = 1000,I0=10))  


packer <- monty::monty_packer(scalar = c("beta", "gamma"),
                              fixed = list(S0 = c(500, 500), I0 = c(2, 2),   
                                           m = matrix(c(
                                             1.0, 0.1,
                                             0.1, 1.0
                                           ), nrow = 2, ncol = 2, byrow = TRUE))
                              )


likelihood <- dust2::dust_likelihood_monty(filter, packer, save_trajectories = TRUE) #save_trajectories = TRUE is optional and requires much memory


#We also have doubts about the priors, vcv, the sampler, monty_sample...
prior <- monty::monty_dsl({
  
  beta ~ Uniform(0.00005, 0.001)
  gamma ~ Uniform(0.00005, 0.001)
  
})


#######################################################################################################
posterior <- likelihood + prior


posterior

#vcv <- matrix(0.005, nrow = 2, ncol = 2)

#diag(vcv) <- 0.01

#testando outras matrizes vcv
vcv <- matrix(0.01, nrow = 2, ncol = 2)

diag(vcv) <- 0.01

print(vcv)

str(vcv)


sampler <- monty::monty_sampler_random_walk(vcv)

n_chains=5
n_iterations=1000

samples <- monty::monty_sample(posterior,sampler,n_iterations,n_chains=n_chains)


######################### PLOT THE LAST TRAJECTORY OF EACH CHAIN############################

library(dplyr)
library(tidyr)
library(ggplot2)

n_chains <- dim(samples$observations$trajectories)[4]
n_iterations <- dim(samples$observations$trajectories)[3]
n_time <- dim(samples$observations$trajectories)[2]


print(n_chains)
print(n_iterations)
print(n_time)

times <- 1:n_time

library(dplyr)
library(tidyr)
library(ggplot2)

# índices de incidência para cada patch no objeto traj
# (ajuste conforme sua ordem real no modelo)
inc_indices <- c(7,8)  

# Criar dataframe longo para cada patch
df_sim_list <- list()
for (patch in 1:2) {
  df_chain_list <- list()
  for (chain in 1:n_chains) {
    # Pegar a última iteração de cada chain
    traj <- samples$observations$trajectories[, , n_iterations, chain]
    
    # extrair incidência do patch correspondente
    inc <- traj[inc_indices[patch], ]
    
    df_chain_list[[chain]] <- data.frame(
      time = times,
      cases = inc,
      chain = paste0("Chain", chain),
      patch = paste0("patch", patch)
    )
  }
  df_sim_list[[patch]] <- bind_rows(df_chain_list)
}

df_sim_all <- bind_rows(df_sim_list)

# Transformar dados reais em formato longo para 4 patches
df_real_long <- data.frame(
  time = data$time,
  patch1 = sapply(data$cases, function(x) x[1]),
  patch2 = sapply(data$cases, function(x) x[2])
) %>%
  pivot_longer(cols = -time, names_to = "patch", values_to = "cases")

# Plot automático para cada patch
ggplot() +
  geom_line(data = df_real_long,
            aes(x = time, y = cases),
            color = "black", linetype = "dashed", size = 1) +
  geom_line(data = df_sim_all,
            aes(x = time, y = cases, color = chain),
            alpha = 0.7,
            size = 0.8) +
  scale_color_manual(values = rainbow(n_chains)) +
  facet_wrap(~patch, scales = "free_y") +
  scale_y_continuous(limits = c(0, 40)) +
  labs(title = "Trajectory vs Real Data - Each Patch",
       x = "Time (week)", y = "Incidence", color = "Chain") +
  theme_minimal()

library(ggplot2)
library(dplyr)

# Cores do rainbow para os 4 patches
patch_colors <- rainbow(2)
names(patch_colors) <- paste0("patch", 1:2)

# Filtrar apenas a chain 2 para simulação
df_sim_chain2 <- df_sim_all %>% filter(chain == "Chain2")

# Criar coluna 'type' para legendas
df_real_long <- df_real_long %>% mutate(type = "Dados Reais")
df_sim_chain2 <- df_sim_chain2 %>% mutate(type = "Chain 2")

# Plot
ggplot() +
  # Dados reais como linha tracejada
  geom_line(data = df_real_long,
            aes(x = time, y = cases, color = patch, linetype = type),
            size = 1.2) +
  # Simulação Chain 2 como linha contínua
  geom_line(data = df_sim_chain2,
            aes(x = time, y = cases, color = patch, linetype = type),
            size = 1.2) +
  scale_color_manual(values = patch_colors) +
  scale_linetype_manual(values = c("Dados Reais" = "twodash",
                                   "Chain 2" = "solid")) +
  scale_y_continuous(limits = c(0, 125)) +
  labs(title = "Simulação vs Dados Reais - 2 Patches",
       x = "Tempo", y = "Incidência",
       color = "Patch", linetype = "Tipo") +
  theme_minimal()

#################################Log posterior probability density################################
par(mfrow = c(1, 1)) # Reseta para 1 plot por vez

matplot(samples$density, type = "l", lty = 1,
        xlab = "Sample", ylab = "Log posterior probability density")

print(samples$density)


matplot(samples$density, type = "l", lty = 1,
        xlab = "Sample", ylab = "Log posterior probability density",
        col = rainbow(ncol(samples$density)))
legend("bottomright", legend = paste("Chain", 1:ncol(samples$density)),
       col = rainbow(ncol(samples$density)), lty = 1)


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
  beta = c(samples$pars["beta", , ]),  
  gamma = c(samples$pars["gamma", , ]), 
  #I0 = c(samples$pars["I0", , ]),  
  #eta = c(samples$pars["eta", , ]),  
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
  scale_y_continuous(limits = c(0.0, 0.001))  


# Trace plot for gamma
ggplot(df_samples, aes(x = iter, y = gamma, color = factor(chain))) +
  geom_line() +
  labs(title = "Evolution of the gamma parameter by chain", color = "Chain") +
  labs(x = "Iteration", y = "Gamma", color = "Chain") +
  theme_minimal() +
  scale_y_continuous(limits = c(0.0, 0.001))  



library(ggplot2)

# Selecione apenas os parâmetros de interesse
samples_subset <- samples_df[, c("beta")]

# Converta para long format
library(tidyr)
samples_long <- pivot_longer(samples_subset, cols = everything(),
                             names_to = "parameter", values_to = "value")

# Plote a densidade
ggplot(samples_long, aes(x = value, fill = parameter, color = parameter)) +
  geom_density(alpha = 0.4) +
  labs(title = "Distribuição a posteriori dos parâmetros",
       x = "Valor", y = "Densidade") +
  theme_minimal()

head(samples_df)

# Selecione apenas os parâmetros de interesse
samples_subset <- samples_df[, c("gamma")]

# Converta para long format
library(tidyr)
samples_long <- pivot_longer(samples_subset, cols = everything(),
                             names_to = "parameter", values_to = "value")

# Plote a densidade
ggplot(samples_long, aes(x = value, fill = parameter, color = parameter)) +
  geom_density(alpha = 0.4) +
  labs(title = "Distribuição a posteriori dos parâmetros",
       x = "Valor", y = "Densidade") +
  theme_minimal()

head(samples_df)

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



