# Read the data
#original_data <- rowSums(read.csv(file = "C:/Users/angel/OneDrive/Desktop/Angelica_R/Newcodes/flu2009.csv", sep = ";"))

#Plot with our data
data<- read.csv("C:/Users/angel/OneDrive/Desktop/Angelica_R/sir_temp_example_epizotias/cases_per_state_pnh_sp_mg_pr.csv", header= TRUE, sep = ",") # this is just example data from the package
#data<- read.csv("fiocruz/data/cases_per_state_pnh_sp_mg_pr.csv", header= TRUE, sep = ",") # this is just example data from the package

data_week_selecionado <- data[, c(1,2)] #para pegar os dados apenas de SP por enquanto

data_week_selecionado$week <- ceiling(data_week_selecionado$day/30)

data<- aggregate(cases_A ~ week, data = data_week_selecionado, sum)
head(data)

library(dplyr)

data <- data %>%
  rename(cases = cases_A)

data <- data %>%
  rename(time = week)



# Create the data frame and add days as time
#data <- data.frame(
#  time = seq(7, by = 7, length.out = length(original_data)),
#  cases = original_data
#)

plot(data$time, data$cases, pch = 19, col = "red", 
     xlab = "Days since start of epidemics",
     ylab = "Official cases count",
     main = "Cases over Time with Holiday Periods",
     ylim = c(0, max(data$cases) * 1.3), # Add padding to the y-axis
     xlim = c(min(data$time), max(data$time)))

library(odin2)
library(dust2)
library(monty)

seir <- odin2::odin({
  # initial conditions
  initial(S) <- (1 - 2 * alpha) * N
  initial(E) <- alpha * N
  initial(I) <- alpha * N
  initial(R) <- 0
  initial(incidence, zero_every = 1) <- 0
  
  # equations
  deriv(S) <- - hol * beta * S * I / N
  deriv(E) <- hol * beta * S * I / N - sigma * E
  deriv(I) <- sigma * E - gamma * I
  deriv(R) <- gamma * I
  deriv(incidence) <- sigma * E
  
  # parameter values
  R_0 <- parameter(1.5)
  L <- 1
  D <- 1.25
  alpha <- parameter(1e-4) # initial proportion
  N <- parameter(1000) # total population
  
  # convert parameters
  hol <- interpolate(h_times, h_values, "constant")
  h_times <- parameter()
  h_values <- parameter()
  dim(h_times) <- parameter(rank = 1)
  dim(h_values) <- length(h_times)
  gamma <- 1 / L
  sigma <- 1 / D
  beta <- R_0 * gamma
  
  # observation model
  rho <- parameter(0.1)
  eta <- parameter(10)
  cases <- data()
  cases ~ NegativeBinomial(size = eta, mu = rho * incidence)
})


hol_t <- c(0,28, 56, 84)         
hol_v <- c(1.0, 0.8, 0.5, 1.0)      



mod <- dust_system_create(seir,
                          pars = list(h_times = hol_t, h_values = hol_v))



pars <- list(
  alpha = 2e-5,
  R_0 = 1.3,
  rho = 0.1,
  eta = 10,
  N = 1000,
  h_times = hol_t,
  h_values= hol_v
)
dust_system_update_pars(sys = mod, pars = pars)

t <- c(0,data$time)
dust_system_set_state_initial(mod)
y <- dust_system_simulate(mod, t)


plot(t, dust_unpack_state(mod, y)$incidence, type = "l", col = "#000088ff")



plot(data$time, data$cases,
     xlab = "Days since start of epidemics",
     ylab = "Official estimates of cases",
     pch = 19,
     col = "red",
     ylim = c(0, max(data$cases) * 1.5))
lines(t, dust_unpack_state(mod, y)$incidence * pars$rho, col = "#000088ff")


filter <- dust_unfilter_create(seir, data = data, time_start = 0)
dust_likelihood_run(filter, pars)


packer <- monty_packer(c("alpha", "R_0", "rho", "eta"),
                       fixed = list(h_times = pars$h_times,
                                    h_values = pars$h_values,
                                    N = pars$N))

likelihood <- dust_likelihood_monty(filter, packer)


prior <- monty_dsl({
  alpha ~ Beta(a = 4e-4, b = 2)
  R_0 ~ Gamma(2, 0.7)
  rho ~ Uniform(0, 1)
  eta ~ Exponential(mean = 1000)
})

posterior <- likelihood + prior

vcv <- diag(c(5e-10, 5e-5, 1e-5, 1))
sampler <- monty_sampler_adaptive(vcv, initial_vcv_weight = 100)

samples <- monty_sample(posterior,
                        sampler,
                        1000,
                        initial = packer$pack(pars),
                        n_chains = 3)


matplot(samples$density, type = "l", lty = 1,
        xlab = "Sample", ylab = "Log posterior probability density")


length(unique(samples$density)) / length(samples$density)


samples_df <- posterior::as_draws_df(samples)
posterior::summarise_draws(samples_df)

n_groups <- 1000
i_rand <- sample(dim(samples$pars)[2], n_groups)
list_par <- lapply(i_rand, FUN = function(i) packer$unpack(samples$pars[, i, sample(3, 1)]))

mod_mult <- dust_system_create(seir,
                               pars = list_par, n_groups = n_groups)

t <- c(0, data$time)
dust_system_set_state_initial(mod_mult)
y <- dust_system_simulate(mod_mult, t)
y <- dust_unpack_state(mod_mult, y)

matplot(t, t(y$incidence), type = "l", lty = 1, col = "#00008822",
        xlab = "Time", ylab = "Infected population")
points(data$time, data$cases / mean(samples$pars["rho", , ]), pch = 19, col = "red")

