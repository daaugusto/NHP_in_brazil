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


sir <- odin({
  deriv(S) <- -beta * S * I / N
  deriv(I) <- beta * S * I / N - gamma * I
  deriv(R) <- gamma * I
  deriv(incidence) <- beta * S * I / N
  
  initial(S) <- N - I0
  initial(I) <- I0
  initial(R) <- 0
  initial(incidence, zero_every = 1) <- 0
  
  N <- parameter(1000)
  I0 <- parameter(10)
  beta <- parameter(0.2)
  gamma <- parameter(0.1)
  
  # observation model
  rho <- parameter(0.1)
  eta <- parameter(10)
  cases <- data()
  cases ~ NegativeBinomial(size = eta, mu = rho * incidence)
})

mod <- dust_system_create(sir,
                          pars = list())



pars <- list(
  rho = 0.1,
  eta = 10,
  N = 1000,
  beta=0.4,
  gamma=0.1
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


filter <- dust_unfilter_create(sir, data = data, time_start = 0)
dust_likelihood_run(filter, pars)


packer <- monty_packer(c("beta", "gamma", "rho", "eta"),
                       fixed = list(I0 =pars$I0,
                                    N = pars$N))

likelihood <- dust_likelihood_monty(filter, packer)


prior <- monty_dsl({
  beta ~ Uniform(min=0.2,max=0.6)
  gamma ~ Uniform(min=0.01,max=0.20)
  
  rho ~ Uniform(0, 1)
  eta ~ Uniform(10, 1000)
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

mod_mult <- dust_system_create(sir,
                               pars = list_par, n_groups = n_groups)

t <- c(0, data$time)
dust_system_set_state_initial(mod_mult)
y <- dust_system_simulate(mod_mult, t)
y <- dust_unpack_state(mod_mult, y)

matplot(t, t(y$incidence), type = "l", lty = 1, col = "#00008822",
        xlab = "Time", ylab = "Infected population")
points(data$time, data$cases / mean(samples$pars["rho", , ]), pch = 19, col = "red")




