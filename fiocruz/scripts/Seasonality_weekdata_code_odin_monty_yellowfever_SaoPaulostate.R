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


# INSTALLATION------------------------------------------------------------------
# Requires devtools package (run install.packages("devtools") if not installed)
devtools::install_github("mrc-ide/monty")
devtools::install_github("mrc-ide/dust2")
devtools::install_github("mrc-ide/odin2")
devtools::install_github("mrc-ide/dust")
# Note the order - odin2 requires dust2 and monty to be installed
library(odin2)
library(dust2)
library(monty)
library(dust)

packageVersion("odin2")
packageVersion("dust2")
packageVersion("monty")

sir <- odin2::odin({
  
  N <- parameter(1000)
  I0 <- parameter(10)
  
  Temp <- interpolate(temp_times, temp_values, "spline")
  temp_times <- parameter()
  temp_values <- parameter()
  dim(temp_times) <- parameter(rank = 1)
  dim(temp_values) <- length(temp_times)
  
  temp_scale <- parameter(0.002)  
  
  beta <- temp_scale * Temp * (Temp - 20.0) * sqrt(30.0 - Temp)  # Brière
  gamma <- parameter(0.1)
  
  p_SI <- max(0, min(1, 1 - exp(-beta * I / N * dt))) # Note parameter dt - this is a time interval set when the model is set up via dust_system_create
  p_IR <- max(0, min(1, 1 - exp(-gamma * dt)))
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
  
  
  # cases <- data()
  # cases ~ Poisson(incidence)
  
  # Negative Binomial
  #rho <- parameter(0.5)
  eta <- parameter(10)
  cases <- data()
  #cases ~ NegativeBinomial(size = eta, mu = incidence) #keith's paper
  # cases ~ NegativeBinomial(max(cases, 0.1)  , mu = incidence) #tinha entendido errado no paper do Keith
  cases ~ NegativeBinomial(size = eta, mu = incidence) #odin and monty book's suggestion
  
  # print("temperatura: {Temp}") 
})

###########################################################################################################


temp_values <- c(
  25.21466909497379, 24.55025728492137, 22.213488667900094, 20.988407724329324, 21.045301900246685,
  22.28036106614246, 23.701607308048104, 25.27081406105458, 26.017891901788467, 26.30525603222325,
  26.336332581714462, 27.15742079093432, 25.985560245143386, 25.760113128276288, 21.821806101603453,
  19.13791724483503, 20.78495991366019, 21.341196423065064, 22.087814138143695, 24.07351410730805,
  25.314309185168053, 26.00403031529448, 26.408926437711994, 27.040500693802034, 25.630931043786617,
  24.02696673604687, 22.99182132670367, 20.449068185322233, 19.825970118331792, 21.31176283533765,
  23.743856961147085, 24.66359032917052, 25.065959374036385, 25.8045405488745, 26.21587843046562,
  25.757415009250693, 26.738706444650017, 24.65603559589886, 22.872687326549492, 21.502746299722478,
  21.66973577705828, 20.648379683163736, 22.815084894387912, 24.28058028831329, 25.34891015263645,
  26.32100389300031, 27.633759250693803, 26.423966042244835, 25.621241905642922, 25.14059127351218,
  23.68020544249152, 21.379808433549183, 20.180022355843356, 21.310100601295098, 23.815164392537774,
  25.29953698350293, 25.740106093894543, 25.403898781991984, 26.192563791242677, 25.50753546099291,
  24.717261216466234, 23.128666551032993, 20.55941643539932, 21.05048373419673, 20.759532550878816,
  20.681542456830094, 24.054145467160037, 24.932667476102374, 25.215348442799876, 25.643612203206906,
  26.436680928152946, 25.77178249306198, 25.726087920135676, 23.18066797718162, 21.9918719164354,
  19.948338247764415, 18.47781881167129, 21.62431342506938, 24.110266342892384, 23.396864400246685,
  25.346385484119644
)

#temp_times <- seq(0, 80) #for month data
temp_times <- c(0, seq(4, 345, by = 4.3))
print(temp_times)
length(temp_times)

plot(temp_times, temp_values, type = "l", col = "red", lwd = 2,
     xlab = "Tempo", ylab = "Temperatura",
     main = "Temperatura ao longo do tempo")

points(temp_times, temp_values, col = "red", pch = 19, cex = 1.2)

# load example data
#data<- read.csv("C:/Users/angel/OneDrive/Desktop/Angelica_R/sir_temp_example_epizotias/cases_per_state_pnh_sp_mg_pr.csv", header= TRUE, sep = ",") # this is just example data from the package
#data <- read.csv(system.file("sir_incidence.csv", package="mcstate")) # this is just example data from the package 
data<- read.csv("fiocruz/data/cases_per_state_pnh_sp_mg_pr.csv", header= TRUE, sep = ",") # this is just example data from the package

data_week_selected <- data[, c(1,2)] #to get data only from SP for now
head(data_week_selected)

plot(data_week_selected[,1], data_week_selected[,2], 
     xlab = "day",  
     ylab = "Number of Cases SP",  
     pch = 16, ylim = c(0,30), cex =0.7)                      

data_week_selected$week <- ceiling(data_week_selected$day/7) #grouping monthly
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

plot(data[,1], data[,2], 
     xlab = "month",  
     ylab = "Number of Cases SP",  
     pch = 16, ylim = c(0,150), cex =0.7)  

pars <- list(
  N = 1000,
  I0 = 10,
  gamma = 0.1,
  eta = 10,
  temp_scale = 0.002,  # valor inicial
  temp_times = temp_times,
  temp_values = temp_values
)

#pars <- list(N = 1000, I0 =10, beta=0.4, gamma = 0.1, eta = 10)
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

t <- c(0,data$time) #mudar de acordo com nossos dados (numero de meses ou semanas)
print(t)

sys <- dust2::dust_system_create(generator = sir, pars = pars, time = 0, dt =1, 
                                 deterministic = FALSE, n_particles = 200, n_threads = 1, seed = 1, preserve_particle_dimension = TRUE)


# The output time sequence has the same units as the calculation time interval dt
# Note that the separation of the output time points needs to be compatible with dt
# For example, if dt = 1, the output time sequence cannot have fractional points
#t <- seq(1, 100, by = 1)  #mudar de acordo com nossos dados (numero de meses ou semanas)


#length(data)
# Initialize model to starting conditions
dust2::dust_system_set_state_initial(sys)

t <- c(0,data$time) #mudar de acordo com nossos dados (numero de meses ou semanas)
print(t)
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
filter <- dust2::dust_filter_create(sir, data = data, time_start = 0, n_threads = 1, n_particles = 1)
#filter <- dust2::dust_filter_create(sir, data = data, time_start = 0, n_particles = 200)

dust2::dust_likelihood_run(filter, pars)

library(ggplot2)

#data created from interpolate when we run using Print inside sir.
time_interpolate = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 
                     39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 
                     75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 
                     109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 
                     138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 
                     167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 
                     196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 
                     225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 
                     254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 
                     283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 
                     312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 
                     341, 342, 343)

temperature_interpolate <- c(25.214669, 25.153037, 25.049617, 24.862620, 24.550257, 24.088789, 23.526664, 22.930377, 22.366427, 21.896630, 21.534500, 21.265054, 
                             21.072941, 20.943057, 20.869880, 20.860337, 20.922185, 21.063180, 21.284156, 21.565396, 21.883383, 22.214600, 22.538512, 22.856641, 
                             23.180391, 23.521214, 23.889800, 24.279286, 24.665261, 25.022552, 25.326038, 25.562250, 25.743708, 25.886451, 26.006515, 26.116578, 
                             26.211135, 26.278565, 26.307241, 26.287605, 26.241042, 26.212761, 26.248581, 26.393700, 26.645087, 26.917967, 27.119643, 27.157421, 
                             26.973212, 26.647354, 26.294791, 26.030467, 25.954451, 26.013247, 26.062763, 25.957738, 25.555679, 24.821806, 23.861263, 22.788540, 
                             21.718124, 20.752920, 19.961435, 19.405816, 19.148208, 19.234285, 19.587965, 20.078596, 20.575268, 20.950889, 21.166175, 21.269648, 
                             21.313646, 21.350465, 21.422759, 21.551673, 21.755443, 22.052306, 22.452675, 22.924656, 23.422118, 23.898916, 24.311168, 24.648795, 
                             24.927742, 25.164624, 25.375954, 25.570483, 25.743797, 25.890209, 26.004030, 26.084894, 26.153708, 26.236702, 26.360105, 26.545180, 
                             26.761931, 26.950122, 27.049128, 26.999448, 26.785442, 26.448438, 26.033574, 25.585982, 25.145113, 24.733531, 24.370682, 24.076009, 
                             23.861723, 23.686579, 23.485374, 23.192791, 22.746648, 22.156848, 21.505383, 20.877378, 20.357821, 20.002964, 19.804937, 19.747198,
                             19.813201, 19.987846, 20.263836, 20.636502, 21.101176, 21.651469, 22.255213, 22.860398, 23.414502, 23.865372, 24.189312, 24.410855,
                             24.559212, 24.663590, 24.750126, 24.832666, 24.921987, 25.028862, 25.163001, 25.323094, 25.501329, 25.689814, 25.880426, 26.056181,
                             26.188580, 26.248359, 26.206267, 26.058113, 25.874149, 25.738392, 25.734861, 25.929094, 26.250075, 26.565581, 26.743101, 26.655096,
                             26.288386, 25.744156, 25.128560, 24.547560, 24.066287, 23.658770, 23.286709, 22.911808, 22.507182, 22.107685, 21.768949, 21.546619,
                             21.490758, 21.567853, 21.680057, 21.727867, 21.612843, 21.318784, 20.968888, 20.699853, 20.648380, 20.908880, 21.406610, 22.024539,
                             22.645635, 23.161724, 23.556082, 23.865934, 24.129203, 24.383282, 24.644982, 24.904375, 25.149746, 25.369389, 25.561185, 25.751511,
                             25.972013, 26.254336, 26.620436, 27.020663, 27.373273, 27.596370, 27.611429, 27.417435, 27.090881, 26.711634, 26.359430, 26.086771,
                             25.885399, 25.738833, 25.630593, 25.543727, 25.458736, 25.355263, 25.212951, 25.012163, 24.744051, 24.408070, 24.003889, 23.531323,
                             23.001523, 22.444854, 21.893541, 21.379808, 20.933491, 20.574859, 20.321791, 20.192166, 20.202212, 20.351132, 20.628079, 21.022080,
                             21.521785, 22.101306, 22.715865, 23.319422, 23.865947, 24.322765, 24.696885, 25.002650, 25.254405, 25.464154, 25.626603, 25.728709,
                             25.757389, 25.701045, 25.582203, 25.457513, 25.385110, 25.422947, 25.591082, 25.825008, 26.048782, 26.186456, 26.177333, 26.043184,
                             25.833538, 25.597941, 25.383476, 25.200311, 25.030200, 24.854169, 24.653208, 24.405883, 24.086649, 23.669558, 23.128667, 22.460519,
                             21.751623, 21.110976, 20.647577, 20.458186, 20.513210, 20.708518, 20.939014, 21.100955, 21.143055, 21.082179, 20.939748, 20.737192,
                             20.512108, 20.350119, 20.345726, 20.593433, 21.166443, 21.980540, 22.880958, 23.712592, 24.325182, 24.679823, 24.848964, 24.909895,
                             24.939789, 24.990916, 25.059986, 25.136188, 25.208714, 25.270728, 25.336889, 25.429090, 25.569229, 25.776247, 26.024840, 26.255645,
                             26.408422, 26.423569, 26.290834, 26.083621, 25.883435, 25.771782, 25.798413, 25.886052, 25.925671, 25.808241, 25.437782, 24.853032,
                             24.172204, 23.514542, 22.997193, 22.655968, 22.421009, 22.215399, 21.962238, 21.606400, 21.157472, 20.637000, 20.066535, 19.474186,
                             18.936580, 18.552083, 18.419164, 18.631946, 19.184540, 19.971047, 20.881221, 21.804886, 22.647274, 23.347979, 23.851249, 24.101330,
                             24.066456, 23.844598, 23.577395, 23.406513, 23.466056, 23.776900, 24.272763, 24.885118)



df_interpolate <- data.frame(time = time_interpolate, temperature = temperature_interpolate)
df_nova <- data.frame(time = temp_times, temperature = temp_values)

# --- Plotar os dados usando ggplot2 ---
p <- ggplot() +
  geom_line(data = df_original, aes(x = time, y = temperature, color = "Série Original"), size = 1.2) +
  geom_point(data = df_original, aes(x = time, y = temperature, color = "Série Original"), size = 2) +
  
  geom_point(data = df_nova, aes(x = time, y = temperature, color = "Nova Série"), size = 2) +
  
  labs(
    title = "Comparação de Duas Séries Temporais de Temperatura",
    x = "Tempo (semana)",
    y = "Temperatura (°C)",
    color = "Série de Dados"
  ) +
  
  scale_color_manual(values = c("Série Original" = "blue", "Nova Série" = "red")) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(p)


# Create packer - divide input parameters into estimated (beta + gamma) and fixed (N + I0)
# Significantly simplified from mcstate!
packer <- monty::monty_packer(c("temp_scale","gamma","eta","I0"),fixed=list(N = 1000,temp_times = temp_times, 
                                                                            temp_values = temp_values))  

likelihood <- dust2::dust_likelihood_monty(filter, packer, save_trajectories = TRUE) #save_trajectories = TRUE is optional and requires much memory

likelihood$properties

# Set prior likelihood distributions for estimated parameters
# Here a simple uniform distribution is used for beta and gamma, with permitted minimum/maximum values
# See https://mrc-ide.github.io/odin2/articles/functions.html#distribution-functions for available distributions
prior <- monty::monty_dsl({
  
  temp_scale ~ Uniform(0.0003, 0.02)
  gamma ~ Uniform(0.1, 5.0)
  
  #beta ~ Exponential(mean = 0.3)
  # gamma ~ Exponential(mean = 0.1)
  
  # beta ~ Normal(0.5, 0.05)
  # gamma ~ Normal(0.2, 0.01)
  
  #beta ~TruncatedNormal(1.5, 1.5, min = 0.1, max = 4)  
  #gamma ~TruncatedNormal(1.5, 1.5, min = 0.1, max = 4) 
  
  eta ~ Exponential(mean = 1000)

  I0 ~ Uniform(2,10)
  
})

prior$properties
prior$domain #Check limits (can adjust manually by adjusting values in prior$domain)

###################################################################################################

n_streams <- 100
r <- monty::monty_rng_create(n_streams = n_streams)
prior_samples <- matrix(monty::monty_model_direct_sample(prior, r), nrow = n_streams)
colnames(prior_samples) <- prior$parameters


#######################################################################################################
posterior <- likelihood + prior

posterior

pars <- list(temp_scale = 0.05, gamma=0.3, eta=10, I0=10)
no_param <- length(pars)

vcv <- matrix(0.005, nrow = 4, ncol = 4)
diag(vcv) <- 0.01
print(vcv)
str(vcv)

print(class(posterior$system))
print(class(posterior$system$methods))


sampler <- monty::monty_sampler_random_walk(vcv)


n_chains=5
n_iterations=10000

#samples <- monty::monty_sample(posterior,sampler,n_iterations,initial=array(rep(c(0.05,0.25),n_chains),
#                                                                  dim=c(2,n_chains)),n_chains=n_chains)
#samples <- monty::monty_sample(posterior,s,n_iterations,n_chains=n_chains)

samples <- monty::monty_sample(posterior,sampler,n_iterations,n_chains=n_chains)
#if we leave it like this, without the initial conditions, it starts with random values within the range of the prior


plot(dexp(1000))

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
length(pars)
dim(vcv)

posterior$parameters
names(pars)
dim(vcv)

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


##################### AREA ##############################

# Número total de cadeias, iterações e tempo
n_chains <- dim(samples$observations$trajectories)[4]
n_iterations <- dim(samples$observations$trajectories)[3]
n_time <- dim(samples$observations$trajectories)[2]

cases_real <- data$cases

# Para guardar resultados por cadeia
results <- data.frame(chain = 1:n_chains,
                      mean_diff_raw = NA,
                      sd_diff_raw = NA)
# mean_diff_scaled = NA,
# sd_diff_scaled = NA)

for (chain in 1:n_chains) {
  diffs_raw <- numeric(n_iterations)
  #diffs_scaled <- numeric(n_iterations)
  
  for (iter in 1:n_iterations) {
    # Trajetória no tempo para essa iteração e cadeia
    trajectory_incidence <- samples$observations$trajectories[4, , iter, chain]
    
    # Soma das diferenças absolutas em relação aos dados reais
    diffs_raw[iter] <- sum(abs(trajectory_incidence - cases_real))
  }
  
  # Guardar média e desvio padrão
  results$mean_diff_raw[chain] <- mean(diffs_raw)
  results$sd_diff_raw[chain] <- sd(diffs_raw)

}

print(results)

# Soma das médias ponderadas pelas iterações (ou só soma direta das médias)
mean_diff_raw_total <- mean(results$mean_diff_raw)
sd_diff_raw_total <- sqrt(sum(results$sd_diff_raw^2)) / length(results$sd_diff_raw)  # Uma aproximação simples

cat("\n--- RESULTADO FINAL (média entre cadeias) ---\n")
cat(sprintf("Média das diferenças absolutas: %.2f ± %.2f\n", mean_diff_raw_total, sd_diff_raw_total))

# Para guardar todos os valores de diferenças absolutas de todas as iterações e cadeias
all_diffs_raw <- c()

for (chain in 1:n_chains) {
  diffs_raw <- numeric(n_iterations)
  
  for (iter in 1:n_iterations) {
    trajectory_incidence <- samples$observations$trajectories[4, , iter, chain]
    diffs_raw[iter] <- sum(abs(trajectory_incidence - cases_real))
  }
  
  # Acumula todos os valores para o cálculo final
  all_diffs_raw <- c(all_diffs_raw, diffs_raw)
  
  # Estatísticas por cadeia (se ainda quiser manter)
  results$mean_diff_raw[chain] <- mean(diffs_raw)
  results$sd_diff_raw[chain] <- sd(diffs_raw)
}

mean_diff_raw_total <- mean(all_diffs_raw)
sd_diff_raw_total <- sd(all_diffs_raw)

cat("\n--- RESULTADO FINAL (todas as cadeias e iterações juntas) ---\n")
cat(sprintf("Média das diferenças absolutas: %.2f ± %.2f\n",
            mean_diff_raw_total, sd_diff_raw_total))


#################################Log posterior probability density################################
par(mfrow = c(1, 1)) # Reseta para 1 plot por vez

matplot(samples$density, type = "l", lty = 1,
        xlab = "Sample", ylab = "Log posterior probability density")

print(samples$density)

start_iteration <- 5000
end_iteration <- 10000

zoomed_density <- samples$density[start_iteration:end_iteration, ]

zoomed_iterations <- start_iteration:end_iteration

matplot(zoomed_iterations, zoomed_density, type = "l", lty = 1,
        xlab = "Iteração", ylab = "Densidade de Probabilidade Log-Posterior",
        main = paste("Zoom na Densidade da Log-Posterior (Iterações", start_iteration, "a", end_iteration, ")"))

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

#######################################################################

samples_df <- posterior::as_draws_df(samples)
samples_df
posterior::summarise_draws(samples_df)


unique(samples_df$.chain)
bayesplot::mcmc_pairs(samples_df)

############## Selecting specific chains after burn-in #######################

library(posterior)
library(dplyr)

# Converta para data.frame se ainda não tiver feito
samples_df <- posterior::as_draws_df(samples)

# Defina as cadeias e o número de iterações para burn-in
selected_chains <- c(5)
burn_in <- 500

# Filtra para cadeias selecionadas e iterações após o burn-in
samples_filtered <- samples_df %>%
  filter(.chain %in% selected_chains, .iteration > burn_in)

# Calcular as estatísticas apenas nas cadeias selecionadas pós burn-in
summary_selected <- posterior::summarise_draws(samples_filtered)
summary_selected




###########################################################################################
library(ggplot2)

n_iter <- dim(samples$pars)[2]  # Number of iterations
n_chain <- dim(samples$pars)[3]  # Number of chains

# data.frame
df_samples <- data.frame(
  iter = rep(1:n_iter, times = n_chain),  
  beta = c(samples$pars["beta", , ]),  
  gamma = c(samples$pars["gamma", , ]), 
  I0 = c(samples$pars["I0", , ]),  
  eta = c(samples$pars["eta", , ]),  
  #rho = c(samples$pars["rho", , ]),  
  chain = rep(1:n_chain, each = n_iter)  # Chain index
)


cor(df_samples$beta, df_samples$gamma)
cor(df_samples$beta, df_samples$I0)
cor(df_samples$gamma, df_samples$I0)


by(df_samples[, c("beta", "gamma")], df_samples$chain,
   function(x) cor(x$beta, x$gamma))

ggplot(df_samples, aes(x = beta, y = gamma, color = factor(chain))) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(color = "Cadeia")

ggplot(df_samples, aes(x = beta, y = gamma, color = factor(chain))) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(color = "Cadeia") +
  facet_wrap(~ chain, scales = "free")  # Cria um gráfico para cada cadeia

########################################
# Trace plot for beta

ggplot(df_samples, aes(x = iter, y = beta, color = factor(chain))) +
  geom_line() +
  labs(title = "Evolution of the beta parameter by chain", color = "Chain") +
  labs(x = "Iteration", y = "Beta", color = "Chain") +
  theme_minimal() +
  scale_y_continuous(limits = c(0.0, 10))  

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


ggplot(df_samples, aes(x = iter, y = beta, color = factor(chain))) +
  geom_line(size = 1.2) +  # linha mais grossa
  labs(
    title = "Evolution of the beta parameter by chain",
    x = "Iteration",
    y = "Beta",
    color = "Chain"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    legend.text = element_text(size = 13)
  )

ggplot(df_samples, aes(x = iter, y = gamma, color = factor(chain))) +
  geom_line(size = 1.5) +
  labs(
    title = "Evolution of the gamma parameter by chain",
    x = "Iteration",
    y = "Gamma",
    color = "Chain"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = c(0.95, 0.55),
    legend.justification = c("right", "bottom"),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 14)
  )


# Trace plot for gamma
ggplot(df_samples, aes(x = iter, y = gamma, color = factor(chain))) +
  geom_line() +
  labs(title = "Evolution of the gamma parameter by chain", color = "Chain") +
  labs(x = "Iteration", y = "Gamma", color = "Chain") +
  theme_minimal() +
  scale_y_continuous(limits = c(0.0, 10))  

# Trace plot for I0
ggplot(df_samples, aes(x = iter, y = I0, color = factor(chain))) +
  geom_line() +
  labs(title = "Evolution of the I0 parameter by chain", color = "Chain") +
  labs(x = "Iteration", y = "I0", color = "Chain") +
  theme_minimal() +
  scale_y_continuous(limits = c(0.0, 10)) 

#####ZOOM NO TRAÇO#######

# Defina seus limites ymin e ymax desejados
my_ymin <- 0.0
my_ymax <- 0.30

ggplot(df_samples_zoomed, aes(x = iter, y = gamma, color = factor(chain))) +
  geom_line() +
  labs(title = paste("Evolução do parâmetro Gamma por cadeia (Iterações", start_iteration, "a", end_iteration, ")"),
       x = "Iteração",
       y = "Gamma",
       color = "Cadeia") +
  theme_minimal() +
  ylim(my_ymin, my_ymax) # Adicione esta linha para definir os limites do eixo Y


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


################DENSIDADE####################

library(ggplot2)

# Selecione apenas os parâmetros de interesse
samples_subset <- samples_df[, c("beta", "gamma")]

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

################# PLOTANDO SEPARADO POR CADEIAS ##############

library(ggplot2)
library(tidyr)
library(dplyr)

# Mantenha a coluna .chain
samples_subset <- samples_df[, c(".chain", "beta", "gamma")]

# Formato longo
samples_long <- pivot_longer(samples_subset, cols = c("beta", "gamma"),
                             names_to = "parameter", values_to = "value")

# Densidade por cadeia
ggplot(samples_long, aes(x = value, color = factor(.chain), fill = factor(.chain))) +
  geom_density(alpha = 0.4) +
  facet_wrap(~parameter, scales = "free") +
  labs(title = "Distribuição dos parâmetros por cadeia",
       x = "Valor", y = "Densidade", color = "Cadeia", fill = "Cadeia") +
  theme_minimal()

###############  BETA ESTIMATED USING SAZONALITY ##################

library(ggplot2)
library(dplyr)
library(tibble)


# 2. Mean Parameters for each chain (after burn-in)
chains <- tibble::tribble(
  ~chain, ~temp_scale,  ~gamma, 
  "Chain 1", 0.0013,    0.209,
  "Chain 2",  0.00541 , 3.07 ,
  "Chain 3",0.00539 ,  4.49, 
  "Chain 4", 0.00546 ,  1.95, 
  "Chain 5", 0.00465, 1.61
)


# 3. Função de Brière
briere <- function(temp,scale) {
  scale * temp * (temp - 20.0) * sqrt(pmax(30.0 - temp, 0))
}


# 4. Calcular β(t) e R0(t) para cada cadeia usando ross-macdonald
df_plot <- chains %>%
  rowwise() %>%
  do({
    chain <- .$chain
    temp_scale <- .$temp_scale
    #sd_scale <- .$sd_temp_scale
    gamma <- .$gamma
    #sd_gamma <- .$sd_gamma
    
    beta <- briere(temperature_interpolate,temp_scale)
    #beta_sd <- ross(temp_values, temp_scale + sd_scale) - beta
    
    R0 <- beta / gamma
    # R0_sd <- R0 * sqrt((beta_sd / beta)^2 + (sd_gamma / gamma)^2)
    
    tibble(
      chain = chain,
      time = time_interpolate,
      beta = beta,
      #  beta_sd = beta_sd,
      R0 = R0,
      #  R0_sd = R0_sd
    )
  }) %>%
  bind_rows()

# 5. Plot: β(t)
ggplot(df_plot, aes(x = time, y = beta, color = chain)) +
  geom_point() +
  #geom_errorbar(aes(ymin = beta - beta_sd, ymax = beta + beta_sd), width = 0.3) +
  labs(x = "Tempo (meses)", y = expression(beta(t))) +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank())+
  scale_y_continuous(limits = c(0.0, 5))

# 6. Plot: R₀(t)
ggplot(df_plot, aes(x = time, y = R0, color = chain)) +
  geom_point() +
  #geom_errorbar(aes(ymin = R0 - R0_sd, ymax = R0 + R0_sd), width = 0.3) +
  labs(x = "Tempo (meses)", y = expression(R[0](t))) +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank())+
  scale_y_continuous(limits = c(0.0, 5.0))

################## TRACES AND DENSITIES ##########################

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

# Trace plot for temp_scale

ggplot(df_samples, aes(x = iter, y = temp_scale, color = factor(chain))) +
  geom_line() +
  labs(title = "Evolution of the temp_scale parameter by chain", color = "Chain") +
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

# Trace plot for I0
ggplot(df_samples, aes(x = iter, y = I0, color = factor(chain))) +
  geom_line() +
  labs(title = "Evolution of the I0 parameter by chain", color = "Chain") +
  labs(x = "Iteration", y = "I0", color = "Chain") +
  theme_minimal() +
  scale_y_continuous(limits = c(0.0, 10)) 


# Trace plot para eta
ggplot(df_samples, aes(x = iter, y = eta, color = factor(chain))) +
  geom_line() +
  labs(title = "Evolution of the eta parameter by chain", color = "Chain") +
  labs(x = "Iteration", y = "Eta", color = "Chain") +
  theme_minimal()


library(ggplot2)

# Selecione apenas os parâmetros de interesse
samples_subset <- samples_df[, c("temp_scale")]

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


########################## SIMULATION ##############################

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
    I0 = as.numeric(samples$pars["I0", n_iterations, chain]),
    beta = as.numeric(samples$pars["beta", n_iterations, chain]), 
    gamma = as.numeric(samples$pars["gamma", n_iterations, chain]),
    eta = as.numeric(samples$pars["eta", n_iterations, chain])
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


