# This function adjust_distributions is designed to analyze a numeric dataset by fitting several
# statistical distributions (Negative Binomial, Binomial, Poisson, and Normal) to the data. It
# visualizes the data with a histogram, provides summaries of the fitted distributions, calculates
# goodness-of-fit statistics, and generates plots for each fitted distribution. Finally, it returns a
# list of the fitted models and their goodness-of-fit statistics for further analysis.

adjust_distributions <- function(data_numeric, titulo = "Histogram of Data") {
  library(fitdistrplus)

  # Certifies that it is numeric
  data_numeric <- as.numeric(data_numeric)

  # Histogram
  hist(data_numeric, breaks = 30, probability = TRUE, 
       main = titulo, col = "lightblue")

  # Fit
  fit_nb <- fitdist(data_numeric, "nbinom")
  fit_binom <- fitdist(data_numeric, "binom", fix.arg = list(size = max(data_numeric)))
  fit_pois <- fitdist(data_numeric, "pois")
  fit_norm <- fitdist(data_numeric, "norm")

  # Summaries
  cat("\n--- Negative Binomial ---\n")
  print(summary(fit_nb))

  cat("\n--- Binomial ---\n")
  print(summary(fit_binom))

  cat("\n--- Poisson ---\n")
  print(summary(fit_pois))

  cat("\n--- Normal ---\n")
  print(summary(fit_norm))

  # Goodness of fit statistics
  gof <- gofstat(list(fit_nb, fit_binom, fit_pois, fit_norm), 
                 fitnames = c("Neg. Binomial", "Binomial", "Poisson", "Normal"))
  print(gof)

  # Plots (em 2x2)
  par(mfrow = c(2, 2))
  plot(fit_nb)
  plot(fit_binom)
  plot(fit_pois)
  plot(fit_norm)
  par(mfrow = c(1, 1))  # back to normal

  # Return the settings if you want to use them later

  return(list(neg_binom = fit_nb,
              binom = fit_binom,
              poisson = fit_pois,
              normal = fit_norm,
              gof = gof))
}
