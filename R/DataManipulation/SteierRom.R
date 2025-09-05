#### SteierRom data #####
simulated_data <- function(n_ages, delta, sigma, seed,  half = 10000){
  ### simulation des données
  set.seed(seed) #fixer la souche pour générer des données
  Y <- rnorm(n_ages, sd = sigma)
  mean <- half  + ((1:n_ages) - (n_ages+1) /2) *delta
  M <-  Y + mean
  depth <- seq(200, 400, length.out = n_ages)

  #######################@
  dd <- data.frame(Depth = depth, Age = mean, std = rep(sigma, n_ages), Measure = M,
                   Names = factor(paste0("A_", 1:n_ages), levels = paste0("A_", 1:n_ages)))
  return(dd)

}

library(ggplot2)

models <<- c("unconstrained", "isotonic", "uo", "nicholls")

paper_color <<-   c("unconstrained" = "blue",
                 "uo" = "red",
                 "nicholls" = "green",
                 "isotonic" = "orange",
                 "True Age" = "purple")

silent <- function(expr) {
  sink(nullfile())
  on.exit(sink())
  eval(substitute(expr), envir = parent.frame())
}


#### Experiment A ####
SteierRomA <- function( Delta, sigma, Nb_sample, seed = 123) {

  #collectiong graph comparaison
  plotList = list()
  for (i in seq_along(Delta)) {

    cli::cat_rule(paste("Starting delta =", Delta[i], "loop"))

    sim <- simulated_data(Nb_sample, Delta[i], sigma, 123)

    Unconstrained <- silent(Compute_AgeS_D(list(M = sim$Measure, sigma= rep(sigma,Nb_sample )), Nb_sample = Nb_sample,
                                           SampleNames = sim$Names, ThetaMatrix = NULL, prior = "classical",
                                           PriorAge = rep(c(1, 40000), Nb_sample),
                                           Iter = 2000, burnin = 50000, t = 10))
    uo <- silent(Compute_AgeS_D(list(M = sim$Measure, sigma= rep(sigma,Nb_sample )), Nb_sample = Nb_sample,
                                SampleNames = sim$Names, ThetaMatrix = NULL, prior = "classicalorder", PriorAge = rep(c(1, 40000),
                                                                                                                      Nb_sample),
                                Iter = 2000, burnin = 50000, t = 10))
    nicholls <- silent(Compute_AgeS_D(list(M = sim$Measure, sigma= rep(sigma,Nb_sample )), Nb_sample = Nb_sample,
                                      SampleNames = sim$Names, ThetaMatrix = NULL, prior = "classicalnicholls",
                                      PriorAge = rep(c(1, 40000), Nb_sample),
                                      Iter = 2000, burnin = 50000, t = 10))
    Isotonic = IsotonicCurve(c(), Unconstrained, interactive = F)

    plotList[[i]] = plotHpd(list(Unconstrained, Isotonic, uo, nicholls), models) +
      ggplot2::geom_point(data = sim, ggplot2::aes(x=Names, y= Age, color = "True Age"), inherit.aes = F) + scale_color_manual(values = paper_color)
      # ggplot2::labs(title = paste("Simulation with Delta = ", Delta[i])) +

  }

  return(list(graphics = plotList, isotonic = Isotonic))
}

sigma = 100
n_sample = 15
delta = c(0,10,20,50,100,150)
outputA = SteierRomA( delta, sigma, n_sample)
file ="../../Isotonic/images/Simulations/"
mapply(ggsave, filename = paste0(file, "expA_delta_", delta, ".png"), plot = outputA$graphics, MoreArgs = list(width = 8, height = 6, dpi = 320))


#### Experiment B ####
SteierRomB <- function(delta, sigma, Nb_sample, seed = 123) {

  #collectiong graph comparaison
  plotList = list()
  for (i in seq_along(Nb_sample)) {

    cli::cat_rule(paste(" Starting delta =", Nb_sample[i], "loop"))

    sim <- simulated_data(Nb_sample[i], delta, sigma, 123)

    Unconstrained <- silent(Compute_AgeS_D(list(M = sim$Measure, sigma= rep(sigma,Nb_sample[i] )), Nb_sample = Nb_sample[i],
                                           SampleNames = sim$Names, ThetaMatrix = NULL, prior = "classical", PriorAge = rep(c(1, 40000), Nb_sample[i]),
                                           Iter = 2000, burnin = 50000, t = 10))
    nicholls = silent(Compute_AgeS_D(list(M = sim$Measure, sigma= rep(sigma,Nb_sample[i] )), Nb_sample = Nb_sample[i],
                                     SampleNames = sim$Names, ThetaMatrix = NULL, prior = "classicalnicholls", PriorAge = rep(c(1, 40000), Nb_sample[i]),
                                     Iter = 2000, burnin = 50000, t = 10))
    uo <- silent(Compute_AgeS_D(list(M = sim$Measure, sigma= rep(sigma,Nb_sample[i] )), Nb_sample = Nb_sample[i],
                                SampleNames = sim$Names, ThetaMatrix = NULL, prior = "classicalorder", PriorAge = rep(c(1, 40000), Nb_sample[i]),
                                Iter = 2000, burnin = 50000, t = 10))
    Isotonic = PlotIsotonicCurve(c(), Unconstrained, interactive = F)

    plotList[[i]] = plotHpd(list(Unconstrained, Isotonic$Iso, uo, nicholls), models) +
      ggplot2::geom_point(data = sim, ggplot2::aes(x=Names, y= Age, color = "True Age"), inherit.aes = F) + ggplot2::scale_color_manual(values = paper_color)
      # ggplot2::labs(title = paste("Simulation with size = ", Nb_sample[i])) +

  }

  return(plotList)
}

sample_size = c(15,30,50)
outputB <- SteierRomB(0, sigma, sample_size)

mapply(ggsave, filename = paste0(file, "ExB_n", sample_size, ".png"), plot = outputB, MoreArgs = list(width = 9, height = 6, dpi = 320))


#### Experiment C ####

generate_sample_pair <- function(delta, sigma, center) {

  #true values
  t1_true <-  center - delta /2
  t2_true <-  center + delta /2
  #add measurement scatter
  x1 <- rnorm(1, mean = 0, sd = sigma)
  x2 <- rnorm(1, mean = 0, sd = sigma)

  # Observed measurements (with scatter)
  m1_obs <- t1_true + x1
  m2_obs <- t2_true + x2

  # Create data structure matching your format
  sim_data <- data.frame(
    Depth = c(200, 400),  # Arbitrary depths
    Age = c(t1_true, t2_true),  # True ages
    std = rep(sigma, 2),
    Measure = c(m1_obs, m2_obs),  # Observed with scatter
    Names = factor(paste0("mu[", 1:2, "]"), levels = paste0("mu[", 1:2, "]"))
  )

  return(list(
    data = sim_data,
    t1_true = t1_true,
    t2_true = t2_true,
    delta = delta
  ))

}

check_coverage <- function(true_ages, object) {

  coverage <- logical(length(true_ages))
  d = dim(object$Ages)[2]
  ci_data <- object$Ages[, c(d, (d-1))]

  for (i in seq_along(true_ages)) {
    lower <- ci_data[i,2]
    upper <- ci_data[i,1]
    coverage[i] <- true_ages[i] >= lower & true_ages[i] <= upper
  }

  return(coverage)


}


check_bias <- function(true_ages, object) {

  mean <- object$Ages$AGE
  bias <- abs(mean - true_ages)

  return(bias)

}

get_variance <- function(true_ages, object) {

  variance = object$Ages$SD**2
  return(variance)
}

run_single_experiment <- function(delta_t, n_pairs) {

  sigma = 100
  Nb_sample = 2
  center = 10000
  results <- list()
  results_bias <- list()

  # cat(paste("Processing delta_t =", delta_t, "- par :"))

  for (i in 1:n_pairs) {
    if (i %% 10 == 0) cli::cat_rule(paste("Start with pair number ", i, ""))

    #generate sample pair
    pair <- generate_sample_pair(delta_t, sigma, center) # list of 3 dataframe data / t1_true, t2_true
    sim <- pair$data
    true_ages <- c(pair$t1_true, pair$t2_true)
    tryCatch({

      #unconstrained model
      unconstrained <- silent(Compute_AgeS_D(
        list(M = sim$Measure, sigma = rep(sigma, Nb_sample)), Nb_sample = Nb_sample, SampleNames = sim$Names,
        ThetaMatrix = NULL, prior = "classical", PriorAge = rep(c(1,40000), Nb_sample),
        Iter = 2000, burnin = 5000, t = 10
      ))

      #order
      uo <- silent(Compute_AgeS_D(
        list(M = sim$Measure, sigma = rep(sigma, Nb_sample)), Nb_sample = Nb_sample, SampleNames = sim$Names,
        ThetaMatrix = NULL, prior = "classicalorder", PriorAge = rep(c(1,40000), Nb_sample),
        Iter = 2000, burnin = 5000, t = 10
      ))

      #Nicholls model
      nicholls <- silent(Compute_AgeS_D(
        list(M = sim$Measure, sigma = rep(sigma, Nb_sample)), Nb_sample = Nb_sample, SampleNames = sim$Names,
        ThetaMatrix = NULL, prior = "classicalnicholls", PriorAge = rep(c(1,40000), Nb_sample),
        Iter = 2000, burnin = 5000, t = 10
      ))
      #isotonic projection
      isotonic <- IsotonicCurve(c(), unconstrained, interactive = F)

    })

    # print(plotHpd(list(isotonic, uo), c("isotonic", "uo")))
    #Check coverage for each model


    uo_coverage <- check_coverage(true_ages, uo)
    nicholls_coverage <- check_coverage(true_ages, nicholls)
    isotonic_coverage <- check_coverage(true_ages, isotonic)

    uo_bias <-  check_bias(true_ages, uo)
    nicholls_bias <-  check_bias(true_ages, nicholls)
    isotonic_bias <-  check_bias(true_ages, isotonic)

    uo_var <-  get_variance(true_ages, uo)
    nicholls_var <-  get_variance(true_ages, nicholls)
    isotonic_var <-  get_variance(true_ages, isotonic)

    results[[i]] <- data.frame(
      pair_id = i,
      delta_t = delta_t,
      sample_1_true = true_ages[1],
      sample_2_true = true_ages[2],

      # Coverage results for each model

      uo_s1_coverage = uo_coverage[1],
      uo_s2_coverage = uo_coverage[2],
      u0_s1_bias = uo_bias[1],
      u0_s2_bias = uo_bias[2],
      u0_s1_var = uo_var[1],
      u0_s2_var = uo_var[2],


      nicholls_s1_coverage = nicholls_coverage[1],
      nicholls_s2_coverage = nicholls_coverage[2],
      nicholls_s1_bias = nicholls_bias[1],
      nicholls_s2_bias = nicholls_bias[2],
      nicholls_s1_var = nicholls_var[1],
      nicholls_s2_var = nicholls_var[2],


      isotonic_s1_coverage = isotonic_coverage[1],
      isotonic_s2_coverage = isotonic_coverage[2],
      isotonic_s1_bias = isotonic_bias[1],
      isotonic_s2_bias = isotonic_bias[2],
      isotonic_s1_var = isotonic_var[1],
      isotonic_s2_var = isotonic_var[2]


    )

  }

  cat("\n")
  return(do.call(rbind, results))
  # return(isotonic)

}

cat("Running Computer Experiment C with models ...")

DELTA_T_VALUES <- c(0, 10, 50, 100, 150)
N_PAIRS <- 500

# Sequential processing (MCMC models don't parallelize well)
all_results <- list()
for (dt in DELTA_T_VALUES) {
  cli::cat_rule(paste("delta = ", dt, ""))
  all_results[[paste0("dt_", dt)]] <- run_single_experiment(dt,N_PAIRS)
}



# Combine results
experiment_results <- do.call(rbind, all_results)

# Calculate summary statistics for each delta_t and model
summary_stats <- experiment_results %>%
  dplyr::filter(!is.na(uo_s1_coverage)) %>%  # Remove failed runs
  dplyr::group_by(delta_t) %>%
  dplyr::summarise(
    n_pairs = dplyr::n(),
    n_samples = 2 * n_pairs,

    # Miss rates for each model (percentage)

    uo_miss_rate_1 = sum(!uo_s1_coverage, na.rm = TRUE) / n_pairs *100,
    uo_miss_rate_2 =  sum(!uo_s2_coverage, na.rm = TRUE) / n_pairs * 100,

    nicholls_miss_rate_1 = sum(!nicholls_s1_coverage, na.rm = TRUE) / n_pairs *100,
    nicholls_miss_rate_2 = sum(!nicholls_s2_coverage, na.rm = TRUE) / n_pairs * 100,

    isotonic_miss_rate_1 = sum(!isotonic_s1_coverage, na.rm = TRUE) / n_pairs *100,
    isotonic_miss_rate_2 = sum(!isotonic_s2_coverage, na.rm = TRUE) / n_pairs * 100,


    ## bias for each model
    uo_bias_1 = mean(u0_s1_bias),
    uo_bias_2 = mean(u0_s2_bias),

    nicholls_bias_1 = mean(nicholls_s1_bias),
    nicholls_bias_2 = mean(nicholls_s2_bias),

    isotonic_bias_1 = mean(isotonic_s1_bias),
    isotonic_bias_2 = mean(isotonic_s2_bias),

    ### variance
    uo_variance_1 = mean(u0_s1_var),
    uo_variance_2 = mean(u0_s2_var),

    nicholls_variance_1 = mean(nicholls_s1_var),
    nicholls_variance_2 = mean(nicholls_s2_var),

    isotonic_variance_1 = mean(isotonic_s1_var),
    isotonic_variance_2 = mean(isotonic_s2_var),

    ## bias vs variance
    uo_bv_1 = mean(u0_s1_bias**2 + u0_s1_var),
    uo_bv_2 = mean(u0_s2_bias**2 + u0_s2_var),

    nicholls_bv_1 = mean(nicholls_s1_bias**2 + nicholls_s1_var),
    nicholls_bv_2 = mean(nicholls_s2_bias**2 + nicholls_s2_var),

    isotonic_bv_1 = mean(isotonic_s1_bias**2 + isotonic_s1_var),
    isotonic_bv_2 = mean(isotonic_s2_bias**2 + isotonic_s2_var),



    .groups = 'drop'
  )


plot_data <- function(summarydata, pattern) {
    summarydata%>% dplyr::select(delta_t, dplyr::ends_with(pattern)) %>%
    tidyr::pivot_longer(cols = -delta_t,
                        names_to = "model",
                        values_to = "value") %>% dplyr::mutate(model = gsub(pattern, "",  model) )

}


plot_experimentC <- function(summarydata, pattern) {

    dataframe <- plot_data(summarydata, pattern)

    plot_type = gsub("^_(.*)_.*", "\\1", pattern)
    subtitle = gsub("^_(.*)_(.*)", "A\\2 \\1", pattern)
    p <- ggplot(dataframe, aes(x = delta_t, y = value, color = model, shape = model))  +

      # Model results
      geom_point(size = 3) +
      geom_line(alpha = 0.7) +

      # Formatting
      scale_x_continuous(name = "True Age Difference Δt (years)") +
      scale_y_continuous(name = paste("Average", plot_type, "values"),
                         limits = c(0, max(dataframe$value, 15))) +
      scale_color_manual(values = c("unconstrained" = "blue",
                                    "uo" = "red",
                                    "nicholls" = "green",
                                    "isotonic" = "orange")) +
      scale_shape_manual(values = c("unconstrained" = 16,
                                    "uo" = 17,
                                    "nicholls" = 18,
                                    "isotonic" = 15)) +
      paper_theme() +
      theme(
        # plot.title = element_text(size = 10, face = "bold"),
        legend.title = element_blank(),
        panel.grid.minor = element_blank()
      )

    if (plot_type == "miss_rate") {

      # Add annotation about the 5% line
        # Theoretical 5% line
     p <-  p + geom_hline(yintercept = 5, linetype = "dashed", color = "black", linewidth = 1) +
        annotate("text", x = max(DELTA_T_VALUES) * .9, y = 5.5,
               label = "Expected 5%", color = "black", size = 5) +
       scale_y_continuous(name = paste("Percentage of 95% CI Missing True Age"),
                          limits = c(0, max(dataframe$value, 15)))

    }
    else if (plot_type == "bv") {
      p <-  p +
        scale_y_continuous(name = paste("Mean Squared Error (MSE)"),
                           limits = c(0, max(dataframe$value, 15))) }

    return(p)

}

all_pattern = paste0("_", c("miss_rate", "bias", "variance", "bv"), "_", c(rep(1, 4), rep(2, 4)))
plots = list()
for (pattern in all_pattern) {
 plots[[pattern]] <- plot_experimentC(summary_stats, pattern)
}
library(patchwork)
(plots$`_variance_1` + plots$`_bias_1`) / (plots$`_miss_rate_1` + plots$`_bv_1`)
ggsave(paste0(file, "ExC_A1.png"), width = 14, height = 8, dpi = 320)
(plots$`_variance_2` + plots$`_bias_2`) / (plots$`_miss_rate_2`+ plots$`_bv_2`)
ggsave(paste0(file, "ExC_A2.png"), width = 14, height = 8, dpi = 320)


# Print results
cat("\n=== EXPERIMENT RESULTS ===\n")
print(summary_stats)

cat("\n=== KEY FINDINGS ===\n")
for (model in paste0( c( "uo", "nicholls", "isotonic"))) {
  miss_rates <- summary_stats[[paste0(model, "_miss_rate_1")]]
  max_miss <- max(miss_rates, na.rm = TRUE)
  max_delta <- summary_stats$delta_t[which.max(miss_rates)]

  cat(sprintf("%s model: Max miss rate %.1f%% at Δt=%d years\n",
              toupper(model), max_miss, max_delta))
}

cat("\n=== INTERPRETATION ===\n")
cat("Models performing well (close to 5% miss rate) are more accurate.\n")
cat("Models with high miss rates for small Δt show the 'spurious precision' problem.\n")

# Save results
save(experiment_results, summary_stats, file = "experiment_c_your_models.RData")
cat("Results saved to 'experiment_c_your_models.RData'\n")





























