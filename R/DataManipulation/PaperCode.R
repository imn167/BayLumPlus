



# Load required libraries
library(mvtnorm)
library(ggplot2)
library(viridis)
library(gridExtra)

# Function k
k <- function(u, mean) {
  m1 <- mean[1]
  m2 <- mean[2]
  return(m1^2 + (2*u - m2)^2 - (m1 - m2 + 2*u)^2/2)
}

# Function to calculate probability
probabilityW <- function(mean, var) {
  m1 <- mean[1]
  m2 <- mean[2]
  sd <- sqrt(var)
  alphaD <- 2 * (m2 - m1) / sd
  p <- pnorm(alphaD)

  result <- list(
    "Interior (x < y)" = round(p, 4),
    "Diagonal (x = y)" = round(1 - p, 4)
  )
  return(result)
}

# Calculate values (assuming postMean and postVar are defined)
# values <- probabilityW(postMean, postVar)
# print(names(values))

# Function g
g <- function(u, mean, var) {
  m1 <- mean[1]
  m2 <- mean[2]
  sd <- sqrt(var)
  alphaD <- 2 * (m2 - m1) / sd
  return(exp(-k(u, mean)/(2*var)) * (1 - pnorm(alphaD)) / (2 * sqrt(2*pi) * sd))
}

# Truncated sampling on cone for Gaussian
truncKone <- function(m_samples, mean, var) {
  n <- length(mean)
  density <- matrix(0, nrow = m_samples, ncol = n)
  k <- 1

  while (k <= m_samples) {
    Astar <- rnorm(n, mean, sqrt(var))
    if (Astar[1] < Astar[2]) {
      density[k, ] <- Astar
      k <- k + 1
    }
  }
  return(density)
}

# Projection sampling on cone
projKone <- function(m_samples, mean, var) {
  n <- length(mean)
  values <- matrix(rnorm(m_samples * n, rep(mean, each = m_samples), sqrt(var)),
                   nrow = m_samples, ncol = n)

  mask <- values[, 1] > values[, 2]
  if (sum(mask) > 0) {
    row_means <- rowMeans(values[mask, , drop = FALSE])
    values[mask, ] <- matrix(rep(row_means, n), ncol = n)
  }

  return(values)
}

# Function to plot gradient/contour
# Function to plot gradient/contour
plotgradient <- function(min_val, max_val, samples = NULL, kernel = "kernel", ...) {
  args <- list(...)

  # Create finer grid for better resolution
  x_seq <- seq(min_val, max_val, by = 0.05)
  y_seq <- seq(min_val, max_val, by = 0.05)
  grid <- expand.grid(x = x_seq, y = y_seq)

  if (kernel == "unconstrained") {
    mean_vec <- args$mean
    var_val <- args$var
    sigma_matrix <- diag(var_val, nrow = 2)
    Z <- dmvnorm(grid, mean = mean_vec, sigma = sigma_matrix)

  } else if (kernel == "truncated") {
    mean_vec <- args$mean
    var_val <- args$var
    alphaD <- (mean_vec[2] - mean_vec[1]) / sqrt(2 * var_val)
    sigma_matrix <- diag(var_val, nrow = 2)

    mask <- grid$x <= grid$y
    prob_norm <- pnorm(alphaD)

    Z <- numeric(nrow(grid))
    Z[mask] <- dmvnorm(grid[mask, ], mean = mean_vec, sigma = sigma_matrix) / prob_norm

  } else if (kernel == "projection") {
    mean_vec <- args$mean
    var_val <- args$var
    alphaD <- (mean_vec[2] - mean_vec[1]) / sqrt(2 * var_val)
    sigma_matrix <- diag(var_val, nrow = 2)

    mask <- grid$x <= grid$y
    equal <- abs(grid$x - grid$y) < 0.001  # Small tolerance for diagonal

    Z <- numeric(nrow(grid))

    # Interior points
    if (sum(mask) > 0) {
      Z[mask] <- dmvnorm(grid[mask, ], mean = mean_vec, sigma = sigma_matrix)
    }

    # Diagonal points
    if (sum(equal) > 0) {
      Z[equal] <- g(grid$x[equal], mean_vec, var_val)
    }
  }

  # Remove zeros and set minimum threshold for better visualization
  Z[Z == 0] <- NA

  # Create the plot
  df <- data.frame(x = grid$x, y = grid$y, z = Z)
  df <- df[!is.na(df$z), ]  # Remove NA values

  # Determine diagonal line parameters
  diag_intercept <- 0
  diag_slope <- 1

  p <- ggplot(df, aes(x = x, y = y, fill = z)) +
    geom_raster(interpolate = TRUE) +
    geom_abline(intercept = diag_intercept, slope = diag_slope,
                color = "black", linetype = "dashed", size = 1.2) +
    annotate("text", x = min_val + (max_val - min_val) * 0.3,
             y = max_val - (max_val - min_val) * 0.2,
             label = "A1 = A2",
             color = "black", size = 7, fontface = "bold") +
    scale_fill_gradientn(colors = c("lightgray", "lightblue", "blue", "purple", "red", "yellow", "white"),
                         name = "Density"
    ) +
    coord_equal() +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, hjust = 0.5),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) +
    labs(title = paste("Kernel:", kernel), x = "x", y = "y")

  print(p)
}

sigma = 6
tau = 3

Atrue = c(10, 11)

mu = c(10,11)
set.seed(123)
M = c(30, 15)
postmean = (sigma**2 * mu + tau**2 * M) / (sigma**2 + tau**2)
postvar = sigma**2 * tau**2 / (sigma**2 + tau**2)

diagcoord = rbind(postmean-2*postvar, postmean+2*postvar)
min(diagcoord)


# Comparing both methods
samplesRejection <- truncKone(8000, postmean, postvar)
sampleProj <- projKone(8000, postmean, postvar)
print(colMeans(samplesRejection))
print(colMeans(sampleProj))

# Kernel density estimation
kde_1 <- density(sampleProj[, 1])
kde_2 <- density(sampleProj[, 2])
kde_3 <- density(samplesRejection[, 1])
kde_4 <- density(samplesRejection[, 2])



# Plot diagonal function
u <- seq(min(diagcoord), max(diagcoord), length.out = 10000)
diagonal_vals <- g(u, postmean, postvar)
df_diag <- data.frame(u = u, diagonal = diagonal_vals)
p_diag <- ggplot(df_diag, aes(x = u, y = diagonal)) +
  geom_line() +
  labs(title = "Diagonal Function", x = "diagonal", y = "density_diagonal") +
  theme_minimal()
print(p_diag)

# Create bar plot for probabilities
values <- probabilityW(postmean, postvar)
df_bar <- data.frame(
  category = names(values),
  probability = unlist(values)
)

p_bar <- ggplot(df_bar, aes(x = category, y = probability)) +
  geom_bar(stat = "identity", fill = c("green", "blue")) +
  geom_text(aes(label = paste0(round(probability * 100, 1), "%")),
            vjust = -0.3, size = 4) +
  ylim(0, 1.1) +
  labs(title = "Probability Mass Distribution", y = "Probability") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10))
print(p_bar)

# Plot truncated kernel
plotgradient(min(diagcoord), max(diagcoord), kernel = "truncated", mean = postmean, var = postvar)

# Plot KDE comparisons
df_kde <- data.frame(
  u = u,
  A1_tilde = approx(kde_1$x, kde_1$y, xout = u)$y,
  A2_tilde = approx(kde_2$x, kde_2$y, xout = u)$y,
  A1 = approx(kde_3$x, kde_3$y, xout = u)$y,
  A2 = approx(kde_4$x, kde_4$y, xout = u)$y
)

# Handle NAs from approximation
df_kde[is.na(df_kde)] <- 0

p_kde <- ggplot(df_kde) +
  geom_line(aes(x = u, y = A1_tilde, color = "Ã1")) +
  geom_line(aes(x = u, y = A2_tilde, color = "Ã2")) +
  geom_line(aes(x = u, y = A1, color = "A1")) +
  geom_line(aes(x = u, y = A2, color = "A2")) +
  labs(title = "KDE Comparison", x = "u", y = "Density", color = "Variable") +
  theme_minimal()
print(p_kde)




library(ggridges)
library(tidyr)

# Your existing KDE dataframe
df_kde <- data.frame(
  u = u,
  A1_tilde = approx(kde_1$x, kde_1$y, xout = u)$y,
  A2_tilde = approx(kde_2$x, kde_2$y, xout = u)$y,
  A1 = approx(kde_3$x, kde_3$y, xout = u)$y,
  A2 = approx(kde_4$x, kde_4$y, xout = u)$y
)

# Transform to long format for ggplot
df_long <- df_kde %>%
  pivot_longer(cols = -u, names_to = "Parameter", values_to = "Density")

# Create the ggridges plot
p <- ggplot(df_long, aes(x = u, y = Parameter, height = Density, fill = Parameter)) +
  geom_ridgeline(alpha = 0.7, scale = 1.5) +
  scale_fill_manual(values = c(
    "A1_tilde" = "steelblue",
    "A2_tilde" = "lightsalmon",
    "A1" = "steelblue",
    "A2" = "lightsalmon"
  )) +
  labs(
    title = "Posterior Density Comparison",
    x = "Value",
    y = "Parameter",
    fill = "Parameter"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",  # Y-axis already shows parameter names
    plot.title = element_text(hjust = 0.5)
  )

print(p)

# Optional: Add MAP estimates as vertical lines
maps <- df_long %>%
  group_by(Parameter) %>%
  summarise(map_estimate = u[which.max(Density)])

p_with_map <- p +
  geom_vline(data = maps, aes(xintercept = map_estimate, color = Parameter),
             linetype = "dashed", size = 1, alpha = 0.8) +
  scale_color_manual(values = c(
    "A1_tilde" = "darkblue",
    "A2_tilde" = "blue",
    "A1" = "darkorange",
    "A2" = "red"
  )) +
  guides(color = "none")

print(p_with_map)


# Add true values as reference lines
true_values <- data.frame(
  Parameter = c("A1", "A1_tilde", "A2", "A2_tilde"),
  true_value = c(10, 10, 11, 11)
)

# Optional: Add MAP estimates as vertical lines
maps <- df_long %>%
  group_by(Parameter) %>%
  summarise(map_estimate = u[which.max(Density)])

p_with_references <- p +
  # Add true values as solid black lines
  geom_vline(data = true_values, aes(xintercept = true_value),
             linetype = "solid", color = "black", size = 1.2, alpha = 0.8) +
  # Add MAP estimates as dashed colored lines
  geom_vline(data = maps, aes(xintercept = map_estimate, color = Parameter),
             linetype = "dashed", size = 1, alpha = 0.8) +
  scale_color_manual(values = c(
    "Projected_A1" = "darkblue",
    "Projected_A2" = "blue",
    "A1" = "darkorange",
    "A2" = "red"
  )) +
  # Add annotation for true values
  annotate("text", x = 10, y = 4.5, label = "True A1",
           hjust = 0, vjust = 0, size = 3.5, color = "black") +
  annotate("text", x = 11, y = 4.5, label = "True A2",
           hjust = 0, vjust = 0, size = 3.5, color = "black") +
  guides(color = "none")

print(p_with_references)



# Function to create comprehensive 4-panel visualization
create_comprehensive_plot <- function(postMean, postVar, DiagCoord) {

  # Calculate probabilities
  values <- probabilityW(postMean, postVar)

  # Create range for plots
  min_coord <- min(DiagCoord)
  max_coord <- max(DiagCoord)

  # Panel 1: Unconstrained posterior
  x_seq <- seq(min_coord, max_coord, by = 0.05)
  y_seq <- seq(min_coord, max_coord, by = 0.05)
  grid <- expand.grid(x = x_seq, y = y_seq)

  sigma_matrix <- diag(postVar, nrow = 2)
  Z_unconstrained <- dmvnorm(grid, mean = postMean, sigma = sigma_matrix)

  df_unconstrained <- data.frame(x = grid$x, y = grid$y, z = Z_unconstrained)

  p1 <- ggplot(df_unconstrained, aes(x = x, y = y, fill = z)) +
    geom_raster(interpolate = TRUE) +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", size = 1) +
    annotate("text", x = min_coord + (max_coord - min_coord) * 0.3,
             y = max_coord - (max_coord - min_coord) * 0.2,
             label = "A1 = A2", color = "black", size = 5, fontface = "bold") +
    scale_fill_gradientn(
      colors = c("lightgray", "lightblue", "blue", "purple", "red", "yellow", "white"),
      name = "Density"
    ) +
    coord_equal() +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 10),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      legend.key.height = unit(0.8, "cm")
    ) +
    labs(x = "x", y = "y")

  # Panel 2: Projected posterior (isotonic regression)
  mask <- grid$x <= grid$y
  Z_projected <- numeric(nrow(grid))
  Z_projected[mask] <- dmvnorm(grid[mask, ], mean = postMean, sigma = sigma_matrix)
  Z_projected[Z_projected == 0] <- NA

  df_projected <- data.frame(x = grid$x, y = grid$y, z = Z_projected)
  df_projected <- df_projected[!is.na(df_projected$z), ]

  p2 <- ggplot(df_projected, aes(x = x, y = y, fill = z)) +
    geom_raster(interpolate = TRUE) +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", size = 1) +
    annotate("text", x = min_coord + (max_coord - min_coord) * 0.3,
             y = max_coord - (max_coord - min_coord) * 0.2,
             label = "A1 = A2", color = "black", size = 5, fontface = "bold") +
    scale_fill_gradientn(
      colors = c("lightgray", "lightblue", "blue", "purple", "red", "orange"),
      name = "Density",
      na.value = "lightgray"
    ) +
    coord_equal() +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 10),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      legend.key.height = unit(0.8, "cm")
    ) +
    labs(x = "x", y = "y")

  # Panel 3: Probability mass distribution
  df_bar <- data.frame(
    category = names(values),
    probability = unlist(values)
  )

  p3 <- ggplot(df_bar, aes(x = category, y = probability)) +
    geom_bar(stat = "identity", fill = c("green", "blue"), width = 0.6) +
    geom_text(aes(label = paste0(round(probability * 100, 2), "%")),
              vjust = -0.3, size = 3.5, fontface = "bold") +
    ylim(0, 1.1) +
    labs(y = "Probability", x = "") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9, face = "bold"),
      axis.text.y = element_text(size = 9),
      axis.title = element_text(size = 10),
      plot.title = element_text(size = 11, hjust = 0.5)
    )

  # Panel 4: Diagonal function
  u <- seq(min_coord, max_coord, length.out = 1000)
  diagonal_vals <- g(u, postMean, postVar)
  df_diag <- data.frame(u = u, diagonal = diagonal_vals)

  p4 <- ggplot(df_diag, aes(x = u, y = diagonal)) +
    geom_line(color = "blue", size = 1) +
    labs(x = "x", y = "Density") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 10),
      axis.text.x = element_text(size = 9, face = "bold"),
      axis.text.y = element_text(size = 9, face = "bold"),
      legend.title = element_text(size = 9)
    ) +
    annotate("text", x = mean(u)-5, y = max(diagonal_vals) * 0.9,
             label = "Diagonal", color = "blue", size = 3.5)

  # Arrange all plots
  library(gridExtra)

  # Create the combined plot
  combined_plot <- grid.arrange(
    p1, p2, p3, p4,
    ncol = 2, nrow = 2,
    heights = c(1, 1),
    widths = c(1, 1)
  )

  return(combined_plot)
}

comprehensive_plot <- create_comprehensive_plot(postmean, postvar, diagcoord)


#### SteierRom data #####
simulated_data <- function(n_ages, delta, sigma, seed,  half = 1000){
  ### simulation des données
  # set.seed(seed) #fixer la souche pour générer des données
  Y <- rnorm(n_ages, sd = sigma)
  mean <- half  + ((1:n_ages) - (n_ages+1) /2) *delta
  M <-  Y + mean
  depth <- seq(200, 400, length.out = n_ages)

  #######################@
  dd <- data.frame(Depth = depth, Age = mean, std = rep(sigma, n_ages), Measure = M,
                   Names = factor(paste0("mu[", 1:n_ages, "]"), levels = paste0("mu[", 1:n_ages, "]")))
  return(dd)

}

silent <- function(expr) {
  sink(nullfile())
  on.exit(sink())
  eval(substitute(expr), envir = parent.frame())
}

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

    plotList[[i]] = plotHpd(list(Unconstrained, Isotonic, uo, nicholls), c("unconstrained", "Iso", "uo", "nicholls")) +
      ggplot2::geom_point(data = sim, ggplot2::aes(x=Names, y= Age, color = "True Age"), inherit.aes = F) +
      ggplot2::labs(title = paste("Simulation with Delta = ", Delta[i]))

  }

  return(plotList)
}



#experiment A
outputA = SteierRomA( c(150), 100, 15)


#experiment B
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

    plotList[[i]] = plotHpd(list(Unconstrained, Isotonic$Iso, uo, nicholls), c("unconstrained", "Iso", "order", "nicholls")) +
      ggplot2::geom_point(data = sim, ggplot2::aes(x=Names, y= Age, color = "True Age"), inherit.aes = F) +
      ggplot2::labs(title = paste("Simulation with size = ", Nb_sample[i]))

  }

  return(plotList)
}

outputB <- SteierRomB(0, 100, c(15,30,50))


#### Experiment C ####

generate_sample_pair <- function(delta, sigma, center) {

  #true values
  t1_true <-  center + delta /2
  t2_true <-  center - delta /2
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



## check coverage
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



run_single_experiment <- function(delta_t, n_pairs) {

  sigma = 100
  Nb_sample = 2
  center = 1000
  results <- list()

  cat(paste("Processing delta_t =", delta_t, "- par :"))

  for (i in 1:n_pairs) {
    if (i %% 10 == 0) cat(paste(i, ""))

    #generate sample pair
    pair <- generate_sample_pair(delta_t, sigma, center)
    sim <- pair$data
    true_ages <- c(pair$t1_true, pair$t1_true)

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

    #Check coverage for each model

    uo_coverage <- check_coverage(true_ages, uo)
    nicholls_coverage <- check_coverage(true_ages, nicholls)
    isotonic_coverage <- check_coverage(true_ages, isotonic)

    results[[i]] <- data.frame(
      pair_id = i,
      delta_t = delta_t,
      sample_1_true = true_ages[1],
      sample_2_true = true_ages[2],

      # Coverage results for each model

      uo_s1_coverage = uo_coverage[1],
      uo_s2_coverage = uo_coverage[2],
      nicholls_s1_coverage = nicholls_coverage[1],
      nicholls_s2_coverage = nicholls_coverage[2],
      isotonic_s1_coverage = isotonic_coverage[1],
      isotonic_s2_coverage = isotonic_coverage[2]
    )

  }

  cat("\n")
  return(do.call(rbind, results))
  # return(isotonic)

}

cat("Running Computer Experiment C with models ...")

DELTA_T_VALUES <- c(0,  150) #c(0, 10, 25, 50, 100, 150)
N_PAIRS <- 10

# Sequential processing (MCMC models don't parallelize well)
all_results <- list()
for (dt in DELTA_T_VALUES) {
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

    uo_miss_rate = (sum(!uo_s1_coverage, na.rm = TRUE) +
                      sum(!uo_s2_coverage, na.rm = TRUE)) / n_samples * 100,

    nicholls_miss_rate = (sum(!nicholls_s1_coverage, na.rm = TRUE) +
                            sum(!nicholls_s2_coverage, na.rm = TRUE)) / n_samples * 100,

    isotonic_miss_rate = (sum(!isotonic_s1_coverage, na.rm = TRUE) +
                            sum(!isotonic_s2_coverage, na.rm = TRUE)) / n_samples * 100,

    .groups = 'drop'
  )


# Create comparison plot
plot_data <- summary_stats %>%
  dplyr::select(delta_t, ends_with("_miss_rate")) %>%
  tidyr::pivot_longer(cols = ends_with("_miss_rate"),
                      names_to = "model",
                      values_to = "miss_rate") %>%
  dplyr::mutate(model = gsub("_miss_rate", "", model))

p <- ggplot(plot_data, aes(x = delta_t, y = miss_rate, color = model, shape = model)) +
  # Theoretical 5% line
  geom_hline(yintercept = 5, linetype = "dashed", color = "black", linewidth = 1) +

  # Model results
  geom_point(size = 3) +
  geom_line(alpha = 0.7) +

  # Formatting
  scale_x_continuous(name = "True Age Difference Δt (years)") +
  scale_y_continuous(name = "Percentage of 95% CI Missing True Age",
                     limits = c(0, max(plot_data$miss_rate, 15))) +
  scale_color_manual(values = c("unconstrained" = "blue",
                                "uo" = "red",
                                "nicholls" = "green",
                                "isotonic" = "orange")) +
  scale_shape_manual(values = c("unconstrained" = 16,
                                "uo" = 17,
                                "nicholls" = 18,
                                "isotonic" = 15)) +

  ggtitle("Computer Experiment C - Your Models",
          subtitle = "Comparison of coverage rates for different Bayesian approaches") +

  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  ) +

  # Add annotation about the 5% line
  annotate("text", x = max(DELTA_T_VALUES) * 0.7, y = 5.5,
           label = "Expected 5%", color = "black", size = 3.5)

print(p)


# Print results
cat("\n=== EXPERIMENT RESULTS ===\n")
print(summary_stats)

cat("\n=== KEY FINDINGS ===\n")
for (model in c( "uo", "nicholls", "isotonic")) {
  miss_rates <- summary_stats[[paste0(model, "_miss_rate")]]
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























































