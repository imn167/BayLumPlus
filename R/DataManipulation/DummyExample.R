##
# Load required libraries
library(mvtnorm)
library(ggplot2)
library(viridis)
library(gridExtra)

## ggplot_theme
paper_theme <- function() {

  tt <-  theme_minimal() +
    theme(axis.text.x = element_text( face = 'bold', color = "black", size = 12),
          axis.text.y = element_text(face = "bold", color = "black", size = 12),
          axis.title.x = element_text(face = 'bold', size = 14, color = "black"),
          axis.title.y = element_text(face = 'bold', size = 14, color = "black"),

          legend.text = element_text(size = 12, face = "bold", color = "black"),
          legend.title = ggplot2::element_text(size = 12, face = "bold", color = "black"),
          legend.key.size = unit(1.5, "cm"))
  return(tt)
}


# Function k
k <- function(u, mean) {
  m1 <- mean[1]
  m2 <- mean[2]
  return(m1^2 + (2*u - m2)^2 - (m1 - m2 + 2*u)^2/2)
}

# Function to calculate probability on the diagonal and in the interior
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

# Function g, density over the diagonal
g <- function(u, mean, var) {
  m1 <- mean[1]
  m2 <- mean[2]
  sd <- sqrt(var)
  alphaD <- 2 * (m2 - m1) / sd
  return(exp(-k(u, mean)/(2*var)) * (1 - pnorm(alphaD)) / (2 * sqrt(2*pi) * sd))
}

# Truncated sampling on cone for Gaussian posterior
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

# Posterior Projection over the cone
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

  } else if (kernel == "truncated") { ## TRUNCATED POSTERIOR
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
             label = "x1 = x2",
             color = "black", size = 7, fontface = "bold") +
    scale_fill_gradientn(colors = c("lightgray", "lightblue", "blue", "purple", "red", "yellow", "white"),
                         name = "Density"
    ) +
    coord_equal() +
    paper_theme() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, hjust = 0.5),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) +
    labs(title = paste( kernel, "posterior"), x = "x", y = "y")

  print(p)
}


## Test
sigma = 6 # error for measures
tau = 3 # error for the prior

Atrue = c(10, 11) #true values

mu = c(10,11)
set.seed(123)
M =   c(15.26369, 9.45099) #rnorm(2, mean = Atrue, sd = sigma)
postmean = (sigma**2 * mu + tau**2 * M) / (sigma**2 + tau**2)
postvar = sigma**2 * tau**2 / (sigma**2 + tau**2)

diagcoord = rbind(postmean-2*postvar, postmean+2*postvar)
min(diagcoord)

## plot of the unconstrained posterior
p1 = plotgradient(min(diagcoord), max(diagcoord), kernel = "unconstrained", mean = postmean, var = postvar)

p2 = plotgradient(min(diagcoord), max(diagcoord), kernel = "truncated", mean = postmean, var = postvar)

usual_case = gridExtra::grid.arrange(p1, p2, ncol = 2, nrow = 1)
ggsave(filename = paste0(file, "UnconstrainedVSTruncated.png"), width = 10, height = 9, dpi = 350, plot = usual_case)

#comprehensive plot for projection
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
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", size = .8) +
    annotate("text", x = min_coord + (max_coord - min_coord) * 0.3,
             y = max_coord - (max_coord - min_coord) * 0.2,
             label = "x = y", color = "black", size = 5, fontface = "bold") +
    scale_fill_gradientn(
      colors = c("lightgray", "lightblue", "blue", "purple", "red", "yellow", "white"),
      name = "Density"
    ) +
    coord_equal() +
    paper_theme() +
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
             label = "x = y", color = "black", size = 5, fontface = "bold") +
    scale_fill_gradientn(
      colors = c("lightgray", "lightblue", "blue", "purple", "red", "orange"),
      name = "Density",
      na.value = "lightgray"
    ) +
    coord_equal() +
    paper_theme() +
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
    paper_theme() +
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
    labs(x = "Diagonal", y = "Density") +
    paper_theme() +
    theme(
      axis.title = element_text(size = 10),
      axis.text.x = element_text(size = 9, face = "bold"),
      axis.text.y = element_text(size = 9, face = "bold"),
      legend.title = element_text(size = 9)
     ) #+
    # annotate("text", x = mean(u)-5, y = max(diagonal_vals) * 0.9,
    #          label = "Diagonal", color = "blue", size = 5)

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

u <- seq(min(diagcoord), max(diagcoord), length.out = 1000)


# Plot KDE comparisons
df_kde <- data.frame(
  u = u,
  Projected_A1 = approx(kde_1$x, kde_1$y, xout = u)$y,
  Projected_A2 = approx(kde_2$x, kde_2$y, xout = u)$y,
  A1 = approx(kde_3$x, kde_3$y, xout = u)$y,
  A2 = approx(kde_4$x, kde_4$y, xout = u)$y
)

# Handle NAs from approximation
df_kde[is.na(df_kde)] <- 0

p_kde <- ggplot(df_kde) +
  geom_line(aes(x = u, y = Projected_A1, color = "isotonic")) +
  geom_line(aes(x = u, y = Projected_A2, color = "isotonic")) +
  geom_line(aes(x = u, y = A1, color = "uo")) +
  geom_line(aes(x = u, y = A2, color = "uo")) +
  geom_vline(aes(xintercept = 10), linetype = "dashed", linewidth = 1) + geom_vline(aes(xintercept = 11), linetype = "dashed", linewidth = 1)+
  labs(title = "KDE Comparison", x = "u", y = "Density", color = "Model") +
  ggrepel::geom_label_repel(data = df_kde %>% tidyr::pivot_longer(-u, , names_to = "Parameter", values_to = "Density") %>% dplyr::group_by(Parameter) %>% dplyr::summarise(x = u[which.max(Density)], y = max(Density)),
            mapping = aes(x = x , y = y, label = Parameter, colour = c( rep("uo",2), rep("isotonic", 2)))) +
  paper_theme() + scale_color_manual(values = paper_color)
print(p_kde)
ggsave(filename = paste0(file, "Marginals.png"), width = 10, height = 8, dpi = 320)
# Transform to long format for ggplot
df_long <- df_kde %>%
  tidyr::pivot_longer(cols = -u, names_to = "Parameter", values_to = "Density")

# Create the ggridges plot
p <- ggplot(df_long, aes(x = u, y = Parameter, height = Density, fill = Parameter)) +
  ggridges::geom_ridgeline(alpha = 0.7, scale = 1.5) +
  scale_fill_manual(values = c(
    "Projected_A1" = "steelblue",
    "Projected_A2" = "steelblue",
    "A1" = "lightsalmon",
    "A2" = "lightsalmon"
  )) +
  labs(
    title = "Marginal's Posterior Density Comparison",
    x = "Value",
    y = "Parameter",
    fill = "Parameter"
  ) +
  paper_theme() +
  theme(
    legend.position = "none",  # Y-axis already shows parameter names
    plot.title = element_text(hjust = 0.5)
  )

print(p)

# Add true values as reference lines
true_values <- data.frame(
  Parameter = c("A1", "A2"),
  true_value = c(10, 11)
)
p + geom_vline(aes(xintercept = true_value), data = true_values, linetype = "dashed", alpha = .8)

maps <- df_long %>%
  dplyr::group_by(Parameter) %>%
  dplyr::summarise(map_estimate = u[which.max(Density)])
maps

p + geom_vline(aes(xintercept = true_value), data = true_values, linetype = "solid", alpha = .8) +
  geom_vline(aes(xintercept = map_estimate, colour = Parameter), data = maps, linetype = "dashed", linewidth = 1) +
  scale_color_manual(
    values = c(
      "Projected_A1" = "steelblue",
      "Projected_A2" = "steelblue",
      "A1" = "lightsalmon",
      "A2" = "lightsalmon"
    )
  ) +
  # Add annotation for true values
  annotate("text", x = 10, y = 3.5, label = "True A1",
           hjust = 1, vjust = 0, size = 4, color = "black") +
  annotate("text", x = 11, y = 4.5, label = "True A2",
           hjust = 0, vjust = 0, size = 4, color = "black") +
  guides(color = "none")




