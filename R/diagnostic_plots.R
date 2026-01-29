#' Create MCMC diagnostic plots for AgeS_Computation / Compute_Ages functions
#'
#' @param echantillon MCMC object from runjags
#' @param SampleNames Character vector of sample names
#' @return List of ggplot/plot objects that can be displayed or saved
#' @export
create_diagnostic_plots <- function(echantillon, SampleNames) {

  plots <- list()

  # Title plot with convergence info
  CV <- coda::gelman.diag(echantillon, multivariate = FALSE)
  convergence_ok <- all(CV$psrf[, 1] < 1.1)
  plot.new()
  title_text <- if(convergence_ok) {
    "MCMC Diagnostics\n(Convergence: OK)"
  } else {
    "MCMC Diagnostics\n(Convergence: CHECK WARNINGS!)"
  }
  # MCMC trace plots
  plots$trace <- tryCatch({
    function() plot_MCMC(echantillon, sample_names = SampleNames)
  }, error = function(e) {
    function() {
      plot(1:10, 1:10, type = "n", axes = FALSE, xlab = "", ylab = "", main = title_text)
      text(5.5, 5.5, paste("Error creating MCMC plot:\n", e$message),
           cex = 1.2, col = "red")
    }
  })



  # Autocorrelation plots
  plots$acf <- tryCatch({
    function() coda::acfplot(echantillon)
  }, error = function(e) {
    function() {
      plot(1:10, 1:10, type = "n", axes = FALSE, xlab = "", ylab = "", main = "Autocorrelation Plot")
      text(5.5, 5.5, paste("Error creating ACF plot:\n", e$message),
           cex = 1.2, col = "red")
    }
  })

  # Add convergence info as attribute
  attr(plots, "converged") <- convergence_ok
  attr(plots, "gelman_rubin") <- CV

  return(plots)
}

#' Display diagnostic plots
#'
#' @param plots List of plot functions from create_diagnostic_plots()
#' @export
display_diagnostic_plots <- function(plots) {
  for (plot_func in plots) {
    plot_func()
  }
}

#' Save diagnostic plots to PDF
#'
#' @param plots List of plot functions from create_diagnostic_plots()
#' @param filepath Full path to output PDF file
#' @export
save_diagnostic_plots <- function(plots, filepath) {

  # Ensure directory exists
  dir.create(dirname(filepath), showWarnings = FALSE, recursive = TRUE)

  # Safer PDF handling
  tryCatch({
    pdf(file = filepath)
    on.exit(dev.off(), add = TRUE)  # Ensures dev.off() is called even if error

    for (plot_func in plots) {
      plot_func()
    }

    message(sprintf("✓ Plots saved to: %s", filepath))
    return(invisible(TRUE))

  }, error = function(e) {
    warning(sprintf("Failed to save plots: %s", e$message))
    return(invisible(FALSE))
  })
}
