## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Title:   Internal Helper Functions for BayLum
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#'@title Create BayLum List
#'
#'@description This function works similar to [list], except the case that it sets a proper
#'class and originator
#'
#'@param ... arguments passed to [list]
#'
#'@param originator [character] (*with default*): argument to set originator manually
#'
#'@author Sebastian Kreutzer, IRAMAT-CRP2A, UMR 5060, CNRS-Université Bordeaux Montaigne (France)
#'
#'@md
#'@noRd
.list_BayLum <- function(..., originator = sys.call(which = -1)[[1]]){
  ## set list
  l <- list(...)

  ## update originators
  attr(l, "class") <- "BayLum.list"
  attr(l, "originator") <- as.character(originator)

  return(l)

}

#' Bayesian Credible Interval
#'
#' Computes the shortest credible interval of the output of the MCMC algorithm
#' for a single parameter
#' @param a_chain Numeric vector containing the output of the MCMC algorithm
#'  for the parameter.
#' @param level Probability corresponding to the level of confidence used for
#'  the credible interval, default = 0.95.
#' @param roundingOfValue Integer indicating the number of decimal places to be
#'  used, default = 0.
#' @details
#'  A \eqn{(100 * level)}\% credible interval is an interval that keeps \eqn{N * (1 - level)}
#'  elements of the sample outside the interval.
#'
#'  The \eqn{(100 * level)}\% credible interval is the shortest of all those intervals.
#' @return
#'  A named vector of values containing the confidence level and the endpoints
#'  of the shortest credible interval in calendar years (BC/AD).
#' @examples
#'   data(Events); attach(Events)
#'   CredibleInterval(Event.1)
#'   CredibleInterval(Event.12, 0.50)
#' @author A. Philippe, M.-A. Vibet
#' @noRd
CredibleInterval <- function(a_chain, level = 0.95, roundingOfValue = 0) {
  sorted_sample <- sort(a_chain) # Ordering the sample
  N <- length(a_chain)           # Calculation of the sample size
  OutSample <- N * (1 - level)   # Calculation of the number of data to be outside the interval

  # Combinasion of all credible intervals
  I = cbind(sorted_sample[1:(OutSample + 1)],
            sorted_sample[(N - OutSample):N])

  l = I[, 2] - I[, 1]   # Length of intervals
  i <- which.min(l)     # Look for the shortest interval

  # Returns the level and the endpoints rounded
  return(
    c(
      "level" = level,
      "Credible.Interval.Inf" = round(I[i, 1], digits = roundingOfValue),
      "Credible.Interval.Sup" = round(I[i, 2], digits = roundingOfValue)
    )
  )
}

#======================= GGplot Theme ======#
#'@description Setting ggplot
#'@author Bouafia Imène (LMJL)
#'@md
#'@noRd

ICAgeTheme <- function(rotation_x = F) {
  tt <-  ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text( face = 'bold', color = "#993355", size = 12),
          axis.text.y = ggplot2::element_text(face = "bold", color = "#993355", size = 12),
          axis.title.x = ggplot2::element_text(face = 'bold', size = 14, color = "black"),
          axis.title.y = ggplot2::element_text(face = 'bold', size = 14, color = "black"),

          legend.text = ggplot2::element_text(size = 12, face = "bold", color = "black"),
          legend.key.size = ggplot2::unit(1.5, "cm"))
  if (rotation_x) {
    tt <- ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text( face = 'bold', color = "#993355",
                                                          size = 12, angle = 90),
            axis.text.y = ggplot2::element_text(face = "bold", color = "#993355", size = 12),
            axis.title.x = ggplot2::element_text(face = 'bold', size = 14, color = "black"),
            axis.title.y = ggplot2::element_text(face = 'bold', size = 14, color = "black"),
            legend.text = ggplot2::element_text(size = 8, face = "bold", color = "black"),
            legend.title = ggplot2::element_text(face = "bold", size = 12, color = "black"))
  }
  return(tt)
}


#=========== HPD Regions ===========#
#'@description This function compute the HPD regions
#'@author Bouafia Imène (LMJL)
#'@md
#'@noRd


HPDRegions <- function(X, level = .95) {
  ### estimation density
  kde = stats::density(X)
  #@ values
  N = length(kde$y)
  quant <- floor(N*(1-level)) #

  #sort probabilities s
  sorted_y = sort(kde$y)
  Kq <- sorted_y[quant]

  density_X <- approx(kde$x, kde$y, xout = X)$y
  ind <-  which(density_X > Kq)
  sim_HPD = X[ind]

  return(sim_HPD)
}

pallet <- c("#FFF0AAAA", "#0000FFA0", "#00AAA0", "#D44D20", "#9DDF3E", "#3BBFDF", "#F3EC5E", "#ED5524")




Transform_HPD <- function(all_hpd) {
  max_length <- max(sapply(all_hpd, length))

  pad_vector <- function(vec, max_length) rep(vec, length.out = max_length)

  D <- lapply(all_hpd, pad_vector, max_length = max_length)
  return(sapply(D, identity))
}

hpd_method <- function(name, chain, level = .95) {
  hpd_output <- apply(chain, 2, arkhe::interval_hdr, level = level)
  if (is.list(hpd_output)) {
    HPD <- data.frame()
    ## list manip
    for (var in names(hpd_output)) {
      d = dim(hpd_output[[var]])[1]
      hpd <- matrix(hpd_output[[var]][, 1:2], nrow = d)
      tab <- data.frame(age = rep(var, d), inf = hpd[, 1], sup = hpd[, 2])
      HPD <- HPD %>% dplyr::bind_rows(tab)
    }
    HPD <- HPD  %>% dplyr::mutate(Models = name)
  }
  else {
    tab = t(hpd_output[1:2, ])
    HPD <- data.frame( Samples = rownames(tab), inf = tab[, 1], sup = tab[, 2])  %>%
      dplyr::mutate(Models = name)
  }
  return(HPD)
}



