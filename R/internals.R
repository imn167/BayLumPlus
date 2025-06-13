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


####
#'@description general Theme
#'@author Bouafia Imène (LMJL)
#'@md
#'@noRd

BayLumTheme <- function() {

  theme <- ggplot2::theme_minimal() + ggplot2::theme(axis.text = ggplot2::element_text(face = "bold", color = "black", size = 12),
                                                     axis.title =  ggplot2::element_text(face = "bold", color = "#342F2E", size = 12),
                                                     legend.text = ggplot2::element_text(size = 8, face = "bold", color = "black"),
                                                     legend.title = ggplot2::element_text(size = 10, face = "bold", color = "black"))

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

#=================================================================================@
AgeApprox <- function(DataMeasures) {
  D = as.numeric(DataMeasures$Measures$D)
  sD = as.numeric(DataMeasures$Measures$sD)
  ddot = as.numeric( DataMeasures$Measures$ddot)
  sddot = as.numeric( DataMeasures$Measures$sddot)
  Ahat = D/ddot
  sdAhat = Ahat * sqrt( (sD/D)**2 + (sddot/ddot)**2 )
  return(list(Ahat = Ahat, sdAhat = sdAhat))

}

#=================================================================================@

strictify_monotonic <-  function(A, min_gap= .5, jitter_strength=1e-2) {
  n = length(A)
  jitter = cumsum(runif(n, min_gap, min_gap + jitter_strength))
  A_strict = A + jitter

  return(A_strict)
}


#### Initialization over the whole study period [T1, T2]
nichollsInit <- function(I, upper, lower) {
  s = runif(1, min = 0, max = (upper-lower))
  first = runif(1, min = lower, max = (upper-s))
  e = rexp((I-1))

  return(list(s =s , e= e, first = first))
}

nichollsBRInit <- function(I, upper, lower) {
  s = runif(1, min = 0, max = (upper-lower))
  first = runif(1, min = lower, max = (upper-s))
  e = rexp((I-1))
  b = rbeta((I-2), .5,(I-2))
  z = rbinom((I-2), 1, .5)

  return(list(s =s , e= e, first = first, b =b, z =z))
}


#=================================================================================@
#### Networks Functions ####

remove_transitive_edges <- function(G) {
  reduced_G = rlang::duplicate(G)
  vertices = igraph::V(G)
  for (u in vertices) {
    message(paste("----------------- traitement sommet", u, "--------------"))
    #get all the neighbors of the vertices u
    u_neighbors = igraph::neighbors(G, u, mode = "out")
    #look for the descendants of each neighbors
    for (nei in u_neighbors) {
      message(paste("traitement voisin", nei))
      childs = igraph::subcomponent(G, nei, mode = "out")[-1]
      for (child in childs) {
        if (igraph::are_adjacent(reduced_G, u, child)) {
          print(igraph::E(reduced_G, c(u, child)))
          reduced_G = igraph::delete_edges(reduced_G, igraph::E(reduced_G, P = c(u,child)))
        }
      }
    }
  }
  return(reduced_G)
}


## Network Visualization

network_vizualization <- function(network, vertices_labels, interactive = FALSE, ...) {

  ##layout for stratigraphic constraints
  layout = igraph::layout_with_sugiyama(network)$layout
  n = length(igraph::V(network))
  if (interactive) {
    nodes <- data.frame(id = 1:n,
      label = vertices_labels,
      size = 25,
      shape = "dot",
      font = list(size = 20, face = "bold"),
      stringsAsFactors = FALSE,
      x = layout[, 1] * 100,
      y = -layout[, 2] * 100   # invert Y for visNetwork
    )

    edges = igraph::as_data_frame(network, what = "edges")
    names(edges)[1:2] <- c("from", "to")

    visNetwork::visNetwork(nodes, edges) %>%
      visNetwork::visEdges(arrows = "to") %>%
      visNetwork::visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visNetwork::visLayout(randomSeed = 123)

  }


  else {
  plot(
    network,
    layout = layout,
    vertex.label = vertices_labels,
    vertex.size = 20,
    vertex.color = adjustcolor("lightblue", alpha.f = 0.6),
    edge.arrow.size = 0.4,  # Smaller arrowheads
    edge.width = 2,
    asp = 0,
    edge.curved = 0.1
  )}
}




#=================================================================================@

#### Isotonic Regression ####

IsotonicRegDAG = function(network, Ahat, weights) {

  n = length(Ahat)
  # cp variables
  A = CVXR::Variable(n, name = "A")

  #graph
  m = length(igraph::E(network))
  ##quadratic expression
  objectif = CVXR::Minimize(CVXR::sum_squares( CVXR::multiply(weights, (A-Ahat)) ))
  #optimization matrix
  M = matrix(0, nrow = m, ncol = n)
  edges_list = igraph::as_edgelist(network, names = F)
  for (i in 1:m) {
    M[i, edges_list[i, ]] = c(1,-1)
  }

  constraints = M %*% A <= 0
  problem = CVXR::Problem(objectif, list(constraints))
  results = CVXR::solve(problem)
  IsoA = results$getValue(A)
  return(list(A=IsoA, solver = results))
}



IsotonicCurve <- function(network, object, level = .95) {
  #get all mcmc samples
  sample = runjags::combine.mcmc(object$Sampling) ## mcmc sample
  SampleNames = object$Ages$SAMPLE

  n = length(SampleNames)

  w = 1/ as.numeric(object$Summary[, 8])^2 #inv of the estimated variance

  ## apply for each age vector the Isotonic Regression
  IsoSamples = as.matrix(pbapply::pbapply(sample, 1, function(Ahat, network, weights) {
    Sys.sleep(.003)
    t(IsotonicRegDAG(network, Ahat, weights )$A)},
                     network = network , weights = w)) # n_ages x n_iter
  IsoSamples = t(IsoSamples)
  HPD = apply(IsoSamples, 2, arkhe::interval_hdr, level = level)
  print(HPD)

  IsoSummary = data.frame(lower = HPD[1, ], upper = HPD[2, ], avg = apply(IsoSamples, 2, mean),
                          Samples = factor(SampleNames, levels = SampleNames), Unit = 1:n )

  return(list(chain = IsoSamples, summary = IsoSummary))

}

PlotIsotonicCurve <- function(network, object, level = .95) {

  df <- IsotonicCurve(network, object, level)[[2]]
  n = dim(df)[1]

  df %>% ggplot2::ggplot(ggplot2::aes(x = Unit, ymin = lower, ymax = upper), fill = "orange") +
    ggplot2::geom_ribbon(alpha = .4) +
    ggplot2::geom_line(ggplot2::aes(y = lower), color = "orange", group = 1) +
    ggplot2::geom_line(ggplot2::aes(y = upper), color = "orange", group = 1) +
    ggplot2::geom_line(ggplot2::aes(y = avg), color = "orange", group = 1, size =1.5) +
    ggplot2::geom_point(ggplot2::aes(x = Unit, y = lower), color = "blue") +
    ggplot2::geom_point(ggplot2::aes(x = Unit, y = upper), color = "red") +
    ggplot2::geom_point(ggplot2::aes(x = Unit, y = avg), color = "black") +
    BayLumTheme() + ggplot2::ylab("IsotonicRegression") +
    ggplot2::scale_x_continuous(breaks = df$Unit, labels = df$Samples) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45))


}







