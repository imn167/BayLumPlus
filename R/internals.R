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
      "Credible.Interval.Inf" = round(I[i, 1], digits = roundingOfValue),
      "Credible.Interval.Sup" = round(I[i, 2], digits = roundingOfValue),
      "level" = level
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
                                                     legend.title = ggplot2::element_text(size = 10, face = "bold", color = "black")
                                                     )

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
  if (is.list(hpd_output)) { # Multi-modal distribution (at leat two interval for the HPD regions)
    HPD <- data.frame()
    ## list manip
    for (var in names(hpd_output)) {
      d = dim(hpd_output[[var]])[1]
      hpd <- matrix(hpd_output[[var]][, 1:2], nrow = d)
      tab <- data.frame(Samples = rep(var, d), inf = hpd[, 1], sup = hpd[, 2])
      HPD <- HPD %>% dplyr::bind_rows(tab)
    }
    HPD <- HPD  %>% dplyr::mutate(Models = name) #If variable with several regions HPD --> repeated columns
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

stratify_pvalue <- function(p_value) {
  stars <- character(length(p_value))
  stars[p_value <= .001] = "***" # not converged
  stars[p_value <= .01 & p_value > .001] = "**" #not converged
  stars[p_value <= .05 & p_value > .01] = "*" #not converged
  stars[p_value < .1 & p_value > .05] = "." #closed to the limit of the test"
  stars[p_value > .1] = "" #converged as the system is stationnary (values of the Z stat between -2, and 2)
  return(stars)
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

C14Init <- function(I, xbound, Sc, unconstrained = FALSE){
  lowers = xbound[(2*(1:I)-1)]
  uppers = xbound[2*(1:I)]

  if (unconstrained) {
    return(list(Age = runif(I, lowers, uppers), invalpha = rgamma(I, shape = 3, rate = 4)))
  }

  A = rep(0, I)
  A[1] = runif(1,lowers[1], uppers[1])
  for (i in 2:I) {
    # l = max(Sc[1:i, i]*c(lowers[i], A[1:(i-1)]))
    A[i] = runif(1,lowers[i], uppers[i])
  }
  network = remove_transitive_edges(Sc)
  edges = igraph::as_edgelist(network, names = F)

  A_order = IsotoneOptimization::solve_isotone_DAG(y = A, w = rep(1, I), Emat = edges)

  list(Age =A_order, invalpha = rgamma(I, shape = 3, rate = 4))
}


#=================================================================================@
#### Networks Functions ####
#Build the network fixed by the matrix of constraints
# Title
# @param StratiConstraints [numeric matrix] or [character] : The stratigraphic relation between samples.
# @param n_samples [numeric] number of samples if a strict order constraints
#
# @returns [igraph] an igraph network that summarize the constraints
#

buildNetwork <- function(StratiConstraints, n_samples = NaN) {
  ##no Strati
  if (length(StratiConstraints) == 0) {
    StratiConstraints <- rbind(rep(1, n_samples),
                               upper.tri(matrix(rep(1), ncol = n_samples, nrow = n_samples))*1)
  }
  ## Strati
  else{
    if (is(StratiConstraints)[1] == "character") {
      SCMatrix <- read.csv(StratiConstraints, sep = sepSC)
      StratiConstraints <- as.matrix(SCMatrix)
    }
  }

  Sc = StratiConstraints[-1, ]
  network = igraph::graph_from_adjacency_matrix(Sc)

  return(network)

}

# EdgePruner / EdgeSweeper / CoreGraph



#=================================================================================@

#### Isotonic Regression ####

## Several solver are available by default ECOS / clarabel (Rust) / OSQP /SCS

IsotonicRegDAG = function(network, Ahat, weights) {

  n = length(Ahat)
  # cp variables
  A = CVXR::Variable(n, name = "A")

  #graph
  m = length(igraph::E(network))
  ##quadratic expression
  objectif = CVXR::Minimize(CVXR::sum_squares( CVXR::multiply(weights, (A-Ahat)) )) # Solver by default ECOS
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


##==================================================================================@
findbound <- function(index,  Sc) {
  Sc = Sc[-1, ] #del lower bound line

  lowerbound_index = which(Sc[, index] == 1)
  if (length(lowerbound_index) > 0) { a = max(lowerbound_index)}

  else { a = 0}

  upperbound_index = which(Sc[index, ] == 1)
  if (length(upperbound_index) > 0 ) {b = min(upperbound_index)}
  else {b = -1}
  return(c(a,b))
}


##==================================================================================@

### C14 Calibration graphs

plotC14_Cal <- function(Nb_sample, PriorAge, TableauCalib, AgePlot95, AgePlotMoy, Data_C14Cal){

couleur=rainbow(Nb_sample)
sup = max(AgePlot95[, 2])
inf = min(AgePlot95[, 1])
par(mfrow=c(1,1),las = 0,mar=c(5,5,2,2))
xl <- c(min(PriorAge[seq(1,(2*Nb_sample-1),2)]),max(PriorAge[seq(2,(2*Nb_sample),2)]))
plot(
  xl,
  xl,
  col = "white",
  xlab = "C-14 Age BP (ka)",
  ylab = "C-14 cal. BP (ka)",
  xaxt = "n",
  yaxt = "n",
  cex.lab = 1.8,
  xlim = c(inf - (sup-inf)* .75, sup + (sup-inf)* .75)
)
axis(2,cex.axis=2)
axis(1,cex.axis=2)
polygon(c(TableauCalib[, 1], rev(TableauCalib[, 1])),
        c(TableauCalib[, 2] + 2 * TableauCalib[, 3],
          rev(TableauCalib[, 2] - 2 * TableauCalib[, 3]))/1000,
        col = "gray", border = "black")
for (i in 1:Nb_sample) {
  lines(c(AgePlot95[i, 1:2]), rep(Data_C14Cal[i], 2)/1000, col = couleur[i],
        lwd = 4)
  lines(AgePlotMoy[i], Data_C14Cal[i]/1000, col = "black", lwd = .5,
        type = "p")
}

}











