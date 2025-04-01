
#=================================================================================@
#'
#'@export
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

#=================================================================================@
#'@export
GibbsDensity <- function(DataMeasures, A, Ai, index, a, b) {
    ## Using the function create_MeasuresDataFrame
  Measures <- DataMeasures$Measures
  D = as.numeric(Measures$D)
  Di = D[index]
  sD = as.numeric(Measures$sD)
  sDi = sD[index]
  ddot = as.numeric(Measures$ddot)
  ddoti = ddot[index]
  sddot = as.numeric(Measures$sddot)
  sddoti = sddot[index]
  sdc = as.numeric(Measures$sddot_shared)
  sdci = sdc[index]
  alpha = as.numeric(Measures$symmetric_error)


  ## getting the modidied
  A[index] = Ai

  detCov <- function(A) {
    1 + sum( (A*alpha*sdc)**2 / ((A*sddot)**2 + sD**2) )
  }


  invdiagvar = ( (Ai*sddoti)**2 + sDi**2 )**(-1)
  centeri = (Di -Ai*ddoti)

  firstExp = exp( (-1/2) * centeri**2 *invdiagvar )
  sumexp = sum( Ai*centeri * invdiagvar* alpha *( (A* (D-A*ddot) * sdc**2)  / ((A*sddot)**2 + sD**2) ) )
  secondExp = exp( (1/2) * detCov(A)**(-1) * sumexp )

  # detCov(A)*invdiagvar *(-.5)  * firstExp * secondExp * (Ai<= b)* (Ai>=a)
  return( detCov(A)**(-1/2) * invdiagvar *firstExp * secondExp * (Ai<= b)* (Ai>=a)/Ai )
}

#=================================================================================@
#'@export
proposal_sd <- function(DataMeasures) {
  Measures = DataMeasures$Measures
  Ahat = as.numeric(Measures$D)/as.numeric(Measures$ddot)
  sdAhat = Ahat *sqrt((as.numeric(Measures$sD / Measures$D))**2 +
                       as.numeric(Measures$sddot/Measures$ddot)**2)
  return(sdAhat)

}
#=================================================================================@
#'@export
logitT <- function(u, bounds) {
  lowerbound = bounds[1]
  upperbound = bounds[2]
  return(
    (exp(u) * upperbound + lowerbound) / (1+ exp(u))
  )
}


#=================================================================================@
#' @export
arctanT <- function(u, bounds) {
  lowerbound = bounds[1]
  upperbound = bounds[2]
  return(
      (upperbound - lowerbound) * tan(u) /pi + (upperbound + lowerbound) / 2
    )
}

#=================================================================================@
#'@export
initialize_SC <- function(Sc, LowerPeriod, UpperPeriod, plotGraph = F, order= T) {
  Sc = Sc[-1, ] #del first line
  n = nrow(Sc)
  rownames(Sc) = colnames(Sc) = paste0("A", 1:n)

  network = igraph::graph_from_adjacency_matrix(Sc)
  if (plotGraph) {
    plot(network)
  }

  if (order) {
    e = rexp((n+1))
    u = cumsum(e[1:n]) / sum(e)
    A = (UpperPeriod-LowerPeriod) * u + LowerPeriod
  }
  else {
    A = rep(NA, n)
  ##Source nodes
  indegree = igraph::V(network)[igraph::degree(network, mode = "in") == 0] #in-degree
  outdegree = igraph::V(network)[igraph::degree(network, mode = "out") == 0] #out-degree

  #shift
  shift = runif(length(indegree), max = (UpperPeriod-LowerPeriod))
  A[indegree] = runif(length(indegree), min = LowerPeriod, max = (UpperPeriod - shift))
  A[outdegree] = A[indegree] + shift

  intermed = igraph::V(network)[igraph::degree(network, mode = "out") > 0 & igraph::degree(network, mode = "in") >0 ]
  print(A)
  topo_order = igraph::topo_sort(network, mode = "out")
  in_topo_order = topo_order[topo_order %in% intermed]
  for (i in in_topo_order) {
    max_bounds = A[igraph::neighbors(network, i, mode = "out")]
    min_bounds = A[igraph::neighbors(network, i, mode = "in")]
    #MAJ max_bound
    if ( any(!is.na(max_bounds))) {
      max_bound = min(max_bounds, na.rm = T)
    }

    #MAJ min_bound
    if (any(!is.na(min_bounds ))) {
      min_bound =  max(min_bounds, na.rm =T)
    }

    print(c(min_bound, max_bound))
    ##simulation
    if (min_bound < max_bound) {
      A[i] = runif(1, min_bound, max_bound)
      print("yes")
    }

    else { A[i] = min_bound}
  }}
  return(A)
}


#=================================================================================@
#'@export
GibbsSampler <- function(DataMeasures, nchain,niter, burnin, Sc,
                      LowerPeriod, UpperPeriod, Transformation = "arctan",
                      lag = 10, plotGraph = T, plotChain = T, plotACF = T,
                      roundingOfValue = 3, proposal_magnitude = 1,
                      ...) {

  #---------------- Pre-seting -----------------#@
  n_ages = DataMeasures$Measures$Nb_sample
  chains = coda::mcmc.list() #list of chains
  As = coda::mcmc.list()
  acceptances = list()
  ### Bounds Index according to SC
  IndexBounds = sapply(1:n_ages, findbound, Sc) + 1 # 2 x n_ages
  IndexBounds[which(IndexBounds == 0, arr.ind = T)] = n_ages + 2

  ## standard deviation for kernel proposal


  for (c in 1:nchain) {
    #initialize
    if (c > 1) {
      plotGraph = F
    }
    init = initialize_SC(Sc, LowerPeriod, UpperPeriod, plotGraph)

    bounds = apply(IndexBounds, 2, function(b) c(LowerPeriod, init, UpperPeriod)[b]) # 2 x n_ages
    #creating array for iteration
    chain <- A <- matrix(NA, nrow = niter+1, ncol = n_ages)
    colnames(A) <- DataMeasures$Measures$SampleNames

    A[1, ] = init
    acceptance = rep(0, n_ages)
    names(acceptance) <- paste("rate", DataMeasures$Measures$SampleNames, sep = "-")


    if (Transformation == "logit") { #Test Sigmoid with parameter lambda
    ##logit transformation for init
    chain[1, ] = log((init-bounds[1,]) / (bounds[2,] - init)) # n_ages
    }

    else if (Transformation == "arctan") {
      chain[1, ] = tan( (pi/(bounds[2, ] - bounds[1, ])) * (init - bounds[1,]) - pi/2 )
    }

      X = A[1, ] #vector of proposition
      for (iter in 1:niter) {
          for (i in sample(1:n_ages)) {
      bounds_i = c(LowerPeriod, X, UpperPeriod)[IndexBounds[, i]] # (bounds for Ai whithin Gibbs)

      # MAJ Ai
      sd = proposal_magnitude * proposal_sd(DataMeasures)
      proposal = chain[iter, i] + rnorm(1, sd = sd[i])

      #proposal depends on the transformation
      if (Transformation == "logit") {
        Aproposal = logitT(proposal, bounds_i)
        top = GibbsDensity(DataMeasures, X, Aproposal, i, bounds_i[1], bounds_i[2]) *
          (bounds_i[2]-bounds_i[1]) * exp(proposal) / ( (1 + exp(proposal))**2)

        bottom = GibbsDensity(DataMeasures, X, A[iter, i], i, bounds_i[1], bounds_i[2]) *
          ((bounds_i[2]-bounds_i[1]) * exp(chain[iter, i]) / (( 1+exp(chain[iter, i]))**2 ) )
      }

      else if (Transformation == "arctan") {
      Aproposal = arctanT(proposal, bounds_i)

      top = GibbsDensity(DataMeasures, X, Aproposal, i, bounds_i[1], bounds_i[2]) *
        (bounds_i[2]-bounds_i[1]) / (pi* (1 + proposal**2))
      bottom = GibbsDensity(DataMeasures, X, A[iter, i], i, bounds_i[1], bounds_i[2]) *
        ((bounds_i[2]-bounds_i[1]) / (pi *(1+chain[iter, i]**2)) )

      }##### END OF ARCTAN TRANSFORMATION
      # print(bounds_i)
      # print(paste("X:", X))
      ratio = top /bottom

      u = runif(1)

      if (u < ratio) {
        chain[(iter+1), i] = proposal
        X[i] = Aproposal
        acceptance[i] = acceptance[i] + 1
      }

      else {
        chain[(iter+1), i] = chain[iter, i]
        }


    }### END OF COORDINATE LOOP

      #MAJ A for (iter+1)
      A[(iter + 1), ] = X

    ## MESSAGE FOR EACH THOUSAND ITER
    if (iter%% 5000 == 0) {
      cat(sprintf("\r -------- Chain %d & Iteration: %d -------", c,  iter))
    }
  }#### END OF ITER LOOP

    seq_lag = seq(burnin, (niter+1), lag)
    chain = chain[seq_lag, ]
    A = A[seq_lag, ]


    #ADD chains, As, acceptances
    chains[[paste0("chain", c)]] <- coda::mcmc(chain)
    As[[paste0("chain", c)]] <- coda::mcmc(A)
    acceptances[[paste0("acceptance", c)]] <- acceptance / niter

  }

  ##plot
  if (plotChain) {
    try(plot(As))
  }
  #---Gelman and Rubin test of convergence of the MCMC ####
  CV <- gelman.diag(As, multivariate = FALSE)
  cat(paste(
    "\n\n>> Results of the Gelman and Rubin criterion of convergence <<\n"
  ))
  for (i in 1:DataMeasures$Measures$Nb_sample) {
    cat("----------------------------------------------\n")
    cat(paste(" Sample name: ", DataMeasures$SampleNames[i], "\n"))
    cat("---------------------\n")
    cat(paste("\t\t", "Point estimate", "Uppers confidence interval\n"))
    cat(paste(
      paste("A_", DataMeasures$Measures$SampleNames[i], sep = ""),
      "\t",
      round(CV$psrf[i, 1], roundingOfValue),
      "\t\t",
      round(CV$psrf[i, 2], roundingOfValue),
      "\n"
    ))

  }

  cat("\n\n ------------------------------------------------------------------------------\n\n")
  cat(paste(
    "\n\n>> Results of sampling autocorrelation <<\n"
  ))
  print(coda::autocorr.diag(As))

  try(plot(coda::acfplot(As)))

  dev.off()

  cat("\n\n ------------------------------------------------------------------------------\n\n")
  message(".   *****  The following information are only valid if the MCMC chains have converged.   ****    ")
  cat("\n\n ------------------------------------------------------------------------------\n\n")

  #---print results ####
  sample = as.data.frame(runjags::combine.mcmc(As))
  rnames <- paste0("A_", DataMeasures$Measures$SampleNames)
  summaryMCMC <- summary(As)
  ##Matrix of results
  rnames <- paste0("A_", DataMeasures$Measures$SampleNames)

  R <- matrix(
    data = NA,
    ncol = 8,
    nrow = DataMeasures$Measures$Nb_sample,
    dimnames = list(rnames,
                    c(
                      "lower bound at 95%",
                      "lower bound at 68%",
                      "Bayes estimate",
                      "upper bound at 68%",
                      "upper bound at 95%",
                      "Convergencies: Point estimate",
                      "Convergencies: uppers confidence interval",
                      "Time Series SE"
                    )
    )
  )


  #Bayes Estimate and credibal iterval
  cat(
    "\n\n>> Bayes estimates of Age, Palaeodose and its dispersion for each sample and credible interval <<\n"
  )

  credible95 <-  apply(sample, 2, CredibleInterval, level = .95)[ 2:3, ]
  credible68 <- apply(sample, 2, CredibleInterval, level = .68)[2:3, ]
  estimate <- apply(sample, 2, mean)

  R[, c(1,5)] <- round(credible95, roundingOfValue)
  R[, c(2,4)] <- round(credible68, roundingOfValue)
  R[, 3] <-   round(estimate, roundingOfValue)

  R[, c(6, 7)] <- round(CV$psrf, roundingOfValue)
  R[, 8] <- round(summaryMCMC$statistics[, 4], roundingOfValue)

  print(R )
  cat("\n----------------------------------------------\n")

  #---Create return object ####
  .list_BayLum <- function(..., originator = sys.call(which = -1)[[1]]){
    ## set list
    l <- list(...)

    ## update originators
    attr(l, "class") <- "BayLum.list"
    attr(l, "originator") <- as.character(originator)

    return(l)

  }

  name_chains = paste(as.character(Transformation), "Sampling", sep ="_")
  output <- .list_BayLum(
    "Ages" = data.frame(
      SAMPLE = DataMeasures$Measures$SampleNames,
      AGE = estimate,
      HPD68.MIN = credible68[ 1,],
      HPD68.MAX = credible68[2, ],
      HPD95.MIN = credible95[1, ],
      HPD95.MAX = credible95[2, ],
      stringsAsFactors = FALSE
    ),
    "Sampling" = As,
    "UpperLowerBounds" = c(UpperPeriod, LowerPeriod),
    "StratiConstraints" = Sc,
    "CovarianceMatrix" = list(DataMeasures$Theta, DataMeasures$covD),
    "model" = "Gibbs",
    "acceptance" = acceptances,
    SummaryMCMC = summaryMCMC,
    name_chains = chains
  )

  BayLum::plot_Ages(object = output, model = "Jeffreys prior")

  return(
    output
  )
}









