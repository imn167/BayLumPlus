
#=================================================================================@


#' Bounds for each age according to the stratigraphic constraints
#'@description This function gives the bounds of each age according to the reduced network (better computational time)
#'it will return a list of list each list within have two vector
#' *index of the ages that preceed the selected age (`lower`)
#'*index of the ages that comes after the selected age (`upper``)
#'The index are shift by one to includ the study period for the younger and older :
#'@export
findbounds <- function(network) {
  vertices = igraph::V(network)
  n = length(vertices)

  verticesTreatement <- function(v, network) {
    neighbors_in = as.numeric(igraph::neighbors(network, v, mode = "in")) #youngers ages
    neighbors_out = as.numeric(igraph::neighbors(network, v , mode = "out")) #older ages

    if (length(neighbors_in)==0) {
      neighbors_in = 0
    }

    if (length(neighbors_out)==0) {
      neighbors_out = n+1
    }

    return(list(upper = neighbors_out+1, lower = neighbors_in+1))

  }

  all_bounds = lapply(vertices, verticesTreatement, network = network)
  all_bounds
}

#=================================================================================@
#'@export
GibbsDensity <- function(DataMeasures, A, Ai, index, a, b, ...) {
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


  ## Likelyhood part
  A[index] = Ai

  detCov <- function(A) {
    1 + sum( (A*alpha*sdc)**2 / ((A*sddot)**2 + sD**2) )
  }


  invdiagvar = ( (Ai*sddoti)**2 + sDi**2 )**(-1)
  centeri = (Di -Ai*ddoti)

  firstExp = exp( (-1/2) * centeri**2 *invdiagvar )
  sumexp = sum( Ai*centeri * invdiagvar* alpha *( (A* (D-A*ddot) * sdc**2)  / ((A*sddot)**2 + sD**2) ) )
  secondExp = exp( (1/2) * detCov(A)**(-1) * sumexp )

  likelyhood = detCov(A)**(-1/2) * invdiagvar *firstExp * secondExp

  ## Priors

  # Simple case of Jeffrey
  prior = (Ai<= b)* (Ai>=a)/Ai

  ## Uniform Order on the partial order
  # prior = (Ai<= b)* (Ai>=a)

  # ## with Nicholls & Jones
  # n = length(D)
  # upper = list(...)$upper
  # lower = list(...)$lower
  # if (index == 1 || index ==  n) {
  #   nicholls = log(A[n] / A[1])**(n-2) * log(upper * A[1] / (lower * A[n]))
  #   dens = detCov(A)**(-1/2) * invdiagvar *firstExp * secondExp * (Ai<= b)* (Ai>=a)/(Ai * nicholls)
  # }
  #
  # else {
  # }

  dens = likelyhood * prior
  return( dens )
}

#=================================================================================@
#'@export
proposal_sd <- function(DataMeasures, lowerPeriod, upperPeriod, indexbounds, lambda) {
  ##
  Measures = DataMeasures$Measures
  D = as.numeric(Measures$D)
  ddot = as.numeric(Measures$ddot)
  sD = as.numeric(Measures$sD)
  sddot = as.numeric(Measures$sddot)
  ## bounds
  Ahat = D/ ddot
  bounds = apply(indexbounds, 2, function(b) c(lowerPeriod, Ahat, upperPeriod)[b])
  upper = bounds[2, ]
  lower = bounds[1, ]
  ##
  const =( Ahat * (upper - lower)  / ( (upper - Ahat) * (Ahat - lower)))
  sdhat = abs(const)*sqrt( (sD/D)**2  + (sddot / ddot)**2 )
  #
  return(sdhat)

}
#=================================================================================@
#'@export
sigmoidT <- function(u, bounds, lambda) {
  lowerbound = bounds[1]
  upperbound = bounds[2]
  return(
    (exp(lambda * u) * upperbound + lowerbound) / (1+ exp(lambda*u))
  )
}

#'@export
logitT <- function(u, bounds, lambda) {
  lowerbound = bounds[1]
  upperbound = bounds[2]
  return(
    log((u-lowerbound) / (upperbound - u)) / lambda
  )
}

#'@export
logitJacobien <- function(u, bounds, lambda) {
  lowerbound = bounds[1]
  upperbound = bounds[2]

  return(
    (upperbound - lowerbound) * lambda * exp(lambda * u) / (1+exp(lambda * u))**2
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
#' Initialization of the Gibbs Algorithm according to the fixed stratigraphy
#' @description
#' In this function we try different way of initializing our Gibbs Sampler.
#' Several limitation can be encoutered :
#' - The constraints need to be respected
#' - The ages need to belongs to a certain period of study [T1, T2] that will be fixed arbitrarily
#' - Initialisation are not only constrained by the Matrix of constraint but also to avoid numerical zeros the interval of initialisation
#' need to be reasonable according to the density of Gibbs Sampler.
#'
#' The most efficient approach would be to considered the approximation Ahat = Dhat / ddot and then make it respect the constraints
#'
#'@export
initialize_SC <- function(Sc, LowerPeriod, UpperPeriod, type, plotGraph = F, ...) {
  Sc = Sc[-1, ] #del first line
  n = nrow(Sc)
  rownames(Sc) = colnames(Sc) = paste0("A", 1:n)

  network = igraph::graph_from_adjacency_matrix(Sc)
  reduced_network = remove_transitive_edges(network)
  if (plotGraph) {
    plot(
      reduced_network,
      layout = igraph::layout_with_sugiyama(reduced_network),
      vertex.size = 11,
      vertex.color = adjustcolor("orange", alpha.f = 0.6),
      edge.arrow.size = 0.3,  # Smaller arrowheads
      edge.width = 1,
      asp = 0,
      edge.curved = 0.1
    )
  }

  if (type == "Exponential") {
    e = rexp((n+1))
    u = cumsum(e[1:n]) / sum(e)
    A = (UpperPeriod-LowerPeriod) * u + LowerPeriod
  }

  else if (type == "isotonic") {
    D = list(...)$D
    sD = list(...)$sD

    ddot = list(...)$ddot
    sddot = list(...)$sddot

    Ahat = D / ddot
    w = 1 / ((sD/D)**2 + (sddot / ddot)**2)
    A = Iso::pava(Ahat, w = w )
    A = strictify_monotonic(A)

  }
  else if (type == "brutal") {
          A = rep(NA, n)
        ##Source nodes
        indegree = igraph::V(network)[igraph::degree(network, mode = "in") == 0] #in-degree
        outdegree = igraph::V(network)[igraph::degree(network, mode = "out") == 0] #out-degree

        #shift
        shift = runif(length(indegree), max = (UpperPeriod-LowerPeriod), min = LowerPeriod)
        u = runif(length(indegree))
        A[indegree] = exp(u * log((UpperPeriod - shift) / LowerPeriod) + log(LowerPeriod) )
        A[outdegree] = A[indegree] + shift

        intermed = igraph::V(network)[igraph::degree(network, mode = "out") > 0 & igraph::degree(network, mode = "in") >0 ]
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

            u = runif(1)
            A[i] = exp(u * log(max_bound / min_bound) + log(min_bound))

          }

          else { A[i] = min_bound}
        }
  }

  else {
    D = list(...)$D
    ddot = list(...)$ddot

    A = sort(D/ddot)

  }

  return(A)
}


#=================================================================================@
#'@export
GibbsSampler <- function(DataMeasures, nchain,niter, burnin, Sc,
                      LowerPeriod, UpperPeriod, Transformation = "arctan",
                      lag = 10, plotGraph = T, plotChain = T, plotACF = T,
                      roundingOfValue = 3, proposal_magnitude = 1, lambda, n_logitStep,
                      ...) {

  #---------------- Pre-setting -----------------#@
  n_ages = DataMeasures$Measures$Nb_sample
    D = as.numeric(DataMeasures$Measures$D)
    sD = as.numeric(DataMeasures$Measures$sD)
    ddot = as.numeric(DataMeasures$Measures$ddot)
    sddot = as.numeric(DataMeasures$Measures$sddot)


  #### Network creation and visualization
  network = igraph::graph_from_adjacency_matrix(Sc, mode = "directed")
  reduced_network = remove_transitive_edges(network)



  ### MCMC lists
  chains = coda::mcmc.list() #list of chains
  As = coda::mcmc.list()
  acceptances = list()
  length_bound = list()


  ### Stratigraphics Bounds
  IndexBounds = sapply(1:n_ages, findbound, Sc) + 1 # 2 x n_ages
  IndexBounds[which(IndexBounds == 0, arr.ind = T)] = n_ages + 2



  for (c in 1:nchain) {
    #initialize
    if (c > 1) {
      plotGraph = F
    }
    init = initialize_SC(Sc, LowerPeriod, UpperPeriod, plotGraph, D = D, ddot = ddot, sddot = sddot, sD = sD, type = list(...)$type)
  print(init)
    bounds = apply(IndexBounds, 2, function(b) c(LowerPeriod, init, UpperPeriod)[b]) # 2 x n_ages

    ## standard deviation for kernel proposal
    sd = rep(proposal_magnitude, DataMeasures$Measures$Nb_sample) * proposal_sd(DataMeasures,
                                            LowerPeriod, UpperPeriod, IndexBounds, lambda)
  print(sd)
    #creating array for iteration
    chain <- A <- lbound <- matrix(NA, nrow = niter+1, ncol = n_ages)

    lbound[1,] = bounds[2, ] - bounds[1, ]
    A[1, ] = init
    acceptance = rep(0, n_ages)
    names(acceptance) <- paste("rate", DataMeasures$Measures$SampleNames, sep = "-")


    if (Transformation == "logit") {
    ##logit transformation for init
    chain[1, ] = log((init-bounds[1,]) / (bounds[2,] - init)) / lambda # n_ages
    }

    else if (Transformation == "arctan") {
      chain[1, ] = tan( (pi/(bounds[2, ] - bounds[1, ])) * (init - bounds[1,]) - pi/2 )
    }
   # sink("R/DataManipulation/output.txt")
   print(data.frame(lower = bounds[1,], A = init, upper = bounds[ 2, ]))

      X = A[1, ] #vector of proposition

      for (iter in 1:niter) {
          for (i in 1:n_ages) {
      bounds_i = c(LowerPeriod, X, UpperPeriod)[IndexBounds[, i]] # (bounds for Ai whithin Gibbs)
      lbound[(iter+1), i] = bounds_i[2]-bounds_i[1]

      logitstep = matrix(, nrow = 1, ncol = (n_logitStep + 1))
      Ak = matrix(, nrow = 1, ncol = (n_logitStep + 1))
      #### Might be an error so re - calculating the logitT
      chain[iter, i] = logitT(X[i], bounds_i, lambda)
      logitstep[1] = chain[iter, i]
      Ak[1] = X[i]
      for (k in 1:n_logitStep) {

        proposal = logitstep[k] + rnorm(1, sd = sd[i])

      #proposal depends on the transformation
      if (Transformation == "logit") {
        Aproposal = sigmoidT(proposal, bounds_i, lambda)
        top = GibbsDensity(DataMeasures, X, Aproposal, i, bounds_i[1], bounds_i[2],
                           upper = UpperPeriod, lower = LowerPeriod) * logitJacobien(proposal, bounds_i, lambda)
        bottom = GibbsDensity(DataMeasures, X, Ak[k], i, bounds_i[1], bounds_i[2], upper = UpperPeriod, lower = LowerPeriod) *
          logitJacobien(logitstep[k], bounds_i, lambda)

      }

      else if (Transformation == "arctan") {
      Aproposal = arctanT(proposal, bounds_i)

      top = GibbsDensity(DataMeasures, X, Aproposal, i, bounds_i[1], bounds_i[2]) *
        (bounds_i[2]-bounds_i[1]) / (pi* (1 + proposal**2))
      bottom = GibbsDensity(DataMeasures, X, Ak[k], i, bounds_i[1], bounds_i[2]) *
        ((bounds_i[2]-bounds_i[1]) / (pi *(1+logitstep[k]**2)) )

      }##### END OF ARCTAN TRANSFORMATION

      ratio = top / bottom

        if (is.na(ratio)) {
      print(data.frame(index = i, theta_star = proposal, A_star = Aproposal, A = Ak[k], theta = logitstep[k],
      lower = bounds_i[1], upper = bounds_i[2], top = top, bottom = bottom, iter = iter ))

        }

        if (runif(1) <= ratio) {
          logitstep[k+1] = proposal
          Ak[k+1] = Aproposal
          acceptance[i] = acceptance[i] + 1
        }
        else {
          logitstep[k+1] = logitstep[k]
          Ak[k+1] = Ak[k]
        }


      }

      # print(bounds_i)
      # print(paste("X:", X))

      X[i] = Ak[(n_logitStep +1)]
      chain[(iter+1), i] = logitstep[(n_logitStep+1)]


      if (iter < 2000 & iter%%100 == 0) {
        if (acceptance[i] > .6) {
          sd[i] = sd[i] * 1.1
        }
        else {
          sd[i] = .9 * sd[i]
          }
      }


    }### END OF COORDINATE LOOP

      #MAJ A for (iter+1)
      A[(iter + 1), ] = X

    ## MESSAGE FOR EACH THOUSAND ITER
    if (iter%% 5000 == 0) {
      cat(sprintf("\r -------- Chain %d & Iteration: %d -------", c,  iter))
    }
  }#### END OF ITER LOOP
    # sink(file = NULL)

    colnames(A) <- DataMeasures$Measures$SampleNames
    seq_lag = seq(burnin, (niter+1), lag)
    chain = chain[seq_lag, ]
    A = A[seq_lag, ]


    #ADD chains, As, acceptances
    chains[[paste0("chain", c)]] <- coda::mcmc(chain)
    As[[paste0("chain", c)]] <- coda::mcmc(A)
    acceptances[[paste0("acceptance", c)]] <- acceptance /(n_logitStep* niter)
    length_bound[[paste0("chain", c)]] <- lbound
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

  # try(plot(coda::acfplot(As)))



  plot(1:n_ages, seq(0,1, length.out = n_ages), ylab = "acceptation rate", xlab = "samples", type = "n")
  a = sapply(acceptances, lines, col = 1:n_ages)
  a = sapply(acceptances, points, col = 1:n_ages)


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
    ncol = 9,
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
                      "Time Series SE",
                      "Bayes sd"
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
  standardError <- apply(sample, 2, sd)


  R[, c(1,5)] <- round(credible95, roundingOfValue)
  R[, c(2,4)] <- round(credible68, roundingOfValue)
  R[, 3] <-   round(estimate, roundingOfValue)
  R[, 9] <-   round(standardError, roundingOfValue)

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
    name_chains = chains,
    length_bound = length_bound,
    "Summary" = R
  )

  plot_Ages(object = output, model = "Jeffreys prior")

  return(
    output
  )
}


######@
#'@export
GibbsSamplerTrunc <- function(DataMeasures, nchain,niter, burnin, Sc,
                              LowerPeriod, UpperPeriod,
                              lag = 10, plotGraph = T, plotChain = T, plotACF = T,
                              roundingOfValue = 3, proposal_magnitude = 1,
                              ...) {

  #### Pre-setting ####

  n_ages = DataMeasures$Measures$Nb_sample
  chains = coda::mcmc.list() #list of chains
  As = coda::mcmc.list()
  acceptances = list()
  length_bound = list()

  ### Network from Sc
  SC = Sc[-1, ] ## KEEP IT UNTIL FUTHER UPDATES
  network = igraph::graph_from_adjacency_matrix()
  reduced_network = remove_transitive_edges(network)

  ## vizualisation 2 modes
  if ("interactive" %in% list(...)) {

  network_vizualization(reduced_network, DataMeasures$Measures$SampleNames, list(...)$interactive)
  }

  network_vizualization(reduced_network, DataMeasures$Measures$SampleNames)


  ### Index of bounds for each age
  IndexBounds = findbounds(reduced_network) # list of list

  # observation's values
  D = as.numeric(DataMeasures$Measures$D)
  sD = as.numeric(DataMeasures$Measures$sD)
  ddot = as.numeric(DataMeasures$Measures$ddot)
  sddot = as.numeric(DataMeasures$Measures$sddot)
  Ahat = D/ddot

  init = initialize_SC(Sc, LowerPeriod, UpperPeriod, plotGraph, D = D, ddot = ddot, sddot = sddot, sD = sD, type = list(...)$type)
  print(init)
  bounds = apply(IndexBounds, 2, function(b) c(LowerPeriod, init, UpperPeriod)[b]) # 2 x n_ages
  sd = Ahat * sqrt( (sD/D)**2 + (sddot/ddot)**2 ) * proposal_magnitude

  # sd[n_ages] = .1*sd[n_ages]
  print(sd)

  print(data.frame(lower = bounds[1,], A = init, upper = bounds[ 2, ]))

  ### CHAINS ###
  for (c in 1:nchain) {

    if (c > 1) {
      plotGraph = F
    }

    chain <- A <- lbound <- matrix(NA, nrow = niter+1, ncol = n_ages)

    lbound[1,] = bounds[2, ] - bounds[1, ]
    A[1, ] = init
    acceptance = rep(0, n_ages)
    names(acceptance) <- paste("rate", DataMeasures$Measures$SampleNames, sep = "-")

    X = A[1,]

    for (iter in 1:niter) {
      for (i in 1:n_ages) {

        bounds_i = c(LowerPeriod, X, UpperPeriod)[IndexBounds[, i]]
        #proposal
        proposal = truncnorm::rtruncnorm(1,bounds_i[1], bounds_i[2], mean = X[i], sd = sd[i]) #runif(1, min = bounds_i[1], max = bounds_i[2])

        top = GibbsDensity(DataMeasures, X, proposal, i, bounds_i[1], bounds_i[2]) #* dunif(X[i], bounds_i[1], bounds_i[2])
        bottom = GibbsDensity(DataMeasures, X, X[i], i, bounds_i[1], bounds_i[2]) #* dunif(proposal, bounds_i[1], bounds_i[2])
        ratio = top /bottom
      # if (i == n_ages) {
      #   print(ratio)
      # }

        ## MAj
        u = runif(1)
        if (u < ratio) {
          X[i] = proposal
          acceptance[i] = acceptance[i] + 1
        }

        ### Acceptance calibration only for the first 2000 ####
        # if (iter < 20000 & iter%%100 == 0) {
        #   if (acceptance[i] > .8) {
        #     sd[i] = sd[i] * 1.1
        #   }
        #   else {
        #     sd[i] = .9 * sd[i]
        #   }
        # }


      } #### End COORDINATE LOOP

      A[(iter +1), ] =X

      ## MESSAGE FOR EACH THOUSAND ITER
      if (iter%% 5000 == 0) {
        cat(sprintf("\r -------- Chain %d & Iteration: %d -------", c,  iter))
      }

    } #### END OF ITER LOOP

    colnames(A) <- DataMeasures$Measures$SampleNames
    seq_lag = seq(burnin, (niter+1), lag)
    chain = chain[seq_lag, ]
    A = A[seq_lag, ]

    #ADD chains, As, acceptances
    chains[[paste0("chain", c)]] <- coda::mcmc(chain)
    As[[paste0("chain", c)]] <- coda::mcmc(A)
    acceptances[[paste0("acceptance", c)]] <- acceptance /(niter)
    length_bound[[paste0("chain", c)]] <- lbound

  } ### END OF CHAIN LOOP ###

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

  # try(plot(coda::acfplot(As)))



  plot(1:n_ages, seq(0,1, length.out = n_ages), ylab = "acceptation rate", xlab = "samples", type = "n")
  a = sapply(acceptances, lines, col = 1:n_ages)
  a = sapply(acceptances, points, col = 1:n_ages)


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
    ncol = 9,
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
                      "Time Series SE",
                      "Bayes sd"
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
  standardError <- apply(sample, 2, sd)


  R[, c(1,5)] <- round(t(credible95), roundingOfValue)
  R[, c(2,4)] <- round(t(credible68), roundingOfValue)
  R[, 3] <-   round(estimate, roundingOfValue)
  R[, 9] <-   round(standardError, roundingOfValue)

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
    length_bound = length_bound,
    "Summary" = R
  )

  plot_Ages(object = output, model = "Jeffreys prior")

  return(
    output
  )



}






