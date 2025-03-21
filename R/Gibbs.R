
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
  return( detCov(A)**(-1/2) * invdiagvar *firstExp * secondExp * (Ai<= b)* (Ai>=a) )
}


#=================================================================================@
#'@export
logitT <- function(A, bounds) {
  lowerbound = bounds[1]
  upperbound = bounds[2]
  return(log((A-lowerbound) / (upperbound-A)))
}


#=================================================================================@
#' @export
arctanT <- function(A, bounds) {
  lowerbound = bounds[1]
  upperbound = bounds[2]
  return( atan( pi* (A-lowerbound) / (upperbound - lowerbound) - pi /2 ) )
}

#=================================================================================@
#'@export
initialize_SC <- function(Sc, LowerPeriod, UpperPeriod, plotGraph = F) {
  Sc = Sc[-1, ] #del. first line
  n = nrow(Sc)
  rownames(Sc) = colnames(Sc) = paste0("A", 1:n)

  network = igraph::graph_from_adjacency_matrix(Sc)
  if (plotGraph) {
    plot(network)
  }
  A = rep(NA, n)
  ##topological order
  topo_order = igraph::topo_sort(network, mode = "out")
  for (i in topo_order) {
    max_bound = UpperPeriod
    min_bound = LowerPeriod
    #MAJ max_bound
    older_ages = A[which(Sc[i, ] == 1 )]
    if ( any(!is.na(older_ages))) {
      max_bound = min(max_bound, min(older_ages, na.rm = T))
    }

    #MAJ min_bound
    younger_ages = A[which(Sc[, i] == 1)]
    if (any(!is.na(younger_ages ))) {
      min_bound = max(min_bound, max(younger_ages, na.rm =T))
    }

    ##simulation
    if (min_bound < max_bound) {
      A[i] = runif(1, min_bound, max_bound)
    }

    else { A[i] = min_bound}
  }
  return(A)
}


#=================================================================================@
#'@export
GibbsSampler <- function(DataMeasures,  sd, nchain, burnin, Sc,
                      LowerPeriod, UpperPeriod, Transformation = "arctan", lag = 10, plotGraph = T) {

  n_ages = DataMeasures$Measures$Nb_sample
  chain = matrix(NA, nrow = nchain+1, ncol = n_ages)
  A = chain #nchains+1 x n_ages
  #initialize
  init = initialize_SC(Sc, LowerPeriod, UpperPeriod, plotGraph)

  A[1, ] = init

  acceptance = 0

  IndexBounds = sapply(1:n_ages, findbound, Sc) + 1 # 2 x n_ages
  IndexBounds[which(IndexBounds == 0, arr.ind = T)] = n_ages + 2
  bounds = apply(IndexBounds, 2, function(b) c(LowerPeriod, init, UpperPeriod)[b])

  if (Transformation == "logit") { #Test Sigmoid with parameter lambda
  ##logit transformation for init
  chain[1, ] = sapply(init, logitT, bounds) # n_ages
  }

  else if (Transformation == "arctan") {
    chain[1, ] = sapply(init, arctanT, bounds)
  }

  for (iter in 1:nchain) {
    for (i in 1:n_ages) {
      bounds_i = c(LowerPeriod, A[iter, ], UpperPeriod)[IndexBounds[, i]] # (a(iter), b(iter))

      #proposal depends on the transformation
      if (Transformation == "arctan") {
      # MAJ Ai
      proposal = chain[iter, i] + truncnorm::rtruncnorm(1, a = -pi/2, b = pi/2, sd = sd[i]) #proposal for the i-th Age in (-pi, pi)/2
      Aproposal = (bounds_i[2] - bounds_i[1]) * tan(proposal)/pi + (bounds_i[2] + bounds_i[1])/2

      top = GibbsDensity(DataMeasures, A[iter, ], Aproposal, i, bounds_i[1], bounds_i[2]) * (1+tan(proposal)**2)
      bottom = GibbsDensity(DataMeasures, A[iter, ], A[iter, i], i, bounds_i[1], bounds_i[2]) * (1+ tan(chain[iter, i])**2)

      log_ratio = log(top) - log(bottom)

      u = log(runif(1))



      if (u < log_ratio) {
        chain[iter+1, i] = proposal
        A[iter+1, i] = Aproposal
        acceptance = acceptance +1
      }

      else {
        chain[iter+1, i] = chain[iter, i]
        A[iter+1, i] = A[iter, i]
        }
      }##### END OF ARCTAN TRANSFORMATION


    }
  ## message
    if (iter%/% 1000) {
      message(paste("Iteration iter =", iter, "done"))
    }
  }
  seq_lag = seq(burnin, (nchain+1), lag)
  chain = chain[seq_lag, ]
  A = A[seq_lag, ]

  ##plot

  return(list(A = A, chain = chain, acceptance = acceptance / (n_ages * nchain)))
}









