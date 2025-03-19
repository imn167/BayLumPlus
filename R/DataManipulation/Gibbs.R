####

#'
#'@export
findbound <- function(Sc, index,  UpperPeriod, LowerPeriod) {

  lowerbound_index = which(Sc[, index] == 1)
  if (length(lowerbound_index) > 0) { a = A[min(lowerbound_index)]}

  else { a = LowerPeriod}

  upperbound_index = which(Sc[index, ] == 1)
  if (length(upperbound_index) > 0 ) {b = A[max(upperbound_index)]}
  else {b = UpperPeriod}
  return(list(a = a, b = b))
}

#'@export
GibbsDensity <- function(DataMeasures, A, Ai, index, SC, LowerPeriod, UpperPeriod) {
    ## Using the function create_MeasuresDataFrame
  Measures <- DataMeasures$Measures
  D = Measures$D
  Di = D[index]
  sD = Measures$Ds
  sDi = sD[index]
  ddot = Measures$ddot
  ddoti = ddot[index]
  sddot = Measures$sddot
  sddoti = sddot[index]
  sdc = Measures$sddot_shared
  sdci = sdc[index]
  alpha = Measures$symmetric_error
  alphai = alpha[index]

  ## getting the modidied
  A[index] = Ai

  detCov <- function(A) {
    1 + sum((A*ddot * alpha)**2 / ((A*sddot)**2 + sD**2))
  }

  invdiagvar = ((Ai*sddoti)**2 + sDi**2)**(-1)
  centeri = (Di-Aiddoti)**2

  firstExp = exp( (-.5) * ( (centeri *invdiagvar ) )
  secondExp = exp(  -.5* detCov(A)**(-1) * sum( (Ai*(centeri) * invdiagvar * (D-A*ddot)**2 / ((A*sddot)**2 + sD**2) )))


  return( detCov(A)**(-.5) * invdiagvar * firstExp * secondExp *(Ai > a) * (Ai < b) / Ai)
}

#'@export
logitT <- function(A, upperbound, lowerbound) {

  return(log((A-lowerbound) / (upperbound-A)))
}

#'@export
GibbsSampler <- function(DataMeasures, init, sd, nchain, burnin, Sc,
                      LowerPeriod, UpperPeriod, lag = 10) {

  n_ages = DataMeasures$Measures$Nb_sample
  chain = matrix(NA, nrow = nchain+1, ncol = n_ages)
  ##logit transformation for init

  chain[1, ] = sapply(init, logitT, upperbound = UpperPeriod, )

  for (iter in 1:nchain) {
    for (i in 1:n_ages) {

    }



  }
}









