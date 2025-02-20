library(BayLum)

#write a function to create DATA for Age computation
#write Models in Nimble or Jags

create_MeasuresDataFrame <- function(
    PalaeodoseObject,
    DATA,
    symetric_error,
    contamination_degree
) {
  MCMCSamples <- runjags::combine.mcmc(PalaeodoseObject$Sampling)
  Obs <- apply(MCMCSamples, 2, median)
  #creating the dataframe for AgeEstimation
  Measures <- list(SampleNames = DATA$SampleNames, Nb_sample = DATA$Nb_sample,
                         ddot = DATA$ddot[, 1], sddot = DATA$ddot[, 2],
                         D = Obs[1:DATA$Nb_sample],
                         sD = Obs[(DATA$Nb_sample+1): (2*DATA$Nb_sample)])
  Theta = diag(Measures$sddot) + (contamination_degree %*% t(contamination_degree)) * symetric_error +
    Diag(Measures$sD)

  return(list(Theta = Theta, Measures = Measures))
}



Compute_AgeS_D <- function(
    DATA, SamplesNames = DATA$SamplesNames,
    Nb_Samples = DATA$Nb_Samples,
    THETA = c(),
    StratiConstraints = c(),
    model = NULL,
    Iter = 10000,
    burnin = 4000,
    adapt = 1000,
    t = 5,
    n.chains = 3,
    prior = "Jeffreys",
    jags_method = "rjags",
    autorun = F,
    quit = F,
    roundingOfValue = 3,
    ...
) {

  ## StratigraphicConstraints
  ##no Strati
  if (length(StratiConstraints) == 0) {
    StratiConstraints <- matrix(
      data = c(rep(1, Nb_sample), rep(0, Nb_sample * Nb_sample)),
      ncol = Nb_sample,
      nrow = (Nb_sample + 1),
      byrow = T
    )
  }
  ## Strati
  else{
    if (is(StratiConstraints)[1] == "character") {
      SCMatrix <- read.csv(StratiConstraints, sep = sepSC)
      StratiConstraints <- as.matrix(SCMatrix)
    }
  }

  ### JagsRun
  ## liste of data
  DataList <— list(
    "I" = Measures$Nb_Sample,
    "Theta" = Theta,
    "ddot" = Measures$ddot,
    "StratiConstraints" = StratiConstraints,
    "xbound" = PriorAges #why this name
  )

  ## select Model
  if (is.null(model)) {
    model <- Model_Prior[[prior]]
  }



}
