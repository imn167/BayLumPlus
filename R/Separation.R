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
    model <- Model_Prior[[prior]] #### CAREFULL THE NAME ISN'T QUITE RIGHT
  }


  ## first way to run jags model : manual
  if (!(autorun)) {
    #write model in tempfile
    temp_file <- tempfile(fileext = ".txt")
    writeLines(model, con = temp_file)

    #run JAGS
    results_runjags <-
      runjags::run.JAGS(
        model = temp_file,
        data = dataList,
        n.chains = n.chains,
        monitor = c("A", "D", "sD"),
        adapt = adapt,
        burnin = burnin,
        sample = Iter,
        silent.jags = quiet,
        method = jags_method,
        thin = t
      )

  }

  #second way to run a jags model : automatic
  if (autorun) {
    ##further settings provided eventually
    process_settings <- modifyList(x = list(
      max.time = Inf,
      interactive = FALSE,
      startburnin = 4000,
      startsample = 10000,
      inits = NA

    ), val = list(...))

    ##a text file is wanted as input, so we have to cheat a little bit
    temp_file <- tempfile(fileext = ".txt")
    writeLines(model, con = temp_file)

    ##run the autoprocessing
    results_runjags <-
      runjags::autorun.JAGS(
        model = temp_file,
        data = dataList,
        n.chains = n.chains,
        monitor = c("A", "D", "sD"),
        adapt = adapt,
        startburnin = process_settings$startburnin,
        startsample = process_settings$startsample,
        silent.jags = quiet,
        method = jags_method,
        thin = t,
        inits = process_settings$inits,
        max.time = process_settings$max.time,
        interactive = process_settings$interactive
      )

  }

  # storing the arguments used for the BayLum-run this way,
  # because it allows us an easy way to code the storage of arguments when extending a JAGS-model.
  results_runjags$args <- list(
    "PriorAge" = Model_Prior,
    "StratiConstraints" = StratiConstraints,
    "CovarianceMatrix" = THETA,
    "model" = model
  )



}
