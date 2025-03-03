#' @title Bayesian Models for Age estimation based on OSL measures
#'
#' @description
#'  This function compute the age (in ka) of at least two samples. \cr
#'  The function is based on the output of the [Palaeodose_Computation] and [create_measuresDataFrame]
#'@export
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
    diag(Measures$sD)

  return(list(Theta = Theta, Measures = Measures))
}

#'@describeIn destination description
#'@export

Compute_AgeS_D <- function(
    DATA, SamplesNames = DATA$SamplesNames,
    Nb_Samples = DATA$Nb_Samples,
    THETA = c(),
    StratiConstraints = c(),
    model = NULL,
    Iter = 10000,
    burnin = 4000,
    adapt = 1000,
    t = 5, #thin
    n.chains = 3,
    prior = "Jeffreys",
    jags_method = "rjags",
    autorun = F,
    quit = F,
    roundingOfValue = 3,
    SavePdf = FALSE,
    OutputFileName = c('MCMCplot', "summary"),
    OutputFilePath = c(""),
    SaveEstimates = FALSE,
    OutputTableName = c("DATA"),
    OutputTablePath = c(''),
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
  DataList = list(
    "I" = Measures$Nb_Sample,
    "Theta" = Theta,
    "ddot" = Measures$ddot,
    "StratiConstraints" = StratiConstraints,
    "xbound" = PriorAges
  )

  ## select Model
  if (is.null(model)) {
    model <- AgePrior[[prior]]
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
        monitor = c("A", "Sigma"),
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
        monitor = c("A", "Sigma"),
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
    "CovarianceMatrix" = Theta,
    "model" = model
  )

  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  # JAGS RUN --------------------- END
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

  #---processing of JAGS results
  ##extract mcmc list from runjags object
  echantillon <- results_runjags$mcmc

  ##remove mcmc-list from runjags output to reduce output object size
  results_runjags$mcmc <- list("MCMC-list is not here. Go to first level -> object named *Sampling*")

  ##combine chains into one data.frame
  sample <- as.data.frame(runjags::combine.mcmc(echantillon))

  ##try makes sure that the function runs
  try(plot_MCMC(echantillon, sample_names = SampleNames))

  if (SavePdf) {
    dev.off()
  }

  #---Gelman and Rubin test of convergence of the MCMC ####
  CV <- gelman.diag(echantillon, multivariate = FALSE)
  cat(paste(
    "\n\n>> Results of the Gelman and Rubin criterion of convergence <<\n"
  ))
  for (i in 1:Nb_sample) {
    cat("----------------------------------------------\n")
    cat(paste(" Sample name: ", SampleNames[i], "\n"))
    cat("---------------------\n")
    cat(paste("\t\t", "Point estimate", "Uppers confidence interval\n"))
    cat(paste(
      paste("A_", SampleNames[i], sep = ""),
      "\t",
      round(CV$psrf[i, 1], roundingOfValue),
      "\t\t",
      round(CV$psrf[i, 2], roundingOfValue),
      "\n"
    ))

  }

  cat(
    "\n\n---------------------------------------------------------------------------------------------------\n"
  )
  cat(
    " *** WARNING: The following information are only valid if the MCMC chains have converged  ***\n"
  )
  cat(
    "---------------------------------------------------------------------------------------------------\n\n"
  )

  #---print results ####
  ##Matrix of results
  rnames <- paste0("A_", SamplesNames)

  R <- matrix(
    data = NA,
    ncol = 8,
    nrow = Nb_Samples,
    dimnames = list(rnames,
                  c(
                      "lower bound at 95%",
                      "lower bound at 68%",
                      "Bayes estimate",
                      "upper bound at 68%",
                      "upper bound at 95%",
                      "",
                      "Convergencies: Point estimate",
                      "Convergencies: uppers confidence interval"
                    )
                  )
  )


  #Bayes Estimate and credibal iterval
  cat(
    "\n\n>> Bayes estimates of Age, Palaeodose and its dispersion for each sample and credible interval <<\n"
  )

  credible95 <- apply(sample, 2, CredibleInterval, level = .95)[, 2:3]
  credible68 <- apply(sampe, 2, CredibleInterval, level = .68)[, 2:3]
  estimate <- apply(sample, 2, mean)
  R[, c(1,5)] <- round(credible95, roundingOfValue)
  R[, c(2,4)] <- round(credible68, roundingOfValue)
  R[, 3] <-   round(estimate, roundingOfValue)

  cat("\n----------------------------------------------\n")
  R[, c(7, 8)] <- round(CV$psrf, roundingOfValue)


  #---print csv table ####
  if (SaveEstimates == TRUE) {
    write.csv(R, file = c(
      paste(
        OutputTablePath,
        "Estimates",
        OutputTableName,
        ".csv",
        sep = ""
      )
    ))
  }

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
      SAMPLE = SampleNames,
      AGE = AgePlotMoy,
      HPD68.MIN = credible68[, 2],
      HPD68.MAX = credible68[, 3],
      HPD95.MIN = credible95[, 2],
      HPD95.MAX = credible95[, 3],
      stringsAsFactors = FALSE
    ),
    "Sampling" = echantillon,
    "PriorAge" = results_runjags$args$PriorAge,
    "StratiConstraints" = results_runjags$args$StratiConstraints,
    "CovarianceMatrix" = results_runjags$args$CovarianceMatrix,
    "model" = results_runjags$model,
    "runjags_object" = results_runjags
  )

  #---Plot ages ####
  BayLum::plot_Ages(object = output, legend.pos = "bottomleft")

  ##TODO: get rid of this ... at some point
  if (SavePdf) {
    dev.print(
      pdf,
      file = paste(OutputFilePath, OutputFileName[3], '.pdf', sep = ""),
      width = 8,
      height = 10
    )
  }

  #---Return output ####
  return(output)


}
