#' @title Create Input for Age Computation
#'
#' @description
#' This function prepares the necessary data frame entries for the function [Computation_AgeS_D].
#' It takes structured input data (as defined in [Generate_DataFile()] and [Generate_DataFile_MG()])
#' and computes relevant parameters for age estimation.
#'
#' @param Data A structured dataset following the format required for OSL sample analysis.
#' @param P A Palaeodose_computation object
#' @param symetric_error The $\alpha$ parameter in Combes & Philippe (2017)
#'@param contamination_degree The common error $\sigma_c$ in Combes & Philippe (2017)
#' @return A list containing:
#' \itemize{
#'   \item \strong{Sigma}: The covariance matrix computed as:
#'   \deqn{\Theta = A \Sigma A + \text{diag}(sD)}
#'   \item \strong{Info}: A list with details including:
#'   \itemize{
#'     \item Number of samples (`Nb_sample`)
#'     \item Sample names (`NamesOfSamples`)
#'     \item Dose rate values (`ddot`)
#'     \item Dose rate uncertainties (`sddot`)
#'     \item Estimated palaeodose values (`D`)
#'     \item Estimated palaeodose uncertainties (`sD`)
#'   }
#' }
#'
#' @seealso [Computation_AgeS_D], [Palaeodose_Computation], [Generate_DataFile], [Generate_DataFile_MG]
#' @export
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
                         ddot = DATA$ddot_env[1, ], sddot = DATA$ddot_env[2, ],
                         D = Obs[1:DATA$Nb_sample],
                         sD = Obs[(DATA$Nb_sample+1): (2*DATA$Nb_sample)],
                   sddot_shared = contamination_degree, symetric_error = symetric_error)
  Theta = diag(Measures$sddot) + (contamination_degree %*% t(contamination_degree)) * symetric_error
  CovD = diag(Measures$sD)

  return(list(Theta = Theta, Measures = Measures, covD = covD))
}
sepSC <- NULL
#' @title Bayesian Models for Age Estimation with Priors from ModelAgePrior Dataset
#'
#' @description
#' This function computes the second part of Combes & Philippe (2017), specifically the age estimation
#' using a Bayesian model. Its behavior is similar to other functions like [AgeS_Computation()],
#' with the primary difference being the first parameter.
#'
#' @param DATAMeasures **(required)** [list]
#' The output of the function [create_MeasuresDataFrame()], containing the necessary input data for computation.
#'
#' @return See the documentation for [AgeS_Computation()] for details on the returned output.
#'
#' @seealso [AgeS_Computation()], [create_MeasuresDataFrame()]
#' @md
#' @export


Compute_AgeS_D <- function(
    DATAMeasures,
    StratiConstraints = c(),
    model = NULL,
    Iter = 10000,
    burnin = 4000,
    adapt = 1000,
    t = 5, #thin
    n.chains = 3,
    prior = "Jeffreys",
    PriorAge = rep(c(0.01, 100), Measures$Nb_sample),
    jags_method = "rjags",
    autorun = F,
    quiet = F,
    roundingOfValue = 3,
    SavePdf = FALSE,
    OutputFileName = c('MCMCplot', "summary"),
    OutputFilePath = c(""),
    SaveEstimates = FALSE,
    OutputTableName = c("DATA"),
    OutputTablePath = c(''),
    ...
) {

  Measures = DATAMeasures$Measures


  ## StratigraphicConstraints
  ##no Strati
  if (length(StratiConstraints) == 0) {
    StratiConstraints <- matrix(
      data = c(rep(1, Measures$Nb_sample), rep(0, Measures$Nb_sample * Measures$Nb_sample)),
      ncol = Measures$Nb_sample,
      nrow = (Measures$Nb_sample + 1),
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
  dataList = list(
    "I" = Measures$Nb_sample,
    "Theta" = DATAMeasures$Theta,
    "covD" = DATAMeasures$covD,
    "ddot" = Measures$ddot,
    "StratiConstraints" = StratiConstraints,
    "xbound" = PriorAge,
    "D" = Measures$D
  )

  ModelAgePrior <- 0
  data(ModelAgePrior, envir = environment())

  ## select Model
  if (is.null(model)) {
    model <- ModelAgePrior[[prior]]
  }

  ## first way to run jags model : manual
  if (!(autorun)) {
    #write model in tempfile
    temp_file <- tempfile(fileext = ".txt")
    writeLines(model, con = temp_file)
    if ( prior == "Jeffreys") {

  inits = list(
    list(u = runif(Measures$Nb_sample)), #chain 1
    list(u = runif(Measures$Nb_sample)), # chain 2
    list(u = runif(Measures$Nb_sample)) #chain 3
  )
    }

  else if ( prior == "StrictOrder") {
    temp_file <- tempfile(fileext = ".txt")
    writeLines(model, con = temp_file)
    inits = list(
     list( e = rexp(Measures$Nb_sample + 1)), #chain 1
     list( e = rexp(Measures$Nb_sample+ 1)), #chain 2
      list(e = rexp(Measures$Nb_sample+ 1)) #chain 3
    )
  }
    #run JAGS
    results_runjags <-
      runjags::run.JAGS(
        model = temp_file,
        data = dataList,
        n.chains = n.chains,
        monitor = c("A"),
        adapt = adapt,
        burnin = burnin,
        sample = Iter,
        silent.jags = quiet,
        method = jags_method,
        thin = t,
        inits = inits
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
    "PriorAge" = PriorAge,
    "StratiConstraints" = StratiConstraints,
    "CovarianceMatrix" = DATAMeasures$Theta,
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
  try(plot_MCMC(echantillon, sample_names = Measures$SampleNames))

  #autocorrelation diagnosis
  cat("\n\n ------------------------------------------------------------------------------\n\n")
  cat(paste(
    "\n\n>> Results of sampling autocorrelation <<\n"
  ))
  print(coda::autocorr.diag(echantillon))

  try(plot(coda::acfplot(echantillon)))

  if (SavePdf) {
    dev.off()
  }

  cat("\n\n ------------------------------------------------------------------------------\n\n")

  #---Gelman and Rubin test of convergence of the MCMC ####
  CV <- gelman.diag(echantillon, multivariate = FALSE)
  cat(paste(
    "\n\n>> Results of the Gelman and Rubin criterion of convergence <<\n"
  ))
  for (i in 1:Measures$Nb_sample) {
    cat("----------------------------------------------\n")
    cat(paste(" Sample name: ", Measures$SampleNames[i], "\n"))
    cat("---------------------\n")
    cat(paste("\t\t", "Point estimate", "Uppers confidence interval\n"))
    cat(paste(
      paste("A_", Measures$SampleNames[i], sep = ""),
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
  rnames <- paste0("A_", Measures$SampleNames)

  R <- matrix(
    data = NA,
    ncol = 7,
    nrow = Measures$Nb_sample,
    dimnames = list(rnames,
                  c(
                      "lower bound at 95%",
                      "lower bound at 68%",
                      "Bayes estimate",
                      "upper bound at 68%",
                      "upper bound at 95%",
                      "Convergencies: Point estimate",
                      "Convergencies: uppers confidence interval"
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

  print(data.frame(R) )
  cat("\n----------------------------------------------\n")


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
      SAMPLE = Measures$SampleNames,
      AGE = estimate,
      HPD68.MIN = credible68[ 1,],
      HPD68.MAX = credible68[2, ],
      HPD95.MIN = credible95[1, ],
      HPD95.MAX = credible95[2, ],
      stringsAsFactors = FALSE
    ),
    "Sampling" = echantillon,
    "PriorAge" = results_runjags$args$PriorAge,
    "StratiConstraints" = results_runjags$args$StratiConstraints,
    "CovarianceMatrix" = results_runjags$args$CovarianceMatrix,
    "model" = results_runjags$model,
    "runjags_object" = results_runjags
  )

  cat("\n ===================================\n")
  print(list(SummaryRunJAGS = summary(results_runjags)))
  cat("\n==============================\n")

  #---Plot ages ####
  BayLum::plot_Ages(object = output, legend.pos = "bottomleft", model = paste("BayLum", prior))

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
