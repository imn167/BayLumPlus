#' @title Create Input for Age Computation
#'@md
#' @description
#' This function prepares the necessary data frame entries for the function [Compute_AgeS_D()].
#' It takes structured input data (as defined in [Generate_DataFile()] and [Generate_DataFile_MG()])
#' and computes relevant parameters for age estimation.
#'
#' @param Data A structured dataset following the format required for OSL sample analysis.
#' @param P A Palaeodose_computation object
#' @param symetric_error The \eqn{\alpha}{ascii} parameter in Combes & Philippe (2017)
#'@param contamination_degree The common error \eqn{\sigma}{ascii} in Combes & Philippe (2017)
#' @return A list containing:
#' * **Sigma**: The covariance matrix computed as:
#'   \deqn{\Theta = A \Sigma A + \text{diag}(sD)}{Theta = A * Sigma * A + diag(sD)}
#' * **Info**: A list with details including:
#'   * Number of samples (`Nb_sample`)
#'   * Sample names (`NamesOfSamples`)
#'   * Dose rate values (`ddot`)
#'   * Dose rate uncertainties (`sddot`)
#'   * Estimated palaeodose values (`D`)
#'   * Estimated palaeodose uncertainties (`sD`)
#' @seealso \code{\link{Computation_AgeS_D}}, \code{\link{Palaeodose_Computation}}, \code{\link{Generate_DataFile}}, \code{\link{Generate_DataFile_MG}}
#'
#'


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
  covD = diag(Measures$sD**2)

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
#' @param DATA **(required)** [list] of objects :
#' - Equivalent doses `D`
#' - dose rates : `ddot`
#' - standard error for D : `sD`
#' The output of the function [create_MeasuresDataFrame()], containing the necessary input data for computation.
#' @param Nb_sample [integer] number of samples
#' @param SampleNames [character] character vector with sample names
#' @param ThetaMatrix [matrix] or [character] input of systematic and individual errors.
#' @param PalaeodoseObject [list] (with default NULL) Output of the [Palaeodose_Computation()]
#' @param StratiConstraints [matrix] or [character] The stratigraphic relation between samples
#' @param model [character] (optional) custom Jags model
#' @param Iter (with default): the number of iterations to run which will be used to assess convergence and ages (see [runjags::run.jags]).
#' @param burnin  [integer] (with default): the number of iterations used to "home in" on the stationary posterior distribution. These are not used for assessing convergence (see [runjags::run.jags]).
#' @param adapt [integer] (with default): the number of iterations used in the adaptive phase of the simulation (see [runjags::run.jags]).
#' @param t [integer] (with default): 1 every \code{t} iterations of the MCMC is considered for sampling the posterior distribution.
#' (for more information see [runjags::run.jags]).
#' @param n.chains [integer] (with default): number of independent chains for the model (for more information see [runjags::run.jags]).
#' @param prior [character] : Character string specifying the name of one of the models
#'   available in the `ModelAgePrior` dataset. Use [extract_Jags_model()] to see all available options
#' @param PriorAge vector (with default): lower and upper bounds for age parameter of each sample (in ka).
#' @param jags_method (with default): select which method to use in order to call JAGS. jags_methods `"rjparallel"`  (the default) and `"rjags"` have been tested. (for more information about these possibilities and others, see [runjags::run.jags])
#' @param autorun [logical] (with default): choose to automate JAGS processing. JAGS model will be automatically extended until convergence is reached (for more information see [runjags::autorun.jags]).
#' @param quiet  [logical] (with default): enables/disables `rjags` messages
#' @param roundingOfValue [integer] (with default):  Integer indicating the number of decimal places to be used, default = 3.
#' @param SavePdf [logical] (with default): if TRUE save graphs in pdf file named `OutputFileName` in folder `OutputFilePath`.
#' @param OutputFileName  OutputFileName [character] (with default): name of the pdf file that will be generated by the function if `SavePdf = TRUE`;
#' \code{length(OutputFileName)=2}, see \bold{PLOT OUTPUT} in \bold{Value} section for more informations.
#' @param OutputFilePath [character] (with default): path to the pdf file that will be generated by the function if \code{SavePdf}=TRUE.
#' If it is not equal to "", it must be terminated by "/".
#' @param SaveEstimates [logical] (with default): if TRUE save Bayes' estimates, credible interval at level 68% and 95%,
#' the result of the Gelman en Rubin test of convergence and the Time Series SE, in a csv table named \code{OutputFileName} in folder \code{OutputFilePath}.
#' @param OutputTableName [character] (with default): name of the table that will be generated by the function if `SaveEstimates = TRUE`.
#' @param OutputTablePath [character] (with default): path to the table that will be generated by the function if `SaveEstimates = TRUE`.
#' If it is not equal to "", it must be terminated by "/".
#' @param ...
#'
#' @return
#' **NUMERICAL OUTPUT**
#'
#'  **A list of type `BayLum.list` containing the following objects:**
#'  \itemize{
#'  \item \bold{Ages} : dataframe containing the Credible interval at 95% and 68%, the bayes mean estimator, the bayes standard deviation estimator and sample names.
#'   \item \bold{Sampling}: that corresponds to a sample of the posterior distributions
#'  of the age (in ka);
#'   \item \bold{prior}: category of prior used. `prior == unconstrained` if no stratigraphic constraints;
#'   \item \bold{PriorAge}: stating the priors used for the age parameter (in ka);
#'   \item \bold{StratiConstraints}: stating the stratigraphic relations between samples considered in the model;
#'   \item \bold{CovarianceMatrix}: stating the covariance matrix of error used in the model, highlighting common errors between samples or not;
#'   \item \bold{model}: returns the model that was used for the Bayesian modelling as a [character];
#'   \item \bold{JAGS model output}: returns the JAGS model with class "runjags";
#'   \item \bold{Summary}: Summary Table of the posterior's MCMC;
#'  }
#'  **PLOT OUTPUT**
#'
#' \enumerate{
#'  \item\bold{MCMC trajectories}: A graph with the MCMC trajectories and posterior distributions of the age.
#' On each line, the plot on the left represents the MCMC trajectories, and the one on the right the posterior distribution of the parameter.
#'  \item \bold{Summary of sample age estimates}: plot credible intervals and Bayes estimate of each sample age on a same graph. \cr
#'  A default `plot_mode` parameter can gives either IC segment or density distribution if *plot_mode = "density"*, see [plot_Ages()]
#' }
#' @details ** Which prior to use regarding the Stratigraphic constraints **
#'If there is a strict order as a stratigraphic constraints, the user would be able to use the following priors :
#'\enumerate{
#'  \item \bold{constrained_Jeffrey} : uniform order over the period of study;
#'  \item \bold{old_BayLum} : old BayLum model : false configuration of the approximated Jeffrey;
#'  \item \bold{StrictNicholls&Jones} : The Nicholls & Jones uniform order applied on ages;
#'  \item \bold{unconstrained_jeffrey} : The approached Jeffrey (see Combes & Philippe 2017) without stratigraphic constraints;
#'
#'}
#'@examples
#'DATA <- list(D = OSLJingbian$D, sD = OSLJingbian$sD, ddot = OSLJingbian$ddot)
#'\dontrun{
#' ### run standard
#'  Output <- Compute_AgeS_D(DATA = DATA,
#'  Nb_sample = OSLJingbian$Nb_Sample,
#'  SampleNames = OSLJingbian$SampleNames,
#'  ThetaMatrix = OSLJingbian$ThetaMatrix,
#'  ‡prior = "unconstrained_Jeffrey",
#'  PriorAge = rep(c(1, 1400), OSLJingbian$Nb_Sample),
#'  Iter = 2000, burnin = 50000, t = 10)
#'
#'
#'}
#'
#' @seealso [AgeS_Computation()]
#' @md
#' @export


Compute_AgeS_D <- function(
    DATA,
    Nb_sample,
    SampleNames,
    ThetaMatrix,
    PalaeodoseObject = NULL,
    StratiConstraints = c(),
    model = NULL,
    Iter = 10000,
    burnin = 4000,
    adapt = 1000,
    t = 5, #thin
    n.chains = 3,
    prior = "unconstrained_jeffrey",
    PriorAge = rep(c(0.01, 100), Nb_sample),
    jags_method = "rjparallel",
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

  ### read csv file of Palaeodose output




  ### JagsRun
  ## liste of data
   if (prior == "constrained_Jeffrey" | prior == "StrictNicholls&Jones" | prior == "unconstrained_Jeffrey" | prior =="Conditional"  ) {
     dataList = list(
       "I" = Nb_sample,
       "Theta" = ThetaMatrix,
       "covD" = diag(DATA$sD**2),
       "ddot" = DATA$ddot,
       "xbound" = PriorAge,
       "D" = DATA$D
     )
   }

  else if (prior == "unconstrained" | prior == "uniform_order" | prior == "Nicholls&Jones") {
    dataList = list(
      "I" = Nb_sample,
      "M" = DATA$M,
      "xbound" = PriorAge,
      "sigma" = DATA$sigma
    )
  }

  else {

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


    dataList = list(
      "I" = Nb_sample,
      "Theta" = ThetaMatrix,
      "covD" = diag(DATA$sD**2),
      "ddot" = DATA$ddot,
      "StratiConstraints" = StratiConstraints,
      "xbound" = PriorAge,
      "D" = DATA$D
    )

  }

  ModelAgePrior <- 0
  data(ModelAgePrior, envir = environment())

  ## select Model
  if (is.null(model)) {
    entry_model = model
    model <- ModelAgePrior[[prior]]
  }

  ## first way to run jags model : manual
  if (!(autorun)) {
    #write model in tempfile
    temp_file <- tempfile(fileext = ".txt")
    writeLines(model, con = temp_file)
    if ( prior == "oldBayLum"  | prior == "unconstrained_Jeffrey"| prior == "unconstrained") {

  inits = replicate(n.chains, list(u = runif(Nb_sample)), simplify = F)
    }

  else if ( prior == "constrained_Jeffrey" | prior == "uniform_order") {

    inits = replicate(n.chains, list( e = rexp(Nb_sample + 1)), simplify = F)
  }

    else if (prior == "StrictNicholls&Jones" | prior == "Nicholls&Jones") {    ######

      inits = replicate(n.chains, nichollsInit(Nb_sample, 1, 0), simplify = F)

    }

    else if (prior == "nichollsBR") {    ######
      inits = replicate(n.chains, nichollsBRInit(Nb_sample, 1, 0), simplify = F)
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
        monitor = c("A"),
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
    "CovarianceMatrix" = ThetaMatrix,
    "model" = model
  )

  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  # JAGS RUN --------------------- END
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

  #---processing of JAGS results
  ##extract mcmc list from runjags object
  echantillon <- results_runjags$mcmc
  SummaryMCMC <- summary(echantillon)

  ##remove mcmc-list from runjags output to reduce output object size
  results_runjags$mcmc <- list("MCMC-list is not here. Go to first level -> object named *Sampling*")

  ##combine chains into one data.frame
  sample <- as.data.frame(runjags::combine.mcmc(echantillon))

  #autocorrelation diagnosis
  cat("\n\n ------------------------------------------------------------------------------\n\n")
  cat(paste(
    "\n\n>> Results of sampling autocorrelation <<\n"
  ))
  print(coda::autocorr.diag(echantillon))


  #---plot MCMC ####
  if (SavePdf) {
    pdf(file = paste(OutputFilePath, OutputFileName[1], '.pdf', sep = ""))
  }

  ##try makes sure that the function runs
  plot(1:10, 1:10, type = "n", axes = FALSE, xlab = "", ylab = "")
  text(5.5, 9, "MCMC Plot", cex = 2, font = 2)  # Large, bold text
  try(plot_MCMC(echantillon, sample_names = SampleNames))
  plot(1:10, 1:10, type = "n", axes = FALSE, xlab = "", ylab = "")
  text(5.5, 9, "Autocorr Plot", cex = 2, font = 2)
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
  rnames <- paste0("A_", SampleNames)

  R <- matrix(
    data = NA,
    ncol = 10,
    nrow = Nb_sample,
    dimnames = list(rnames,
                  c(
                      "lower bound at 95%",
                      "lower bound at 68%",
                      "Bayes estimate",
                      "upper bound at 68%",
                      "upper bound at 95%",
                      "Bayes sd",
                      "Convergencies: Point estimate",
                      "Convergencies: uppers confidence interval",
                      "Time Series SE",
                      "Geweke Criteria pvalue"
                    )
                  )
  )


  #Bayes Estimate and credibal iterval
  cat(
    "\n\n>> Bayes estimates of Age, Palaeodose and its dispersion for each sample and credible interval <<\n"
  )

  credible95 <-  apply(sample, 2, CredibleInterval, level = .95)[ 1:2, ]
  credible68 <- apply(sample, 2, CredibleInterval, level = .68)[1:2, ]
  estimate <- apply(sample, 2, mean)
  standardError <- apply(sample, 2, sd)
  geweke = geweke.diag(as.mcmc(sample))$z
  pvalue = 2*apply(rbind(pnorm(geweke), pnorm(geweke, lower.tail = F)), 2, min)

  R[, c(1,5)] <- round(t(credible95), roundingOfValue)
  R[, c(2,4)] <- round(t(credible68), roundingOfValue)
  R[, 3] <-   round(estimate, roundingOfValue)

  R[, c(7, 8)] <- round(CV$psrf, roundingOfValue)
  R[, 6] <- round(standardError, roundingOfValue)
  R[, 9] <- round(SummaryMCMC$statistics[, 4], roundingOfValue)
  R[, 10] <- pvalue

  R <- dplyr::as_tibble(R)
  R <- R %>% dplyr::mutate(GewekeStars = stratify_pvalue(pvalue))


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
      SAMPLE = SampleNames,
      AGE = estimate,
      SD = standardError,
      HPD68.MIN = credible68[ 1,],
      HPD68.MAX = credible68[2, ],
      HPD95.MIN = credible95[1, ],
      HPD95.MAX = credible95[2, ],
      stringsAsFactors = FALSE
    ),
    "Sampling" = echantillon,
    "prior" = if(is.null(entry_model)) prior else NaN,
    "PriorAge" = results_runjags$args$PriorAge,
    "StratiConstraints" = results_runjags$args$StratiConstraints,
    "CovarianceMatrix" = results_runjags$args$CovarianceMatrix,
    "model" = results_runjags$model,
    "runjags_object" = results_runjags,
    "Summary" = R
  )

  cat("\n ===================================\n")
  print(list(SummaryRunJAGS = summary(results_runjags)))
  cat("\n==============================\n")

  #---Plot ages ####
  plot_Ages(object = output, legend.pos = "bottomleft", model = paste("BayLum", prior))

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
