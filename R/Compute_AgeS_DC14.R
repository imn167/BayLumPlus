
#' Bayesian Models for Age Estimation with Priors from Model_DC14 Dataset
#'
#' @param **(required)** [list] of objects :
#' - Vector of Equivalent doses and C14 Calibrated Ages in the right order according to stratigraphy
#' - dose rates : `ddot`
#' - standard error for D : `sD`
#' - Error of the C14 Calibrated age : `SigmaC14Cal`
#' @param Nb_sample [integer] number of samples
#' @param SampleNames [character] character vector with sample names
#' @param encoding
#' @param ThetaMatrix [matrix] or [character] input of systematic and individual errors.
#' @param sepTHETA
#' @param prior [character] : Character string specifying the name of one of the models
#'   available in the `Model_DC14` dataset. Use [extract_Jags_model()] to see all available options
#' @param model [character] (optional) custom Jags model
#' @param sepSC
#' @param StratiConstraints [matrix] or [character] The stratigraphic relation between samples
#' @param PriorAge vector (with default): lower and upper bounds for age parameter of each sample (in ka)
#' @param SavePdf
#' @param OutputFileName
#' @param OutputFilePath
#' @param SaveEstimates
#' @param OutputTableName
#' @param OutputTablePath
#' @param Model_C14
#' @param CalibrationCurve
#' @param MCMC_plots
#' @param Iter (with default): the number of iterations to run which will be used to assess convergence and ages (see [runjags::run.jags]).
#' @param burnin  [integer] (with default): the number of iterations used to "home in" on the stationary posterior distribution. These are not used for assessing convergence (see [runjags::run.jags]).
#' @param adapt [integer] (with default): the number of iterations used in the adaptive phase of the simulation (see [runjags::run.jags]).
#' @param t [integer] (with default): 1 every \code{t} iterations of the MCMC is considered for sampling the posterior distribution.
#' (for more information see [runjags::run.jags]).
#' @param n.chains [integer] (with default): number of independent chains for the model (for more information see [runjags::run.jags]).
#' @param jags_method
#' @param autorun
#' @param quiet
#' @param roundingOfValue
#' @param ...
#'
#' @returns
#' @export
#'
#' @examples
Compute_AgeS_DC14 <- function(
    DATA, #Data_C14Cal, Data_SigmaC14Cal
    Nb_sample ,
    SampleNames,
    encoding,
    ThetaMatrix = c(),
    sepTHETA = c(','),
    prior = "Unconstrained",
    model = NULL,
    sepSC = c(','),
    StratiConstraints = c(),
    PriorAge = rep(c(1, 60), Nb_sample),
    SavePdf = FALSE,
    OutputFileName = c('MCMCplot', 'HPD_Cal14CCurve', "summary"),
    OutputFilePath = c(""),
    SaveEstimates = FALSE,
    OutputTableName = c("DATA"),
    OutputTablePath = c(''),
    Model_C14 = c("full"),
    CalibrationCurve = c("IntCal20"),
    MCMC_plots = F,
    Iter = 10000,
    burnin = 4000,
    adapt = 1000,
    t = 5,
    n.chains = 3,
    jags_method = "rjparallel",
    autorun = FALSE,
    quiet = FALSE,
    roundingOfValue = 3,
    ...
) {

  ## sigma_D For Equivalent doses D when mixing
  encode_sD = rep(0, Nb_sample)
  encode_sD[encoding ==1] = DATA$sD**2

  # ## Measure vector for both Equivalent doses and calibrated C14 data
  # M = rep(0, Nb_sample)
  # M[encoding ==1] = DATA$D
  # M[encoding == 0] = DATA$Data_C14_Cal

  ## index in C14 references only
  order_osl = 1:sum(encoding)
  order_c14 = 1:sum(1-encoding)

  CS_osl = cumsum(encoding)
  CS_c14 = cumsum(1-encoding)

  index_osl = which(encoding == 1)
  index_c14 = which(encoding == 0)

  # grid = expand.grid(index_osl, index_osl)
  Theta_Tilde = matrix(0, Nb_sample, Nb_sample)
  Theta_Tilde[index_osl, index_osl] = ThetaMatrix

  #--- Calibration curve ####
  TableauCalib = c()
  if (CalibrationCurve %in%  c(
    "IntCal13", "IntCal20", "Marine13", "Marine20", "SHCal13" , "SHCal20")) {
    TableauCalib = get(data(list = CalibrationCurve, envir = environment()))
  } else {
    TableauCalib = read.csv(file = CalibrationCurve,
                            sep = ",",
                            dec = ".")
  }

  AgeBP = rev(TableauCalib[, 1])
  CalC14 = rev(TableauCalib[, 2])
  SigmaCalC14 = rev(TableauCalib[, 3])

  if (prior == "Constrained" || prior == "NichollsJones") {

    dataList = list(
      "I" = Nb_sample,
      "Theta" = Theta_Tilde,
      "encoded_covD" = diag(encode_sD),
      "ddot" = DATA$ddot,
      "xbound" = PriorAge,
      "M" = DATA$M,
      "sigma" = DATA$sigmaC14Cal,
      "order_osl" = order_osl,
      "order_c14" = order_c14,
      "index_osl" = index_osl,
      "index_c14" = index_c14,
      # "StratiConstraints" = StratiConstraints,
      # "CS_osl" = CS_osl,
      # "CS_c14" = CS_c14,
      "xTableauCalib" = AgeBP,
      "yTableauCalib" = CalC14,
      "zTableauCalib" = SigmaCalC14
    )
  }


  else if ( prior == "Unconstrained") {

  dataList = list(
    "Theta" = Theta_Tilde,
    "encoded_covD" = diag(encode_sD),
    "ddot" = DATA$ddot,
    "xbound" = PriorAge,
    "M" = DATA$M,
    "sigma" = DATA$sigmaC14Cal,
    "order_osl" = order_osl,
    "order_c14" = order_c14,
    "index_osl" = index_osl,
    "index_c14" = index_c14,
    # "StratiConstraints" = StratiConstraints,
    # "CS_osl" = CS_osl,
    # "CS_c14" = CS_c14,
    "xTableauCalib" = AgeBP,
    "yTableauCalib" = CalC14,
    "zTableauCalib" = SigmaCalC14
  )
  }


    Model_DC14 <- 0
    data(Model_DC14, envir = environment())
    if (is.null(model)) {
      entry_model = model
      model <- Model_DC14[[prior]]
    }


    #
    if (!(autorun)) {
      #write model in tempfile
      temp_file <- tempfile(fileext = ".txt")
      writeLines(model, con = temp_file)
      if ( prior == "Unconstrained") {

        inits = replicate(n.chains, list(u = runif(Nb_sample)), simplify = F)
      }

      else if (  prior == "Constrained") {

        inits = replicate(n.chains, list( e = rexp(Nb_sample + 1)), simplify = F)
      }

      else if (prior == "NichollsJones") {    ######

        inits = replicate(n.chains, nichollsInit(Nb_sample, 1, 0), simplify = F)

      }

      results_runjags <- runjags::run.JAGS(
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

      # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
      # JAGS RUN --------------------- END
      # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

      #---processing of JAGS results
      ##extract mcmc list from runjags object
      echantillon <- results_runjags$mcmc
      SummaryMCMC <- summary(echantillon)

      ##combine all of them
      sample = as.data.frame(runjags::combine.mcmc(echantillon))

      # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
      #                           START OF MCMC'S DIAGNOSIS
      # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

      diagnostics_plots = create_diagnostic_plots(echantillon, SampleNames)
      CV <- gelman.diag(echantillon, multivariate = FALSE)
      cat(paste(
        "\n\n>> Results of the Gelman and Rubin criterion of convergence <<\n"
      ))
      print(CV)

      #optional display for plots
      if (MCMC_plots) {
        display_diagnostic_plots(diagnostic_plots)
      }

      if (SavePdf) {
        pdf_path = paste0(OutputFilePath, OutputTableName[1], ".pdf")
        save_diagnostic_plots(diagnostic_plots, pdf_path)
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
        # "prior" = if(is.null(entry_model)) prior else NaN,
        "PriorAge" = PriorAge,
        "StratiConstraints" = StratiConstraints,
        "CovarianceMatrix" = ThetaTilde,
        "model" = model,
        "diagnostics_plots" = diagnostics_plots,
        # "runjags_object" = results_runjags,
        "Summary" = R
      )


      #---Plot ages ####
      plot_Ages(object = output, legend.pos = "bottomleft", model = paste("BayLum", prior), plot_mode = "density")

      return(output)

}
