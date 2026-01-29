

Compute_AgeS_DC14 <- function(
    DATA, #Data_C14Cal, Data_SigmaC14Cal
    Nb_sample ,
    SampleNames,
    encoding,
    ThetaMatrix = c(),
    sepTHETA = c(','),
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
    Iter = 10000,
    burnin = 4000,
    adapt = 1000,
    t = 5,
    n.chains = 3,
    jags_method = "rjags",
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
  index_c14 = which(encodinf == 0)

  grid = expand.grid(index_osl, index_osl)
  Theta_Tilde = matrix(0, Nb_sample, Nb_sample)
  Theta_Tilde[grid$Var1, grid$Var2] = ThetaMatrix

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

  dataList = list(
    "I" = Nb_sample,
    "Theta" = Theta_Tilde,
    "encoded_covD" = diag(encode_sD),
    "ddot" = DATA$ddot,
    "xbound" = PriorAge,
    "M" = DATA$M,
    "sigma" = DATA$Data_SigmaC14Cal,
    "order_osl" = order_osl,
    "order_c14" = order_c14,
    "index_osl" = index_osl,
    "index_c14" = index_c14,
    "StratiConstraints" = StratiConstraints
    "CS_osl" = CS_osl,
    "CS_c14" = CS_c14,
    "xTableauCalib" = AgeBP,
    "yTableauCalib" = CalC14,
    "zTableauCalib" = SigmaCalC14
  )


    ModelAgePrior <- 0
    data(ModelAgePrior, envir = environment())

    if ()






}
