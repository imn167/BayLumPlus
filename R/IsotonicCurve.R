#'@export
IsotonicCurve <- function(StratiConstraints, object, level = .95, method = "SBM", plotNetwork) {
  #get all mcmc samples
  sample = runjags::combine.mcmc(object$Sampling) ## mcmc sample
  SampleNames = object$Ages$SAMPLE

  n = length(SampleNames)

  w = 1/ as.numeric(object$Summary[, 8])^2 #inv of the estimated variance

  #Build The Network

  ## StratigraphicConstraints
  ##no Strati
  if (length(StratiConstraints) == 0) {
    StratiConstraints <- rbind(rep(1, n), upper.tri(matrix(rep(1), ncol = n, nrow = n))*1)

    IsoSamples = as.matrix(pbapply::pbapply(sample, 1, function(Ahat, weights) {
      Sys.sleep(.003)
      t(Iso::pava(Ahat, w = weights))
      }, weights = w))

  }
  ## Strati
  else{
    if (is(StratiConstraints)[1] == "character") {
      SCMatrix <- read.csv(StratiConstraints, sep = NULL)
      StratiConstraints <- as.matrix(SCMatrix)
    }
  }

  network = remove_transitive_edges(buildNetwork(StratiConstraints))

  print(network_vizualization(network, paste0("A", 1:n), interactive = T))


  ## apply for each age vector the Isotonic Regression
  if (method == "SBM") {
    edges = igraph::as_edgelist(network, names = F)
    IsoSamples = as.matrix(pbapply::pbapply(sample, 1, function(Ahat, edges, weights) {
      Sys.sleep(.003)
      t(IsotoneOptimization::solve_isotone_DAG(Ahat, w = weights, Emat = edges))},
      edges = edges , weights = w))
  }

  else {IsoSamples = as.matrix(pbapply::pbapply(sample, 1, function(Ahat, network, weights) {
    Sys.sleep(.003)
    t(IsotonicRegDAG(network, Ahat, weights )$A)},
    network = network , weights = w))} # n_ages x n_iter
  IsoSamples = t(IsoSamples)
  HPD = apply(IsoSamples, 2, arkhe::interval_hdr, level = level)

  if (is.list(HPD)) {
    #divided_hpd <- sapply(HPD,length)
    HPD = apply(IsoSamples, 2, CredibleInterval, level = level, roundingOfValue = 3)
  }


  IsoSummary = data.frame(lower = HPD[1, ], upper = HPD[2, ], avg = apply(IsoSamples, 2, mean),
                          Samples = factor(SampleNames, levels = SampleNames), Unit = 1:n )

  return(list(chain = IsoSamples, summary = IsoSummary))

}

#'@export
PlotIsotonicCurve <- function(StratiConstraints, object, level = .95, method = "SBM") {

  df <- IsotonicCurve(StratiConstraints, object, level, method)[[2]]
  n = dim(df)[1]

  curve = df %>% ggplot2::ggplot(ggplot2::aes(x = Unit, ymin = lower, ymax = upper), fill = "orange") +
    ggplot2::geom_ribbon(alpha = .4) +
    ggplot2::geom_line(ggplot2::aes(y = lower), color = "orange", group = 1) +
    ggplot2::geom_line(ggplot2::aes(y = upper), color = "orange", group = 1) +
    ggplot2::geom_line(ggplot2::aes(y = avg), color = "orange", group = 1, size =1.5) +
    ggplot2::geom_point(ggplot2::aes(x = Unit, y = lower), color = "blue") +
    ggplot2::geom_point(ggplot2::aes(x = Unit, y = upper), color = "red") +
    ggplot2::geom_point(ggplot2::aes(x = Unit, y = avg), color = "black") +
    BayLumTheme() + ggplot2::ylab("IsotonicRegression") +
    ggplot2::scale_x_continuous(breaks = df$Unit, labels = df$Samples) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45))

  print(curve)
}
