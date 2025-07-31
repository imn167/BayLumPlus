#' Compute Isotonic Regression for different stratigraphic constraints
#' @description This function compute the isotonic distorsion of a posterior distrubution obtained by `Compute_AgeS_D()`
#' @param StratiConstraints : The stratigraphic matrix, it can be stored in a file or given directly. If there is none, then the model suppose that it is a strict order
#' @param object : the returned object of the age model function `Compute_AgeS_D()`
#' @param level = 0.95 by default for the level of High Posterior Densities regions
#'
#'@export
IsotonicCurve <- function(StratiConstraints, object, level = .95, method = "SBM", graphPath = file.path(tempdir(), "graph.html"), interactive) {
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


    network = remove_transitive_edges(buildNetwork(StratiConstraints))
    visplot = network_vizualization(network, SampleNames, interactive)

    if (interactive) {

      visNetwork::visSave(visplot, file = graphPath)

      utils::browseURL(graphPath)
    }

    else{

      print(visplot)
    }

    IsoSamples = as.matrix(pbapply::pbapply(sample, 1, function(Ahat, weights) {
      Sys.sleep(.003)
      t(Iso::pava(Ahat, w = weights))
      }, weights = w))
    HPD = apply(IsoSamples, 1, arkhe::interval_hdr, level = level)
    if (is.list(HPD)) {
      message("\t  \t Multiples HPD Regions -- Using Credible Interval Instead \t \t ")
      HPD = apply(IsoSamples, 2, CredibleInterval, level = level, roundingOfValue = 3)
      HPD = HPD[-1, ]
    }
    IsoSamples <- t(IsoSamples)
    colnames(IsoSamples) <-  SampleNames
    IsoSummary = data.frame(lower = HPD[1, ], upper = HPD[2, ], AGE = apply(IsoSamples, 2, mean),
                            SAMPLE = factor(SampleNames, levels = SampleNames), Unit = 1:n )

    return(.list_BayLum(Sampling = coda::as.mcmc.list(coda::as.mcmc(IsoSamples)), Ages = IsoSummary, network= network))
  }
  ## Strati

  else if (is(StratiConstraints)[1] == "character") {
      SCMatrix <- read.csv(StratiConstraints, sep = NULL)
      StratiConstraints <- as.matrix(SCMatrix)
    }



  network = remove_transitive_edges(buildNetwork(StratiConstraints))
  visplot = network_vizualization(network, SampleNames, interactive)

  if (interactive) {

  visNetwork::visSave(visplot, file = graphPath)

  utils::browseURL(graphPath)
  }

  else{

    print(visplot)
  }

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
    message("\t  \t Multiples HPD Regions -- Using Credible Interval Instead \t \t ")
    HPD = apply(IsoSamples, 2, CredibleInterval, level = level, roundingOfValue = 3)
    HPD = HPD[-1, ]
  }
  colnames(IsoSamples) <- SampleNames
  IsoSummary = data.frame(lower = HPD[1, ], upper = HPD[2, ], AGE = apply(IsoSamples, 2, mean),
                          SAMPLE = factor(SampleNames, levels = SampleNames), Unit = 1:n )

  return(.list_BayLum(Sampling = coda::as.mcmc.list(coda::as.mcmc(IsoSamples)), Ages = IsoSummary, network= network))

}

#'@export
PlotIsotonicCurve <- function(StratiConstraints, object, level = .95, method = "SBM", ...) {
  arg = list(...)
  if (!is.null(arg$interactive)) {
    Iso <- IsotonicCurve(StratiConstraints, object, level, method, interactive = arg$interactive)
  }
  else {
  Iso <- IsotonicCurve(StratiConstraints, object, level, method, interactive = T)
  }
  df <- Iso[[2]]
  network <- Iso[[3]]
  n = dim(df)[1]

  curve = df %>% ggplot2::ggplot(ggplot2::aes(x = Unit, ymin = lower, ymax = upper), fill = "orange") +
    ggplot2::geom_ribbon(alpha = .4) +
    ggplot2::geom_line(ggplot2::aes(y = lower), color = "orange", group = 1) +
    ggplot2::geom_line(ggplot2::aes(y = upper), color = "orange", group = 1) +
    ggplot2::geom_line(ggplot2::aes(y = AGE), color = "orange", group = 1, size =1.5) +
    ggplot2::geom_point(ggplot2::aes(x = Unit, y = lower), color = "blue") +
    ggplot2::geom_point(ggplot2::aes(x = Unit, y = upper), color = "red") +
    ggplot2::geom_point(ggplot2::aes(x = Unit, y = AGE), color = "black") +
    BayLumTheme() + ggplot2::ylab("IsotonicRegression") +
    ggplot2::scale_x_continuous(breaks = df$Unit, labels = df$SAMPLE) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

  print(curve)
  SAMPLE = object$Ages$SAMPLE
  tg <- tidygraph::as_tbl_graph(network)
  tg <- tg %>% tidygraph::activate(nodes) %>% tidygraph::mutate(SAMPLE = SAMPLE)
  tg <- tg %>% tidygraph::activate(nodes) %>% tidygraph::left_join(df, by = "SAMPLE") %>%
    tidygraph::mutate(translation = (upper-lower)/2)
  layout <- ggraph::create_layout(tg, layout = "sugiyama") %>% dplyr::mutate(x1 = x-translation, x2 = x + translation, y = -(AGE ))

  dag <- ggraph::ggraph(layout) +
    ggplot2::geom_segment(data = layout, ggplot2::aes(x = x1, xend = x2, y = y, color = AGE), linewidth = 2) + ggraph::theme_graph() +
    ggrepel::geom_text_repel(ggplot2::aes(x = x, y = y, label = SAMPLE), size = 5, max.overlaps = Inf) +
    ggplot2::scale_color_viridis_c(name = "Ages") +
    ggraph::geom_edge_link(arrow = grid::arrow(length = grid::unit(.8, 'mm')),end_cap = ggraph::circle(3, 'mm'), alpha = 0.1)

  print(dag)
  return(Iso)
}






