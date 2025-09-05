#' Compute Isotonic Regression for different stratigraphic constraints
#' @description This function compute the isotonic distorsion of a posterior distrubution obtained by `Compute_AgeS_D()`.\cr
#'The efficient algorithm considered for the Isotnic Regression (IR) is the ** Sequetial Block Merging (SBM)**. \cr
#'If the `StratiConstraints` input is a null matrix then the model suppose a strict order
#'
#' @param StratiConstraints [numeric matrix] or [character] : The stratigraphic relation between samples.
#' @param object [BayLum.list]: output of the function [Compute_AgeS_D()] when using the prior **no_strat**. If the `StratiConstraints` input is a null matrix then the model suppose a strict order
#' @param level [numeric] c(0.95, 0.68) by default for the level of High Posterior Densities (HPD) regions. If the HPD region is composed of more than 1 interval then the model return the Credible Interval
#'  at the level indicated.
#' @param graphPath [character] (Default) path for saving the `StratiConstraits`'s DAG
#' @param interactive [logical] (Default) Indicating wether we want an interactive html file or an plot image for the DAG structure
#'
#'@return ** NUMERICAL OUTPUT : A list of type `BayLum.list` containing the following objects**
#'\enumerate{
#'  \item \bold{Sampling} : Samples from the posterior distribution after distorsion by the isotonic regression;
#'  \item \bold{network} : The DAG constructed from the `StratiConstraints` matrix
#'  \item \bold{Ages} : data frame containing the Credible interval at 95% and 68%, the bayes mean estimator, the bayes standard deviation estimator, the MAP and sample names.
#'}

#' @seealso [Compute_AgeS_D(), PlotIsotonicCurve()]
#'@export


IsotonicCurve <- function(StratiConstraints, object, levels = c(.65,.95), path = tempdir(), interactive) {

  ##### VARIABLES #######@
  #get all mcmc samples
  sample = runjags::combine.mcmc(object$Sampling) ## mcmc sample
  SampleNames = object$Ages$SAMPLE

  n = length(SampleNames)

  w = 1/ as.numeric(object$Summary[, 8])^2 #inv of the estimated variance

  # html path
  graphPath = file.path(path, "graph.html")

  # function to get the map a posterior
    get_map <- function(x) {
      kde = density(x)
      m = min(x)
      M = max(x)
      u = seq(m, M, length.out = 1000)
      values = approx(kde$x, kde$y, xout = u)$y

      map = x[which.max(values)]
    }
  ## Stratigraph in a csv file
  if (is(StratiConstraints)[1] == "character") {
    SCMatrix <- read.csv(StratiConstraints, sep = NULL)
    StratiConstraints <- as.matrix(SCMatrix)
  }


    ## if empty StratiConstraints matrix --> strict order assumption
    if (length(StratiConstraints) == 0) {
      StratiConstraints <- rbind(rep(1, n), upper.tri(matrix(rep(1), ncol = n, nrow = n))*1)
    }

    ## Network vizualisation before

    network = remove_transitive_edges(buildNetwork(StratiConstraints))
    visplot = network_vizualization(network, SampleNames, interactive)

    if (interactive) {

      visNetwork::visSave(visplot, file = graphPath)

      utils::browseURL(graphPath)
    }

    else{

      print(visplot)
    }

    ## init dataframe for HPD regions
    HPD_frame <- data.frame()
    ## Isotonic Regression
    edges = igraph::as_edgelist(network, names = F)
    IsoSamples = as.matrix(pbapply::pbapply(sample, 1, function(Ahat, edges, weights) {
      Sys.sleep(.003)
      t(IsotoneOptimization::solve_isotone_DAG(Ahat, w = weights, Emat = edges))},
      edges = edges , weights = w))
    IsoSamples = t(IsoSamples)

    for (elt in levels) {
      HPD = apply(IsoSamples, 2, arkhe::interval_hdr, level = level)
      if (is.list(HPD)) {
        message(paste("\t  \t Multiples HPD Regions -- Using", elt*100, " Credible Interval Instead \t \t "))
        HPD = apply(IsoSamples, 2, CredibleInterval, level = level, roundingOfValue = 3)
        HPD = HPD[-1, ]
      }

      #create column names
      min_col <- paste0("HPD.", elt, "MIN")
      max_col <- paste0("HPD", elt, "MAX")

      if (is.matrix(HPD) && nrow(HPD) == 2) {
        HPD_frame[[min_col]] <- HPD[1,]
        HPD_frame[[max_col]] <- HPD[2,]}
    }
    colnames(IsoSamples) <-  SampleNames
    IsoSummary = data.frame( AGE = apply(IsoSamples, 2, mean), SD = apply(IsoSamples, 2, sd), MAP = apply(IsoSamples, 2, get_map),
                             SAMPLE = factor(SampleNames, levels = SampleNames), Unit = 1:n )
    IsoSummary = cbind(IsoSummary, HPD_frame)

    return(.list_BayLum(Sampling = coda::as.mcmc.list(coda::as.mcmc(IsoSamples)), Ages = IsoSummary, network= network))

    }



#' Plot graphs build to vizualise istonic distorsion of the function [IsotonicCurve()]
#' @param StratiConstraints [numeric matrix] or [character] : The stratigraphic relation between samples.
#'
#' @param object [BayLum.list]: output of the function [Compute_AgeS_D()] when using the prior **no_strat**. If the `StratiConstraints` input is a null matrix then the model suppose a strict order
#' @param level [numeric] c(0.95, 0.68) by default for the level of High Posterior Densities (HPD) regions. If the HPD region is composed of more than 1 interval then the model return the Credible Interval
#'  at the level indicated.
#' @param ...
#'
#'@return ** NUMERICAL OUTPUT : A list of type `BayLum.list` containing the following objects**
#'\enumerate{
#'  \item \bold{ribbon} : Ribbon of the Credible Interval and solid orange line representing the Bayes mean estimator;
#'  \item \bold{dag} : The DAG constructed from the `StratiConstraints` matrix with segments length represent IC(level) and gradient color the value of ages for each samples
#'  \item \bold{Iso} : returned value of the function [IsotonicCurve()]
#'}
#'
#' @importFrom rlang .data
#'@export
#'
PlotIsotonicCurve <- function(StratiConstraints, object, level = .95,  ...) {
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

  curve = df %>% ggplot2::ggplot(ggplot2::aes(x = .data$Unit, ymin = .data$lower, ymax = .data$upper), fill = "orange") +
    ggplot2::geom_ribbon(alpha = .4) +
    ggplot2::geom_line(ggplot2::aes(y = .data$lower), color = "orange", group = 1) +
    ggplot2::geom_line(ggplot2::aes(y = .data$upper), color = "orange", group = 1) +
    ggplot2::geom_line(ggplot2::aes(y = .data$AGE), color = "orange", group = 1, size =1.5) +
    ggplot2::geom_point(ggplot2::aes(x = .data$Unit, y = .data$lower), color = "blue") +
    ggplot2::geom_point(ggplot2::aes(x = .data$Unit, y = .data$upper), color = "red") +
    ggplot2::geom_point(ggplot2::aes(x = .data$Unit, y = .data$AGE), color = "black") +
    BayLumTheme() + ggplot2::ylab("IsotonicRegression") + ggplot2::xlab("Samples")
    ggplot2::scale_x_continuous(breaks = df$Unit, labels = df$SAMPLE)  +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) + ggplot2::coord_flip()

  print(curve)
  SAMPLE = object$Ages$SAMPLE
  tg <- tidygraph::as_tbl_graph(network)
  tg <- tg %>% tidygraph::activate(nodes) %>% tidygraph::mutate(SAMPLE = SAMPLE)
  tg <- tg %>% tidygraph::activate(nodes) %>% tidygraph::left_join(df, by = "SAMPLE") %>%
    tidygraph::mutate(translation = (.data$upper-.data$lower)/2)
  layout <- ggraph::create_layout(tg, layout = "sugiyama") %>% dplyr::mutate(x1 = .data$x-.data$translation, x2 = .data$x + .data$translation, y = -.data$AGE )

  dag <- ggraph::ggraph(layout) +
    ggplot2::geom_segment(data = layout, ggplot2::aes(x = .data$x1, xend = .data$x2, y = .data$y, color = .data$AGE), linewidth = 2) + ggraph::theme_graph() +
    ggrepel::geom_text_repel(ggplot2::aes(x = .data$x, y = .data$y, label = .data$SAMPLE), size = 5, max.overlaps = Inf) +
    ggplot2::scale_color_viridis_c(name = "Ages") +
    ggraph::geom_edge_link(arrow = grid::arrow(length = grid::unit(.8, 'mm')),end_cap = ggraph::circle(3, 'mm'), alpha = 0.1)

  print(dag)
  return(.list_BayLum(Iso = Iso, DAG = dag, ribbon = curve))
}






