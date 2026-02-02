#' Compute Isotonic Regression for different stratigraphic constraints
#' @description This function compute the isotonic distorsion of a posterior distrubution obtained by [Compute_AgeS_D()].\cr
#'The efficient algorithm considered for the Isotnic Regression (IR) is the **Sequetial Block Merging (SBM)**. \cr
#'If the `StratiConstraints` input is a null matrix then the model suppose a strict order
#'
#' @param StratiConstraints [matrix] or [character] : The stratigraphic relation between samples.
#' @param object [list]: output of the function [Compute_AgeS_D()] when using the prior **no_strat**. If the `StratiConstraints` input is a null matrix then the model suppose a strict order
#' @param interactive [logical] (Default) Indicating wether we want an interactive html file or a plot image for the DAG structure
#' @param levels [numeric] c(0.95, 0.68) by default for the level of High Posterior Densities (HPD) regions. If the HPD region is composed of more than 1 interval then the model return the Credible Interval
#'  at the level indicated.
#' @param path [character] (Default) path for saving the `StratiConstraints`'s DAG
#' @param rounding_digits [integer] (Default) digits for rounding estimated values.
#'
#'@return **NUMERICAL OUTPUT : A list of type `BayLum.list` containing the following objects**
#'\enumerate{
#'  \item \bold{Sampling} : Samples from the posterior distribution after distorsion by the isotonic regression;
#'  \item \bold{network} : The DAG constructed from the `StratiConstraints` matrix
#'  \item \bold{Ages} : data frame containing the Credible interval at 95% and 68% , the bayes mean estimator, the bayes standard deviation estimator, the MAP and sample names.
#'}

#' @seealso [Compute_AgeS_D()] [PlotIsotonicCurve()]
#'@export


IsotonicCurve <- function(StratiConstraints, object, interactive, levels = c(.68,.95), path = tempdir(), rounding_digits = 3) {

  ##### VARIABLES #######@
  ## Forbid user from using anything other than unconstrained_Jeffrey
  if(!is.nan(object$prior)) {
    if(!startsWith(tolower(object$prior), "unconstrained")) stop(paste("The Posterior should be computed with the prior unconstrained_Jeffrey not", object$prior))
  }

  #get all mcmc samples
  sample = runjags::combine.mcmc(object$Sampling) ## mcmc sample
  SampleNames = object$Ages$SAMPLE

  n = length(SampleNames)

  w = 1/ as.numeric(object$Summary$`Bayes sd`)^2 #inv of the estimated variance
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
    SCMatrix <- read.csv(StratiConstraints, sep = ",")
    StratiConstraints <- as.matrix(SCMatrix)
  }


    ## if empty StratiConstraints matrix --> strict order assumption
    if (length(StratiConstraints) == 0) {
      StratiConstraints <- rbind(rep(1, n), upper.tri(matrix(rep(1), ncol = n, nrow = n))*1)
    }

    ## Network vizualisation before

    network = remove_transitive_edges(StratiConstraints)
    visplot = network_vizualization(network, SampleNames, interactive)

    if (interactive) {

      visNetwork::visSave(visplot, file = graphPath)

      utils::browseURL(graphPath)
    }

    else{

      print(visplot)
    }

    ## init dataframe for HPD regions
    HPD_frame <- data.frame(SAMPLE = SampleNames)
    ## Isotonic Regression
    edges = igraph::as_edgelist(network, names = F)
    IsoSamples = as.matrix(pbapply::pbapply(sample, 1, function(Ahat, edges, weights) {
      Sys.sleep(.003)
      t(IsotoneOptimization::solve_isotone_DAG(Ahat, w = weights, Emat = edges))},
      edges = edges , weights = w))
    IsoSamples = t(IsoSamples)

    for (elt in levels) {
      HPD = apply(IsoSamples, 2, arkhe::interval_hdr, level = elt)
      if (is.list(HPD)) {
        message(paste("\t  \t Multiples HPD Regions -- Using", elt*100, " Credible Interval Instead \t \t "))
        HPD = apply(IsoSamples, 2, CredibleInterval, level = elt, roundingOfValue = 12)
      }
      HPD = HPD[-3, ]
      HPD = round(HPD, digits = rounding_digits)

      #create column names
      min_col <- paste0("HPD", elt*100, ".MIN")
      max_col <- paste0("HPD", elt*100, ".MAX")

      if (is.matrix(HPD) && nrow(HPD) == 2) {

        HPD_frame[[min_col]] <- HPD[1,]
        HPD_frame[[max_col]] <- HPD[2,]}
    }

    colnames(IsoSamples) <-  SampleNames
    IsoSummary = data.frame( AGE = round(apply(IsoSamples, 2, mean), digits = rounding_digits),
                             SD = round(apply(IsoSamples, 2, sd), digits = rounding_digits), MAP = round(apply(IsoSamples, 2, get_map), digits = rounding_digits),
                             SAMPLE = factor(SampleNames, levels = SampleNames), Unit = 1:n )
    IsoSummary = dplyr::left_join(IsoSummary, HPD_frame) #cbind(IsoSummary, HPD_frame)

    return(.list_BayLum(Sampling = coda::as.mcmc.list(coda::as.mcmc(IsoSamples)), Ages = IsoSummary, network= network))

    }



#' Plot graphs build to vizualise isotonic distorsion of the function [IsotonicCurve()].
#' The user needs to note that only one HPD regions / Credible Interval level can be allowed and vizualized.
#' @param object [list]: output of the function [IsotonicCurve()] when using the prior **no_strat**.
#' @param level [numeric] 0.95 by default. One of the computed level given in [IsotonicCurve()]  Notice that only one level can be considered.
#'
#'@return ** NUMERICAL OUTPUT : A list of type `BayLum.list` containing the following objects**
#'\enumerate{
#'  \item \bold{ribbon} : Ribbon of the Credible Interval and solid orange line representing the Bayes mean estimator;
#'  \item \bold{dag} : The DAG constructed from the `StratiConstraints` matrix with segments length represent IC(level) and gradient color the value of ages for each samples
#'  \item \bold{Iso} : Data frame used for generating curve and DAG plots. This data frame is extracted from the **Ages** element of the list object returned by [IsotonicCurve()]
#'}
#'
#' @importFrom rlang .data
#'@export
#'
PlotIsotonicCurve <- function( object, level = .95) {
  # arg = list(...)
  # if (!is.null(arg$interactive)) {
  #   Iso <- IsotonicCurve(StratiConstraints, object, level, interactive = arg$interactive)
  # }
  # else {
  # Iso <- IsotonicCurve(StratiConstraints, object, level, interactive = T)
  # }
  #

  if (length(level)>1) {
    stop("This function accepts only a single HPD level. Please provide one credible interval level instead of multiple.", call. = FALSE)
  }

  #check if the HPD columns exist
  hpd_pattern <-  paste0("HPD", level*100)
  has_hpd <- any(grepl(hpd_pattern, colnames(object$Ages)))

  if (!has_hpd) {
    stop(paste0("No HPD regions / Credible Interval found with the level", level *100, "%"))
  }

  df <- object[[2]] %>% dplyr::select(Unit, SAMPLE, AGE,dplyr::starts_with(paste0("HPD", level*100))) %>%
    dplyr::rename_with(~c("lower", "upper"), dplyr::starts_with("HPD"))
  network <- object[[3]]
  n = dim(df)[1]
  ## reversing the order for the ribbon and Graph plot
  df <- df %>% dplyr::mutate(SAMPLE = factor(.data$SAMPLE, levels = rev(.data$SAMPLE)), Unit = rev(.data$Unit))

  curve = df  %>% ggplot2::ggplot(ggplot2::aes(x = .data$Unit, ymin = .data$lower, ymax = .data$upper)) +
    ggplot2::geom_ribbon(alpha = .4, fill = "orange") +
    ggplot2::geom_line(ggplot2::aes(y = .data$lower), color = "orange", group = 1) +
    ggplot2::geom_line(ggplot2::aes(y = .data$upper), color = "orange", group = 1) +
    ggplot2::geom_line(ggplot2::aes(y = .data$AGE), color = "orange", group = 1, linewidth =1.5) +
    ggplot2::geom_point(ggplot2::aes(x = .data$Unit, y = .data$lower), color = "blue") +
    ggplot2::geom_point(ggplot2::aes(x = .data$Unit, y = .data$upper), color = "red") +
    ggplot2::geom_point(ggplot2::aes(x = .data$Unit, y = .data$AGE), color = "black") +
    BayLumTheme() + ggplot2::ylab("IsotonicRegression") + ggplot2::xlab("Samples") +
    ggplot2::scale_x_continuous(breaks = df$Unit, labels = df$SAMPLE)  +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) +
    ggplot2::coord_flip() ## Flipped coordinates but needed to change the ordering ?

  # print(curve)

  # SAMPLE = object$Ages$SAMPLE
  tg <- tidygraph::as_tbl_graph(network)
  tg <- tg %>% tidygraph::activate(nodes) %>% tidygraph::mutate(SAMPLE = df$SAMPLE)
  tg <- tg %>% tidygraph::activate(nodes) %>% tidygraph::left_join(df, by = "SAMPLE") %>%
    tidygraph::mutate(translation = (.data$upper-.data$lower)/2)
  layout <- ggraph::create_layout(tg, layout = "sugiyama") %>% dplyr::mutate(x1 = .data$x-.data$translation, x2 = .data$x + .data$translation, y = -.data$AGE )

  dag <- ggraph::ggraph(layout) +
    ggplot2::geom_segment(data = layout, ggplot2::aes(x = .data$x1, xend = .data$x2, y = .data$y, color = .data$AGE), linewidth = 2) + ggraph::theme_graph() +
    ggrepel::geom_text_repel(ggplot2::aes(x = .data$x, y = .data$y, label = .data$SAMPLE), size = 5, max.overlaps = Inf) +
    ggplot2::scale_color_viridis_c(name = "Ages") +
    ggraph::geom_edge_link(arrow = grid::arrow(length = grid::unit(.8, 'mm')),end_cap = ggraph::circle(3, 'mm'), alpha = 0.1)

  # print(dag)

  if (interactive()) {
    curve
    dag
  }

  return(.list_BayLum(Iso = df, DAG = dag, ribbon = curve))
}






