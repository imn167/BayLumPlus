#### HPDPlot ####
#'@description
#'@author Imène Bouafia (LMJL)
#'@md
#'@export

plotHpd <- function(
    list_object,
    ModelNames,
    level = .95
                    ) {

  if (length(ModelNames) != length(list_object)) {
    stop("[plotHPD()] Wrong input. 'ModelNames' should have the same length as 'list_object'", call. = FALSE )
  }

  ## dataframe of all points for each model
  DtHpd <- data.frame()
  sampleNames <- list_object[[1]]$Ages$SAMPLE

  for (i in seq_along(list_object)) {
  object = list_object[[i]]
  if (is.null(attributes(object)$class) || attributes(object)$class != "BayLum.list")
    stop("[plotHPD()] Wrong input, only objects of type 'BayLum.list' are allowed. Please check the manual!",
         call. = FALSE
    )
  samples = runjags::combine.mcmc(object$Sampling)
  Estimate = object$Age$AGE ## return estimate (mean)

  #if user wants other names
  colnames(samples) <- sampleNames
  DtHpd <- DtHpd %>% dplyr::bind_rows(hpd_method(ModelNames[i], samples))

  }

  ##

   DtHpd <- DtHpd %>% dplyr::mutate(Samples = factor(Samples, levels = sampleNames))
  #
  plotting <- DtHpd %>% ggplot2::ggplot(ggplot2::aes(x = Samples, ymin= inf, ymax = sup, colour = Models)) +
    ggplot2::geom_linerange(position = ggplot2::position_dodge(.5), ) + ICAgeTheme(rotation_x = T)

  return(plotting)
}









