library(ArchaeoPhases)


Path=c("R/DataManipulation/")

########### sans micocouliers et sans UCIAMS-109995 (reconnu par BayLum comme un outlier)


C14_SampleNames = c("OxA-21261","OxA-21262","UCIAMS-98210","UCIAMS-103134","P-782","UCIAMS-103138","OxA-23247",
                    "OxA-9774","OxA-9946","UCIAMS-98208","UCIAMS-98209","OxA-9947","OxA-9775","OxA-9948","OxA-23523",
                    "UCIAMS-103135","UCIAMS_103136","UCIAMS-103137",'OxA-23248',"UCIAMS-103141","UCIAMS-103140",
                    "OxA-27087","OxA-9949","UCIAMS-103139","OxA-23249","OxA-9950","UCIAMS-109991",
                    "OxA-23250","UCIAMS-109992","OxA-9776","OxA-23251","UCIAMS-109993","OxA-9892",
                    "OxA-9777","UCIAMS-109994","OxA-9778","OxA-23252",
                    "OxA-9893","PL-980525A","AA-27982")

C14_Nb_sample = length(C14_SampleNames)

C14ages = c(8033,7955,7940,7955,8092,7920,8027,7935,7980,7965,7990,7985,8090,8090,7931,7970,7940,7965,8082,7980,
            8025,8000,8050,7970,8024,8030,8035,8085,8030,7985,8137,8160,8150,8160,8210,8240,
            8199,8155,8390,8195)

C14agesEr = c(39,40,30,25,98,25,37,50,55,25,25,50,55,50,38,25,25,30,37,25,25,50,40,25,35,50,30,36,30,55,36,
              30,50,50,30,55,36,50,90,80)

AC14_WithStratiWithout109995= AgeC14_Computation(Data_C14Cal=C14ages, Data_SigmaC14Cal=C14agesEr,
                                                SampleNames=C14_SampleNames, Nb_sample = C14_Nb_sample,
                                                PriorAge = rep(c(7, 13), C14_Nb_sample), SavePdf
                                                = TRUE,
                                                monitors = c("Age"),
                                                OutputFileName = c("MCMCplot",
                                                                   "HPD_CalC-14Curve", "summary"),
                                                OutputFilePath = Path, SaveEstimates = TRUE,
                                                OutputTableName = c("AllC14"), OutputTablePath
                                                = c(""),
                                                StratiConstraints = c(), sepSC = c(","), Model       #contraintes stratigraphiques
                                                = c("full"),
                                                CalibrationCurve = c("IntCal20"), Iter = 50000,
                                                t = 5,
                                                n.chains = 3, quiet = FALSE)


network <-igraph::graph_from_adjacency_matrix( Sc)
reduced_network = remove_transitive_edges(network)

sample = runjags::combine.mcmc(AC14_WithStratiWithout109995$Sampling)[1:200, ]

iso = as.matrix(pbapply::pbapply(sample, 1, function(Ahat, network, weights) {
  Sys.sleep(.003)
  t(IsotonicRegDAG(network, Ahat, weights )$A)},
  network = reduced_network , weights = rep(1, 40)))

iso = t(iso)

hpd = apply(iso, 2, CredibleInterval, level = .68, roundingOfValue = 3)
AC14_WithStratiWithout109995$Ages$SAMPLE

IsoC14 <- IsotonicCurve(reduced_network, object = AC14_WithStratiWithout109995, level = .68)
IsoC14$summary


## creation du graphe
tg = tidygraph::as_tbl_graph(reduced_network)
tg <- tg %>% tidygraph::activate(nodes) %>% tidygraph::mutate(Samples = C14_SampleNames)

tg <- tg %>% tidygraph::activate(nodes) %>% tidygraph::left_join(IsoC14$summary, by = "Samples") %>%
  tidygraph::mutate(translation = (upper-lower)/2)
layout <- ggraph::create_layout(tg, layout = "sugiyama") %>% dplyr::mutate(x1 = x-translation, x2 = x + translation)

ggraph::ggraph(layout) + ggraph::geom_edge_link(arrow = grid::arrow(length = grid::unit(.8, 'mm')),end_cap = ggraph::circle(3, 'mm'), alpha = 0.2) +
  ggplot2::geom_segment(data = layout, ggplot2::aes(x = x1, xend = x2, y = y, color = avg), size = 2) + ggraph::theme_graph() +
  ggplot2::geom_text(ggplot2::aes(x = x, y = y, label = Samples), vjust = -1, size = 3)



creape_subPlot <- function(x1, x2) {

  ggplot2::ggplot() + ggplot2::geom_segment(ggplot2::aes(x = .5, xend = .5, y = x1, yend = x2), size = 2) +
    ggplot2::theme_void()
}

layout


ggraph::ggraph(tg, layout = "sugiyama") + ggraph::geom_node_point() + ggraph::theme_graph() ##ggreppel

edges <- cbind(1:(C14_Nb_sample-1), 2:C14_Nb_sample)
G = igraph::graph_from_edgelist(edges)

plot(G,
     vertex.label = igraph::V(G)$name,
     vertex.color = "lightblue",
     vertex.size = 5,
     edge.arrow.size = 0.5,
     layout = igraph::layout_with_sugiyama(G))


PlotIsotonicCurve(G, AC14_WithStratiWithout109995, level = .68)
AC14_WithStratiWithout109995$


