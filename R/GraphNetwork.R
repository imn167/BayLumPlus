
#' The Edge Pruner Algorithm (EPA). This algorithm handles the redundancy of the constraints matrix (see [create_ThetaMatrix()]) by deleting all unecessary edges.
#'The result is
#' @param Sc  [numeric matrix] or [character] : The stratigraphic relation between samples.
#' @param n_samples [numeric] NaN by default. Optional size of samples to create a strict order DAG without constructing `Sc`
#'
#' @returns [igraph]  a non transitive DAG that can also be referred to as a Hasse Diagram.
#' @author Imène Bouafia
#' @export
#'
#' @examples
#' #Example for strict order
#' network <- remove_transitive_edges(c(), n_samples = 5)
#' print(network)
remove_transitive_edges <- function(Sc, n_samples = NaN) {
  #build network
  G = buildNetwork(Sc, n_samples)

  #EdgePruner Algorithm
  reduced_G = rlang::duplicate(G)
  vertices = igraph::V(G)
  for (u in vertices) {
    # message(paste("----------------- traitement sommet", u, "--------------"))
    #get all the neighbors of the vertices u
    u_neighbors = igraph::neighbors(G, u, mode = "out")
    #look for the descendants of each neighbors
    for (nei in u_neighbors) {
      # message(paste("treatement for the neighbours", nei))
      childs = igraph::subcomponent(G, nei, mode = "out")[-1]
      for (child in childs) {
        if (igraph::are_adjacent(reduced_G, u, child)) {
          # print(igraph::E(reduced_G, c(u, child)))
          reduced_G = igraph::delete_edges(reduced_G, igraph::E(reduced_G, P = c(u,child)))
        }
      }
    }
  }
  return(reduced_G)
}




#' Plot a Graph Network according to sugiyama layout
#'
#' @param network [igraph network] An igraph object
#' @param vertices_labels [character] Samples Names for each node
#' @param interactive [bool] TRUE by default, An optional parameter whether to use the interactive visualization with VisNetwork
#'
#' @returns \enumerate{
#' \item If *interactive = TRUE* : `visNetwork` plot
#' \item If *interactive = FALSE* : `ggraph` plot
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' network = remove_transitive_edges(c(), n_samples = 5)
#' network_vizualization(network, paste0("A", 1:5))
#' }
network_vizualization <- function(network, vertices_labels, interactive = TRUE) {

  ##layout for stratigraphic constraints
  layout = igraph::layout_with_sugiyama(network)$layout
  n = length(igraph::V(network))
  if (interactive) {
    nodes <- data.frame(id = 1:n,
                        label = vertices_labels,
                        size = 25,
                        shape = "dot",
                        font = list(size = 20, face = "bold"),
                        stringsAsFactors = FALSE,
                        x = layout[, 1] * 100,
                        y = -layout[, 2] * 100   # invert Y for visNetwork
    )

    edges = igraph::as_data_frame(network, what = "edges")
    names(edges)[1:2] <- c("from", "to")

    visNetwork::visNetwork(nodes, edges, width = "100%", height = "90vh") %>%
      visNetwork::visNodes(font = list(size=20, align = "center")) %>%
      visNetwork::visEdges(arrows = "to") %>%
      visNetwork::visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visNetwork::visLayout(randomSeed = 123)

  }


  else {
    tg <- tidygraph::as_tbl_graph(network)
    tg <- tg %>% tidygraph::activate(nodes) %>% tidygraph::mutate(Samples = vertices_labels)

    layout <- ggraph::create_layout(tg, layout = "sugiyama")

    ggraph::ggraph(layout) + ggraph::geom_edge_link(arrow = grid::arrow(length = grid::unit(.8, 'mm')),end_cap = ggraph::circle(3, 'mm'), alpha = 0.5, edge_colour = "red") +
      ggraph::theme_graph() + ggraph::geom_node_circle(ggplot2::aes(r = .05), fill = "lightyellow", color = "blue", linewidth = 1)  +
      ggrepel::geom_text_repel(ggplot2::aes(x = .data$x, y = .data$y, label = .data$Samples), size = 3.5, max.overlaps = Inf)
    # plot(
    #   network,
    #   layout = layout,
    #   vertex.label = vertices_labels,
    #   vertex.size = 10,
    #   vertex.color = adjustcolor("lightblue", alpha.f = 0.6),
    #   edge.arrow.size = 0.4,  # Smaller arrowheads
    #   edge.width = 2,
    #   edge.arrow.length = 10,
    #   asp = 0,
    #   edge.curved = 0.1
    # )
  }
}

