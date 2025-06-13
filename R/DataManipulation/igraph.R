library(igraph)
## matrice des contraintes Catalhöyük
SC = matrix(data=0,ncol=40,nrow=41)
SC[1,]=rep(1,40)
SC[2,]=c(rep(0,1),rep(1,1),rep(0,3),rep(1,35)) #4853 (OxA -2161)
SC[3,]=c(rep(0,5),rep(1,35)) #4861 (OxA-2162)
SC[4,]=c(rep(0,5),rep(1,35)) #4555 (UCIAMS-98210)
SC[5,]=c(rep(0,5),rep(1,35)) #4779 (UCIAMS-103134)
SC[6,]=rep(0,40) #P782 After (Hearth Level 10)
SC[7,]=c(rep(0,9),rep(1,31)) #4850
SC[8,]=c(rep(0,9),rep(1,31)) #4850
SC[9,]=rep(0,40) #4715 After
SC[10,]=rep(0,40) #4715 After
SC[11,]=c(rep(0,17),rep(1,23)) #4517
SC[12,]=c(rep(0,17),rep(1,23)) #4517
SC[13,]=rep(0,40) #4822 After
SC[14,]=rep(0,40) #4826 After
SC[15,]=rep(0,40) #4826 After
SC[16,]=c(rep(0,16),rep(1,24)) #4828
SC[17,]=c(rep(0,16),rep(1,24)) #4828
SC[18,]=c(rep(0,17),rep(1,1),rep(0,4),rep(1,1),rep(0,1),rep(1,16)) #4836
SC[19,]=c(rep(0,22),rep(1,1),rep(0,1),rep(1,16)) #4837
SC[20,]=c(rep(0,20),rep(1,2),rep(0,1),rep(1,17)) #4869
SC[21,]=c(rep(0,20),rep(1,2),rep(0,1),rep(1,17)) #4869
SC[22,]=c(rep(0,22),rep(1,18)) #4867
SC[23,]=c(rep(0,22),rep(1,18)) #4867
SC[24,]=c(rep(0,24),rep(1,16)) #4848
SC[25,]=c(rep(0,24),rep(1,16)) #4865
SC[26,]=c(rep(0,25),rep(1,15)) #4878
SC[27,]=c(rep(0,26),rep(1,14)) #5276
SC[28,]=c(rep(0,27),rep(1,13)) #5279
SC[29,]=c(rep(0,29),rep(1,11)) #5283
SC[30,]=c(rep(0,29),rep(1,11)) #5283
SC[31,]=c(rep(0,30),rep(1,10)) #5292
SC[32,]=c(rep(0,31),rep(1,9)) #5308
SC[33,]=c(rep(0,32),rep(1,8)) #5316
SC[34,]=c(rep(0,33),rep(1,7)) #5317
SC[35,]=c(rep(0,35),rep(1,5)) #5323
SC[36,]=c(rep(0,35),rep(1,5)) #5323
SC[37,]=c(rep(0,36),rep(1,4)) #5324
SC[38,]=c(rep(0,37),rep(1,3)) #5328
SC[39,]=c(rep(0,38),rep(1,2)) #5329
SC[40,]=rep(0,40) #pre XII After
SC[41,]=rep(0,40) #pre XII After

Sc = SC[-1,]
Sc = Sc[1:8, 1:8]
network <-igraph::graph_from_adjacency_matrix( Sc)
E(network)
plot(
  network,
  layout = layout_with_sugiyama(network),
  vertex.size = 8,
  vertex.color = adjustcolor("orange", alpha.f = 0.6),
  edge.arrow.size = 0.3,  # Smaller arrowheads
  edge.width = 1,
  asp = 0,
  edge.curved = 0.1
)



reduced_network = remove_transitive_edges(network)
V(reduced_network)$name <- paste0("A", 1:40)


plot(
  reduced_network,
  layout = layout_with_sugiyama(reduced_network),
  vertex.size = 10,
  vertex.color = adjustcolor("orange", alpha.f = 0.6),
  edge.arrow.size = 0.3,  # Smaller arrowheads
  edge.width = 1,
  asp = 0,
  edge.curved = 0.1
)
levels = distances(reduced_network)
depth = apply(levels, 2, function(x) ifelse(all(is.infinite(x)), 0, max(x[!is.infinite(x)])))
depth
C14ages = c(8033,7955,7940,7955,8092,7920,8027,7935,7980,7965,7990,7985,8090,
            8090,7931,7970,7940,7965,8082,7980,8025,8000,8050,7970,8024,8030,8035,8085,8030,7985,8137,8160,8150,8160,8210,8240,8199,8155,8390,8195)
C14agesEr = c(39,40,30,25,98,25,37,50,55,25,25,50,55,50,38,25,25,30,37,25,25,50,40,25,35,50,30,36,30,55,36,30,50,50,30,55,36,50,90,80)

CalAges = rcarbon::calibrate(x = C14ages, errors = C14agesEr, calCurves = 'intcal20', ids = 1:40)

plot(CalAges, 6, HPD = T, credMass = 0.68,calendar = "BCAD",
     xlab = "Years BC/AD")
sum=summary(CalAges, prob = .68)
sum
l =strsplit(sum$p_0.68_BP_1, "to")
ll =sapply(l, as.numeric)
w = 1 / ((ll[1,]-ll[2,]) /2)**2

Agehat = sapply(sum$MedianBP, function(x) (1950-x)*(x<1950) + (x-1949)*(x>=1950))
Agehat

results = IsotonicRegDAG(network = reduced_network, Ahat = Agehat, w)
V(reduced_network)$estimate <- results$A
V(reduced_network)$label <- paste0("Node", V(reduced_network), "\n",round( V(reduced_network)$estimate,2))

layout = layout_with_sugiyama(reduced_network)
layout = layout$layout *2.5
plot(
  reduced_network,
  layout =  layout,
  vertex.size = 20,
  vertex.color = adjustcolor("orange", alpha.f = 0.6),
  edge.arrow.size = 0.3,  # Smaller arrowheads
  edge.width = 1,
  asp = 0,
  vertex.label.dist = 0,
  edge.curved = 0.1,
  vertex.label.degree = -pi/4)

tg = tidygraph::as_tbl_graph(reduced_network)
library(ggraph)
ggraph(tg, layout = "sugiyama") +
  geom_node_point(aes(color = estimate), size = 5) +
  geom_node_text(aes(label = paste0("Node", name, "\n", round(estimate, 2))),
                 repel = TRUE, size = 3) +  # Enable repel to reduce overlap
  scale_color_viridis_c() +
  theme_void() + theme(legend.position = "right")


plot(c(0,41), c(min(results$A), max(results$A)), type = "n", ylab = "C(Age) BP", xlab = "TopologicalOrder")
lines(1:40,results$A, col = "red")
points(1:40, results$A, pch =12, col = "red")
points(1:40, Agehat, pch =15, col = "green")



#### données Anne en bloc
G <- graph_from_adjacency_matrix(Sc, mode = "directed")
E(G)
plot(
  G,
  layout = layout_with_sugiyama(G),
  vertex.size = 11,
  vertex.color = adjustcolor("orange", alpha.f = 0.6),
  edge.arrow.size = 0.3,  # Smaller arrowheads
  edge.width = 1,
  asp = 0,
  edge.curved = 0.1
)

reduced_G <- remove_transitive_edges(G)
plot(
  reduced_G,
  layout = layout_with_sugiyama(reduced_G),
  vertex.size = 11,
  vertex.color = adjustcolor("orange", alpha.f = 0.6),
  edge.arrow.size = 0.3,  # Smaller arrowheads
  edge.width = 1,
  asp = 0,
  edge.curved = 0.1
)
E(reduced_G)

Approx = AgeApprox(DtMeasures)

IsoBloc <- IsotonicRegDAG(reduced_G, Approx$Ahat, (1/Approx$sdAhat**2))

plot(c(1,13), c(min(Approx$Ahat)*.5, max(Approx$Ahat)*1.5), type = "n",
     xlab = "samples", ylab = "Ages")
points(1:13, IsoBloc$A, pch = 12, col = "red")
points(1:13, Approx$Ahat, pch=13, col = "green")

### representation with VisNetwork
library(visNetwork)

layout <- igraph::layout_with_sugiyama(reduced_network)$layout
nodes <- data.frame(
  label =1:5,
  stringsAsFactors = FALSE,
  x = layout[, 1] * 100,
  y = -layout[, 2] * 100   # invert Y for visNetwork
)
edges = igraph::as_data_frame(reduced_network, what = "edges")
names(edges)[1:2] <- c("from", "to")


visNetwork(nodes, edges) %>%
  visNodes(shape = "dot", size = 10) %>%
  visEdges(arrows = "to") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visLayout(randomSeed = 123)


library(htmlwidgets)
# Create a folder to store HTML frames (if recording)
dir.create("frames", showWarnings = FALSE)
edges_all <- rlang::duplicate(edges)
nodes_all <- rlang::duplicate(nodes)
# Reveal edges one-by-one
for (i in 1:nrow(edges_all)) {

  edges <- edges_all[1:i, ]
  node_ids <- unique(c(edges$from, edges$to))
  nodes <- nodes_all[nodes_all$id %in% node_ids, ]

  vis <- visNetwork(nodes, edges) %>%
    visNodes(shape = "dot", size = 20) %>%
    visEdges(arrows = "to") %>%
    visOptions(highlightNearest = TRUE) %>%
    visLayout(randomSeed = 123)

  # Save each frame
  saveWidget(vis, file = sprintf("frames/frame_%02d.html", i), selfcontained = TRUE)

  Sys.sleep(0.1)  # pause for preview (optional)
}


#### Example for the paper ####
Sc = rbind(c(0,1,0,1,1), c(0,0,0,1,1), c(0,1,0,1,1), rep(0,5), rep(0,5))
network <-graph_from_adjacency_matrix( Sc)
E(network)
plot(
  network,
  layout = layout_with_sugiyama(network),
  vertex.label = paste0("A", 1:5),
  vertex.size = 20,
  vertex.color = adjustcolor("lightblue", alpha.f = 0.6),
  edge.arrow.size = 0.4,  # Smaller arrowheads
  edge.width = 2,
  asp = 0,
  edge.curved = 0.1
)

reduced_network = remove_transitive_edges(network)
plot(
  reduced_network,
  layout = layout_with_sugiyama(reduced_network),
  vertex.label = paste0("A", 1:5),
  vertex.size = 20,
  vertex.color = adjustcolor("lightblue", alpha.f = 0.6),
  edge.arrow.size = 0.4,  # Smaller arrowheads
  edge.width = 2,
  asp = 0,
  edge.curved = 0.1
)

layout = igraph::layout_with_sugiyama(reduced_network)$layout
edges = igraph::as_data_frame(reduced_network, what = "edges")

node = data.frame(id = 1:5, label = 1:5, title = paste0("A", 1:5), x = layout[,1],
                  y= -layout[,2])
node
visNetwork::visNetwork(node, edges)  %>%
  visNetwork::visNodes(shape = "dot", size = 20) %>%
  visNetwork::visEdges(arrows = "to")

#### Clustering ####
opt <- igraph::cluster_optimal(reduced_network)
plot(opt, reduced_network, layout = igraph::layout_with_sugiyama(reduced_network))
hc_bet <- cluster_edge_betweenness(reduced_network)
plot(hc_bet, reduced_network, layout = layout_with_sugiyama(reduced_network))




Sc = matrix(0, 6,6)
Sc[1, ] = c(rep(0,2), rep(1,2), rep(0,2))
Sc[2, ] = c(rep(0,2), 1,0, rep(1,2))
Sc[5,]= c(rep(0,5),1)


G = igraph::graph_from_adjacency_matrix(Sc)
plot(G, layout = igraph::layout_with_sugiyama(G))
reduced_G = remove_transitive_edges(G)

plot(
  reduced_G,
  layout = layout_with_sugiyama(reduced_G),
  vertex.size = 11,
  vertex.color = adjustcolor("orange", alpha.f = 0.6),
  edge.arrow.size = 0.3,  # Smaller arrowheads
  edge.width = 1,
  asp = 0,
  edge.curved = 0.1
)
opt = cluster_optimal(reduced_G)
plot(opt, reduced_G, layout = layout_with_sugiyama(reduced_G))
hc_bet = cluster_edge_betweenness(reduced_G)
plot(hc_bet, reduced_G, layout = layout_with_sugiyama(reduced_G))




#### Find All Topo ####

# Function to initialize graph metadata: adjacency list and in-degree
build_graph_data <- function(g) {
  n <- vcount(g)
  adj_list <- vector("list", n)
  in_deg <- rep(0, n)

  for (v in 1:n) {
    neighbors_v <- as.integer(neighbors(g, v, mode = "out"))
    adj_list[[v]] <- neighbors_v
    for (u in neighbors_v) {
      in_deg[u] <- in_deg[u] + 1
    }
  }
  list(adj_list = adj_list, in_deg = in_deg)
}

# Recursive backtracking function
find_all_topo_orders <- function(adj_list, in_deg, path, discovered, n, result) {
  any_found <- FALSE

  for (v in 1:n) {
    if (in_deg[v] == 0 && !discovered[v]) {
      # Mark node as visited
      discovered[v] <- TRUE
      path <- c(path, v)

      # Decrease in-degree of neighbors
      for (u in adj_list[[v]]) {
        in_deg[u] <- in_deg[u] - 1
      }

      # Recurse
      result <- find_all_topo_orders(adj_list, in_deg, path, discovered, n, result)

      # Backtrack
      for (u in adj_list[[v]]) {
        in_deg[u] <- in_deg[u] + 1
      }
      discovered[v] <- FALSE
      path <- path[-length(path)]

      any_found <- TRUE
    }
  }

  # If all nodes are in path, store this order
  if (!any_found && length(path) == n) {
    result[[length(result) + 1]] <- path
  }

  return(result)
}

# Wrapper function
all_topo_orders <- function(g) {
  data <- build_graph_data(g)
  n <- vcount(g)
  discovered <- rep(FALSE, n)
  path <- integer(0)
  result <- list()

  result <- find_all_topo_orders(
    adj_list = data$adj_list,
    in_deg = data$in_deg,
    path = path,
    discovered = discovered,
    n = n,
    result = result
  )

  return(result)
}

alltopo = all_topo_orders(reduced_network)






#### coding another function for the
all_neighbors = list()
for (v in igraph::V(reduced_G)) {
  nei_v =  igraph::neighbors(reduced_G, v, mode = "total")
  all_neighbors[[v]] = nei_v
}
all_neighbors

igraph::neighbors(reduced_G, 2, mode = "in")

(neighbors(reduced_network, 2, mode = "out"))
as.numeric(neighbors(reduced_network, 1, mode = "in"))

findbounds <- function(network) {
  vertices = igraph::V(network)
  n = length(vertices)

  verticesTreatement <- function(v, network) {
    neighbors_in = as.numeric(igraph::neighbors(network, v, mode = "in")) #youngers ages
    neighbors_out = as.numeric(igraph::neighbors(network, v , mode = "out")) #older ages

    if (length(neighbors_in)==0) {
      neighbors_in = 0
    }

    if (length(neighbors_out)==0) {
      neighbors_out = n+1
    }

    return(list(upper = neighbors_out+1, lower = neighbors_in+1))

  }

    all_bounds = lapply(vertices, verticesTreatement, network = network)
  all_bounds
}

hey= findbounds(reduced_network)
hey[[4]]






library(igraph)
library(visNetwork)

# Sample graph
g <- make_ring(5)

# Layout
layout <- layout_with_fr(g)

# Labels
vertices_labels <- LETTERS[1:5]

# Nodes
nodes <- data.frame(
  id = 1:5,
  label = vertices_labels,
  shape = "dot",
  size = 30,
  font = list(color = "white", size = 20, face = "bold"),
  x = layout[, 1] * 100,
  y = -layout[, 2] * 100,
  stringsAsFactors = FALSE
)

# Edges
edges <- as_data_frame(g, what = "edges")
names(edges)[1:2] <- c("from", "to")

# Visualize
visNetwork(nodes, edges) %>%
  visEdges(arrows = "to") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visLayout(randomSeed = 123) %>% visPhysics(enabled = F)









