# Plot a Graph Network according to sugiyama layout

Plot a Graph Network according to sugiyama layout

## Usage

``` r
network_vizualization(network, vertices_labels, interactive = TRUE)
```

## Arguments

- network:

  igraph An igraph object

- vertices_labels:

  [character](https://rdrr.io/r/base/character.html) Samples Names for
  each node

- interactive:

  [logical](https://rdrr.io/r/base/logical.html) TRUE by default, An
  optional parameter whether to use the interactive visualization with
  VisNetwork

## Value

1.  If *interactive = TRUE* : `visNetwork` plot

2.  If *interactive = FALSE* : `ggraph` plot

## Examples

``` r
if (FALSE) { # \dontrun{
network = remove_transitive_edges(c(), n_samples = 5)
network_vizualization(network, paste0("A", 1:5))
} # }
```
