# The Edge Pruner Algorithm (EPA). This algorithm handles the redundancy of the constraints matrix (see [`create_ThetaMatrix()`](https://imn167.github.io/BayLumPlus/reference/create_ThetaMatrix.md)) by deleting all unecessary edges. The result is

The Edge Pruner Algorithm (EPA). This algorithm handles the redundancy
of the constraints matrix (see
[`create_ThetaMatrix()`](https://imn167.github.io/BayLumPlus/reference/create_ThetaMatrix.md))
by deleting all unecessary edges. The result is

## Usage

``` r
remove_transitive_edges(Sc, n_samples = NaN)
```

## Arguments

- Sc:

  [matrix](https://rdrr.io/r/base/matrix.html) or
  [character](https://rdrr.io/r/base/character.html) : The stratigraphic
  relation between samples.

- n_samples:

  [numeric](https://rdrr.io/r/base/numeric.html) NaN by default.
  Optional size of samples to create a strict order DAG without
  constructing `Sc`

## Value

igraph a non transitive DAG that can also be referred to as a Hasse
Diagram.

## Author

Imène Bouafia

## Examples

``` r
#Example for strict order
network <- remove_transitive_edges(c(), n_samples = 5)
print(network)
#> IGRAPH f8f679b D--- 5 4 -- 
#> + edges from f8f679b:
#> [1] 1->2 2->3 3->4 4->5
```
