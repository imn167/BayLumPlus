# BayLumPlus ![](reference/figures/logo_BayLum.png)

An R package for chronological **Bay**esian models integrated for
Optically Stimulated (OSL) **Lum**inescence Dating

**BayLumPlus** is a a refined fork of the original
[`BayLum`](https://crp2a.github.io/BayLum/) R package. This update
introduces several new capabilities, particularly for the Age
processing, where we aim to test different prior assumptions without
incurring the computational cost associated with the Palaeodose Model.

``` R
To cite the R package 'BayLumPlus' please cite the R package itself and the following article:

  Bouafia I, Christophe C, Philippe A, Kreutzer S, Guérin G, Baumgarten F, Frerebeau N (2024). _BayLumPlus: Chronological
  Bayesian Models Integrating Optically Stimulated Luminescence and Radiocarbon Age Dating_. R package version 1.0.0,
  <https://imn167.github.io/BayLumPlus/>.

  Philippe A, Guerin G, Kreutzer S (2019). "BayLum - An R package for Bayesian analysis of OSL ages: An introduction."
  _Quaternary Geochronology_, *49*, 16-24. doi:10.1016/j.quageo.2018.05.009 <https://doi.org/10.1016/j.quageo.2018.05.009>.
```

## Installation

**You need to have [JAGS](https://mcmc-jags.sourceforge.io) installed on
your computer.**

The package **BayLumPlus** has only a development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("imn167/BayLumPLus")
```

Please note that development versions may change day by day.

## About this Extension

### Graph Theory for Stratigraphic Constraints

To support the modeling of stratigraphic relationships, we introduce
graph-based functions:

- [`network_vizualization()`](https://imn167.github.io/BayLumPlus/reference/network_vizualization.md)
  – visualizes stratigraphic constraints as a graph  
- [`remove_transitive_edges()`](https://imn167.github.io/BayLumPlus/reference/remove_transitive_edges.md)
  – called Edge Pruner Algorithm (EPA), it removes redundant edges
  implied by transitivity

These tools help simplify and explore complex stratigraphic
relationships more intuitively.

------------------------------------------------------------------------

### Bayesian Age Modeling with `Compute_AgeS_D()`

We introduce a new modeling function:
[`Compute_AgeS_D()`](https://imn167.github.io/BayLumPlus/reference/Compute_AgeS_D.md),
focused on Bayesian age estimation under various prior structures.

**Supported priors include:**

- Bayesian prior for OSL datasets:

  - `constrained_Jeffrey`: Uniform order on the log-scale (chain-like
    constraints)  
  - `StrictNicholls&Jones`: Based on the Uniform Order prior, from the
    original `BayLum`  
  - `unconstrained_Jeffrey`: For unstructured or weakly constrained
    stratigraphy

- Bayesian prior for the simple approach where the likelihood is
  $M_{i} \sim \mathcal{N}\left( A_{i},\sigma_{i}^{2} \right)\quad\forall i$:

  - `unconstrained`
  - `uniform_order`
  - `Nicholls&Jones`

*All priors for age processing are stored in the data `ModelAgePrior` to
get an easy access please use `extract_Jags_Model()`interactive
function.*

### Isotonic Distortion Framework

`BayLumPlus` introduces a new Bayesian modeling strategy called
**Isotonic Distortion**, implemented via:

- [`IsotonicCurve()`](https://imn167.github.io/BayLumPlus/reference/IsotonicCurve.md)
  – fits the isotonic model  
- [`PlotIsotonicCurve()`](https://imn167.github.io/BayLumPlus/reference/PlotIsotonicCurve.md)
  – visualizes the results with ggplots

This framework supports partial order constraints in ages, offering a
flexible and interpretable alternative to traditional priors.

------------------------------------------------------------------------

### Comparing Priors with `plotHpd()`

To facilitate the comparison of different modeling choices and prior
structures, the
[`plotHpd()`](https://imn167.github.io/BayLumPlus/reference/plotHpd.md)
function allows visual comparison of **Highest Posterior Density (HPD)**
intervals under different model settings.

This is particularly useful for:

- Sensitivity analysis  
- Model comparison  
- Reporting credible intervals with varying assumptions

## License

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or any later
version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

## Acknowledgements

The development of **BayLumPlus** received a european financial support
by the European Research Center ERC through the grant QuinaWorld.

The development of **BayLum** received a state financial support managed
by the Agence Nationale de la Recherche (France) through the program
*Investissements d’avenir* (ref. ANR-10-LABX-52) between 2015 and 2018.
