# Changelog

## BayLumPlus 1.0.0

**BayLumPlus** is a refined fork of the original *BayLum* R package.
This update introduces several new capabilities, particularly for the
Age processing, where we aim to test different prior assumptions without
incurring the computational cost associated with the Palaeodose Model.

### Major New Features

- Graph-based Stratigraphic Tools: New functions to build, visualize,
  and simplify stratigraphic relationships using graph theory.

- Advanced Bayesian Age Modeling: A new function, Compute_AgeS_D(), for
  more flexible age estimation with multiple prior structures.

- Isotonic Distortion Framework: A novel Bayesian modeling strategy for
  structured monotonic trends in ages.

- Enhanced Visualizations: New plotting functions, including plotHpd(),
  to visually compare Highest Posterior Density (HPD) intervals across
  different model settings.

### New Functions

`buildNetwork()`: Builds a network from stratigraphic input.

[`network_vizualization()`](https://imn167.github.io/BayLumPlus/reference/network_vizualization.md):
Visualizes stratigraphic constraints as a graph.

[`remove_transitive_edges()`](https://imn167.github.io/BayLumPlus/reference/remove_transitive_edges.md):
Removes redundant edges implied by transitivity.

[`Compute_AgeS_D()`](https://imn167.github.io/BayLumPlus/reference/Compute_AgeS_D.md):
A new Bayesian age modeling function that supports StrictOrder,
StrictNicholls, and Independence priors.

[`IsotonicCurve()`](https://imn167.github.io/BayLumPlus/reference/IsotonicCurve.md):
Fits the isotonic distortion model.

[`PlotIsotonicCurve()`](https://imn167.github.io/BayLumPlus/reference/PlotIsotonicCurve.md):
Visualizes the results of the isotonic distortion model.

[`plotHpd()`](https://imn167.github.io/BayLumPlus/reference/plotHpd.md):
Compares Highest Posterior Density intervals for different model
settings.

### Breaking Changes from BayLum

Some function names have been updated for clarity and to better reflect
their functionality within the new architecture.

The
[`Compute_AgeS_D()`](https://imn167.github.io/BayLumPlus/reference/Compute_AgeS_D.md)
function adds two new outputs (standard deviation and time-series
standard error), which may affect scripts expecting the original output
format.

### Migration from BayLum

Users migrating from BayLum should:

- Review the function documentation for API changes, especially for
  functions related to age and palaeodose computation.

- Update code to use the new
  [`Compute_AgeS_D()`](https://imn167.github.io/BayLumPlus/reference/Compute_AgeS_D.md)
  function for advanced age modeling.

- Familiarize yourself with the new graph theory and isotonic distortion
  functions to leverage the package’s full capabilities.
