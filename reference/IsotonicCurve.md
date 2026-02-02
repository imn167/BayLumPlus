# Compute Isotonic Regression for different stratigraphic constraints

This function compute the isotonic distorsion of a posterior
distrubution obtained by
[`Compute_AgeS_D()`](https://imn167.github.io/BayLumPlus/reference/Compute_AgeS_D.md).  
The efficient algorithm considered for the Isotnic Regression (IR) is
the **Sequetial Block Merging (SBM)**.  
If the `StratiConstraints` input is a null matrix then the model suppose
a strict order

## Usage

``` r
IsotonicCurve(
  StratiConstraints,
  object,
  interactive,
  levels = c(0.68, 0.95),
  path = tempdir(),
  rounding_digits = 3
)
```

## Arguments

- StratiConstraints:

  [matrix](https://rdrr.io/r/base/matrix.html) or
  [character](https://rdrr.io/r/base/character.html) : The stratigraphic
  relation between samples.

- object:

  [list](https://rdrr.io/r/base/list.html): output of the function
  [`Compute_AgeS_D()`](https://imn167.github.io/BayLumPlus/reference/Compute_AgeS_D.md)
  when using the prior **no_strat**. If the `StratiConstraints` input is
  a null matrix then the model suppose a strict order

- interactive:

  [logical](https://rdrr.io/r/base/logical.html) (Default) Indicating
  wether we want an interactive html file or a plot image for the DAG
  structure

- levels:

  [numeric](https://rdrr.io/r/base/numeric.html) c(0.95, 0.68) by
  default for the level of High Posterior Densities (HPD) regions. If
  the HPD region is composed of more than 1 interval then the model
  return the Credible Interval at the level indicated.

- path:

  [character](https://rdrr.io/r/base/character.html) (Default) path for
  saving the `StratiConstraints`'s DAG

- rounding_digits:

  [integer](https://rdrr.io/r/base/integer.html) (Default) digits for
  rounding estimated values.

## Value

**NUMERICAL OUTPUT : A list of type `BayLum.list` containing the
following objects**

1.  **Sampling** : Samples from the posterior distribution after
    distorsion by the isotonic regression;

2.  **network** : The DAG constructed from the `StratiConstraints`
    matrix

3.  **Ages** : data frame containing the Credible interval at 95% and
    68% , the bayes mean estimator, the bayes standard deviation
    estimator, the MAP and sample names.

## See also

[`Compute_AgeS_D()`](https://imn167.github.io/BayLumPlus/reference/Compute_AgeS_D.md)
[`PlotIsotonicCurve()`](https://imn167.github.io/BayLumPlus/reference/PlotIsotonicCurve.md)
