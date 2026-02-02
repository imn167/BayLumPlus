# Plot graphs build to vizualise isotonic distorsion of the function [`IsotonicCurve()`](https://imn167.github.io/BayLumPlus/reference/IsotonicCurve.md). The user needs to note that only one HPD regions / Credible Interval level can be allowed and vizualized.

Plot graphs build to vizualise isotonic distorsion of the function
[`IsotonicCurve()`](https://imn167.github.io/BayLumPlus/reference/IsotonicCurve.md).
The user needs to note that only one HPD regions / Credible Interval
level can be allowed and vizualized.

## Usage

``` r
PlotIsotonicCurve(object, level = 0.95)
```

## Arguments

- object:

  [list](https://rdrr.io/r/base/list.html): output of the function
  [`IsotonicCurve()`](https://imn167.github.io/BayLumPlus/reference/IsotonicCurve.md)
  when using the prior **no_strat**.

- level:

  [numeric](https://rdrr.io/r/base/numeric.html) 0.95 by default. One of
  the computed level given in
  [`IsotonicCurve()`](https://imn167.github.io/BayLumPlus/reference/IsotonicCurve.md)
  Notice that only one level can be considered.

## Value

\*\* NUMERICAL OUTPUT : A list of type `BayLum.list` containing the
following objects\*\*

1.  **ribbon** : Ribbon of the Credible Interval and solid orange line
    representing the Bayes mean estimator;

2.  **dag** : The DAG constructed from the `StratiConstraints` matrix
    with segments length represent IC(level) and gradient color the
    value of ages for each samples

3.  **Iso** : Data frame used for generating curve and DAG plots. This
    data frame is extracted from the **Ages** element of the list object
    returned by
    [`IsotonicCurve()`](https://imn167.github.io/BayLumPlus/reference/IsotonicCurve.md)
