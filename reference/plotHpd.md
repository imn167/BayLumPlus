# Compute High Posterior Density (HPD) Interval

This function compute the HPD intervals for a set of posterior law

## Usage

``` r
plotHpd(list_object, ModelNames, level = 0.95)
```

## Arguments

- list_object:

  : A list of returned objects from the used methods for example
  [`Compute_AgeS_D()`](https://imn167.github.io/BayLumPlus/reference/Compute_AgeS_D.md)

- ModelNames:

  : the names given to each method for the legend

- level:

  : Level of the desired HPD region. Default is 0.95

## Author

Imène Bouafia (LMJL)
