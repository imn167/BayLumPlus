# Create Age Plot

Create Age Plot

## Usage

``` r
plot_Ages(
  object,
  sample_names = NULL,
  sample_order = NULL,
  plot_mode = "ages",
  model = "Baylum",
  ...
)
```

## Arguments

- object:

  [list](https://rdrr.io/r/base/list.html) or
  [data.frame](https://rdrr.io/r/base/data.frame.html) (**required**):
  Output as created by functions like
  [AgeC14_Computation](https://imn167.github.io/BayLumPlus/reference/AgeC14_Computation.md),
  which is a list of class `BayLum.list`. Alternatively the function
  supports a [data.frame](https://rdrr.io/r/base/data.frame.html) as
  input, however, in such a case the
  [data.frame](https://rdrr.io/r/base/data.frame.html) must resemble the
  ages [data.frame](https://rdrr.io/r/base/data.frame.html) created by
  the computation functions otherwise the input will be silently
  ignored.

- sample_names:

  [character](https://rdrr.io/r/base/character.html) (optional):
  alternative sample names used for the plotting. If the length of the
  provided [character](https://rdrr.io/r/base/character.html) vector is
  shorter than the real number of samples, the names are recycled.

- sample_order:

  [numeric](https://rdrr.io/r/base/numeric.html) (optional): argument to
  rearrange the sample order, e.g., `sample_order = c(4:1)` plots the
  last sample first.

- plot_mode:

  [character](https://rdrr.io/r/base/character.html) (*with default*):
  allows to switch from displaying ages as points with lines (`"ages"`)
  for the credible intervals to densities (`"density"`)

- model:

  [character](https://rdrr.io/r/base/character.html) "BayLum" by
  default. Name of the used Bayesian Age model (see `ModelAgePrior`).

- ...:

  further arguments to control the plot output, standard arguments are:
  `cex`, `xlim`, `main`, `xlab`, `col` further (non-standard) arguments
  are: `grid` (`TRUE`/`FALSE`), `legend` (`TRUE`/`FALSE`), `legend.text`
  ([character](https://rdrr.io/r/base/character.html) input needed),
  `legend.pos`
  [graphics::legend](https://rdrr.io/r/graphics/legend.html),
  `legend.cex`. Additional arguments: `d_scale` (scales density plots),
  `show_ages` (add ages to density plots)

## Value

The function returns a plot and the
[data.frame](https://rdrr.io/r/base/data.frame.html) used to display the
data

## Details

This function creates an age plot showing the mean ages along with the
credible intervals. The function provides various arguments to modify
the plot output, however, for an ultimate control the function returns
the [data.frame](https://rdrr.io/r/base/data.frame.html) extracted from
the input object for own plots.

## Function version

0.1.5

## See also

[AgeC14_Computation](https://imn167.github.io/BayLumPlus/reference/AgeC14_Computation.md),
[AgeS_Computation](https://imn167.github.io/BayLumPlus/reference/AgeS_Computation.md)

## Author

Sebastian Kreutzer, Institute of Geography, Ruprecht-Karl-University of
Heidelberg (Germany), based on code written by Claire Christophe

## Examples

``` r
## load data
data(DATA_C14,envir = environment())
C14Cal <- DATA_C14$C14[,1]
SigmaC14Cal <- DATA_C14$C14[,2]
Names <- DATA_C14$Names
nb_sample <- length(Names)

## Age computation
Age <- AgeC14_Computation(
   Data_C14Cal = C14Cal,
   Data_SigmaC14Cal = SigmaC14Cal,
   SampleNames = Names,
   Nb_sample = nb_sample,
   PriorAge = rep(c(20,60),nb_sample),
   Iter = 500,
   quiet = TRUE)
#> Finished running the simulation
#> Warning: [plot_MCMC()] Argument 'n.iter' out of range, reset to number of observations in mcmc.list / mcmc




#> 
#> 
#> >> MCMC Convergence of Age parameters <<
#> ----------------------------------------------
#> Sample name   Bayes estimate   Uppers credible interval
#> A_S-EVA-26510     1       1 
#> A_S-EVA-26506     1.009       1.023 
#> A_S-EVA-26507     1.002       1.007 
#> A_S-EVA-26508     1.003       1.01 
#> 
#> 
#> ---------------------------------------------------------------------------------------------------
#>  *** WARNING: The following information are only valid if the MCMC chains have converged  ***
#> ---------------------------------------------------------------------------------------------------
#> 
#> 
#> 
#> >> Bayes estimates of Age for each sample and credible interval <<
#> ------------------------------------------------------
#> Sample name   Bayes estimate  Credible interval: 
#> A_S-EVA-26510     41.9614143424911 
#>                       lower bound     upper bound
#>               at level 95%    41.473          42.356 
#>               at level 68%    41.827          42.24 
#> ------------------------------------------------------
#> Sample name   Bayes estimate  Credible interval: 
#> A_S-EVA-26506     45.712113182689 
#>                       lower bound     upper bound
#>               at level 95%    45.035          46.207 
#>               at level 68%    45.492          46.024 
#> ------------------------------------------------------
#> Sample name   Bayes estimate  Credible interval: 
#> A_S-EVA-26507     44.8905482327229 
#>                       lower bound     upper bound
#>               at level 95%    43.719          45.866 
#>               at level 68%    44.455          45.394 
#> ------------------------------------------------------
#> Sample name   Bayes estimate  Credible interval: 
#> A_S-EVA-26508     45.0256760593151 
#>                       lower bound     upper bound
#>               at level 95%    44.037          46.22 
#>               at level 68%    44.577          45.576 
#> 
#> ------------------------------------------------------



## plot output
plot_Ages(Age)

#>        SAMPLE      AGE HPD68.MIN HPD68.MAX HPD95.MIN HPD95.MAX ALT_SAMPLE_NAME
#> 1 S-EVA-26510 41.96141    41.827    42.240    41.473    42.356              NA
#> 2 S-EVA-26506 45.71211    45.492    46.024    45.035    46.207              NA
#> 3 S-EVA-26507 44.89055    44.455    45.394    43.719    45.866              NA
#> 4 S-EVA-26508 45.02568    44.577    45.576    44.037    46.220              NA
#>   AT
#> 1  4
#> 2  3
#> 3  2
#> 4  1

## plot output
plot_Ages(Age, plot_mode = "density", legend.pos = "topright")

#>        SAMPLE      AGE HPD68.MIN HPD68.MAX HPD95.MIN HPD95.MAX ALT_SAMPLE_NAME
#> 1 S-EVA-26510 41.96141    41.827    42.240    41.473    42.356              NA
#> 2 S-EVA-26506 45.71211    45.492    46.024    45.035    46.207              NA
#> 3 S-EVA-26507 44.89055    44.455    45.394    43.719    45.866              NA
#> 4 S-EVA-26508 45.02568    44.577    45.576    44.037    46.220              NA
#>   AT
#> 1  4
#> 2  3
#> 3  2
#> 4  1
```
