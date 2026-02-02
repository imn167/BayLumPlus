# Bayesian Models for Age Estimation with Priors from ModelAgePrior Dataset

This function computes the second part of Combes & Philippe (2017),
specifically the age estimation using a Bayesian model. Its behavior is
similar to other functions like
[`AgeS_Computation()`](https://imn167.github.io/BayLumPlus/reference/AgeS_Computation.md),
with the primary difference being the first parameter.

## Usage

``` r
Compute_AgeS_D(
  DATA,
  Nb_sample,
  SampleNames,
  ThetaMatrix,
  PalaeodoseObject = NULL,
  StratiConstraints = c(),
  model = NULL,
  Iter = 10000,
  burnin = 4000,
  adapt = 1000,
  t = 5,
  n.chains = 3,
  prior = "unconstrained_jeffrey",
  PriorAge = rep(c(0.01, 100), Nb_sample),
  jags_method = "rjparallel",
  autorun = F,
  quiet = F,
  roundingOfValue = 3,
  display_plots = F,
  SavePdf = FALSE,
  OutputFileName = c("MCMCplot", "summary"),
  OutputFilePath = c(""),
  SaveEstimates = FALSE,
  OutputTableName = c("DATA"),
  OutputTablePath = c(""),
  ...
)
```

## Arguments

- DATA:

  **(required)** [list](https://rdrr.io/r/base/list.html) of objects :

  - Equivalent doses `D`

  - dose rates : `ddot`

  - standard error for D : `sD` The output of the function
    [`create_MeasuresDataFrame()`](https://imn167.github.io/BayLumPlus/reference/create_MeasuresDataFrame.md),
    containing the necessary input data for computation.

- Nb_sample:

  [integer](https://rdrr.io/r/base/integer.html) number of samples

- SampleNames:

  [character](https://rdrr.io/r/base/character.html) character vector
  with sample names

- ThetaMatrix:

  [matrix](https://rdrr.io/r/base/matrix.html) or
  [character](https://rdrr.io/r/base/character.html) input of systematic
  and individual errors.

- PalaeodoseObject:

  [list](https://rdrr.io/r/base/list.html) (with default NULL) Output of
  the
  [`Palaeodose_Computation()`](https://imn167.github.io/BayLumPlus/reference/Palaeodose_Computation.md)

- StratiConstraints:

  [matrix](https://rdrr.io/r/base/matrix.html) or
  [character](https://rdrr.io/r/base/character.html) The stratigraphic
  relation between samples

- model:

  [character](https://rdrr.io/r/base/character.html) (optional) custom
  Jags model

- Iter:

  (with default): the number of iterations to run which will be used to
  assess convergence and ages (see
  [runjags::run.jags](https://rdrr.io/pkg/runjags/man/run.jags.html)).

- burnin:

  [integer](https://rdrr.io/r/base/integer.html) (with default): the
  number of iterations used to "home in" on the stationary posterior
  distribution. These are not used for assessing convergence (see
  [runjags::run.jags](https://rdrr.io/pkg/runjags/man/run.jags.html)).

- adapt:

  [integer](https://rdrr.io/r/base/integer.html) (with default): the
  number of iterations used in the adaptive phase of the simulation (see
  [runjags::run.jags](https://rdrr.io/pkg/runjags/man/run.jags.html)).

- t:

  [integer](https://rdrr.io/r/base/integer.html) (with default): 1 every
  `t` iterations of the MCMC is considered for sampling the posterior
  distribution. (for more information see
  [runjags::run.jags](https://rdrr.io/pkg/runjags/man/run.jags.html)).

- n.chains:

  [integer](https://rdrr.io/r/base/integer.html) (with default): number
  of independent chains for the model (for more information see
  [runjags::run.jags](https://rdrr.io/pkg/runjags/man/run.jags.html)).

- prior:

  [character](https://rdrr.io/r/base/character.html) : Character string
  specifying the name of one of the models available in the
  `ModelAgePrior` dataset. Use
  [`extract_Jags_model()`](https://imn167.github.io/BayLumPlus/reference/extract_Jags_model.md)
  to see all available options

- PriorAge:

  vector (with default): lower and upper bounds for age parameter of
  each sample (in ka).

- jags_method:

  (with default): select which method to use in order to call JAGS.
  jags_methods `"rjparallel"` (the default) and `"rjags"` have been
  tested. (for more information about these possibilities and others,
  see
  [runjags::run.jags](https://rdrr.io/pkg/runjags/man/run.jags.html))

- autorun:

  [logical](https://rdrr.io/r/base/logical.html) (with default): choose
  to automate JAGS processing. JAGS model will be automatically extended
  until convergence is reached (for more information see
  [runjags::autorun.jags](https://rdrr.io/pkg/runjags/man/autorun.jags.html)).

- quiet:

  [logical](https://rdrr.io/r/base/logical.html) (with default):
  enables/disables `rjags` messages

- roundingOfValue:

  [integer](https://rdrr.io/r/base/integer.html) (with default): Integer
  indicating the number of decimal places to be used, default = 3.

- display_plots:

  [logical](https://rdrr.io/r/base/logical.html) (with default):
  enable/disable MCMC and ACF plots.

- SavePdf:

  [logical](https://rdrr.io/r/base/logical.html) (with default): if TRUE
  save graphs in pdf file named `OutputFileName` in folder
  `OutputFilePath`.

- OutputFileName:

  OutputFileName [character](https://rdrr.io/r/base/character.html)
  (with default): name of the pdf file that will be generated by the
  function if `SavePdf = TRUE`; `length(OutputFileName)=2`, see **PLOT
  OUTPUT** in **Value** section for more informations.

- OutputFilePath:

  [character](https://rdrr.io/r/base/character.html) (with default):
  path to the pdf file that will be generated by the function if
  `SavePdf`=TRUE. If it is not equal to "", it must be terminated by
  "/".

- SaveEstimates:

  [logical](https://rdrr.io/r/base/logical.html) (with default): if TRUE
  save Bayes' estimates, credible interval at level 68% and 95%, the
  result of the Gelman en Rubin test of convergence and the Time Series
  SE, in a csv table named `OutputFileName` in folder `OutputFilePath`.

- OutputTableName:

  [character](https://rdrr.io/r/base/character.html) (with default):
  name of the table that will be generated by the function if
  `SaveEstimates = TRUE`.

- OutputTablePath:

  [character](https://rdrr.io/r/base/character.html) (with default):
  path to the table that will be generated by the function if
  `SaveEstimates = TRUE`. If it is not equal to "", it must be
  terminated by "/".

- ...:

## Value

**NUMERICAL OUTPUT**

**A list of type `BayLum.list` containing the following objects:**

- **Ages** : dataframe containing the Credible interval at 95% and 68%,
  the bayes mean estimator, the bayes standard deviation estimator and
  sample names.

- **Sampling**: that corresponds to a sample of the posterior
  distributions of the age (in ka);

- **prior**: category of prior used. `prior == unconstrained` if no
  stratigraphic constraints;

- **PriorAge**: stating the priors used for the age parameter (in ka);

- **StratiConstraints**: stating the stratigraphic relations between
  samples considered in the model;

- **CovarianceMatrix**: stating the covariance matrix of error used in
  the model, highlighting common errors between samples or not;

- **model**: returns the model that was used for the Bayesian modelling
  as a [character](https://rdrr.io/r/base/character.html);

- **diagnostics_plots**: List of MCMC and ACF plots for the diagnostics
  of convergence – check attributes for Gelman CV;

- **Summary**: Summary Table of the posterior's MCMC;

**PLOT OUTPUT**

1.  **MCMC trajectories**: A graph with the MCMC trajectories and
    posterior distributions of the age. On each line, the plot on the
    left represents the MCMC trajectories, and the one on the right the
    posterior distribution of the parameter.

2.  **Summary of sample age estimates**: plot credible intervals and
    Bayes estimate of each sample age on a same graph.  
    A default `plot_mode` parameter can gives either IC segment or
    density distribution if *plot_mode = "density"*, see
    [`plot_Ages()`](https://imn167.github.io/BayLumPlus/reference/plot_Ages.md)

## Details

\*\* Which prior to use regarding the Stratigraphic constraints \*\* If
there is a strict order as a stratigraphic constraints, the user would
be able to use the following priors :

1.  **constrained_Jeffrey** : uniform order over the period of study;

2.  **old_BayLum** : old BayLum model : false configuration of the
    approximated Jeffrey;

3.  **StrictNicholls&Jones** : The Nicholls & Jones uniform order
    applied on ages;

4.  **unconstrained_jeffrey** : The approached Jeffrey (see Combes &
    Philippe 2017) without stratigraphic constraints;

## See also

[`AgeS_Computation()`](https://imn167.github.io/BayLumPlus/reference/AgeS_Computation.md)

## Examples
