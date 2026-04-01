# JAGS Models for OSL Age Estimation in [`Compute_AgeS_DC14`](https://imn167.github.io/BayLumPlus/reference/Compute_AgeS_DC14.md)

JAGS models used to estimate true ages based on data obtained from OSL
measures and C14 measures. Both kind of measurements are used to
computed a common stratigraphy.

## Usage

``` r
Model_DC14
```

## Format

Unconstrained : Model with log-uniform settings for OSL, plain uniform
for C14 and UTh.

Constrained : Model with log-uniform order settings for OSL and plain
uniform order for C14.

NichollsJones : Nicholls&Jones' Age Model applied on ages directly

## Details

The models are designed to refine age estimation by integrating these
measurements into a Bayesian framework.

## References

To cite this package, please use: citation("BayLumPlus")
