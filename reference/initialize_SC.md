# Initialization of the Gibbs Algorithm according to the fixed stratigraphy

In this function we try different way of initializing our Gibbs Sampler.
Several limitation can be encoutered :

- The constraints need to be respected

- The ages need to belongs to a certain period of study (T1, T2) that
  will be fixed arbitrarily

- Initialisation are not only constrained by the Matrix of constraint
  but also to avoid numerical zeros the interval of initialisation need
  to be reasonable according to the density of Gibbs Sampler.

The most efficient approach would be to considered the approximation
Ahat = Dhat / ddot and then make it respect the constraints

## Usage

``` r
initialize_SC(Sc, LowerPeriod, UpperPeriod, type, plotGraph = F, ...)
```
