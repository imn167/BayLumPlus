# JAGS Models for OSL Age Estimation in [`Compute_AgeS_D`](https://imn167.github.io/BayLumPlus/reference/Compute_AgeS_D.md)

JAGS models used to estimate true OSL ages based on data obtained from
the Bayesian OSL analysis performed by the function
[`Palaeodose_Computation`](https://imn167.github.io/BayLumPlus/reference/Palaeodose_Computation.md).
Note that *unconstrained*, *uniform_order*, *Nicholls&Jones* are model
base on the likelihood \\M_i \sim \mathcal{N}(A_i, \sigma_i^2)\\ which
is different from the *OSL* model.

## Usage

``` r
ModelAgePrior
```

## Format

oldBaylum : BayLum's Old Age Model (wrong vector law)

constrained_Jeffrey : Jeffrey's Age Model with log-uniform order
settings

StrcitNicholls&Jones : Nicholls&Jones' Age Model applied on ages
directly

unconstrained_jeffrey : unconstrained Age Model (no stratigraphic
constraints)

unconstrained : Age uniformly distributed over the period of study

uniform_order : Age uniformly ordered over the period of study

Nicholls&Jones : Nicholls&Jones prior applied directly to Ages. The
duration \\A_n-A_1 \sim \mathcal{U}(0, T_2-T1)\\

## Details

These models take as input the estimated dose response (D) from
[`Palaeodose_Computation`](https://imn167.github.io/BayLumPlus/reference/Palaeodose_Computation.md)
along with the structured data matrix computed by
[`create_MeasuresDataFrame`](https://imn167.github.io/BayLumPlus/reference/create_MeasuresDataFrame.md).
The models are designed to refine age estimation by integrating these
measurements into a Bayesian framework.

## References

To cite this package, please use: citation("BayLumPlus")
