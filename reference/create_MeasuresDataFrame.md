# Create Input for Age Computation

This function prepares the necessary data frame entries for the function
[`Compute_AgeS_D()`](https://imn167.github.io/BayLumPlus/reference/Compute_AgeS_D.md).
It takes structured input data (as defined in
[`Generate_DataFile()`](https://imn167.github.io/BayLumPlus/reference/Generate_DataFile-deprecated.md)
and
[`Generate_DataFile_MG()`](https://imn167.github.io/BayLumPlus/reference/Generate_DataFile_MG-deprecated.md))
and computes relevant parameters for age estimation.

## Usage

``` r
create_MeasuresDataFrame(
  PalaeodoseObject,
  DATA,
  symetric_error,
  contamination_degree
)
```

## Arguments

- symetric_error:

  The \\\alpha\\ parameter in Combes & Philippe (2017)

- contamination_degree:

  The common error \\\sigma\\ in Combes & Philippe (2017)

- Data:

  A structured dataset following the format required for OSL sample
  analysis.

- P:

  A Palaeodose_computation object

## Value

A list containing:

- **Sigma**: The covariance matrix computed as: \$\$\Theta = A \Sigma
  A + \text{diag}(sD)\$\$

- **Info**: A list with details including:

  - Number of samples (`Nb_sample`)

  - Sample names (`NamesOfSamples`)

  - Dose rate values (`ddot`)

  - Dose rate uncertainties (`sddot`)

  - Estimated palaeodose values (`D`)

  - Estimated palaeodose uncertainties (`sD`)

## See also

`Computation_AgeS_D`,
[`Palaeodose_Computation`](https://imn167.github.io/BayLumPlus/reference/Palaeodose_Computation.md),
[`Generate_DataFile`](https://imn167.github.io/BayLumPlus/reference/Generate_DataFile-deprecated.md),
[`Generate_DataFile_MG`](https://imn167.github.io/BayLumPlus/reference/Generate_DataFile_MG-deprecated.md)
