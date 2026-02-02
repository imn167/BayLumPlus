# Chronological Bayesian Models Integrating Optically Stimulated Luminescence and C-14 Dating

A collection of various R functions for Bayesian analysis of
luminescence data and C-14 age estimates. This includes, amongst others,
data import, export, application of age and palaeodose models.

## Details

This package is based on the functions:
[`Generate_DataFile()`](https://imn167.github.io/BayLumPlus/reference/Generate_DataFile-deprecated.md)
and
[`Generate_DataFile_MG()`](https://imn167.github.io/BayLumPlus/reference/Generate_DataFile_MG-deprecated.md)
to import luminescence data. These functions create a list containing
all informations to compute age of single-grain OSL measurements for the
first function and multi-grain OSL measurements for the second.

The functions:
[`Age_Computation()`](https://imn167.github.io/BayLumPlus/reference/Age_Computation.md)
and
[`AgeS_Computation()`](https://imn167.github.io/BayLumPlus/reference/AgeS_Computation.md)
use Bayesian analysis for OSL age estimation for one or various samples
according to difference models (e.g. different dose-response curves and
different equivalent dose distributions around the palaeodose).

It is possible to consider various BIN/BINX-files per sample, to compute
ages of samples in stratigraphic constraints and to integrate systematic
errors.

It is possible to calibrate C-14 age with the function
[`AgeC14_Computation()`](https://imn167.github.io/BayLumPlus/reference/AgeC14_Computation.md).
We can also estimate chronology containing 14C age and OSL samples with
the function
[`Age_OSLC14()`](https://imn167.github.io/BayLumPlus/reference/Age_OSLC14.md).

## Note

This work received a state financial support managed by the Agence
Nationale de la Recherche (France) through the program *Investissements
d'avenir* (ref. ANR-10-LABX-52).

## References

Philippe, A., Guérin, G., Kreutzer, S., 2019. BayLum - An R package for
Bayesian analysis of OSL ages: An introduction. *Quaternary
Geochronology* 49, 16-24.
[doi:10.1016/j.quageo.2018.05.009](https://doi.org/10.1016/j.quageo.2018.05.009)

## See also

Useful links:

- <https://imn167.github.io/BayLumPlus/>

- <https://github.com/imn167/BayLumPlus>

- Report bugs at <https://github.com/imn167/BayLumPlus/issues>

## Author

**Maintainer**: Anne Philippe <anne.philippe@univ-nantes.fr>
([ORCID](https://orcid.org/0000-0002-5331-5087))

Authors:

- Imene Bouafia

- Claire Christophe

- Sebastian Kreutzer ([ORCID](https://orcid.org/0000-0002-0734-2199))

- Guillaume Guérin ([ORCID](https://orcid.org/0000-0001-6298-5579))

- Frederik Harly Baumgarten
  ([ORCID](https://orcid.org/0000-0002-4374-5948))

- Nicolas Frerebeau ([ORCID](https://orcid.org/0000-0001-5759-4944))

Other contributors:

- ERC (Institutional contributor) \[copyright holder, funder\]

- CNRS (Institutional contributor) \[funder\]

- QuinaWorld (Institutional contributor) \[funder\]

- LMJL (Institutional contributor) \[copyright holder, funder\]

- Université de Rennes (Institutional contributor) \[copyright holder,
  funder\]

- Université Bordeaux Montaigne (Institutional contributor) \[copyright
  holder, funder\]

- LabEx Sciences archéologiques de Bordeaux (Institutional contributor)
  \[funder\]
