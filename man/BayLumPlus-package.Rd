% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BayLum-package.R
\docType{package}
\name{BayLumPlus-package}
\alias{BayLumPlus}
\alias{BayLumPlus-package}
\alias{BayLum}
\title{Chronological Bayesian Models Integrating Optically Stimulated Luminescence
and C-14 Dating}
\description{
A collection of various R functions for Bayesian analysis of luminescence
data and C-14 age estimates. This includes, amongst others, data import,
export, application of age and palaeodose models.
}
\details{
This package is based on the functions: \code{\link[=Generate_DataFile]{Generate_DataFile()}} and
\code{\link[=Generate_DataFile_MG]{Generate_DataFile_MG()}} to import luminescence data. These functions create
a list containing all informations to compute age of single-grain OSL
measurements for the first function and multi-grain OSL measurements for
the second.

The functions: \code{\link[=Age_Computation]{Age_Computation()}} and \code{\link[=AgeS_Computation]{AgeS_Computation()}} use Bayesian
analysis for OSL age estimation for one or various samples according to
difference models (e.g. different dose-response curves and different
equivalent dose distributions around the palaeodose).

It is possible to consider various BIN/BINX-files per sample, to compute
ages of samples in stratigraphic constraints and to integrate systematic
errors.

It is possible to calibrate C-14 age with the function
\code{\link[=AgeC14_Computation]{AgeC14_Computation()}}. We can also estimate chronology containing 14C age
and OSL samples with the function \code{\link[=Age_OSLC14]{Age_OSLC14()}}.
}
\note{
This work received a state financial support managed by the Agence Nationale
de la Recherche (France) through the program \emph{Investissements d'avenir}
(ref. ANR-10-LABX-52).
}
\references{
Philippe, A., Guérin, G., Kreutzer, S., 2019. BayLum - An R package for
Bayesian analysis of OSL ages: An introduction. \emph{Quaternary Geochronology}
49, 16-24. \doi{10.1016/j.quageo.2018.05.009}
}
\seealso{
Useful links:
\itemize{
  \item \url{https://CRAN.r-project.org/package=BayLum}
  \item \url{https://crp2a.github.io/BayLum/}
  \item Report bugs at \url{https://github.com/crp2a/BayLum/issues}
}

}
\author{
\strong{Maintainer}: Anne Philippe \email{anne.philippe@univ-nantes.fr} (\href{https://orcid.org/0000-0002-5331-5087}{ORCID})

Authors:
\itemize{
  \item Claire Christophe
  \item Sebastian Kreutzer (\href{https://orcid.org/0000-0002-0734-2199}{ORCID})
  \item Guillaume Guérin (\href{https://orcid.org/0000-0001-6298-5579}{ORCID})
  \item Frederik Harly Baumgarten (\href{https://orcid.org/0000-0002-4374-5948}{ORCID})
  \item Nicolas Frerebeau (\href{https://orcid.org/0000-0001-5759-4944}{ORCID})
  \item Imene Bouafia
}

Other contributors:
\itemize{
  \item Université Bordeaux Montaigne (Institutional contributor) [copyright holder, funder]
  \item CNRS (Institutional contributor) [funder]
  \item LabEx Sciences archéologiques de Bordeaux (Institutional contributor) [funder]
}

}
\keyword{internal}
