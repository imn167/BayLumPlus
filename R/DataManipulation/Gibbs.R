####
GibbsDensity <- function(DataMeasures, A, Ai, index) {
    ## Using the function create_MeasuresDataFrame
  Measures <- DataMeasures$Measures
  D = Measures$D
  Ds = Measures$Ds
  ddot = Measures$ddot
  sddot = Measures$sddot

  ## getting the modidied
  A[index] = Ai

  detfunc <- function(A) {
    1 + sum((A*ddot)**2)
  }
}
