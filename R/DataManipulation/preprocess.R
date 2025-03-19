### data used in paper
SamplesNames <- c("OSL8", "OSL7" ,"OSL6", "OSL4", "OSL1")
ddot = c( 1.40, 1.30, 1.25, 1.33, 1.13)
sddot <- rep(.03, 4)
sdc <- c(0.08, 0.09, 0.09, 0.09, 0.10)
alpha <- rep(1, 4)

DATA4 <-  combine_DataFiles(DATA1, DATA2, DATA3)
DATA4$SampleNames <-  c("GDB3", "GDB5", "FER1")

P <- Palaeodose_Computation(DATA4, DATA4$SampleNames, DATA4$Nb_sample)

DATA4$ddot_env
DtMeasures <- create_MeasuresDataFrame(P, DATA4,alpha[1], sdc[1:3])
DtMeasures$Measures$D
Sc = matrix(c(rep(1, 3), rep(0,3), c(1,0,1), c(1,0,0)), nrow = 4, byrow = T)

AgeAsBayLum <-Compute_AgeS_D(DtMeasures, Sc, ModelAgePrior)
