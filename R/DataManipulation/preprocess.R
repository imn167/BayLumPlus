### data used in paper
SamplesNames <- c("OSL8", "OSL7" ,"OSL6", "OSL4", "OSL1")
ddot = c( 1.40, 1.30, 1.25, 1.33, 1.13)
sddot <- rep(.03, 4)
sdc <- c(0.08, 0.09, 0.09, 0.09, 0.10)
alpha <- rep(1, 4)

DATA4 <-  combine_DataFiles(DATA2, DATA3 , DATA1)
DATA4$SampleNames <-  c( "GDB5", "FER1","GDB3")

P <- Palaeodose_Computation(DATA4, DATA4$SampleNames, DATA4$Nb_sample)


DtMeasures <- create_MeasuresDataFrame(P, DATA4,alpha[1], sdc[1:3])
DtMeasures$Measures$D
Sc = matrix(c(rep(1, 3), c(0,1,1), c(0,0,1), rep(0,3)), nrow = 4, byrow = T) #change it !!

AgeAsBayLum <-Compute_AgeS_D(DtMeasures, Sc, ModelAgePrior$Jeffreys)

AgeCorrected <- Compute_AgeS_D(DtMeasures, Sc, prior = "StrictOrder")
AgeAsBayLum$Ages
AgeCorrected$Ages


