### Data simulated
Measures <- list(SampleNames = c("OSL-001", "OSL-002", "OSL-003", "OSL-004", "OSL-005"),
ddot = c(2.50, 2.30, 2.10, 1.90, 1.75),  # Environmental dose rate (Gy/ka)
sddot = c(0.10, 0.12, 0.15, 0.10, 0.08), # Error on dose rate (Gy/ka)
D = c(75.0, 85.0, 92.0, 105.0, 115.0),  # OSL Dose (Gy)
sD = c(3.0, 3.5, 4.0, 5.0, 4.5),  # Error on dose (Gy)
sddot_shared = c(1.5, 1.6, 1.7, 1.9, 2.0),  # Common error term (Gy)
Nb_sample = 5)


### List of the Measures
Theta = diag(as.numeric(Measures$sddot)) +
  as.numeric(Measures$sddot_shared) %*% t(as.numeric(Measures$sddot_shared))
covD = diag(as.numeric(Measures$sD))

DtMeasures <- list(Theta = Theta, Measures = Measures, covD = covD)

## RW step
Ahat = as.numeric(Measures$D)/as.numeric(Measures$ddot)
vAhat = Ahat *sqrt((as.numeric(Measures$sD / Measures$D))**2 +
                     as.numeric(Measures$sddot/Measures$ddot)**2)

#### Try inverse
invVhat <- vAhat
invVhat[c(2,3)] = invVhat[c(3,2)]

invMeasures = lapply(Measures[2:6], function(l){
  l[c(2,3)] = l[c(3,2)]

  return(l)} )
invMeasures$SampleNames <- Measures$SampleNames
invMeasures$Nb_sample <- Measures$Nb_sample
invTheta = diag(as.numeric(invMeasures$sddot)) +
  as.numeric(invMeasures$sddot_shared) %*% t(as.numeric(invMeasures$sddot_shared))
invcovD = diag(as.numeric(invMeasures$sD))
invDtMeasures <- list(Theta = invTheta, Measures = invMeasures, covD = invcovD)

##-------------------------- True Ages -------------------#@
Age = c(30.0, 37.0, 43.8, 55.3, 65.7)  # Calculated age (ka)
##-------------------------- Constraints Matrix -------------------#@
Sc = rbind(rep(1, 5), upper.tri(matrix(rep(1), ncol = 5, nrow = 5))*1)

AgeAsBayLum <-Compute_AgeS_D(DtMeasures, Sc, ModelAgePrior$Jeffreys, Iter = 50000, burnin = 30000,
                             PriorAge = rep(c(1, 150),  Measures$Nb_sample))

invAgeAsBayLum <- Compute_AgeS_D(invDtMeasures, Sc, ModelAgePrior$Jeffreys, Iter = 50000, burnin = 30000,
                                 PriorAge = rep(c(1, 150),  Measures$Nb_sample))

AgeCorrected <- Compute_AgeS_D(DtMeasures, Sc, prior = "StrictOrder", Iter = 70000, burnin = 50000,
                               PriorAge = rep(c(1, 150),  Measures$Nb_sample))


init_dist <- replicate(1000, initialize_SC(Sc, 1, 150))
ldensity <- apply(init_dist, 1, function(l) plot(density(l)))


GibbsOutput = GibbsSampler(DtMeasures,  vAhat, 50000, 30000,Sc, 1, 150, "logit")
invGibbsOutput = GibbsSampler(invDtMeasures,  invVhat, 50000, 30000,Sc, 1, 150, "logit")

plot(coda::as.mcmc.list(coda::as.mcmc(GibbsOutput$A)))





###-------------------------------------------------------------------------------@
# DATA4 <-  combine_DataFiles(DATA2, DATA3 , DATA1)
# DATA4$SampleNames <-  c( "GDB5", "FER1","GDB3")
#
# P <- Palaeodose_Computation(DATA4, DATA4$SampleNames, DATA4$Nb_sample)
#
#
# DtMeasures <- create_MeasuresDataFrame(P, DATA4,alpha[1], sdc[1:3])
# DtMeasures$Measures$D
C14_SampleNames = c("OxA-21261","OxA-21262","UCIAMS-98210","UCIAMS-103134","P-782","UCIAMS-103138","OxA-23247","OxA-9774",
                    "OxA-9946","UCIAMS-98208","UCIAMS-98209",
                    "OxA-9947","OxA-9775","OxA-9948","OxA-23523","UCIAMS-103135","UCIAMS_103136","UCIAMS-103137",
                    'OxA-23248',"UCIAMS-103141","UCIAMS-103140",
                    "OxA-27087","OxA-9949","UCIAMS-103139","OxA-23249","OxA-9950","UCIAMS-109991","OxA-23250","UCIAMS-109992","OxA-9776",
                    "OxA-23251","UCIAMS-109993","OxA-9892","OxA-9777","UCIAMS-109994","OxA-9778","OxA-23252","OxA-9893","PL-980525A","AA-27982")

SC = matrix(data=0,ncol=40,nrow=41) ### matrix to account for stratigraphic constraints
SC[1,]=rep(1,40)
SC[2,]=c(rep(0,1),rep(1,1),rep(0,3),rep(1,35)) #4853
SC[3,]=c(rep(0,5),rep(1,35)) #4861
SC[4,]=c(rep(0,5),rep(1,35)) #4555
SC[5,]=c(rep(0,5),rep(1,35)) #4779
SC[6,]=rep(0,40) #P782 After
SC[7,]=c(rep(0,9),rep(1,31)) #4850
SC[8,]=c(rep(0,9),rep(1,31)) #4850
SC[9,]=rep(0,40) #4715 (After)
SC[10,]=rep(0,40) #4715 (After)
SC[11,]=c(rep(0,18),rep(1,4),rep(0,1),rep(1,17)) #4517
SC[12,]=c(rep(0,18),rep(1,4),rep(0,1),rep(1,17)) #4517
SC[13,]=rep(0,40)#4822 (After)
SC[14,]=rep(0,40) #4826 (After)
SC[15,]=rep(0,40) #4826 (After)
SC[16,]=c(rep(0,16),rep(1,2),rep(0,4),rep(1,1),rep(0,1),rep(1,16)) #4828
SC[17,]=c(rep(0,16),rep(1,2),rep(0,4),rep(1,1),rep(0,1),rep(1,16)) #4828
SC[18,]=c(rep(0,17),rep(1,1),rep(0,4),rep(1,1),rep(0,1),rep(1,16)) #4836
SC[19,]=c(rep(0,22),rep(1,1),rep(0,1),rep(1,16)) #4837
SC[20,]=c(rep(0,19),rep(1,2),rep(0,1),rep(1,18)) #4869
SC[21,]=c(rep(0,19),rep(1,2),rep(0,1),rep(1,18)) #4869
SC[22,]=c(rep(0,22),rep(1,18)) #4867
SC[23,]=c(rep(0,22),rep(1,18)) #4867
SC[24,]=c(rep(0,24),rep(1,16)) #4848
SC[25,]=c(rep(0,24),rep(1,16)) #4865
SC[26,]=c(rep(0,25),rep(1,15)) #4878
SC[27,]=c(rep(0,26),rep(1,14)) #5276
SC[28,]=c(rep(0,27),rep(1,13)) #5279
SC[29,]=c(rep(0,29),rep(1,11)) #5283
SC[30,]=c(rep(0,29),rep(1,11)) #5283
SC[31,]=c(rep(0,30),rep(1,10)) #5292
SC[32,]=c(rep(0,31),rep(1,9)) #5308
SC[33,]=c(rep(0,32),rep(1,8)) #5316
SC[34,]=c(rep(0,33),rep(1,7)) #5317
SC[35,]=c(rep(0,35),rep(1,5)) #5323
SC[36,]=c(rep(0,35),rep(1,5)) #5323
SC[37,]=c(rep(0,36),rep(1,4)) #5324
SC[38,]=c(rep(0,37),rep(1,3)) #5328
SC[39,]=c(rep(0,38),rep(1,2)) #5329
SC[40,]=rep(0,40) #pre XII (After)
SC[41,]=rep(0,40) #pre XII (After)

initialize_SC(SC, 1,100,T)



