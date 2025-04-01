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


GibbsOutput = GibbsSampler(DtMeasures,  3, 50000, 30000,Sc, 1, 150, "logit")
invGibbsOutput = GibbsSampler(invDtMeasures,  invVhat, 50000, 30000,Sc, 1, 150, "logit")

plot(coda::as.mcmc.list(coda::as.mcmc(GibbsOutput$A)))

traceplot(mcmc(GibbsOutput$A))
densplot(mcmc(GibbsOutput$A))
summary(mcmc(GibbsOutput$A))
autocorr.diag(mcmc(GibbsOutput$A))
acfplot(mcmc(GibbsOutput$A))
pairs(GibbsOutput$A)
gelman.diag(mcmc(GibbsOutput$A))


###-------------------------------------------------------------------------------@
datasets1 <- readxl::read_excel("R/DataManipulation/dataOSL.xlsx")
colnames(datasets1)[3:6] <- c("ddot", "sddot", "D", "sD")
Sc = rbind(rep(1, 3), upper.tri(matrix(rep(1), ncol = 3, nrow = 3))*1)

Measures = list(ddot = datasets1$ddot, sddot = datasets1$sddot, D = datasets1$D, sD = datasets1$sD,
                sddot_shared = rep(0,3), Nb_sample = 3, SampleNames = datasets1$Sample)

Theta = diag(as.numeric(Measures$sddot)) +
  as.numeric(Measures$sddot_shared) %*% t(as.numeric(Measures$sddot_shared))
covD = diag(as.numeric(Measures$sD))

DtMeasures <- list(Theta = Theta, Measures = Measures, covD = covD)

GibbsOutput <-  GibbsSampler(DtMeasures, 3, 70000, 50000, Sc, 1, 500, "logit", proposal_magnitude = 1)
pdf("R/DataManipulation/dataOSL.pdf", width = 18)
plot(GibbsOutput$Sampling)
coda::acfplot(GibbsOutput$Sampling)
dev.off()
AgeAsBayLum <-Compute_AgeS_D(DtMeasures, Sc, ModelAgePrior$Jeffreys, Iter = 2000, burnin = 50000, t = 10,adapt = 0,
                             PriorAge = rep(c(1, 500),  Measures$Nb_sample))

plotHpd(list(GibbsOutput, AgeAsBayLum), c("Gibbs", "BayLum"))


###-------------------------------------------------------------------------------@
datasets2 <- readxl::read_excel("R/DataManipulation/dataOSL_3sites.xlsx")
colnames(datasets2)[3:6] <- c("ddot", "sddot", "D", "sD")

BlockLength <- datasets2 %>% dplyr::mutate(block = stringr::str_extract(Sample, "Dhab\\d+")) %>%
  dplyr::group_by(block) %>%
  dplyr::summarise(blockLength = dplyr::n()) %>% dplyr::select(blockLength)

Sc = as.matrix(Matrix::bdiag(sapply((BlockLength$blockLength), function(n) upper.tri(matrix(1, nrow = n, ncol = n)))))
Sc = as.matrix( rbind(rep(1, dim(datasets2)[1]), Sc))

Measures = list(ddot = datasets2$ddot, sddot = datasets2$sddot, D = datasets2$D, sD = datasets2$sD,
                sddot_shared = rep(0,dim(datasets2)[1]), Nb_sample = dim(datasets2)[1],
                SampleNames = datasets2$Sample)


Theta = diag(as.numeric(Measures$sddot)) +
  as.numeric(Measures$sddot_shared) %*% t(as.numeric(Measures$sddot_shared))

covD = diag(as.numeric(Measures$sD))

DtMeasures <- list(Theta = Theta, Measures = Measures, covD = covD)

GibbsOutput <-  GibbsSampler(DtMeasures, 3, 90000, 70000, Sc, 1, 1400, "logit", proposal_magnitude = 1,
                             plotlegend.pos = locator(1))
pdf("R/DataManipulation/blockdataset.pdf", width = 17)
plot(GibbsOutput$Sampling)
coda::acfplot(GibbsOutput$Sampling)
dev.off()
AgeAsBayLum <-Compute_AgeS_D(DtMeasures, Sc, ModelAgePrior$Jeffreys, Iter = 2000, burnin = 50000, t = 10, adapt = 0,
                             PriorAge = rep(c(1, 1400),  Measures$Nb_sample))


plotHpd(list(GibbsOutput, AgeAsBayLum), c("Gibbs", "JagsBayLum"))

## ============ Elaine data sets ======== ##

extractElaine <- function() {
    #strati
  SampleNames <- readxl::read_xlsx("R/DataManipulation/All_Qz_Grains_BayLum_doses.xlsx", sheet = 3)[, 1]
  OSLestimate <- as.matrix(readxl::read_xlsx("R/DataManipulation/All_Qz_Grains_BayLum_doses.xlsx", sheet = 3)[, c(4,7)])
  Theta <- as.matrix(readxl::read_xlsx("R/DataManipulation/All_Qz_Grains_BayLum_doses.xlsx", sheet = 5))
  Sc <- as.matrix(readxl::read_xlsx("R/DataManipulation/All_Qz_Grains_BayLum_doses.xlsx", sheet = 6))
  n = dim(Theta)[1]
  ddot =  as.numeric(readxl::read_xlsx("R/DataManipulation/All_Qz_Grains_BayLum_doses.xlsx", sheet = 4)[3,c(12,13)])

  Measures = list(SampleNames = SampleNames$Layer, D = OSLestimate[, 1], sD = OSLestimate[, 2], Nb_sample = n ,
                  ddot = rep(ddot[1], n), sddot = rep(ddot[2], n), sddot_shared = rep(0, n))

  return( list(Measures = Measures, Theta = Theta, covD = diag(Measures$sD), Sc = Sc))

}

DtMeasures <- extractElaine()
GibbsOutput <- GibbsSampler(DtMeasures,3,  90000, 70000, DtMeasures$Sc, 5, 500, "logit", proposal_magnitude = 2, order = T)
initialize_SC(DtMeasures$Sc, 5, 300)
pdf("R/DataManipulation/Elainedataset.pdf", width = 17)
plot(GibbsOutput$Sampling)
coda::acfplot(GibbsOutput$Sampling)
dev.off()

AgeAsBayLum <-Compute_AgeS_D(DtMeasures, DtMeasures$Sc, ModelAgePrior$Jeffreys, Iter = 2000, burnin = 50000, t = 10,
                             PriorAge = rep(c(1, 500),  DtMeasures$Measures$Nb_sample))

AgeCorrected <- Compute_AgeS_D(DtMeasures, DtMeasures$Sc, prior = "StrictOrder", Iter = 2000, burnin = 50000, t = 10, adapt = 0,
                               PriorAge = rep(c(1, 500),  DtMeasures$Measures$Nb_sample))




