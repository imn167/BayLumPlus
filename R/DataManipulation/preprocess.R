### Data simulated
Measures <- list(SampleNames = c("OSL-001", "OSL-002", "OSL-003", "OSL-004", "OSL-005"),
ddot = c(2.50, 2.30, 2.10, 1.90, 1.75),  # Environmental dose rate (Gy/ka)
sddot = c(0.10, 0.12, 0.15, 0.10, 0.08), # Error on dose rate (Gy/ka)
D = c(75.0, 85.0, 92.0, 105.0, 115.0),  # OSL Dose (Gy)
sD = c(3.0, 3.5, 4.0, 5.0, 4.5),  # Error on dose (Gy)
sddot_shared = c(1.5, 1.6, 1.7, 1.9, 2.0),  # Common error term (Gy)
Nb_sample = 5)


### List of the Measures
Theta = diag(as.numeric(Measures$sddot)**2) +
  as.numeric(Measures$sddot_shared) %*% t(as.numeric(Measures$sddot_shared))
covD = diag(as.numeric(Measures$sD)**2)

DtMeasures <- list(Theta = Theta, Measures = Measures, covD = covD)


##-------------------------- True Ages -------------------#@
Age = c(30.0, 37.0, 43.8, 55.3, 65.7)  # Calculated age (ka)
##-------------------------- Constraints Matrix -------------------#@
Sc = rbind(rep(1, 5), upper.tri(matrix(rep(1), ncol = 5, nrow = 5))*1)

AgeAsBayLum <-Compute_AgeS_D(DtMeasures, Sc, prior = "Jeffreys", Iter = 50000, burnin = 30000,
                             PriorAge = rep(c(1, 1000),  Measures$Nb_sample))

AgeCorrected <- Compute_AgeS_D(DtMeasures, Sc, prior = "StrictOrder", Iter = 2000, burnin = 50000,
                               PriorAge = rep(c(1, 1000),  Measures$Nb_sample))
AgeCorrected <- Compute_AgeS_D(DtMeasures, Sc, prior = "StrictNicholls", Iter = 2000, burnin = 50000,
                               PriorAge = rep(c(1, 1000),  Measures$Nb_sample))


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

Theta = diag(as.numeric(Measures$sddot)**2) +
  as.numeric(Measures$sddot_shared) %*% t(as.numeric(Measures$sddot_shared))
covD = diag(as.numeric(Measures$sD)**2)

DtMeasures <- list(Theta = Theta, Measures = Measures, covD = covD)

GibbsOutput <-  GibbsSampler(DtMeasures, 3, 70000, 50000, Sc, 1, 500, "logit", proposal_magnitude = 1, lambda = 1, lag = 10, n_logitStep = 5,
                             type = "other")
pdf("R/DataManipulation/dataOSL.pdf", width = 18)
plot(GibbsOutput$Sampling)
coda::acfplot(GibbsOutput$Sampling)
plot(GibbsOutput$name_chains)
coda::acfplot(GibbsOutput$name_chains)
dev.off()
AgeAsBayLum <-Compute_AgeS_D(DtMeasures, Sc, prior = "Jeffreys", Iter = 2000, burnin = 50000, t = 10,
                             PriorAge = rep(c(1, 1400),  Measures$Nb_sample))

AgeCorrected <- Compute_AgeS_D(DtMeasures, Sc, prior = "StrictOrder", Iter = 2000, burnin = 50000,t = 10,
                               PriorAge = rep(c(1, 1400),  Measures$Nb_sample))

AgeNicholls <- Compute_AgeS_D(DtMeasures, DtMeasures$Sc, prior = "StrictNicholls", Iter = 2000, burnin = 50000, t = 10,
                              PriorAge = rep(c(1, 1400),  DtMeasures$Measures$Nb_sample))


plotHpd(list(AgeNicholls, AgeAsBayLum, AgeCorrected), c("Nicholls", "BayLum", "Exponential")) +
  ggplot2::geom_point(mapping = ggplot2::aes(SAMPLE, AGE), data = AgeNicholls$Ages, inherit.aes = F, color = "blue") +
  ggplot2::geom_point(ggplot2::aes(SAMPLE, AGE ), data = AgeCorrected$Ages, inherit.aes = F, color = "green") +
  ggplot2::geom_point(ggplot2::aes(SAMPLE, AGE ), data = AgeAsBayLum$Ages, inherit.aes = F, color = "red")


###-------------------------------------------------------------------------------@
datasets2 <- readxl::read_excel("R/DataManipulation/dataOSL_3sites.xlsx")
colnames(datasets2)[3:6] <- c("ddot", "sddot", "D", "sD")

BlockLength <- datasets2 %>% dplyr::mutate(block = stringr::str_extract(Sample, "Dhab\\d+")) %>%
  dplyr::group_by(block) %>%
  dplyr::summarise(blockLength = dplyr::n()) %>% dplyr::select(blockLength)

Sc = as.matrix(Matrix::bdiag(sapply((BlockLength$blockLength), function(n) upper.tri(matrix(1, nrow = n, ncol = n))))) *1
Sc = as.matrix( rbind(rep(1, dim(datasets2)[1]), Sc))

network = igraph::graph_from_adjacency_matrix(Sc[-1, ], mode = "directed")
plot(
  network,
  layout = igraph::layout_with_sugiyama(network),
  vertex.size = 16,
  vertex.color = adjustcolor("lightblue", alpha.f = 0.8),
  edge.arrow.size = 0.4,  # Smaller arrowheads
  edge.width = 1,
  asp = 0,
  edge.curved = 0.1
)

reduced_network = remove_transitive_edges(network)
V(reduced_network)$name <- paste0("A", 1:13)


plot(
  reduced_network,
  layout = igraph::layout_with_sugiyama(reduced_network),
  vertex.size = 16,
  vertex.color = adjustcolor("lightblue", alpha.f = 0.8),
  edge.arrow.size = 0.4,  # Smaller arrowheads
  edge.width = 1,
  asp = 0,
  edge.curved = 0.1
)

Measures = list(ddot = datasets2$ddot, sddot = datasets2$sddot, D = datasets2$D, sD = datasets2$sD,
                sddot_shared = rep(0,dim(datasets2)[1]), Nb_sample = dim(datasets2)[1],
                SampleNames = datasets2$Sample)


Theta = diag(as.numeric(Measures$sddot)**2) +
  as.numeric(Measures$sddot_shared) %*% t(as.numeric(Measures$sddot_shared))

covD = diag(as.numeric(Measures$sD)**2)

DtMeasures <- list(Theta = Theta, Measures = Measures, covD = covD)

GibbsOutput <-  GibbsSampler(DtMeasures, 3, 70000, 50000, Sc, 1, 1400, "logit", proposal_magnitude = 2, lambda = 1, n_logitStep = 5,
                             type = "other")
pdf("R/DataManipulation/blockdataset.pdf", width = 17)
plot(GibbsOutput$Sampling)
coda::acfplot(GibbsOutput$Sampling)
plot(GibbsOutput$name_chains)
coda::acfplot(GibbsOutput$name_chains)
dev.off()
AgeAsBayLum <-Compute_AgeS_D(DtMeasures, Sc, ModelAgePrior$Jeffreys, Iter = 2000, burnin = 50000, t = 10, adapt = 1000,
                             PriorAge = rep(c(1, 1400),  Measures$Nb_sample))


plotHpd(list(GibbsOutput, AgeAsBayLum), c("Gibbs", "JagsBayLum")) +  ggplot2::geom_point(mapping = ggplot2::aes(SAMPLE, AGE), data = GibbsOutput$Ages, inherit.aes = F,
                                                                                         color = "red") +
  ggplot2::geom_point(mapping = ggplot2::aes(SAMPLE, AGE), data = AgeAsBayLum$Ages, inherit.aes = F,
                      color = "blue")

## ============ Elaine data sets ======== ##
library(igraph)
extractElaine <- function(name) {
  path = paste0("R/DataManipulation/", name)
    #strati
  SampleNames <- readxl::read_xlsx(path, sheet = 2)[, 1] #Genral Strat
  # Depth <- as.matrix(readxl::read_xlsx(path, sheet = 2)[, 2])[,1]
  OSLestimate <- as.matrix(readxl::read_xlsx(path, sheet = 3)[, 5]) # Central Dose and Their uncertainties with No Strati
  lower95 =  as.matrix(readxl::read_xlsx(path, sheet = 3)[16:30, 3])
  upper95 = as.matrix(readxl::read_xlsx(path, sheet = 3)[16:30, 7])
  sD = (upper95-lower95)/ (2*1.96)
  Theta <- as.matrix(readxl::read_xlsx(path, sheet = 6))
  Sc <- as.matrix(readxl::read_xlsx(path, sheet = 7))
  n = dim(Theta)[1]
  ddot =  as.numeric(readxl::read_xlsx(path, sheet = 5)[3,c(12,13)])

  Measures = list(SampleNames = SampleNames$Layer, D = OSLestimate[16:30], sD = sD, Nb_sample = n ,
                  ddot = rep(ddot[1], n), sddot = rep(ddot[2], n), sddot_shared = rep(0, n), Depth = sort(Depth, decreasing = F))

  return( list(Measures = Measures, Theta = Theta, covD = diag(as.numeric(Measures$sD)**2), Sc = Sc))

}

extractElaine <- function(name) {
  path = paste0("R/DataManipulation/", name)
  #strati
  SampleNames <- readxl::read_xlsx(path, sheet = 2)[, 1]
  OSLestimate <- as.matrix(readxl::read_xlsx(path, sheet = 3)[, c(4,7)]) # Central Dose and Their uncertainties with No Strati

  Theta <- as.matrix(readxl::read_xlsx(path, sheet = 5))
  Sc <- as.matrix(readxl::read_xlsx(path, sheet = 6))
  n = dim(Theta)[1]
  ddot =  as.numeric(readxl::read_xlsx(path, sheet = 4)[3,c(12,13)])

  Measures = list(SampleNames = SampleNames$Layer, D = OSLestimate[,1], sD = OSLestimate[,2], Nb_sample = n ,
                  ddot = rep(ddot[1], n), sddot = rep(ddot[2], n), sddot_shared = rep(0, n))

  return( list(Measures = Measures, Theta = Theta, covD = diag(as.numeric(Measures$sD)**2), Sc = Sc))

}

DtMeasures <- extractElaine("All_Qz_Grains_BayLum_doses.xlsx")
DtMeasures
# depth = read.csv("R/DataManipulation/Layer_DepthsElaine.csv")
# colnames(depth) <- c("samples", "depth")

# DtMeasures$Sc = rbind(rep(1, Measures$Nb_sample), matrix(rep(0), ncol = DtMeasures$Measures$Nb_sample, nrow = DtMeasures$Measures$Nb_sample))

data.frame(Samples = DtMeasures$Measures$SampleNames, D = DtMeasures$Measures$D, sD = DtMeasures$Measures$sD)

ElaineApprox = AgeApprox(DtMeasures)
G = igraph::graph_from_adjacency_matrix(DtMeasures$Sc[-1, ])
plot(
  G,
  layout = layout_with_sugiyama(G),
  vertex.size = 11,
  vertex.color = adjustcolor("orange", alpha.f = 0.6),
  edge.arrow.size = 0.3,  # Smaller arrowheads
  edge.width = 1,
  edge.length = DtMeasures$Measures$Depth[2:15]-DtMeasures$Measures$Depth[1:14],
  asp = 0,
  edge.curved = 0.1
)

reduced_G <- remove_transitive_edges(G)
plot(
  reduced_G,
  layout = layout_with_sugiyama(reduced_G, weights = DtMeasures$Measures$Depth[2:15]-DtMeasures$Measures$Depth[1:14]),
  vertex.size = 11,
  vertex.color = adjustcolor("orange", alpha.f = 0.6),
  edge.arrow.size = 0.3,  # Smaller arrowheads
  edge.width = 1,
  asp = 0,
  edge.length = DtMeasures$Measures$Depth[2:15]-DtMeasures$Measures$Depth[1:14],
  edge.curved = 0.1
)

IsoBloc <- IsotonicRegDAG(reduced_G, ElaineApprox$Ahat, (1/ElaineApprox$sdAhat**2))

plot(c(min(DtMeasures$Measures$Depth),max(DtMeasures$Measures$Depth)), c(min(ElaineApprox$Ahat)*.5, max(ElaineApprox$Ahat)*1.5), type = "n",
     xlab = "depth", ylab = "Ages")
points(DtMeasures$Measures$Depth, IsoBloc$A, pch = 12, col = "red")
lines(DtMeasures$Measures$Depth,IsoBloc$A,  col = "red")
points(DtMeasures$Measures$Depth, ElaineApprox$Ahat, pch=13, col = "green")

legend("topright", c("IsotonicReg", "OSLApprox"), col = c("red", "green"), pch = c(12, 13))



GibbsOutput <- GibbsSampler(DtMeasures,3,  50000, 30000, DtMeasures$Sc, 1, 1400, "logit", proposal_magnitude = 2,
                            lambda = .1, n_logitStep = 5, type = "isotonic", lag = 20)
GibbsTrunc <- GibbsSamplerTrunc(DtMeasures, 3, 70000, 50000, DtMeasures$Sc, 1, 1400, type = "isotonic", proposal_magnitude = 2)

pdf("R/DataManipulation/Elainedataset.pdf", width = 17)
plot(GibbsOutput$Sampling)
plot(GibbsTrunc$Sampling)
coda::acfplot(GibbsOutput$Sampling)
coda::acfplot(GibbsTrunc$Sampling)
plot(GibbsOutput$name_chains)
coda::acfplot(GibbsOutput$name_chains)
dev.off()

plotHpd(list(GibbsOutput, GibbsTrunc), c("WithTrans", "WithoutTrans"))


AgeAsBayLum <-Compute_AgeS_D(DtMeasures, DtMeasures$Sc, prior = "Jeffreys", Iter = 2000, burnin = 30000, t = 10,
                             PriorAge = rep(c(1, 1400),  DtMeasures$Measures$Nb_sample))

plotHpd(list(AgeAsBayLum, GibbsTrunc), c("WithTrans", "WithoutTrans"))


IndepSc = rlang::duplicate(DtMeasures$Sc)
IndepSc[2:16, ] = IndepSc[2:16, ] *0
IndepSc
Independant  <- GibbsSamplerTrunc(DtMeasures, 3, 70000, 50000, IndepSc, 1, 1400, type = "isotonic")

DtIso <- Independant$Summary[, c(1,3, 5, 9)]
lower = IsotonicRegDAG(reduced_G, Ahat = DtIso[, 1], weights = 1/ (DtIso[, 4])**2)$A
upper = IsotonicRegDAG(reduced_G, Ahat = DtIso[, 3], weights = 1/ (DtIso[, 4])**2)$A
data.frame(lower = lower, upper = upper, depth = DtMeasures$Measures$Depth) %>%
  ggplot2::ggplot(ggplot2::aes(x = depth, ymin = lower, ymax = upper), fill = "orange") +
  ggplot2::geom_ribbon(alpha = .4) +
  ggplot2::geom_line(ggplot2::aes(y = lower), color = "orange", group = 1) +
  ggplot2::geom_line(ggplot2::aes(y = upper), color = "orange", group = 1)  + ggplot2::theme_minimal()+
  ggplot2::ylab("IsotonicRegression") + ggplot2::geom_point(mapping = ggplot2::aes(depth, upper), inherit.aes = F, color = "blue") +
  ggplot2::geom_point(mapping = ggplot2::aes(depth, lower), inherit.aes = F, color = "red")

AgeCorrected <- Compute_AgeS_D(DtMeasures, DtMeasures$Sc, prior = "StrictOrder", Iter = 2000, burnin = 50000, t = 10,
                               PriorAge = rep(c(1, 1400),  DtMeasures$Measures$Nb_sample))

AgeNicholls <- Compute_AgeS_D(DtMeasures, DtMeasures$Sc, prior = "StrictNicholls", Iter = 2000, burnin = 30000, t = 10,
                               PriorAge = rep(c(1, 1400),  DtMeasures$Measures$Nb_sample))

IsoData = rlang::duplicate(AgeNicholls$Ages[, 1:2]) %>% dplyr::mutate(AGE = IsoBloc$A, depth = DtMeasures$Measures$Depth)

methods = list(AgeCorrected, AgeNicholls, Independant)

my_list <- lapply(methods, function(inner_list){
  inner_list$Ages$depth = DtMeasures$Measures$Depth
  inner_list
})

plot_Ages(AgeAsBayLum)
my_list[[1]]$Ages

plotHpd(methods,
        c( "BayLum", "Nicholls", "Independant")) +
  ggplot2::geom_point(mapping = ggplot2::aes(SAMPLE, AGE), data = methods[[3]]$Ages, inherit.aes = F,
                                                                                   color = "green") +
  ggplot2::geom_point(ggplot2::aes(SAMPLE, AGE ), data = methods[[1]]$Ages, inherit.aes = F, color = "red") +
  ggplot2::geom_point(ggplot2::aes(SAMPLE, AGE ), data = methods[[2]]$Ages, inherit.aes = F, color = "blue") +
  ggplot2::geom_point(ggplot2::aes(SAMPLE, AGE ), data = IsoData, inherit.aes = F, color = "black", shape = 18,
                      size = 3)
##-----------------------------------------------##
simulated = read.csv("R/DataManipulation/osl_constrained_data.csv")
commonError = .03
colnames(simulated)[2:7] <- c("D", "sD", "ddot", "sddot", "Age", "sAge" )
simulated
Measures <- list(SampleNames = simulated$Sample,
                 ddot = simulated$ddot,  # Environmental dose rate (Gy/ka)
                 sddot = simulated$sddot, # Error on dose rate (Gy/ka)
                 D = simulated$D,  # OSL Dose (Gy)
                 sD = simulated$sD,  # Error on dose (Gy)
                 sddot_shared = rep(commonError, 15),  # Common error term (Gy)
                 Nb_sample = 15)

Sc = rbind(rep(1, Measures$Nb_sample), upper.tri(matrix(rep(1), ncol = Measures$Nb_sample, nrow = Measures$Nb_sample))*1)
# Sc = rbind(rep(1, Measures$Nb_sample), matrix(rep(0), ncol = Measures$Nb_sample, nrow = Measures$Nb_sample))
G = graph_from_adjacency_matrix(Sc[-1, ], mode = "directed")
plot(
  G,
  layout = layout_with_sugiyama(G),
  vertex.size = 11,
  vertex.color = adjustcolor("orange", alpha.f = 0.6),
  edge.arrow.size = 0.3,  # Smaller arrowheads
  edge.width = 1,
  asp = 0,
  edge.curved = 0.1
)

reduced_G <- remove_transitive_edges(G)

plot(
  reduced_G,
  layout = layout_with_sugiyama(reduced_G),
  vertex.size = 11,
  vertex.color = adjustcolor("orange", alpha.f = 0.6),
  edge.arrow.size = 0.3,  # Smaller arrowheads
  edge.width = 1,
  asp = 0,
  edge.curved = 0.1
)

Theta = diag(as.numeric(Measures$sddot)**2) +
  as.numeric(Measures$sddot_shared) %*% t(as.numeric(Measures$sddot_shared))
covD = diag(as.numeric(Measures$sD)**2)

DtMeasures <- list(Theta = Theta, Measures = Measures, covD = covD)
GibbsOutput = GibbsSampler(DtMeasures, 3, 70000, 50000, Sc, 1, 200, "logit", lambda = .5, n_logitStep = 5, type = "other")
pdf("R/DataManipulation/simulated15.pdf", width = 18)
plot(GibbsOutput$Sampling)
coda::acfplot(GibbsOutput$Sampling)
plot(GibbsOutput$name_chains)
coda::acfplot(GibbsOutput$name_chains)
dev.off()




AgeAsBayLum <-Compute_AgeS_D(DtMeasures, Sc, prior = "Jeffreys", Iter = 2000, burnin = 50000, t = 10,
                             PriorAge = rep(c(1, 200),  DtMeasures$Measures$Nb_sample))
AgeCorrected <- Compute_AgeS_D(DtMeasures, Sc, prior = "StrictOrder", Iter = 2000, burnin = 50000,t = 10,
                               PriorAge = rep(c(1, 200),  Measures$Nb_sample))
AgeConditionnal <- Compute_AgeS_D(DtMeasures, Sc, prior = "Conditional", Iter = 2000, burnin = 50000, t = 10,
                                 PriorAge = rep(c(1, 200),  Measures$Nb_sample))

AgeNicholls <- Compute_AgeS_D(DtMeasures,Sc, prior = "StrictNicholls", Iter = 2000, burnin = 30000, t = 10,
                              PriorAge = rep(c(1, 1400),  DtMeasures$Measures$Nb_sample))

plotHpd(list(GibbsOutput, AgeConditionnal, AgeCorrected, AgeAsBayLum), c("Gibbs", "Conditionnal", "Exponential", "BayLum")) +
  ggplot2::geom_point(mapping = ggplot2::aes(SAMPLE, AGE), data = AgeConditionnal$Ages, inherit.aes = F,
                                                                                    color = "green") +
  ggplot2::geom_point(ggplot2::aes(SAMPLE, AGE), data = GibbsOutput$Ages, inherit.aes = F, color = "purple") +
  ggplot2::geom_point(ggplot2::aes(SAMPLE, AGE), data = AgeAsBayLum$Ages, inherit.aes = F, color = "red") +
  ggplot2::geom_point(ggplot2::aes(Sample,Age ), data = simulated, inherit.aes = F, alpha = .6)


##======================================================================================#
#### Données Guillaume ####
# DEMANDER A GUILLAUME DE PARTAGER COMMENT SONT CALCULÉES LES VALEURS DE LA MATRICE DE COV
# IMPOSSIBLE D'UTILISER LE GIBBS A LA MAIN POUR L'INSTANT CAR IL FAUT AVOIR CERTAINES DONNÉES NON PRÉCISÉES ICI (VOIR LATEX)
#### POUR TOUS LES ECHANTILLONS
dt <- readxl::read_xlsx("R/DataManipulation/DoseAndDRestimates-forModelling.xlsx", skip = 1)[, -8]
colnames(dt) <- c("unit", "samples", "lower95", "lower68", "estimate", "upper68", "upper95",
                  "ddot", "varddot")
dt
Measures <- list(
  SampleNames = dt$samples,
  ddot = dt$ddot,  # Environmental dose rate (Gy/ka)
  sddot = sqrt(dt$varddot), # Error on dose rate (Gy/ka)
  D = dt$estimate,  # OSL Dose (Gy)
  sD = (dt$upper95-dt$lower95) /(2*1.96),  # Error on dose (Gy)
  Nb_sample = length(dt$unit)
)

Measures


Theta <- as.matrix(read.csv("R/DataManipulation/CovarianceMatrix_1313.csv", sep = ";",
                            col.names = paste0("A", 1:Measures$Nb_sample)))

covD = diag(Measures$sD**2)
DtMeasures <- list(Measures = Measures, Theta = Theta, covD = covD)
Sc = rbind(rep(1, Measures$Nb_sample), upper.tri(matrix(rep(1), ncol = Measures$Nb_sample, nrow = Measures$Nb_sample))*1)

G = igraph::graph_from_adjacency_matrix(Sc[-1, ])
plot(
  G,
  layout = igraph::layout_with_sugiyama(G),
  vertex.size = 11,
  vertex.color = adjustcolor("orange", alpha.f = 0.6),
  edge.arrow.size = 0.3,  # Smaller arrowheads
  edge.width = 1,
  asp = 0,
  edge.curved = 0.1
)


reduced_G <- remove_transitive_edges(G)
plot(
  reduced_G,
  layout = igraph::layout_with_sugiyama(reduced_G),
  vertex.size = 11,
  vertex.color = adjustcolor("orange", alpha.f = 0.6),
  edge.arrow.size = 0.3,  # Smaller arrowheads
  edge.width = 1,
  asp = 0,
  edge.curved = 0.1
)

AgeCorrected <- Compute_AgeS_D(DtMeasures, Sc, prior = "StrictOrder", Iter = 2000, burnin = 50000, t = 10,
                               PriorAge = rep(c(1, 1400),  DtMeasures$Measures$Nb_sample))

AgeNicholls <- Compute_AgeS_D(DtMeasures,Sc, prior = "StrictNicholls", Iter = 2000, burnin = 50000, t = 10,
                              PriorAge = rep(c(1, 1400),  DtMeasures$Measures$Nb_sample))
AgeAsBayLum <- Compute_AgeS_D(DtMeasures, Sc, prior = "Jeffreys", Iter = 2000, burnin = 50000, t = 10,
                              PriorAge = rep(c(1, 1400),  DtMeasures$Measures$Nb_sample))

Sc = rbind(rep(1, Measures$Nb_sample), upper.tri(matrix(rep(0), ncol = Measures$Nb_sample, nrow = Measures$Nb_sample))*0)

Independant <-Compute_AgeS_D(DtMeasures, , prior = "Jeffreys", Iter = 2000, burnin = 50000, t = 10,
                             PriorAge = rep(c(1, 1400),  DtMeasures$Measures$Nb_sample))

sample = as.matrix(runjags::combine.mcmc(Independant$Sampling))
w
IsotonicRegDAG(reduced_G, sample[1,], 1/w^2)$A
apply(sample, 1, function(Ahat, network, weights) IsotonicRegDAG(network, Ahat, weights )$A,
      network = reduced_G , weights = 1/w^2)

IsoSample = ISotonicCurve(reduced_G, Independant)
IsoSample = t(IsoSample)
colnames(IsoSample) <- DtMeasures$Measures$SampleNames
HPDISO = hpd_method("IsotonicReg",IsoSample ) %>% dplyr::mutate(Samples = factor(Samples, levels = DtMeasures$Measures$SampleNames))

HPDISO <- HPDISO %>% dplyr::mutate(avg = apply(IsoSample, 2, mean))
DtIso <- Independant$Summary[, c(1,3, 5, 9)]
lower = IsotonicRegDAG(G, Ahat = DtIso[, 1], weights = 1/ (DtIso[, 4])**2)$A
upper = IsotonicRegDAG(G, Ahat = DtIso[, 3], weights = 1/ (DtIso[, 4])**2)$A

data.frame(lower = HPDISO$inf, upper = HPDISO$sup, Unit = 1:13, Samples =HPDISO$Samples, avg = HPDISO$avg) %>%
  ggplot2::ggplot(ggplot2::aes(x = Unit, ymin = lower, ymax = upper), fill = "orange") +
  ggplot2::geom_ribbon(alpha = .4) +
  ggplot2::geom_line(ggplot2::aes(y = lower), color = "orange", group = 1) +
  ggplot2::geom_line(ggplot2::aes(y = upper), color = "orange", group = 1)  + BayLumTheme() +
  ggplot2::geom_line(ggplot2::aes(y = avg), color = "orange", group = 1, size =1.5) +
  ggplot2::ylab("IsotonicRegression") + ggplot2::scale_x_continuous(breaks = unique(dt$unit)) +
  ggplot2::geom_point(ggplot2::aes(x = Unit, y = lower), color = "blue") +
  ggplot2::geom_point(ggplot2::aes(x = Unit, y = upper), color = "red") +
  ggplot2::geom_point(ggplot2::aes(x = Unit, y = avg), color = "black") +
  ggplot2::scale_x_continuous(breaks = 1:13, labels = HPDISO$Samples) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45))

plotHpd(list(AgeCorrected, AgeNicholls), c("LogUniform", "Nicholls")) +
  ggplot2::geom_linerange(ggplot2::aes(x = Samples, ymin= inf, ymax = sup, colour = Models),
        data = HPDISO,
       position = ggplot2::position_dodge(.5))


## comparing estimate points
Approx = AgeApprox(DtMeasures)
estimateMethod = data.frame(IsoReg = Iso::pava(Approx$Ahat, 1/(Approx$sdAhat**2)), Nicholls = AgeNicholls$Ages$AGE,
                            BayLum = AgeAsBayLum$Ages$AGE, Jeffrey = AgeCorrected$Ages$AGE,
                            Independant = Independant$Ages$AGE, Unit = dt$unit)
estimateMethod %>% tidyr::pivot_longer(-Unit) %>% dplyr::arrange(name) %>%
  ggplot2::ggplot( ggplot2::aes(x = Unit, y =value, color = name, shape = name )) + ggplot2::geom_point(size = 2) +
  ggplot2::geom_line() +
  ggplot2::theme_minimal() + ggplot2::ylab("Age")

ggplot2::ggplot(estimateMethod, ggplot2::aes(y = Independant, x = Unit)) + ggplot2::geom_line() +ggplot2::geom_point() +
  ggplot2::ylab("Age")

#### POUR SEULEMENT LES 11 SELECTIONÉS
Measures <- list(
  SampleNames = dt$samples[-c(8,11)],
  ddot = dt$ddot[-c(8,11)],  # Environmental dose rate (Gy/ka)
  sddot = sqrt(dt$varddot)[-c(8,11)], # Error on dose rate (Gy/ka)
  D = dt$estimate[-c(8,11)],  # OSL Dose (Gy)
  sD = (dt$upper95-dt$lower95)[-c(8,11)] /(2*1.96),  # Error on dose (Gy)
  Nb_sample = length(dt$unit)-2
)

Theta <- as.matrix(read.csv("R/DataManipulation/CovarianceMatrix_1111.csv", sep = ";",
                            col.names = paste0("A", 1:Measures$Nb_sample)))
covD = diag(Measures$sD**2)

DtMeasures <- list(Measures = Measures, Theta = Theta, covD = covD)
Sc = rbind(rep(1, Measures$Nb_sample), upper.tri(matrix(rep(1), ncol = Measures$Nb_sample, nrow = Measures$Nb_sample))*1)

AgeCorrected <- Compute_AgeS_D(DtMeasures, Sc, prior = "StrictOrder", Iter = 2000, burnin = 50000, t = 10,
                               PriorAge = rep(c(1, 1400),  DtMeasures$Measures$Nb_sample))

AgeNicholls <- Compute_AgeS_D(DtMeasures,Sc, prior = "StrictNicholls", Iter = 2000, burnin = 50000, t = 10,
                              PriorAge = rep(c(1, 1400),  DtMeasures$Measures$Nb_sample))

AgeNichollsBR <- Compute_AgeS_D(DtMeasures,Sc, prior = "nichollsBR", Iter = 2000, burnin = 50000, t = 10,
                                PriorAge = rep(c(1, 1400),  DtMeasures$Measures$Nb_sample))

AgeAsBayLum <- Compute_AgeS_D(DtMeasures, Sc, prior = "Jeffreys", Iter = 2000, burnin = 50000, t = 10,
                              PriorAge = rep(c(1, 1400),  DtMeasures$Measures$Nb_sample))

Sc = rbind(rep(1, Measures$Nb_sample), upper.tri(matrix(rep(0), ncol = Measures$Nb_sample, nrow = Measures$Nb_sample))*0)

Independant <-Compute_AgeS_D(DtMeasures, , prior = "Jeffreys", Iter = 2000, burnin = 50000, t = 10,
                             PriorAge = rep(c(1, 1400),  DtMeasures$Measures$Nb_sample))

plotHpd(list(AgeCorrected, AgeNicholls, AgeAsBayLum, Independant,
             AgeNichollsBR), c("JeffreysOrder", "Nicholls", "BayLum", "Independant", "nichollsBR"))

Approx = AgeApprox(DtMeasures)
estimateMethod = data.frame(IsoReg = Iso::pava(Approx$Ahat, 1/(Approx$sdAhat**2)), Nicholls = AgeNicholls$Ages$AGE,
                            BayLum = AgeAsBayLum$Ages$AGE, Jeffrey = AgeCorrected$Ages$AGE, NichollsBR = AgeNichollsBR$Ages$AGE,
                            Independant = Independant$Ages$AGE, Unit = dt$unit[-c(8,11)])
estimateMethod %>% tidyr::pivot_longer(-Unit) %>% dplyr::arrange(name) %>%
  ggplot2::ggplot( ggplot2::aes(x = Unit, y =value, color = name, shape = name )) + ggplot2::geom_point(size = 2) +
  ggplot2::geom_line() +
  ggplot2::theme_minimal() + ggplot2::ylab("Age")

data.frame(IsoReg = Iso::pava(Approx$Ahat, 1/(Approx$sdAhat**2)), NichollsBR = AgeNichollsBR$Ages$AGE,
           Independant = Independant$Ages$AGE, Unit = dt$unit[-c(8,11)])  %>% tidyr::pivot_longer(-Unit) %>% dplyr::arrange(name) %>%
  ggplot2::ggplot( ggplot2::aes(x = Unit, y =value, color = name, shape = name )) + ggplot2::geom_point(size = 2) +
  ggplot2::geom_line() +
  ggplot2::theme_minimal() + ggplot2::ylab("Age")

##======================================================================================#
#### Jingbian loess ####
OslJingbian <- readxl::read_xlsx("R/DataManipulation/OSL_EastMound.xlsx", sheet = "Sheet1")
colnames(OslJingbian)[c(4, 7)] <- c("Depth", "std")
dim_data <- dim(OslJingbian)
OslJingbian <- OslJingbian %>% dplyr::select(Depth, Age, std)
IsoJingbian = Iso::pava(OslJingbian$Age, 1/(OslJingbian$std)**2)
OslJingbian <- OslJingbian %>% dplyr::mutate(iso = IsoJingbian)
ggplot2::ggplot(data = OslJingbian, ggplot2::aes(x = Depth, y = Age, color = "Measures")) + ggplot2::geom_point(size = 2, shape = 15) +
  ggplot2::geom_point(ggplot2::aes(x = Depth, y= iso, color = "IsotonicReg"), size = 2) + BayLumTheme() +ggplot2::labs(colour = "Method")

## apply the uncertainty after only using the independant bayesian model


##======================================================================================#

#### Graph Clustering ####
Sc = matrix(0, 6,6)
Sc[1, ] = c(rep(0,2), rep(1,2), rep(0,2))
Sc[2, ] = c(rep(0,2), 1,0, rep(1,2))
Sc[5,]= c(rep(0,5),1)


G = igraph::graph_from_adjacency_matrix(Sc)
plot(G, layout = igraph::layout_with_sugiyama(G))
### avec la matrice d'adj
hc_basic <- hclust(as.dist(1-Sc), method = "ward.D2")
plot(hc_basic, labels = F)

igraph::topo_sort(G, mode = "out")
igraph::dfs(G, 1)
igraph::bfs(G,2)

Dpath = igraph::distances(G)
Dpath[which(Dpath == Inf, arr.ind = T)] = -2
corrplot::corrplot(Dpath, tl.pos = 'n', is.corr = FALSE)
hc_path <- hclust(as.dist(Dpath), method = "ward.D2")
plot(hc_path, labels = F)

# hc_modularity <- igraph::cluster_fast_greedy(G) ONLY WORKS FOR UNDIRECT GRAPH NOT THE CASE HERE
hc_bet = igraph::cluster_edge_betweenness(G)
plot(hc_bet, G, layout = igraph::layout_with_sugiyama(G))
igraph::membership(hc_bet)
degree = igraph::degree(G)
vertices = igraph::V(G)
n_nodes = igraph::gorder(G)
n_edges = igraph::gsize(G)



network = igraph::graph_from_adjacency_matrix(SC[-1,])
plot(network, layout = igraph::layout_with_sugiyama(network))
degree = igraph::degree(network)
vertices = igraph::V(network)
n_nodes = igraph::gorder(network)
n_edges = igraph::gsize(network)
hc_bet = igraph::cluster_edge_betweenness(network)
plot(hc_bet, network, layout = igraph::layout_with_sugiyama(network))
igraph::tkplot(network, layout = igraph::layout_with_sugiyama(network))

L <- igraph::laplacian_matrix(network)
spec_L <- eigen(L)
practical_zero <- 1e-12
lambda  <- min(spec_L$values[spec_L$values>practical_zero])
fiedler <- spec_L$vectors[, which(spec_L$values == lambda)]

ggplot2::qplot(y = fiedler) +
  viridis::scale_color_viridis(discrete = TRUE) + ggplot2::theme_minimal()


opt = igraph::cluster_optimal(G)
plot(opt,G)



SC = matrix(data=0,ncol=40,nrow=41)
SC[1,]=rep(1,40)
SC[2,]=c(rep(0,1),rep(1,1),rep(0,3),rep(1,35)) #4853 (OxA -2161)
SC[3,]=c(rep(0,5),rep(1,35)) #4861 (OxA-2162)
SC[4,]=c(rep(0,5),rep(1,35)) #4555 (UCIAMS-98210)
SC[5,]=c(rep(0,5),rep(1,35)) #4779 (UCIAMS-103134)
SC[6,]=rep(0,40) #P782 After (Hearth Level 10)
SC[7,]=c(rep(0,9),rep(1,31)) #4850
SC[8,]=c(rep(0,9),rep(1,31)) #4850
SC[9,]=rep(0,40) #4715 After
SC[10,]=rep(0,40) #4715 After
SC[11,]=c(rep(0,17),rep(1,23)) #4517
SC[12,]=c(rep(0,17),rep(1,23)) #4517
SC[13,]=rep(0,40) #4822 After
SC[14,]=rep(0,40) #4826 After
SC[15,]=rep(0,40) #4826 After
SC[16,]=c(rep(0,16),rep(1,24)) #4828
SC[17,]=c(rep(0,16),rep(1,24)) #4828
SC[18,]=c(rep(0,17),rep(1,1),rep(0,4),rep(1,1),rep(0,1),rep(1,16)) #4836
SC[19,]=c(rep(0,22),rep(1,1),rep(0,1),rep(1,16)) #4837
SC[20,]=c(rep(0,20),rep(1,2),rep(0,1),rep(1,17)) #4869
SC[21,]=c(rep(0,20),rep(1,2),rep(0,1),rep(1,17)) #4869
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
SC[40,]=rep(0,40) #pre XII After
SC[41,]=rep(0,40) #pre XII After


Sc = SC[-1,]
library(igraph)

sparsingSc <- function(col) {
  n = length(col)
  index = which(col == 1)
  if (length(index) == 0) {
    return(col)
  }
  else {

  m = min(index)
  new_col = rep(0, n)
  new_col[m] = 1
  return(new_col)
  }
}

bound <- function(index,  Sc) {
  Sc = Sc[-1, ] #del lower bound line

  lowerbound_index = which(Sc[, index] == 1)
  if (length(lowerbound_index) > 0) { a = max(lowerbound_index)}

  else { a = 0}

  upperbound_index = which(Sc[index, ] == 1)
  if (length(upperbound_index) > 0 ) {b = min(upperbound_index)}
  else {b = -1}
  return(c(a,b))
}

bounds = sapply(1:40, bound, SC)
bounds
Sc = matrix(0, 40, 40)
for (age in 1:40) {
  inf = bounds[1, age]
  sup = bounds[2,age]
  if (inf == 0 & sup ==-1) {
    next
  }
  else if (inf == 0) {
    Sc[age, sup] = 1
  }
  else if (sup == -1) {
    Sc[inf, age] = 1
  }
  else {
    Sc[inf, sup ] = 1
  }
}
Sc

M = SC[-1, ]

for (age in 1:40) {
  if (sum(M[, age] == age)) {
    M[, (age-1)]
  }
}


Sc = SC[-1,]
library(igraph)
network <- igraph::graph_from_adjacency_matrix(Sc[1:10,1:10], mode = "directed")
uu= V(network)
for (u in uu) {
  message(paste("check des voisins du sommet", u))
  nei = neighbors(network, u, mode = "out")
  for (v in nei) {

    print(subcomponent(network, v, mode = "out"))
  }
}

uu[1]
v = neighbors(network, uu[1], mode = "out")[1]
des = subcomponent(network, 2, mode = "out")
des[1]

remove_transitive_edge <- function(G) {
  reduced_G = G
  for (u in igraph::V(G)) {
    for (v in igraph::neighbors(G, u, mode = "out")) {
      descendants <- igraph::subcomponent(G, v, mode = "out")
      ndex = length(descendants)
      for (w in descendants[-1]) {
        print(E(reduced_G, P = c(u,w)))
        #if w is a direct descendant of v and not u then we delete the edges between w and u
        if (igraph::are_adjacent(G, u, w)) {
          reduced_G = igraph::delete_edges(reduced_G, igraph::E(reduced_G, P = c(u,w)))
        }
      }
    }
  }
  return(reduced_G)
}

reduced_network = remove_transitive_edge(network)

layout <- igraph::layout_with_sugiyama(reduced_network, hgap = 1)
plot(reduced_network, layout = igraph::layout_with_sugiyama(reduced_network, hgap = 5, vgap = 5),
     vertex.size = 8,
     vertex.color = adjustcolor("orange", alpha.f = 0.6),
     edge.arrow.size = 0.3,  # Smaller arrowheads
     edge.width = 1,
    asp = 0,
     edge.curved = 0.1)
igraph::tkplot(network)

igraph::topo_sort(network)
igraph::degree(network, igraph::topo_sort(network))

opt = igraph::cluster_optimal(network)
plot(opt,network,  layout = igraph::layout_with_sugiyama(network, hgap = 5, vgap = 5))





# Adjust node numbering since row 1 was removed
edges$from <- as.character(edges$from)
edges$to <- as.character(edges$to)

# STEP 3: Create nodes data frame
unique_nodes <- unique(c(edges$from, edges$to))
nodes <- data.frame(id = unique_nodes, label = paste0("Node ", unique_nodes))

# STEP 4: Visualize with stable layout (hierarchical)
library(visNetwork)

visNetwork(nodes, edges) %>%
  visEdges(arrows = "to") %>%
  visHierarchicalLayout(direction = "Rl") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)




############## testing findbounds


all_bounds = findbounds(reduced_network)
all_bounds[[2]]










