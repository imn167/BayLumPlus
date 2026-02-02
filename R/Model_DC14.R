#' @title JAGS Models for OSL Age Estimation in \code{\link{Compute_AgeS_DC14}}
#'
#' @description
#' JAGS models used to estimate true ages based on data obtained from OSL measures and C14 measures. Both kind of measurements are used to computed a common stratigraphy.
#'
#'@md
#' @details
#' The models are designed to refine age estimation by integrating these measurements into a Bayesian framework.
#' @references
#' To cite this package, please use: citation("BayLumPlus")
#
#'@format Unconstrained : Model with log-uniform settings for OSL and plain uniform for C14.
# Unconstrained <- "model {
#               ###### Likelihood ####
#               ## Initialize Sigma matrix with zeros (off-diagonal for mixed types)
#               for (i in index_osl) {
#                   for (j in index_c14) {
#                       Sigma[i, j] <- 0
#                       Sigma[j,i] <- 0
#                   }
#               }
#               ## init off-diag for C14
#               for (i in 1:(length(order_c14) - 1)) {
#                 for (j in (i + 1):length(order_c14)) {
#                   Sigma[index_c14[i], index_c14[j]] <- 0
#                   Sigma[index_c14[j], index_c14[i]] <- 0
#                 }
#               }
#               ## OSL filling
#               for (i in order_osl) {
#                   for (j in order_osl){
#                     Sigma[index_osl[i], index_osl[j]] <-  Atilde[index_osl[i]]*Atilde[index_osl[j]] * Theta[index_osl[i], index_osl[j]]
#                   }
#                 mu[index_osl[i]] <-  Atilde[index_osl[i]] * ddot[i]
#               }
#               ###C14 filling
#               for (i in order_c14) {
#                   mu[index_c14[i]] <-  interp.lin(Atilde[index_c14[i]], xTableauCalib, yTableauCalib)
#                   Z[i] ~ dcat(c(0.1, 0.9))
#                   err[i] <- interp.lin(Atilde[index_c14[i]], xTableauCalib, zTableauCalib)
#                   sigma_C14[i] <-  (alpha[i]^(-Z[i]+2)*(pow(sigma[i],2)+pow(err[i],2)))
#                   Sigma[index_c14[i], index_c14[i]] <-  sigma_C14[i]
#               }
#
#               invSigma  <-  inverse(Sigma + encoded_covD)
#               M ~ dmnorm(mu, invSigma)
#
#
#
#
#               ### prior
#               for (i in order_osl) {
#                   u[index_osl[i]] ~ dunif(0,1)
#                   Atilde[index_osl[i]] <- exp(u[index_osl[i]] * log(xbound[2*index_osl[i]] / xbound[2*index_osl[i]-1])) + xbound[2*index_osl[i]-1]
#                   A[index_osl[i]] = Atilde[index_osl[i]] -0.075
#               }
#
#               for (i in order_c14) {
#                 u[index_c14[i]] ~dunif(0,1)
#                 Atilde[index_c14[i]] = u[index_c14[i]] * (xbound[2*index_c14[i]] - xbound[2*index_c14[i]-1]) + xbound[2*index_c14[i]-1]
#                 A[index_c14[i]] = Atilde[index_c14[i]]
#                 invalpha[i] ~ dgamma(3,4)
#                 alpha[i] = 1/invalpha[i]
#               }
#
#
#
#
# }  "
#
#'@format Constrained : Model with log-uniform order settings for OSL and plain uniform order for C14.
# Constrained <-  "model {
#               ###### Likelihood ####
#               ## Initialize Sigma matrix with zeros (off-diagonal for mixed types)
#               for (i in index_osl) {
#                   for (j in index_c14) {
#                       Sigma[i, j] <- 0
#                       Sigma[j,i] <- 0
#                   }
#               }
#
#               ## init off-diag for C14
#               for (i in 1:(length(order_c14) - 1)) {
#                 for (j in (i + 1):length(order_c14)) {
#                   Sigma[index_c14[i], index_c14[j]] <- 0
#                   Sigma[index_c14[j], index_c14[i]] <- 0
#                 }
#               }
#
#               ## OSL filling
#               for (i in order_osl) {
#                   for (j in order_osl){
#                     Sigma[index_osl[i], index_osl[j]] <-  Atilde[index_osl[i]]*Atilde[index_osl[j]] * Theta[index_osl[i], index_osl[j]]
#                   }
#                 mu[index_osl[i]] <-  Atilde[index_osl[i]] * ddot[i]
#               }
#               ###C14 filling
#               for (i in order_c14) {
#                   mu[index_c14[i]] <-  interp.lin(Atilde[index_c14[i]], xTableauCalib, yTableauCalib)
#                   Z[i] ~ dcat(c(0.1, 0.9))
#                   err[i] <- interp.lin(Atilde[index_c14[i]], xTableauCalib, zTableauCalib)
#                   sigma_C14[i] <-  (alpha[i]^(-Z[i]+2)*(pow(sigma[i],2)+pow(err[i],2)))
#                   Sigma[index_c14[i], index_c14[i]] <-  sigma_C14[i]
#               }
#
#               invSigma  <-  inverse(Sigma + encoded_covD)
#               M ~ dmnorm(mu, invSigma)
#
#
#
#
#               ### prior
#               for (i in 1:(I+1)) {
#                 e[i] ~dexp(1)
#               }
#
#
#               for (i in order_osl) {
#                   u[index_osl[i]]  <-  sum(e[1:index_osl[i]]) / sum(e[1:(I+1)])
#                   Atilde[index_osl[i]] <- exp(u[index_osl[i]] * log(xbound[2*index_osl[i]] / xbound[2*index_osl[i]-1])) + xbound[2*index_osl[i]-1]
#                   A[index_osl[i]] = Atilde[index_osl[i]] -0.075
#               }
#
#               for (i in order_c14) {
#                 u[index_c14[i]] <-  sum(e[1:index_c14[i]]) / sum(e[1:(I+1)])
#                 Atilde[index_c14[i]] = u[index_c14[i]] * (xbound[2*index_c14[i]] - xbound[2*index_c14[i]-1]) + xbound[2*index_c14[i]-1]
#                 A[index_c14[i]] = Atilde[index_c14[i]]
#                 invalpha[i] ~ dgamma(3,4)
#                 alpha[i] = 1/invalpha[i]
#               }
#
#
#
#
# }  "
#'@format NichollsJones : Nicholls&Jones' Age Model applied on ages directly
# StrictNichollsJones <- "model {
#               ###### Likelihood ####
#               ## Initialize Sigma matrix with zeros (off-diagonal for mixed types)
#               for (i in index_osl) {
#                   for (j in index_c14) {
#                       Sigma[i, j] <- 0
#                       Sigma[j,i] <- 0
#                   }
#               }
#
#               ## init off-diag for C14
#               for (i in 1:(length(order_c14) - 1)) {
#                 for (j in (i + 1):length(order_c14)) {
#                   Sigma[index_c14[i], index_c14[j]] <- 0
#                   Sigma[index_c14[j], index_c14[i]] <- 0
#                 }
#               }
#
#               ## OSL filling
#               for (i in order_osl) {
#                   for (j in order_osl){
#                     Sigma[index_osl[i], index_osl[j]] <-  Atilde[index_osl[i]]*Atilde[index_osl[j]] * Theta[index_osl[i], index_osl[j]]
#                   }
#                 mu[index_osl[i]] <-  Atilde[index_osl[i]] * ddot[i]
#               }
#               ###C14 filling
#               for (i in order_c14) {
#                   mu[index_c14[i]] <-  interp.lin(Atilde[index_c14[i]], xTableauCalib, yTableauCalib)
#                   Z[i] ~ dcat(c(0.1, 0.9))
#                   err[i] <- interp.lin(Atilde[index_c14[i]], xTableauCalib, zTableauCalib)
#                   sigma_C14[i] <-  (alpha[i]^(-Z[i]+2)*(pow(sigma[i],2)+pow(err[i],2)))
#                   Sigma[index_c14[i], index_c14[i]] <-  sigma_C14[i]
#               }
#
#               invSigma  <-  inverse(Sigma + encoded_covD)
#               M ~ dmnorm(mu, invSigma)
#
#
#
#
#               ### prior
#               for (i in 1:(I-1)) {
#                 e[i] ~dexp(1)
#               }
#
#               s ~dunif(0,1)
#               first ~ dunif(0, (1-s))
#               u[1] = first
#               u[I] = s + first
#
#               for (i in 2:(I-1)) {
#                 u[i] = (sum(e[1:(i-1)]) / sum(e[1:(I-1)])) *s + first
#               }
#
#
#               for (i in order_osl) {
#                   Atilde[index_osl[i]] <- exp(u[index_osl[i]] * log(xbound[2*index_osl[i]] / xbound[2*index_osl[i]-1])) + xbound[2*index_osl[i]-1]
#                   A[index_osl[i]] = Atilde[index_osl[i]] -0.075
#               }
#
#               for (i in order_c14) {
#                 Atilde[index_c14[i]] = u[index_c14[i]] * (xbound[2*index_c14[i]] - xbound[2*index_c14[i]-1]) + xbound[2*index_c14[i]-1]
#                 A[index_c14[i]] = Atilde[index_c14[i]]
#                 invalpha[i] ~ dgamma(3,4)
#                 alpha[i] = 1/invalpha[i]
#               }
#
# }  "
#
# Model_DC14 = list()
#
# Model_DC14$Unconstrained <- Unconstrained
# Model_DC14$Constrained <- Constrained
# Model_DC14$NichollsJones <- StrictNichollsJones
# #


"Model_DC14"





