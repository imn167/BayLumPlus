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
#'@format Unconstrained : Model with log-uniform settings for OSL, plain uniform for C14 and UTh.
# Unconstrained_OSL_C14_UTh <- "model {
#               ###### Likelihood ####
#               ## Initialize Sigma matrix with zeros (off-diagonal for mixed types)
#               for (i in index_osl) {
#                   for (j in index_c14) {
#                       Sigma[i, j] <- 0
#                       Sigma[j, i] <- 0
#                   }
#               }
#
#               for (i in index_osl) {
#                   for (k in index_uth) {
#                       Sigma[i, k] <- 0
#                       Sigma[k, i] <- 0
#                   }
#               }
#
#               for (j in index_c14) {
#                   for (k in index_uth) {
#                       Sigma[j, k] <- 0
#                       Sigma[k, j] <- 0
#                   }
#               }
#               ## init off-diag for C14 (only diag)
#               for (i in 1:(length(order_c14) - 1)) {
#                 for (j in (i + 1):length(order_c14)) {
#                   Sigma[index_c14[i], index_c14[j]] <- 0
#                   Sigma[index_c14[j], index_c14[i]] <- 0
#                 }
#               }
#
#               ## init off diag for UTh (only diag)
#               for (i in 1:(length(order_uth) - 1)) {
#                 for (j in (i + 1):length(order_uth)) {
#                   Sigma[index_uth[i], index_uth[j]] <- 0
#                   Sigma[index_uth[j], index_uth[i]] <- 0
#                 }
#               }
#
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
#                   err[i] <- interp.lin(Atilde[index_c14[i]], xTableauCalib, zTableauCalib)
#                   sigma_C14[i] <-  (pow(sigma[i],2)+pow(err[i],2))
#                   Sigma[index_c14[i], index_c14[i]] <-  sigma_C14[i]
#               }
#
#               ### UTh filling
#               for (i in order_uth) {
#                 mu[index_uth[i]] <- Atilde[index_uth[i]] #translated mean -0.075
#                 sigma_UTh[i] <- pow(sigma_uth[i], 2) #variance
#                 Sigma[index_uth[i], index_uth[i]] <- sigma_UTh[i]
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
#                   Atilde[index_osl[i]] <- exp(u[index_osl[i]] * log(xbound[2*index_osl[i]] / xbound[2*index_osl[i]-1])) * xbound[2*index_osl[i]-1]
#                   A[index_osl[i]] = Atilde[index_osl[i]] -0.075
#               }
#
#               for (i in order_c14) {
#                 u[index_c14[i]] ~dunif(0,1)
#                 Atilde[index_c14[i]] = u[index_c14[i]] * (xbound[2*index_c14[i]] - xbound[2*index_c14[i]-1]) + xbound[2*index_c14[i]-1]
#                 A[index_c14[i]] = Atilde[index_c14[i]]
#               }
#
#               for (i in order_uth) {
#                 u[index_uth[i]] ~ dunif(0,1)
#                 Atilde[index_uth[i]] <-  u[index_uth[i]] * (xbound[2* index_uth[i]] - xbound[2* index_uth[i]-1]) + xbound[2* index_uth[i]-1]
#                 A[index_uth[i]] <- Atilde[index_uth[i]] - 0.075
#               }
#
#
#
#
# }  "
#
#
#
# Unconstrained_OSL_C14 <- "model {
#               ###### Likelihood ####
#               ## Initialize Sigma matrix with zeros (off-diagonal for mixed types)
#               for (i in index_osl) {
#                   for (j in index_c14) {
#                       Sigma[i, j] <- 0
#                       Sigma[j,i] <- 0
#                   }
#               }
#               ## init off-diag for C14 (only diag)
#               for (i in 1:(length(order_c14) - 1)) {
#                 for (j in (i + 1):length(order_c14)) {
#                   Sigma[index_c14[i], index_c14[j]] <- 0
#                   Sigma[index_c14[j], index_c14[i]] <- 0
#                 }
#               }
#
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
#                   err[i] <- interp.lin(Atilde[index_c14[i]], xTableauCalib, zTableauCalib)
#                   sigma_C14[i] <-  (pow(sigma[i],2)+pow(err[i],2))
#                   Sigma[index_c14[i], index_c14[i]] <-  sigma_C14[i]
#               }
#
#
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
#                   Atilde[index_osl[i]] <- exp(u[index_osl[i]] * log(xbound[2*index_osl[i]] / xbound[2*index_osl[i]-1])) * xbound[2*index_osl[i]-1]
#                   A[index_osl[i]] = Atilde[index_osl[i]] -0.075
#               }
#
#               for (i in order_c14) {
#                 u[index_c14[i]] ~dunif(0,1)
#                 Atilde[index_c14[i]] = u[index_c14[i]] * (xbound[2*index_c14[i]] - xbound[2*index_c14[i]-1]) + xbound[2*index_c14[i]-1]
#                 A[index_c14[i]] = Atilde[index_c14[i]]
#               }
#
#
#
#
# }  "
#
# Unconstrained_OSL_UTh <- "model {
#               ###### Likelihood ####
#               ## Initialize Sigma matrix with zeros (off-diagonal for mixed types)
#               for (i in index_osl) {
#                       for (k in index_uth) {
#                       Sigma[i, k] <- 0
#                       Sigma[k,i] <- 0
#                       }
#               }
#
#               ## init off diag for UTh (only diag)
#               for (i in 1:(length(order_uth) - 1)) {
#                 for (j in (i + 1):length(order_uth)) {
#                   Sigma[index_uth[i], index_uth[j]] <- 0
#                   Sigma[index_uth[j], index_uth[i]] <- 0
#                 }
#               }
#
#
#               ## OSL filling
#               for (i in order_osl) {
#                   for (j in order_osl){
#                     Sigma[index_osl[i], index_osl[j]] <-  A[index_osl[i]]*A[index_osl[j]] * Theta[index_osl[i], index_osl[j]]
#                   }
#                 mu[index_osl[i]] <-  A[index_osl[i]] * ddot[i]
#               }
#
#               ### UTh filling
#               for (i in order_uth) {
#                 mu[index_uth[i]] <- A[index_uth[i]]
#                 sigma_UTh[i] <- pow(sigma_uth[i], 2) #variance
#                 Sigma[index_uth[i], index_uth[i]] <- sigma_UTh[i]
#               }
#
#               invSigma  <-  inverse(Sigma + encoded_covD)
#               M ~ dmnorm(mu, invSigma)
#
#               ### prior
#               for (i in order_osl) {
#                   u[index_osl[i]] ~ dunif(0,1)
#                   A[index_osl[i]] <- exp(u[index_osl[i]] * log(xbound[2*index_osl[i]] / xbound[2*index_osl[i]-1])) * xbound[2*index_osl[i]-1]
#               }
#
#               for (i in order_uth) {
#                 u[index_uth[i]] ~ dunif(0,1)
#                 A[index_uth[i]] <-  u[index_uth[i]] * (xbound[2* index_uth[i]] - xbound[2* index_uth[i]-1]) + xbound[2* index_uth[i]-1]
#               }
#
#
# }  "
#
#
# Unconstrained_C14_UTh <- "model {
#               ###### Likelihood ####
#
#
#               ###C14 filling
#               for (i in order_c14) {
#                   mu[index_c14[i]] <-  interp.lin(Atilde[index_c14[i]], xTableauCalib, yTableauCalib)
#                   err[i] <- interp.lin(Atilde[index_c14[i]], xTableauCalib, zTableauCalib)
#                   precision_C14[i] <-  1/ (pow(sigma[i],2)+pow(err[i],2))
#                   M[index_c14[i]] ~ dnorm(mu[index_c14[i]], precision_C14[i])
#               }
#
#               ### UTh filling
#               for (i in order_uth) {
#                 mu[index_uth[i]] <- Atilde[index_uth[i]]
#                 precision_uth[i] <- 1/ pow(sigma_uth[i], 2) #variance
#                 M[index_uth[i]] ~ dnorm(mu[index_uth[i]], precision_uth[i])
#               }
#
#
#
#
#               ### prior
#
#               for (i in order_c14) {
#                 u[index_c14[i]] ~ dunif(0,1)
#                 Atilde[index_c14[i]] = u[index_c14[i]] * (xbound[2*index_c14[i]] - xbound[2*index_c14[i]-1]) + xbound[2*index_c14[i]-1]
#                 A[index_c14[i]] = Atilde[index_c14[i]]
#               }
#
#               for (i in order_uth) {
#                 u[index_uth[i]] ~ dunif(0,1)
#                 Atilde[index_uth[i]] <-  u[index_uth[i]] * (xbound[2* index_uth[i]] - xbound[2* index_uth[i]-1]) + xbound[2* index_uth[i]-1]
#                 A[index_uth[i]] <- Atilde[index_uth[i]] - 0.075 # origin 1950
#               }
#
#
# }  "

#

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
# #
# Model_DC14$Unconstrained_C14_UTh <- Unconstrained_C14_UTh
# Model_DC14$Unconstrained_OSL_UTh <- Unconstrained_OSL_UTh
# Model_DC14$Unconstrained_OSL_C14 <- Unconstrained_OSL_C14
# Model_DC14$Unconstrained_OSL_C14_UTh <- Unconstrained_OSL_C14_UTh


# Model_DC14$Constrained <- Constrained
# Model_DC14$NichollsJones <- StrictNichollsJones
# #


"Model_DC14"





