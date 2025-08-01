#' @title JAGS Models for OSL Age Estimation in \code{\link{Compute_AgeS_D}}
#'
#' @description
#' JAGS models used to estimate true OSL ages based on data obtained from the Bayesian OSL analysis performed
#' by the function  \code{\link{Palaeodose_Computation}}.
#'
#' @details
#' These models take as input the estimated dose response ($D$) from  \code{\link{Palaeodose_Computation}}
#' along with the structured data matrix computed by  \code{\link{create_MeasuresDataFrame}}.
#' The models are designed to refine age estimation by integrating these measurements into a Bayesian framework.
#'@md
#' @references
#' To cite this package, please use: citation("BayLum")

#'@format BayLum's Old Age Model (wrong vector law)

# Jeffreys <- " model {
#   ###### Likelyhood ####
#   for (i in 1:I) {
#     for (j in 1:I) {
#       Sigma[i, j] = A[i]* A[j] * Theta[i,j]
#     }
#     mu[i] = A[i] * ddot[i]
#   }
#   invSigma = inverse(Sigma + covD) #DP and symmetric ?
#   D ~ dmnorm(mu, invSigma)
#
#   ###### Prior #####
#   ## borne inf de la strati, soit la + proche de 0
#   Atemp[1]=xbound[1]
#   # i0=2
#   u[1]~dunif(0,1)
#   CS[1]=xbound[1]
#   # alpha1 et beta1 : borne pour l'age 1
#   # la loi de Atem2 est proportionnel à celle de jeffrey sur l'intervalle (alpha1, beta1)
#   Atemp[2]=exp( u[1]*log(xbound[2]/CS[1]) +log(CS[1]) )
#   # i0>2
#   for(i0 in 3:(I+1)){
#     u[i0-1]~dunif(0,1)
#     #On choisit le marjorant des bornes inf soit entre alpha_(i-1) et A1 ..... A(i-1)
#     CS[i0-1]=max( StratiConstraints[(1:(i0-1)),(i0-1)] *c( xbound[(2*(i0-1)-1)],Atemp[2:(i0-1)]) )
#     Atemp[i0]=exp(u[(i0-1)]*log(xbound[2*(i0-1)]/CS[i0-1])+log(CS[i0-1]))
#   }
#   A=Atemp[2:(I+1)] #A1, A2, ...., An
#
# }"

## saving the model to data in rda form


# ModelAgePrior$Jeffreys <- Jeffreys



#'@format Jeffrey's Age Model with log-uniform order settings
# StrictOrder <- " model {
#   ###### Likelyhood ####
#   for (i in 1:I) {
#     for (j in 1:I) {
#       Sigma[i, j] = A[i]* A[j] * Theta[i,j]
#     }
#     mu[i] = A[i] * ddot[i]
#   }
#   invSigma = inverse(Sigma + covD) #DP and symmetric ?
#   D ~ dmnorm(mu, invSigma)
#
#   #bounds
#   alpha = xbound[1]
#   beta = xbound[2]
# ###### Prior #####
#   for (i in 1: (I+1)) {
#   e[i] ~ dexp(1)
#   }
#
#   for (i in 1:I) {
#   u[i] = sum(e[1:i]) / sum(e[1:(I+1)])
#   A[i] = exp(u[i] * (log(beta) - log(alpha)) + log(alpha))
#   }
#
# }"


#'@format Nicholls' Age Model applied on ages directly
# Strict_nicholls <- "model {
#   ###### Likelyhood ####
#   for (i in 1:I) {
#     for (j in 1:I) {
#       Sigma[i, j] = A[i]* A[j] * Theta[i,j]
#     }
#     mu[i] = A[i] * ddot[i]
#   }
#   invSigma = inverse(Sigma + covD) #DP and symmetric ?
#   D ~ dmnorm(mu, invSigma)
#
#   #bounds
#   alpha = xbound[1]
#   beta = xbound[2]
# ###### Prior #####
#   for (i in 1: (I-1)) {
#   e[i] ~ dexp(1)
#   }
#
#   s ~ dunif(0, 1)
#   first ~ dunif(0, (1-s))
#   u[1]= first
#   u[I]= s + first
#
#   for (i in 2:(I-1)) {
#     u[i] = (sum(e[1:(i-1)]) / sum(e[1:(I-1)])) *s + first
#   }
#   for (i in 1:I) {
#     A[i] = u[i] * (beta - alpha) + alpha
#
#   }
#
# }"

# nichollsBR <- "model {
#   ###### Likelyhood ####
#   for (i in 1:I) {
#     for (j in 1:I) {
#       Sigma[i, j] = A[i]* A[j] * Theta[i,j]
#     }
#     mu[i] = A[i] * ddot[i]
#   }
#   invSigma = inverse(Sigma + covD) #DP and symmetric ?
#   D ~ dmnorm(mu, invSigma)
#
#   #bounds
#   alpha = xbound[1]
#   beta = xbound[2]
# ###### Prior #####
#   for (i in 1: (I-1)) {
#   e[i] ~ dexp(1)
#   }
#
#   s ~ dunif(0, 1)
#   first ~ dunif(0, (1-s))
#   A[1]= first * (beta-alpha) + alpha
#   A[I]= (s + first) * (beta-alpha) + alpha
#
#   for (i in 1:(I-2)) {
#     u[i] = (sum(e[1:(i)]) / sum(e[1:(I-1)])) *s + first
#     z[i] ~ dbinom(.5, 1)
#     b[i] ~ dbeta(1, (I-2))
#   }
#   for (i in 2:(I-1)) {
#     A[i] = (u[i-1] + z[i-1]*b[i-1]) * (beta - alpha) + alpha
#   }
#
# }"

#'@format Jeffrey's Age Model with conditional setting
# Conditional <- " model {
#   ###### Likelyhood ####
#   for (i in 1:I) {
#     for (j in 1:I) {
#       Sigma[i, j] = A[i]* A[j] * Theta[i,j]
#     }
#     mu[i] = A[i] * ddot[i]
#   }
#   invSigma = inverse(Sigma + covD) #DP and symmetric
#   D ~ dmnorm(mu, invSigma)
#
#   ###### Prior #####
#   T1=xbound[1]
#   T2 = xbound[2]
#
#   #i = 1
#   u[1]~dunif(0,1)
#   A[1]=exp( log(T2)- (u[1])^(1/I) * log(T2/T1) ) # simulation ~ pi(A1)
#
#   # i>2
#   for(i in 2:I){
#     u[i]~dunif(0,1)
#
#     A[i]=exp((u[i])^(1 / (I-i+1))*log(A[i-1]/T2)+log(T2)) #simulation ~ pi(Ai | A1...A(i-1))
#   }
#
#
# }"

#'@format Independance Age Model

# Independance <-  "model {
#   #### Likelyhood #####
#   for (i in 1:I) {
#     for (j in 1:I) {
#
#       Sigma[i, j] = A[i]* A[j] * Theta[i,j]
#
#     }
#     mu[i] = A[i] * ddot[i]
#   }
#   invSigma = inverse(Sigma + covD)
#
#   D ~dmnorm(mu, invSigma)
#
#   T1 = xbound[1]
#   T2 = xbound[2]
#   ### Prior ####
#
#   for(i in 1:I) {
#     u[i] ~ dunif(0,1)
#     A[i] = exp(u[i] * log(T2/T1) + log(T1))
#   }
#
# }"


"ModelAgePrior"

#
# ModelAgePrior$StrictOrder <- StrictOrder
# ModelAgePrior$StrictNicholls <- Strict_nicholls
# ModelAgePrior$nichollsBR <- nichollsBR
# ModelAgePrior$Independance <- Independance
# usethis::use_data(ModelAgePrior, overwrite = T)
# ModelAgePrior$Conditional <- Conditional
