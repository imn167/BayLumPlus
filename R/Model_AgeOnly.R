
### Model with Jeffreys prior
ModelAgePrior <- list()

Jeffreys <- " model {
  ###### Likelyhood ####
  for (i in 1:I) {
    for (j in 1:I) {
      Sigma[i, j] = A[i]* A[j] * Theta[i,j]
    }
    mu[i] <— A[i] * ddot[i]
  }
  invSigma <- inverse(Sigma) #DP and symmetric ?
  D ~ dnorm(mu, invSigma)

  ###### Prior #####
  ## borne inf de la strati, soit la + proche de 0
  Atemp[1]=xbound[1]
  # i0=2
  u[1]~dunif(0,1)
  CS[1]=xbound[1]
  # alpha1 et beta1 : borne pour l'age 1
  # la loi de Atem2 est proportionnel à celle de jeffrey sur l'intervalle (alpha1, beta1)
  Atemp[2]<-exp( u[1]*log(xbound[2]/CS[1]) +log(CS[1]) )
  # i0>2
  for(i0 in 3:(I+1)){
    u[i0-1]~dunif(0,1)
    #On choisit le marjorant des bornes inf soit entre alpha_(i-1) et A1 ..... A(i-1)
    CS[i0-1]=max( StratiConstraints[(1:(i0-1)),(i0-1)] *c( xbound[(2*(i0-1)-1)],Atemp[2:(i0-1)]) )
    Atemp[i0]<-exp(u[(i0-1)]*log(xbound[2*(i0-1)]/CS[i0-1])+log(CS[i0-1]))
  }
  A=Atemp[2:(I+1)] #A1, A2, ...., An

}"

## saving the model to data in rda form


ModelAgePrior$Jeffreys <- Jeffreys ## see how to save the Model in the data of the package ?
usethis::use_data(ModelAgePrior, overwrite = T)
