

Model_DC14 <- "model {
              ###### Likelihood ####
              ## Initialize Sigma matrix with zeros (off-diagonal for mixed types)
              for (i in 1:I) {
                  for (j in 1:I) {
                      Sigma[i, j] <- 0
                  }
              }
              ## OSL filling
              for (i in order_osl) {
                  for (j in order_osl){
                    Sigma[index_osl[i], index_osl[j]] <-  A[index_osl[i]]*A[index_osl[j]] * Theta[index_osl[i], index_osl[j]]
                  }
                mu[index_osl[i]] <-  A[index_osl[i]] * ddot[i]
              }
              ###C14 filling
              for (i in order_c14) {
                  mu[index_c14[i]] <-  interp.lin(A[index_c14[i]], xTableauCalib, yTableauCalib)
                  Z[i] ~ dcat(c(0.1, 0.9))
                  err[i] <- interp.lin(A[index_c14[i]], xTableauCalib, zTableauCalib)
                  sigma_C14[i] <-  (alpha[i]^(-Z[i]+2)*(pow(sigma[i],2)+pow(err[i],2)))
                  Sigma[index_c14[i], index_c14[i]] <-  sigma_C14[i]
              }

              invSigma  <-  inverse(Sigma + encoded_covD)
              M ~ dmnorm(mu, invSigma)


              #### prior ###
              ## starting with a OSL
              u[1] ~ dunif(0,1)
              if (equals(encoding[1], 1)) {
                A[1] <- exp(u[1] * log(x_bound[2] / x_bound[1])) + x_bound[1]
                Atilde[1]  <-  A[1]-0.075
              }
              else {
              A[1] <- u[1] * (x_bound[2]-x_bound[1]) + x_bound[1]
              invalpha[1] <-  dgamma(3,4)
              alpha[1] <- 1/invalpha[1]
              Atilde[1] <- A[1]
              }

              for (i in 2:Nb_sample) {
                  u[i] ~ dunif(0,1)
                  if (equals(encoding[i], 1)) {
                      CS[i] <- max(StratiConstraints[1:i, i] * c(x_bound[2*(i-1) +1], A[1:(i-1)] * encoding[1:(i-1)] + (A[1:(i-1)] + 0.075) * (1-encoding[1:(i-1)])))
                      A[i] <- exp(u[i] * log(x_bound[2*i] / CS[i])) + CS[i]
                      Atilde[i] <-  A[i] - 0.075
                  }

                  else {
                      CS[i] <- max(StratiConstraints[1:i, i] * c(x_bound[2*(i-1) +1], A[1:(i-1)] * (1-encoding[1:(i-1)]) + (Atilde[1:(i-1)]) * encoding[1:(i-1)] ))
                      A[i] <- u[i] * (x_bound[2*i] - CS[i]) + CS[i]
                      Atilde[i] = A[i]
                      invalpha[CS_osl[i]] ~ dgamma(3,4)
                      alpha <- 1/inv_alpha[i]
                  }


              }



}  "



