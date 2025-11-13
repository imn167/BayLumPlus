#
#
# cat(ModelPrior$C14)
# ModelPrior$Sample1_C14 <- "
#                           u[1] ~ dunif(0, 1)
#                           A[1] = u[1] * (xbound[2]-xbound[1]) + xbound[1]
#                           invalpha[1]~dgamma(3,4)
#                           alpha[1]<-1/invalpha[1]
#                           Atilde[1] = A[1]
#                           "
#
# ModelPrior$Sample1_OSL <- "
#                 u[1]~dunif(0,1)
#                 CS[1]<-xbound[1]
#                 A[1]<-exp(u[1]*log(xbound[2]/CS[1])+log(CS[1]))
#                 Atilde[1] = A[1] -0.75
#                 "
#
# ModelPrior$OSL_C14 <- "
#             for(k in 1:q){
#                   # donnee OSL
#                   for(i in (ind_change[2*k-1]+1):ind_change[2*k]){
#                   u[i]~dunif(0,1)
#                   CS[i]<-max(StratiConstraints[(1:i),i]*c(xbound[(2*i-1)],A[1:(i-1)]))
#                   A[i]<-exp(u[i]*log(xbound[2*(i-1)]/CS[i])+log(CS[i]))
#                   Atilde[i] = A[i]-0.75
#                   }
#                   # donnee C14
#                   for(i in (ind_change[2*k]+1):ind_change[2*k+1]){
#                   CS[i]<-max(StratiConstraints[(1:i),i]*c(xbound[(2*(i-1)+1)],Atilde[1:(i-1)]))
#                   u[i] ~ dunif(0,1)
#                   A[i] = u[i]* (xbound[2*i]-CS[i]) +CS[i]
#                   invalpha[CS_C14[i]]~dgamma(3,4)
#                   alpha[CS_C14[i]]<-1/invalpha[CS_C14[i]]
#                   Atilde[i] = A[i]
#                   }
#             }"
#
# ModelPrior$C14_OSL <- "
#               for(k in 1:q){
#                   # donnee C14
#                   for(i in (ind_change[2*k-1]+1):ind_change[2*k]){
#                     CS[i]<-max(StratiConstraints[(1:i),i]*c(xbound[(2*(i-1)+1)],Atilde[1:(i-1)]))
#                     u[i] ~ dunif(0, 1)
#                     A[i] = u[i] * (xbound[2*i] - CS[i]) + CS[i]
#                     invalpha[CS_C14[i]]~dgamma(3,4)
#                     alpha[CS_C14[i]]<-1/invalpha[CS_C14[i]]
#                     Atilde[i] = A[i]
#                   }
#                   # donne OSL
#                   for(i in (ind_change[2*k]+1):ind_change[2*k+1]){
#                     u[i]~dunif(0,1)
#                     CS[i]<-max(StratiConstraints[(1:i),i]*c(xbound[(2*i-1)],A[1:(i-1)]))
#                     A[i]<-exp(u[i]*log(xbound[2*(i-1)]/CS[i])+log(CS[i]))
#                     Atilde[i] = A[i]-0.75
#                   }
#             }"
#
# ModelPrior$OSL <- "
#     # donnee OSL
#       kk<-k+1
#       for(i in (ind_change[2*q+1]+1):ind_change[2*(q+1)]){
#         u[i]~dunif(0,1)
#         CS[i]<-max(StratiConstraints[(1:i),i]*c(xbound[(2*i-1)],A[1:(i-1)]))
#         A[i]<-exp(u[i]*log(xbound[2*(i-1)]/CS[i])+log(CS[i]))
#         Atilde[i] = A[i]-0.75
#       }"
#
#
# ModelPrior$C14 <- "
#         # donnee C14
#         for(i in (ind_change[2*q+1]+1):ind_change[2*(q+1)]){
#             CS[i]<-max(StratiConstraints[(1:i),i]*c(xbound[(2*(i-1)+1)],Atilde[1:(i-1)]))
#             u[i] ~dunif(0, 1)
#             A[i] = u[i] * (xbound[2*i]- CS[i]) + CS[i]
#             invalpha[CS_C14[i]]~dgamma(3,4)
#             alpha[CS_C14[i]]<-1/invalpha[CS_C14[i]]
#             Atilde[i] = A[i]
#         }"




