#Pregunta 1

A1<-cbind(c(.7,0,.9),c(.1,.4,0),c(0,.1,.8))
A2<-cbind(c(-.2,0,0),c(0,.1,0),c(0,.1,0))

Sigmau<-cbind(c(.26,.03,0),c(.03,.09,0),c(0,0,.81))
SIGMAU<-rbind(cbind(Sigmau,diag(3)-diag(3)),cbind(diag(3)-diag(3),diag(3)-diag(3)))
VecU<-cbind(c(c(.26,.03,0),c(0,0,0),c(.03,.09,0),c(0,0,0),c(0,0,.81),c(0,0,0),diag(3)-diag(3),diag(3)-diag(3)))

A<-rbind(cbind(A1,A2),cbind(diag(3),diag(3)-diag(3)))

VecGamma0<-solve((diag(36)-kronecker(A,A)))%*%VecU

Gamma0<-cbind(VecGamma0[1:3],VecGamma0[7:9],VecGamma0[13:15])
Gamma1<-cbind(VecGamma0[19:21],VecGamma0[25:27],VecGamma0[31:33])

Gamma2<-A1%*%Gamma1+A2%*%Gamma0
Gamma3<-A1%*%Gamma2+A2%*%Gamma1

D<-diag(c(sqrt(Gamma0[1,1]),sqrt(Gamma0[2,2]),sqrt(Gamma0[3,3])))
invD<-solve(D)

R0<-invD%*%Gamma0%*%invD
R1<-invD%*%Gamma1%*%invD
R2<-invD%*%Gamma2%*%invD
R3<-invD%*%Gamma3%*%invD

par(mfrow=c(3,3),title='Correlaciones')

plot(c(0,1,2,3),c(R0[1,1],R1[1,1],R2[1,1],R3[1,1]),type='overplotted',col='4',xlab="h",ylab="Corr(GNP_t,GNP_{t-h})")
plot(c(0,1,2,3),c(R0[1,2],R1[1,2],R2[1,2],R3[1,2]),type='overplotted',col='4',xlab="h",ylab="Corr(GNP_t,M2_{t-h})")
plot(c(0,1,2,3),c(R0[1,3],R1[1,3],R2[1,3],R3[1,3]),type='overplotted',col='4',xlab="h",ylab="Corr(GNP_t,IR_{t-h})")
plot(c(0,1,2,3),c(R0[2,1],R1[2,1],R2[2,1],R3[2,1]),type='overplotted',col='4',xlab="h",ylab="Corr(M2_t,GNP_{t-h})")
plot(c(0,1,2,3),c(R0[2,2],R1[2,2],R2[2,2],R3[2,2]),type='overplotted',col='4',xlab="h",ylab="Corr(M2_t,M2_{t-h})")
plot(c(0,1,2,3),c(R0[2,3],R1[2,3],R2[2,3],R3[2,3]),type='overplotted',col='4',xlab="h",ylab="Corr(M2_t,IR_{t-h})")
plot(c(0,1,2,3),c(R0[3,1],R1[3,1],R2[3,1],R3[3,1]),type='overplotted',col='4',xlab="h",ylab="Corr(IR_t,GNP_{t-h})")
plot(c(0,1,2,3),c(R0[3,2],R1[3,2],R2[3,2],R3[3,2]),type='overplotted',col='4',xlab="h",ylab="Corr(IR_t,M2_{t-h})")
plot(c(0,1,2,3),c(R0[3,3],R1[3,3],R2[3,3],R3[3,3]),type='overplotted',col='4',xlab="h",ylab="Corr(IR_t,IR_{t-h})")

#Pregunta 2

#a)

x2000<-c(.7,1,1.5)
x1999<-c(1,1.5,3)

nu<-c(2,1,0)
mu<-solve((diag(3)-A1-A2))%*%nu

y2000<-x2000-mu
y1999<-x1999-mu

Y2000<-cbind(c(y2000,y1999))

Y2001<-A%*%Y2000
Y2002<-A%*%Y2001
Y2003<-A%*%Y2002

#b)

J<-cbind(diag(3),diag(3)-diag(3))

phi0<-diag(3)
phi1<-J%*%A%*%t(J)
phi2<-J%*%A%*%A%*%t(J)

MSE1<-phi0%*%Sigmau%*%t(phi0)
MSE2<- MSE1+phi1%*%Sigmau%*%t(phi1)
MSE3<- MSE2+phi2%*%Sigmau%*%t(phi2)

#c

#Para t=2001 al 90%

Y2001[1]+qnorm(.05)*sqrt(MSE1[1,1])+mu[1]
Y2001[1]-qnorm(.05)*sqrt(MSE1[1,1])+mu[1]

Y2001[2]+qnorm(.05)*sqrt(MSE1[2,2])+mu[2]
Y2001[2]-qnorm(.05)*sqrt(MSE1[2,2])+mu[2]

Y2001[3]+qnorm(.05)*sqrt(MSE1[3,3])+mu[3]
Y2001[3]-qnorm(.05)*sqrt(MSE1[3,3])+mu[3]

#Para t=2002 al 90%

Y2002[1]+qnorm(.05)*sqrt(MSE2[1,1])+mu[1]
Y2002[1]-qnorm(.05)*sqrt(MSE2[1,1])+mu[1]

Y2002[2]+qnorm(.05)*sqrt(MSE2[2,2])+mu[2]
Y2002[2]-qnorm(.05)*sqrt(MSE2[2,2])+mu[2]

Y2002[3]+qnorm(.05)*sqrt(MSE2[3,3])+mu[3]
Y2002[3]-qnorm(.05)*sqrt(MSE2[3,3])+mu[3]

#Para t=2003 al 90%

Y2003[1]+qnorm(.05)*sqrt(MSE3[1,1])+mu[1]
Y2003[1]-qnorm(.05)*sqrt(MSE3[1,1])+mu[1]

Y2003[2]+qnorm(.05)*sqrt(MSE3[2,2])+mu[2]
Y2003[2]-qnorm(.05)*sqrt(MSE3[2,2])+mu[2]

Y2003[3]+qnorm(.05)*sqrt(MSE3[3,3])+mu[3]
Y2003[3]-qnorm(.05)*sqrt(MSE3[3,3])+mu[3]

#Para t=2001 al 95%

Y2001[1]+qnorm(.025)*sqrt(MSE1[1,1])+mu[1]
Y2001[1]-qnorm(.025)*sqrt(MSE1[1,1])+mu[1]

Y2001[2]+qnorm(.025)*sqrt(MSE1[2,2])+mu[2]
Y2001[2]-qnorm(.025)*sqrt(MSE1[2,2])+mu[2]

Y2001[3]+qnorm(.025)*sqrt(MSE1[3,3])+mu[3]
Y2001[3]-qnorm(.025)*sqrt(MSE1[3,3])+mu[3]

#Para t=2002 al 95%

Y2002[1]+qnorm(.025)*sqrt(MSE2[1,1])+mu[1]
Y2002[1]-qnorm(.025)*sqrt(MSE2[1,1])+mu[1]

Y2002[2]+qnorm(.025)*sqrt(MSE2[2,2])+mu[2]
Y2002[2]-qnorm(.025)*sqrt(MSE2[2,2])+mu[2]

Y2002[3]+qnorm(.025)*sqrt(MSE2[3,3])+mu[3]
Y2002[3]-qnorm(.025)*sqrt(MSE2[3,3])+mu[3]

#Para t=2003 al 95%

Y2003[1]+qnorm(.025)*sqrt(MSE3[1,1])+mu[1]
Y2003[1]-qnorm(.025)*sqrt(MSE3[1,1])+mu[1]

Y2003[2]+qnorm(.025)*sqrt(MSE3[2,2])+mu[2]
Y2003[2]-qnorm(.025)*sqrt(MSE3[2,2])+mu[2]

Y2003[3]+qnorm(.025)*sqrt(MSE3[3,3])+mu[3]
Y2003[3]-qnorm(.025)*sqrt(MSE3[3,3])+mu[3]

#d)

#Como alpha=0.03, debemos crear intervalos de confianza del 99% para cada GNP_k

Y2001[1]+qnorm(.005)*sqrt(MSE1[1,1])+mu[1]
Y2001[1]-qnorm(.005)*sqrt(MSE1[1,1])+mu[1]

Y2002[1]+qnorm(.005)*sqrt(MSE2[1,1])+mu[1]
Y2002[1]-qnorm(.005)*sqrt(MSE2[1,1])+mu[1]

Y2003[1]+qnorm(.005)*sqrt(MSE3[1,1])+mu[1]
Y2003[1]-qnorm(.005)*sqrt(MSE3[1,1])+mu[1]