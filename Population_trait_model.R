# Packages required for ODE and network metrics
library(deSolve)

Urand<-function(mu_y,var_y,m,n){
  a<-mu_y
  b<-sqrt((3*var_y)/(mu_y**2))

  Urand<-a*(1+b*(2*matrix(runif(m*n,0,1),m,n)-1))
  return(Urand)
}
varP<-0.1
varA<-0.01
M=3# No. of plant
N=1# No. of animal
#set.seed(123)
{xP0=runif(M,1,1.5)
  xA0=runif(N,1,3.5)
  zP0=runif(M,0.5,1.5)
  zA0=runif(N,1,3.5)
  xm=2
  ym=3
  sd_k=exp(-1)
  sd_c=exp(-5)

  #Maximum carrying capacity
  K=10
  # Trait elongation cost
  lambdaP=0.08
  lambdaA=0.08
  # Intrinsic growth rate
  rP=matrix(runif(M,0,1), nr=M)
  rA=matrix(runif(N,0,1), nr=N)
  # floral resource production rate
  alpha=1*matrix(runif(M,0,1), nr=M)
  #floral resource decay rate
  w=1*matrix(runif(M,0,1), nr=M)

  mu=0.1
  s=0.5
  # Benefit scaling parameter
  c=1
  # Strength of assortative
  a=0.3

  #Initial foraging effort
  bet0<-matrix(1/M,M,N)
  #Initial density
  XP0<-matrix(runif(M,0,1), nr=M)
  XA0<-matrix(runif(N,0,1), nr=N)
  #Animal adaptation rate
  G=20*runif(N,1,2)
  #Initial densities
  X0=rbind(XP0,XA0)}

# Mutualistic model function
lotka<-function(t,y,parameters){
  with(as.list(c(y,parameters)),{
    Xx<-y
    XP<-Xx[1:M]
    XA<-Xx[(M+1):(M+N)]
    bet<-matrix(Xx[(M+N+1):(M*N+M+N)],M,N)
    xP<-Xx[(M*N+M+N+1):(M*N+2*M+N)]
    zP<-Xx[(M*N+2*M+N+1):(M*N+3*M+N)]
    xA<-Xx[(M*N+3*M+N+1):(M*N+3*M+2*N)]
    zA<-Xx[(M*N+3*M+2*N+1):(M*N+3*M+3*N)]


    kP=K*exp(-((xP-xm)^2)/(2*(sd_k)^2))
    kA=K*exp(-((xA-ym)^2)/(2*(sd_k)^2))
    disMatP=c(matrix(xP,M,1))-matrix(xP,M,M, byrow = T)
    disMatA=c(matrix(xA,N,1))-matrix(xA,N,N, byrow = T)
    disMat=c(matrix(zP,M,1))-matrix(zA,M,N, byrow = T)
    cP=1*exp(-(disMatP^2)/(2*(sd_c)^2))
    cA=.1*exp(-(disMatA^2)/(2*(sd_c)^2))

    sigma_P=c*((1/(1+exp(-a*disMat)))-(exp(lambdaP*zP)-1))
    sigma_A=c*((1/(1+exp(a*disMat)))-matrix((exp(lambdaA*zA)-1),M,N, byrow = TRUE))



    dXP<-(rP-(rP*(cP%*%XP))/kP+(((sigma_P*bet)%*%XA)*(alpha))/(w+((sigma_A*bet)%*%XA)))

    dXA<-(rA-(rA*(cA%*%XA))/kA+(t(sigma_A*bet)%*%((alpha*XP)/(w+((sigma_A*bet)%*%XA)))))

    dX<-(Xx[1:(M+N)]*rbind(dXP,dXA))

    dbet<-((bet%*%diag(G,ncol = N))*sweep((sigma_A*c(alpha)*c(XP))/c(w+((sigma_A*bet)%*%XA)),2,
                                    colSums((bet*sigma_A*c(alpha)*c(XP))/c(w+((sigma_A*bet)%*%XA))),"-"))


    #  Selection gradient of plant competitive trait
    gx=c(rP*XP/kP)*rowSums(cP*(-(xP-xm)/sd_k^2+disMatP/sd_c^2))

    # Selection gradient of plant mutualistic trait
    g_zi=(a*alpha*(((bet^2)*(sigma_P)*(sigma_A^2)*exp(a*disMat))%*%(XA^2)))/((((bet*sigma_A)%*%XA)+w)^2)+
      ((c(alpha)*bet*(a*sigma_P^2*exp(-a*disMat)-lambdaP*exp(lambdaP*zP)))%*%XA)/(((bet*sigma_A)%*%XA)+w)

    # Selection gradient of animal competitive trait
    gy=colSums((cA*(matrix((-(xA-ym)/sd_k^2),N,N,byrow = TRUE)+disMatA/sd_c^2))%*%diag(c((rA*XA/kA)),ncol = N))

    # Selection gradient of animal mutualistic trait
    #g_zj=colSums((c(alpha*XP)*bet*(a*(sigma_A^2)*exp(a*disMat)-matrix(lambdaA*exp(lambdaA*zA),M,N, byrow = TRUE)))/c(((bet*sigma_A)%*%XA)+w))

    g_zj=colSums((((c(-alpha*XP)*(bet^2)*(sigma_A)*(a*sigma_A^2*exp(a*disMat)-matrix(lambdaA*exp(lambdaA*zA),M,N, byrow = TRUE)))%*%diag(XA,ncol = N))/c(((bet*sigma_A)%*%XA)+w)^2)+
                   (c(alpha*XP)*bet*(a*(sigma_A^2)*exp(a*disMat)-matrix(lambdaA*exp(lambdaA*zA),M,N, byrow = TRUE)))/c(((bet*sigma_A)%*%XA)+w))

    dxP=0.5*mu*s^2*XP*gx
    dzP=0.5*mu*s^2*XP*g_zi
    dxA=0.5*mu*s^2*XA*gy
    dzA=0.5*mu*s^2*XA*g_zj

    return(list(c(dX,dbet,dxP,dzP,dxA,dzA)))
  })
}


# Simulation of the model equation
{yini=c(c(X0),c(bet0),xP0,zP0,xA0,zA0)
times=seq(0,50000,1)
parameters=list(rP=rP,rA=rA,alpha=alpha,w=w, G=G,c=c, a=a)
solution<-ode(y=yini, times=times, func=lotka, parms=parameters)}

#Extraction of each state variables

{X<-as.matrix(solution[,2:(M+N+1)])
  ForEffMatA<-as.matrix(solution[,(M+N+2):(M*N+M+N+1)])
  x_trait<-as.matrix(solution[,(M*N+M+N+2):(M*N+2*M+N+1)])
  zP_trait<-as.matrix(solution[,(M*N+2*M+N+2):(M*N+3*M+N+1)])
  y_trait<-as.matrix(solution[,(M*N+3*M+N+2):(M*N+3*M+2*N+1)])
  zA_trait<-as.matrix(solution[,(M*N+3*M+2*N+2):(M*N+3*M+3*N+1)])
}

tail(X,5)

layout(matrix(1:8, ncol = 2), widths = 1, heights = c(1,0.85,1), respect = FALSE)

{
  par(mar = c(0.5, 4.5, 4.1, 2.0))

  matplot(ForEffMatA, type = "l", lwd=2,lty ="solid" ,pch = 1, col = rep(1:N, each=M),
          main=NA ,ylab = "Foraging effort ",xlab=NULL,xaxt="n",
          cex.lab=2.0,cex.axis=2.0)
  par(mar = c(0.5, 4.5, 1, 2.0))
  matplot(X[,1:M], type = "l",lwd=2,lty = "solid" , pch = 1, col=1:M,
          main=NA,
          ylab = "Plant density", xlab = NA,xaxt="n",cex.lab=2.0,cex.axis=2.0)

  par(mar = c(0.5, 4.5, 1, 2.0))
  matplot(x_trait, type = "l",lwd=2,lty = "solid" , pch = 1, col=1:M,
          main=NA,
          ylab = "x trait", xlab = NA,xaxt="n",cex.lab=2.0,cex.axis=2.0)
  par(mar = c(4.1, 4.5, 1, 2.0))
  matplot(zP_trait, lwd=2, type = "l",lty = "solid" , pch = 1, col=1:M,
          main=NA, ylab = "zP trait", xlab ="Time", cex.lab=2.0,cex.axis=2.0)

  net<-matrix(tail(ForEffMatA,1),M,N)*(tail(X,1)[1:M]%*%t(tail(X,1)[(M+1):(M+N)]))
  net[net<0.00001]<-0
  bipartite::plotweb(net,method='normal')
  par(mar = c(0.5, 4.5, 1, 2.0))
  matplot(X[,(M+1):(M+N)], type = "l",lwd=2,lty = "solid" , pch=1,col = 1:N,
          main=NA, ylab = "Animal density", xlab = NA,xaxt="n",cex.lab=2.0,cex.axis=2.0)

  par(mar = c(0.5, 4.5, 1, 2.0))
  matplot(y_trait, type = "l",lwd=2,lty = "solid" , pch=1,col = 1:N,
          main=NA, ylab = "y trait", xlab = "Time",xaxt="n",cex.lab=2.0,cex.axis=2.0)

  par(mar = c(4.1, 4.5, 1, 2.1))
  matplot(zA_trait, type = "l", lwd=2,lty ="solid" ,pch = 1, col = 1:N,
          main=NA ,ylab = "zA trait ",xlab="Time",
          cex.lab=2.0,cex.axis=2.0)
}



tail(X,5)
tail(x_trait,1)
tail(zP_trait,1)
tail(y_trait,1)
tail(zA_trait,1)


matrix(tail(ForEffMatA,1),M,N)

net

dev.copy(jpeg,"trait_mediated.jpeg",width = 300, height = 300,units = "mm", res = 400)

dev.off()





zA=zP=seq.int(1,10,1)
y1=exp(-(x^2)/(2*(exp(1))^2))
y2=1/(1+exp(0.6*x))
plot(x,y1,type = "l",lty = "solid" ,lwd=2, pch = 1)
plot(x,y2,type = "l",lty = "solid" ,lwd=2, pch = 1,col="blue")


x=seq(-5,5)
y=1/(1+exp(0.7*x))
plot(x,y,type = "l",lty = "solid" ,lwd=2, pch = 1)


disMat=c(matrix(zP,M,1))-matrix(zA,M,N, byrow = T)

sigma_P=c*((1/(1+exp(-a*disMat)))-(exp(0.1*zP)-1))
sigma_A=c*((1/(1+exp(a*disMat)))-matrix((exp(0.1*zA)-1),M,N, byrow = TRUE))
plot(disMat,sigma_P,col=1:10)
plot(disMat,sigma_A,col=1:10)



{xP0=c(1.56827,1.15267)
  xA0=c(2.6836,1.38599)
  zP0=c(0.82854,1.23288)
  zA0=c(2.43836,2.23784)
  xm=2
  ym=3
  sd_k=exp(-1)
  sd_c=exp(-2)

  #Maximum carrying capacity
  K=100
  # Trait elongation cost
  lambdaP=0.05
  lambdaA=0.05
  # Intrinsic growth rate
  rP=c(0.767,0.3717)
  rA=c(0.3111,0.3437)
  # floral resource production rate
  alpha=c(1,1)
  #floral resource decay rate
  w=c(1,1)

  mu=0.01
  s=0.05
  # Benefit scaling parameter
  c=.3
  # Strength of assortative
  a=0.3

  #Initial foraging effort
  bet0<-matrix(1/M,M,N)
  #Initial density
  XP0<-matrix(c(0.6761,0.727), nr=M)
  XA0<-matrix(c(0.978,0.1953), nr=N)
  #Animal adaptation rate
  G=matrix(c(10,20))
  #Initial densities
  X0=rbind(XP0,XA0)}

{xP0=runif(M,1,2.5)
  xA0=runif(N,1,3.5)
  zP0=runif(M,1,2.5)
  zA0=runif(N,1,3.5)
  xm=2
  ym=3
  sd_k=exp(-1)
  sd_c=exp(-2)

  # Intrinsic growth rate
  rP=matrix(runif(M,0,1), nr=M)
  rA=matrix(runif(N,0,1), nr=N)
  # floral resource production rate
  alpha=10*matrix(runif(M,0,1), nr=M)
  #floral resource decay rate
  w=10*matrix(runif(M,0,.1), nr=M)

  mu=0.01
  s=0.05
  #conversion rate
  c=30

  a=0.5

  #Initial foraging effort
  bet0<-matrix(1/M,M,N)
  # Initial floral resource abundance
  #Fi0<-matrix(runif(M,0,1), nr=M)
  #Initial density
  XP0<-matrix(runif(M,0,1), nr=M)
  XA0<-matrix(runif(N,0,1), nr=N)
  #Animal adaptation rate
  G=30*matrix(runif(N,1,2), nr=N)
  #Initial densities
  X0=rbind(XP0,XA0)}




gx=-(rP*XP*(xP-xm))/((sd_k^2)*kP)
g_zi=(a*alpha*(((bet^2)*(sigma_P)*(sigma_A^2)*exp(a*disMat))%*%(XA^2)))/((((bet*sigma_A)%*%XA)+w)^2)+
  (a*alpha*((bet*(sigma_P^2)*exp(-a*disMat))%*%XA))/(((bet*sigma_A)%*%XA)+w)
gy=-(rA*XA*(xA-ym))/((sd_k^2)*kA)
g_zj=colSums((((c(-a*alpha*XP)*(bet^2)*(sigma_A^3)*exp(a*disMat))%*%(XA))/c(((bet*sigma_A)%*%XA)+w)^2)+
               (c(a*alpha*XP)*bet*(sigma_A^2)*exp(a*disMat))/c(((bet*sigma_A)%*%XA)+w))





kP=100*exp(-((xP0-xm)^2)/(2*(sd_k)^2))
kA=100*exp(-((xA0-ym)^2)/(2*(sd_k)^2))
disMatP=c(matrix(xP0,M,1))-matrix(xP0,M,M, byrow = T)
disMatA=c(matrix(xA0,N,1))-matrix(xA0,N,N, byrow = T)
disMat=c(matrix(zP0,M,1))-matrix(zA0,M,N, byrow = T)
cP=1*exp(-(disMatP^2)/(2*(sd_c)^2))
cA=1*exp(-(disMatA^2)/(2*(sd_c)^2))

# sigma_P=c*((1/(1+exp(-a*disMat))))
# sigma_A=c*((1/(1+exp(a*disMat))))
sigma_P=c*((1/(1+exp(-a*disMat)))-(exp(lambdaP*zP0)-1))
sigma_A=c*((1/(1+exp(a*disMat)))-matrix((exp(lambdaA*zA0)-1),M,N, byrow = TRUE))


(c(rP*XP0/kP)*rowSums(cP))*(-(xP0-xm)/sd_k^2+rowSums(disMatP)/sd_c^2)
