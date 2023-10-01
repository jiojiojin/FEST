# Function
g<-function(u){
  return(1/(1+exp(-u)))
}
a<-function(u){
  return(log(1+exp(u)))
}

data_generation_c1<-function(rho_z,N,ns,Ns,nt,Nt,p,q){
  mu<-rep(0,p+q)
  r<-rho_z^seq(0,p+q-1)
  Sigma<-toeplitz(r)
  # generate Z based on latent W
  W<-mvrnorm(N,mu,Sigma)
  Z<-floor(a(W))
  # separate X and S
  X<-Z[,1:p]
  S<-Z[,(p+1):(p+q)]
  S[,1]<-floor(a(rowSums(W[,1:5])+W[,(1+p)]))
  S[,2]<-floor(a(rowSums(W[,3:7])+W[,(2+p)]))
  S[,3]<-floor(a(rowSums(W[,5:10])+W[,(3+p)]))
  # data normalization
  X_std<-scale(X)
  S_std<-scale(S)
  # generate Y based on Z_std (m(z) is correctly specified)
  Z<-cbind(rep(1,N),X_std,S_std)
  tempy<-0.2-X_std[,1]+X_std[,2]-X_std[,3]+2*(S_std[,1]-S_std[,2]+S_std[,3])
  u<-runif(N)
  Y<-ifelse(u<g(tempy),1,0)
  # generate R to separate the source and target populations
  tempr<-X_std[,1]-X_std[,2]-X_std[,3]+S_std[,1]
  uu<-runif(N)
  R<-ifelse(uu<g(tempr),1,0)
  # separate the data to four subsets
  # first separate data from target (R=0) and source (R=1)
  Ys<-Y[R==1]
  Yt<-Y[R==0]
  Xs<-X_std[R==1,]
  Xt<-X_std[R==0,]
  Ss<-S_std[R==1,]
  St<-S_std[R==0,]
  # then create labelled and unlabelled data in the two populations 
  Ysl<-Ys[1:ns]
  Ytl<-Yt[1:nt]
  Xsl<-Xs[1:ns,]
  Xtl<-Xt[1:nt,]
  Ssl<-Ss[1:ns,]
  Stl<-St[1:nt,]
  Ysu<-Ys[(ns+1):(Ns+ns)]
  Ytu<-Yt[(nt+1):(Nt+nt)]
  Xsu<-Xs[(ns+1):(Ns+ns),]
  Xtu<-Xt[(nt+1):(Nt+nt),]
  Ssu<-Ss[(ns+1):(Ns+ns),]
  Stu<-St[(nt+1):(Nt+nt),]
  return(list(Ysl=Ysl, Xsl=Xsl, Ssl=Ssl, Ytl=Ytl, Xtl=Xtl, Stl=Stl,
              Ysu=Ysu, Xsu=Xsu, Ssu=Ssu, Ytu=Ytu, Xtu=Xtu, Stu=Stu))
}

data_generation_c2<-function(rho_z,N,ns,Ns,nt,Nt,p,q){
  mu<-rep(0,p+q)
  r<-rho_z^seq(0,p+q-1)
  Sigma<-toeplitz(r)
  # generate Z based on latent W
  W<-mvrnorm(N,mu,Sigma)
  Z<-floor(a(W))
  # separate X and S
  X<-Z[,1:p]
  S<-Z[,(p+1):(p+q)]
  S[,1]<-floor(a(rowSums(W[,1:5])+W[,(1+p)]))
  S[,2]<-floor(a(rowSums(W[,3:7])+W[,(2+p)]))
  S[,3]<-floor(a(rowSums(W[,5:10])+W[,(3+p)]))
  # data normalization
  X_std<-scale(X)
  S_std<-scale(S)
  # generate Y based on Z_std (m(z) is correctly specified)
  Z<-cbind(rep(1,N),X_std,S_std)
  tempy<--X_std[,1]+X_std[,2]-X_std[,3]^2-2*S_std[,2]+2*S_std[,3]+2*S_std[,1]/(1+exp(-X_std[,1]*S_std[,1]^2))
  u<-runif(N)
  Y<-ifelse(u<g(tempy),1,0)
  # generate R to separate the source and target populations
  tempr<-X_std[,1]-X_std[,2]-X_std[,3]+S_std[,1]
  uu<-runif(N)
  R<-ifelse(uu<g(tempr),1,0)
  # separate the data to four subsets
  # first separate data from target (R=0) and source (R=1)
  Ys<-Y[R==1]
  Yt<-Y[R==0]
  Xs<-X_std[R==1,]
  Xt<-X_std[R==0,]
  Ss<-S_std[R==1,]
  St<-S_std[R==0,]
  # then create labelled and unlabelled data in the two populations 
  Ysl<-Ys[1:ns]
  Ytl<-Yt[1:nt]
  Xsl<-Xs[1:ns,]
  Xtl<-Xt[1:nt,]
  Ssl<-Ss[1:ns,]
  Stl<-St[1:nt,]
  Ysu<-Ys[(ns+1):(Ns+ns)]
  Ytu<-Yt[(nt+1):(Nt+nt)]
  Xsu<-Xs[(ns+1):(Ns+ns),]
  Xtu<-Xt[(nt+1):(Nt+nt),]
  Ssu<-Ss[(ns+1):(Ns+ns),]
  Stu<-St[(nt+1):(Nt+nt),]
  return(list(Ysl=Ysl, Xsl=Xsl, Ssl=Ssl, Ytl=Ytl, Xtl=Xtl, Stl=Stl,
              Ysu=Ysu, Xsu=Xsu, Ssu=Ssu, Ytu=Ytu, Xtu=Xtu, Stu=Stu))
}

data_generation_c3<-function(rho_z,N,ns,Ns,nt,Nt,p,q){
  mu<-rep(0,p+q)
  r<-rho_z^seq(0,p+q-1)
  Sigma<-toeplitz(r)
  # generate Z based on latent W
  W<-mvrnorm(N,mu,Sigma)
  Z<-floor(a(W))
  # separate X and S
  X<-Z[,1:p]
  S<-Z[,(p+1):(p+q)]
  S[,1]<-floor(a(rowSums(W[,1:5])+W[,(1+p)]))
  S[,2]<-floor(a(rowSums(W[,3:7])+W[,(2+p)]))
  S[,3]<-floor(a(rowSums(W[,5:10])+W[,(3+p)]))
  # data normalization
  X_std<-scale(X)
  S_std<-scale(S)
  # generate Y based on Z_std (m(z) is correctly specified)
  Z<-cbind(rep(1,N),X_std,S_std)
  tempy<-0.2-X_std[,1]+X_std[,2]-X_std[,3]+2*(S_std[,1]-S_std[,2]+S_std[,3])
  u<-runif(N)
  Y<-ifelse(u<g(tempy),1,0)
  # generate R to separate the source and target populations
  tempr<-1.8*S_std[,2]-2*S_std[,3]+S_std[,1]*(1+(X_std[,1]+X_std[,2])/(1+exp(-rowSums(X_std[,3:10]))))
  uu<-runif(N)
  R<-ifelse(uu<g(tempr),1,0)
  # separate the data to four subsets
  # first separate data from target (R=0) and source (R=1)
  Ys<-Y[R==1]
  Yt<-Y[R==0]
  Xs<-X_std[R==1,]
  Xt<-X_std[R==0,]
  Ss<-S_std[R==1,]
  St<-S_std[R==0,]
  # then create labelled and unlabelled data in the two populations 
  Ysl<-Ys[1:ns]
  Ytl<-Yt[1:nt]
  Xsl<-Xs[1:ns,]
  Xtl<-Xt[1:nt,]
  Ssl<-Ss[1:ns,]
  Stl<-St[1:nt,]
  Ysu<-Ys[(ns+1):(Ns+ns)]
  Ytu<-Yt[(nt+1):(Nt+nt)]
  Xsu<-Xs[(ns+1):(Ns+ns),]
  Xtu<-Xt[(nt+1):(Nt+nt),]
  Ssu<-Ss[(ns+1):(Ns+ns),]
  Stu<-St[(nt+1):(Nt+nt),]
  return(list(Ysl=Ysl, Xsl=Xsl, Ssl=Ssl, Ytl=Ytl, Xtl=Xtl, Stl=Stl,
              Ysu=Ysu, Xsu=Xsu, Ssu=Ssu, Ytu=Ytu, Xtu=Xtu, Stu=Stu))
}
