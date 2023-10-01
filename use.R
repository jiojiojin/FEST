
logit<-function(u){
  return(1/(1+exp(-u)))
}
d_logit<-function(u){
  return(logit(u)*(1-logit(u)))
}

# evaluation funtion
Error_auc<-function(Y,Z,theta){
  theta<-matrix(theta,ncol=1)
  y_hat<-logit(Z%*%theta)
  roc0<-roc(Y,y_hat)
  return(auc(roc0))
}


model_w<-function(Zu,EZ){
  mm<-function(xi){
    xi<-matrix(xi,ncol=1)
    ff<-t(Zu)%*%exp(Zu%*%xi)/nrow(Zu)-EZ
    return(ff)
  }
  return(mm)
}


Loss_beta<-function(X,Y,beta_0,lambda){
  beta_0<-matrix(beta_0,ncol=1)
  Q1<--matrix(rep(1,nrow(X)),nrow=1)%*%(X*c(Y-logit(X%*%beta_0))) 
  Q2<-t(X)%*%(X*c(d_logit(X%*%beta_0)))
  Func<-function(beta){
    beta<-matrix(beta,ncol=1)
    f<-(Q1%*%beta+t(beta-beta_0)%*%Q2%*%(beta-beta_0)/2)/nrow(X)+lambda*norm(beta)
    return(f)
  }
  return(Func)
}

# For TL_TU
Loss_mt_sas<-function(Ztl,Ytl,site_nt,theta_0,lambda,m){
  theta_0<-matrix(theta_0,ncol=1)
  Q1<-matrix(rep(0,ncol(Ztl)),nrow=1)
  Q2<-matrix(rep(0,ncol(Ztl)^2),ncol=ncol(Ztl))
  for (i in 1:m) {
    Ztl_i<-Ztl[site_nt==i,]
    Ytl_i<-Ytl[site_nt==i]
    Q1<-Q1-matrix(rep(1,nrow(Ztl_i)),nrow=1)%*%(Ztl_i*c(Ytl_i-logit(Ztl_i%*%theta_0))) 
    Q2<-Q2+t(Ztl_i)%*%(Ztl_i*c(d_logit(Ztl_i%*%theta_0)))
  }
  Func<-function(theta){
    theta<-matrix(theta,ncol=1)
    f<-(Q1%*%theta+t(theta-theta_0)%*%Q2%*%(theta-theta_0)/2)/nrow(Ztl)+lambda*norm(theta)
    return(f)
  }
  return(Func)
}

# For CS_st
Loss_ms<-function(Zsl,Ysl,site_ns,xi,theta_s,lambda){
  Q1<-matrix(rep(0,ncol(Zsl)),nrow=1)
  Q2<-matrix(rep(0,ncol(Zsl)^2),ncol=ncol(Zsl))
  for (ii in 1:m){
    Zsl_i<-Zsl[site_ns==ii,]
    Ysl_i<-Ysl[site_ns==ii]
    weight_i<-exp(Zsl_i%*%matrix(xi[,ii],ncol=1))
    Q1<-Q1-t(weight_i)%*%(Zsl_i*c(Ysl_i-logit(Zsl_i%*%theta_s)))
    Q2<-Q2+t(Zsl_i*c(weight_i))%*%(Zsl_i*c(d_logit(Zsl_i%*%theta_s)))
  }
  Func<-function(theta){
    theta<-matrix(theta,ncol=1)
    f<-(Q1%*%theta+t(theta-theta_s)%*%Q2%*%(theta-theta_s)/2)/nrow(Zsl)+lambda*norm(theta)
    return(f)
  }
  return(Func)
}

Loss_mt<-function(Ztl,Ytl,site_nt,theta_s,delta_0,lambda){
  delta_0<-matrix(delta_0,ncol=1)
  theta_s<-matrix(theta_s,ncol=1)
  Q1<-matrix(rep(0,ncol(Ztl)),nrow=1)
  Q2<-matrix(rep(0,ncol(Ztl)^2),ncol=ncol(Ztl))
  for (ii in 1:m) {
    Ztl_i<-Ztl[site_nt==ii,]
    Ytl_i<-Ytl[site_nt==ii]
    Q1<-Q1-matrix(rep(1,nrow(Ztl_i)),nrow=1)%*%(Ztl_i*c(Ytl_i-logit(Ztl_i%*%(theta_s+delta_0)))) 
    Q2<-Q2+t(Ztl_i)%*%(Ztl_i*c(d_logit(Ztl_i%*%(theta_s+delta_0))))
  }
  Func<-function(delta){
    delta<-matrix(delta,ncol=1)
    f<-(Q1%*%delta+t(delta-delta_0)%*%Q2%*%(delta-delta_0)/2)/nrow(Ztl)+lambda*norm(delta)
    return(f)
  }
  return(Func)
}

# For CS_free & CS_st
Loss_pool<-function(Zsl,Ysl,site_ns,xi,Ztl,Ytl,site_nt,theta_t,delta_0,lambda){
  Q1<-matrix(rep(0,ncol(Ztl)),nrow=1)
  Q2<-matrix(rep(0,ncol(Ztl)^2),ncol=ncol(Ztl))
  for (ii in 1:m) {
    Ztl_i<-Ztl[site_nt==ii,]
    Ytl_i<-Ytl[site_nt==ii]
    Zsl_i<-Zsl[site_ns==ii,]
    Ysl_i<-Ysl[site_ns==ii]
    weight_i<-exp(Zsl_i%*%matrix(xi[,ii],ncol=1))
    Q1<-Q1-t(weight_i)%*%(Zsl_i*c(Ysl_i-logit(Zsl_i%*%(theta_t-delta_0))))
    Q2<-Q2+t(Zsl_i*c(weight_i))%*%(Zsl_i*c(d_logit(Zsl_i%*%(theta_t-delta_0))))
    Q1<-Q1-matrix(rep(1,nrow(Ztl_i)),nrow=1)%*%(Ztl_i*c(Ytl_i-logit(Ztl_i%*%theta_t))) 
    Q2<-Q2+t(Ztl_i)%*%(Ztl_i*c(d_logit(Ztl_i%*%theta_t)))
  }
  Func<-function(theta){
    theta<-matrix(theta,ncol=1)
    f<-(Q1%*%theta+t(theta-theta_t)%*%Q2%*%(theta-theta_t)/2)/(nrow(Zsl)+nrow(Ztl))+lambda*norm(theta)
    return(f)
  }
  return(Func)
}


# FEST
Loss_mt_d<-function(Ztl,Ytl,theta_s_m,delta_0,lambda){
  delta_0<-matrix(delta_0,ncol=1)
  theta_s_m<-matrix(theta_s_m,ncol=1)
  Q1<--matrix(rep(1,nrow(Ztl)),nrow=1)%*%(Ztl*c(Ytl-logit(Ztl%*%(theta_s_m+delta_0)))) 
  Q2<-t(Ztl)%*%(Ztl*c(d_logit(Ztl%*%(theta_s_m+delta_0))))
  Func<-function(delta){
    delta<-matrix(delta,ncol=1)
    f<-(Q1%*%delta+t(delta-delta_0)%*%Q2%*%(delta-delta_0)/2)/nrow(Ztl)+lambda*norm(delta)
    return(f)
  }
  return(Func)
}

Loss_pool_d<-function(Zsl,Ysl,Ztl,Ytl,theta_t,xi,site_ns,site_nt,delta,lambda){
  Q1<-matrix(rep(0,ncol(Ztl)),nrow=1)
  Q2<-matrix(rep(0,ncol(Ztl)^2),ncol=ncol(Ztl))
  for (ii in site_nt) {
    Ztl_i<-Ztl[site_nt==ii,]
    Ytl_i<-Ytl[site_nt==ii]
    Zsl_i<-Zsl[site_ns==ii,]
    Ysl_i<-Ysl[site_ns==ii]
    weight_i<-exp(Zsl_i%*%matrix(xi[,ii],ncol=1))
    Q1<-Q1-t(weight_i)%*%(Zsl_i*c(Ysl_i-logit(Zsl_i%*%(theta_t-delta[,ii]))))
    Q2<-Q2+t(Zsl_i*c(weight_i))%*%(Zsl_i*c(d_logit(Zsl_i%*%(theta_t-delta[,ii]))))
    Q1<-Q1-matrix(rep(1,nrow(Ztl_i)),nrow=1)%*%(Ztl_i*c(Ytl_i-logit(Ztl_i%*%theta_t))) 
    Q2<-Q2+t(Ztl_i)%*%(Ztl_i*c(d_logit(Ztl_i%*%theta_t)))
  }
  Func<-function(theta){
    theta<-matrix(theta,ncol=1)
    f<-(Q1%*%theta+t(theta-theta_t)%*%Q2%*%(theta-theta_t)/2)/(nrow(Zsl)+nrow(Ztl))+lambda*norm(theta)
    return(f)
  }
  return(Func)
}

# For CS_st
Loss_ms_m<-function(Zsl_i,Ysl_i,weight_i,theta_s,lambda){
  Q1<--t(weight_i)%*%(Zsl_i*c(Ysl_i-logit(Zsl_i%*%theta_s)))
  Q2<-t(Zsl_i*c(weight_i))%*%(Zsl_i*c(d_logit(Zsl_i%*%theta_s)))
  Func<-function(theta){
    theta<-matrix(theta,ncol=1)
    f<-(Q1%*%theta+t(theta-theta_s)%*%Q2%*%(theta-theta_s)/2)/nrow(Zsl_i)+lambda*norm(theta)
    return(f)
  }
  return(Func)
}

# SUPTrans

Loss_delta<-function(X,Y,beta_s,delta_0,lambda){
  delta_0<-matrix(delta_0,ncol=1)
  Q1<--matrix(rep(1,nrow(X)),nrow=1)%*%(X*c(Y-logit(X%*%(beta_s+delta_0)))) 
  Q2<-t(X)%*%(X*c(d_logit(X%*%(beta_s+delta_0))))
  Func<-function(delta){
    delta<-matrix(delta,ncol=1)
    f<-(Q1%*%delta+t(delta-delta_0)%*%Q2%*%(delta-delta_0)/2)/nrow(X)+lambda*norm(delta)
    return(f)
  }
  return(Func)
}

Loss_beta_st<-function(Xt,Yt,Xs,Ys,delta,beta_0,lambda){
  beta_0<-matrix(beta_0,ncol=1)
  delta<-matrix(delta,ncol=1)
  Q1s<--matrix(rep(1,nrow(Xs)),nrow=1)%*%(Xs*c(Ys-logit(Xs%*%(beta_0-delta))))
  Q2s<-t(Xs)%*%(Xs*c(d_logit(Xs%*%(beta_0-delta))))
  Q1t<--matrix(rep(1,nrow(Xt)),nrow=1)%*%(Xt*c(Yt-logit(Xt%*%(beta_0)))) 
  Q2t<-t(Xt)%*%(Xt*c(d_logit(Xt%*%(beta_0))))
  Q1<-Q1s+Q1t
  Q2<-Q2s+Q2t
  Func<-function(beta){
    beta<-matrix(beta,ncol=1)
    f<-(Q1%*%beta+t(beta-beta_0)%*%Q2%*%(beta-beta_0)/2)/(nrow(Xt)+nrow(Xs))+lambda*norm(beta)
    return(f)
  }
  return(Func)
}
