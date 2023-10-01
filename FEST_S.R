FEST_S<-function(t0,t1,t2,t3){
  theta_s_hat<-theta_s<-matrix(rep(0,(d+1)*m),ncol=m) # record the final theta
  for (ms in 1:m) {
    Zsl_i<-Zsl[site_ns==ms,]
    Ysl_i<-Ysl[site_ns==ms]
    weight_i<-c(exp(Zsl_i%*%matrix(xi[,ms],ncol=1)))
    cv.out<-cv.glmnet(Zsl_i[,2:(d+1)],Ysl_i,family='binomial',alpha=1,seed=1,weights = weight_i)
    best.lambda<-cv.out$lambda.min
    ini.mod<-glmnet(Zsl_i[,2:(d+1)],Ysl_i,family='binomial',alpha=1,
                    lambda=best.lambda,seed=1,weights = weight_i)
    theta_s[,ms]<-matrix(coef(ini.mod),ncol=1)
  }
  auc_theta_s<-apply(theta_s,2,function(x)Error_auc(Yt_vali,Zt_vali,x))
  #theta_s_hat<-matrix(c(theta_s,rep(NA,(d+1)*t1)),ncol=t1+1) # record the final theta
  grid<-exp(-seq(-0.5,3,0.5))
  for (ms in 1:m) {
    Zsl_i<-Zsl[site_ns==ms,]
    Ysl_i<-Ysl[site_ns==ms]
    weight_i<-c(exp(Zsl_i%*%matrix(xi[,ms],ncol=1)))
    theta_s_m<-matrix(c(theta_s[,ms],rep(0,(d+1)*t0)),ncol=t0+1)
    for (j in 1:t0) {
      AUC<-c()
      for (i in 1:length(grid)) {
        model_ms<-nlm(Loss_ms_m(Zsl_i,Ysl_i,weight_i,theta_s_m[,j],grid[i]),theta_s_m[,j],iterlim = 500)
        para<-model_ms$estimate
        AUC[i]<-Error_auc(Yt_vali,Zt_vali,para) 
      }
      model_ms<-nlm(Loss_ms_m(Zsl_i,Ysl_i,weight_i,theta_s_m[,j],grid[which.max(AUC)]),theta_s_m[,j],iterlim = 500)
      theta_s_m[,j+1]<-matrix(model_ms$estimate,ncol=1)
    }
    auc_s_m<-apply(theta_s_m,2,function(x)Error_auc(Yt_vali,Zt_vali,x))
    theta_s_hat[,ms]<-theta_s_m[,which.max(auc_s_m)]
  }
  ## initialization of \theta
  theta_tT<-delta_T<-matrix(rep(NA,m*(1+d)),ncol=m)
  for (k in 1:m) {
    delta_0<-rep(0,d+1)
    theta_s_k<-theta_s_hat[,k]
    theta_t_hat<-matrix(c(theta_s_k,rep(NA,(d+1)*t1)),ncol=t1+1) # record the final theta
    for (j in 1:t1) {
      AUC<-c()
      for (i in 1:length(grid)) {
        model_mt<-nlm(Loss_mt_d(Ztl,Ytl,theta_s_k,delta_0,grid[i]),delta_0,iterlim = 500)
        para<-model_mt$estimate
        AUC[i]<-Error_auc(Yt_vali,Zt_vali,para+theta_s_k) 
      }
      model_mt<-nlm(Loss_mt_d(Ztl,Ytl,theta_s_k,delta_0,grid[which.max(AUC)]),delta_0,iterlim = 500)
      delta_0<-matrix(model_mt$estimate,ncol=1)
      theta_t_hat[,j+1]<-theta_t_hat[,1]+delta_0
    }
    auc_theta_t<-apply(theta_t_hat,2,function(x)Error_auc(Yt_vali,Zt_vali,x))
    theta_tT[,k]<-theta_t_hat[,which.max(auc_theta_t)] #theta_T
    delta_T[,k]<-theta_tT[,k]-theta_s_k
  }
  auc_tT<-apply(theta_tT,2,function(x)Error_auc(Yt_vali,Zt_vali,x))
  theta_ini<-theta_tT[,(which.max(auc_tT))]
  
  theta_t<-matrix(rep(0,(d+1)*m),ncol=m) # record the final theta
  for (ii in 1:m) {
    Ztl_i<-Ztl[site_nt==ii,]
    Ytl_i<-Ytl[site_nt==ii]
    cv.out<-cv.glmnet(Ztl_i[,2:(d+1)],Ytl_i,family='binomial',alpha=1,seed=1)
    best.lambda<-cv.out$lambda.min
    ini.mod<-glmnet(Ztl_i[,2:(d+1)],Ytl_i,family='binomial',alpha=1,
                    lambda=best.lambda,seed=1)
    theta_t[,ii]<-matrix(coef(ini.mod),ncol=1)
  }
  auc_t<-apply(theta_t,2,function(x)Error_auc(Yt_vali,Zt_vali,x))
  
  aaa<-seq(1,m)
  site_use<-aaa[auc_tT>=min(auc_t)]
  site_ns_new<-site_ns[site_ns %in% site_use]
  Zsl_new<-Zsl[(site_ns %in% site_use),]
  Ysl_new<-Ysl[site_ns %in% site_use]
  
  theta_hat<-matrix(c(theta_ini,rep(NA,(d+1)*t2)),ncol=t2+1) # record the final theta
  grid<-c(seq(5,1,-0.5),exp(-seq(0,5)))
  for (j in 1:t2) {
    AUC<-c()
    for (i in 1:length(grid)) {
      model_mt<-nlm(Loss_pool_d(Zsl_new,Ysl_new,Ztl,Ytl,theta_hat[,j],xi,site_ns_new,site_nt,
                                delta_T,grid[i]),theta_hat[,j],iterlim = 500)
      para<-model_mt$estimate
      AUC[i]<-Error_auc(Yt_vali,Zt_vali,para) 
    }
    model_mt<-nlm(Loss_pool_d(Zsl_new,Ysl_new,Ztl,Ytl,theta_hat[,j],xi,site_ns_new,site_nt,
                              delta_T,grid[which.max(AUC)]),theta_hat[,j],iterlim = 500)
    theta_hat[,j+1]<-matrix(model_mt$estimate,ncol=1)
  }
  auc_pool<-apply(theta_hat,2,function(x)Error_auc(Yt_vali,Zt_vali,x))
  # estimate beta
  theta_final<-theta_hat[,which.max(auc_pool)] #theta_T
  Y_t_bar<-logit(Ztu%*%theta_final)
  #theta_tT<-theta_t_hat[,which.max(auc_theta_t)] #theta_T
  # estimate beta
  Y_t_bar<-c(data$Ytl,Y_t_bar)
  X_t<-cbind(rep(1,nt+Nt),rbind(data$Xtl,data$Xtu))
  
  beta_m<-matrix(rep(NA,(1+p)*m),ncol=m)
  for (m0 in 1:m) {
    Xtl_m0<-Ztl[site_nt==m0,1:(1+p)]
    Ytl_m0<-Ytl[site_nt==m0]
    cv.out<-cv.glmnet(Xtl_m0[,2:(1+p)],Ytl_m0,family='binomial',alpha=1,set.seed(1))
    best.lambda<-cv.out$lambda.min
    ini.mod<-glmnet(Xtl_m0[,2:(1+p)],Ytl_m0,family='binomial',alpha=1,lambda=best.lambda,set.seed(1))
    beta_m[,m0]<-matrix(coef(ini.mod),ncol=1)
  }
  auc_beta_m<-apply(beta_m,2,function(x)Error_auc(Yt_vali,Xt_vali,x))
  beta_ini<-beta_m[,(which.max(auc_beta_m))]
  
  beta<-matrix(c(beta_ini,rep(NA,t3*(p+1))),ncol=t3+1)
  for (j in 1:t3){
    AUC<-c()
    for (i in 1:length(grid)) {
      model_mb<-nlm(Loss_beta(X_t,Y_t_bar,beta[,j],grid[i]),beta[,j],iterlim = 500)
      para<-model_mb$estimate
      AUC[i]<-Error_auc(Yt_vali,Xt_vali,para) 
    }
    model_mb<-nlm(Loss_beta(X_t,Y_t_bar,beta[,j],grid[which.max(AUC)]),beta[,j],iterlim = 500)
    beta[,j+1]<-model_mb$estimate
  }
  auc_beta<-apply(beta,2,function(x)Error_auc(Yt_test,Xt_test,x))
  abs_beta<-apply(beta,2,function(x)norm(matrix(x-beta_true,ncol=1)))
  l2_beta<-apply(beta,2,function(x)norm(matrix(x-beta_true,ncol=1),type = '2'))
  
  return(list(auc_beta=auc_beta,abs_beta=abs_beta,l2_beta=l2_beta))
}

