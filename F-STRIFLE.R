F_STRIFLE<-function(t1,t2,t3){
  theta_s<-matrix(rep(0,(d+1)*m),ncol=m) # record the final theta
  for (ms in 1:m) {
    Zsl_i<-Zsl[site_ns==ms,]
    Ysl_i<-Ysl[site_ns==ms]
    cv.out<-cv.glmnet(Zsl_i[,2:(d+1)],Ysl_i,family='binomial',alpha=1,seed=1)
    best.lambda<-cv.out$lambda.min
    ini.mod<-glmnet(Zsl_i[,2:(d+1)],Ysl_i,family='binomial',alpha=1,
                    lambda=best.lambda,seed=1)
    theta_s[,ms]<-matrix(coef(ini.mod),ncol=1)
  }
  auc_s<-apply(theta_s,2,function(x)Error_auc(Yt_vali,Zt_vali,x))
  theta_ini<-theta_s[,which.max(auc_s)]
  
  # Step 2
  lam<-3*sqrt(log(p+q)/ns)
  grid<-lam*1.15^c(0:5)
  theta_s_hat<-matrix(c(theta_ini,rep(NA,(d+1)*t1)),ncol=t1+1) # record the final theta
  for (j in 1:t1) {
    AUC<-c()
    for (i in 1:length(grid)) {
      model_ms<-nlm(Loss_ms(Zsl,Ysl,site_ns,xi,theta_s_hat[,j],grid[i]),theta_s_hat[,j],iterlim=500)
      para<-model_ms$estimate
      AUC[i]<-Error_auc(Yt_vali,Zt_vali,para) # CV on the source labeled data of site_m0
    }
    model_ms<-nlm(Loss_ms(Zsl,Ysl,site_ns,xi,theta_s_hat[,j],grid[which.max(AUC)]),theta_s_hat[,j],iterlim=500)
    theta_s_hat[,j+1]<-model_ms$estimate
  }
  auc_theta_s<-apply(theta_s_hat,2,function(x)Error_auc(Yt_vali,Zt_vali,x))
  # Step 3
  lam<-3*sqrt(log(p+q)/nt)
  grid<-lam*1.15^c(0:5)
  delta_0<-rep(0,d+1)
  #theta_t_hat<-matrix(c(theta_s_hat[,which.max(auc_theta_s[-1])+1],rep(NA,(d+1)*t2)),ncol=t2+1) # record the final theta
  theta_t_hat<-matrix(c(theta_s_hat[,which.max(auc_theta_s)],rep(NA,(d+1)*t2)),ncol=t2+1) # record the final theta
  
  for (j in 1:t2) {
    AUC_t<-c()
    for (i in 1:length(grid)) {
      model_mt<-nlm(Loss_mt(Ztl,Ytl,site_nt,theta_t_hat[,1],delta_0,grid[i]),delta_0,iterlim = 500)
      para<-model_mt$estimate
      AUC_t[i]<-Error_auc(Yt_vali,Zt_vali,para) 
    }
    model_mt<-nlm(Loss_mt(Ztl,Ytl,site_nt,theta_t_hat[,1],delta_0,grid[which.max(AUC_t)]),delta_0,iterlim = 500)
    delta_0<-matrix(model_mt$estimate,ncol=1)
    theta_t_hat[,j+1]<-theta_t_hat[,1]+delta_0
  }
  auc_theta_t<-apply(theta_t_hat,2,function(x)Error_auc(Yt_vali,Zt_vali,x))
  theta_tT<-theta_t_hat[,which.max(auc_theta_t)] #theta_T
  # estimate beta
  Y_t_bar<-c(Ytl,logit(Ztu%*%theta_tT))
  X_t<-cbind(rep(1,nt+Nt),rbind(Xtl,Xtu))
  
  ### initialization of beta
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
  
  lam<-3*sqrt(log(p)/(nt+Nt))
  grid<-lam*1.15^c(0:5)
  beta<-matrix(c(beta_ini,rep(NA,t3*(p+1))),ncol=t3+1)
  for (j in 1:t3){
    AUC_b<-c()
    for (i in 1:length(grid)) {
      model_mb<-nlm(Loss_beta(X_t,Y_t_bar,beta[,j],grid[i]),beta[,j],iterlim = 500)
      para<-model_mb$estimate
      AUC_b[i]<-Error_auc(Yt_vali,Xt_vali,para) 
    }
    model_mb<-nlm(Loss_beta(X_t,Y_t_bar,beta[,j],grid[which.max(AUC_b)]),beta[,j],iterlim = 500)
    beta[,j+1]<-model_mb$estimate
  }
  auc_beta<-apply(beta,2,function(x)Error_auc(Yt_test,Xt_test,x))
  abs_beta<-apply(beta,2,function(x)norm(matrix(x-beta_true,ncol=1)))
  l2_beta<-apply(beta,2,function(x)norm(matrix(x-beta_true,ncol=1),type = '2'))
  
  return(list(auc_beta=auc_beta,abs_beta=abs_beta,l2_beta=l2_beta))
}
