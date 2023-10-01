TSEM<-function(t2,t3){
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
  auc_t_m<-apply(theta_t,2,function(x)Error_auc(Yt_vali,Zt_vali,x))

  theta_ini<-theta_t[,which.max(auc_t_m)]
  
  grid<-exp(-seq(0,3,1))
  theta_t_hat<-matrix(c(theta_ini,rep(NA,(d+1)*t2)),ncol=t2+1) # record the final theta
  for (j in 1:t2) {
    AUC_t<-c()
    for (i in 1:length(grid)) {
      model_mt<-nlm(Loss_mt_sas(Ztl,Ytl,site_nt,theta_t_hat[,j],grid[i],m),
                    theta_t_hat[,j],iterlim = 500)
      para<-model_mt$estimate
      AUC_t[i]<-Error_auc(Yt_vali,Zt_vali,para) 
    }
    model_mt<-nlm(Loss_mt_sas(Ztl,Ytl,site_nt,theta_t_hat[,j],grid[which.max(AUC_t)],m),theta_t_hat[,j],iterlim = 500)
    theta_t_hat[,j+1]<-matrix(model_mt$estimate,ncol=1)
  }
  auc_theta_t<-apply(theta_t_hat,2,function(x)Error_auc(Yt_vali,Zt_vali,x))
  # estimate beta
  theta_tT<-theta_t_hat[,which.max(auc_theta_t)]
  
  Y_t_bar<-c(data$Ytl,logit(Ztu%*%theta_tT))
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
  
  # initialization 
  beta<-matrix(c(beta_ini,rep(NA,t3*(p+1))),ncol=t3+1)
  max_b<-c()
  for (j in 1:t3){
    AUC_b<-c()
    for (i in 1:length(grid)) {
      model_mb<-nlm(Loss_beta(X_t,Y_t_bar,beta[,j],grid[i]),beta[,j],iterlim = 500)
      para<-model_mb$estimate
      AUC_b[i]<-Error_auc(Yt_vali,Xt_vali,para) 
    }
    max_b[j]<-which.max(AUC_b)
    model_mb<-nlm(Loss_beta(X_t,Y_t_bar,beta[,j],grid[which.max(AUC_b)]),beta[,j],iterlim = 500)
    beta[,j+1]<-model_mb$estimate
  }
  auc_beta<-apply(beta,2,function(x)Error_auc(Yt_test,Xt_test,x))
  abs_beta<-apply(beta,2,function(x)norm(matrix(x-beta_true,ncol=1)))
  l2_beta<-apply(beta,2,function(x)norm(matrix(x-beta_true,ncol=1),type = '2'))
  return(list(auc_beta=auc_beta,abs_beta=abs_beta,l2_beta=l2_beta))
}
