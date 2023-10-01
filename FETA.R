FETA<-function(t1,t2,t3){
  # initialization
  m0<-which.max(table(site_ns))
  Zsl_m0<-Zsl[site_ns==m0,]
  Ysl_m0<-Ysl[site_ns==m0]
  cv.out<-cv.glmnet(Zsl_m0[,2:(1+p)],Ysl_m0,family='binomial',alpha=1)
  best.lambda<-cv.out$lambda.min
  ini.mod<-glmnet(Zsl_m0[,2:(1+p)],Ysl_m0,family='binomial',alpha=1,lambda=best.lambda)
  beta_ini<-matrix(coef(ini.mod),ncol=1)
  
  # estimate beta_s
  grid<-exp(-seq(0,5,1))
  beta_s<-matrix(c(beta_ini,rep(NA,(p+1)*t1)),ncol=t1+1) # record the final theta
  for (j in 1:t1) {
    AUC<-c()
    for (i in 1:length(grid)) {
      model_om<-nlm(Loss_beta(Zsl[,1:(1+p)],Ysl,beta_s[,j],grid[i]),beta_s[,j],iterlim=500)
      para<-model_om$estimate
      AUC[i]<-Error_auc(Ys_vali,Xs_vali,para) 
    }
    model_om<-nlm(Loss_beta(Zsl[,1:(1+p)],Ysl,beta_s[,j],grid[which.max(AUC)]),beta_s[,j],iterlim=500)
    beta_s[,j+1]<-model_om$estimate
  }
  auc_beta_s<-apply(beta_s,2,function(x)Error_auc(Ys_vali,Xs_vali,x))
  
  # estimate difference between source and target
  delta_0<-rep(0,p+1)
  delta_hat<-matrix(c(delta_0,rep(NA,(p+1)*t2)),ncol=t2+1)
  ws<-beta_s[,which.max(auc_beta_s)]
  for (j in 1:t2) {
    AUC_t<-c()
    for (i in 1:length(grid)) {
      model_t<-nlm(Loss_delta(Ztl[,1:(1+p)],Ytl,ws,delta_hat[,j],grid[i]),delta_hat[,j],iterlim=500)
      para<-model_t$estimate
      AUC_t[i]<-Error_auc(Yt_vali,Xt_vali,para)
    }
    model_t<-nlm(Loss_delta(Ztl[,1:(1+p)],Ytl,ws,delta_hat[,j],grid[which.max(AUC_t)]),delta_hat[,j],iterlim=500)
    delta_hat[,j+1]<-matrix(model_t$estimate,ncol=1)
  }
  auc_beta_t<-apply(delta_hat,2,function(x)Error_auc(Yt_vali,Xt_vali,x+ws))
  # estimate beta
  # initialization 
  delta_T<-delta_hat[,which.max(auc_beta_t)]
  beta<-matrix(c(delta_T+ws,rep(NA,t3*(p+1))),ncol=t3+1)
  for (j in 1:t3){
    AUC_b<-c()
    for (i in 1:length(grid)) {
      model_mb<-nlm(Loss_beta_st(Ztl[,1:(1+p)],Ytl,Zsl[,1:(1+p)],Ysl,delta_T,beta[,j],grid[i]),beta[,j],iterlim = 500)
      para<-model_mb$estimate
      AUC_b[i]<-Error_auc(Yt_vali,Xt_vali,para) 
    }
    model_mb<-nlm(Loss_beta_st(Ztl[,1:(1+p)],Ytl,Zsl[,1:(1+p)],Ysl,delta_T,beta[,j],grid[which.max(AUC_b)]),beta[,j],iterlim = 500)
    beta[,j+1]<-model_mb$estimate
  }
  auc_beta<-apply(beta,2,function(x)Error_auc(Yt_test,Xt_test,x))
  abs_beta<-apply(beta,2,function(x)norm(matrix(x-beta_true,ncol=1)))
  l2_beta<-apply(beta,2,function(x)norm(matrix(x-beta_true,ncol=1),type = '2'))
  return(list(auc_beta=auc_beta,abs_beta=abs_beta,l2_beta=l2_beta))
}
