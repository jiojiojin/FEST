
TSUP<-function(t3){
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
  grid<-exp(-seq(-0.5,3,0.5))
  # initialization 
  beta<-matrix(c(beta_ini,rep(NA,t3*(p+1))),ncol=t3+1)
  for (j in 1:t3){
    AUC_b<-c()
    for (i in 1:length(grid)) {
      model_mb<-nlm(Loss_beta(Ztl[,1:(1+p)],Ytl,beta[,j],grid[i]),beta[,j],iterlim = 500)
      para<-model_mb$estimate
      AUC_b[i]<-Error_auc(Yt_vali,Xt_vali,para) 
    }
    model_mb<-nlm(Loss_beta(Ztl[,1:(1+p)],Ytl,beta[,j],grid[which.max(AUC_b)]),beta[,j],iterlim = 500)
    beta[,j+1]<-model_mb$estimate
  }
  auc_beta<-apply(beta,2,function(x)Error_auc(Yt_test,Xt_test,x))
  abs_beta<-apply(beta,2,function(x)norm(matrix(x-beta_true,ncol=1)))
  l2_beta<-apply(beta,2,function(x)norm(matrix(x-beta_true,ncol=1),type = '2'))
  return(list(auc_beta=auc_beta,abs_beta=abs_beta,l2_beta=l2_beta))
}
