library(doParallel)
library(tidyverse)
library(foreach)
library(grf)
library(latex2exp)

myCluster <- makeCluster(detectCores()-1)
registerDoParallel(myCluster)

dtrue = 2
nonlin = function(z){1/(1+exp(-z))}
baselinex = function(x){(x[,1]<0)*(-1+2*(x[,2]>0))}
effectx   = function(x){2*(x[,1]>=0)*(-1+2*(x[,2]>0))}
baseline0 = 0
effect0   = 0
mu = function(x,a,beta){nonlin(beta*(baselinex(x)+(a-0.5)*effectx(x))+baseline0+effect0)}

### true bounds

ntest = 1000000
Xtest = matrix(data=rnorm(ntest*dtrue),nrow=ntest,ncol=dtrue)
betas = seq(0,8,.1)
truebounds = foreach(beta=betas, .combine=rbind) %dopar% {
  mu1test = mu(Xtest,0,beta)
  mu0test = mu(Xtest,1,beta)
  ub0=min(mean(mu0test),1-mean(mu1test))
  lb0=-min(mean(mu1test)-mean(mu0test),0)
  ubx=mean(mu0test+(1-mu1test<mu0test)*(1-mu1test-mu0test))
  lbx=-mean((mu1test-mu0test<0)*(mu1test-mu0test))
  data.frame(lb0=lb0,ub0=ub0,lbx=lbx,ubx=ubx,beta=beta)
}
truebounds = truebounds %>% mutate(opt0 = ub0-lb0, optx = ubx-lbx)
truebounds %>% pivot_longer(!beta) %>% ggplot() + aes(x=beta, y=value, color=name) + geom_line()


### actual sim

# splits data into K folds
make.cvgroup = function(n, K, right = TRUE) {
  split     = runif(n)
  return(as.numeric(cut(split, quantile(split, probs = seq(0, 1, 1/K)), include.lowest = TRUE, right = right)))
}

# conditions on fraction of 1s to 0s being similar across folds
make.cvgroup.balanced = function(data, K, form_t) {
  cvgroup = numeric(nrow(data))
  cvgroup[data[[form_t]]==1] = make.cvgroup(sum(data[[form_t]]==1), K, right = TRUE)
  cvgroup[data[[form_t]]==0] = make.cvgroup(sum(data[[form_t]]==0), K, right = FALSE)
  return(cvgroup)
}

propensity = function(x){1/(1+exp((1/2-(x[,3]>0))-(1/2-(x[,4]>0))/2))}
dextra = 5
d = dtrue+dextra

dorun = function(beta,n,seed) {
  set.seed(seed)
  
  Xorg = matrix(data=rnorm(n*d),nrow=n,ncol=d)
  X = Xorg
  ee  = propensity(X)
  A = runif(n)<ee
  muax = A*mu(X,1,beta) + (1-A)*mu(X,0,beta)
  Y = runif(n)<muax
  
  X[,1:2] = X[,1:2]*(2*apply(X[,(dtrue+1):d]>0,1,function(x){reduce(x,xor)})-1)
  
  data = data.frame(X=X, A=A, Y=Y)
  
  K = 5
  cvgroup = data %>% make.cvgroup.balanced(., K, 'A')
  
  mu0.pred  = numeric(nrow(data))
  mu1.pred  = numeric(nrow(data))
  e.pred   = numeric(nrow(data))
  for (k in 1:K) {
    fmu0 = regression_forest(X[cvgroup!=k & data$A==0,], Y[cvgroup!=k & data$A==0])
    fmu1 = regression_forest(X[cvgroup!=k & data$A==1,], Y[cvgroup!=k & data$A==1])
    fehat = regression_forest(X[cvgroup!=k,], A[cvgroup!=k])
    mu0.pred[cvgroup==k] = predict(fmu0, X[cvgroup==k,])$predictions
    mu1.pred[cvgroup==k] = predict(fmu1, X[cvgroup==k,])$predictions
    e.pred[cvgroup==k] = predict(fehat, X[cvgroup==k,])$predictions
  }
  data = data %>% mutate(
    mu0 = mu0.pred,
    mu1 = mu1.pred,
    e = e.pred
  )
  tau.pred  = numeric(nrow(data))
  sau.pred  = numeric(nrow(data))
  Yalt = A*Y - (1-A)*Y
  for (k in 1:K) {
    ftau = causal_forest(X[cvgroup!=k,], Y[cvgroup!=k], A[cvgroup!=k])
    fsau = causal_forest(X[cvgroup!=k,], Yalt[cvgroup!=k], A[cvgroup!=k])
    tau.pred[cvgroup==k] = predict(ftau, X[cvgroup==k,])$predictions
    sau.pred[cvgroup==k] = predict(fsau, X[cvgroup==k,])$predictions
  }
  data = data %>% mutate(
    tau = tau.pred,
    sau = sau.pred
  )
  
  return(rbind(
    data %>% mutate(eif = -(tau<=0)*(mu1-mu0+A*(Y-mu1)/e-(1-A)*(Y-mu0)/(1-e))) %>% summarise(m=mean(eif), s=sd(eif)/sqrt(n())) %>% mutate(type = 'lb.dr', beta=beta, n=n, seed=seed)
    ,
    data %>% mutate(eif = -(mu1-mu0<=0)*(mu1-mu0)) %>% summarise(m=mean(eif), s=sd(eif)/sqrt(n())) %>% mutate(type = 'lb.pi', beta=beta, n=n, seed=seed)
    ,
    data %>% mutate(eif = (1-sau<=0)*(1-mu1-A*(Y-mu1)/e)+(1-sau>0)*(mu0+(1-A)*(Y-mu0)/(1-e))) %>% summarise(m=mean(eif), s=sd(eif)/sqrt(n())) %>% mutate(type = 'ub.dr', beta=beta, n=n, seed=seed)
    ,
    data %>% mutate(eif = (1-mu1-mu0<=0)*(1-mu1)+(1-mu1-mu0>0)*(mu0)) %>% summarise(m=mean(eif), s=sd(eif)/sqrt(n())) %>% mutate(type = 'ub.pi', beta=beta, n=n, seed=seed)
    ,
    data %>% mutate(eif = (tau<=0)*(mu1-mu0+A*(Y-mu1)/e-(1-A)*(Y-mu0)/(1-e))+(1-sau<=0)*(1-mu1-A*(Y-mu1)/e)+(1-sau>0)*(mu0+(1-A)*(Y-mu0)/(1-e))) %>% summarise(m=mean(eif), s=sd(eif)/sqrt(n())) %>% mutate(type = 'opt.dr', beta=beta, n=n, seed=seed)
    ,
    data %>% mutate(eif = (mu1-mu0<=0)*(mu1-mu0)+(1-mu1-mu0<=0)*(1-mu1)+(1-mu1-mu0>0)*(mu0)) %>% summarise(m=mean(eif), s=sd(eif)/sqrt(n())) %>% mutate(type = 'opt.pi', beta=beta, n=n, seed=seed)
  ))
}

betas = c(3)
ns = c(200,400,800,1600,3200,6400,12800)
seeds.per.run = 4
runs = 25

simresultslist = list()

for (j in 1:runs) {
  print(paste('run',j,Sys.time()))
  simresultslist[[j]] = foreach(seed = (1+(j-1)*seeds.per.run):(j*seeds.per.run), .combine=rbind, .inorder=FALSE) %:%
    foreach(beta = betas, .combine=rbind, .inorder=FALSE) %:%
      foreach(n = ns, .combine=rbind, .inorder=FALSE, .packages=c('foreach','grf','tidyverse')) %dopar% {
        dorun(beta,n,seed)
      }
}

sims = do.call(rbind,simresultslist)

simsjoin = sims %>% separate(type, c('bound', 'type')) %>% left_join(
  truebounds %>% select(!c(lb0,ub0,opt0)) %>% pivot_longer(!beta) %>% mutate(bound=substr(name,0,nchar(name)-1)) %>% select(!name) %>% rename(true = value)
  ,
  by = c('bound', 'beta')
)

estimand.names = unname(TeX(c("$\\mathrm{FNA}^-$", "$\\mathrm{FNA}^+$", "$\\mathrm{FNA}^+_{1-\\pi^*\\rightarrow\\pi^*}$")))
simsjoin %>% group_by(bound,type,beta,n) %>% summarise(rmse = sqrt(mean((m-true)^2)), rmsese = sd((m-true)^2)/2./rmse/sqrt(n()), cov80 = mean(abs(m-true)<=1.282*s), cov90 = mean(abs(m-true)<=1.645*s), cov95 = mean(abs(m-true)<=1.960*s), sermse = sqrt(mean((s-rmse)^2))) %>% 
  filter(beta==3) %>% filter(n>300) %>%
  mutate(Estimator=recode_factor(type,dr='Alg. 1',pi='Plugin')) %>%
  mutate(Estimand=recode_factor(bound,lb='FNA^+',ub='FNA^-',opt='Opt')) %>%
  ggplot() + aes(x=n,y=rmse,ymax=rmse+1.96*rmsese,ymin=rmse-1.96*rmsese,color=Estimand,fill=Estimand,linetype=Estimator,shape=Estimator) + geom_line() + geom_ribbon(alpha=0.5,color=NA) + geom_point() +
  scale_color_discrete(labels = estimand.names) +
  scale_fill_discrete(labels = estimand.names)  + theme_minimal() + 
  theme(legend.box = "horizontal", legend.direction = "vertical", legend.position = c(0.7, 0.8), axis.title.y=element_blank(), axis.title.x = element_text(hjust=.9975,vjust=6), legend.text.align = 0) + xlab(TeX('$n$'))

ggsave('simsrmse.pdf',plot=last_plot()    + theme(plot.margin=grid::unit(c(0,0,0,0), "mm")), dpi = 300, height = 4, width = 4)

simsjoin %>% group_by(bound,type,beta,n) %>% summarise(cov95 = mean(abs(m-true)<=1.960*s), cov95se = sd(abs(m-true)<=1.960*s)/sqrt(n())) %>% 
  filter(beta==3) %>% filter(n>300) %>%
  mutate(Estimator=recode_factor(type,dr='Alg. 1',pi='Plugin')) %>%
  mutate(Estimand=recode_factor(bound,lb='FNA^+',ub='FNA^-',opt='Opt')) %>%
  ggplot() + aes(x=n,y=cov95,ymax=cov95+1.96*cov95se,ymin=cov95-1.96*cov95se,color=Estimand,fill=Estimand,linetype=Estimator,shape=Estimator) + geom_line() + geom_ribbon(alpha=0.5,color=NA) + geom_point() +
  scale_color_discrete(labels = estimand.names) + 
  scale_fill_discrete(labels = estimand.names)  + theme_minimal() + 
  theme(legend.box = "horizontal", legend.direction = "vertical", legend.position = c(0.7, 0.3), axis.title.y=element_blank(), axis.title.x = element_text(hjust=.9975,vjust=6), legend.text.align = 0) + xlab(TeX('$n$'))

ggsave('simscov.pdf',plot=last_plot()    + theme(plot.margin=grid::unit(c(0,0,0,0), "mm")), dpi = 300, height = 4, width = 4)

betas = seq(0,8,.1)
runsplits = 3
runsplit = make.cvgroup(length(betas),runsplits)
cisimslist = list()
foreach(j = 1:runsplits) %do% {
  print(paste('run',j,Sys.time()))
  cisimslist[[j]] = foreach(beta = betas[runsplit==j], .combine=rbind, .inorder=FALSE, .packages=c('foreach','grf','tidyverse')) %dopar% {
    dorun(beta,12800,0)
  }
}

cisims = do.call(rbind,cisimslist)

estimand.names = unname(TeX(c("$\\mathrm{FNA}^-$", "$\\mathrm{FNA}^+$", "$\\mathrm{FNA}^+_{1-\\pi^*\\rightarrow\\pi^*}$")))
rbind(
cisims %>% separate(type, c('bound', 'type')) %>% filter(type=='dr') %>% select(m,s,bound,beta,type),
truebounds %>% select(!c(lb0,ub0,opt0)) %>% pivot_longer(!beta) %>% mutate(bound=substr(name,0,nchar(name)-1),s=NA,type='true') %>% select(!name) %>% rename(m = value)
) %>% 
  mutate(Bound=recode_factor(type,dr='Alg. 1',true='True')) %>%
  mutate(Estimand=recode_factor(bound,lb='FNA^+',ub='FNA^-',opt='Opt')) %>%
  ggplot() + aes(x=beta, y=m, ymax=m+1.96*s, ymin=m-1.96*s, color=Estimand, fill=Estimand, linetype=Bound) + geom_line() + geom_ribbon(alpha=.5, color=NA) +
  scale_color_discrete(labels = estimand.names) +
  scale_fill_discrete(labels = estimand.names)  + theme_minimal() + scale_x_continuous(breaks=c(0,2,4,6)) +
  theme(legend.box = "horizontal", legend.direction = "vertical", legend.position = c(0.7, 0.8), axis.title.y=element_blank(), axis.title.x = element_text(hjust=.9975,vjust=6), legend.text.align = 0) + xlab(TeX('$\\beta$'))

ggsave('cisims.pdf',plot=last_plot()    + theme(plot.margin=grid::unit(c(0,0,0,0), "mm")), dpi = 300, height = 4, width = 4)





