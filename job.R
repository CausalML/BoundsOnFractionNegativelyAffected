library(tidyverse)
library(foreach)
library(latex2exp)
library(gbm)
library(cmna)

job_org = read.csv('behanghel_processed.csv')
job = job_org %>% filter(POIDS_PZ_6MOIS>0) %>% mutate(A_standard = CLA, A_private = OPP, A_public = CVE, Y = EMPLOI_6MOIS) %>% mutate(sw = POIDS_PZ_6MOIS, ipw = 1. / (A_standard*mean(sw*A_standard) + A_private*mean(sw*A_private) + A_public*mean(sw*A_public)))

# Prep covariates

Xall = c(
  'College_education',
  'nivetude2',
  'Vocational',
  'High_school_dropout',
  'Manager',
  'Technician',
  'Skilled_clerical_worker',
  'Unskilled_clerical_worker',
  'Skilled_blue_colar',
  'Unskilled_blue_colar',
  'Woman',
  'Married',
  'French',
  'African',
  'Other_Nationality',
  'Paris_region',
  'North',
  'Other_regions',
  'Employment_component_level_1',
  'Employment_component_level_2',
  'Employment_component_missing',
  'Economic_Layoff',
  'Personnal_Layoff',
  'End_of_Fixed_Term_Contract',
  'End_of_Temporary_Work',
  'Other_reasons_of_unemployment',
  'Statistical_risk_level_2',
  'Statistical_risk_level_3',
  'Other_Statistical_risk',
  'Search_for_a_full_time_position',
  'Sensitive_suburban_area',
  'Insertion',
  'Interim',
  'Conseil',
  'No_child',
  'One_child',
  'More_than_one_child',
  'Wage_target_1200_1349_euros',
  'Wage_target_1350_1549_euros',
  'Wage_target_1550_1799_euros',
  'Wage_target_1800_2200_euros',
  'Wage_target_2200_euros',
  'age',
  'exper',
  'mois_saisie_occ',
  'ndem',
  'ndem1'
)

job = job%>%mutate(exper = factor(ntile(exper, 10)), age = factor(ntile(age, 10)), ndem1 = (ndem>1), mois_saisie_occ = factor(ntile(mois_saisie_occ,10)))


AA = 'A_private'
BB = 'A_public'

arms = c('A_standard','A_private','A_public')

results = foreach(AA = arms, .combine = rbind) %:% foreach(BB = arms, .combine = rbind) %do% {
if(AA==BB){data.frame()}else{
job_binary = job %>% filter((!!as.symbol(AA)) == 1 | (!!as.symbol(BB)) == 1) %>% mutate(sw = sw/mean(sw)) %>% mutate(A = (!!as.symbol(BB)), ipw = 1 / (A_standard*mean(sw*A_standard) + A_private*mean(sw*A_private) + A_public*mean(sw*A_public)))

mu1.pred = predict(glm(formula=as.formula(paste('Y ~ (',paste(Xall,collapse = ' + '),')')), data=job_binary[job_binary$A==1,], family = "binomial", weights=(job_binary$sw*job_binary$ipw)[job_binary$A==1]),job_binary,type="response")
mu0.pred = predict(glm(formula=as.formula(paste('Y ~ (',paste(Xall,collapse = ' + '),')')), data=job_binary[job_binary$A==0,], family = "binomial", weights=(job_binary$sw*job_binary$ipw)[job_binary$A==0]),job_binary,type="response")
job_binary = job_binary %>% mutate(mu1 = mu1.pred, mu0 = mu0.pred, mu = A*mu1+(1-A)*mu0)
tau.pred = predict(lm(formula=as.formula(paste('mu1-mu0+(2*A-1)*ipw*(Y-mu) ~ (',paste(Xall,collapse = ' + '),')')), data=job_binary,weights=job_binary$sw))
sau.pred = predict(lm(formula=as.formula(paste('mu1+mu0+ipw*(Y-mu) ~ (',paste(Xall,collapse = ' + '),')')), data=job_binary,weights=job_binary$sw))
job_binary = job_binary %>% mutate(tau = tau.pred, sau = sau.pred)
ATE = job_binary %>% mutate(IF = (mu1-mu0+(2*A-1)*ipw*(Y-mu))) %>% summarise(m = mean(sw*IF,na.rm=T), se = sd(sw*IF,na.rm=T)/sqrt(n()))
A1 = job_binary %>% mutate(IF = (mu1+A*ipw*(Y-mu))) %>% summarise(m = mean(sw*IF,na.rm=T), se = sd(sw*IF,na.rm=T)/sqrt(n()))
A0 = job_binary %>% mutate(IF = (mu0+(1-A)*ipw*(Y-mu))) %>% summarise(m = mean(sw*IF,na.rm=T), se = sd(sw*IF,na.rm=T)/sqrt(n()))

rbind(
job_binary %>% mutate(IF = -(tau<0)*(mu1-mu0+(2*A-1)*ipw*(Y-mu))) %>% summarise(m = mean(sw*IF,na.rm=T), se = sd(sw*IF,na.rm=T)/sqrt(n())) %>% mutate(type='wx',bound='lb'),
job_binary %>% mutate(IF = mu0+(1-A)*ipw*(Y-mu)+(sau>=1)*(1-mu1-mu0-ipw*(Y-mu))) %>% summarise(m = mean(sw*IF,na.rm=T), se = sd(sw*IF,na.rm=T)/sqrt(n())) %>% mutate(type='wx',bound='ub'),
job_binary %>% mutate(IF = (tau<0)*((2*A-1)*ipw*(Y))+(1-A)*ipw*(Y)+(sau>=1)*(1-ipw*Y)) %>% summarise(m = mean(sw*IF,na.rm=T), se = sd(sw*IF,na.rm=T)/sqrt(n())) %>% mutate(type='wx',bound='opt'),
data.frame(m=-min(ATE$m,0), se=if(ATE$m<0){ATE$se}else{0}, type='nx',bound='lb'),
data.frame(m=min(A0$m,1-A1$m), se=if(A0$m<1-A1$m){A0$se}else{A1$se}, type='nx',bound='ub'),
data.frame(m=min(A0$m,A1$m,1-A0$m,1-A1$m), se=if(min(A0$m,1-A0$m)<min(A1$m,1-A1$m)){A0$se}else{A1$se}, type='nx',bound='opt')
) %>% mutate(A=AA,B=BB)
}
}

plot.lb = results %>% mutate(A=substr(A,3,5), B=substr(B,3,5)) %>% unite(AB, c("A", "B"), sep='>') %>% mutate(Type=recode_factor(type,wx='With X',nx='No X')) %>% filter(bound=='lb') %>% ggplot() +
  aes(x=AB,y=m,ymin=m-1.96*se,ymax=m+1.96*se, color=Type, linetype=Type, shape=Type) + theme_minimal() + geom_errorbar(size=1) + geom_point(size=5) + theme(axis.title.y=element_blank(), axis.title.x=element_blank())


plot.ub = results %>% mutate(A=substr(A,3,5), B=substr(B,3,5)) %>% unite(AB, c("A", "B"), sep='>') %>% mutate(Type=recode_factor(type,wx='With X',nx='No X')) %>% filter(bound=='ub') %>% ggplot() +
  aes(x=AB,y=m,ymin=m-1.96*se,ymax=m+1.96*se, color=Type, linetype=Type, shape=Type) + theme_minimal() + geom_errorbar(size=1) + geom_point(size=5) + theme(axis.title.y=element_blank(), axis.title.x=element_blank()) + scale_y_continuous(breaks=c(0.22,0.24,0.26))

plot.opt = results %>% mutate(A=substr(A,3,5), B=substr(B,3,5)) %>% unite(AB, c("A", "B"), sep='-') %>% mutate(Type=recode_factor(type,wx='With X',nx='No X')) %>% filter(bound=='opt', AB %in% c('pri-pub','pri-sta','pub-sta')) %>% ggplot() +
  aes(x=AB,y=m,ymin=m-1.96*se,ymax=m+1.96*se, color=Type, linetype=Type, shape=Type) + theme_minimal() + geom_errorbar(size=1) + geom_point(size=5) + theme(axis.title.y=element_blank(), axis.title.x=element_blank()) + scale_y_continuous(breaks=c(0.20,0.22,0.24))

library(gridExtra)
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

ggsave('beh_leg.pdf',plot=g_legend(plot.lb+guides(color = guide_legend(label.position = "top",label.vjust=-5))+theme(legend.title=element_blank(),legend.key.size = unit(1.25, 'cm'))))
ggsave('beh_lb.pdf' ,plot=plot.lb   + theme(legend.position="none", plot.margin=grid::unit(c(0,0,0,0), "mm")), dpi = 300, height = 1.33, width = 4)
ggsave('beh_ub.pdf' ,plot=plot.ub   + theme(legend.position="none", plot.margin=grid::unit(c(0,0,0,0), "mm")), dpi = 300, height = 1.33, width = 4)
ggsave('beh_opt.pdf',plot=plot.opt  + theme(legend.position="none", plot.margin=grid::unit(c(0,0,0,0), "mm")), dpi = 300, height = 1.33, width = 2.25)
