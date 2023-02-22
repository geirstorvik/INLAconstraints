
# Script taken from Rachel Lowe (2021) and altered by Aanes and Storvik. Only parts of the scripts are used.
# https://github.com/drrachellowe/hydromet_dengue

# Step 0: load packages and pre-processed data
# Step 1: formulate a baseline model on a subset of the data, including spatiotemporal random effects
# Step 2: fit model on smaller dataset, 8 years of data using data for month 6 only

#NB: We remove one region to get a fully connected graph.

# Step 0: load packages and pre-processed data
rm(list=ls())
library(INLA)
INLA::inla.setOption("pardiso.license","~/sys/licences/pardiso.lic")
INLA::inla.pardiso.check()
library(SpaceTimePaper)
library(Matrix)
library(sf)
library(spdep)
library(data.table)
library(dlnm)
library(ggplot2)
library(xtable)


nt=6

objects_to_be_used2=dengue_data_set_up(nt=nt,month=6)

objects_to_be_used = readRDS("exec/objects_to_be_used.RDS")

#We fit a simplified model first on a smaller subset of the data, using six years of data.
df_check=objects_to_be_used$data_check
df= objects_to_be_used$data_check
df_month=objects_to_be_used$df_month
Q_s1=objects_to_be_used$Q_s1
ns=nrow(Q_s1)


Q_RW2_before=GMRF_RW(n=nt,order=2)
#Scale model
A_t=rbind(matrix(rep(1,nt),nrow=1),seq(1,nt))
Q_RW2 = inla.scale.model(Q_RW2_before,list(A=A_t,e=rep(0,2)))

#Q_RW2_scaled=inla.scale.model(GMRF_RW(order=2,n=10),constr = list(A=rbind(rep(1,nt),seq(1,nt)),e=rep(0,2)))
Q_ICAR=inla.scale.model(Q_s1,constr=list(A=matrix(rep(1,ns),nrow=1),e=0))

Q_st=kronecker(Q_RW2,Q_ICAR)
source("R/SpaceTimeProjConstr.R")
PC = SpaceTimeProjConstr(ns,nt)

df_month$T3 = df_month$T2
df_month$S2 = df_month$S1

eps = 1e-05
kap = 1e06
baseformula =Y~offset(log(E))+Vu+
  f(T2,model="generic0",Cmatrix=Q_RW2+diag(nt)*eps,constr=T)+
  #f(T3,model="iid")+
  f(S1,model="generic0",Cmatrix=Q_ICAR+diag(ns)*eps,constr=T)+
  #f(S2,model="iid")+
  f(S1T1,model="generic0",Cmatrix=Q_st+diag(ns*nt)*eps,constr=F,
    extraconstr=list(A=as.matrix(PC$A),e=rep(0,nrow(PC$A))))


baseformula.proj=Y~offset(log(E))+Vu+
  f(T2,model="generic0",Cmatrix=Q_RW2+diag(nt)*eps,constr=T)+
  #f(T3,model="iid")+
  f(S1,model="generic0",Cmatrix=Q_ICAR+diag(ns)*eps,constr = T) +
  #f(S2,model="iid")+
  f(S1T1,Cmatrix=Q_st+diag(ns*nt)*eps,precision=kap,model="z",Z=as.matrix(PC$P),constr=F,
    extraconstr=list(A=as.matrix(PC$A2),e=rep(0,nrow(PC$A2))))

run.goicoa=FALSE
if(!file.exists("dengue.goicoa.proj.RDS") | !file.exists("dengue.goicoa.RDS"))
  run.goicoa = TRUE
if(run.goicoa)
{
  dengue.goicoa.proj=inla(baseformula.proj, family = "nbinomial",data =df_month,num.threads =4,inla.mode="experimental",
                         control.fixed = list(
                           prec.intercept =0.01),verbose=T,control.inla=list(strategy="gaussian" ))
  saveRDS(dengue.goicoa.proj,file="dengue.goicoa.proj.RDS")
  
  dengue.goicoa=inla(baseformula, family = "nbinomial",data =df_month,num.threads =4,inla.mode="experimental",
                  control.fixed = list(
                    prec.intercept =0.01),verbose=T,control.inla=list(strategy="gaussian" ))
  saveRDS(dengue.goicoa,file="dengue.goicoa.RDS")
}
if(!run.goicoa)
{
  dengue.goicoa.proj=readRDS("dengue.goicoa.proj.RDS")
  dengue.goicoa=readRDS("dengue.goicoa.RDS")
}

show(c(dengue.goicoa$cpu.used[4],dengue.goicoa.proj$cpu.used[4]))

#Table for latex file
tab.goicoa = rbind(dengue.goicoa$summary.fixed[,1:2],
                 dengue.goicoa$summary.hyperpar[,1:2])
tab.goicoa.proj = rbind(dengue.goicoa.proj$summary.fixed[,1:2],
                 dengue.goicoa.proj$summary.hyperpar[,1:2])
tab.goicoa2 = cbind(tab.goicoa,tab.goicoa.proj)
#rownames(tab.goicoa2) = c("mu","Vu","Disp","tau_alpha","tau_T_iid","tau_theta","tau_S_iid","tau_delta")
rownames(tab.goicoa2) = c("mu","Vu","Disp","tau_alpha","tau_theta","tau_delta")
xtable(tab.goicoa2,digits=3)
saveRDS(tab.goicia2,file="dengue_tab_goicoa2.RDS")

#Plotting interaction terms, mean and standard deviations
plotData=data.table::data.table(StandardModelE=dengue.goicoa$summary.random$S1T2$mean,
                                NewParametrizationE=dengue.goicoa.proj$summary.random$S1T2$mean[1:(ns*nt)],
                                StandardModelSD=dengue.goicoa$summary.random$S1T2$sd,
                                NewParametrizationSD=dengue.goicoa.proj$summary.random$S1T2$sd[1:(ns*nt)])

ggplot(data=plotData)+geom_point(aes(y=StandardModelE,x=NewParametrizationE),colour="red",size=1.25)+xlab("Estimated mean new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+
  ggtitle("Estimated mean standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))
ggsave("Dengue_EstimatedMeanSimulatedDateGoicoa.pdf",height=5,width=5)
ggplot(data=plotData)+geom_point(aes(y=StandardModelSD,x=NewParametrizationSD),colour="red",size=1.25)+xlab("Estimated sd new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+
  ggtitle("Estimated sd standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))#+
ggsave("Dengue_EstimatedSDSimulatedDataGoicoa.pdf",height=5,width=5)

#Comparison of marginal likelihoods
#Need to correct for the standard method

C2 = forceSymmetric(as(Q_st+diag(ns*nt)*eps,"sparseMatrix"))
A = PC$A
QA = solve(C2,t(A))
AQA = forceSymmetric(solve(t(QA)%*%t(A)))
Sig.cond = forceSymmetric(solve(C2)-QA%*%AQA%*%t(QA))
val = eigen(Sig.cond)$val
np = diff(dim(A))
logdet = sum(log(val[1:np])) 
show(c(dengue.goicoa$mlik[1,1]-0.5*logdet,dengue.goicoa.proj$mlik[1,1]))

