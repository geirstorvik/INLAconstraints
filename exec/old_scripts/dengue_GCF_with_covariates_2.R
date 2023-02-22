
# Script taken from Rachel Lowe (2021) and altered by Aanes and Storvik. Only parts of the scripts are used.
# https://github.com/drrachellowe/hydromet_dengue

# Step 0: load packages and pre-processed data
# Step 1: formulate a baseline model on a subset of the data, including spatiotemporal random effects
# Step 2: fit model on smaller dataset, 8 years of data using data for month 6 only

#NB: We remove one region to get a fully connected graph.

# Step 0: load packages and pre-processed data
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
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


# Script taken from Rachel Lowe (2021) and altered by Aanes and Storvik. Only parts of the scripts are used.
# https://github.com/drrachellowe/hydromet_dengue

# Step 0: load packages and pre-processed data
# Step 1: formulate a baseline model on a subset of the data, including spatiotemporal random effects
# Step 2: fit model on smaller dataset, 7 years of data using data for all months

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


nt=7

#objects_to_be_used=dengue_data_set_up(nt=nt)
objects_to_be_used = readRDS("../data/objects_to_be_used.RDS")

#objects_to_be_used = readRDS("exec/objects_to_be_used.RDS")

#We fit a simplified model first on a smaller subset of the data, using six years of data.
df=objects_to_be_used$df
basis_tmin=objects_to_be_used$basis_tmin
basis_pdsi=objects_to_be_used$basis_pdsi
urban_basis1_pdsi=objects_to_be_used$urban_basis1_pdsi
Q_s1=objects_to_be_used$Q_s1
ns=nrow(Q_s1)
df=df[df$T2<=nt,]
df$S1T2_iid=df$S1T2

Q_RW2_before=GMRF_RW(n=nt,order=2)
#Scale model
A_t=rbind(matrix(rep(1,nt),nrow=1),seq(1,nt))
Q_RW2 = inla.scale.model(Q_RW2_before,list(A=A_t,e=rep(0,2)))

#Q_RW2_scaled=inla.scale.model(GMRF_RW(order=2,n=10),constr = list(A=rbind(rep(1,nt),seq(1,nt)),e=rep(0,2)))
Q_ICAR=inla.scale.model(Q_s1,constr=list(A=matrix(rep(1,ns),nrow=1),e=0))

Q_st=kronecker(Q_RW2,Q_ICAR)
source("../R/SpaceTimeProjConstr.R")

PC = SpaceTimeProjConstr(ns,nt,type ="GCF")

#Make graphs for the temporal components
library(Matrix)
graph.ar = function(nt)
{
  Adj = diag(c(1,rep(2,nt-2),1))
  for(i in 2:nt)
  {
    Adj[i,i-1] = -1
    Adj[i-1,i] = -1
  }
  as(Adj,"sparseMatrix")
}
graph.T1 = graph.ar(12)
graph.T2 = graph.ar(nt)
Q_s1 = as(Q_s1,"sparseMatrix")



baseformula <- Y ~ offset(log(E)) + basis_tmin + basis_pdsi+ urban_basis1_pdsi + Vu+
  #f(T1,model="bym2",constr=TRUE,graph=graph.T1,scale.model=T)+
  f(T2,model="bym2",diagonal=eps,graph=graph.T2,constr=T) +
  f(S1,model="bym2",diagonal=eps,graph=Q_s1,constr=T)+
  f(S1T2,model="generic0",Cmatrix = Q_st+diag(ns*nt)*eps,constr=F,
    extraconstr = list(A=as.matrix(PC$A),e=rep(0,nrow(PC$A))))#+f(S1T2_iid,model="iid")


baseformula.proj <- Y~ offset(log(E)) + basis_tmin + basis_pdsi+ urban_basis1_pdsi + Vu+
  #f(T1,model="bym2",constr=TRUE,graph=graph.T1,scale.model=T)+
  f(T2,model="bym2",diagonal=eps,graph=graph.T2,constr=T) +
  f(S1,model="bym2",diagonal=eps,graph=Q_s1,constr=T)+
  f(S1T2,model="z",Z=as.matrix(PC$P),precision=kap,Cmatrix=Q_st+diag(ns*nt)*eps,constr=F,
    extraconstr=list(A=as.matrix(PC$A2),e=rep(0,nrow(PC$A2))))#+f(S1T2_iid,model="iid")



eps = 1e-05
kap = 1e06
run.goicoa=FALSE
if(!file.exists("dengue.goicoa2.proj.RDS") | !file.exists("dengue.goicoa2.RDS"))
  run.goicoa = TRUE
if(run.goicoa)
{
  dengue.goicoa2.proj=inla(baseformula.proj, family = "poisson",num.threads =4,inla.mode="experimental",
                         control.fixed = list(
                           prec.intercept =0.01),verbose=T,
                         #control.inla=list(strategy="gaussian" ),
                         data=df)
  saveRDS(dengue.goicoa2.proj,file="dengue.goicoa2.proj.RDS")

  dengue.goicoa2=inla(baseformula, family = "poisson",data =df,num.threads =4,inla.mode="experimental",
                  control.fixed = list(
                    prec.intercept =0.01),verbose=T
                  #,control.inla=list(strategy="gaussian" )
                  )
  saveRDS(dengue.goicoa2,file="dengue.goicoa2.RDS")
}
if(!run.goicoa)
{
  dengue.goicoa2.proj=readRDS("dengue.goicoa2.proj.RDS")
  dengue.goicoa2=readRDS("dengue.goicoa2.RDS")
}

show(c(dengue.goicoa2$cpu.used[4],dengue.goicoa2.proj$cpu.used[4]))

#Table for latex file
tab.goicoa = rbind(dengue.goicoa2$summary.fixed[,1:2],
                 dengue.goicoa2$summary.hyperpar[,1:2])
tab.goicoa.proj = rbind(dengue.goicoa2.proj$summary.fixed[,1:2],
                 dengue.goicoa2.proj$summary.hyperpar[,1:2])
tab.goicoa2 = cbind(tab.goicoa,tab.goicoa.proj)
rownames(tab.goicoa2) = c("mu","Vu","Disp","tau_alpha","tau_T_iid","tau_theta","tau_S_iid","tau_delta")
rownames(tab.goicoa2) = c("mu","Vu","Disp","tau_alpha","tau_theta","tau_delta")
xtable(tab.goicoa2,digits=3)
saveRDS(tab.goicia2,file="dengue_tab_goicoa2.RDS")

#Plotting interaction terms, mean and standard deviations
plotData=data.table::data.table(StandardModelE=dengue.goicoa2$summary.random$S1T2$mean,
                                NewParametrizationE=dengue.goicoa2.proj$summary.random$S1T2$mean[1:(ns*nt)],
                                StandardModelSD=dengue.goicoa2$summary.random$S1T2$sd,
                                NewParametrizationSD=dengue.goicoa2.proj$summary.random$S1T2$sd[1:(ns*nt)])

ggplot(data=plotData)+geom_point(aes(y=StandardModelE,x=NewParametrizationE),colour="red",size=1.25)+xlab("Estimated mean new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+
  ggtitle("Estimated mean standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))
  ggsave("Dengue_EstimatedMeanSimulatedDateGoicoa2.pdf")
ggplot(data=plotData)+geom_point(aes(y=StandardModelSD,x=NewParametrizationSD),colour="red",size=1.25)+xlab("Estimated sd new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+
  ggtitle("Estimated sd standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))#+
ggsave("Dengue_EstimatedSDSimulatedDataGoicoa2.pdf")



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
show(c(dengue.goicoa2$mlik[1,1]-0.5*logdet,dengue.goicoa2.proj$mlik[1,1]))

