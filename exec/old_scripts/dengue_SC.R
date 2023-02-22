
# Script taken from Rachel Lowe (2021) and altered by Aanes and Storvik. Only parts of the scripts are used.
# https://github.com/drrachellowe/hydromet_dengue

# Step 0: load packages and pre-processed data
# Step 1: formulate a baseline model on a subset of the data, including spatiotemporal random effects
# Step 2: fit model on smaller dataset, 8 years of data using data for month 6 only

#NB: We remove one region to get a fully connected graph.

# Step 0: load packages and pre-processed data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
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

#objects_to_be_used=dengue_data_set_up(nt=nt)

objects_to_be_used = readRDS("../data/objects_to_be_used.RDS")


#We fit a simplified model first on a smaller subset of the data, using six years of data.
df_check=objects_to_be_used$data_check
df= objects_to_be_used$data_check
df_month=objects_to_be_used$df_month
Q_s1=objects_to_be_used$Q_s1
ns=nrow(Q_s1)

df$basis_tmin=objects_to_be_used$basis_tmin
df$basis_pdsi=objects_to_be_used$basis_pdsi
df$urban_basis1_pdsi=objects_to_be_used$urban_basis1_pdsi



Q_RW2_before=GMRF_RW(n=nt,order=2)
#Scale model
A_t=rbind(matrix(rep(1,nt),nrow=1),seq(1,nt))
Q_RW2 = inla.scale.model(Q_RW2_before,list(A=A_t,e=rep(0,2)))

#Q_RW2_scaled=inla.scale.model(GMRF_RW(order=2,n=10),constr = list(A=rbind(rep(1,nt),seq(1,nt)),e=rep(0,2)))
Q_ICAR=inla.scale.model(Q_s1,constr=list(A=matrix(rep(1,ns),nrow=1),e=0))

Q_st=kronecker(Q_RW2,Q_ICAR)
source("../R/SpaceTimeProjConstr.R")

#Knorr-Held:

PC2 = SpaceTimeProjConstr(ns,nt,type="SC")

prior.fixed=list(prec.intercept=0.001)

eps = 1e-05
kap = 1e06
baseformula<- Y ~offset(log(E)) + Vu + 
  f(S1,model="generic0",diagonal=eps,Cmatrix = Q_s1,constr=T)+
  f(T2,model="generic0",Cmatrix=Q_RW2+diag(nt)*eps,constr=T) +
  f(S1T1,model="generic0",Cmatrix = Q_st+diag(ns*nt)*eps,constr=F,
    extraconstr = list(A=as.matrix(PC2$A),e=rep(0,nrow(PC2$A))))

baseformula.proj <- Y ~ offset(log(E)) + Vu + 
  f(S1,model="generic0",diagonal=eps,Cmatrix=Q_s1,constr=T) +
  f(T2,model="generic0",Cmatrix=Q_RW2+diag(nt)*eps,constr=T) +
  f(S1T1,model="z",Z=as.matrix(PC2$P),precision=kap,Cmatrix=Q_st+diag(ns*nt)*eps,constr=F,
    extraconstr=list(A=as.matrix(PC2$A2),e=rep(0,nrow(PC2$A2))))

run.knorr = FALSE
if(!file.exists("dengue.proj.nb.knorr.RDS") | !file.exists("dengue.nb.knorr.RDS"))
  run.knorr = TRUE
if(run.knorr)
{
 dengue.proj.nb.knorr <- inla(baseformula.proj, family = "nbinomial",data =df_month,num.threads =4,inla.mode="experimental",control.fixed = list(
   prec.intercept =0.001),verbose=T,control.inla=list(cmin=0 ))
 saveRDS(dengue.proj.nb.knorr,"dengue.proj.nb.knorr.RDS")

 dengue.nb.knorr <- inla(baseformula, family = "nbinomial",data=df_month,num.threads = 4,control.fixed = list(
   prec.intercept =0.001),control.inla=list(cmin=0),verbose=T,inla.mode="experimental")
 saveRDS(dengue.nb.knorr,"dengue.nb.knorr.RDS")
}
if(!run.knorr)
{
  dengue.proj.nb.knorr = readRDS("dengue.proj.nb.knorr.RDS")
  dengue.nb.knorr = readRDS("dengue.nb.knorr.RDS")
}

show(c(dengue.nb.knorr$cpu.used[4],dengue.proj.nb.knorr$cpu.used[4]))

#Table for latex file
tab.knorr = rbind(dengue.nb.knorr$summary.fixed[,1:2],
                   dengue.nb.knorr$summary.hyperpar[,1:2])
tab.knorr.proj = rbind(dengue.proj.nb.knorr$summary.fixed[,1:2],
                        dengue.proj.nb.knorr$summary.hyperpar[,1:2])
tab.knorr2 = cbind(tab.knorr,tab.knorr.proj)
rownames(tab.knorr2) = c("mu","Vu","Disp","tau_alpha","tau_theta","tau_delta")
xtable(tab.knorr2,digits=3)
saveRDS(tab.knorr2,file="dengue_tab_knorr2.RDS")

tab = cbind(tab.goicoa2,tab.knorr2)
xtable(tab,digits=3)

#Plotting interaction terms, mean and standard deviations
plotData=data.table::data.table(StandardModelE=dengue.nb.knorr$summary.random$S1T1$mean,
                                NewParametrizationE=dengue.proj.nb.knorr$summary.random$S1T1$mean[1:(ns*nt)],
                                StandardModelSD=dengue.nb.knorr$summary.random$S1T1$sd,
                                NewParametrizationSD=dengue.proj.nb.knorr$summary.random$S1T1$sd[1:(ns*nt)])

ggplot(data=plotData)+geom_point(aes(y=StandardModelE,x=NewParametrizationE),colour="red",size=1.25)+xlab("Estimated mean new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+ggtitle("Estimated mean standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))+
  coord_fixed()+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+labs(x=NULL, y=NULL)
ggsave("Dengue_EstimatedMeanSimulatedData_Knorr.pdf")
ggplot(data=plotData)+geom_point(aes(y=StandardModelSD,x=NewParametrizationSD),colour="red",size=1.25)+xlab("Estimated sd new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+ggtitle("Estimated sd standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))+
  coord_fixed()+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
ggsave("Dengue_EstimatedSDSimulatedData_Knorr.pdf")

#Comparison of marginal likelihoods
#Need to correct for the standard method

C2 = forceSymmetric(as(Q_st+diag(ns*nt)*eps,"sparseMatrix"))
A = PC2$A
QA = solve(C2,t(A))
AQA = forceSymmetric(solve(t(QA)%*%t(A)))
Sig.cond = forceSymmetric(solve(C2)-QA%*%AQA%*%t(QA))
val = eigen(Sig.cond)$val
np = diff(dim(A))
logdet = sum(log(val[1:np])) 
#dengue.nb.knorr$mlik= dengue.nb.knorr$mlik
show(c(dengue.nb.knorr$mlik[1,1]-0.5*logdet,dengue.proj.nb.knorr$mlik[1,1]))
