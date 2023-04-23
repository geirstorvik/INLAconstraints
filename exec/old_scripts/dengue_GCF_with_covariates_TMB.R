
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

objects_to_be_used=dengue_data_set_up(nt=nt)

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
source("R/SpaceTimeProjConstr.R")

PC2 = SpaceTimeProjConstr(ns,nt,type ="GCF")

eps=1e-05
Z_st=create_designmatrix_random_effect(Q_st,PC2$A,eps=eps)
Z_s=create_designmatrix_random_effect(Q_ICAR,matrix(1,ncol=ns,nrow=1),eps=eps)
Z_t=create_designmatrix_random_effect(Q_RW2,matrix(1,ncol=nt,nrow=1),eps=eps)


Z_s=kronecker(rep(1,nt),diag(ns))%*%Z_s

Z_t=kronecker(diag(nt),rep(1,ns))%*%Z_t

#Z_st = diag(ns*nt)*Z_st

rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

Z_s_repeated=NULL
for(j in 1:12){
  print(j)
  Z_s_repeated=rbind(Z_s_repeated,Z_s)
}


Z_t_repeated=NULL
for(j in 1:12){
  print(j)

  Z_t_repeated=rbind(Z_t_repeated,Z_t)}


Z_st_repeated=NULL
for(j in 1:12){
  print(j)
  Z_st_repeated=rbind(Z_st_repeated,Z_st)
}

data=list(y=df$Y,
          Q_t=Q_RW2+diag(nt)*eps,
          Q_s=Q_ICAR+diag(ns)*eps,
          Q_st=Q_st+diag(ns*nt)*eps,
          Z_s=as(Z_s_repeated,"dgTMatrix"),
          Z_t=as(Z_t_repeated,"dgTMatrix"),
          Z_st=Z_st_repeated,
          X=matrix(1,nrow=nrow(df),ncol=1),
          m_s=ns-1,
          m_t=nt-2,
          m_st=(ns-1)*(nt-2),
          n_data=nrow(df),
          shape=1,
          scale=1,
          fixed_sd=1)
parameters=list(xs=rep(0,ns),
                xt=rep(0,nt),
                xst=rep(0,ns*nt),
                log_tau_s=0,
                log_tau_st=0,
                log_tau_t=0,
                beta_fixed=rep(0,ncol(data$X)))


library(TMB)
dyn.load(dynlib("C:/Users/fredr/OneDrive/Dokumenter/Paper/INLAconstraints/exec/SpatioTemporalModelTMB_sparse"))

obj = MakeADFun(data,parameters,random=c("xst","xs","xt"),DLL="SpatioTemporalModelTMB_sparse")
opt = nlminb(obj$par, obj$fn,obj$gr)

