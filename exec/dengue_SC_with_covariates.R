# Script taken from Rachel Lowe (2021) and altered by Aanes and Storvik. Only parts of the scripts are used.
# https://github.com/drrachellowe/hydromet_dengue

# Step 0: load packages and pre-processed data
# Step 1: formulate a baseline model on a subset of the data, including spatiotemporal random effects
# Step 2: fit model on smaller dataset, 7 years of data using data for all months

#NB: We remove one region to get a fully connected graph.

# Step 0: load packages and pre-processed data
#rm(list=ls())
library(INLA)
INLA::inla.setOption("pardiso.license","~/sys/licences/pardiso.lic")
INLA::inla.pardiso.check()
library(INLAconstraints)
library(Matrix)
library(sf)
library(spdep)
library(data.table)
library(dlnm)
library(ggplot2)
library(xtable)


nt=7

data(denguedata)


#We fit a simplified model first on a smaller subset of the data, using six years of data.
df=denguedata$df
basis_tmin=denguedata$basis_tmin
basis_pdsi=denguedata$basis_pdsi
urban_basis1_pdsi=denguedata$urban_basis1_pdsi
Q_ICAR=denguedata$Q_s1
ns=nrow(Q_ICAR)
df=df[df$T2<=nt,]

order = 1
Q_RW=GMRF_RW(n=nt,order=order)
SCALED = FALSE
if(SCALED)
{
  #Scale model
  A_t=rbind(matrix(rep(1,nt),nrow=1),seq(1,nt))
  Q_RW = inla.scale.model(Q_RW,list(A=A_t,e=rep(0,2)))
  Q_ICAR=inla.scale.model(Q_ICAR,constr=list(A=matrix(rep(1,ns),nrow=1),e=0))
  
}

Q_st=kronecker(Q_RW,Q_ICAR)

PC = SpaceTimeProjConstr(ns,nt,type ="SC")

df2 = cbind(df,basis_tmin)

eps = 1e-05
kap = 1e06

baseformula <- Y ~ offset(log(E)) + Vu + 
  basis_tmin.v1.l1 + basis_tmin.v1.l2 + basis_tmin.v1.l3+
  basis_tmin.v2.l1 + basis_tmin.v2.l2 + basis_tmin.v2.l3+
  basis_tmin.v3.l1 + basis_tmin.v3.l2 + basis_tmin.v3.l3+
  f(T2,model="generic0",Cmatrix=Q_RW+diag(nt)*eps) +
  #f(T2,model="rw1",constr=T,cyclic=TRUE) +
  f(S1,model="generic0",Cmatrix = Q_ICAR+diag(ns)*eps,constr=T)+
  f(S1T2,model="generic0",Cmatrix = Q_st+diag(ns*nt)*eps,constr=F,
    extraconstr = list(A=as.matrix(PC$A),e=rep(0,nrow(PC$A))))


baseformula.proj <- Y~ offset(log(E)) + Vu + 
  basis_tmin.v1.l1 + basis_tmin.v1.l2 + basis_tmin.v1.l3+
  basis_tmin.v2.l1 + basis_tmin.v2.l2 + basis_tmin.v2.l3+
  basis_tmin.v3.l1 + basis_tmin.v3.l2 + basis_tmin.v3.l3+
  f(T2,model="generic0",Cmatrix=Q_RW+diag(nt)*eps) +
  f(S1,model="generic0",Cmatrix = Q_ICAR+diag(ns)*eps,constr=T)+
  f(S1T2,model="z",Z=as.matrix(PC$P),precision=kap,Cmatrix=Q_st+diag(ns*nt)*eps,constr=F,
    extraconstr=list(A=as.matrix(PC$A2),e=rep(0,nrow(PC$A2))))


run.sc=FALSE
if(!file.exists("dengue.sc.proj.RDS") | !file.exists("dengue.sc.RDS"))
  run.sc = TRUE
if(run.sc)
{
  dengue.sc.proj=inla(baseformula.proj, family = "nbinomial",data =df2,num.threads =10, 
                      inla.mode="experimental",
                      control.fixed = list(
                      prec.intercept =0.01),verbose=F,
                      control.inla=list(strategy="gaussian" ))
   saveRDS(dengue.sc.proj,file="dengue.sc.proj.RDS")

  dengue.sc=inla(baseformula, family = "nbinomial",data =df2,num.threads =10, 
                 inla.mode="experimental",
                 control.fixed = list(
                 prec.intercept =0.01),verbose=F,
                 control.inla=list(strategy="gaussian" ))
   saveRDS(dengue.sc,file="dengue.sc.RDS")
} 
if(!run.sc)
{
  dengue.sc.proj=readRDS("dengue.sc.proj.RDS")
  dengue.sc=readRDS("dengue.sc.RDS")
}

show(c(dengue.sc$cpu.used[4],dengue.sc.proj$cpu.used[4]))

#Plotting interaction terms, mean and standard deviations
plotData=data.table::data.table(StandardModelE=dengue.sc$summary.random$S1T2$mean,
                                NewParametrizationE=dengue.sc.proj$summary.random$S1T2$mean[1:(ns*nt)],
                                StandardModelSD=dengue.sc$summary.random$S1T2$sd,
                                NewParametrizationSD=dengue.sc.proj$summary.random$S1T2$sd[1:(ns*nt)])

plotData=data.table::data.table(StandardModelE=dengue.sc$summary.random$T2$mean,
                                NewParametrizationE=dengue.sc.proj$summary.random$T2$mean,
                                StandardModelSD=dengue.sc$summary.random$T2$sd,
                                NewParametrizationSD=dengue.sc.proj$summary.random$T2$sd)

ggplot(data=plotData)+geom_point(aes(y=StandardModelE,x=NewParametrizationE),colour="red",size=1.25)+xlab("Estimated mean new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+
  ggtitle("Estimated mean standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))
  ggsave("Dengue_EstimatedMeanSimulatedDatesc.pdf",height=5,width=5)
ggplot(data=plotData)+geom_point(aes(y=StandardModelSD,x=NewParametrizationSD),colour="red",size=1.25)+xlab("Estimated sd new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+
  ggtitle("Estimated sd standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))#+
ggsave("Dengue_EstimatedSDSimulatedDatasc.pdf",height=5,width=5)

#Comparison of marginal likelihoods
#Need to correct for the standard method
mcorS = LikCorrectGeneric0(Q_ICAR,matrix(rep(1,ns),nrow=1),eps)
mcorT = LikCorrectGeneric0(Q_RW,matrix(rep(1,nt),nrow=1),eps)
mcorST = LikCorrectGeneric0(Q_st,as.matrix(PC$A),eps)
show(c(dengue.sc$mlik[1,1]+mcorS+mcorT+mcorST,
       dengue.sc.proj$mlik[1,1]+mcorS+mcorT)/nrow(df))

#Table for latex file
tab.sc = rbind(dengue.sc$summary.fixed[,1:2],
               dengue.sc$summary.hyperpar[,1:2])
tab.sc.proj = rbind(dengue.sc.proj$summary.fixed[,1:2],
                    dengue.sc.proj$summary.hyperpar[,1:2])
tab.sc2 = cbind(tab.sc,tab.sc.proj)
show(round(tab.sc2,3))
rownames(tab.sc2) = c("mu","Vu","Disp","tau_alpha","tau_T_iid","tau_theta","tau_S_iid","tau_delta")
rownames(tab.sc2) = c("mu","Vu","Disp","tau_alpha","tau_theta","tau_delta")
xtable(tab.sc2,digits=3)
saveRDS(tab.sc2,file="dengue_tab_sc2.RDS")

foo1 = dengue.sc$marginals.hyperpar$`Precision for T2`
foo2 = dengue.sc.proj$marginals.hyperpar$`Precision for T2`
matplot(cbind(foo1[,1],foo2[,1]),cbind(foo1[,2],foo2[,2]),type="l")


#Table for latex file
tab.sc = rbind(dengue.sc$summary.fixed[,1:2],
               dengue.sc2$summary.fixed[,1:2])
tab.sc.proj = rbind(dengue.sc.proj$summary.fixed[,1:2],
                    dengue.sc.proj2$summary.fixed[,1:2])
tab.sc2 = cbind(tab.sc,tab.sc.proj)

