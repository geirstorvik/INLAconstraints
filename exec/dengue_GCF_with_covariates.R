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
library(SpaceTimePaper)
library(Matrix)
library(sf)
library(spdep)
library(data.table)
library(dlnm)
library(ggplot2)
library(xtable)

#We fit a simplified model first on a smaller subset of the data, using seven years of data.
nt=7

data(denguedata)

df=denguedata$df
basis_tmin=denguedata$basis_tmin
basis_pdsi=denguedata$basis_pdsi
urban_basis1_pdsi=denguedata$urban_basis1_pdsi
Q_ICAR=denguedata$Q_s1
ns=nrow(Q_ICAR)
df=df[df$T2<=nt,]

#RW(2) for time
Q_RW2=GMRF_RW(n=nt,order=2)

SCALED = FALSE
if(SCALED)
{
  #Scale model
  A_t=rbind(matrix(rep(1,nt),nrow=1),seq(1,nt))
  Q_RW2 = inla.scale.model(Q_RW2,list(A=A_t,e=rep(0,2)))
  Q_ICAR=inla.scale.model(Q_ICAR,constr=list(A=matrix(rep(1,ns),nrow=1),e=0))
  
}
Q_st=kronecker(Q_RW2,Q_ICAR)

PC = SpaceTimeProjConstr(ns,nt,type ="GCF")

baseformula <- Y ~ offset(log(E)) + basis_tmin + basis_pdsi+ urban_basis1_pdsi + Vu+
  f(T2,model="generic0",Cmatrix=Q_RW2+diag(nt)*eps,constr=T) +
  f(S1,model="generic0",diagonal=eps,Cmatrix = Q_ICAR,constr=T)+
  f(S1T2,model="generic0",Cmatrix = Q_st+diag(ns*nt)*eps,constr=F,
    extraconstr = list(A=as.matrix(PC$A),e=rep(0,nrow(PC$A))))


baseformula.proj <- Y~ offset(log(E)) + basis_tmin + basis_pdsi+ urban_basis1_pdsi + Vu+
  f(T2,model="generic0",Cmatrix=Q_RW2+diag(nt)*eps,constr=T) +
  f(S1,model="generic0",diagonal=eps,Cmatrix = Q_ICAR,constr=T)+
  f(S1T2,model="z",Z=as.matrix(PC$P),precision=kap,Cmatrix=Q_st+diag(ns*nt)*eps,constr=F,
    extraconstr=list(A=as.matrix(PC$A2),e=rep(0,nrow(PC$A2))))

eps = 1e-05
kap = 1e06
dengue.goicoa.proj=inla(baseformula.proj, family = "nbinomial",data=df,num.threads =10, inla.mode="experimental",
                         control.fixed = list(prec.intercept =0.01),verbose=F,
                         control.inla=list(strategy="gaussian"))
  saveRDS(dengue.goicoa.proj,file="dengue.goicoa.proj.RDS")

dengue.goicoa=inla(baseformula, family = "nbinomial",data =df,num.threads =10, inla.mode="experimental",
                         control.fixed = list(prec.intercept =0.01),verbose=F,
                     control.inla=list(strategy="gaussian"))
  saveRDS(dengue.goicoa,file="dengue.goicoa.RDS")
  
show(c(dengue.goicoa$cpu.used[4],dengue.goicoa.proj$cpu.used[4]))

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
mcorS = LikCorrectGeneric0(Q_ICAR,matrix(rep(1,ns),nrow=1),eps)
mcorT = LikCorrectGeneric0(Q_RW2,matrix(rep(1,nt),nrow=1),eps)
mcorST = LikCorrectGeneric0(Q_st,as.matrix(PC$A),eps)
show(c(dengue.goicoa$mlik[1,1]+mcorS+mcorT+mcorST,
       dengue.goicoa.proj$mlik[1,1]+mcorS+mcorT)/nrow(df))

#Table for latex file
tab.goicoa = rbind(dengue.goicoa$summary.fixed[,1:2],
                   dengue.goicoa$summary.hyperpar[,1:2])
tab.goicoa.proj = rbind(dengue.goicoa.proj$summary.fixed[,1:2],
                        dengue.goicoa.proj$summary.hyperpar[,1:2])
tab.goicoa2 = cbind(tab.goicoa,tab.goicoa.proj)
rownames(tab.goicoa2) = c("mu","Vu","Disp","tau_alpha","tau_T_iid","tau_theta","tau_S_iid","tau_delta")
rownames(tab.goicoa2) = c("mu","Vu","Disp","tau_alpha","tau_theta","tau_delta")
xtable(tab.goicoa2,digits=3)
saveRDS(tab.goicia2,file="dengue_tab_goicoa2.RDS")


