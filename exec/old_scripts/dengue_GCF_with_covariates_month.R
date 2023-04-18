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
data(denguedata)

df=denguedata$df
basis_tmin=denguedata$basis_tmin
basis_pdsi=denguedata$basis_pdsi
urban_basis1_pdsi=denguedata$urban_basis1_pdsi
Q_ICAR=denguedata$Q_s1
ns=nrow(Q_ICAR)
df=df[df$T2==1,]
df=df[df$T1<7,]
nt = max(df$T1)

#RW(2) for time
Q_RW2=GMRF_RW(n=nt,order=1)

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

df2 = df
df2 = cbind(df2,basis_tmin[1:nrow(df2),])

eps = 1e-05
kap = 1e06

baseformula <- Y ~ offset(log(E)) + Vu + #basis_tmin + basis_pdsi+ urban_basis1_pdsi + Vu+
  basis_tmin.v1.l1 + basis_tmin.v1.l2 + basis_tmin.v1.l3+
  basis_tmin.v2.l1 + basis_tmin.v2.l2 + basis_tmin.v2.l3+
  basis_tmin.v3.l1 + basis_tmin.v3.l2 + basis_tmin.v3.l3+
  f(T1,model="generic0",Cmatrix=Q_RW2+diag(nt)*eps,constr=T) +
  f(S1,model="generic0",Cmatrix = Q_ICAR+diag(ns)*eps,constr=T)+
  f(S1T1,model="generic0",Cmatrix = Q_st+diag(ns*nt)*eps,constr=F,
    extraconstr = list(A=as.matrix(PC$A),e=rep(0,nrow(PC$A))))


baseformula.proj <- Y~ offset(log(E)) + Vu + #basis_tmin + basis_pdsi+ urban_basis1_pdsi + Vu+
  basis_tmin.v1.l1 + basis_tmin.v1.l2 + basis_tmin.v1.l3+
  basis_tmin.v2.l1 + basis_tmin.v2.l2 + basis_tmin.v2.l3+
  basis_tmin.v3.l1 + basis_tmin.v3.l2 + basis_tmin.v3.l3+
 f(T1,model="generic0",Cmatrix=Q_RW2+diag(nt)*eps,constr=T) +
  f(S1,model="generic0",Cmatrix = Q_ICAR+diag(ns)*eps,constr=T)+
  f(S1T1,model="z",Z=as.matrix(PC$P),precision=kap,Cmatrix=Q_st+diag(ns*nt)*eps,constr=F,
    extraconstr=list(A=as.matrix(PC$A2),e=rep(0,nrow(PC$A2))))


dengue.goicoa.proj=inla(baseformula.proj, family = "nbinomial",data=df2,num.threads =10, inla.mode="experimental",
                         control.fixed = list(prec.intercept =0.01),verbose=T,
                         control.inla=list(strategy="gaussian"),
                        control.compute=list(config=TRUE))
saveRDS(dengue.goicoa.proj,file="dengue.goicoa.proj.RDS")

dengue.goicoa=inla(baseformula, family = "nbinomial",data =df2,num.threads =10, inla.mode="experimental",
                         control.fixed = list(prec.intercept =0.01),verbose=T,
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
show(round(tab.goicoa2,3))
rownames(tab.goicoa2) = c("mu","Vu","Disp","tau_alpha","tau_T_iid","tau_theta","tau_S_iid","tau_delta")
rownames(tab.goicoa2) = c("mu","Vu","Disp","tau_alpha","tau_theta","tau_delta")
xtable(tab.goicoa2,digits=3)
saveRDS(tab.goicia2,file="dengue_tab_goicoa2.RDS")

plotData=data.table::data.table(StandardModelE=dengue.goicoa$summary.random$T2$mean,
                                NewParametrizationE=dengue.goicoa.proj$summary.random$T2$mean,
                                StandardModelSD=dengue.goicoa$summary.random$T2$sd,
                                NewParametrizationSD=dengue.goicoa.proj$summary.random$T2$sd)

ggplot(data=plotData)+geom_point(aes(y=StandardModelE,x=NewParametrizationE),colour="red",size=1.25)+xlab("Estimated mean new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+
  ggtitle("Estimated mean standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))

foo = inla.posterior.sample(n=1,dengue.goicoa.proj)[[1]]$latent
foo2 = foo[grep("T2",rownames(foo)),,drop=FALSE]
x3 = foo2[1:nt,]

tau = seq(0.01,4,length=43)
l = rep(NA,length(tau))
x = matrix((dengue.goicoa$summary.random$T2$mean),ncol=1)
l2 = rep(NA,length(tau))
x2 = matrix((dengue.goicoa.proj$summary.random$T2$mean),ncol=1)
l3 = rep(NA,length(tau))
for(i in 1:length(tau))
{
  l[i] = exp(0.5*(nt-2)*log(tau[i])-0.5*tau[i]*t(x)%*%Q_RW2%*%x)
  l2[i] = exp(0.5*(nt-2)*log(tau[i])-0.5*tau[i]*t(x2)%*%Q_RW2%*%x2)
  l3[i] = exp(0.5*(nt-2)*log(tau[i])-0.5*tau[i]*t(x3)%*%Q_RW2%*%x3)
}
d = diff(tau)[1]
l = l/(sum(l)*d)
l2 = l2/(sum(l2)*d)
l3 = l3/(sum(l3)*d)
foo1 = dengue.goicoa$marginals.hyperpar$`Precision for T2`
foo2 = dengue.goicoa.proj$marginals.hyperpar$`Precision for T2`
matplot(cbind(tau,tau,tau,foo1[,1],foo2[,1]),
        cbind(l,l2,l3,foo1[,2],foo2[,2]),type="l",lty=1,col=1:5,lwd=1.5)
legend("topright",c("l1","l2","l3","stand","proj"),lty=1,col=1:5,lwd=1.5)

foo2 = foo[grep("S1T2",rownames(foo)),,drop=FALSE]
x3 = foo2[1:(ns*nt),]
tau = seq(0.1,0.6,length=43)
l = rep(NA,length(tau))
x = matrix(dengue.goicoa$summary.random$S1T2$mean[1:(ns*nt)],ncol=1)
l2 = rep(NA,length(tau))
x2 = matrix(dengue.goicoa.proj$summary.random$S1T2$mean[1:(ns*nt)],ncol=1)
l3 = rep(NA,length(tau))
for(i in 1:length(tau))
{
  l[i] = 0.5*(ns*nt-(2*ns+nt-2))*log(tau[i])-0.5*tau[i]*t(x)%*%Q_st%*%x
  l2[i] = 0.5*(ns*nt-(2*ns+nt-2))*log(tau[i])-0.5*tau[i]*t(x2)%*%Q_st%*%x2
  l3[i] = 0.5*(ns*nt-(2*ns*nt-2))*log(tau[i])-0.5*tau[i]*t(x3)%*%Q_st%*%x3
}
l = exp(l-max(l))
l2 = exp(l2-max(l2))
l3 = exp(l3-max(l3))
d = diff(tau)[1]
l = l/(sum(l)*d)
l2 = l2/(sum(l2)*d)
l3 = l3/(sum(l3)*d)
foo1 = dengue.goicoa$marginals.hyperpar$`Precision for S1T2`
foo2 = dengue.goicoa.proj$marginals.hyperpar$`Precision for S1T2`
matplot(cbind(tau,tau,tau,foo1[,1],foo2[,1]),
        cbind(l,l2,l3,foo1[,2],foo2[,2]),type="l",lty=1,col=1:5,lwd=1.5)
legend("topright",c("l1","l2","l3","stand","proj"),lty=1,col=1:5,lwd=1.5)

tab.goicoa = cbind(rbind(dengue.goicoa$summary.fixed[,1:2],
                         dengue.goicoa$summary.hyperpar[,1:2]),
                   rbind(dengue.goicoa2$summary.fixed[,1:2],
                         dengue.goicoa2$summary.hyperpar[,1:2]))
tab.goicoa.proj = cbind(rbind(dengue.goicoa.proj$summary.fixed[,1:2],
                              dengue.goicoa.proj$summary.hyperpar[,1:2]),
                        rbind(dengue.goicoa.proj2$summary.fixed[,1:2],
                              dengue.goicoa.proj2$summary.hyperpar[,1:2]))
tab.goicoa2 = cbind(tab.goicoa,tab.goicoa.proj)

