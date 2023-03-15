##INLA analysis of covid data
library(INLA)
INLA::inla.setOption("pardiso.license","~/sys/licences/pardiso.lic")
INLA::inla.pardiso.check()
library(SpaceTimePaper)
library(Matrix)
library(sf)
library(spdep)
library(data.table)
library(dlnm)
library(xtable)

rm(list=ls())


#Read data and make indices for temporal and interaction terms
d = readRDS("data/covidUK.RDS")
df = d$data
nt = max(df$T1)
ns = max(df$S1)

#Reduce number of time points
nt = 250
df = df[T1<(nt+1)]

df$S1T1 = (df$T1-1)*ns+df$S1


#Make precision matrices
Qs = -d$adj
for(i in 1:ns)
  Qs[i,i] = -sum(Qs[i,-i])
Q_RW2_before=GMRF_RW(n=nt,order=2)

#Scale model
A_t=rbind(matrix(rep(1,nt),nrow=1),seq(1,nt))
Q_RW2 = inla.scale.model(Q_RW2_before,list(A=A_t,e=rep(0,2)))
Q_RW2 = Q_RW2_before
#Q_RW2 = inla.scale.model(Q_RW2_before)

#Q_RW2_scaled=inla.scale.model(GMRF_RW(order=2,n=10),constr = list(A=rbind(rep(1,nt),seq(1,nt)),e=rep(0,2)))
Q_ICAR=inla.scale.model(Qs,constr=list(A=matrix(rep(1,ns),nrow=1),e=0))
Q_ICAR = Qs
#Q_ICAR=inla.scale.model(Qs)




Q_st=kronecker(Q_RW2,Q_ICAR)
#Q_st = inla.scale.model(Q_st,list(A=D.g,e=rep(0,nrow(D.g))))
source("R/SpaceTimeProjConstr.R")
PC = SpaceTimeProjConstr(ns,nt,type="time")

#PC = SpaceTimeProjConstr(ns,nt,type="space")

df$S2 = df$S1
df$T2 = df$T1

eps = 1e-05
kap = 1e06
df$E = log(d$pop[df$S1])
baseformula =cases~offset(E)+
  f(T1,model="generic0",Cmatrix=Q_RW2+diag(nt)*eps,constr=T)+
  f(T2,model="iid")+
  f(S1,model="generic0",Cmatrix=Q_ICAR+diag(ns)*eps,constr=T) +
  f(S2,model="iid")+
   f(S1T1,model="generic0",Cmatrix=Q_st+diag(nrow(Q_st))*eps,constr=F,
     extraconstr=list(A=as.matrix(PC$C),e=rep(0,nrow(PC$C))))

baseformula.proj=cases~offset(E)+
  f(T1,model="generic0",Cmatrix=Q_RW2+diag(nt)*eps,constr=T)+
  f(T2,model="iid")+
  f(S1,model="generic0",Cmatrix=Q_ICAR+diag(ns)*eps,constr = T) +
  f(S2,model="iid")+
  f(S1T1,Cmatrix=Q_st+diag(nrow(Q_st))*eps,precision=kap,model="z",Z=as.matrix(PC$P),constr=F,
    extraconstr=list(A=as.matrix(PC$C1),e=rep(0,nrow(PC$C1))))

#df = df[order(S1T1)]

run.goicoa=FALSE
if(!file.exists("covid.goicoa.proj.RDS") | !file.exists("covid.goicoa.RDS"))
  run.goicoa = TRUE
if(run.goicoa)
{
  covid.goicoa.proj=inla(baseformula.proj, family = "poisson",data =df,num.threads =4,inla.mode="experimental",
                         control.fixed = list(
                           prec.intercept =0.01),verbose=T,control.inla=list(strategy="laplace" ))
  saveRDS(covid.goicoa.proj,file="covid.goicoa.proj.RDS")
  
  covid.goicoa=inla(baseformula, family = "poisson",data =df,num.threads =4,inla.mode="experimental",
                    control.fixed = list(
                      prec.intercept =0.01),verbose=T,control.inla=list(strategy="laplace" ))
  saveRDS(covid.goicoa,file="covid.goicoa.RDS")
}
if(!run.goicoa)
{
  covid.goicoa.proj=readRDS("covid.goicoa.proj.RDS")
  covid.goicoa=readRDS("covid.goicoa.RDS")
}

show(c(covid.goicoa$cpu.used[4],covid.goicoa.proj$cpu.used[4]))

#Table for latex file
tab.goicoa = rbind(covid.goicoa$summary.fixed[,1:2],
                   covid.goicoa$summary.hyperpar[,1:2])
tab.goicoa.proj = rbind(covid.goicoa.proj$summary.fixed[,1:2],
                        covid.goicoa.proj$summary.hyperpar[,1:2])
tab.goicoa2 = cbind(tab.goicoa,tab.goicoa.proj)
rownames(tab.goicoa2) = c("mu","Disp","tau_alpha","tau_theta","tau_delta")
xtable(tab.goicoa2,digits=3)

#Plotting interaction terms, mean and standard deviations
plotData=data.table::data.table(StandardModelE=covid.goicoa$summary.random$S1T1$mean,
                                NewParametrizationE=covid.goicoa.proj$summary.random$S1T1$mean[1:(ns*nt)],
                                StandardModelSD=covid.goicoa$summary.random$S1T1$sd,
                                NewParametrizationSD=covid.goicoa.proj$summary.random$S1T1$sd[1:(ns*nt)])
plot(plotData[,1:2],col=df$S1)
i = 1850
matplot(cbind(covid.goicoa$marginals.random$S1T1[[i]][,1],covid.goicoa.proj$marginals.random$S1T1[[i]][,1]),
        cbind(covid.goicoa$marginals.random$S1T1[[i]][,2],covid.goicoa.proj$marginals.random$S1T1[[i]][,2]),type="l",lty=1,col=1:2)
kl = rep(NA,ns*nt)
for(i in 1:(ns*nt))
{
 # print(i)
 kl[i] = kldivsym(covid.goicoa$marginals.random$S1T1[[i]],covid.goicoa.proj$marginals.random$S1T1[[i]])
}

ggplot(data=plotData)+geom_point(aes(y=StandardModelE,x=NewParametrizationE),colour="red",size=1.25)+xlab("Estimated mean new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+
  ggtitle("Estimated mean standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))
ggsave("Dengue_EstimatedMeanSimulatedDateGoicoa.pdf")
ggplot(data=plotData)+geom_point(aes(y=StandardModelSD,x=NewParametrizationSD),colour="red",size=1.25)+xlab("Estimated sd new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+
  ggtitle("Estimated sd standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))#+
ggsave("Dengue_EstimatedSDSimulatedDataGoicoa.pdf")

plotData=data.table::data.table(StandardModelE=covid.goicoa$summary.random$T1$mean,
                                NewParametrizationE=covid.goicoa.proj$summary.random$T1$mean,
                                StandardModelSD=covid.goicoa$summary.random$T1$sd,
                                NewParametrizationSD=covid.goicoa.proj$summary.random$T1$sd)
ggplot(data=plotData)+geom_point(aes(y=StandardModelE,x=NewParametrizationE),colour="red",size=1.25)+xlab("Estimated mean new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+
  ggtitle("Estimated mean standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))

plotData=data.table::data.table(StandardModelE=covid.goicoa$summary.random$S1$mean,
                                NewParametrizationE=covid.goicoa.proj$summary.random$S1$mean,
                                StandardModelSD=covid.goicoa$summary.random$S1$sd,
                                NewParametrizationSD=covid.goicoa.proj$summary.random$S1$sd)
ggplot(data=plotData)+geom_point(aes(y=StandardModelE,x=NewParametrizationE),colour="red",size=1.25)+xlab("Estimated mean new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+
  ggtitle("Estimated mean standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))



#Comparison of marginal likelihoods
#Need to correct for the standard method


C2 = forceSymmetric(as(Q_st+diag(nrow(Q_st))*eps,"sparseMatrix"))
#A1=kronecker(Diagonal(ns),matrix(rep(1,nt),nrow=1))  #C1
#A2=kronecker(matrix(rep(1,ns),nrow=1),Diagonal(nt))  #A
A = PC$C
A = as(rbind(A2,A1[-1,]),"sparseMatrix")
QA = solve(C2,t(A))
AQA = forceSymmetric(solve(t(QA)%*%t(A)))
Sig.cond = forceSymmetric(solve(C2)-QA%*%AQA%*%t(QA))
val = eigen(Sig.cond)$val
np = diff(dim(A))
logdet = sum(log(val[1:np])) 
covid.goicoa$mlik= covid.goicoa$mlik-0.5*logdet
show(c(covid.goicoa$mlik[1,1],covid.goicoa.proj$mlik[1,1]))


#Knorr-Held:

Projector=create_projector_matrix_interactionIV(ns,nt)

NewConstraints=kronecker(seq(1,nt),diag(ns))
NewConstraintsScaled=kronecker(seq(1,nt),diag(ns))
#NewConstraintsScaled=NewConstraintsScaled-mean(NewConstraintsScaled[which(NewConstraintsScaled[,1]>0),1])
NewConstraintsScaled=apply(NewConstraintsScaled,2,function(x){
  x[which(x>0)]=x[which(x>0)]-mean(x[which(x>0)])
  x=x/sqrt(sum(x^2))
  #x=x/sqrt(sum(x^2))
  return(x)})

#NewConstraintsScaledProjectorMatrix=diag(ns*nt)-(NewConstraintsScaled)%*%t(NewConstraintsScaled)
NewConstraintsScaledProjectorMatrix=-(NewConstraintsScaled)%*%t(NewConstraintsScaled)


Projector2=Projector+NewConstraintsScaledProjectorMatrix


A=rbind(A1,A2[-nrow(A2),],t(NewConstraints))
A=A[-nrow(A),]
D=t(kronecker(diag(nt),rep(1,ns)))
D=D[-1,]

A_extra=cbind(D,0*D)
#A_extra=cbind(0*D,D)

prior.fixed=list(prec.intercept=0.001)

#f(Vw,model="rw1",scale.model=T,constr=TRUE)+f(Vu,model="rw1",scale.model=TRUE,constr=TRUE)
baseformula<- Y ~offset(log(E)) +f(S1,model="generic0",diagonal=eps,Cmatrix = Q_s1,constr=T)+
  f(T2,model="generic0",Cmatrix=Q_RW2+diag(nt)*eps,constr=T) +
  f(S1T1,model="generic0",Cmatrix = Q_st+diag(nrow(Q_st))*eps,constr=F,extraconstr = list(A=A,e=rep(0,nrow(A))))+Vu
#+f(Vw,model="rw1",scale.model=T,constr=TRUE)+f(Vu,model="rw1",scale.model=TRUE,constr=TRUE)+
baseformula.proj <- Y ~ offset(log(E)) +f(S1,model="generic0",diagonal=eps,Cmatrix=Q_s1,constr=T) +
  f(T2,model="generic0",Cmatrix=Q_RW2+diag(nt)*eps,constr=T) +
  f(S1T1,model="z",Z=Projector2, precision=kap,Cmatrix = Q_st+diag(nrow(Q_st))*eps,constr=F,extraconstr = list(A=cbind(D,0*D),e=rep(0,nrow(D))))+Vu

run.knorr = FALSE
if(!file.exists("covid.proj.nb.knorr.RDS") | !file.exists("covid.nb.knorr.RDS"))
  run.knorr = TRUE
if(run.knorr)
{
  covid.proj.nb.knorr <- inla(baseformula.proj, family = "nbinomial",data =df_month,num.threads =4,inla.mode="experimental",control.fixed = list(
    prec.intercept =0.001),verbose=T,control.inla=list(cmin=0 ))
  saveRDS(covid.proj.nb.knorr,"covid.proj.nb.knorr.RDS")
  
  covid.nb.knorr <- inla(baseformula, family = "nbinomial",data=df_month,num.threads = 4,control.fixed = list(
    prec.intercept =0.001),control.inla=list(cmin=0),verbose=T,inla.mode="experimental")
  saveRDS(covid.nb.knorr,"covid.nb.knorr.RDS")
}
if(!run.knorr)
{
  covid.proj.nb.knorr = readRDS("covid.proj.nb.knorr.RDS")
  covid.nb.knorr = readRDS("covid.nb.knorr.RDS")
}

show(c(covid.nb.knorr$cpu.used[4],covid.proj.nb.knorr$cpu.used[4]))

#Table for latex file
tab.knorr = rbind(covid.nb.knorr$summary.fixed[,1:2],
                  covid.nb.knorr$summary.hyperpar[,1:2])
tab.knorr.proj = rbind(covid.proj.nb.knorr$summary.fixed[,1:2],
                       covid.proj.nb.knorr$summary.hyperpar[,1:2])
tab.knorr2 = cbind(tab.knorr,tab.knorr.proj)
rownames(tab.knorr2) = c("mu","Vu","Disp","tau_alpha","tau_theta","tau_delta")
xtable(tab.knorr2,digits=3)

tab = cbind(tab.goicoa2,tab.knorr2)
xtable(tab,digits=3)

#Plotting interaction terms, mean and standard deviations
plotData=data.table::data.table(StandardModelE=covid.nb.knorr$summary.random$S1T1$mean,NewParametrizationE=covid.proj.nb.knorr$summary.random$S1T1$mean[1:(ns*nt)],StandardModelSD=covid.nb.knorr$summary.random$S1T1$sd,NewParametrizationSD=covid.proj.nb.knorr$summary.random$S1T1$sd[1:(ns*nt)])

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

Cmatrix = Q_st+diag(nrow(Q_st))*eps
C2 = forceSymmetric(as(Cmatrix,"sparseMatrix"))
A = as(list(A=A,e=rep(0,nrow(A)))$A,"sparseMatrix")
QA = solve(C2,t(A))
AQA = forceSymmetric(solve(t(QA)%*%t(A)))
Sig.cond = forceSymmetric(solve(C2)-QA%*%AQA%*%t(QA))
val = eigen(Sig.cond)$val
np = diff(dim(A))
logdet = sum(log(val[1:np])) 
covid.nb.knorr$mlik= covid.nb.knorr$mlik-0.5*logdet
show(c(covid.nb.knorr$mlik[1,1],covid.proj.nb.knorr$mlik[1,1]))
