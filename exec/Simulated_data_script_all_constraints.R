library(Matrix)
library(SpaceTimePaper)
library(INLA)
library(data.table)
library(ggplot2)
library(xtable)
rm(list=ls())
source("../R/SpaceTimeProjConstr.R")

dataDir <- system.file("extdata", package = "spatioTemporalIndices")
dataDir <- system.file("extdata", package = "SpaceTimePaper")
data.sum.zero=readr::read_rds(paste0(dataDir,"/SpatioTemporalDatasetNew.RDS"))
graph=system.file("demodata/germany.graph", package="INLA")

Q_ICAR=INLA::inla.graph2matrix(graph)
diag(Q_ICAR)=0
Q_ICAR=-Q_ICAR
diag(Q_ICAR)=-rowSums(Q_ICAR)
ns=nrow(Q_ICAR)

nt=10
Q_RW2=GMRF_RW(n=nt,order=2)

Projector=create_projector_matrix_interactionIV(ns,nt)

tstar = seq(1,nt)
tstar = tstar-mean(tstar)
tstar = tstar/sqrt(sum(tstar))
NewConstraints=kronecker(seq(1,nt),diag(ns))
NewConstraints=kronecker(tstar,diag(ns))
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

A=t(kronecker(rep(1,nt),diag(ns)))
A=rbind(A,t(kronecker(diag(nt),rep(1,ns))))
A=A[-nrow(A),]
A=rbind(A,t(NewConstraints))
A=A[-nrow(A),]
D=t(kronecker(diag(nt),rep(1,ns)))
D=D[-nrow(D),]

A_extra=cbind(D,0*D)
#A_extra=cbind(0*D,D)

#Scale model
SCALE=FALSE
if(SCALE)
{
 A_t=rbind(matrix(rep(1,nt),nrow=1),seq(1,nt))
 Q_RW2 = inla.scale.model(Q_RW2,list(A=A_t,e=rep(0,2)))
 Q_ICAR=inla.scale.model(Q_ICAR,constr=list(A=matrix(rep(1,ns),nrow=1),e=0))
}
Q_st=kronecker(Q_RW2,Q_ICAR)

prior.fixed=list(prec.intercept=0.001)

eps=1e-05
kap=1e06

#GCF-constraints:

PC = SpaceTimeProjConstr(ns,nt,type ="GCF")

#resINLA2.full.GCF.Constr.Projector=
#  inla(Y~f(main_temporal,model="generic0",Cmatrix=Q_RW2+diag(nt)*1e-05,constr=T)+
#         f(main_spatial,model="generic0",Cmatrix=Q_ICAR+diag(nrow(Q_ICAR))*1e-05,constr=T)+
#         f(interaction,model="z",precision=1e06,Z=Projector2,Cmatrix = Q_st+diag(ns*nt)*1e-05,constr=F,
#           extraconstr = list(A=A_extra,e=rep(0,nrow(A_extra)))),data=data.sum.zero,verbose=T,family="poisson",control.fixed=prior.fixed,num.threads=4,control.predictor=list(compute=TRUE),inla.mode="experimental")

resINLA2.full.GCF.Constr.Projector=
  inla(Y~f(main_temporal,model="generic0",Cmatrix=Q_RW2+diag(nt)*eps,constr=T)+
         f(main_spatial,model="generic0",Cmatrix=Q_ICAR+diag(nrow(Q_ICAR))*eps,constr=T)+
         f(interaction,model="z",precision=kap,Z=PC$P,Cmatrix = Q_st+diag(ns*nt)*eps,constr=F,
           extraconstr = list(A=as.matrix(PC$A2),e=rep(0,nrow(A_extra)))),
       data=data.sum.zero,verbose=T,family="poisson",control.fixed=prior.fixed,num.threads=4,
       control.predictor=list(compute=TRUE),inla.mode="experimental")

#resINLA2.full.GCF.Constr.ordinary=
#  inla(Y~f(main_temporal,model="generic0",Cmatrix=Q_RW2+diag(nt)*1e-05,constr=T)+
#         f(main_spatial,model="generic0",Cmatrix=Q_ICAR+diag(nrow(Q_ICAR))*1e-05,constr=T)+
#         f(interaction,model="generic0",Cmatrix = Q_st+diag(ns*nt)*1e-05,constr=F,extraconstr = list(A=A,e=rep(0,nrow(A)))),data=data.sum.zero,verbose=T,family="poisson",control.fixed=prior.fixed,num.threads=4,control.predictor=list(compute=TRUE),inla.mode="experimental")

resINLA2.full.GCF.Constr.ordinary=
  inla(Y~f(main_temporal,model="generic0",Cmatrix=Q_RW2+diag(nt)*eps,constr=T)+
         f(main_spatial,model="generic0",Cmatrix=Q_ICAR+diag(nrow(Q_ICAR))*eps,constr=T)+
         f(interaction,model="generic0",Cmatrix = Q_st+diag(ns*nt)*eps,constr=F,
           extraconstr = list(A=as.matrix(PC$A),e=rep(0,nrow(PC$A)))),
       data=data.sum.zero,verbose=T,family="poisson",control.fixed=prior.fixed,num.threads=4,control.predictor=list(compute=TRUE),inla.mode="experimental")

show(c(resINLA2.full.GCF.Constr.ordinary$cpu.used[4],resINLA2.full.GCF.Constr.Projector$cpu.used[4]))

plotData=data.table::data.table(StandardModelE=resINLA2.full.GCF.Constr.ordinary$summary.random$interaction$mean,NewParametrizationE=resINLA2.full.GCF.Constr.Projector$summary.random$interaction$mean[1:(ns*nt)],StandardModelSD=resINLA2.full.GCF.Constr.ordinary$summary.random$interaction$sd,NewParametrizationSD=resINLA2.full.GCF.Constr.Projector$summary.random$interaction$sd[1:(ns*nt)])

ggplot(data=plotData)+geom_point(aes(y=StandardModelE,x=NewParametrizationE),colour="skyblue",size=1.25)+xlab("Estimated mean new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+ggtitle("Estimated mean standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))+
  coord_fixed()+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+labs(x=NULL, y=NULL)
#ggsave("~/Documents/SpatioTemporalModels/EstimatedMeanSimulatedData.pdf")
ggplot(data=plotData)+geom_point(aes(y=StandardModelSD,x=NewParametrizationSD),colour="red",size=1.25)+xlab("Estimated sd new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+ggtitle("Estimated sd standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))+
  coord_fixed()+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
# ggsave("~/Documents/SpatioTemporalModels/EstimatedSDSimulatedData.pdf")

C2 = forceSymmetric(as(Q_st+diag(ns*nt)*eps,"sparseMatrix"))
A = PC$A
QA = solve(C2,t(A))
AQA = forceSymmetric(solve(t(QA)%*%t(A)))
Sig.cond = forceSymmetric(solve(C2)-QA%*%AQA%*%t(QA))
val = eigen(Sig.cond)$val
np = diff(dim(A))
logdet.GCF = sum(log(val[1:np]))
show(c(resINLA2.full.GCF.Constr.ordinary$mlik[1,1]-0.5*logdet.GCF,
       resINLA2.full.GCF.Constr.Projector$mlik[1,1])/nrow(data.sum.zero))

tab.gcf = rbind(resINLA2.full.GCF.Constr.ordinary$summary.fixed[,1:2],
                resINLA2.full.GCF.Constr.ordinary$summary.hyperpar[,1:2])
tab.gcf.proj = rbind(resINLA2.full.GCF.Constr.Projector$summary.fixed[,1:2],
                     resINLA2.full.GCF.Constr.Projector$summary.hyperpar[,1:2])
tab.gcf2 = cbind(tab.gcf,tab.gcf.proj)


#SC-constraints:


PC = SpaceTimeProjConstr(ns,nt,type ="SC")

#resINLA2.full.GCF.Constr.Projector=
#  inla(Y~f(main_temporal,model="generic0",Cmatrix=Q_RW2+diag(nt)*1e-05,constr=T)+
#         f(main_spatial,model="generic0",Cmatrix=Q_ICAR+diag(nrow(Q_ICAR))*1e-05,constr=T)+
#         f(interaction,model="z",precision=1e06,Z=Projector2,Cmatrix = Q_st+diag(ns*nt)*1e-05,constr=F,
#           extraconstr = list(A=A_extra,e=rep(0,nrow(A_extra)))),data=data.sum.zero,verbose=T,family="poisson",control.fixed=prior.fixed,num.threads=4,control.predictor=list(compute=TRUE),inla.mode="experimental")

resINLA2.full.SC.Constr.Projector=
  inla(Y~f(main_temporal,model="generic0",Cmatrix=Q_RW2+diag(nt)*eps,constr=T)+
         f(main_spatial,model="generic0",Cmatrix=Q_ICAR+diag(nrow(Q_ICAR))*eps,constr=T)+
         f(interaction,model="z",precision=kap,Z=PC$P,Cmatrix = Q_st+diag(ns*nt)*eps,constr=F,
           extraconstr = list(A=as.matrix(PC$A2),e=rep(0,nrow(A_extra)))),
       data=data.sum.zero,verbose=T,family="poisson",control.fixed=prior.fixed,num.threads=4,
       control.predictor=list(compute=TRUE),inla.mode="experimental")

#resINLA2.full.SC.Constr.ordinary=
#  inla(Y~f(main_temporal,model="generic0",Cmatrix=Q_RW2+diag(nt)*1e-05,constr=T)+
#         f(main_spatial,model="generic0",Cmatrix=Q_ICAR+diag(nrow(Q_ICAR))*1e-05,constr=T)+
#         f(interaction,model="generic0",Cmatrix = Q_st+diag(ns*nt)*1e-05,constr=F,extraconstr = list(A=A,e=rep(0,nrow(A)))),data=data.sum.zero,verbose=T,family="poisson",control.fixed=prior.fixed,num.threads=4,control.predictor=list(compute=TRUE),inla.mode="experimental")

resINLA2.full.SC.Constr.ordinary=
  inla(Y~f(main_temporal,model="generic0",Cmatrix=Q_RW2+diag(nt)*eps,constr=T)+
         f(main_spatial,model="generic0",Cmatrix=Q_ICAR+diag(nrow(Q_ICAR))*eps,constr=T)+
         f(interaction,model="generic0",Cmatrix = Q_st+diag(ns*nt)*eps,constr=F,
           extraconstr = list(A=as.matrix(PC$A),e=rep(0,nrow(PC$A)))),
       data=data.sum.zero,verbose=T,family="poisson",control.fixed=prior.fixed,num.threads=4,control.predictor=list(compute=TRUE),inla.mode="experimental")

show(c(resINLA2.full.SC.Constr.ordinary$cpu.used[4],resINLA2.full.SC.Constr.Projector$cpu.used[4]))

plotData=data.table::data.table(StandardModelE=resINLA2.full.SC.Constr.ordinary$summary.random$interaction$mean,NewParametrizationE=resINLA2.full.SC.Constr.Projector$summary.random$interaction$mean[1:(ns*nt)],StandardModelSD=resINLA2.full.SC.Constr.ordinary$summary.random$interaction$sd,NewParametrizationSD=resINLA2.full.SC.Constr.Projector$summary.random$interaction$sd[1:(ns*nt)])

ggplot(data=plotData)+geom_point(aes(y=StandardModelE,x=NewParametrizationE),colour="skyblue",size=1.25)+xlab("Estimated mean new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+ggtitle("Estimated mean standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))+
  coord_fixed()+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+labs(x=NULL, y=NULL)
#ggsave("~/Documents/SpatioTemporalModels/EstimatedMeanSimulatedData.pdf")
ggplot(data=plotData)+geom_point(aes(y=StandardModelSD,x=NewParametrizationSD),colour="red",size=1.25)+xlab("Estimated sd new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+ggtitle("Estimated sd standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))+
  coord_fixed()+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
# ggsave("~/Documents/SpatioTemporalModels/EstimatedSDSimulatedData.pdf")


C2 = forceSymmetric(as(Q_st+diag(ns*nt)*eps,"sparseMatrix"))
A = PC$A
QA = solve(C2,t(A))
AQA = forceSymmetric(solve(t(QA)%*%t(A)))
Sig.cond = forceSymmetric(solve(C2)-QA%*%AQA%*%t(QA))
val = eigen(Sig.cond)$val
np = diff(dim(A))
logdet.SC = sum(log(val[1:np]))
show(c(resINLA2.full.SC.Constr.ordinary$mlik[1,1]-0.5*logdet.SC,
       resINLA2.full.SC.Constr.Projector$mlik[1,1])/nrow(data.sum.zero))

tab.sc = rbind(resINLA2.full.SC.Constr.ordinary$summary.fixed[,1:2],
                resINLA2.full.SC.Constr.ordinary$summary.hyperpar[,1:2])
tab.sc.proj = rbind(resINLA2.full.SC.Constr.Projector$summary.fixed[,1:2],
                     resINLA2.full.SC.Constr.Projector$summary.hyperpar[,1:2])
tab.sc2 = cbind(tab.sc,tab.sc.proj)

xtable(cbind(tab.gcf2,tab.sc2),3)

if(0)
{
A=t(kronecker(rep(1,nt),diag(ns)))
A=rbind(A,t(kronecker(diag(nt),rep(1,ns))))
A=A[-nrow(A),]

D=t(kronecker(diag(nt),rep(1,ns)))
D=D[-1,]
A_extra=cbind(D,0*D)

resINLA2.full.goicoa.Constr.Projector=inla(Y~f(main_temporal,model="generic0",Cmatrix=Q_RW2+diag(nt)*1e-05,constr=T)+f(main_spatial,model="generic0",Cmatrix=Q_ICAR+diag(nrow(Q_ICAR))*1e-05,constr=T)+f(interaction,model="z",precision=1e06,Z=Projector,Cmatrix = Q_st+diag(ns*nt)*1e-05,constr=F,extraconstr = list(A=A_extra,e=rep(0,nrow(A_extra)))),data=data.sum.zero,verbose=T,family="poisson",control.fixed=prior.fixed,num.threads=4,control.predictor=list(compute=TRUE),inla.mode="experimental")
resINLA2.full.goicoa.Constr.ordinary=inla(Y~f(main_temporal,model="generic0",Cmatrix=Q_RW2+diag(nt)*1e-05,constr=T)+f(main_spatial,model="generic0",Cmatrix=Q_ICAR+diag(nrow(Q_ICAR))*1e-05,constr=T)+f(interaction,model="generic0",Cmatrix = Q_st+diag(ns*nt)*1e-05,constr=F,extraconstr = list(A=A,e=rep(0,nrow(A)))),data=data.sum.zero,verbose=T,family="poisson",control.fixed=prior.fixed,num.threads=4,control.predictor=list(compute=TRUE),inla.mode="experimental")



plotData=data.table::data.table(StandardModelE=resINLA2.full.goicoa.Constr.ordinary$summary.random$interaction$mean,NewParametrizationE=resINLA2.full.goicoa.Constr.Projector$summary.random$interaction$mean[1:(ns*nt)],StandardModelSD=resINLA2.full.goicoa.Constr.ordinary$summary.random$interaction$sd,NewParametrizationSD=resINLA2.full.goicoa.Constr.Projector$summary.random$interaction$sd[1:(ns*nt)])

library(ggplot2)
ggplot(data=plotData)+geom_point(aes(y=StandardModelE,x=NewParametrizationE),colour="skyblue",size=1.25)+xlab("Estimated mean new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+ggtitle("Estimated mean standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))+
  coord_fixed()+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+labs(x=NULL, y=NULL)
#ggsave("~/Documents/SpatioTemporalModels/EstimatedMeanSimulatedData.pdf")
ggplot(data=plotData)+geom_point(aes(y=StandardModelSD,x=NewParametrizationSD),colour="red",size=1.25)+xlab("Estimated sd new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+ggtitle("Estimated sd standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))+
  coord_fixed()+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
# ggsave("~/Documents/SpatioTemporalModels/EstimatedSDSimulatedData.pdf")

}
