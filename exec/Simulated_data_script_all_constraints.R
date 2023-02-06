library(Matrix)
library(SpaceTimePaper)
dataDir <- system.file("extdata", package = "spatioTemporalIndices")
library(INLA)
n=10

graph=system.file("demodata/germany.graph", package="INLA")
I=INLA::inla.graph2matrix(graph)
diag(I)=0
I=-I
diag(I)=-rowSums(I)

I2=I
epsilon=1e-05



library(data.table)
dataDir <- system.file("extdata", package = "SpaceTimePaper")
data.sum.zero=readr::read_rds(paste0(dataDir,"/SpatioTemporalDatasetNew.RDS"))

ns=nrow(I2)
nt=nrow(GMRF_RW(n=10,order=2))

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

A=t(kronecker(rep(1,nt),diag(ns)))
A=rbind(A,t(kronecker(diag(nt),rep(1,ns))))
A=A[-nrow(A),]
A=rbind(A,t(NewConstraints))
A=A[-nrow(A),]
D=t(kronecker(diag(nt),rep(1,ns)))
D=D[-1,]

A_extra=cbind(D,0*D)
#A_extra=cbind(0*D,D)

Q_RW2_before=GMRF_RW(n=nt,order=2)
#Scale model
A_t=rbind(matrix(rep(1,nt),nrow=1),seq(1,nt))
Q_RW2 = inla.scale.model(Q_RW2_before,list(A=A_t,e=rep(0,2)))

#Q_RW2_scaled=inla.scale.model(GMRF_RW(order=2,n=10),constr = list(A=rbind(rep(1,nt),seq(1,nt)),e=rep(0,2)))
Q_ICAR=inla.scale.model(I2,constr=list(A=matrix(rep(1,ns),nrow=1),e=0))
Q_st=kronecker(Q_RW2,Q_ICAR)

prior.fixed=list(prec.intercept=0.001)

resINLA2.full.Rue.Constr.Projector=inla(Y~f(main_temporal,model="generic0",Cmatrix=Q_RW2+diag(10)*1e-05,constr=T)+f(main_spatial,model="generic0",Cmatrix=Q_ICAR+diag(nrow(I2))*1e-05,constr=T)+f(interaction,model="z",precision=1e06,Z=Projector2,Cmatrix = Q_st+diag(ns*nt)*1e-05,constr=F,extraconstr = list(A=A_extra,e=rep(0,nrow(A_extra)))),data=data.sum.zero,verbose=T,family="poisson",control.fixed=prior.fixed,num.threads=4,control.predictor=list(compute=TRUE),inla.mode="experimental")
resINLA2.full.Rue.Constr.ordinary=inla(Y~f(main_temporal,model="generic0",Cmatrix=Q_RW2+diag(10)*1e-05,constr=T)+f(main_spatial,model="generic0",Cmatrix=Q_ICAR+diag(nrow(I2))*1e-05,constr=T)+f(interaction,model="generic0",Cmatrix = Q_st+diag(ns*nt)*1e-05,constr=F,extraconstr = list(A=A,e=rep(0,nrow(A)))),data=data.sum.zero,verbose=T,family="poisson",control.fixed=prior.fixed,num.threads=4,control.predictor=list(compute=TRUE),inla.mode="experimental")


plotData=data.table::data.table(StandardModelE=resINLA2.full.Rue.Constr.ordinary$summary.random$interaction$mean,NewParametrizationE=resINLA2.full.Rue.Constr.Projector$summary.random$interaction$mean[1:(ns*nt)],StandardModelSD=resINLA2.full.Rue.Constr.ordinary$summary.random$interaction$sd,NewParametrizationSD=resINLA2.full.Rue.Constr.Projector$summary.random$interaction$sd[1:(ns*nt)])

library(ggplot2)
ggplot(data=plotData)+geom_point(aes(y=StandardModelE,x=NewParametrizationE),colour="skyblue",size=1.25)+xlab("Estimated mean new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+ggtitle("Estimated mean standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))+
  coord_fixed()+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+labs(x=NULL, y=NULL)
#ggsave("~/Documents/SpatioTemporalModels/EstimatedMeanSimulatedData.pdf")
ggplot(data=plotData)+geom_point(aes(y=StandardModelSD,x=NewParametrizationSD),colour="red",size=1.25)+xlab("Estimated sd new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+ggtitle("Estimated sd standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))+
  coord_fixed()+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
# ggsave("~/Documents/SpatioTemporalModels/EstimatedSDSimulatedData.pdf")


#Goicoa-constraints:


A=t(kronecker(rep(1,nt),diag(ns)))
A=rbind(A,t(kronecker(diag(nt),rep(1,ns))))
A=A[-nrow(A),]

D=t(kronecker(diag(nt),rep(1,ns)))
D=D[-1,]
A_extra=cbind(D,0*D)

resINLA2.full.goicoa.Constr.Projector=inla(Y~f(main_temporal,model="generic0",Cmatrix=Q_RW2+diag(10)*1e-05,constr=T)+f(main_spatial,model="generic0",Cmatrix=Q_ICAR+diag(nrow(I2))*1e-05,constr=T)+f(interaction,model="z",precision=1e06,Z=Projector,Cmatrix = Q_st+diag(ns*nt)*1e-05,constr=F,extraconstr = list(A=A_extra,e=rep(0,nrow(A_extra)))),data=data.sum.zero,verbose=T,family="poisson",control.fixed=prior.fixed,num.threads=4,control.predictor=list(compute=TRUE),inla.mode="experimental")
resINLA2.full.goicoa.Constr.ordinary=inla(Y~f(main_temporal,model="generic0",Cmatrix=Q_RW2+diag(10)*1e-05,constr=T)+f(main_spatial,model="generic0",Cmatrix=Q_ICAR+diag(nrow(I2))*1e-05,constr=T)+f(interaction,model="generic0",Cmatrix = Q_st+diag(ns*nt)*1e-05,constr=F,extraconstr = list(A=A,e=rep(0,nrow(A)))),data=data.sum.zero,verbose=T,family="poisson",control.fixed=prior.fixed,num.threads=4,control.predictor=list(compute=TRUE),inla.mode="experimental")



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


