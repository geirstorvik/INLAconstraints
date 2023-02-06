#INLA::"/nr/samba/user/fredrik"
INLA::inla.setOption("pardiso.license","~/Documents/licences/pardiso.lic")
INLA::inla.pardiso.check()
library(Matrix)
library(SpaceTimePaper)
dataDir <- system.file("extdata", package = "spatioTemporalIndices")
library(INLA)
n=10
Q_RW2_before=GMRF_RW(n=n,order=2)
#Q_RW2 = inla.scale.model(Q_RW2_before,list(A=A_t,e=0))


graph=system.file("demodata/germany.graph", package="INLA")
I=INLA::inla.graph2matrix(graph)
diag(I)=0
I=-I
diag(I)=-rowSums(I)

I2=I
epsilon=1e-05



library(data.table)
dataDir <- system.file("extdata", package = "SpaceTimePaper")
data.sum.zero=readr::read_rds(paste0(dataDir,"/SimulatedData.RDS"))
#data=data[!duplicated_rows,]
Q_st_t_ord=kronecker(GMRF_RW(n=10,order=2),I2)

ns=nrow(I2)
nt=nrow(GMRF_RW(n=10,order=2))

B=matrix(0,ns*nt,ns*nt)
for(i in 1:nrow(I2)){
  out=rep(0,ns)
  out[i]=1
  tot=rep(out,nt)
  tot=tot/sqrt(sum(tot^2))
 # print(tot[1:ns])
  B=B+(tot)%*%t(tot)
}
C=as(B,"dgTMatrix")

Projector=diag(nrow(C))-C

A=t(kronecker(rep(1,nt),diag(ns)))
A=rbind(A,t(kronecker(diag(nt),rep(1,ns))))
A=A[-1,]
a.qr=qr(t(A))
new_vectors=NULL
new_vectors$Q=qr.Q(a.qr,complete=F)
A_extra_new=new_vectors$Q[,(nrow(I2)+1):ncol(new_vectors$Q)]
A_extra=cbind(t(A_extra_new),0*t(A_extra_new))





prior.fixed=list(prec.intercept=0.0001)





resINLASum2.full.2=inla(Y~f(time_group,model="generic0",Cmatrix=GMRF_RW(order=2,n=10)+diag(10)*1e-05,constr=T)+f(spatial_group,model="generic0",Cmatrix=I2+diag(nrow(I2))*1e-05,constr=T)+f(interaction,model="generic0",Cmatrix = Q_st_t_ord+diag(nrow(Q_st_t_ord))*1e-05,constr=F,extraconstr = list(A=as.matrix(A),e=rep(0,nrow((A))))),data=data.sum.zero[1:150000,],verbose=T,family="poisson",control.fixed=prior.fixed,num.threads=4)
resINLA2.full=inla(Y~f(time_group,model="generic0",Cmatrix=GMRF_RW(order=2,n=10)+diag(10)*1e-05,constr=T)+f(spatial_group,model="generic0",Cmatrix=I2+diag(nrow(I2))*1e-05,constr=T)+f(interaction,model="z",precision=1000,Z=Projector,Cmatrix = Q_st_t_ord+diag(nrow(Q_st_t_ord))*1e-05,constr=F,extraconstr = list(A=as.matrix(A_extra),e=rep(0,nrow(A_extra)))),data=data.sum.zero[1:150000,],verbose=T,family="poisson",control.fixed=prior.fixed,num.threads=4)
#resINLASum2.full.alt2.normk=inla(Y~f(time_group,model="generic0",Cmatrix=GMRF_RW(order=2,n=10)+diag(10)*1e-05,constr=T)+f(spatial_group,model="generic0",Cmatrix=I2+diag(nrow(I2))*1e-05,constr=T)+f(interaction,model="generic0",Cmatrix = Q_st_t_ord+diag(nrow(Q_st_t_ord))*1e-05,constr=F,extraconstr = list(A=as.matrix(A.norm),e=rep(0,nrow((A_sum_constr))))),data=data.sum.zero[1:no,],verbose=T,family="poisson",control.fixed=prior.fixed,num.threads=4)p
#saveRDS(resINLASum2.full.2,"~/Documents/SpatioTemporalModels/SimulatedDataINLAResSumZero.RDS")
#saveRDS(resINLA2.full,"~/Documents/SpatioTemporalModels/SimulatedDataINLAResNewPar.RDS")

plotData=data.table::data.table(StandardModelE=resINLASum2.full.2$summary.random$interaction$mean,NewParametrizationE=resINLA2.full$summary.random$interaction$mean[1:5440],StandardModelSD=resINLASum2.full.2$summary.random$interaction$sd,NewParametrizationSD=resINLA2.full$summary.random$interaction$sd[1:5440])
library(ggplot2)
ggplot(data=plotData)+geom_point(aes(y=StandardModelE,x=NewParametrizationE),colour="skyblue",size=1.25)+xlab("Estimated mean new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+ggtitle("Estimated mean standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))+
  coord_fixed()+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+labs(x=NULL, y=NULL)
#ggsave("~/Documents/SpatioTemporalModels/EstimatedMeanSimulatedData.pdf")
ggplot(data=plotData)+geom_point(aes(y=StandardModelSD,x=NewParametrizationSD),colour="red",size=1.25)+xlab("Estimated sd new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+ggtitle("Estimated sd standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))+
  coord_fixed()+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
# ggsave("~/Documents/SpatioTemporalModels/EstimatedSDSimulatedData.pdf")
