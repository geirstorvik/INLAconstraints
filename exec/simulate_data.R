library(Matrix)
library(INLA)
library(sparseMVN)
library(SpaceTimePaper)
library(data.table)
#Number of time points



set.seed(20)
nt=10
graph=system.file("demodata/germany.graph", package="INLA")
I=INLA::inla.graph2matrix(graph)
diag(I)=0
I=-I
diag(I)=-rowSums(I)

I2=I
ns=nrow(I2)
#Small value, so that we get proper precision matrix
#We first sample from a proper distribution with too large variances in some directions, and later correct these variances.
epsilon=1e-05
#Set up constraints
A_t = matrix(1,1,nt)
A_t=rbind(A_t,seq(1,nt))
A_s=matrix(1,ncol=ns,nrow=1)
A_st=cbind(kronecker(rep(1,nt),diag(ns)),kronecker(diag(nt),rep(1,ns))[,-1],kronecker(seq(1,nt),diag(ns))[,-1])
A_st=t(A_st)
#Set up precision matrix
Q_RW2_before=GMRF_RW(n=nt,order=2)
#Scale model
Q_RW2 = inla.scale.model(Q_RW2_before,list(A=A_t,e=rep(0,2)))

#Q_RW2_scaled=inla.scale.model(GMRF_RW(order=2,n=10),constr = list(A=rbind(rep(1,nt),seq(1,nt)),e=rep(0,2)))
Q_ICAR=inla.scale.model(I2,constr=list(A=matrix(rep(1,ns),nrow=1),e=0))
Q_st=kronecker(Q_RW2,Q_ICAR)

#Use INLA to create samples
set.seed(20)
x_s_main= 0.25*INLA::inla.qsample(Q=Q_ICAR+diag(ns)*1e-05,n=1,mu=rep(0,ns),constr=list(A=matrix(rep(1,ns),nrow=1),e=0))
x_t_main=0.2*INLA::inla.qsample(Q=Q_RW2+diag(nt)*1e-05,n=1,mu=rep(0,nt),constr=list(A=A_t,e=rep(0,2)))
x_st_main=0.25*INLA::inla.qsample(Q=Q_st+diag(nt*ns)*1e-05,n=1,mu=rep(0,nt*ns),constr=list(A=A_st,e=rep(0,nrow(A_st))))

# #NB: we must add the slope to get a simple sum to zero constraint.
# slope = rnorm(1,mean=0,sd=1)
# slopes_st_spatial=rnorm(2*ns,mean=0,sd=1)
# slopes_st_temp=rnorm(nt,mean=0,sd=1)
main=seq(1,ns*nt)%%ns
main[main==0]=ns
temporal=ceiling(seq(1,ns*nt)/ns)

data=data.table(main_spatial=main,main_temporal=temporal,interaction=seq(1,ns*nt))
dataF=copy(data)
dataF$main_spatial=factor(dataF$main_spatial)
dataF$main_temporal=factor(dataF$main_temporal)
dataF$interaction=factor(dataF$interaction)

DesignMatrixMainSpatial=model.matrix(~main_spatial-1,dataF)
DesignMatrixMainTemporal=model.matrix(~main_temporal-1,dataF)
DesignMatrixMainInteraction=model.matrix(~interaction-1,dataF)



temporal_C=t(kronecker(diag(nt),rep(1,ns))/sqrt(ns))
spatial_C=t(kronecker(rep(1,nt),diag(ns))/sqrt(nt))
#V=(seq(1,nt)-mean(seq(1,nt)))/norm((seq(1,nt)-mean(seq(1,nt))),"2")
#spatial_trend_C=t(kronecker(V,diag(ns)))
temporal_main_C=rep(seq(1,nt)/norm(seq(1,nt),"2"),ns)
mean_pred=DesignMatrixMainTemporal%*%as.vector(x_t_main)+DesignMatrixMainSpatial%*%as.vector(x_s_main)+DesignMatrixMainInteraction%*%as.vector(x_st_main)
mean_pred_model=mean_pred+1.5

combined_pred=do.call(rbind,replicate(30,mean_pred_model,simplify=FALSE))
data=do.call(rbind,replicate(30,data,simplify = FALSE))

set.seed(20)
mu = 1
data$Y=rpois(n=nrow(combined_pred),lambda=exp(combined_pred))
show(range(data$Y))
SpatioTemporalDataset = data
save(SpatioTemporalDataset,file="data/SpatioTemporalDataset.rda")

#saveRDS(data,file="data/SpatioTemporalDataset.RDS")
m= nt

#min_limit_y=min(c(x_t_main[,1]+slope*seq(1,m)/m,x_t_main[,1]))
#max_limit_y=max(c(x_t_main[,1]+slope*seq(1,m)/m,x_t_main[,1]))


#Create plot, using base graphics. Black: simple sum to zero, while blue is "fully" constrained field.
#plot(x_t_main+slope*seq(1,m)/m,ylim=c(-1+min_limit_y,max_limit_y+1),type="l")
#points(x_t_main,col="blue",type="l")




mean_pred=mean_pred-0.75*temporal_main_C
for(i in 1:(nt)){
  print(i)
  mean_pred=mean_pred-0.5*temporal_C[i,]
}

for(i in 1:(ns)){
  print(i)
  mean_pred=mean_pred+.25*spatial_C[i,]
  # mean_pred=mean_pred+slopes_st_spatial[i+ns]*spatial_trend_C[i,]

}

