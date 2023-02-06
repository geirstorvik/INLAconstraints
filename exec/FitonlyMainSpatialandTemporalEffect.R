library(INLA)
library(Matrix)
library(INLA)


ns = 10000
nt=10

I2=GMRF_RW(n=ns,order=1)

Q_st_t_ord=kronecker(GMRF_RW(n=nt,order=2),I2)

Projector=create_projector_matrix_interactionIV(ns,nt)

library(data.table)

prior.fixed=list(prec.intercept=0.0001)

KeepData=data.table::data.table(interaction=1:(ns*nt),Y=rnorm(ns*nt,0,3))
KeepData$spatial_group=KeepData$interaction%%ns
KeepData$spatial_group[which(KeepData$spatial_group==0)]=ns


#Inla needs sparse matrices!
I2INLA=INLA::inla.as.sparse(I2+Matrix::Diagonal(x=1*1e-05,ns))
Q_st_t_ordINLA= INLA::inla.as.sparse(Q_st_t_ord+Matrix::Diagonal(x=1*1e-05,ns*nt))

newINLA=inla(Y~f(spatial_group,model="generic0",Cmatrix=I2INLA,constr=T)+f(interaction,model="z",Z=Projector,Cmatrix =Q_st_t_ordINLA,constr=F,precision=1e08),data=KeepData,verbose=T,family="gaussian",control.fixed=prior.fixed,num.threads=1,control.inla=list(int.strategy="eb"))
