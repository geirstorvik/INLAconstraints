
# Script taken from Rachel Lowe (2021) and altered by Aanes and Storvik. Only parts of the scripts are used.
# https://github.com/drrachellowe/hydromet_dengue

# Step 0: load packages and pre-processed data
# Step 1: formulate a baseline model on a subset of the data, including spatiotemporal random effects
# Step 2: fit model on smaller dataset, 8 years of data using data for month 6 only

#NB: We remove one region to get a fully connected graph.

# Step 0: load packages and pre-processed data
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
library(SpaceTimePaper)
library(INLA)
library(Matrix)
nt=6

objects_to_be_used=dengue_data_set_up(nt=nt)



#We fit a simplified model first on a smaller subset of the data, using six years of data.
df_check=objects_to_be_used$data_check
df= objects_to_be_used$data_check
df_month=objects_to_be_used$df_month
Q_s1=objects_to_be_used$Q_s1
graph=objects_to_be_used$graph
ns=nrow(Q_s1)
#Create projector matrix for the projection method.
Projector=create_projector_matrix_interactionIV(ns=nrow(Q_s1),nt)

#Create precision matrix for the space time term for the first six years.

#Set up constraints for space time term
A1=t(kronecker(diag(1,nt),rep(1,nrow(Q_s1))))
A2=t(kronecker(rep(1,nt),diag(nrow(Q_s1))))
D.g=rbind(A1,A2[-nrow(A2),])

Q_RW2_before=GMRF_RW(n=nt,order=2)
#Scale model
A_t=rbind(matrix(rep(1,nt),nrow=1),seq(1,nt))
Q_RW2 = inla.scale.model(Q_RW2_before,list(A=A_t,e=rep(0,2)))

#Q_RW2_scaled=inla.scale.model(GMRF_RW(order=2,n=10),constr = list(A=rbind(rep(1,nt),seq(1,nt)),e=rep(0,2)))
Q_ICAR=inla.scale.model(Q_s1,constr=list(A=matrix(rep(1,ns),nrow=1),e=0))
Q_st=kronecker(Q_RW2,Q_ICAR)


baseformula = 
  Y ~ f(T2, model = "generic0", Cmatrix = Q_RW2 + diag(nt) * 1e-05, constr = T) + 
  f(S1, model = "generic0", Cmatrix = Q_ICAR +diag(ns) * 1e-05, constr = T) + f(S1T1, model = "generic0",
                                                                                 Cmatrix = Q_st + diag(nrow(Q_st)) * 1e-05, constr = F, extraconstr = list(A = rbind(A1,
                                                                                                                                                                     A2[-1, ]), e = rep(0, -1 + nrow(A2) + nrow(A1)))) + offset(log(E))
baseformula.proj=Y ~f(T2, model = "generic0", Cmatrix = Q_RW2 + diag(nt) *
                                        1e-05, constr = T)+ f(S1, model = "generic0", Cmatrix = Q_ICAR +
                        diag(ns) * 1e-05, constr = T) + f(S1T1,
                                                          Cmatrix = Q_st + diag(nrow(Q_st)) * 1e-05,precision=1e06, model="z",Z=Projector,extraconstr=list(A=cbind(A1[-1,],0*A1[-1,]),e=rep(0,nt-1))) +offset(log(E))
model.goicoa.proj=inla(baseformula.proj, family = "poisson",data =df_month,num.threads =4,inla.mode="experimental",control.fixed = list(
  prec.intercept =0.01),verbose=T,control.inla=list(cmin=0,strategy="gaussian" ))

model.goicoa=inla(baseformula, family = "poisson",data =df_month,num.threads =4,inla.mode="experimental",control.fixed = list(
  prec.intercept =0.01),verbose=T,control.inla=list(cmin=0,strategy="gaussian" ))


plotData=data.table::data.table(StandardModelE=model.goicoia$summary.random$S1T1$mean,NewParametrizationE=model.goicoa.proj$summary.random$S1T1$mean[1:(ns*nt)],StandardModelSD=model.goicoa$summary.random$S1T1$sd,NewParametrizationSD=model.goicoa.proj$summary.random$S1T1$sd[1:(ns*nt)])

library(ggplot2)
ggplot(data=plotData)+geom_point(aes(y=StandardModelE,x=NewParametrizationE),colour="skyblue",size=1.25)+xlab("Estimated mean new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+ggtitle("Estimated mean standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))
  #ggsave("~/EstimatedMeanSimulatedDatedGoicoa2.pdf")
ggplot(data=plotData)+geom_point(aes(y=StandardModelSD,x=NewParametrizationSD),colour="red",size=1.25)+xlab("Estimated sd new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+ggtitle("Estimated sd standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))+
  coord_fixed()+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
 # ggsave("~/EstimatedSDSimulatedData.pdf")




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
baseformula<- Y ~offset(log(E)) +f(S1,model="generic0",diagonal=1e-05,Cmatrix = Q_ICAR,constr=T)+f(T2,model="generic0",Cmatrix=Q_RW2+diag(nt)*1e-05,constr=T) +
  f(S1T1,model="generic0",Cmatrix = Q_st+diag(nrow(Q_st))*1e-05,constr=F,extraconstr = list(A=A,e=rep(0,nrow(A))))
#+f(Vw,model="rw1",scale.model=T,constr=TRUE)+f(Vu,model="rw1",scale.model=TRUE,constr=TRUE)+
baseformula.proj <- Y ~ offset(log(E)) +f(S1,model="generic0",diagonal=1e-05,Cmatrix=Q_ICAR,constr=T) +f(T2,model="generic0",Cmatrix=Q_RW2+diag(nt)*1e-05,constr=T) +
  f(S1T1,model="z",Z=Projector2, precision=1e6,Cmatrix = Q_st+diag(nrow(Q_st))*1e-05,constr=F,extraconstr = list(A=cbind(D,0*D),e=rep(0,nrow(D))))


model.proj.nb.knorr <- inla(baseformula.proj, family = "poisson",data =df_month,num.threads =4,inla.mode="experimental",control.fixed = list(
  prec.intercept =0.001),verbose=T,control.inla=list(cmin=0,strategy="gaussian" ))
model.nb.knorr <- inla(baseformula, family = "poisson",data=df_month,num.threads = 4,control.fixed = list(
  prec.intercept =0.001),control.inla=list(cmin=0,strategy="gaussian" ),verbose=T,inla.mode="experimental")


plotData=data.table::data.table(StandardModelE=model.nb.knorr$summary.random$S1T1$mean,NewParametrizationE=model.proj.nb.knorr$summary.random$S1T1$mean[1:(ns*nt)],StandardModelSD=model.nb.knorr$summary.random$S1T1$sd,NewParametrizationSD=model.proj.nb.knorr$summary.random$S1T1$sd[1:(ns*nt)])


library(ggplot2)
ggplot(data=plotData)+geom_point(aes(y=StandardModelE,x=NewParametrizationE),colour="skyblue",size=1.25)+xlab("Estimated mean new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+ggtitle("Estimated mean standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))+
  coord_fixed()+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+labs(x=NULL, y=NULL)
#ggsave("~/Documents/SpatioTemporalModels/EstimatedMeanSimulatedData.pdf")
ggplot(data=plotData)+geom_point(aes(y=StandardModelSD,x=NewParametrizationSD),colour="red",size=1.25)+xlab("Estimated sd new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+ggtitle("Estimated sd standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))+
  coord_fixed()+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
# ggsave("~/Documents/SpatioTemporalModels/EstimatedSDSimulatedData.pdf")

