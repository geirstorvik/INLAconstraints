##INLA analysis of covid data
library(INLA)
#Can turn on pardiso if available
#INLA::inla.setOption("pardiso.license","~/sys/licences/pardiso.lic")
#INLA::inla.pardiso.check()
library(INLAconstraints)
library(Matrix)
<<<<<<< HEAD
=======
#library(sf)
#library(spdep)
library(data.table)
library(dlnm)
>>>>>>> 077e27a (Update from MI machine)
library(ggplot2)
library(xtable)

rm(list=ls())

#Read data and make indices for temporal and interaction terms
data(coviddata)
df = coviddata$data
df$weekday = as.factor(weekdays(df$date))
df = df[df$date>as.Date("2020-10-01"),]
df$T1 = as.numeric(df$date-df$date[1]+1)
nt = max(df$T1)
df$county = as.numeric(as.factor(df$location_code))
ns = max(df$county)

#Reduce number of time points
nt = 500
df = df[df$T1<(nt+1),]

df$S1T1 = (df$T1-1)*ns+df$county


#Make precision matrices
Q_ICAR = -coviddata$adj
for(i in 1:ns)
  Q_ICAR[i,i] = -sum(Q_ICAR[i,-i])
Q_RW2=GMRF_RW(n=nt,order=1)

eps = 1e-05
kap = 1e06


SCALE = FALSE
if(SCALE)
{
  #Scale model
  A_t=matrix(rep(1,nt),nrow=1)
  Q_RW2 = inla.scale.model(Q_RW2+diag(nt)*eps,list(A=A_t,e=0))
  Q_ICAR=inla.scale.model(Q_ICAR+diag(ns)*eps,constr=list(A=matrix(rep(1,ns),nrow=1),e=0))
}


Q_st=kronecker(Q_RW2,Q_ICAR)

PC = SpaceTimeProjConstr(ns,nt,dim="time",type="SC")

df$E = log(df$pop)
baseformula.proj=cases~offset(E)+#weekday+
  f(T1,model="generic0",Cmatrix=Q_RW2+diag(nt)*eps,constr=T)+
   f(county,model="generic0",Cmatrix=Q_ICAR+diag(ns)*eps,constr = T) +
 f(S1T1,Cmatrix=Q_st+diag(nrow(Q_st))*eps,precision=kap,model="z",Z=as.matrix(PC$P),constr=F,
    extraconstr=list(A=as.matrix(PC$A2),e=rep(0,nrow(PC$A2))))

baseformula =cases~offset(E)+#weekday+
  f(T1,model="generic0",Cmatrix=Q_RW2+diag(nt)*eps,constr=T)+
  f(county,model="generic0",Cmatrix=Q_ICAR+diag(ns)*eps,constr=T) +
  f(S1T1,model="generic0",Cmatrix=Q_st+diag(ns*nt)*eps,constr=F,
     extraconstr=list(A=as.matrix(PC$A),e=rep(0,nrow(PC$A))))



covid.SC.proj=inla(baseformula.proj, family = "poisson",data =df,num.threads ="6:2",inla.mode="experimental",
                         control.fixed = list(
                           prec.intercept =0.01),verbose=T,control.inla=list(strategy="gaussian" ))
#saveRDS(covid.SC.proj,file="covid.SC.proj.RDS")
  
covid.SC=inla(baseformula, family = "poisson",data =df,num.threads ="6:2",inla.mode="experimental",
                    control.fixed = list(
                      prec.intercept =0.01),verbose=T,control.inla=list(strategy="gaussian" ))
#saveRDS(covid.SC,file="covid.SC.RDS")

show(c(covid.SC$cpu.used[4],covid.SC.proj$cpu.used[4]))

#Plotting interaction terms, mean and standard deviations
plotData=data.table::data.table(StandardModelE=covid.SC$summary.random$S1T1$mean,
                                NewParametrizationE=covid.SC.proj$summary.random$S1T1$mean[1:(ns*nt)],
                                StandardModelSD=covid.SC$summary.random$S1T1$sd,
                                NewParametrizationSD=covid.SC.proj$summary.random$S1T1$sd[1:(ns*nt)])
ggplot(data=plotData)+geom_point(aes(y=StandardModelE,x=NewParametrizationE),colour="red",size=1.25)+xlab("Estimated mean new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+
  ggtitle("Estimated mean standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))
ggsave("Covid_EstimatedMeanSimulatedDateSC.pdf",height=5,width=5)
ggplot(data=plotData)+geom_point(aes(y=StandardModelSD,x=NewParametrizationSD),colour="red",size=1.25)+xlab("Estimated sd new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+
  ggtitle("Estimated sd standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))#+
ggsave("Covid_EstimatedSDSimulatedDataSC.pdf",height=5,width=5)

#Comparison of marginal likelihoods
#Need to correct for the standard method
mcorS = LikCorrectGeneric0(Q_ICAR,matrix(rep(1,ns),nrow=1),eps)
mcorT = LikCorrectGeneric0(Q_RW2,matrix(rep(1,nt),nrow=1),eps)
mcorST = LikCorrectGeneric0(Q_st,as.matrix(PC$A),eps)
show(c(covid.SC$mlik[1,1]+mcorS+mcorT+mcorST,covid.SC.proj$mlik[1,1]+mcorS+mcorT)/nrow(df))

#Table for latex file
tab.SC = rbind(covid.SC$summary.fixed[,1:2],
               covid.SC$summary.hyperpar[,1:2])
tab.SC.proj = rbind(covid.SC.proj$summary.fixed[,1:2],
                    covid.SC.proj$summary.hyperpar[,1:2])
tab.SC2 = cbind(tab.SC,tab.SC.proj)
xtable(tab.SC2,digits=3)
