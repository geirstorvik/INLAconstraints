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
library(ggplot2)
library(xtable)

#rm(list=ls())

#Read data and make indices for temporal and interaction terms
d = readRDS("../data/coviddata.RDS")
df = d$data
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
Qs = -d$adj
for(i in 1:ns)
  Qs[i,i] = -sum(Qs[i,-i])
Q_RW2_before=GMRF_RW(n=nt,order=2)

eps = 1e-05
kap = 1e06


#Scale model
A_t=rbind(matrix(rep(1,nt)/sqrt(nt),nrow=1),seq(1,nt)/sqrt(sum(seq(1,nt)^2)))
A_t=rbind(matrix(rep(1,nt),nrow=1),seq(1,nt))
Q_RW2 = inla.scale.model(Q_RW2_before,list(A=A_t,e=rep(0,2)))
A_t=matrix(rep(1,nt),nrow=1)
Q_RW2 = inla.scale.model(Q_RW2_before+diag(nt)*eps,list(A=A_t,e=0))
Q_RW2 = Q_RW2_before

Q_ICAR=inla.scale.model(Qs+diag(ns)*eps,constr=list(A=matrix(rep(1,ns),nrow=1),e=0))
Q_ICAR = Qs


Q_st=kronecker(Q_RW2,Q_ICAR)
#Q_st = inla.scale.model(Q_st,list(A=D.g,e=rep(0,nrow(D.g))))
source("../R/SpaceTimeProjConstr.R")
PC = SpaceTimeProjConstr(ns,nt,dim="time")

df$T2 = df$T1
df$county2 = df$county

scT = Q_RW2[1,1]/Q_RW2_before[1,1]
scS = Q_ICAR[1,1]/Qs[1,1]
scS = scT = 1

df$E = log(df$pop)
baseformula.proj=cases~offset(E)+#weekday+
  f(T1,model="generic0",Cmatrix=Q_RW2+diag(nt)*eps*scT,constr=T)+
  #f(T2,model="iid",constr=TRUE)+
  f(county,model="generic0",Cmatrix=Q_ICAR+diag(ns)*eps*scS,constr = T) +
  #f(county2,model="iid",constr=TRUE)+
  f(S1T1,Cmatrix=Q_st+diag(nrow(Q_st))*eps*(scT*scS),precision=kap,model="z",Z=as.matrix(PC$P),constr=F,
    extraconstr=list(A=as.matrix(PC$A2),e=rep(0,nrow(PC$A2))))

baseformula =cases~offset(E)+#weekday+
  f(T1,model="generic0",Cmatrix=Q_RW2+diag(nt)*eps*scT,constr=T)+
  #f(T2,model="iid",constr=TRUE)+
  f(county,model="generic0",Cmatrix=Q_ICAR+diag(ns)*eps*scS,constr=T) +
  #f(county2,model="iid",constr=TRUE)+
  f(S1T1,model="generic0",Cmatrix=Q_st+diag(ns*nt)*eps*(scT*scS),constr=F,
     extraconstr=list(A=as.matrix(PC$A),e=rep(0,nrow(PC$A))))



covid.goicoa.proj=inla(baseformula.proj, family = "poisson",data =df,num.threads =10,#inla.mode="experimental",
                         control.fixed = list(
                           prec.intercept =0.01),verbose=T,control.inla=list(strategy="gaussian" ))
#saveRDS(covid.goicoa.proj,file="covid.goicoa.proj.RDS")
  
covid.goicoa=inla(baseformula, family = "poisson",data =df,num.threads =10,#inla.mode="experimental",
                    control.fixed = list(
                      prec.intercept =0.01),verbose=T,control.inla=list(strategy="gaussian" ))
#saveRDS(covid.goicoa,file="covid.goicoa.RDS")

show(c(covid.goicoa$cpu.used[4],covid.goicoa.proj$cpu.used[4]))

#Table for latex file
tab.goicoa = rbind(covid.goicoa$summary.fixed[,1:2],
                   covid.goicoa$summary.hyperpar[,1:2])
tab.goicoa.proj = rbind(covid.goicoa.proj$summary.fixed[,1:2],
                        covid.goicoa.proj$summary.hyperpar[,1:2])
tab.goicoa2 = cbind(tab.goicoa,tab.goicoa.proj)
#rownames(tab.goicoa2) = c("mu","tau_alpha","tau_theta","tau_delta")
xtable(tab.goicoa2,digits=3)

#Plotting interaction terms, mean and standard deviations
plotData=data.table::data.table(StandardModelE=covid.goicoa$summary.random$S1T1$mean,
                                NewParametrizationE=covid.goicoa.proj$summary.random$S1T1$mean[1:(ns*nt)],
                                StandardModelSD=covid.goicoa$summary.random$S1T1$sd,
                                NewParametrizationSD=covid.goicoa.proj$summary.random$S1T1$sd[1:(ns*nt)])
ggplot(data=plotData)+geom_point(aes(y=StandardModelE,x=NewParametrizationE),colour="red",size=1.25)+xlab("Estimated mean new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+
  ggtitle("Estimated mean standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))
ggsave("Covid_EstimatedMeanSimulatedDateGoicoa.pdf",height=5,width=5)
ggplot(data=plotData)+geom_point(aes(y=StandardModelSD,x=NewParametrizationSD),colour="red",size=1.25)+xlab("Estimated sd new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+
  ggtitle("Estimated sd standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))#+
ggsave("Covid_EstimatedSDSimulatedDataGoicoa.pdf",height=5,width=5)

plotData=data.table::data.table(StandardModelE=covid.goicoa$summary.random$T1$mean,
                                NewParametrizationE=covid.goicoa.proj$summary.random$T1$mean,
                                StandardModelSD=covid.goicoa$summary.random$T1$sd,
                                NewParametrizationSD=covid.goicoa.proj$summary.random$T1$sd)
ggplot(data=plotData)+geom_point(aes(y=StandardModelE,x=NewParametrizationE),colour="red",size=1.25)+xlab("Estimated mean new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+
  ggtitle("Estimated mean standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))

plotData=data.table::data.table(StandardModelE=covid.goicoa$summary.random$county$mean,
                                NewParametrizationE=covid.goicoa.proj$summary.random$county$mean,
                                StandardModelSD=covid.goicoa$summary.random$county$sd,
                                NewParametrizationSD=covid.goicoa.proj$summary.random$county$sd)
ggplot(data=plotData)+geom_point(aes(y=StandardModelE,x=NewParametrizationE),colour="red",size=1.25)+xlab("Estimated mean new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+
  ggtitle("Estimated mean standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))



#Comparison of marginal likelihoods
#Need to correct for the standard method
source('../R/LikCorrectGeneric0.R')
mcorS = LikCorrectGeneric0(Q_ICAR,matrix(rep(1,ns),nrow=1),eps*scS)
mcorT = LikCorrectGeneric0(Q_RW2,matrix(rep(1,nt),nrow=1),eps*scT)
mcorST = LikCorrectGeneric0(Q_st,PC$A,eps*scT*scS)
mcorST = LikCorrectGeneric0(Q_st,PC$A,eps)
show(c(covid.goicoa$mlik[1,1]+mcorS+mcorT+mcorST,covid.goicoa.proj$mlik[1,1]+mcorS+mcorT)/nrow(df))

dates = df$date[df$location_code=="county03"]
xx = data.frame(date=dates,value=c(covid.goicoa.proj$summary.random$T1$mean,covid.goicoa.proj$summary.random$T2$mean),
                variable=c(rep(1,nt),rep(2,nt)))
xx$variable = as.character(xx$variable)
ggplot(xx, aes(x = date, y = value)) + 
  geom_line(aes(color = variable), size = 1) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal()
plot.ts(cbind(covid.goicoa.proj$summary.random$T1$mean,covid.goicoa.proj$summary.random$T2$mean),plot.type="single",col=1:2)
