
Tran = rep(NA,ns*nt)
Tran2 = rep(NA,ns*nt)
for(s in 1:ns)
  for(t in 1:nt)
  {
    i = (t-1)*ns+s
    j = (s-1)*nt+t
    Tran[i] = j
    Tran2[j] = i
  }


covid.goicoa.proj=readRDS("covid.goicoa.proj.RDS")
covid.goicoa=readRDS("covid.goicoa.RDS")

covid2.goicoa.proj=readRDS("covid2.goicoa.proj.RDS")
covid2.goicoa=readRDS("covid2.goicoa.RDS")

plotData=data.table::data.table(StandardModelE=covid.goicoa$summary.random$S1T1$mean[1:(ns*nt)],
                                NewParametrizationE=covid2.goicoa$summary.random$S1T1$mean[Tran],
                                StandardModelSD=covid.goicoa$summary.random$S1T1$sd[1:(ns*nt)],
                                NewParametrizationSD=covid2.goicoa$summary.random$S1T1$sd[Tran])

ggplot(data=plotData)+geom_point(aes(y=StandardModelE,x=NewParametrizationE),colour="red",size=1.25)+xlab("Estimated mean new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+
  ggtitle("Estimated mean standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))



plotData=data.table::data.table(StandardModelE=covid.goicoa.proj$summary.random$S1T1$mean[1:(ns*nt)],
                                NewParametrizationE=covid2.goicoa.proj$summary.random$S1T1$mean[Tran],
                                StandardModelSD=covid.goicoa.proj$summary.random$S1T1$sd[1:(ns*nt)],
                                NewParametrizationSD=covid2.goicoa.proj$summary.random$S1T1$sd[Tran])
plotData=plotData[df2$county==1,]
plot(plotData[,1:2],col=df2$county)

ggplot(data=plotData)+geom_point(aes(y=StandardModelE,x=NewParametrizationE),colour="red",size=1.25)+xlab("Estimated mean new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+
  ggtitle("Estimated mean standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))

df2 = df
df2$meandelta = covid.goicoa.proj$summary.random$S1T1$mean[1:(ns*nt)]
df2$mean2delta = covid2.goicoa.proj$summary.random$S1T1$mean[Tran]
df2$diffdelta = covid.goicoa.proj$summary.random$S1T1$mean[1:(ns*nt)]-covid2.goicoa.proj$summary.random$S1T1$mean[Tran]
  
library(mgcv)
fit = gam(diffdelta~s(county)+s(T1),data=df2)

df2$county = as.factor(df2$county)
df2$T1 = as.factor(df2$T1)
fit = lm(diffdelta~county+T1,data=df2)
