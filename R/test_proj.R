library(graph4lg)
ns = 11
nt = 200
PC = SpaceTimeProjConstr(ns,nt,type="time")

PC2 = SpaceTimeProjConstr(nt,ns,type="space")

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

PC2$Ptran = PC2$P[Tran,Tran]
show(range(PC$P-PC2$Ptran))
PC2$Ctran = PC2$C[,Tran]
show(range(PC$C-PC2$C))
PC2$C1tran = PC2$C1[,c(Tran,ns*nt+Tran)]
show(range(PC$C1-PC2$C1))

covid.goicoa.proj=readRDS("covid.goicoa.proj.RDS")
covid.goicoa=readRDS("covid.goicoa.RDS")

covid2.goicoa.proj=readRDS("covid2.goicoa.proj.RDS")
covid2.goicoa=readRDS("covid2.goicoa.RDS")

plotData=data.table::data.table(StandardModelE=covid.goicoa$summary.random$S1T1$mean[Tran],
                                NewParametrizationE=covid2.goicoa$summary.random$S1T1$mean[1:(ns*nt)],
                                StandardModelSD=covid.goicoa$summary.random$S1T1$sd[Tran],
                                NewParametrizationSD=covid2.goicoa$summary.random$S1T1$sd[1:(ns*nt)])

ggplot(data=plotData)+geom_point(aes(y=StandardModelE,x=NewParametrizationE),colour="red",size=1.25)+xlab("Estimated mean new parametrization")+
  ylab("Estimated mean standard parametrization")+geom_abline(intercept=0,slope=1,size=0.5)+
  ggtitle("Estimated mean standard vs new parametrization")+ theme(plot.title=element_text(hjust=0.5))
