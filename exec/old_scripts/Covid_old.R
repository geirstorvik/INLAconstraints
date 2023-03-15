library(data.table)
d = data.table(read.csv("data/municipality.csv"))
d2 = d[,sum(cases),by=list(date,fylke_name)]
d3 = d2[fylke_name=="Oslo"]
plot.ts(d3$V1)
d3$cases = c(NA,diff(d3$V1))
plot.ts(d3$cases)
d4 = d[,sum(cases),by=list(date)]
d4$cases = c(NA,diff(d4$V1))


a = data.table(read.csv("data/admissions.csv"))



county = list(Agder = c("Sørlandet sykehus HF"),
Innlandet = c("Sykehuset Innlandet HF"),
More_og_Romsdal = c("Helse Møre og Romsdal"),
Nordland = c("Helgelandssykehuset HF","Nordlandssykehuset HF"),
Oslo = c("Oslo universitetssykehus HF",
         "Private med avtale med Helse Sør-Øst i Oslo"),
Rogaland = c("Helse Stavanger HF"),
Vestfold_og_Telemark = c("Sykehuset i Vestfold HF","Sykehuset Telemark HF"),
Troms_og_Finnmark=c("Finnmarkssykehuset HF","Universitetssykehuset Nord-Norge HF"),
Trondelag = c("Helse Nord-Trøndelag","St. Olavs hospital "),
Vestland = c("Helse Bergen HF + Privat med avtale med Helse Vest i Bergen",
             "Helse Fonna HF","Helse Førde HF"),
Viken = c("Akershus universitetssykehus HF","Sunnaas Sykehus HF",
          "Sykehuset Østfold HF","Vestre Viken HF"))

a$county = NA
for(f in 1:11)
{
  for(i in 1:length(county[[f]]))
    a$county[a$health_org_name==county[[f]][i]] = f
}

acounty = a[,.(cases=sum(admissions)),by=list(date,county)]

a$county = NA
for(f in 1:11)
{
 for(i in 1:length(county[[f]]))
  a$county[a$health_org_name==county[[f]][i]] = f
}


##Making adjency matrix
Adj = matrix(0,nrow=11,ncol=11)
Adj[1,c(6,7)] = 1
Adj[2,c(11,10,3,9)] = 1
Adj[3,c(10,2,9)] = 1
Adj[4,c(9,8)] = 1
Adj[5,c(11)] = 1
Adj[6,c(1,10)] = 1
Adj[7,c(1,10,11)] = 1
Adj[8,c(4)] = 1
Adj[9,c(2,3,4)] = 1
Adj[10,c(6,7,11,2,3)] = 1
Adj[11,c(7,10,2,5)] = 1

acounty$date = as.Date(acounty$date)

coviddata = list(data=acounty,adj=Adj)
saveRDS(coviddata,file="data/coviddata.RDS")





