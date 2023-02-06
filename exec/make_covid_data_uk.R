##Covid data from UK
##Link:  https://api.coronavirus.data.gov.uk/v2/data?areaType=region&metric=cumDailyNsoDeathsByDeathDate&format=csv
d = read.csv("data/region_2022-12-22.csv")
d$date = as.Date(d$date)
cname = unique(d$areaName)
ns = length(cname)
for(i in 1:ns)
  d$S1[d$areaName==cname[i]] = i
nt = as.numeric(max(d$date)-min(d$date)+1)
d$T1 = as.numeric(d$date-min(d$date))+1
d$S1T1 = (d$T1-1)*ns+d$S1
d = d[,-c(1,3)]
if(1)
{
  foo = pmatch(seq(1,ns*nt),d$S1T1)
  foo2 = data.frame(S1T1=seq(1,ns*nt)[is.na(foo)])
  foo2$S1 = (foo2$S1T1-1)%%ns+1
  foo2$T1 = floor((foo2$S1T1-1)/ns)+1
  foo2$cumDailyNsoDeathsByDeathDate = 0
  foo2$date = min(d$date)+foo2$T1-1
  foo2$areaName = cname[foo2$S1]
  d2 = rbind(d,foo2[,pmatch(names(d),names(foo2))])
  rm(foo,foo2)
}
d2 = data.table(d2)
d2 = d2[order(S1T1)]
d2$cases = NA
for(i in 1:ns)
{
  cases = d2$cumDailyNsoDeathsByDeathDate[d2$S1==i]
  d2$cases[d2$S1==i] = c(cases[1],diff(cases))
}
pop = c(2646772,7422295,5481431,4880094,5964240,6348096,8796628,9294023,5712840)

##Making adjency matrix
Adj = matrix(0,nrow=ns,ncol=ns)
Adj[1,c(2,3)] = 1
Adj[2,c(1,3,4,5)] = 1
Adj[3,c(1,2,4)] = 1
Adj[4,c(2,3,4,5,6,8)] = 1
Adj[5,c(2,4,8,9)] = 1
Adj[6,c(4,7,8)] = 1
Adj[7,c(6,8)] = 1
Adj[8,c(4,5,6,7,9)] = 1
Adj[9,c(5,8)] = 1

coviduk = list(data=d2,adj=Adj,pop=pop)
saveRDS(coviduk,file="data/covidUK.RDS")
