library(data.table)
#acounty = data.table(read.csv("data/hosp_fylke.csv"))
#acounty = read.csv("data/covid_msis_latest.csv",sep=",")

require(RCurl)


start_date = as.Date("2020-06-03")
end_date = as.Date("2022-11-11")
n = as.numeric(end_date-start_date)
dir = "https://raw.githubusercontent.com/folkehelseinstituttet/surveillance_data/master/covid19/data_covid19_msis_by_location_"
d = NULL
i=1
file = paste(dir,as.character(start_date+i),".csv",sep="")
foo1 = read.csv(file,sep=",",nrows=11)[,c("location_code","n","pop","location_name")]
foo1$date = start_date+i
d = foo1
foo1$n = NA

for(i in 1:n)
{
  file = paste(dir,as.character(start_date+i),".csv",sep="")
  show(start_date+i)
  foo = foo1
  if(url.exists(file))
  {
   foo = read.csv(file,sep=",",nrows=11)[,c("location_code","n","pop","location_name")]
  }
  foo$date = start_date+i
  d = rbind(d,foo)
}
##Fill in gaps
counties = unique(d$location_code)
dfill = NULL
for(i in 1:length(counties))
{
 d2 = d[d$location_code==counties[i],]
 ind = c(1:nrow(d2))[is.na(d2$n)]
 while(length(ind)>0)
 {
  ind1 = ind[1]
  n=1
  while((n<length(ind)) & (ind[n+1]==(ind1+n)))
   n=n+1
  a=d2$n[ind1-1]
  b=d2$n[ind1+n]
  for(j in 1:n)
  {
   d2$n[ind1-1+j] = as.integer((j*b+(n-j+1)*a)/(n+1))
  }
  ind = ind[-c(1:n)]
 }
 d2$cases = c(0,diff(d2$n))
 dfill = rbind(dfill,d2)
}

#Ad hoc, putting negative cases to zero!!
dfill$cases[dfill$cases<0]=0


saveRDS(d,file="Covid_Norway_test_by_county_and_date.RDS")
saveRDS(dfill,file="Covid_Norway_test_by_county_and_date_impute.RDS")

d3= d[d$location_code==counties[3],]
d3$diff = c(0,diff(d3$n))
c(1:nrow(d3))[d3$diff<0]
