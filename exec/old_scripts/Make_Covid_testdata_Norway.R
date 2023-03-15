library(data.table)
#acounty = data.table(read.csv("data/hosp_fylke.csv"))
acounty = readRDS("data/Covid_Norway_test_by_county_and_date_impute.RDS")

cname = c("Oslo","Rogaland","Møre og Romsdal","Nordland","Viken","Innlandet",
          "Vestfold og Telemark","Agder","Vestland","Trøndelag","Troms og Finnmark")
county = paste("county",c(3,11,15,18,30,34,38,42,46,50,54),sep="")
county[1] = "county03"



##Making adjency matrix
Adj = matrix(0,nrow=11,ncol=11)
Adj[1,c(5)] = 1
Adj[2,c(8,9)] = 1
Adj[3,c(6,9,10)] = 1
Adj[4,c(10,11)] = 1
Adj[5,c(1,6,7,9)] = 1
Adj[6,c(3,5,9,10)] = 1
Adj[7,c(5,8,9)] = 1
Adj[8,c(2,7)] = 1
Adj[9,c(2,3,5,6,7)] = 1
Adj[10,c(3,4,6)] = 1
Adj[11,c(4)] = 1

acounty = acounty[,c("date","location_code","cases","pop")]

coviddata = list(data=acounty,adj=Adj)
saveRDS(coviddata,file="data/coviddata.RDS")





