#take hourly files and convert them into daily for summer months of 2008
#did not adjust for local time zone

library(reshape2)
library(plyr)
setwd("M:/NET MyDocuments/Projects/FireResearch/SAMSI")


#summer months only
#datadir = "HYSPLIThourlydataSummer2008/2008_RAW_HourlyData_ugm3/ascii/"

datadir = "2008_RAW_HourlyData_ugm3/ascii/"

filenms = list.files(datadir)
dataoutdir = "HYSPLIT_fire_daily_2008/"
#dataoutdir = "HYSPLIT_fire_daily_summer2008/"

dt = NULL
for ( i in 1:length (filenms)){
   dt = read.csv(paste(datadir,filenms[i],sep=""))
   dt$PM2500100 = dt$PM2500100*1000000000
   dtm = melt(dt, id.vars=c("LAT","LON","HR","DA","YEAR", "MO"))
   davg = dcast(dtm,LAT+LON ~ variable, fun.aggregate= mean)
   
   write.csv(davg, file=paste(dataoutdir, "hysplit_pm_2008-",
                   substr(filenms[i], 20,21), "-",
                    substr(filenms[i], 22,23), ".csv", sep=""))
   print(i)
   
}

stop()

