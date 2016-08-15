#take out the extreme observations on the slope >40

library(rgdal)
w.dir <- "/data/shang009/big/soildepth"
setwd(w.dir)

load("./profs/wells.depth.rda")
slp <- readGDAL("./profs/slop30.tif")
proj4string(slp)<-proj4string(wells.depth)

tmp <- subset(wells.depth, wells.depth$BDRICM >200)
tmp$slp <- over(tmp,slp )$band1
sum(tmp$slp>0,na.rm=T)
tmp<-subset(tmp, is.na(tmp$slp))
tmp$slp<-NULL

wells.depth<- rbind(tmp,subset(wells.depth, wells.depth$BDRICM <=200))
save(wells.depth, file = "./profs/wells.depth2.rda")

