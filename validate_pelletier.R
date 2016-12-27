library(psych)

load("~/big/soildepth/profs/wells.depth.rda")
load("~/big/soildepth/profs/sprofs.depth.rda")
ptif <- readGDAL("../average_soil_and_sedimentary-deposit_thickness.tif")
proj4string(ptif) <- proj4string(wells.depth)
dtb0 <- wells.depth$BDRICM
dtb0[dtb0>=5000] <-NA
dtb <- over(wells.depth, ptif)$band1
dtb[dtb==255] <-NA
dtb[dtb==50] <-NA
range(dtb,na.rm=T)
dtb <- dtb*100
dtb2 <- cbind(dtb,dtb0)
dtb2 <- dtb2[complete.cases(dtb2),]
dim(dtb2)
cor(dtb2[,1],dtb2[,2])^2


#take out the interpolation area
usa.m <- map('state', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(usa.m$names, ":"), function(x) x[1])
state <- map2SpatialPolygons(usa.m, IDs=IDs)
proj4string(state) <-proj4string(wells.depth)
wells.depth$s <- over(wells.depth,state)
for(i in 1:49)
{
  print(i)
  print(state[i]@polygons[[1]]@ID)
}
#13indiana,16kentucky,31new york,37pennsylvania,

wells.sub1 <- subset(wells.depth, wells.depth$s %in% c(13,16,31,37))
wells.sub2 <- subset(wells.depth, !(wells.depth$s %in% c(13,16,31,37)))



dtb0 <- wells.sub1$BDRICM
dtb0[dtb0>=5000] <-NA
dtb <- over(wells.sub1, ptif)$band1
dtb[dtb==255] <-NA
dtb[dtb==50] <-NA
range(dtb,na.rm=T)
dtb <- dtb*100
dtb2 <- cbind(dtb,dtb0)
dtb2 <- dtb2[complete.cases(dtb2),]
dim(dtb2)
cor(dtb2[,1],dtb2[,2])^2
mean(log1p(dtb2[,1]) - log1p(dtb2[,2]), na.rm=TRUE)
scatter.hist(dtb2[,2],dtb2[,1], xlab="Observation", ylab="Pelletier et al. (2016)", pch=19, col="lightblue", cex=1.5)


dtb0 <- wells.sub2$BDRICM
dtb0[dtb0>=5000] <-NA
dtb <- over(wells.sub2, ptif)$band1
dtb[dtb==255] <-NA
dtb[dtb==50] <-NA
range(dtb,na.rm=T)
dtb <- dtb*100
dtb2 <- cbind(dtb,dtb0)
dtb2 <- dtb2[complete.cases(dtb2),]
dim(dtb2)
cor(dtb2[,1],dtb2[,2])^2
mean(log1p(dtb2[,1]) - log1p(dtb2[,2]), na.rm=TRUE)
