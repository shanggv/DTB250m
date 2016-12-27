## Compare BDTICM based on regional studies and SoilGrids predictions

library(rgdal)
library(raster)
library(plotKML)
library(maps)
library(maptools)
library(psych)
library(scales)
usa.m <- map('state', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(usa.m$names, ":"), function(x) x[1])
prj0 ="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "
prj = "+proj=utm +zone=15 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
state = as(map2SpatialPolygons(usa.m, IDs=IDs), "SpatialLines")
proj4string(state) = "+proj=longlat +datum=WGS84"
state <- spTransform(state, CRS(prj))
state.n <- sapply(state@lines, function(x){slot(x, "ID")})

if(.Platform$OS.type == "windows"){
  gdal.dir <- shortPathName("C:/Program files/GDAL")
  gdal.dir <- shortPathName("D:\\GDAL\\bin\\gdal\\apps")
  gdal_translate <- paste0(gdal.dir, "/gdal_translate.exe")
  gdalbuildvrt <- paste0(gdal.dir, "/gdalbuildvrt.exe")
  gdalwarp <- paste0(gdal.dir, "/gdalwarp.exe")
} else {
  gdalwarp = "/usr/local/bin/gdalwarp"
  gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
  gdal_translate =  "/usr/local/bin/gdal_translate"
}
load("/data/shang009/big/soildepth/profs/wells.depth.rda")

setwd("/data/shang009/big/soildepth2")
## Regional maps:
grd.lst <- c("Ken.tif", "Pen.tif")
ext <- rbind(c(-90, -81.5, 36.25, 39.5 ), c(-81, -74.5, 39.5, 43))
ext
## SoilGrids-based BDTICM (predictions at 250 m):
for(j in 1:length(grd.lst)){
  out_GSRS = gsub(".tif", "_GSRS.tif", grd.lst[j])
  out_utm = gsub(".tif", "_utm.tif", grd.lst[j])
  if(!file.exists(grd.lst[j])){
    unlink(grd.lst[j])
    system(paste0(gdalwarp, ' BDTICM_M_1km_ll.tif ', grd.lst[j], ' -te ', paste(ext[j,][c(1,3,2,4)], collapse=" "), ' -tr 0.008333333 0.008333333 -t_srs \"', prj0, '\" -co \"COMPRESS=DEFLATE\"'))
    system(paste0(gdalwarp, ' ', grd.lst[j], ' ', out_utm,  ' -tr 500 500 -t_srs \"', prj, '\" -co \"COMPRESS=DEFLATE\"'))
    r <- raster(out_utm)
    system(paste0(gdalwarp, ' average_soil_and_sedimentary-deposit_thickness.tif ', out_GSRS, ' -tr ', res(r)[1], ' ', res(r)[2], ' -t_srs \"', proj4string(r), '\" -co \"COMPRESS=DEFLATE\" -te ', paste(extent(r)[c(1,3,2,4)], collapse=" ")))
    }
}

Ken_p <- subset(wells.depth, wells.depth@coords[,1] > ext[1,1] & wells.depth@coords[,1] < ext[1,2] & wells.depth@coords[,2] > ext[1,3] & wells.depth@coords[,2] < ext[1,4])
Ken_p  <- spTransform(Ken_p , CRS(prj))
Ken_p$BDRICM2<-Ken_p$BDRICM
Ken_p$BDRICM2[Ken_p$BDRICM>5000] <-5000
Ken_p$BDRICM2[Ken_p$BDRICM<100] <-100
Ken_p$BD <- log1p(Ken_p$BDRICM2)

Pen_p <- subset(wells.depth, wells.depth@coords[,1] > ext[2,1] & wells.depth@coords[,1] < ext[2,2] & wells.depth@coords[,2] > ext[2,3] & wells.depth@coords[,2] < ext[2,4])
Pen_p  <- spTransform(Pen_p , CRS(prj))
Pen_p$BDRICM2<-Pen_p$BDRICM
Pen_p$BDRICM2[Pen_p$BDRICM>5000] <-5000
Pen_p$BDRICM2[Pen_p$BDRICM<100] <-100
Pen_p$BD <- log1p(Pen_p$BDRICM2)


## plot next to each other:
Ken <- stack(c("Ken_GSRS.tif","Ken_utm.tif"))
Ken <- as(as(Ken, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
proj4string(Ken)  <- prj
Ken$Ken_GSRS[Ken$Ken_GSRS==255] <-NA
Ken$Ken_GSRS <- Ken$Ken_GSRS*100
rn = quantile(c( Ken$Ken_utm, Ken$Ken_GSRS), c(.01, .99), na.rm=TRUE)
rx = rev(as.character(round(c(round(rn[1], 0), NA, round(mean(rn), 0), NA, round(rn[2], 0)), 2)))
Ken$Ken_utmf <- ifelse(Ken$Ken_utm<rn[1], rn[1], ifelse(Ken$Ken_utm>rn[2], rn[2], Ken$Ken_utm))
Ken$Ken_GSRSf <- ifelse(Ken$Ken_GSRS<rn[1], rn[1], ifelse(Ken$Ken_GSRS>rn[2], rn[2], Ken$Ken_GSRS))
jpeg(file="Fig_Ken_comparison0.jpg", res=150, width=900, height=800)
par(mai=c(0.8,0.8,0.5,0), oma=c(0,0,0,0),xaxs='i', yaxs='i')
spplot(Ken_p, zcol = "BD",col.regions=SAGA_pal[[1]],main="kentucky", cex=0.4,colorkey = FALSE,  sp.layout=list(state[which(state.n=="kentucky"),], col = 'black', lwd = 1), auto.key=FALSE, scales=list(draw=T), xlab="Easting (m)", ylab="Northing (m)" ,add =T)
dev.off()
jpeg(file="Fig_Ken_comparison.jpg", res=150, width=750, height=750*2)
#spplot(Ken, col.regions=SAGA_pal[[1]])
par(mfrow=c(2,1))
par(mai=c(0.8,0.8,0.5,0), oma=c(0,0,0,0),xaxs='i', yaxs='i')
image(log1p(raster(Ken["Ken_utmf"])), col=SAGA_pal[[1]], zlim=log1p(rn), main="SoilGrids", asp=1, axes=T, xlab="Easting (m)", ylab="Northing (m)") # , cex.lab=.7, cex.axis=.7
lines(state[which(state.n=="kentucky"),], lwd=2)
scalebar(50000, type="bar", xy=c(1300000, 3950000), divs=2, cex=.8, cex.axis=.8, below="m")
legend(x=1380000, y=4050000, rx, fill=rev(SAGA_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, bty="n", cex=.8)
image(log1p(raster(Ken["Ken_GSRSf"])), col=SAGA_pal[[1]], zlim=log1p(rn), main="Pelletier et al. (2016)", asp=1, axes=T, xlab="Easting (m)", ylab="Northing (m)") # cex.lab=.7, cex.axis=.7)
lines(state[which(state.n=="kentucky"),], lwd=2)
scalebar(50000, type="bar", xy=c(1300000, 3950000), divs=2, cex=.8, cex.axis=.8, below="m")
dev.off()


## Scatter plot histograms:
rn[2]<-5000

Ken_p2 <- Ken_p[sample.int(length(Ken_p),20000),] 
Ken_p2$BD_utm <- over( Ken_p2,Ken)$Ken_utm
Ken_p2$BD_GSRS <- over( Ken_p2,Ken)$Ken_GSRS
Ken_p2$BD_GSRS2 <- ifelse(Ken_p2$BD_GSRS>rn[2], NA, Ken_p2$BD_GSRS)
Ken_p2$BDRICM2 <- ifelse(Ken_p2$BDRICM>=rn[2], NA, Ken_p2$BDRICM2)
jpeg(file="Fig_Ken_scatterplot1.jpg", res=150, width=750, height=750)
#spplot(Ken, col.regions=SAGA_pal[[1]])
par(mfrow=c(2,1))
par(mai=c(0.8,0.8,0.5,0), oma=c(0,0,0,0),xaxs='i', yaxs='i')
with(Ken_p2@data,scatter.hist(BDRICM,BD_utm, xlab="Observation", ylab="SoilGrids", pch=19, col="lightblue", cex=1.5))
dev.off()
jpeg(file="Fig_Ken_scatterplot2.jpg", res=150, width=750, height=750)
#spplot(Ken, col.regions=SAGA_pal[[1]])
par(mfrow=c(2,1))
par(mai=c(0.8,0.8,0.5,0), oma=c(0,0,0,0),xaxs='i', yaxs='i')
with(Ken_p2@data,scatter.hist(BDRICM2,BD_GSRS2, xlab="Observation", ylab="Pelletier et al. (2016)", pch=19, col="lightblue", cex=1.5),correl=T)
dev.off()






## plot next to each other:
Pen <- stack(c("Pen_GSRS.tif","Pen_utm.tif"))
Pen <- as(as(Pen, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
proj4string(Pen)  <- prj
Pen$Pen_GSRS[Pen$Pen_GSRS==255] <-NA
Pen$Pen_GSRS <- Pen$Pen_GSRS*100
rn = quantile(c(Pen$Pen_GSRS, Pen$Pen_utm), c(.01, .99), na.rm=TRUE)
rx = rev(as.character(round(c(round(rn[1], 0), NA, round(mean(rn), 0), NA, round(rn[2], 0)), 2)))
Pen$Pen_utmf <- ifelse(Pen$Pen_utm<rn[1], rn[1], ifelse(Pen$Pen_utm>rn[2], rn[2], Pen$Pen_utm))
Pen$Pen_GSRSf <- ifelse(Pen$Pen_GSRS<rn[1], rn[1], ifelse(Pen$Pen_GSRS>rn[2], rn[2], Pen$Pen_GSRS))
jpeg(file="Fig_Pen_comparison0.jpg", res=150, width=850, height=800)
par(mai=c(0.8,0.8,0.5,0), oma=c(0,0,0,0),xaxs='i', yaxs='i')
spplot(Pen_p, zcol = "BD",col.regions=SAGA_pal[[1]],main="Pennsylvania", cex=0.4,colorkey = FALSE,  sp.layout=list(state[which(state.n=="pennsylvania"),], col = 'black', lwd = 1), auto.key=FALSE, scales=list(draw=T), xlab="Easting (m)", ylab="Northing (m)" ,add =T)
dev.off()

jpeg(file="Fig_Pen_comparison.jpg", res=150, width=750, height=800*2)
#spplot(Pen, col.regions=SAGA_pal[[1]])
par(mfrow=c(2,1))
par(mai=c(0.8,0.8,0.5,0), oma=c(0,0,0,0),xaxs='i', yaxs='i')
image(log1p(raster(Pen["Pen_utmf"])), col=SAGA_pal[[1]], zlim=log1p(rn), main="SoilGrids", asp=1, axes=T, xlab="Easting (m)", ylab="Northing (m)") # , cex.lab=.7, cex.axis=.7
lines(state[which(state.n=="pennsylvania"),], lwd=2)
scalebar(50000, type="bar", xy=c(1800000, 4415000), divs=2, cex=.8, cex.axis=.8, below="m")
legend(x=2000000, y=4540000, rx, fill=rev(SAGA_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, bty="n", cex=.8)
image(log1p(raster(Pen["Pen_GSRSf"])), col=SAGA_pal[[1]], zlim=log1p(rn), main="Pelletier et al. (2016)", asp=1, axes=T, xlab="Easting (m)", ylab="Northing (m)") # cex.lab=.7, cex.axis=.7)
lines(state[which(state.n=="pennsylvania"),], lwd=2)
scalebar(50000, type="bar", xy=c(1800000, 4415000), divs=2, cex=.8, cex.axis=.8, below="m")
dev.off()

## Scatter plot histograms:
Pen_p2 <- Pen_p[sample.int(length(Pen_p),20000),] 
Pen_p2$BD_utm <- over( Pen_p2,Pen)$Pen_utm
Pen_p2$BD_GSRS <- over( Pen_p2,Pen)$Pen_GSRS
Pen_p2$BD_GSRS2 <- ifelse(Pen_p2$BD_GSRS>rn[2], NA, Pen_p2$BD_GSRS)
Pen_p2$BDRICM2 <- ifelse(Pen_p2$BDRICM>=rn[2], NA, Pen_p2$BDRICM2)
jpeg(file="Fig_Pen_scatterplot1.jpg", res=150, width=750, height=750)
#spplot(Pen, col.regions=SAGA_pal[[1]])
par(mfrow=c(2,1))
par(mai=c(0.8,0.8,0.5,0), oma=c(0,0,0,0),xaxs='i', yaxs='i')
with(Pen_p2@data,scatter.hist(BDRICM,BD_utm, xlab="Observation", ylab="SoilGrids", pch=19, col="lightblue", cex=1.5))
dev.off()
jpeg(file="Fig_Pen_scatterplot2.jpg", res=150, width=750, height=750)
#spplot(Pen, col.regions=SAGA_pal[[1]])
par(mfrow=c(2,1))
par(mai=c(0.8,0.8,0.5,0), oma=c(0,0,0,0),xaxs='i', yaxs='i')
with(Pen_p2@data,scatter.hist(BDRICM2,BD_GSRS2, xlab="Observation", ylab="Pelletier et al. (2016)", pch=19, col="lightblue", cex=1.5),correl=T)
dev.off()

