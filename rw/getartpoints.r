# title         : getartpoints.R
# purpose       : get  points from outcrop map and slope map;
# reference     :
# producer      : Prepared by W. Shangguan 
# address       : In Wageningen, NL.
# inputs        : outcrop map and slope map
# outputs       : outcrop points and slope points;
# remarks 1     : Takes ca ??? mins to run with 10 cpus in use;

rm(list = ls(all = TRUE))
library(rgdal)
library(plotKML)
library(gdalUtils)
library(plyr)
library(foreign)
library(ggplot2)
library(snowfall)

w.dir <- "E:/data/soildata/depth"
wg.dir <- "/home/shang009/big/worldgrids"
#wg.dir <- "D:/worldgrids/1km"
dir_f <- "E:\\data\\soildata\\depth\\points\\codegsifb\\head\\"
den.area <- 1e7  # get one point within den.area, square meter
source(paste(dir_f, "functions.r", sep = ""))



kmloutcrop <- function(tmp,tname,i)
{
    kml(tmp, colour = 1,
        raster_name = paste0(tname, i, ".png"),
        colour_scale = SAGA_pal[[1]],
        folder.name =   paste0(tname, i, ".kml" ),
        file.name = paste0(tname, i, ".kml"))
    flist <-  list.files(pattern = tname)
    file.copy(from=flist, to = paste0(w.dir, "/outcropkml/"), overwrite = T)
    file.remove(flist)
}

getsp <- function(obj, den.area, tname, i){
    p.area <- sapply(slot(obj, "polygons"), function(x) sapply(slot(x, "Polygons"), function(x)
            {
                p1 <- 0
                if(slot(x, "hole")== FALSE)
                {
                    p1 <- slot(x, "area")
                }
                return(p1)
            }))

    p.area <- unlist(p.area)
    p.area <- sum(p.area)
    p.num  <- round(p.area/den.area)
    ret <- spsample(obj,  p.num , "random")
    ret <- reproject(ret)
    source <- paste0(tname, i)
    ret <- data.frame(source,ret)
    return(ret)
}
getsp2 <- function(obj, den.area, tname, i){
    p.area <- obj@grid@cellsize[1]*obj@grid@cellsize[2] * length(obj)
    p.num  <- round(p.area/den.area)
    ret <- spsample(obj,  p.num , "random")
    ret <- reproject(ret)
    source <- paste0(tname, i)
    ret <- data.frame(source,ret)
    return(ret)
}
obj<-tmp

get.equalarea <- function(bbox){

    listWorld <- data.frame(projection=c(
        #//Equal-area world map projections with poles represented as points
        "Mollweide",
        "Hammer (or Hammer-Aitoff)",
        "Boggs Eumorphic",
        "Sinusoidal",
        #//Equal-area world map projections with poles represented as lines
        "Eckert IV",
        "Wagner IV (or Putnins P2')",
        "Wagner VII (or Hammer-Wagner)",
        "McBryde-Thomas flat-polar quartic",
        "Eckert VI",
         #//Equal-area interrupted projections for world maps with poles represented as points
         "Mollweide",
         "Boggs Emorphic",
         "Goode homolosine",
         "Sinusoidal",
         #//Equal-area interrupted projections for world maps with poles represented as lines
         "McBryde S3"
    ),PROJ4 = c(
        "+proj=moll",
        "+proj=hammer",
        "+proj=boggs",
        "+proj=sinu",
        "+proj=eck4",
        "+proj=wag4",
        "+proj=wag7",
        "+proj=mbtfpq",
        "+proj=eck6",
        "nocode",
        "nocode",
        "+proj=igh",
        "nocode",
        "nocode"
    ))
    lon.min <- bbox[1,1]
    lon.max <- bbox[1,2]
    lat.min <- bbox[2,1]
    lat.max <- bbox[2,2]
    lon.cen <- mean(bbox[1,])
    lat.cen <- mean(bbox[2,])
    scale <- 360/(lon.max - lon.min) * 180/ (lat.max - lat.min)
    if (scale < 1.5){ #World (small-scale) map
        ret <- NULL
        #//Equal-area world map projections with poles represented as points
        for (i in 1:4) {
			ret <- c(ret, paste0(listWorld[i,2], " +lon_0=", center.lng))
		}
        #//Equal-area world map projections with poles represented as lines
		for (i in 5:9) {
			ret <- c(ret, paste0(listWorld[i,2], " +lon_0=", center.lng))
		}
        #//Equal-area interrupted projections for world maps with poles represented as points
		for (i in 10:13) {
			ret <- c(ret, listWorld[i,2])
		}
        #//Equal-area interrupted projections for world maps with poles represented as lines
        ret <- c(ret, listWorld[13,2])
    }else if ((scale < 6) || ((abs(lon.max - lon.min) > 200) && (abs(lat.cen) > 15)) ){
        #/***PRINTING HEMISPHERE MAP PROJECTIONS***/
        lon <- round(lon.cen * 10) / 10
    	if (lat.cen > 85) {
    		lat = 90
    	} else if (lat.cen < -85) {
    		lat <- -90
    	} else {
    		lat <- round(lat.cen * 10) / 10
    	}
        centerPROJ4 <- paste0(" +lat_0=", lat, " +lon_0=", lon)
        ret <- paste0("+proj=aea", centerPROJ4)
    }else{
        #/***PRINTING LARGE-SCALE MAP PROJECTIONS***/
        ratio <- (lat.max - lat.min) / (lon.max - lon.min)
        if (ratio > 1.25){ #/*Funcion for regional maps with an north-south extent*/
                ret <- paste0("+proj=tcea +lon_0=", lon.cen)
        }else if (ratio < 0.8){#/*Funcion for regional maps with an east-west extent*/
            if ((lat.cen > 67.5 && scale < 8) || lat.cen > 70){
                ret <- paste0("+proj=aea +lat_0=90.0 +lon_0=", lon.cen)
            }else if ((lat.cen < -67.5 && scale < 8) || lat.cen < -70) {
                ret <- paste0("+proj=aea +lat_0=-90.0 +lon_0=", lon.cen)
            }else if (lat.cen > -15. && lat.cen < 15.){
                interval <- (lat.max - lat.min)/4
    			latS1 <- lat.cen + interval
                latS2 <- lat.cen - interval
                if ((latS1 > 0 && latS2 > 0) || (latS1 < 0 && latS2 < 0))
    				latS <- max(abs(lat.max), abs(lat.min))/2
    			else
    				latS <- 0
                ret <- paste0("+proj=cea +lat_ts=" + latS + " +lon_0=", lon.cen)
            }else{
                interval <- (lat.max - lat.min)/6
                ret <- paste0("+proj=aea +lat_1=", (lat.min + interval), " +lat_2=", (lat.max - interval), " +lon_0=", lon.cen)
            }
        }else { #/*Funcion for regional maps in square format*/
            centerPROJ4 <- paste0(" +lat_0=", lat.cen, " +lon_0=", lon.cen)
            if(lat.cen > 75){
                ret <- paste0("+proj=aea +lat_0=90.0 +lon_0=", lon.cen)
            }else if (lat.cen < -75){
                ret <- paste0("+proj=aea +lat_0=-90.0 +lon_0=", lon.cen)
            }else if (lat.cen > -15 && lat.cen < 15) {		
        		ret <- paste0("+proj=aea +lon_0=", lon.cen)		
        	}else{
                ret <- paste0("+proj=aea", centerPROJ4)
            }
        }
    }

    return(ret)
}




outpoint <- NULL
#newyork
setwd(paste0(w.dir, "/outcrop/newyork"))
do_unzip()
tname <- "newyork_crop"
flist <- list.files(path = "./tmp", pattern = ".SHP")
flist <- substr(flist,1,6)
for(i in 1:5)
{
    tmp <- readOGR("./tmp", flist[i])
    tmp <- subset(tmp, tmp$MATERIAL == "r")
    proj4string(tmp) <- "+proj=utm +zone=18 +ellps=clrk66 +units=m +no_defs"
    #kmloutcrop(tmp,tname,i)
    outpoint <- rbind(outpoint,getsp(tmp, den.area, tname, i))

}
del_unzip()
# not captured by the prediction



#vermont
setwd(paste0(w.dir, "/outcrop/vermont"))
do_unzip()
tname <- "vermont_crop"
tmp <- readOGR("./tmp", "Geologic_SURFICIAL62K_poly")
tmp <- subset(tmp, tmp$LITHCODE == "r")
i <- 1
#kmloutcrop(tmp,tname,i)
outpoint <- rbind(outpoint,getsp(tmp, den.area, tname, i))
del_unzip()
# not captured by the prediction


#Alaska
setwd(paste0(w.dir, "/outcrop/alaska"))
do_unzip()
tname <- "Alaska_crop"
#ogr2ogr("gaar_geo", "bbb") # not working, converted to shp by ArcGIS
tmp <- readOGR("./tmp", "gaar_geo")
tmp <- subset(tmp, tmp$polygon_UN == "BEDROCK SURFACE FORMS")
i <- 1
outpoint <- rbind(outpoint,getsp(tmp, den.area, tname, i))
#kmloutcrop(tmp,tname,i)
del_unzip()

#Alberta
setwd(paste0(w.dir, "/outcrop/alberta"))
do_unzip()
tname <- "Alberta_crop"
tmp <- readOGR("./tmp", "surf1m_py_ll")
tmp <- subset(tmp, tmp$MAP_LABEL == "R")
i <- 1
#kmloutcrop(tmp,tname,i)
tmp <- reproject(tmp, get.equalarea(tmp@bbox))
outpoint <- rbind(outpoint,getsp(tmp, den.area, tname, i))
del_unzip()
#tmp$gid <- 1:length(tmp)
#tmp <- subset(tmp, tmp$gid==1093)
#df.sp <- data.frame(sp)
#county <- tmp
#county@data$id <- rownames(county@data)
#county.points <- fortify(county, coords="id")
#county.df <- join(county.points, county@data, by="id")
#ggplot()+
#geom_polygon(county.df,mapping=aes(long,lat,group=group, fill="hole"), color="white" ) +
#geom_path(county.df, mapping=aes(long,lat,group=group), color="white") +
#coord_equal() +
#geom_point(df.sp, mapping=aes(x=x,y=y), size=3)



#manitoba
setwd(paste0(w.dir, "/outcrop/manitoba"))
do_unzip()
tname <- "Manitoba_crop"
tmp <- readOGR("./tmp", "SGCMS_1_1Million")
tmp2 <- read.dbf("./tmp/LegLink.dbf")
tmp2 <- subset(tmp2, tmp2$MATERIAL == "Rock")
tmp  <- subset(tmp, tmp$Link %in% tmp2$LINK)
i <- 1
#kmloutcrop(tmp,tname,i)
outpoint <- rbind(outpoint,getsp(tmp, den.area, tname, i))
del_unzip()



#Newfoundland
setwd(paste0(w.dir, "/outcrop/Newfoundland"))
#"drift poor":greater then 80% bedrock; includes areas of till and other surficial materials generally less than 1m thick and discontinuous;
do_unzip()
tname <- "Newfoundland_crop"
tmp <- readOGR("./tmp/zipfolder", "Regional_Surficial_Geology")
tmp <- subset(tmp, tmp$GENETIC1MA == "exposed bedrock" | tmp$GENETIC1MA == "drift poor")
i <- 1
#kmloutcrop(tmp,tname,i)
tmp <- reproject(tmp, get.equalarea(tmp@bbox))
outpoint <- rbind(outpoint,getsp(tmp, den.area, tname, i))
del_unzip()

#NRCan
setwd(paste0(w.dir, "/outcrop/NRCan"))
do_unzip()
tname <- "NRCan_crop"
flist <-  list.files(path= "./tmp", pattern = ".shp")
flist <- strtrim(flist,nchar(flist)-4)
for(i in 1 :10)
{
    tmp <- readOGR("./tmp", flist[i])
    tmp <- subset(tmp, tmp$gin_lithol %in% c("R","R1", "R3"))
    proj4string(tmp) <- "+proj=longlat +datum=WGS84"
    #kmloutcrop(tmp,tname,i)
    tmp <- reproject(tmp,get.equalarea(tmp@bbox))
    outpoint <- rbind(outpoint,getsp(tmp, den.area, tname, i))
}
del_unzip()

#Australia
setwd(paste0(w.dir, "/Australia/ngis/regolithmap"))
tname <- "Australia"
tmp <- readGDAL("DER_1km.tif")
tmp <- subset(tmp, tmp$band1 < 0.2) #less than 1 meter
tmp <- reproject(tmp,get.equalarea(tmp@bbox))
gdalinfo("DER_1km.tif")
writeGDAL(tmp, fname = "oucrop.tif", drivername = "GTiff", type = "Float32", mvFlag = -9999)
i <- 1
tmp <- readGDAL("oucrop.tif")
tmp <- as(tmp,  "SpatialPixelsDataFrame")
outpoint <- subset(outpoint, as.character(outpoint$source) != paste0(tname,i))
outpoint <- as.data.frame(outpoint)
outpoint <- rbind(outpoint,getsp2(tmp, den.area, tname, i))
count(outpoint, "source")
coordinates(outpoint) <- ~ x + y
proj4string(outpoint) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
save(outpoint, file = paste0(w.dir, "/outcrop/outpoint.rda"))
#Slop map
setwd(wg.dir)
tmp  <- readGDAL("SLPSRT3a.tif")
outpoint$SLPSRT3a <- over(outpoint, tmp)$band1
save(outpoint, file =  "../soildepth/outpoint.rda")
#tmp2 <- readGDAL("SMKISR3a.tif")
#tmp  <- subset(tmp, tmp$band1 > (30*255/90) & (tmp2$band1==1|tmp2$band1==2))
#writeGDAL(tmp, fname = "../soildepth/slop30.tif", drivername = "GTiff", type = "Byte", mvFlag = 255)
#rm(tmp)
#tm(tmp2)


wrapper.sample <- function(i)
{
    slppoint <- NULL
    tmp <- readGDAL("../soildepth/slop30.tif", offset = c(0,(i-1)*108), region.dim = c(21600, 108))
    if(sum(!is.na(tmp$band1)) > 1){
        tmp <- as.data.frame(tmp, row.names = NULL)
        p.num <- round(dim(tmp)[1]*1e6/den.area)  # assuming each grid is 1km*1km,
        spn <- sample(dim(tmp)[1], p.num)
        if(length(spn) > 0) slppoint <- tmp[spn,]
    }
    save(slppoint,file=paste0("../soildepth/slppoint", i, ".rda"))
}


nRuns <- 1:400
nSlaves <- 10


sfInit(parallel=TRUE, cpus = nSlaves, slaveOutfile="~/errorwg.log")
sfLibrary(sp)
sfLibrary(rgdal)
sfExport("den.area")
ptm <- proc.time()
sfClusterApplyLB(nRuns, wrapper.sample)
proc.time() - ptm
sfStop()
tmp <- NULL
for(i in nRuns)
{
    if(file.exists(paste0("../soildepth/slppoint", i, ".rda")))
    {
        load(paste0("../soildepth/slppoint", i, ".rda"))
        tmp <- rbind(tmp, slppoint)
    }
}


slppoint <- tmp
coordinates(slppoint) <- ~ x + y
proj4string(slppoint) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
names(slppoint)[1] <- "SLPSRT3a"
slppoint$source <- "slope"
save(slppoint, file =  "../soildepth/slppoint.rda")


for(i in nRuns)
{
    file.remove(paste0("../soildepth/slppoint", i, ".rda"))
}
