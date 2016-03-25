rm(list = ls(all = TRUE))

library(sp)
library(rgdal)
library(plyr)
#define global, change if needed
dir_f <- "E:\\data\\soildata\\depth\\points\\codegsifb\\head\\"
source(paste(dir_f, "functions.r", sep = ""))
dir1  <- "E:\\data\\soildata\\depth\\points\\europe\\"
c_names  <- c("Source", "Long", "Lat", "D_BR", "D_water", "D_well", "Accu_xy")
setwd(dir1)
dirs   <- shell("dir * /B", intern = T)
states <- cbind(1:length(dirs), dirs)
colnames(states) <- c("s_code", "states")
states
dir_out <- "E:\\data\\soildata\\depth\\points\\"#      s_code states


rec.lst <- as.list(rep(NA, length(dirs)))
names(rec.lst) <- dirs

#-------------------------------Ireland, meter
#0 are  kept,
s_n <- 1#number of source
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
tmp <- readOGR("./tmp", "Geotech_Borehole_Bedrock_Met")
tmp <- spTransform(tmp, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
tmp <- as.data.frame(tmp)
tmp <- cbind(tmp$x, tmp$y, tmp$DTB_Bdk_M, NA, tmp$DEPTH_HOLE, "")
rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
#"108 out of 4250 is 0"
del_unzip()
#-------------------------------Sweden, meter
#0 are  kept,
s_n <- 2#number of source
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
flist <- list.files("./tmp", pattern = ".csv")
vnames <- c("BRUNNS_ID",  "E", "N", "DJUP_TILL_BERG", "GRUNDVATTENNIVA", "TOTALDJUP", "LAGESNOGGRANNHET" , "TJ")
tmp <- NULL
for(i in 1:length(flist))
{
    tmp2 <- read.csv(paste0("./tmp/",flist[i] ), sep = ";")[vnames]
    print(i)
    print(dim(tmp2)[1])
    colnames(tmp2) <- c(c_names, "TJ")
    tmp2 <- tmp2[(!is.na(tmp2[,"Long"]) & tmp2[,"Long"] != "")&
        (!is.na(tmp2[,"Lat"]) & tmp2[,"Lat"] != "") &
         (!is.na(tmp2[,"D_BR"]) &tmp2[,"D_BR"] != "") &
         (is.na(tmp2[,"TJ"]) | tmp2[,"TJ"] ==""), ]
    print(dim(tmp2)[1])
    tmp <- rbind(tmp, tmp2)
}
tmp <- tmp[, c( "Long", "Lat", "D_BR", "D_water", "D_well", "Accu_xy")]
tmp <- tmp[tmp$D_BR> -0.001,]
tmp$Long <- as.integer(as.character(tmp$Long))
tmp$Lat <- as.integer(as.character(tmp$Lat))
tmp <- tmp[!is.na(tmp$Long),]
tmp <- tmp[!is.na(tmp$Lat),]

tmp2 <- tmp
coordinates(tmp2) <- ~ Long +Lat
proj4string(tmp2) <- "+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
tmp2 <- spTransform(tmp2, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
tmp2 <- as.data.frame(tmp2)
tmp  <- tmp2[, c( "x", "y", "D_BR", "D_water", "D_well", "Accu_xy")]
rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
#"14222 out of 320451 is 0"
del_unzip()



#-------------------------------merge rec.lst
wells <- rec.lst[[1]]
for(i in 2:length(rec.lst))
{
    wells <- rbind(wells,rec.lst[[i]])

}
setwd(dir1)
#-------------------------------check values
for(i in 1:length(rec.lst))
{
    png(filename = paste0("..\\pic\\", names(rec.lst)[i], ".png"))
    hist(rec.lst[[i]][, 4], breaks = 200,
        xlab = names(rec.lst)[i])
    dev.off()
    png(filename = paste0("..\\pic\\", 2,names(rec.lst)[i], ".png"))
    hist(rec.lst[[i]][, 4],
        breaks = seq(from = -1, to = max(rec.lst[[i]][, 4])+10, by = 1 ),
        xlim = c(0,200), xlab = names(rec.lst)[i])
    dev.off()
    print(names(rec.lst)[i])
    print_0(rec.lst[[i]][, 4])
}

for( i in 4:6)
{
    print(colnames(wells[i]))
    print(max(wells[i], na.rm = TRUE))
    print(min(wells[i], na.rm = TRUE))
}
#wells[wells[,3]<20,]
wells <- wells[wells[, 4] < 4000, ]
print_0(wells[,4])
#"2718 out of 661970 is 0"

#-------------------------------get the list of accuracy of position
#levels(as.factor(wells[, 7]))
#levels(as.factor(wells[wells[, 7] != "", 1]))
#Alaska
s_n = 1
tmp  <- c("0", "1", "2", "3", "8" )
tmp2 <- c("<100m", "<250m", "unknown", "not controlled", "no information")
acc <- cbind(s_n, tmp, tmp2)
colnames(acc) <- c("source", "flag", "descrip")

#-------------------------------output
#rm(rec.lst)
write.table(cbind(states,count(wells$Source)[2]), paste(dir_out, "states_eu.txt", sep = ""),
    row.names = FALSE, sep = "\t")
write.table(acc, paste(dir_out, "acc_eu.txt", sep = ""),
    row.names = FALSE, sep = "\t")
write.table(wells, paste(dir_out, "wells_eu.txt", sep = ""),
    row.names = FALSE, sep = "\t")
save.image(paste(dir_out, "eu.RData", sep = ""))
