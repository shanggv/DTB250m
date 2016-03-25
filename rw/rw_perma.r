library(stringr)
library(maptools)
#define global, change if needed
options(stringsAsFactors = FALSE) # not usful for readShapePoints, 
                                # fix s_n == 10
# change dir* if needed
dir_f <- "E:\\data\\soildata\\depth\\points\\codegsifb\\head\\"
source(paste(dir_f, "functions.r", sep = ""))
dir1  <- "E:\\data\\soildata\\depth\\points\\permafrost\\"
c_names  <- c("Source", "Long", "Lat", "D_BR", "D_water", "D_well", "Accu_xy") 
setwd(dir1)
dirs   <- shell("dir * /B", intern = T)
states <- cbind(1:length(dirs), dirs)
colnames(states) <- c("s_code", "states")
states
dir_out <- "E:\\data\\soildata\\depth\\points\\"#      s_code states  

to_deg <- function(t)
{
  
    t <- unlist(str_split(t, " "))
    t <- str_replace(t, pattern = "¡å", replacement = "")
    t <- as.numeric(t[1:3])
    t <- t[1] + t[2]/60 + t[3]/60/60
    t
}
to_deg2 <- function(t)
{
  t <- str_trim(t)
  t <- str_replace(t, pattern = "'", replacement = "")  
  t <- unlist(str_split(t, "   "))
  t <- as.numeric(t)
  t <- t[1] + t[2]/60 
  t
}

#search within the next 35 lines 
search_dat <- function(n, lines, targ)
{
    out <- rep(NA,5)    
    for(j in 1:length(targ))
    #for( j in 1:1)
    {
         i <- n
        while(i < n + 35)
        {
            if(!is.na(str_locate(lines[i], targ[j])[1]))
            {
                tmp <- str_sub(lines[i], str_locate(lines[i], targ[j])[2])
                if(!is.na(str_locate(tmp, "[>0-9]")[1]))
                {
                    out[j] <- str_sub(tmp, str_locate(tmp, "[>0-9]")[1])
                    break
                }else if(!is.na(str_locate(lines[i+1], "[>0-9]")[1]))
                {
                    out[j] <- str_sub(lines[i+1], 
                        str_locate(lines[i+1], "[>0-9]")[1])
                    break
                }
            }        
            i <- i + 1
        }        
    }
    out <- matrix(out, nrow=1)
    colnames(out) <- c("lat", "long", "Depth", "unlithified", "bottomPF")
    out
}

next_rec <- function(n,lines)
{
    i <- n+20
    while(i < n+50)
    {
        if(!is.na(str_locate(lines[i], "Number of the borehole")[1])) break
        i <- i+1
    }
    
    i
}


rec.lst <- as.list(rep(NA, length(dirs)))
names(rec.lst) <- dirs
#-------------------------------Alaska, meter
#0 are kept
s_n <- 1 #number of source
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
tmp <- read.csv(".\\tmp\\perma.txt", sep = "\t")
tmp$Latitude  <- sapply(tmp$Latitude, to_deg)
tmp$Longitude <- - sapply(tmp$Longitude, to_deg)
tmp$pf_depth  <- as.numeric(str_replace(tmp$pf_depth, fixed(" *"), ""))
tmp <- cbind(tmp[, c("Longitude", "Latitude", "pf_depth")], 
    as.numeric(NA), as.numeric(NA), "")
rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
#"0 out of 47 is 0"
del_unzip()

#-------------------------------Canada, meter, not use first
#0 are kept
s_n <- 2 #number of source
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
tmp <- readShapePoints(".\\tmp\\sumtbl2_lam")
tmp <- cbind(tmp@data[, 8:9], as.numeric(as.character(tmp@data$DEPTH_OF_Z)), 
    as.numeric(NA), as.numeric(NA), "")
tmp <- tmp[!is.na(tmp[, 3]), ]
rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
#"1 out of 188 is 0"
del_unzip()

#-------------------------------Northwest Territories and northern Yukon, meter, not use first
#0 are kept


#-------------------------------russia_mong, meter
#0 are kept
s_n <- 4 #number of source
setwd(paste(dir1, dirs[s_n], sep = ""))
do_unzip()
tmp  <- readLines(".\\tmp\\brhcat.dat")
tmp2 <- c("latitude", "longitude", "Depth of the borehole", "unlithified", 
        "bottom of PF")
#c("lat", "long", "Depth", "unlithified", "bottomPF")
i <- 10 
tmp3 <- NULL
while(i < length(tmp))
{

    #search within the next 35 lines 
    tmp3 <- rbind(tmp3, search_dat(i, tmp, tmp2))
    i <- next_rec(i, tmp)

}
tmp3[, 4] <- as.numeric(tmp3[, 4]) 
tmp3 <- tmp3[!is.na(tmp3[, 4]),]
tmp3[, 1]  <- sapply(tmp3[, 1], to_deg2)
tmp3[, 2]  <- sapply(tmp3[, 2], to_deg2)
tmp <- cbind(apply(tmp3[, c(2,1,4)], 2, as.numeric), 
    as.numeric(NA), as.numeric(NA), "")
rec.lst[[s_n]] <- form_rec(tmp, s_n)
print_0(rec.lst[[s_n]][, 4])
#"0 out of 26 is 0"
del_unzip()
tmp<- rec.lst[[s_n]]
write.table(tmp, paste(dir_out, "wells_pr.txt", sep = ""), 
            row.names = FALSE, sep = "\t")


library(sp)
coordinates(tmp) <- ~Long+Lat
tmp$tmp<- log1p(tmp$D_BR)
spplot(tmp["tmp"])


