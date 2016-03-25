rm(list = ls(all = TRUE))
dir_f <- "E:\\data\\soildata\\depth\\code\\head\\"
source(paste(dir_f, "functions.r", sep = ""))
dir1 <- "E:\\data\\soildata\\depth\\points\\Southamerica\\brazil"
c_names  <- c("Source", "Long", "Lat", "D_BR", "D_water", "D_well", "Accu_xy") 
setwd(dir1)
dir_out <- "E:\\data\\soildata\\depth\\points\\profs\\well\\"

####siagas
tmp <- read.table("siagas100.txt", sep ="\t", head =T)
tmp[,2] <- tmp[,2]/10000
tmp[,3] <- -tmp[,3]/10000
tmp <- tmp[, c(1, 3, 2, 4)]
tmp <- cbind(tmp, NA,NA,"")
colnames(tmp) <- c_names
print_0(tmp[, 4])
#"293 out of 2004 is 0"
hist(tmp[, 4])
hist(log1p(tmp[, 4]))

write.table(tmp, paste(dir_out, "wells_br.txt", sep = ""), 
            row.names = FALSE, sep = "\t")

library(sp)
coordinates(tmp) <- ~Long+Lat
tmp$tmp<- log1p(tmp$D_BR)
spplot(tmp["tmp"])
