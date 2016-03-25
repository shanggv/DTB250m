dir_f <- "E:\\data\\soildata\\depth\\code\\head\\"
source(paste(dir_f, "functions.r", sep = ""))
dir1 <- "E:\\data\\soildata\\depth\\points\\china\\"
c_names  <- c("Source", "Long", "Lat", "D_BR", "D_water", "D_well", "Accu_xy") 
setwd(dir1)
dir_out <- "E:\\data\\soildata\\depth\\points\\profs\\well\\"

lit <- read.table("china_lit.txt", header = T, as.is = T, sep = "\t")[1:5]


lit[,4] <- apply(lit, 1, function(x){
  if(!is.na(x[5]))
  {
    if(x[5]==">" && x[4] <50) return(NA) else return(x[4])
  }else return(x[4])
  
})

lit[,4] <- as.numeric(lit[,4])
tmp <- cbind(lit[1:4], NA,NA,"")
tmp <- tmp[!is.na(tmp[,4]),]
colnames(tmp) <- c_names
print_0(tmp[, 4])
#"122 out of 598 is 0"
hist(tmp[, 4])
hist(log1p(tmp[, 4]))

write.table(tmp, paste(dir_out, "wells_cn.txt", sep = ""), 
            row.names = FALSE, sep = "\t")

library(sp)
coordinates(tmp) <- ~Long+Lat
tmp$tmp<- log1p(tmp$D_BR)
spplot(tmp["tmp"])
