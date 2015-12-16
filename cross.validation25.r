# title         : cross.validation25.R
# purpose       : cross validation25, call by cv.n25.r;
# reference     :
# producer      : Prepared by W. Shangguan 
# address       : In Beijing.
# inputs        : points with overlayed covariates 
# outputs       : compressed RDA files ;
# remarks 1     : Takes ca 1.2 hrs to run with 5 cpus for randomforest in the defualt setting in use 

##Load regression matrices:
#suba.sp <- read.csv("ov.BDR_SoilGrids250m.csv")
suba.sp <- read.csv("ov.BDR_SoilGrids250m.csv", col.names=names(read.csv("ov.BDR_SoilGrids250m.csv",nrows=1)), nrows=1000, skip=1000000)

if(soil.flag == 0) 
{
    
    suba.sp <- subset(suba.sp, suba.sp$SOURCEDB=="Wells" |suba.sp$SOURCEDB=="Simulated")     
}
if(arti.flag == 0) suba.sp <- subset(suba.sp, suba.sp$SOURCEDB=="Wells")

## Target vars:
tvar.lst <- c("BDTICM", "BDRICM" , "BDRLOG")

## Summary table:
tbl <- data.frame(ATTRIBUTE_LABEL=tvar.lst, PR_DEPTHS=rep(1, length(tvar.lst)))
tbl$N_OBSERVATIONS <- NA
tbl$OBSERVED_RANGE_MIN <- NA
tbl$OBSERVED_RANGE_MAX <- NA
tbl$rf.SIGMA_EXPLAINED.log <- NA
if(m.flag == 1)
{
    tbl$rk.SIGMA_EXPLAINED.log <- NA
    tbl$dl.SIGMA_EXPLAINED.log <- NA 
    tbl$var.test.H1 <- NA
    tbl$t.test.H1 <- NA
    tbl$var.test.H2 <- NA
    tbl$t.test.H2 <- NA
}
tbl$rf.ME.log <- NA
tbl$rf.RMSE.log <- NA
if(m.flag == 1)
{
    tbl$rk.ME.log <- NA
    tbl$dl.ME.log <- NA
    tbl$rk.RMSE.log <- NA
    tbl$dl.RMSE.log <- NA
}
tbl$rf.ME <- NA
tbl$rf.RMSE <- NA
if(m.flag == 1)
{
    tbl$rk.ME<- NA
    tbl$dl.ME <- NA
    tbl$rk.RMSE <- NA
    tbl$dl.RMSE <- NA
}

## Prepare list to store CV output
cv.lst <- as.list(rep(NA, nrow(tbl)))
names(cv.lst) <- tvar.lst

j<-1
## Run cross-validation and var.test and write the results to a table
## Subset to sub.N samples to speed things up!

for(j in 1:length(tbl$ATTRIBUTE_LABEL)){
#for(j in 1:1){
   tvar <- paste(tbl$ATTRIBUTE_LABEL[j])
   #merge
   rmat <- suba.sp[,c(tvar.lst[1],"LATWGS84","LATWGS84","LONWGS84",as.character(pr.lst))]
  
     set.seed(10002)
   #### subset just as test!!!!!
   rmat <- rmat[sample(1:dim(rmat)[1],200), ]
   sel <- !is.na(rmat[, tvar])
   tbl$N_OBSERVATIONS[j] <- sum(sel)
   if(tvar == "BDRLOG"){
    rv <- range(as.character(rmat[, tvar]))
    fm.g <- as.formula(paste(tvar, ' ~ LATWGS84 + ', paste(pr.lst, collapse="+")))
   } else{
    rv <- quantile(rmat[,tvar], c(.005,.995), na.rm=TRUE)
    fm.g <- as.formula(paste0("log1p(", tvar, ") ~ LATWGS84 + ", paste(pr.lst, collapse="+")))
   }
   tbl$OBSERVED_RANGE_MIN[j] <- rv[1]
   tbl$OBSERVED_RANGE_MAX[j] <- rv[2]
   nfold.x <- kfold(rmat, nfold)
   ####n-fold in parrelal
   #sfInit(parallel=TRUE, cpus=nfold)
   sfInit(parallel=F)
   sfLibrary(h2o)
   sfLibrary(GSIF)
   sfLibrary(randomForestSRC)
   sfLibrary(sp)
   sfExport("PC.flag", "pr.lst", "arti.flag", "soil.flag", "fit.name","m.flag" )
   cv.tvar <- sfClusterApplyLB(1 : nfold, cv_nfold, rmat=rmat,
    nfold.x=nfold.x, fm.g=fm.g, tvar=tvar, sub.N=sub.N, cell.size=cell.size , n=n)
   sfStop()
   gc()
   cv.lst[[j]] <- do.call(rbind, cv.tvar)
}
rm(suba.sp)

for(j in 1:length(tbl$ATTRIBUTE_LABEL))
{
    
  tvar <- paste(tbl$ATTRIBUTE_LABEL[j])
    if(tvar %in% c("BDRICM", "SAPICM2")){
        ## derive ME & RMSE & r2 on log prediction:
        try( tbl$rf.ME.log[j] <- round( mean(cv.lst[[j]]$rf.pred-cv.lst[[j]]$meas, na.rm=TRUE), 3) )
        try( tbl$rf.RMSE.log[j] <- signif( sqrt(mean((cv.lst[[j]]$rf.pred-cv.lst[[j]]$meas)^2, na.rm=TRUE)), 3) )    
        try( tbl$rf.SIGMA_EXPLAINED.log[j] <- round( (1-(var(cv.lst[[j]]$rf.pred-cv.lst[[j]]$meas, na.rm=TRUE)/var(cv.lst[[j]]$meas, na.rm=TRUE)))*100, 1) )
        ## derive ME & RMSE & r2:        
        try( tbl$rf.ME[j] <- round( mean(expm1(cv.lst[[j]]$rf.pred)-expm1(cv.lst[[j]]$meas), na.rm=TRUE), 3) )
        try( tbl$rf.RMSE[j] <- signif( sqrt(mean((expm1(cv.lst[[j]]$rf.pred)-expm1(cv.lst[[j]]$meas))^2, na.rm=TRUE)), 3) )
        try( tbl$rf.SIGMA_EXPLAINED[j] <- round( (1-(var(expm1(cv.lst[[j]]$rf.pred)-expm1(cv.lst[[j]]$meas), na.rm=TRUE)/var(expm1(cv.lst[[j]]$meas), na.rm=TRUE)))*100, 1) )       
        if(m.flag==1){
        ## derive ME & RMSE & r2 on log prediction:
        try( tbl$rk.ME.log[j] <- round( mean(cv.lst[[j]]$rk.pred-cv.lst[[j]]$meas, na.rm=TRUE), 3) )
        try( tbl$rk.RMSE.log[j] <- signif( sqrt(mean((cv.lst[[j]]$rk.pred-cv.lst[[j]]$meas)^2, na.rm=TRUE)), 3) )
        try( tbl$dl.ME.log[j] <- round( mean(cv.lst[[j]]$gbm.pred-cv.lst[[j]]$meas, na.rm=TRUE), 3) )
        try( tbl$dl.RMSE.log[j] <- signif( sqrt(mean((cv.lst[[j]]$dl.pred-cv.lst[[j]]$meas)^2, na.rm=TRUE)), 3) )
        try( tbl$rk.SIGMA_EXPLAINED.log[j] <- round( (1-(var((cv.lst[[j]]$rk.pred)-(cv.lst[[j]]$meas), na.rm=TRUE)/var((cv.lst[[j]]$meas), na.rm=TRUE)))*100, 1) )
        try( tbl$dl.SIGMA_EXPLAINED.log[j] <- round( (1-(var((cv.lst[[j]]$dl.pred)-(cv.lst[[j]]$meas), na.rm=TRUE)/var((cv.lst[[j]]$meas), na.rm=TRUE)))*100, 1) )
        ## derive ME & RMSE & r2:
        try( tbl$rk.ME[j] <- round( mean(expm1(cv.lst[[j]]$rk.pred)-expm1(cv.lst[[j]]$meas), na.rm=TRUE), 3) )
        try( tbl$rk.RMSE[j] <- signif( sqrt(mean((expm1(cv.lst[[j]]$rk.pred)-expm1(cv.lst[[j]]$meas))^2, na.rm=TRUE)), 3) )
        try( tbl$dl.ME[j] <- round( mean(expm1(cv.lst[[j]]$dl.pred)-expm1(cv.lst[[j]]$meas), na.rm=TRUE), 3) )
        try( tbl$dl.RMSE[j] <- signif( sqrt(mean((expm1(cv.lst[[j]]$dl.pred)-expm1(cv.lst[[j]]$meas))^2, na.rm=TRUE)), 3) )
        try( tbl$rk.SIGMA_EXPLAINED[j] <- round( (1-(var(expm1(cv.lst[[j]]$rk.pred)-expm1(cv.lst[[j]]$meas), na.rm=TRUE)/var(expm1(cv.lst[[j]]$meas), na.rm=TRUE)))*100, 1) )
        try( tbl$dl.SIGMA_EXPLAINED[j] <- round( (1-(var(expm1(cv.lst[[j]]$dl.pred)-expm1(cv.lst[[j]]$meas), na.rm=TRUE)/var(expm1(cv.lst[[j]]$meas), na.rm=TRUE)))*100, 1) )
         #rk vs rf on log prediction
        try( tbl$var.test.H1[j] <- round(var.test((cv.lst[[j]]$rk.pred)-(cv.lst[[j]]$meas), (cv.lst[[j]]$rf.pred)-(cv.lst[[j]]$meas), alternative="greater")$p.value, 11) )
        try( tbl$t.test.H1[j] <- round(t.test(abs(cv.lst[[j]]$rk.pred-cv.lst[[j]]$meas), abs(cv.lst[[j]]$rf.pred-cv.lst[[j]]$meas),paired=T, alternative="greater")$p.value, 11) )
        #dl vs rf
        try( tbl$var.test.H2[j] <- round(var.test((cv.lst[[j]]$dl.pred)-(cv.lst[[j]]$meas), (cv.lst[[j]]$rf.pred)-(cv.lst[[j]]$meas), alternative="greater")$p.value, 11) )
        try( tbl$t.test.H2[j] <- round(t.test(abs(cv.lst[[j]]$dl.pred-cv.lst[[j]]$meas), abs(cv.lst[[j]]$rf.pred-cv.lst[[j]]$meas),paired=T, alternative="greater")$p.value, 11) )
         }
    }
    if(tvar %in% c("BDRLOG")){
        #get confucitonmatrix
        c_m <- list(NULL)
        c_m[[1]] <- getconfusion(cv.lst[[j]]$meas,cv.lst[[j]]$rf.pred)
        names(c_m) <- "rf"
        if(m.flag == 1){            
            c_m[[2]] <- getconfusion(cv.lst[[j]]$meas,cv.lst[[j]]$svm.pred)
            names(c_m)[2] <- "svm"
        }                   
      }
}

save(cv.lst,tbl,c_m, file = paste0("./cv/cv_p", PC.flag, "_a", arti.flag, "_s", soil.flag, fit.name, ".RData"))
rm(cv.tvar,cv.lst)
