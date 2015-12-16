## Function to fit models using nfold CV:
cv_nfold <- function(nf, rmat, nfold.x, fm.g, tvar, sub.N, cell.size, n){
   set.seed(10002)
   #nf <-2
   t.rmat <- rmat[!(nfold.x==nf),]
   v.rmat <- rmat[nfold.x==nf,]
   coordinates(t.rmat) <- ~LONWGS84 +LATWGS84.1
   proj4string(t.rmat) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
   #get rid of spatial clustering
   t.rmat <-  as.data.frame(GSIF::sample.grid(t.rmat, cell.size = c(cell.size, cell.size), n = n)$sub)
   ## To speed things up, select only fixed random sub-sample:
   if(sub.N > nrow(t.rmat)){ sub.N <- nrow(t.rmat) }
   sel.s  <- sample(1:nrow(t.rmat), sub.N)
   t.rmat <- t.rmat[sel.s, ]

   options(rf.cores=2, mc.cores=2)
   vs.rfs <- var.select(fm.g, t.rmat, ntree =100, nsplit = 10)
   #save(vs.rfs, file = paste0("./cv/", tvar, "m.rf_", nf, "_p", PC.flag, "_a", arti.flag, "_s", soil.flag,  fit.name, ".Rda"))
   #load(paste0("./cv/", tvar, "m.rf_", nf, "_p", PC.flag, "_a", arti.flag, "_s", soil.flag,  fit.name, ".Rda"))  
   if(tvar == "BDRLOG"){
     fm.g2 <- as.formula(paste0(tvar," ~ ",  paste(vs.rfs$topvars, collapse="+")))
   }else fm.g2 <- as.formula(paste0("log1p(", tvar, ") ~ ",  paste(vs.rfs$topvars, collapse="+")))
   ####Random forests modelling: mtry, 0.8 hour for BDRICM
   localH2O = h2o.init(nthreads = -1)
   dfs <- t.rmat[,all.vars(fm.g2)]
   dfs.hex <- as.h2o(dfs[complete.cases(dfs),], destination_frame = "dfs.hex")
   try( m.rf <- h2o.randomForest(y=1, x=2:length(all.vars(fm.g2)), training_frame=dfs.hex) ) ## TAKES ONLY 2-3 MINS 
   if(!class(.Last.value)[1]=="try-error"&!is.null(m.rf)){
     save(vs.rfs, m.rf, file = paste0("./cv/", tvar, "m.rf_", nf, "_p", PC.flag, "_a", arti.flag, "_s", soil.flag,  fit.name, ".Rda"))
     m.grid <- as.h2o(v.rmat,  destination_frame="m.grid")
     rf.pred <- as.data.frame(h2o.predict(m.rf, m.grid, na.action=na.pass))$predict
   } else {
     rf.pred <- rep(NA, nrow(v.rmat))
   }  
   #load(paste0("./cv/", tvar, "m.rf_", nf, "_p", PC.flag, "_a", arti.flag, "_s", soil.flag, fit.name, ".Rda"))
   #rf.pred <- predict(m.rf, v.rmat)    
    rm(m.rf)
    gc()
 if(m.flag == 1)
 {   
    if(tvar == "BDRLOG")
    {
        rk.pred <- NA 
        dl.pred <- NA 
    }else {
        # Fit linear regression model using stepwise selection
        m.rk <- lm(fm.g, t.rmat)
        m.rk <- step(m.rk)
        save(m.rk, file = paste0("./cv/",tvar, "m.rk_", nf, "_p", PC.flag, "_a", arti.flag, "_s", soil.flag,  fit.name, ".Rda"))
#        load( paste0("./cv/",tvar, "m.rk_", nf, "_p", PC.flag, "_a", arti.flag, "_s", soil.flag,  fit.name, ".Rda"))
        try(rk.pred <- predict(m.rk, v.rmat))  
        rm(m.rk)      
        # Stochastic Gradient Boosting  trees: interacion.depth, n.trees, shrinkage, 2 mins
        ###data should be t.rmat[,c(tvar,pr.lst)], or error happens
        try(m.dl <- h2o.deeplearning(y=1, x=2:length(all.vars(fm.g2)), training_frame=dfs.hex))
        if(!class(.Last.value)[1]=="try-error"&!is.null(m.dl)){
         save(m.dl, file = paste0("./cv/", tvar, "m.dl_", nf, "_p", PC.flag, "_a", arti.flag, "_s", soil.flag,  fit.name, ".Rda"))
         dl.pred <- as.data.frame(h2o.predict(m.dl,m.grid, na.action=na.pass))$predict       
       } else {
         dl.pred <- rep(NA, nrow(v.rmat))
       }
#        load(paste0("./cv/", tvar, "m.dl_", nf, "_p", PC.flag, "_a", arti.flag, "_s", soil.flag, fit.name,  ".Rda"))
#        dl.pred <- predict(m.dl, v.rmat) 
        rm(m.dl)  
        gc()
        }

 } 
 h2o.shutdown(prompt = F)
    out.df <- data.frame(  rf.pred = rf.pred,
    meas=eval(fm.g[[2]], v.rmat), longitude=v.rmat$LONWGS84, latitude=v.rmat$LATWGS84)
    if(m.flag == 1)
    {
        out.df$rk.pred  = rk.pred
        out.df$dl.pred  = dl.pred
    } 
   return(out.df)
}
