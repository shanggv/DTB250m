rm(list = ls(all = TRUE))
try(pkgs <- names(sessionInfo()$otherPkgs))
try(pkgs <- paste('package:', pkgs, sep = ""))
try(lapply(pkgs, detach, character.only = TRUE, unload = TRUE))
library(h2o)
library(dismo)
library(GSIF)
library(randomForestSRC)
library(snowfall)

# global define
nfold <- 5 ## Set CV fold
sub.N= 60000
cell.size = 0.1
n = 50
#sub.N= 35000
#cell.size = 0.1
#n = 20
a.dir <- "/home/shang009/big/"# dir of the project
w.dir <- paste0(a.dir, "/worldgrids")
m.dir <- paste0(a.dir, "/soildepth2")
## Cross-validation of models:
setwd(m.dir)
setwd("E:/data/soildata/depth/code/250m")

#source(paste0(a.dir, "soildepth/code/head/cv_nfold25.r"))
source( "cv_nfold25.r")

des <- read.csv("SoilGrids250m_COVS250m.csv")
pr.lst <- des$WORLDGRIDS_CODE
paste(' ~ LATWGS84 + ', paste(pr.lst, collapse="+"))
fit.name <- "all" 
####default global # too few for us
PC.flag <- 0 # 0: not use the PC as predictors; 1: not use
arti.flag <- 1  #1: add artificial points; 0: not
soil.flag <- 0 # 0: without soil profiles; 1: add soil profiles
m.flag <- 1  # 0: with only randomfores; 1: with all models 
source(paste0(a.dir, "soildepth/code/1km/cross.validation25.r"))
