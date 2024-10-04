## Synchrony in spread-climate Master File
## Author: Clare Rodenberg

#clear workspace
rm(list=ls())

#load relevant packages
library(wsyn)
library(tidyverse)
library(data.table)
library(tidyr)
library(reshape2)

#this package and associated setting speeds up computation time (may not be relevant for smaller datasets)
library(parallel)
numCores <- 6
options(mc.cores = numCores)

#set directory
work_dir<-("C:/Users/Clare/Desktop/synchrony_spread/code/")
setwd(work_dir)

#There are two types of source code, DataClean and DataAnalyses, with one of each code for each ecoregion. 
#DataClean produces cleaned time series from input data for use in wavelet procedures (cleaned time series are denoted as
#variable_name.cln)
#DataAnalyses contains code to run all wavelet analyses supporting this manuscript and produces the individual
#panels of each figure
#When running the below scripts, only run one ecoregion at a time, clearing the workspace in between ecoregions. 

#Southeastern USA Plains (SUP) ecoregion
system.time(source("DataCleanSUP.R"))
system.time(source("DataAnalysesSUP.R"))

#Appalachian Forest (AF) ecoregion
system.time(source("DataCleanAF.R"))
system.time(source("DataAnalysesAF.R"))

#Central USA Plains (CUP) ecoregion
system.time(source("DataCleanCUP.R"))
system.time(source("DataAnalysesCUP.R"))

#Mixed Wood Plains (MWP) ecoregion
system.time(source("DataCleanMWP.R"))
system.time(source("DataAnalysesMWP.R"))

#Mixed Wood Shield (MWS) ecoregion
system.time(source("DataCleanMWS.R"))
system.time(source("DataAnalysesMWS.R"))