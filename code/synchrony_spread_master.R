## Synchrony in spread-climate Master File
## Author: Clare Rodenberg

#load relevant packages
library(wsyn)
library(tidyverse)
library(data.table)
library(tidyr)
library(reshape2)
library(here)

#where am i?
here()

#this package and associated setting speeds up computation time (may not be relevant for smaller datasets)
library(parallel)
numCores <- 6
options(mc.cores = numCores)

#There are two types of source code, DataClean and DataAnalyses, with one of each code for each ecoregion. 
#DataClean produces cleaned time series from input data for use in wavelet procedures (cleaned time series are denoted as
#variable_name.cln)
#DataAnalyses contains code to run all wavelet analyses supporting this manuscript and produces the individual
#panels of each figure

#Run the climate index compilation code. This code creates compiled datasets for each climate index by season 
#and ecoregion. The code only compiles a climate index dataset for the ecoregions where we found that synchrony is spread was related to synchrony in a climate variable. 
#Those ecoregions are the MWS, MWP, and SUP.
system.time(source(here("code", "ClimateIndexManipulation.R")))

#When running the below scripts, only run one ecoregion at a time, clearing the workspace in between ecoregions. 

#Southeastern USA Plains (SUP) ecoregion
system.time(source(here("code", "DataCleanSUP.R")))
system.time(source(here("code", "DataAnalysesSUP.R")))

#Appalachian Forest (AF) ecoregion
system.time(source(here("code", "DataCleanAF.R")))
system.time(source(here("code", "DataAnalysesAF.R")))

#Central USA Plains (CUP) ecoregion
system.time(source(here("code", "DataCleanCUP.R")))
system.time(source(here("code", "DataAnalysesCUP.R")))

#Mixed Wood Plains (MWP) ecoregion
system.time(source(here("code", "DataCleanMWP.R")))
system.time(source(here("code", "DataAnalysesMWP.R")))

#Mixed Wood Shield (MWS) ecoregion
system.time(source(here("code", "DataCleanMWS.R")))
system.time(source(here("code", "DataAnalysesMWS.R")))
