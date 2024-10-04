#CREATE and CLEAN time series for raw spread rate, climate variable (tmean, tmin, ppt, and snow_depth),
#and climate index data. These time series are the input for the wavelet analyses methods employed in this
#manuscript.

#Read in combined spread rate + associated climate variable data (tmean, tmin, ppt, and snow_depth) for the 
#Southern USA Plains (SUP) ecoregion
in.data.t <- read.csv("data/sup_spread_climate.csv", header = TRUE, sep = ",")

#exclude time series <20 years in length because wavelet-based approaches require relatively long time series
in.data.t <- setDT(in.data.t)[, grp := cumsum(c(0, diff(year_t)) > 1), by = bearing2
][, if (.N > 20) .SD, by = .(bearing2, grp)][, grp := NULL][]

#create separate datasets from the master dataset by variable --> spread rate and each of the 4 climate variables
#with each separate dataset we will subset and reshape the dataframe to create separate time series of each input 
#we want each column to represent a year of data for all locations, with rows representing the 'time series'
#a location

#create a time series of the spread rate data, preserving bearing, year, and spread rate
in.data.spread <- in.data.t[,c("bearing2", "year_t", "spread_rate")]
in.data.spread <- reshape(in.data.spread, idvar = "bearing2", timevar = "year_t", direction = "wide")
#remove column 'bearing' as it is no longer needed
in.data.spread <- in.data.spread[, -c(1)] 
#remove years with too many NA values (i.e., too few locations). what is 'too few' is based on analyst judgement.
#for this manuscript we never retained a year that had fewer than five locations
in.data.spread <- in.data.spread[, -c(1:3)] 
#remove any rows with NA & convert to matrix
in.data.spread <- in.data.spread[complete.cases(in.data.spread),]
in.data.spread <- as.matrix(in.data.spread)

#assign 'time' element to the data based on years available in final in.data.spread matrix 
time <- c(1993:2020) 

#create a time series of the mean temperature data, preserving bearing, year, and mean temperature
in.data.tmean <- in.data.t[,c("bearing2", "year_t", "tmean")]
in.data.tmean <- reshape(in.data.tmean, idvar = "bearing2", timevar = "year_t", direction = "wide")
in.data.tmean <- in.data.tmean[, -c(1)] 
in.data.tmean <- in.data.tmean[, -c(1:3)] 
in.data.tmean <- in.data.tmean[complete.cases(in.data.tmean),]
in.data.tmean <- as.matrix(in.data.tmean)

#create a time series of the minimum temperature data, preserving bearing, year, and minimum temperature
in.data.tmin <- in.data.t[,c("bearing2", "year_t", "tmin")]
in.data.tmin <- reshape(in.data.tmin, idvar = "bearing2", timevar = "year_t", direction = "wide")
in.data.tmin <- in.data.tmin[, -c(1)] 
in.data.tmin <- in.data.tmin[, -c(1:3)] 
in.data.tmin <- in.data.tmin[complete.cases(in.data.tmin),]
in.data.tmin <- as.matrix(in.data.tmin)

#create a time series of the precipitation data, preserving bearing, year, and precipitation
in.data.ppt <- in.data.t[,c("bearing2", "year_t", "ppt")]
in.data.ppt <- reshape(in.data.ppt, idvar = "bearing2", timevar = "year_t", direction = "wide")
in.data.ppt <- in.data.ppt[, -c(1)]
in.data.ppt <- in.data.ppt[, -c(1:3)]
in.data.ppt <- in.data.ppt[complete.cases(in.data.ppt),]
in.data.ppt <- as.matrix(in.data.ppt)

#create a time series of the snow depth data, preserving bearing, year, and snow depth
in.data.snow_depth <- in.data.t[,c("bearing2", "year_t", "snow_depth")]
in.data.snow_depth <- reshape(in.data.snow_depth, idvar = "bearing2", timevar = "year_t", direction = "wide")
in.data.snow_depth <- in.data.snow_depth[, -c(1)] 
in.data.snow_depth <- in.data.snow_depth[, -c(1:3)] 
in.data.snow_depth <- in.data.snow_depth[complete.cases(in.data.snow_depth),]
in.data.snow_depth <- as.matrix(in.data.snow_depth)

#clean all time series - process involves Box-Cox transformation, linear detrend, and standardizing data to 
#a variance of one. look at the help files of function 'cleandat' for further information.
#help(cleandat) 
spread.cln<-cleandat(in.data.spread,time,clev=5)$cdat #prepare data for analysis
#plot(time,spread.cln[1,],type="l") #take a look at 1 time series

tmean.cln<-cleandat(in.data.tmean,time,clev=5)$cdat
#plot(time,tmean.cln[1,],type="l") 
tmin.cln<-cleandat(in.data.tmin,time,clev=5)$cdat
#plot(time,tmin.cln[1,],type="l") 
ppt.cln<-cleandat(in.data.ppt,time,clev=5)$cdat
#plot(time,ppt.cln[1,],type="l") 
snow_depth.cln<-cleandat(in.data.snow_depth,time,clev=5)$cdat
#plot(time,snow_depth.cln[1,],type="l")

##Lines 87-130 are the same data processes as those that were applied to create time series for spread rate and 
#the 4 climate variables, but are for the climate index data. The climate index datasets must be read in and joined
#to the original input dataframe in.data.t  Note that the climate index data are only relevant for those ecoregions for 
#which we identified syncrhony in spread was explained by a climate variable (the MWS and SUP ecoregions) as we were 
#interested in understanding whether teleconnections explained synchrony in the climate variable that drove synchrony in 
#spread. Therefore, **these lines are only relevant to the MWS and SUP ecoregions.**

#read in climate index data
in.data.climate.indices <- read.csv("data/clim.indices.sup.spring.csv", header = TRUE, sep = ",")

#join main input data (in.data.t) to climate index data by YEAR
in.data.joined <- left_join(in.data.t, in.data.climate.indices, by = c("year_t" = "year"), keep = TRUE)
#omit any NA values from the joined dataset
in.data.joined <- na.omit(in.data.joined)
#rename dataset 
in.data.t <- in.data.joined

# PDO time series for SUP
in.data.pdo <- in.data.t[,c("bearing2", "year_t", "pdo.spring.mean")]
in.data.pdo <- reshape(in.data.pdo, idvar = c("bearing2"), timevar = "year_t", direction = "wide")
in.data.pdo <- in.data.pdo[, -c(1)]
# remove any rows with NA & convert to matrix
in.data.pdo <- in.data.pdo[complete.cases(in.data.pdo),]
in.data.pdo <- as.matrix(in.data.pdo)

# ENSO (MEI index) time series for SUP
in.data.enso <- in.data.t[,c("bearing2", "year_t", "enso.spring.mean")]
in.data.enso <- reshape(in.data.enso, idvar = c("bearing2"), timevar = "year_t", direction = "wide")
in.data.enso <- in.data.enso[, -c(1)]
# remove any rows with NA & convert to matrix
in.data.enso <- in.data.enso[complete.cases(in.data.enso),]
in.data.enso <- as.matrix(in.data.enso)

# NAO time series for SUP
in.data.nao <- in.data.t[,c("bearing2", "year_t", "nao.spring.mean")]
in.data.nao <- reshape(in.data.nao, idvar = c("bearing2"), timevar = "year_t", direction = "wide")
in.data.nao <- in.data.nao[, -c(1)]
# remove any rows with NA & convert to matrix
in.data.nao <- in.data.nao[complete.cases(in.data.nao),]
in.data.nao <- as.matrix(in.data.nao)

#clean climate index data
nao.cln<-cleandat(in.data.nao,time,clev=5)$cdat
#plot(time,nao.cln[1,],type="l") 
pdo.cln<-cleandat(in.data.pdo,time,clev=5)$cdat
#plot(time,pdo.cln[1,],type="l") 
enso.cln<-cleandat(in.data.enso,time,clev=5)$cdat
#plot(time,pdo.cln[1,],type="l") 

