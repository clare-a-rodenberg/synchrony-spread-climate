#Create compiled datasets for each climate index by season and ecoregion. 
#Climate Indices are as follows: Pacific Decadal Oscillation (PDO), North Atlantic Oscillation (NAO),
#and the El Nino Southern Oscillation - ENSO (MEI).

#This code only compiles a climate index dataset for the ecoregions where we found that synchrony in spread was related to synchrony 
#in a climate variable. Those ecoregions were the Mixed Wood Shield (MWS), Mixed Wood Plains (MWP), and the Southeastern USA Plains (SUP). 
#For these ecoregions, this code only compiles the climate index data for the seasons that corresponded to the climate 
#variable that was related to synchrony in spread.
#SUP = spring, MWP = spring, MWS = winter#

#Note that this script MUST be run prior to running the DataClean and DataAnalyses codes for these three ecoregions. 

# read in the csv of interest
pdo <- read.csv(here("data", "ERSST PDO Index.csv"), header = TRUE, sep = ",")
enso <- read.csv(here("data", "MEI ENSO.csv"), header = TRUE, sep = ",")
nao <- read.csv(here("data", "NAO.csv"), header = TRUE, sep = ",")

enso.long <- melt(enso[,1:13],id.vars=c("year"),variable.name="MONTHS",value.name="MEI")
enso.long$Index <- rep(1:12, each=45)

#Per the methodology in the manuscript, spring months vary by ecoregion. Winter months are the 
#same, regardless of ecoregion (December through February)

###Sotheastern USA Plains (SUP)###
# # subset years to match number of years used in wavelet analysis
pdo.sup <- pdo[-c(1:139,168:170),]

# SPRING = march through june (SUP)
# PDO,spring #
# subset months
pdo.spring.sup <- pdo.sup[c(1,4:7)]
# calculate mean index for each year
pdo.spring.sup$mean <- rowMeans(pdo.spring.sup[,2:5])
pdo.spring.sup <- pdo.spring.sup[,c("Year", "mean")]
names(pdo.spring.sup)[names(pdo.spring.sup) == 'mean'] <- 'pdo.spring.mean'

# NAO,spring #
# subset years and months
nao.sup <- nao[-c(1:505,853:878),]
nao.sup <- reshape(nao.sup, idvar = "year", timevar = "month", direction = "wide")
nao.spring.sup <- nao.sup[-c(1),c(1,3:6)]
# calculate mean index for each year
nao.spring.sup$mean <- rowMeans(nao.spring.sup[,2:5])
nao.spring.sup <- nao.spring.sup[,c("year", "mean")]
names(nao.spring.sup)[names(nao.spring.sup) == 'mean'] <- 'nao.spring.mean'

# ENSO,spring #
# subset years #
enso.spring.sup <- rep(NA, length(1993:2020))
years<-unique(enso.long$year[enso.long$year%in%c(1993:2020)])
# subset months #
for(year in years){
  tmp<-enso.long[enso.long$year==year & enso.long$Index%in%c(3:6),]
  enso.spring.sup[years %in% year]<-mean(tmp$MEI,na.rm=T)
}
enso.spring.sup<-as.data.frame(enso.spring.sup)
enso.spring.sup$year <- 1993:2020
enso.spring.sup <- enso.spring.sup[,c(2,1)]
names(enso.spring.sup)[names(enso.spring.sup) == 'enso.spring.sup'] <- 'enso.spring.mean'

# write csv with all climate indices for spring
nao.pdo.enso.sup.spring <- cbind(nao.spring.sup, pdo.spring.sup, enso.spring.sup)
nao.pdo.enso.sup.spring <- nao.pdo.enso.sup.spring[,-c(3,5)]
write.csv(nao.pdo.enso.sup.spring, here("data", "clim.indices.sup.spring.csv"), row.names = FALSE)

###Mixed Wood Shield (MWS)### 
# 2000-2020 (start time series one year earlier to account for shift of Dec data)

# PDO,winter #
# # subset years to match number of years used in wavelet analysis
pdo.mws <- pdo[-c(1:145,168:170),]

# subset for winter months (Dec, Jan, Feb)
pdo.winter.mws <- pdo.mws[c(1:3,13)]
# shift Dec column down one so that year t-1 Dec aligns with year t 
pdo.winter.mws['Dec'] <- c(-1.92, head(pdo.winter.mws['Dec'], dim(pdo.winter.mws)[1]-1)[[1]])
# remove the first year that now remains in dataframe (was used to do the shift above)
pdo.winter.mws <- pdo.winter.mws[-c(1),]
pdo.winter.mws$mean <- rowMeans(pdo.winter.mws[,2:4])
pdo.winter.mws <- pdo.winter.mws[,c("Year", "mean")]
names(pdo.winter.mws)[names(pdo.winter.mws) == 'mean'] <- 'pdo.winter.mean'

# NAO,winter #
# subset years and months
nao.mws <- nao[-c(1:588,853:878),]
nao.mws <- reshape(nao.mws, idvar = "year", timevar = "month", direction = "wide")

nao.winter.mws <- nao.mws[c(1:3,13)]
nao.winter.mws['index.12'] <- c(1.6100, head(nao.winter.mws['index.12'], dim(nao.winter.mws)[1]-1)[[1]])
nao.winter.mws <- nao.winter.mws[-c(1),]
nao.winter.mws$mean <- rowMeans(nao.winter.mws[,2:4])
nao.winter.mws <- nao.winter.mws[,c("year", "mean")]
names(nao.winter.mws)[names(nao.winter.mws) == 'mean'] <- 'nao.winter.mean'

# ENSO,winter #
# subset years
enso.winter.mws <- rep(NA, length(2000:2020))
years<-unique(enso.long$year[enso.long$year%in%c(2000:2020)])
# subset months
for(year in years){
  tmp<-enso.long[enso.long$year==year & enso.long$Index<=2,]
  tmp<-rbind(tmp, enso.long[enso.long$year==(year-1) & enso.long$Index>=12,])
  enso.winter.mws[years %in% year]<-mean(tmp$MEI,na.rm=T)
}
enso.winter.mws<-as.data.frame(enso.winter.mws)
enso.winter.mws$year <- 2000:2020
enso.winter.mws <- enso.winter.mws[,c(2,1)]
names(enso.winter.mws)[names(enso.winter.mws) == 'enso.winter.mws'] <- 'enso.winter.mean'

# write csv with all climate indices for winter
nao.pdo.enso.mws.winter <- cbind(nao.winter.mws, pdo.winter.mws, enso.winter.mws)
nao.pdo.enso.mws.winter <- nao.pdo.enso.mws.winter[,-c(3,5)]
write.csv(nao.pdo.enso.mws.winter, here("data", "clim.indices.mws.winter.csv"), row.names = FALSE)

###Mixed Wood Plains (MWP)###
# # subset years of dataframe to match number of years used in wavelet analysis
pdo.mwp <- pdo[-c(1:145,168:170),]

# SPRING = april through july (MWP)
# PDO,spring #
# subset months
pdo.spring.mwp <- pdo.mwp[c(1,5:8)]
# calculate mean index for each year
pdo.spring.mwp$mean <- rowMeans(pdo.spring.mwp[,2:5])
pdo.spring.mwp <- pdo.spring.mwp[,c("Year", "mean")]
names(pdo.spring.mwp)[names(pdo.spring.mwp) == 'mean'] <- 'pdo.spring.mean'

# NAO,spring #
# subset years and months
nao.mwp <- nao[-c(1:588,853:878),]
nao.mwp <- reshape(nao.mwp, idvar = "year", timevar = "month", direction = "wide")
nao.spring.mwp <- nao.mwp[c(1,5:8)]
# calculate mean index for each year
nao.spring.mwp$mean <- rowMeans(nao.spring.mwp[,2:5])
nao.spring.mwp <- nao.spring.mwp[,c("year", "mean")]
names(nao.spring.mwp)[names(nao.spring.mwp) == 'mean'] <- 'nao.spring.mean'

# ENSO,spring #
# subset years
enso.spring.mwp <- rep(NA, length(1999:2020))
years<-unique(enso.long$year[enso.long$year%in%c(1999:2020)])
# subset months
for(year in years){
  tmp<-enso.long[enso.long$year==year & enso.long$Index%in%c(4:7),]
  enso.spring.mwp[years %in% year]<-mean(tmp$MEI,na.rm=T)
}

enso.spring.mwp<-as.data.frame(enso.spring.mwp)
enso.spring.mwp$year <- 1999:2020
enso.spring.mwp <- enso.spring.mwp[,c(2,1)]
names(enso.spring.mwp)[names(enso.spring.mwp) == 'enso.spring.mwp'] <- 'enso.spring.mean'

# write csv with all climate indices for spring
nao.pdo.enso.mwp.spring <- cbind(nao.spring.mwp, pdo.spring.mwp, enso.spring.mwp)
nao.pdo.enso.mwp.spring <- nao.pdo.enso.mwp.spring[,-c(3,5)]
write.csv(nao.pdo.enso.mwp.spring, here("data", "clim.indices.mwp.spring.csv"), row.names = FALSE)

