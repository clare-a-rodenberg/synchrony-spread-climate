#Calculate wavelet mean field (WMF) for each variable
#check arguments in help file if needed
# help(wmf) 
#take the WMF
wmf.spread<-wmf(spread.cln, time) 
wmf.tmean<-wmf(tmean.cln, time) 
wmf.tmin<-wmf(tmin.cln, time) 
wmf.ppt<-wmf(ppt.cln, time) 
wmf.snow_depth<-wmf(snow_depth.cln, time) 

#Plot WMF for spread rate, creating Figure 2b 
png(here("results", "figures", "Fig2b.png"))
plotmag(wmf.spread)
segments(x0=1995,x1=2018,y0=log2(4),y1=log2(4),col='white',lty=1,lwd=5)
dev.off()
#Plot WMF for the four climate variables, creating Figure 3b,g,l,q 
png(here("results", "figures", "Fig3b.png"))
plotmag(wmf.tmean)
segments(x0=1995,x1=2018,y0=log2(4),y1=log2(4),col='white',lty=1,lwd=5)
dev.off()
png(here("results", "figures", "Fig3g.png"))
plotmag(wmf.ppt)
segments(x0=1995,x1=2018,y0=log2(4),y1=log2(4),col='white',lty=1,lwd=5)
dev.off()
png(here("results", "figures", "Fig3l.png"))
plotmag(wmf.tmin)
segments(x0=1995,x1=2018,y0=log2(4),y1=log2(4),col='white',lty=1,lwd=5)
dev.off()
png(here("results", "figures", "Fig3q.png"))
plotmag(wmf.snow_depth)
segments(x0=1995,x1=2018,y0=log2(4),y1=log2(4),col='white',lty=1,lwd=5)
dev.off()

#Calculate wavelet phasor mean field (WPMF) for each variable using the sig. thresh. of P<0.001
#check arguments in help file if needed
# help(wpmf)
wpmf.spread<-wpmf(spread.cln, time, sigmethod = "fft")
wpmf.tmean<-wpmf(tmean.cln, time, sigmethod = "fft")
wpmf.tmin<-wpmf(tmin.cln, time, sigmethod = "fft")
wpmf.ppt<-wpmf(ppt.cln, time, sigmethod = "fft")
wpmf.snow_depth<-wpmf(snow_depth.cln, time, sigmethod = "fft")

#Plot WPMF for spread rate, creating Figure 2g
png(here("results", "figures", "Fig2g.png"))
plotmag(wpmf.spread, sigthresh = 0.999)
segments(x0=1995,x1=2018,y0=log2(4),y1=log2(4),col='white',lty=1,lwd=5)
dev.off()

#The below plots are for exploratory purposes only. 
# plotmag(wpmf.tmean, sigthresh = 0.999)
# plotmag(wpmf.tmin, sigthresh = 0.999)
# plotmag(wpmf.ppt, sigthresh = 0.999)
# plotmag(wpmf.snow_depth, sigthresh = 0.999)

#SPATIAL COHERENCE#
#remind ourselves of the arguments to this function
# help(coh)

#Set short and long time-scales
#short timescale band, 2-4 year period lengths. the short timescale band is the same for all ecoregions
bshort<-c(2,4)
#long timescale band, > 4 year periods. the long timescale band varies by ecoregion. 
blong<-c(4,11)

#Run the coherence function beetween spread rate and each climate variable (4) for both short and long timescales
spcoh.spread.tmean<-coh(spread.cln, tmean.cln, time, norm="powall", sigmethod="fftsurrog12", nrand=2000)
spcoh.spread.tmin<-coh(spread.cln, tmin.cln, time, norm="powall", sigmethod="fftsurrog12", nrand=2000)
spcoh.spread.ppt<-coh(spread.cln, ppt.cln, time, norm="powall", sigmethod="fftsurrog12", nrand=2000)
spcoh.spread.snow_depth<-coh(spread.cln, snow_depth.cln, time, norm="powall", sigmethod="fftsurrog12", nrand=2000)

#Assign timescale bands to the spatial coherence output
#short timescale band
spcoh.spread.tmean<-bandtest(spcoh.spread.tmean, bshort)
spcoh.spread.tmin<-bandtest(spcoh.spread.tmin, bshort)
spcoh.spread.ppt<-bandtest(spcoh.spread.ppt, bshort)
spcoh.spread.snow_depth<-bandtest(spcoh.spread.snow_depth, bshort)
#long timescale band
spcoh.spread.tmean<-bandtest(spcoh.spread.tmean, blong)
spcoh.spread.tmin<-bandtest(spcoh.spread.tmin, blong)
spcoh.spread.ppt<-bandtest(spcoh.spread.ppt, blong)
spcoh.spread.snow_depth<-bandtest(spcoh.spread.snow_depth, blong)

#Retrieves p-values and mean phase for spatial coherence - note the climate variables that have significant 
#spatial coherence with spread rate, these climate variables will be used in further analyses
get_bandp(spcoh.spread.tmean)
saveRDS(spcoh.spread.tmean,file=here("results", "spcoh_band", "af", "spread_tmean_af.rds"))
get_bandp(spcoh.spread.tmin)
saveRDS(spcoh.spread.tmin,file=here("results", "spcoh_band","af",  "spread_tmin_af.rds"))
get_bandp(spcoh.spread.ppt)
saveRDS(spcoh.spread.ppt,file=here("results", "spcoh_band", "af", "spread_ppt_af.rds"))
get_bandp(spcoh.spread.snow_depth)
saveRDS(spcoh.spread.snow_depth,file=here("results", "spcoh_band", "af", "spread_snow-depth_af.rds"))

#This is essentially a visual depiction of output from the 'get_bandp' function
# dev.off()
# plotmag(spcoh.spread.tmean)
# plotmag(spcoh.spread.tmin)
# plotmag(spcoh.spread.ppt)
# plotmag(spcoh.spread.snow_depth)

#These plots indicate phase differences - information on the temporal lag between oscillations of two variables 
#For information purposes only. 
png(here("results", "phase_dif", "af",  "spread-tmean_af.png"))
plotphase(spcoh.spread.tmean)
dev.off()
png(here("results", "phase_dif", "af", "spread-tmin_af.png"))
plotphase(spcoh.spread.tmin)
dev.off()
png(here("results", "phase_dif", "af", "spread-ppt_af.png"))
plotphase(spcoh.spread.ppt)
dev.off()
png(here("results", "phase_dif", "af", "spread-snow-depth_af.png"))
plotphase(spcoh.spread.snow_depth)
dev.off()