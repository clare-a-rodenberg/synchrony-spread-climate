#calculate wavelet mean field (WMF) for each variable
#check arguments in help file if needed
#help(wmf) 
#take the WMF
wmf.spread<-wmf(spread.cln, time) #take the wavelet mean field
wmf.tmean<-wmf(tmean.cln, time) 
wmf.tmin<-wmf(tmin.cln, time) 
wmf.ppt<-wmf(ppt.cln, time) 
wmf.snow_depth<-wmf(snow_depth.cln, time) 

#plot WMF for spread rate, creating Figure 3d 
png("results/Fig3d.png")
plotmag(wmf.spread)
abline(h=log2(4),lty=2)
dev.off()
#plot WMF for the four climate variables, creating Figure 4m-p 
png("results/Fig4m.png")
plotmag(wmf.tmean)
abline(h=log2(4),lty=2)
dev.off()
png("results/Fig4n.png")
plotmag(wmf.ppt)
abline(h=log2(4),lty=2)
dev.off()
png("results/Fig4o.png")
plotmag(wmf.tmin)
abline(h=log2(4),lty=2)
dev.off()
png("results/Fig4p.png")
plotmag(wmf.snow_depth)
abline(h=log2(4),lty=2)
dev.off()

#calculate wavelet phasor mean field (WPMF) for each variable using the sig. thresh. of P<0.001
#check arguments in help file if needed
#help(wpmf)
wpmf.spread<-wpmf(spread.cln, time, sigmethod = "fft")
wpmf.tmean<-wpmf(tmean.cln, time, sigmethod = "fft")
wpmf.tmin<-wpmf(tmin.cln, time, sigmethod = "fft")
wpmf.ppt<-wpmf(ppt.cln, time, sigmethod = "fft")
wpmf.snow_depth<-wpmf(snow_depth.cln, time, sigmethod = "fft")

#plot WPMF for spread rate, creating Figure 3i
png("results/Fig3i.png")
plotmag(wpmf.spread, sigthresh = 0.999)
abline(h=log2(4),lty=2)
dev.off()
#The below plots are for exploratory purposes only. 
plotmag(wpmf.tmean, sigthresh = 0.999)
plotmag(wpmf.tmin, sigthresh = 0.999)
plotmag(wpmf.ppt, sigthresh = 0.999)
plotmag(wpmf.snow_depth, sigthresh = 0.999)

#spatial coherence
#remind ourselves of the arguments to this function
#help(coh)

#set short and long time-scales
#short timescale band, 2-4 year period lengths. the short timescale band is the same for all ecoregions
bshort<-c(2,4)
#long timescale band, > 4 year periods. the long timescale band varies by ecoregion. 
blong<-c(4,8)

#run the coherence function beetween spread rate and each climate variable (4) for both short and long timescales
spcoh.spread.tmean<-coh(spread.cln, tmean.cln, time, norm="powall", sigmethod="fftsurrog12", nrand=2000)
spcoh.spread.tmin<-coh(spread.cln, tmin.cln, time, norm="powall", sigmethod="fftsurrog12", nrand=2000)
spcoh.spread.ppt<-coh(spread.cln, ppt.cln, time, norm="powall", sigmethod="fftsurrog12", nrand=2000)
spcoh.spread.snow_depth<-coh(spread.cln, snow_depth.cln, time, norm="powall", sigmethod="fftsurrog12", nrand=2000)

#assign timescale bands to the spatial coherence output
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

#retrieves p-values and mean phase for spatial coherence - note the climate variables that have significant 
#spatial coherence with spread rate, these climate variables will be used in further analyses
get_bandp(spcoh.spread.tmean)
get_bandp(spcoh.spread.tmin)
get_bandp(spcoh.spread.ppt)
get_bandp(spcoh.spread.snow_depth)

#this is essentially a visual depiction of output from the 'get_bandp' function
plotmag(spcoh.spread.tmean)
plotmag(spcoh.spread.tmin)
plotmag(spcoh.spread.ppt)
plotmag(spcoh.spread.snow_depth)

#these plots indicate phase differences - information on the temporal lag between oscillations of two variables 
png("results/phase_dif_spread-tmean_mwp.png")
plotphase(spcoh.spread.tmean)
dev.off()
png("results/phase_dif_spread-tmin_mwp.png")
plotphase(spcoh.spread.tmin)
dev.off()
png("results/phase_dif_spread-ppt_mwp.png")
plotphase(spcoh.spread.ppt)
dev.off()
png("results/phase_dif_spread-snow_depth_mwp.png")
plotphase(spcoh.spread.snow_depth)
dev.off()