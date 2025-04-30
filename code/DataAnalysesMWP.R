#Calculate wavelet mean field (WMF) for each variable
#check arguments in help file if needed
# help(wmf) 
#take the WMF
wmf.spread<-wmf(spread.cln, time) #take the wavelet mean field
wmf.tmean<-wmf(tmean.cln, time) 
wmf.tmin<-wmf(tmin.cln, time) 
wmf.ppt<-wmf(ppt.cln, time) 
wmf.snow_depth<-wmf(snow_depth.cln, time) 
wmf.nao<-wmf(nao.cln, time)
wmf.pdo<-wmf(pdo.cln, time)
wmf.enso<-wmf(enso.cln, time)

#Plot WMF for spread rate, creating Figure 2d 
png(here("results", "figures", "Fig2d.png"))
plotmag(wmf.spread)
segments(x0=2000,x1=2018,y0=log2(4),y1=log2(4),col='white',lty=1,lwd=5)
dev.off()

#Plot WMF for the four climate variables, creating Figure 4d,i,n,s 
png(here("results", "figures", "Fig4d.png"))
plotmag(wmf.tmean)
segments(x0=2000,x1=2018,y0=log2(4),y1=log2(4),col='white',lty=1,lwd=5)
dev.off()
png(here("results", "figures", "Fig4i.png"))
plotmag(wmf.ppt)
segments(x0=2000,x1=2018,y0=log2(4),y1=log2(4),col='white',lty=1,lwd=5)
dev.off()
png(here("results", "figures", "Fig4n.png"))
plotmag(wmf.tmin)
segments(x0=2000,x1=2018,y0=log2(4),y1=log2(4),col='white',lty=1,lwd=5)
dev.off()
png(here("results", "figures", "Fig4s.png"))
plotmag(wmf.snow_depth)
segments(x0=2000,x1=2018,y0=log2(4),y1=log2(4),col='white',lty=1,lwd=5)
dev.off()
#these plots are exploratory - WMFs of climate indices
# plotmag(wmf.nao)
# plotmag(wmf.pdo)
# plotmag(wmf.enso)

#Calculate wavelet phasor mean field (WPMF) for each variable using the sig. thresh. of P<0.001
#check arguments in help file if needed
# help(wpmf)
wpmf.spread<-wpmf(spread.cln, time, sigmethod = "fft")
wpmf.tmean<-wpmf(tmean.cln, time, sigmethod = "fft")
wpmf.tmin<-wpmf(tmin.cln, time, sigmethod = "fft")
wpmf.ppt<-wpmf(ppt.cln, time, sigmethod = "fft")
wpmf.snow_depth<-wpmf(snow_depth.cln, time, sigmethod = "fft")

#Plot WPMF for spread rate, creating Figure 2i
png(here("results", "figures", "Fig2i.png"))
plotmag(wpmf.spread, sigthresh = 0.999)
segments(x0=2000,x1=2018,y0=log2(4),y1=log2(4),col='white',lty=1,lwd=5)
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
blong<-c(4,8)

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
saveRDS(spcoh.spread.tmean,file=here("results", "spcoh_band", "mwp", "spread_tmean_mwp.rds"))
get_bandp(spcoh.spread.tmin)
saveRDS(spcoh.spread.tmin,file=here("results", "spcoh_band","mwp",  "spread_tmin_mwp.rds"))
get_bandp(spcoh.spread.ppt)
saveRDS(spcoh.spread.ppt,file=here("results", "spcoh_band", "mwp", "spread_ppt_mwp.rds"))
get_bandp(spcoh.spread.snow_depth)
saveRDS(spcoh.spread.snow_depth,file=here("results", "spcoh_band", "mwp", "spread_snow-depth_mwp.rds"))


#This is essentially a visual depiction of output from the 'get_bandp' function
# dev.off()
# plotmag(spcoh.spread.tmean)
# plotmag(spcoh.spread.tmin)
# plotmag(spcoh.spread.ppt)
# plotmag(spcoh.spread.snow_depth)

#These plots indicate phase differences - information on the temporal lag between oscillations of two variables 
#For information purposes only. 
png(here("results", "phase_dif", "mwp", "spread-tmean_mwp.png"))
plotphase(spcoh.spread.tmean)
dev.off()
png(here("results", "phase_dif", "mwp", "spread-tmin_mwp.png"))
plotphase(spcoh.spread.tmin)
dev.off()
png(here("results", "phase_dif", "mwp", "spread-ppt_mwp.png"))
plotphase(spcoh.spread.ppt)
dev.off()
png(here("results", "phase_dif", "mwp", "spread-snow_depth_mwp.png"))
plotphase(spcoh.spread.snow_depth)
dev.off()

#Invoke wavelet Moran theorem (with 'wlm' function) to quantify the percentage of synchrony in spread that can be 
#explained by synchronous,multi-annual climatic fluctuations. Further calculate cross-terms, a diagnostic of an 
#independence assumption of the wavelet Moran theorem.

#Put predictor and response data matrices into a list - choose predictor based on spatial coherence results 
#(statistical significance). 
dlist.mwp<-list(ppt.cln, spread.cln) 
wlm.mwp<-wlm(dlist.mwp, time, resp=2, pred=1, norm="powall") #invoke the wavelet Moran theorem
# print(wlm.mwp$coefs)

#Determine the average synchrony explained by timescale (short and long), calculate cross-terms
se.mwp<-syncexpl(wlm.mwp)
#print(se.mwp)
se_short.mwp<-se.mwp[se.mwp$timescales>=bshort[1]&se.mwp$timescales<=bshort[2],]
round(100*colMeans(se_short.mwp[,3:6])/mean(se_short.mwp$sync),4)
saveRDS(se_short.mwp,file=here("results", "synch_expl", "mwp", "se_short_mwp.rds"))
se_long.mwp<-se.mwp[se.mwp$timescales>=blong[1]&se.mwp$timescales<=blong[2],]
round(100*colMeans(se_long.mwp[,3:6])/mean(se_long.mwp$sync),4)
saveRDS(se_short.mwp,file=here("results", "synch_expl", "mwp", "se_long_mwp.rds"))

#Extract information from the WMFs for spread (created on line 5) and the model prediction (created on line 147-48) to 
#create a visual representation of the WMFs for observed (Figure 4g) and predicted (Figure 4h) synchrony in spread rate 
#except timescale-specifc synchrony is averaged across all years, representing the mean squared synchrony
wmf.values <- Mod(get_values(wmf.spread))
wmf.timescales <- Mod(get_timescales(wmf.spread))
wmf.values <- colMeans(wmf.values, na.rm = TRUE)
#add information to dataframe
wmf.df <- as.data.frame(cbind(wmf.timescales,wmf.values))
#get predictions from wlm
preds.mwp<-predsync(wlm.mwp)
pred.values <- Mod(get_values(preds.mwp))
pred.values <- colMeans(pred.values, na.rm = TRUE)
#add information to dataframe
wmf.df$pred.values <- pred.values

#Visually compare predicted synchrony to actual synchrony in spread rate
#Figure 4 panels d,e
png(here("results", "figures", "Fig4d.png"))
plotmag(wmf.spread) #plot the WMF for spread rate
segments(x0=2000,x1=2018,y0=log2(4),y1=log2(4),col='white',lty=1,lwd=5)
dev.off()
png(here("results", "figures", "Fig4e.png"))
plotmag(preds.mwp) #plot the predicted WMF
segments(x0=2000,x1=2018,y0=log2(4),y1=log2(4),col='white',lty=1,lwd=5)
dev.off()

#Figure 4 panel f 
test.data.long <- melt(wmf.df, id="wmf.timescales")
png(here("results", "figures", "Fig4f.png"))
ggplot(data=test.data.long, aes(x=wmf.timescales, y=value, colour=variable)) +
  geom_line(linewidth=1) +
  labs(x = "Timescale (years)", y = "Time-averaged \n synchrony") +
  theme(axis.title.x.bottom = element_text(size=10), axis.title.y.left = element_text(size=10), legend.title = element_text(size=0)) +
  scale_color_manual(labels = c("Synchrony", "Pred. sync."),values = c("#1E88E5", "#D81B60")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.key = element_blank())+
  geom_vline(xintercept=4, lwd=1) +
  scale_x_continuous(breaks = c(2,4,6,8,10,12))+
  scale_y_continuous(limits = c(0, 1))
dev.off()

#Run the coherence function between each climate index and the climate variable(s) for which we invoked the 
#wavelet Moran theorem (i.e., climate variables that had significant spatial coherence with spread rate)
#Relationahips were not significant for this ecoregion. Outputs for information only. 

#MWP ecoregion - ppt
spcoh.ppt.nao<-coh(ppt.cln, nao.cln, time, norm="powall", sigmethod="fftsurrog12", nrand=2000)
spcoh.ppt.pdo<-coh(ppt.cln, pdo.cln, time, norm="powall", sigmethod="fftsurrog12", nrand=2000)
spcoh.ppt.enso<-coh(ppt.cln, enso.cln, time, norm="powall", sigmethod="fftsurrog12", nrand=2000)

##
spcoh.ppt.nao<-bandtest(spcoh.ppt.nao, bshort)
spcoh.ppt.nao<-bandtest(spcoh.ppt.nao, blong)
spcoh.ppt.nao.band <- get_bandp(spcoh.ppt.nao)
saveRDS(spcoh.ppt.nao.band,file=here("results", "spcoh_band", "mwp", "ppt_nao_mwp.rds"))
# plotmag(spcoh.ppt.nao)
png(here("results", "phase_dif", "mwp", "ppt-nao_mwp.png"))
plotphase(spcoh.ppt.nao)
dev.off()

##
spcoh.ppt.pdo<-bandtest(spcoh.ppt.pdo, bshort)
spcoh.ppt.pdo<-bandtest(spcoh.ppt.pdo, blong)
spcoh.ppt.pdo.band <- get_bandp(spcoh.ppt.pdo)
saveRDS(spcoh.ppt.pdo.band,file=here("results", "spcoh_band", "mwp", "ppt_pdo_mwp.rds"))
# plotmag(spcoh.ppt.pdo)
png(here("results", "phase_dif", "mwp", "ppt-pdo_mwp.png"))
plotphase(spcoh.ppt.pdo)
dev.off()

##
spcoh.ppt.enso<-bandtest(spcoh.ppt.enso, bshort)
spcoh.ppt.enso<-bandtest(spcoh.ppt.enso, blong)
saveRDS(spcoh.ppt.enso.band,file=here("results", "spcoh_band", "mwp", "ppt_enso_mwp.rds"))
# plotmag(spcoh.ppt.enso)
png(here("results", "phase_dif", "mwp", "ppt-enso_mwp.png"))
plotphase(spcoh.ppt.enso)
dev.off()
