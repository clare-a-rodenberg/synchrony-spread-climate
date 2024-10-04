#calculate wavelet mean field (WMF) for each variable
#check arguments in help file if needed
#help(wmf) 
#take the WMF
wmf.spread<-wmf(spread.cln, time) #take the wavelet mean field
wmf.tmean<-wmf(tmean.cln, time) 
wmf.tmin<-wmf(tmin.cln, time) 
wmf.ppt<-wmf(ppt.cln, time) 
wmf.snow_depth<-wmf(snow_depth.cln, time) 
wmf.nao<-wmf(nao.cln, time)
wmf.pdo<-wmf(pdo.cln, time)
wmf.enso<-wmf(enso.cln, time)

#plot WMF for spread rate, creating Figure 3a 
png("results/Fig3a.png")
plotmag(wmf.spread)
abline(h=log2(4),lty=2)
dev.off()
#plot WMF for the four climate variables, creating Figure 4a-d 
png("results/Fig4a.png")
plotmag(wmf.tmean)
abline(h=log2(4),lty=2)
dev.off()
png("results/Fig4b.png")
plotmag(wmf.ppt)
abline(h=log2(4),lty=2)
dev.off()
png("results/Fig4c.png")
plotmag(wmf.tmin)
abline(h=log2(4),lty=2)
dev.off()
png("results/Fig4d.png")
plotmag(wmf.snow_depth)
abline(h=log2(4),lty=2)
dev.off()
#these plots are exploratory
plotmag(wmf.nao)
plotmag(wmf.pdo)
plotmag(wmf.enso)

#calculate wavelet phasor mean field (WPMF) for each variable using the sig. thresh. of P<0.001
#check arguments in help file if needed
#help(wpmf)
wpmf.spread<-wpmf(spread.cln, time, sigmethod = "fft")
wpmf.tmean<-wpmf(tmean.cln, time, sigmethod = "fft")
wpmf.tmin<-wpmf(tmin.cln, time, sigmethod = "fft")
wpmf.ppt<-wpmf(ppt.cln, time, sigmethod = "fft")
wpmf.snow_depth<-wpmf(snow_depth.cln, time, sigmethod = "fft")
wpmf.nao<-wpmf(nao.cln, time, sigmethod = "fft")
wpmf.pdo<-wpmf(pdo.cln, time, sigmethod = "fft")
wpmf.enso<-wpmf(enso.cln, time, sigmethod = "fft")

#plot WPMF for spread rate, creating Figure 3f
png("results/Fig3f.png")
plotmag(wpmf.spread, sigthresh = 0.999)
abline(h=log2(4),lty=2)
dev.off()
#The below plots are for exploratory purposes only. 
plotmag(wpmf.tmean, sigthresh = 0.999)
plotmag(wpmf.tmin, sigthresh = 0.999)
plotmag(wpmf.ppt, sigthresh = 0.999)
plotmag(wpmf.snow_depth, sigthresh = 0.999)
plotmag(wpmf.nao, sigthresh = 0.999)
plotmag(wpmf.pdo, sigthresh = 0.999)
plotmag(wpmf.enso, sigthresh = 0.999)

#spatial coherence
#remind ourselves of the arguments to this function
#help(coh)

#set short and long time-scales
#short timescale band, 2-4 year period lengths. the short timescale band is the same for all ecoregions
bshort<-c(2,4)
#long timescale band, > 4 year periods. the long timescale band varies by ecoregion. 
blong<-c(4,11)

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
png("results/phase_dif_spread-tmean_sup.png")
plotphase(spcoh.spread.tmean)
dev.off()
png("results/phase_dif_spread-tmin_sup.png")
plotphase(spcoh.spread.tmin)
dev.off()
png("results/phase_dif_spread-ppt_sup.png")
plotphase(spcoh.spread.ppt)
dev.off()
png("results/phase_dif_spread-snow_depth_sup.png")
plotphase(spcoh.spread.snow_depth)
dev.off()

##for the analyses supporting this manuscript, lines 124-206 only apply to the MWS and SUP ecoregions. 

#Invoke wavelet Moran theorem (with 'wlm' function) to quantify the percentage of synchrony in spread that can be 
#explained by synchronous,multi-annual climatic fluctuations. Further calculate cross-terms, a diagnostic of an 
#independence assumption of the wavelet Moran theorem.

#put predictor and response data matrices into a list - choose predictor based on spatial coherence results 
#(statistical significance). 
dlist.sup<-list(ppt.cln, spread.cln) 
wlm.sup<-wlm(dlist.sup, time, resp=2, pred=1, norm="powall") #invoke the wavelet Moran theorem
#print(wlm.sup$coefs)

#Determine the average synchrony explained by timescale (short and long), calculate cross-terms
se.sup<-syncexpl(wlm.sup)
#print(se.sup)
se_short.sup<-se.sup[se.sup$timescales>=bshort[1]&se.sup$timescales<=bshort[2],]
round(100*colMeans(se_short.sup[,3:6])/mean(se_short.sup$sync),4)
saveRDS(se_short.sup,file="results/se_short_sup.rds")
se_long.sup<-se.sup[se.sup$timescales>=blong[1]&se.sup$timescales<=blong[2],]
round(100*colMeans(se_long.sup[,3:6])/mean(se_long.sup$sync),4)
saveRDS(se_short.sup,file="results/se_long_sup.rds")

#extract information from the WMFs for spread (created on line 5) and the model prediction (created on line 123) to 
#create a diagram of essentially the same information as the WMFs for observed (Figure 5 d) and predicted (Figure 5 e) 
#synchrony in spread rate except timescale-specifc synchrony is averaged across all years, which represents the mean 
#squared synchrony
wmf.values <- Mod(get_values(wmf.spread))
wmf.timescales <- Mod(get_timescales(wmf.spread))
wmf.values <- colMeans(wmf.values, na.rm = TRUE)
#add information to dataframe
wmf.df <- as.data.frame(cbind(wmf.timescales,wmf.values))
#get predictions from wlm
preds.mws<-predsync(wlm.mws)
pred.values <- Mod(get_values(preds.mws))
pred.values <- colMeans(pred.values, na.rm = TRUE)
#add information to dataframe
wmf.df$pred.values <- pred.values

#visually compare predicted synchrony to actual synchrony in spread rate
#Figure 5 panels a,b
png("results/Fig5a.png")
plotmag(wmf.spread) #plot the WMF for spread rate
abline(h=log2(4),lty=2)
dev.off()
png("results/Fig5b.png")
plotmag(preds.mws) #plot the predicted WMF
abline(h=log2(4),lty=2)
dev.off()

#Figure 5 panel c 
test.data.long <- melt(wmf.df, id="wmf.timescales")
png("results/Fig5c.png")
ggplot(data=test.data.long, aes(x=wmf.timescales, y=value, colour=variable)) +
  geom_line(linewidth=1) +
  labs(x = "Timescale (years)", y = "Time-averaged \n synchrony") +
  theme(axis.title.x.bottom = element_text(size=10), axis.title.y.left = element_text(size=10), legend.title = element_text(size=0)) +
  scale_color_manual(labels = c("Synchrony", "Pred. sync."),values = c("#1E88E5", "#D81B60")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.key = element_blank())+
  geom_vline(xintercept = 4, size = 1)+
  scale_x_continuous(breaks = c(2,4,6,8,10,12))
dev.off()

#run the coherence function between each climate index and the climate variable(s) for which we invoked the 
#wavelet Moran theorem (i.e., climate variables that had significant spatial coherence with spread rate)

#SUP ecoregion - ppt
spcoh.ppt.nao<-coh(ppt.cln, nao.cln, time, norm="powall", sigmethod="fftsurrog12", nrand=2000)
spcoh.ppt.pdo<-coh(ppt.cln, pdo.cln, time, norm="powall", sigmethod="fftsurrog12", nrand=2000)
spcoh.ppt.enso<-coh(ppt.cln, enso.cln, time, norm="powall", sigmethod="fftsurrog12", nrand=2000)
##
spcoh.ppt.nao<-bandtest(spcoh.ppt.nao, bshort)
spcoh.ppt.nao<-bandtest(spcoh.ppt.nao, blong)
plotmag(spcoh.ppt.nao)
plotphase(spcoh.ppt.nao)
##
spcoh.ppt.pdo<-bandtest(spcoh.ppt.pdo, bshort)
spcoh.ppt.pdo<-bandtest(spcoh.ppt.pdo, blong)
plotmag(spcoh.ppt.pdo)
plotphase(spcoh.ppt.pdo)
##
spcoh.ppt.enso<-bandtest(spcoh.ppt.enso, bshort)
spcoh.ppt.enso<-bandtest(spcoh.ppt.enso, blong)
plotmag(spcoh.ppt.enso)
plotphase(spcoh.ppt.enso)