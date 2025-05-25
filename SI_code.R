# R-code to run and post process models

# Workflow
# 1. import and pre-process data to select functional group of interest (analysis demonstrated for perennial forbs, but was the same approach)
# 2. prepare and run OZAB model (code from Irvine et al 2016 DOI: 10.1007/s13253-016-0265-2)
# 3. extract posterior predictions and analyse
# 4. plot outputs

# import data ####
dat <- read.csv(file=file.path(data.dir,'terricks.entire.csv'))
trt <- read.csv(file=file.path(data.dir,'traits2.csv'), row.names=1, header=TRUE)
env <- read.csv(file=file.path(data.dir,'env_data_all.csv'))
xy <- read.csv(file=file.path(data.dir, 'new way points TTNP.csv'))

#**************####
# 1. preprocess ####
# create unique ID for each sample
rownom <- paste(substr(dat[,1],start=1,stop=1),dat[,2],dat[,3],dat[,4],dat[,5],sep=".")
rownames(dat) <- rownom

# break into veg and metadata
mdat <- dat[,c(1:5)]
veg <- dat[,c(6:ncol(dat))]

# subset to 'e' samples only 
e.idx <- which(dat$position == "e")
mdat.e <- mdat[e.idx,]
veg.e <- veg[e.idx, ]

# check missing spp
mis.idx <- which(apply(veg.e, 2, sum)==0)
vege <- veg.e[,-mis.idx]
mdate <- mdat.e

# match trait data
trt.i <- match(colnames(vege), rownames(trt))
trt.sort <- trt[trt.i,]

# index to select perennial forbs
pf.ix <- which(trt.sort$origin == "native" & trt.sort$type=="perennial" &  trt.sort$lifeform == "forb"|
                 trt.sort$origin == "native" & trt.sort$type=="perennial" & trt.sort$lifeform == "geophyte")

# select
pf.nom <- rownames(trt.sort[pf.ix,]) # names
perf.all <- vege[,which(names(vege)%in% pf.nom)] # data subset

# add xy data
mdate$y <- xy[match(mdate$site, xy$site),"lat"]
mdate$x <- xy[match(mdate$site, xy$site),"lon"]

#******************************#####
# select perennial forbs to model ####
occa <- apply(perf.all > 0, 2, sum) # names of all per forbs
o97 <- apply(perf.all[1:55,]>0, 2, sum) # 1998 samples in 1:55, 2021 in 56:110
o21 <- apply(perf.all[56:110,]>0, 2, sum)

# select non-singletons in both surveys (as with exotics)
nspp.ns <- names(which(o97 >= 1 & o21>= 1))

nspp.ns
# [1] "Arthropodium.fimbriatum"           "Arthropodium.minus"                "Asperula.conferta"                
# [4] "Bulbine.bulbosa"                   "Calotis.scabiosifolia"             "Chrysocephalum.apiculatum"        
# [7] "Convolvulus.spp."                  "Cressa.australis"                  "Euphorbia.dallachyana"            
# [10] "Leptorhynchos.squamatus"           "Linum.marginale"                   "Maireana.enchylaenoides"          
# [13] "Maireana.excavata"                 "Maireana.pentagona"                "Oxalis.perennans"                 
# [16] "Ptilotus.macrocephalus"            "Ptilotus.semilanatus"              "Ptilotus.spathulatus"             
# [19] "Swainsona.plagiotropis"            "Swainsona.procumbens"              "Teucrium.racemosum"               
# [22] "Vittadinia.cuneata.s.l."           "Wurmbea.latifolia.subsp..vanessae"
 
# create spp codes
pf.spp <- c("Art_fim", "Art_min", "Asp_con", "Bul_bul", 
            "Cal_sca", "Chr_api", "Con_spp", "Cre_aus", "Eup_dal",
            "Lep_squ", "Lin_mar", "Mai_enc", "Mai_exc", "Mai_pen", "Oxa_per", 
            "Pti_mac", "Pti_sem", "Pti_spa", "Swa_pla", "Swa_pro", "Teu_rac", 
            "Vit_cun", "Wur_lat")

# collate
natpf.nom <- colnames(veg)

pf.df <- data.frame(site=mdate$site, community= mdate$community, paddock=mdate$paddock,
                     year = mdate$year, gx = mdate$x, gy=mdate$y,
                     vege[,which(colnames(vege)%in% nspp.ns)])

# need in long form
library(reshape2)
dat.pf <- melt(pf.df, id.var= c( "site","community","paddock","year","gx","gy"), variable="species", value.name ="braun")
head(dat.pf)

# recode braun to integers

lutab <- data.frame(bb = c(0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0), 
                    coverclass = c(0, 1, 2, 3, 4, 5, 6))
# 0 = 0.0       = 0
# r = 0.5 = 0-2 = 1
# + = 0.5 = 0-2 = 1
# 1 = 2-5 =     = 2
# 2 = 5-25      = 3
# 3 = 25-50     = 4
# 4 = 50-75     = 5
# 5 = 75-100    = 6

dat.pf$cover <- lutab$coverclass[match(dat.pf$braun, lutab$bb)]

#add detection indicator variable
dat.pf$detect<- ifelse(dat.pf$braun>0,1,0)

#need to switch RE for allotment to species group-level instead
dat.pf$allot <- as.numeric(as.factor(dat.pf$species))

table(dat.pf$species,dat.pf$allot)

# group id for detection/non-det
grp.id <- dat.pf$allot

# prep code for OZAB model (per Irvine et al 2016)
#*************************####
# 2. data prep for OZAB in stan ####
# 1997 
dat.pf97 <- dat.pf[which(dat.pf$year == 1997),]

# group id for detection/non-det
grp.id.97 <- dat.pf97$allot

nspp <- length(unique(dat.pf$species))
# number of groups
J <- nspp # 35
N1.97 <- length(dat.pf97$detect)

#detection/non-detection
z.97 <-dat.pf97$detect

#number of plots with a detection
N2_new.97 <- sum(z.97)

#corresponding ordinal value for first visit to each plot
y_int2.97 <-dat.pf97$cover

#group level id for only plots with non-zeros
grp.id2.97 <- dat.pf97$allot[which(y_int2.97 > 0)]

#removing the zeros from the ordinal variable
y_int2.97 <- y_int2.97[which(y_int2.97 > 0)]

#number of categories
K <- 6
#cutpoints on % cover scale
cuts <- c(0.02, 0.05, 0.25, 0.50, 0.75)

# Stan likes the data to be in the correct order...
pf.spp.data97 <- list("N1" = N1.97,  # number of detections
                     "z" = as.integer(z.97),    # occupancy
                     "N2" = as.integer(N2_new.97),  # number of occupied plots
                     "K" = as.integer(K),        # number of cover classes
                     "cuts" = cuts,  # thresholds 
                     "y_int" = y_int2.97, # non-zero cover scores
                     "J" = as.integer(J),          # number of groups (here species)
                     "grp_id" = as.integer(grp.id.97),  # group for spp detection
                     "grp_id2" = as.integer(grp.id2.97)) # group for non-zero cover 

#*****************************
# 2021 
dat.pf21 <- dat.pf[which(dat.pf$year == 2021),]

# group id for detection/non-det
grp.id.21 <- dat.pf21$allot

# number of groups
J <- nspp
N1.21 <- length(dat.pf21$detect)

#detection/non-detection
z.21 <-dat.pf21$detect

#number of plots with a detection
N2_new.21 <- sum(z.21)

#corresponding ordinal value for first visit to each plot
y_int2.21 <-dat.pf21$cover

#group level id for only plots with non-zeros
grp.id2.21 <- dat.pf21$allot[which(y_int2.21 > 0)]

#removing the zeros from the ordinal variable
y_int2.21 <- y_int2.21[which(y_int2.21 > 0)]

pf.spp.data21 <- list("N1" = N1.21,  # number of detections
                      "z" = as.integer(z.21),    # occupancy
                      "N2" = as.integer(N2_new.21),  # number of occupied plots
                      "K" = as.integer(K),        # number of cover classes
                      "cuts" = cuts,  # thresholds 
                      "y_int" = y_int2.21, # non-zero cover scores
                      "J" = as.integer(J),          # number of groups (here species)
                      "grp_id" = as.integer(grp.id.21),  # group for spp detection
                      "grp_id2" = as.integer(grp.id2.21)) # group for non-zero cover 


#**************************#####
# modelling 

library(rstan)
##Set rstan to work on multiple cores
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#load Stan models from local folder
mod_ozab_mspp <- stan_model('ozab_multispp.stan')

# 1997
pf.spp.fit97 <- sampling(mod_ozab_mspp, pf.spp.data97, chains = 4, iter = 2000, 
                         control = list(adapt_delta=0.999, max_treedepth = 20), seed=108)

# save(pf.spp.fit97, file=file.path(model.dir, "pf.spp.fit97_geqs"))

# 2021 
pf.spp.fit21 <- sampling(mod_ozab_mspp, pf.spp.data21, control = list(adapt_delta=0.999, max_treedepth=20),
                         chains = 4, iter = 2000, seed = 108)
# save(pf.spp.fit21, file=file.path(model.dir, "pf.spp.fit21_geqs"))

#**********************#####
# change in posterior ####
# 1. mean cover in native per forbs
# load(file=file.path(model.dir, "pf.spp.fit21_geqs"))
# load(file=file.path(model.dir, "pf.spp.fit97_geqs"))
# 
pp97mu <- rstan::extract(pf.spp.fit97, pars = "mean_mu")
pp21mu <- rstan::extract(pf.spp.fit21, pars = "mean_mu")

# perforb cover 97
quantile(plogis(pp97mu$mean_mu), prob=c(0.025, 0.5, 0.975))
#  2.5%          50%        97.5%
# 0.0003911702 0.0060655262 0.0181739719 

quantile(plogis(pp21mu$mean_mu), prob=c(0.025, 0.5, 0.975))
#       2.5%         50%       97.5%
#    0.00729718 0.01796206 0.03315902 

# difference in posterior distributions - basis for inference
quantile(plogis(pp21mu$mean_mu) - plogis(pp97mu$mean_mu), prob=c(0.025, 0.5, 0.975))
#      2.5%          50%        97.5%
# -0.003928173  0.011607711  0.027954170 

# wrap up
pdperf <- data.frame(surv2= boot::inv.logit(pp21mu$mean_mu), surv1 = boot::inv.logit(pp97mu$mean_mu))

datmu.pf <- data.frame(ppd = c(surv2= boot::inv.logit(pp21mu$mean_mu),
                          surv1= boot::inv.logit(pp97mu$mean_mu)),
                  year = c(rep("2021",4000), rep("1998", 4000)),
                  var = "Cover")

brms::hypothesis(pdperf, "surv2 > surv1")
#             Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob Star
# 1 (surv2)-(surv1) > 0     0.01      0.01        0     0.02       12.7      0.93     

# numerical increase in median cover
0.011607711/0.0060655262  # 1.913719 = 191%

# 2. occupancy
pp97psi <- rstan::extract(pf.spp.fit97, pars = "mean_psi")
pp21psi <- rstan::extract(pf.spp.fit21, pars = "mean_psi")

hist(boot::inv.logit(pp97psi$mean_psi))
hist(boot::inv.logit(pp21psi$mean_psi), col=2,add=TRUE)

quantile(boot::inv.logit(pp97psi$mean_psi), prob=c(0.025, 0.5, 0.975))
#     2.5%        50%      97.5%
# 0.05941886 0.10230659 0.16772505 

quantile(boot::inv.logit(pp21psi$mean_psi) , prob=c(0.025, 0.5, 0.975))
#    2.5%       50%     97.5%
# 0.07564207 0.12579331 0.19469436 

quantile(boot::inv.logit(pp21psi$mean_psi) - boot::inv.logit(pp97psi$mean_psi), prob=c(0.025, 0.5, 0.975))
#         2.5%         50%       97.5%
#    -0.06004849  0.02383627  0.10471682 

pfoc <- data.frame(surv2 = boot::inv.logit(pp21psi$mean_psi), surv1 = boot::inv.logit(pp97psi$mean_psi))

datocc.pf <- data.frame(ppd = c(boot::inv.logit(pp21psi$mean_psi), boot::inv.logit(pp97psi$mean_psi)),
                  year = c(rep("2021",4000), rep("1998", 4000)),
                  var = "Occupancy")

brms::hypothesis(pfoc, "surv2 > surv1")
# Hypothesis Tests for class :
#             Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob Star
# 1 (surv2)-(surv1) > 0     0.02      0.04    -0.04     0.09       2.61      0.72     

# 72% probability occupancy has increased
0.02383627/0.10230659 # 23% increase in occupancy

# *********************************************************************************
# package change data
ppdmu <-plogis(pp21mu$mean_mu) - plogis(pp97mu$mean_mu)
ppdocc <- plogis(pp21psi$mean_psi) - plogis(pp97psi$mean_psi)

var <- rep(c("mu","occ"), each=4000)

pfnat.ppd <- data.frame(grp = "pfnat", var=var, val = c(ppdmu, ppdocc))
# save(pfnat.ppd, file=file.path(data.dir, "pfnat_ppd_ns.RData"))

dmu <-quantile(plogis(pp21mu$mean_mu) - plogis(pp97mu$mean_mu), prob=c(0.025, 0.5, 0.975))

docc <- quantile(boot::inv.logit(pp21psi$mean_psi) - boot::inv.logit(pp97psi$mean_psi), prob=c(0.025, 0.5, 0.975))

dmat <- matrix(c(dmu,docc), nrow=2, byrow=TRUE)
colnames(dmat) <-   c("loci"  ,"median" ,"hici")
d_pfnat <- data.frame(grp = rep("natperforb",2), var = c("mu","occ"), dmat)
# save(d_pfnat,file=file.path(data.dir,"d_pfnat_ns.RData"))


#************************####
# post process ##### 

source("fn_post_proc.R")
# custom function to post process data and return an object including plots and hypothesis tests

nspp = length(pf.spp)
measx = c(rep("mu", nspp + 1),rep("psi", nspp + 1)) # to add the mean values 


gx <- ppdat(mod1 = pf.spp.fit97, mod2 = pf.spp.fit21, spp.nom = pf.spp, meas= measx)

pppf <- gx

# names(gx)

grid.newpage(); grid.draw(gx[[3]]) # top ten occupancy changing species
grid.newpage(); grid.draw(gx[[2]]) # top ten cover changing species
grid.newpage(); grid.draw(gx[[1]]) # all perennial forb species changes

# plot mean density  ####
library(ggplot2)
library(grid)

df.delt <- gx[[14]] # difference in posteriors for ddensity plots

df.delt$var <- factor(df.delt$var, levels= c("Occupancy", "Cover"),
                      labels = c("a) Mean occupancy", "b) Mean cover"))
# 1. mean overall 
delt.pf <- ggplot(df.delt, aes(x=ppd, fill= year, group = year)) +
  geom_density(alpha=0.4)+
  theme_classic()+
  theme(strip.background = element_blank(),
        #strip.text = element_blank(),
        strip.text.x = element_text(hjust=0),
        legend.position = c(0.8, 0.2),
        #legend.title=element_blank(),
        #panel.border = element_rect(colour = "black",fill = NA),
        panel.grid = element_blank())+
  labs(x= "Prevalence", y = "Predicted posterior density", fill = "Survey\nyear")+
  facet_wrap(var~., nrow=2, scales="free")+
  scale_x_continuous(expand = c(0,0))+
  scale_fill_manual(values=c("#999999", "#E69F00")); print(delt.pf)

names(gx)
gx$ave_change
#        qmed   loci  hici
# cov97 0.006  0.000 0.018
# cov21 0.018  0.007 0.033
# dcov  0.012 -0.004 0.028 # change in cover over 
# occ97 0.102  0.059 0.168 
# occ21 0.126  0.076 0.195
# docc  0.024 -0.060 0.105 # change in occupancy

gx$hyp_cov
#            Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob Star
# 1 (surv2)-(surv1) > 0     0.01      0.01        0     0.02       12.7      0.93     

gx$hyp_occ
#            Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob Star
# 1 (surv2)-(surv1) > 0     0.02      0.04    -0.04     0.09       2.61      0.72 

gx$wl_cov # no individual species had non-overlapping CI
# NULL
gx$wl_occ # two per forbs increased occ 
#          dq50 dq025 dq975 q50.98 q025.98 q975.98 q50.21 q025.21 q975_21
# Eup_dal 0.345 0.217 0.481  0.031   0.006   0.089  0.379   0.257   0.510
# Vit_cun 0.318 0.181 0.458  0.072   0.028   0.150  0.396   0.275   0.525
