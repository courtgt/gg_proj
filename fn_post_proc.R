ppdat <- function(mod1, mod2, spp.nom, meas){
  # mod1, mod2 = model objects
  # spp.nom    = list of species names (2 less than plotting points)
  
  # Function extracts summary stats for each species and the mean across all for cover (mu) and occupancy (psi)
  # accepts models fit in Stan using Irvine et al OZAB model framework
  nspp <- length(spp.nom)

  s97x <- rstan::summary(mod1, pars = c("mu", "mean_mu", "psi", "mean_psi"), probs=c(0.025, 0.5, 0.975), seed=108)$summary
  s97x <- data.frame(s97x, year=1998, species = rep(c(spp.nom,"Mean"), 2), meas = meas)
  s97x["mean_mu",1:6] <- plogis(as.numeric(s97x["mean_mu",1:6]))
  s97x["mean_psi",1:6] <- plogis(as.numeric(s97x["mean_psi",1:6]))
  s97x$species = factor(s97x$species)
  s97x$species <- relevel(s97x$species, "Mean")
  
  s21x <- rstan::summary(mod2, pars = c("mu", "mean_mu", "psi", "mean_psi"), probs=c(0.025, 0.5, 0.975), seed=108)$summary
  s21x <- data.frame(s21x, year = 2021, species = rep(c(spp.nom,"Mean"), 2), meas = meas)
  s21x["mean_mu",1:6] <- plogis(as.numeric(s21x["mean_mu",1:6]))
  s21x["mean_psi",1:6] <- plogis(as.numeric(s21x["mean_psi",1:6]))
  s21x$species = factor(s21x$species)
  s21x$species <- relevel(s21x$species, "Mean")
  
  out <- rbind(s97x, s21x)
  out$lw <- ifelse(out$species == "Mean", 1, 0.2)

  mx <- tapply(out$X97.5., out$meas, max)
  labs <- data.frame(x = c(1,1), y = mx, meas = names(mx), labs= c("(a)", "(b)"), year = c(1997,1997))
  library(ggplot2)
  library(grid)
  
  oka <- c("#999999", "#E69F00") # okabe ito colours
  iplx <- ggplot(data= out, aes(x = species, y = mean, group = as.factor(year), color= as.factor(year)))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, vjust=0.2),
          strip.background = element_blank(),
          strip.text = element_blank(),
          legend.position = c(0.9, 0.95),
          legend.title = element_blank(),
          panel.border = element_rect(colour = "black",fill = NA))+#,
    #panel.grid = element_blank())+
    scale_color_discrete(type = oka, name = "Survey year")+
    xlab("")+
    ylab("")+
    geom_point(position = position_dodge(0.5), size=2)+
    geom_errorbar(aes(ymin = X2.5., ymax=X97.5.), lwd = out$lw, width=.2, position=position_dodge(.5))+
    geom_text(aes(x=x, y=y, label= labs), data=labs, hjust=0,col=1)+
    facet_grid(meas~., scales="free"); # print(iplx)
  
  gx <- ggplotGrob(iplx)
  yax <- which(gx$layout$name=="ylab-l")
  gx[["grobs"]][[yax]]$children[[1]]$label <- c("Occupancy", "Cover")
  gx[["grobs"]][[yax]]$children[[1]]$y <- grid::unit(seq(0.25, 0.75, length=2), "npc")
  
  # vector of row numbers
  mix97 <- 1:nspp # first mu vector
  mix21 <- (nspp*2+3):(nspp*2+2+nspp) # second mu vector
  ocix97 <- (nspp+2):(nspp*2+1)
  ocix21 <- (nspp*3 + 4):(nspp*4+3)
  
  # top ten species 
  dmed.cov.raw <- out$X50.[mix21] - out$X50.[mix97]; 
  dmed.cov <- dmed.cov.raw/max(dmed.cov.raw)   # difference between 2021 and 1998 cov 
  dmed.occ.raw <- out$X50.[ocix21] - out$X50.[ocix97]; 
  dmed.occ <- dmed.occ.raw/max(dmed.occ.raw)   # difference between 2021 and 1998 occ
  names(dmed.occ.raw) <- names(dmed.cov.raw) <- out$species[1:nspp]
  
  ttixc <- sort(abs(dmed.cov), dec=TRUE)[10] # tt = topten changes abs
  ttixo <- sort(abs(dmed.occ), dec=TRUE)[10] # tt = topten changes abs
  
  ttspc <- as.character(out$species[which(abs(dmed.cov) >= ttixc)])
  ttspo <- as.character(out$species[which(abs(dmed.occ) >= ttixo)])
  ttspc <- c(ttspc, "Mean")
  ttspo <- c(ttspo, "Mean")
  
  out.ttc <- out[which(out$species %in% ttspc), ]
  out.tto <- out[which(out$species %in% ttspo), ]
  
  mxtc <- tapply(out.ttc$X97.5., out.ttc$meas, max)
  labttc <- data.frame(x = c(1,1), y = mxtc, meas = names(mxtc), labs= c("(a)", "(b)"), year = c(1998,1998))
  
  oka <- c("#999999", "#E69F00") # okabe ito colours
  
  out.ttc$species <- factor(out.ttc$species)
  iplttc <- ggplot(data= out.ttc, aes(x = species, y = mean, group = as.factor(year), color= as.factor(year)))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, vjust=0.2),
          strip.background = element_blank(),
          strip.text = element_blank(),
          legend.position = c(0.9, 0.95),
          legend.title = element_blank(),
          panel.border = element_rect(colour = "black",fill = NA))+#,
    #panel.grid = element_blank())+
    scale_color_discrete(type = oka, name = "Survey year")+
    xlab("")+
    ylab("")+
    geom_point(position = position_dodge(0.5), size=2)+
    geom_errorbar(aes(ymin = X2.5., ymax=X97.5.), lwd = out.ttc$lw, width=.2, position=position_dodge(.5))+
    geom_text(aes(x=x, y=y, label= labs), data=labttc, hjust=0,col=1)+
    facet_grid(meas~., scales="free"); # print(iplttc)
  
  gxt <- ggplotGrob(iplttc)
  yaxt <- which(gxt$layout$name=="ylab-l")
  gxt[["grobs"]][[yaxt]]$children[[1]]$label <- c("Occupancy", "Cover")
  gxt[["grobs"]][[yaxt]]$children[[1]]$y <- grid::unit(seq(0.25, 0.75, length=2), "npc")
  # grid.draw(gxt)
  
  mxto <- tapply(out.tto$X97.5., out.tto$meas, max)
  labtto <- data.frame(x = c(1,1), y = mxto, meas = names(mxto), labs= c("(a)", "(b)"), year = c(1998,1998))
  out.tto$species <- factor(out.tto$species)
  
  oka <- c("#999999", "#E69F00") # okabe ito colours
  ipltto <- ggplot(data= out.tto, aes(x = species, y = mean, group = as.factor(year), color= as.factor(year)))+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, vjust=0.2),
          strip.background = element_blank(),
          strip.text = element_blank(),
          legend.position = c(0.9, 0.95),
          legend.title = element_blank(),
          panel.border = element_rect(colour = "black",fill = NA))+#,
    #panel.grid = element_blank())+
    scale_color_discrete(type = oka, name = "Survey year")+
    xlab("")+
    ylab("")+
    geom_point(position = position_dodge(0.5), size=2)+
    geom_errorbar(aes(ymin = X2.5., ymax=X97.5.), lwd = out.tto$lw, width=.2, position=position_dodge(.5))+
    geom_text(aes(x=x, y=y, label= labs), data=labtto, hjust=0,col=1)+
    facet_grid(meas~., scales="free"); # print(ipltto)
  
  gxto <- ggplotGrob(ipltto)
  yaxt <- which(gxto$layout$name=="ylab-l")
  gxto[["grobs"]][[yaxt]]$children[[1]]$label <- c("Occupancy", "Cover")
  gxto[["grobs"]][[yaxt]]$children[[1]]$y <- grid::unit(seq(0.25, 0.75, length=2), "npc")
  # grid.draw(gxto)
  
  #*********************************####
  # calculate winners and losers ####
  
  wlcov_mat <- wlocc_mat <- NULL
  
  # 1. cover 
  dcov <- out$X97.5[mix97] - out$X2.5.[mix21] # difference between hici 1998 and loci 2021; negative values = winners
  dcovl<- out$X2.5[mix97] - out$X97.5.[mix21] # difference between loci 1998 and hici 2021; positive values = losers
  
  wlcov <- c(spp.nom[which(dcov < 0)], spp.nom[which(dcovl > 0)])
  
  if(length(wlcov) > 0) { 
    wcspp98 <- rownames(out)[which(out$species %in% wlcov & out$year == 1998 & out$meas == "mu")]#, c(5, 4,6)]
    cspp98 <- out[which(out$species %in% wlcov & out$year == 1998 & out$meas == "mu"), c(5, 4,6)]
    
    wcspp21 <- rownames(out)[which(out$species %in% wlcov & out$year == 2021 & out$meas == "mu")]#, c(5, 4,6)]
    cspp21 <- out[which(out$species %in% wlcov & out$year == 2021 & out$meas == "mu"), c(5, 4,6)]
    
    wlcov_mat <- matrix(NA, nrow= length(wlcov), ncol = 9)
    colnames(wlcov_mat) <- c("dq50", "dq025", "dq975", "q50.98", "q025.98", "q975.98","q50.21", "q025.21", "q975_21")
    rownames(wlcov_mat) <- wlcov
    
    for(i in 1:length(wcspp98)){
      ex_97c <- rstan::extract(mod1, pars = wcspp98[i])
      ex_21c <- rstan::extract(mod2, pars = wcspp98[i])
      difc <- unlist(ex_21c) - unlist(ex_97c)
      dc <- quantile(difc, prob = c(0.5, 0.025, 0.975))
      
      wlcov_mat[i,] <- unlist(c(round(dc,3), round(cspp98[i,],3), round(cspp21[i, ],3)))
    }
  }
  
  # 2. occupancy 
  docc <- out$X97.5[ocix97] - out$X2.5.[ocix21] # difference between hici 1998 and loci 2021 
  doccl <- out$X2.5[ocix97] - out$X97.5.[ocix21] # difference between hici 1998 and loci 2021 
  
  # test 95% CI overlap
  wlocc <- c(spp.nom[which(docc < 0)] , spp.nom[which(doccl > 0)] )# 

  if(length(wlocc) > 0) { 
    wospp98 <- rownames(out)[which(out$species %in% wlocc & out$year == 1998 & out$meas == "psi")]#, c(5, 4,6)]
    spp98 <- out[which(out$species %in% wlocc & out$year == 1998 & out$meas == "psi"), c(5, 4,6)]
    
    wospp21 <- rownames(out)[which(out$species %in% wlocc & out$year == 2021 & out$meas == "psi")]#, c(5, 4,6)]
    spp21 <- out[which(out$species %in% wlocc & out$year == 2021 & out$meas == "psi"), c(5, 4,6)]
    
    wlocc_mat <- matrix(NA, nrow= length(wlocc), ncol = 9)
    colnames(wlocc_mat) <- c("dq50", "dq025", "dq975", "q50.98", "q025.98", "q975.98","q50.21", "q025.21", "q975_21")
    rownames(wlocc_mat) <- wlocc
    
    for(i in 1:length(wospp98)){
      ex_97o <- rstan::extract(mod1, pars = wospp98[i])
      ex_21o <- rstan::extract(mod2, pars = wospp98[i])
      difo <- unlist(ex_21o) - unlist(ex_97o)
      do <- quantile(difo, prob = c(0.5, 0.025, 0.975))
      
      wlocc_mat[i,] <- unlist(c(round(do,3), round(spp98[i,],3), round(spp21[i, ],3)))
      }
  }
  
  # change in PPD ####
  # 1. cover 
  pp97mu <- rstan::extract(mod1, pars = "mean_mu")
  pp21mu <- rstan::extract(mod2, pars = "mean_mu")
  
  chng_mat <- matrix(NA, nrow= 6, ncol = 3)
  rownames(chng_mat) <- c("cov97","cov21","dcov", "occ97","occ21","docc")
  colnames(chng_mat) <- c("qmed","loci","hici")
  
  # perforb cover 97
  chng_mat[1,] <- quantile(plogis(pp97mu$mean_mu), prob=c(0.5,0.025, 0.975))
  chng_mat[2,] <- quantile(plogis(pp21mu$mean_mu), prob=c(0.5, 0.025, 0.975))
  chng_mat[3,] <- quantile(plogis(pp21mu$mean_mu) - plogis(pp97mu$mean_mu), prob=c(0.5, 0.025,0.975))

  ppcov <- data.frame(surv2= plogis(pp21mu$mean_mu), surv1 = plogis(pp97mu$mean_mu))
  
  covdf <- data.frame(ppd = c(plogis(pp21mu$mean_mu), plogis(pp97mu$mean_mu)),
                    year = c(rep("2021",4000), rep("1998", 4000)),
                    var = "Cover")
  
  hyp_cov <- brms::hypothesis(ppcov, "surv2 > surv1")

  
  # 2. occupancy
  pp97psi <- rstan::extract(mod1, pars = "mean_psi")
  pp21psi <- rstan::extract(mod2, pars = "mean_psi")
  
  chng_mat[4,] <- quantile(plogis(pp97psi$mean_psi), prob=c(0.5, 0.025, 0.975))
  chng_mat[5,] <- quantile(plogis(pp21psi$mean_psi) , prob=c(0.5, 0.025, 0.975))
  chng_mat[6,] <- quantile(plogis(pp21psi$mean_psi) - plogis(pp97psi$mean_psi), prob=c(0.5,0.025, 0.975))

  occdf <- data.frame(ppd = c(plogis(pp21psi$mean_psi), plogis(pp97psi$mean_psi)),
                    year = c(rep("2021",4000), rep("1998", 4000)),
                    var = "Occupancy")
  
  ppocc <- data.frame(surv2 = plogis(pp21psi$mean_psi), 
                       surv1 = plogis(pp97psi$mean_psi))
  
  hyp_occ <- brms::hypothesis(ppocc, "surv2 > surv1")

  df.delt <- rbind(covdf, occdf)
  df.delt$var <- factor(df.delt$var, levels= c("Occupancy", "Cover"), 
                        labels = c("a) Mean occupancy", "b) Mean cover"))
  
  # change data
  ppdmu <-plogis(pp21mu$mean_mu) - plogis(pp97mu$mean_mu)
  ppdocc <- plogis(pp21psi$mean_psi) - plogis(pp97psi$mean_psi)
  
  var <- rep(c("mean_mu","mean_occ"), each=4000)
  
  ppds <- data.frame(var=var, val = c(ppdmu, ppdocc))
  
  dmu <-quantile(plogis(pp21mu$mean_mu) - plogis(pp97mu$mean_mu), prob=c(0.5, 0.025, 0.975))
  docc <- quantile(plogis(pp21psi$mean_psi) - plogis(pp97psi$mean_psi), prob=c(0.5,0.025, 0.975))
  
  dmat <- matrix(c(dmu,docc), nrow=2, byrow=TRUE)
  colnames(dmat) <-   c("median", "loci", "hici")
  dmat_out <- data.frame(var = c("mean_mu","mean_occ"), dmat)

  return(list(gg_all_spp = gx, gg_toptenc = gxt, gg_topteno = gxto, # three ggplot objects all spp, top ten by cover change, and occ change
              dat_tt_cov = out.ttc, dat_tt_occ = out.tto, dat_all = out,  # ... and their data frames: three data sources  
              wl_cov = wlcov_mat, wl_occ = wlocc_mat, # two winner-loser objects
              dmed.cov = dmed.cov.raw, dmed.occ= dmed.occ.raw, # two change in median objects (cov and occ)
              ave_change = round(chng_mat,3), hyp_cov = hyp_cov, hyp_occ = hyp_occ, df.delt = df.delt,
              ppd_dif = ppds))
  }




