# load packages ####
library(tidyverse)
library(ggplot2)
library(cowplot)
library(glmmTMB)
library(performance)
library(ggeffects)
library(ggpubr)
library(DHARMa)
library(broom.mixed)

#get data ####
setwd("")
seeddat <- read.csv("seed_exp_r.csv")

#data prep ####
#all blanks are true 0s
seeddat[is.na(seeddat)] <- 0

# other stuff ####
### histogram of final counts (counts at time 4 and remove controls)
t4 <- filter(seeddat, treat != "EC", treat != "GC")
hist(t4$t4)

sum <- summarise(t4, t4=sum(t4))
#total count of germinants across whole expriment at time 4 is 201 seedlings 

sum2 <- summarise(t4, t4=sum(density))
#total number of seeds sown is 14,400

# remove controls from df, controls all 0 at t4
# create new colum of germination fraction to remove density as factor
seed <- pivot_longer(seeddat, c(7:10), names_to = "time", values_to = "count") %>% 
  filter(treat != "EC", treat != "GC") %>% 
  mutate(frac = (count/density)*100) %>% 
  mutate(Week = case_when(time == "t1" ~ "5",
                          time == "t2" ~ "9",
                          time == "t3" ~ "15",
                          time == "t4" ~ "21"))

#summarise with germ frac and se
seed.summ <- group_by(seed, Week, microsite, species, treat) %>% 
  summarise(germ_frac = mean(frac), sd = sd(frac), n = n()) %>% 
  mutate(se = sd/sqrt(n))

seed.summ$Week <- as.numeric(seed.summ$Week)

#summary plot of all time points using germ frac
ggplot(seed.summ, aes(Week, germ_frac, shape = microsite))+
  geom_line(aes(Week, germ_frac, group = microsite))+
  geom_pointrange(aes(Week, germ_frac,
                      ymin=germ_frac-se, ymax = germ_frac+se),
                    fill = "white")+
  scale_shape_manual(name = "Microsite", 
                     labels=c("Crust", "Digging", "Litter"),
                     values=c(21, 24, 22))+
  facet_grid(rows = vars(species), cols = vars(treat))+
  geom_vline(data = subset(seed.summ, treat=="G"),
             aes(xintercept = 12),
             colour = "red", linetype = 2)+
  scale_x_continuous(breaks=c(5,9,15, 21))+
  xlab("Weeks post-sowing")+
  ylab("Germination fraction (%)")+
  theme(text = element_text(size = 20))


# models ####
## data prep ####
#remove controls
seed.fin <- filter(seeddat, treat != "EC", treat != "GC")

#set all independent variables as factors - 
#density treatment is technically an interger but it
#doesn't seem right to treat it as continuous where there
#are only two levels (20,60)
seed.fin$density <- as.ordered(seed.fin$density)
seed.fin$treat <- as.factor(seed.fin$treat)
seed.fin$microsite <- as.factor(seed.fin$microsite)
seed.fin$species <- as.factor(seed.fin$species)

#df for modelling without microsite litter
seed.fil <- filter(seed.fin, microsite!="litt")
seed.fil$Crust <- as.numeric(seed.fil$microsite == "crust")
seed.fil$Soil_disturbance <- as.numeric(seed.fil$microsite == "dig")

## different models####
glmm.pos <-glmmTMB(t4 ~ (1|species) + microsite+ treat + density,
                  family = poisson, data = seed.fil)
summary(glmm.pos)
check_model(glmm.pos)

glmm.nb <- glmmTMB(t4 ~ microsite+ treat + density + (1|species),
                 family = nbinom1, data = seed.fil)
summary(glmm.nb)
check_model(glmm.nb)

glmm.nb2 <- glmmTMB(t4 ~ microsite+ treat + density + (1|species),
                   family = nbinom2, data = seed.fil)
summary(glmm.nb2)
check_model(glmm.nb2)

glmm.nb.offset <-glmmTMB(t4 ~ microsite+ treat + density + (1|species), offset = t1,
                  family = nbinom1, data = seed.fil)
summary(glmm.nb.offset)
check_model(glmm.nb.offset)

glmm.nb.inter <-glmmTMB(t4 ~ (1|species) + microsite * treat * density,
                  family = nbinom1, data = seed.fil)
summary(glmm.nb.inter)
check_model(glmm.nb.inter)

glmm.zi.nb <-glmmTMB(t4 ~ (1|species) + microsite+ treat + density,
                 family = nbinom1, ziformula = ~microsite, data = seed.fil)
summary(glmm.zi.nb)
check_model(glmm.zi.nb)

glmm.zi.pos <-glmmTMB(t4 ~ (1|species) + microsite+ treat + density,
                       family = poisson, ziformula = ~microsite, data = seed.fil)
summary(glmm.zi.pos)

#check which model has the best AIC
AIC(glmm.pos, glmm.nb, glmm.nb2, glmm.nb.offset, glmm.nb.inter, glmm.zi.nb, glmm.zi.pos)

#check fit ####
#diagnosing 
fixef(glmm.nb)
diagnose(glmm.nb)

#model checks using dharma
#test for zero inflation and dispersion
testZeroInflation(glmm.nb) #looks good
testDispersion(glmm.nb) #looks good
simulationOutput <- simulateResiduals(fittedModel = glmm.nb)
plot(simulationOutput)
testResiduals(simulationOutput)

#extract predictions and coefs ####
# get predictions from model for specific treatments
nb_int <- ggpredict(glmm.nb , terms = c("microsite", "treat","density"), bias_correction = TRUE)
nb_int

predicted.counts <- as.data.frame(nb_int)

#extract coefficients from the model
coef <- as.data.frame(broom.mixed::tidy(glmm.nb, conf.int = TRUE))
coef <- filter(coef, effect == "fixed")
coef$term <- as.factor(coef$term)

# plots #####
### plot predicted counts 
partial.dep <- ggplot(predicted.counts, aes(y=predicted, x = facet, fill = group, shape = x)) + 
  geom_pointrange(aes(ymin=conf.low, ymax=conf.high), position=position_dodge(0.5), size = .7)+
  scale_fill_manual(name = " ", labels=c("Ungrazed", "Grazed"), values=c("white", "black"))+
  scale_shape_manual(name = "Microsite", labels=c("Crust", "Dig"), values=c(21,23))+
  #geom_vline(xintercept = 0, size = .3, linetype = 2)+
  theme_cowplot(font_size = 12)+
 ylab("Seedling count")+
  xlab("Sowing Density")+
  guides(fill = guide_legend(override.aes = list(shape = 21)))+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.5))

#coefficient plot
coef.plt<-ggplot(coef, aes(x=estimate, y=term)) + 
  geom_pointrange(aes(xmin=conf.low, xmax=conf.high))+
  scale_y_discrete(labels = c("treatG" = "Grazing", 
                              "micrositedig" = "Microsite: Digging",
                              "density.L" = "Sowing density (60)"), limits=rev)+
  geom_vline(xintercept = 0, size = .3, linetype = 2)+
  theme_cowplot(font_size = 12)+
  ylab("Parameter")+
  xlab("Estimate")+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=.5))

#put the two plots together
ggarrange(coef.plt, partial.dep, labels = "auto", widths = c(1.5, 1.2))




