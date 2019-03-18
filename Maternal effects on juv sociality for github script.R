######
# Code related to the AJP DP paper on maternal effects on social behavior in blue monkeys
######


library(gamlss)
library(tidyverse)
library(RCurl)
library(lubridate)

# (A) load data sets ------

# Data for first year models and average lifetime models (parts 1-3)
mjdf_ajp <- read.csv(text=getURL("https://raw.githubusercontent.com/Gavago/Maternal-effects-on-social-behavior-in-a-Cercopithecoid-monkey/master/maternal-effects-blue-monkeys-social.csv"),stringsAsFactors = F, header=T) %>%
  mutate(mum = as.factor(mum))  
dim(mjdf_ajp) #41 20


# Data for mum juvenile association over time (part 4)
url <- getURL("https://raw.githubusercontent.com/Gavago/Maternal-effects-on-social-behavior-in-a-Cercopithecoid-monkey/master/maternal-effects-blue-monkeys-social-groom-and-near-data-by-year.csv")
mjdf_year_ajp <- read.csv(text = url, stringsAsFactors = F, header=T) %>%
  mutate(mum = as.factor(mum)) %>%
  filter(!end %in% c(ymd("2015-08-01"), ymd("2010-10-30"), ymd("2012-04-09")))
dim(mjdf_year_ajp) #183 8

z. <- function(x) scale(x)


# (1) First-year models (Mother's behavior in subjects' first year of life) -------

# (1a) Grooming behavior:
# Predict juvenile grooming with any partner (not mum)
  #by mothers grooming with any groupmates (not juv)
fy.gg.all<-gamlss(j.p.gg ~ z.(fy.m.p.gg_not_juv) + sex + z.(age)  + z.(m.rank.sp) +  z.(fy.m.p.gg_not_juv)*sex+ random(mum) , family="BE",data=mjdf_ajp)
  #by mothers grooming with juvenile as infant
fy.mj.gg.all<-gamlss(j.p.gg ~ z.(fy.mj.gg) + sex + z.(age)  + z.(m.rank.sp) +  z.(fy.mj.gg)*sex + random(mum) , family="BE",data=mjdf_ajp)

# Predict juvenile grooming with peers (other juveniles)
  # by mothers grooming with her peers (adult females)
fy.gg.prs<-gamlss(p.gg.j ~ z.(fy.p.gg.af) + sex + z.(age)  + z.(m.rank.sp) +  z.(fy.p.gg.af)*sex+ random(mum) , family="BE",data=mjdf_ajp)
  # by mother's groomin with juvenile as infant
fy.mj.gg.prs<-gamlss(p.gg.j ~ z.(fy.mj.gg) + sex + z.(age)  + z.(m.rank.sp) +  z.(fy.mj.gg)*sex + random(mum) , family="BE",data=mjdf_ajp)


# (1b) Time near:
# Predict juvenile time near any partner (not mum)
  #by mothers grooming with any groupmates (not juv)
fy.rc.all<-gamlss(j.p.rc ~ z.(fy.m.p.rc_not_juv) + sex + z.(age)  + z.(m.rank.sp) +  z.(fy.m.p.rc_not_juv)*sex+ random(mum) , family="BE",data=mjdf_ajp)
  #by mothers grooming with juvenile as infant
fy.mj.rc.all<-gamlss(j.p.rc ~ z.(fy.mj.rc) + sex + z.(age)  + z.(m.rank.sp) +  z.(fy.mj.rc)*sex + random(mum) , family="BE",data=mjdf_ajp)

# Predict juvenile grooming with peers (other juveniles)
  # by mothers grooming with her peers (adult females)
fy.rc.prs<-gamlss(p.rc.j ~ z.(fy.p.rc.af) + sex + z.(age)  + z.(m.rank.sp) +  z.(fy.p.rc.af)*sex+ random(mum) , family="BE",data=mjdf_ajp)
  # by mother's groomin with juvenile as infant
fy.mj.rc.prs<-gamlss(p.rc.j ~ z.(fy.mj.rc) + sex + z.(age)  + z.(m.rank.sp) +  z.(fy.mj.rc)*sex + random(mum) , family="BE",data=mjdf_ajp)


# (1c) Playing:
# Predict juvenile time playing with any partner (not mum)
  #by mothers grooming with any groupmates (not juv)
fy.pl.all<-gamlss(j.p.pl ~ z.(fy.m.p.gg_not_juv) + sex + z.(age)  + z.(m.rank.sp) +  z.(fy.m.p.gg_not_juv)*sex+ random(mum) , family="BE",data=mjdf_ajp)
  #by mothers grooming with juvenile as infant
fy.mj.pl.all<-gamlss(j.p.pl ~ z.(fy.mj.gg) + sex + z.(age)  + z.(m.rank.sp) +  z.(fy.mj.gg)*sex + random(mum) , family="BE",data=mjdf_ajp)

# Predict juvenile playing with peers (other juveniles)
  # by mothers grooming with her peers (adult females)
fy.pl.prs<-gamlss(p.pl.j ~ z.(fy.p.gg.af) + sex + z.(age)  + z.(m.rank.sp) +  z.(fy.p.gg.af)*sex+ random(mum) , family="BEZI",data=mjdf_ajp, method=RS(50))
  # by mother's groomin with juvenile as infant
fy.mj.pl.prs<-gamlss(p.pl.j ~ z.(fy.mj.gg) + sex + z.(age)  + z.(m.rank.sp) +  z.(fy.mj.gg)*sex + random(mum) , family="BE",data=mjdf_ajp)




# (2) Significant interactions and Confidence intervals adjusted for false discovery rate ----

# identify all model names
gams <- Filter( function(x) 'gamlss' %in% class(get(x) ), ls())
length(gams) 

# gather model objects
mods<-list()
for(i in seq(gams)){
  mods[[i]]<-get(gams[i])  
}
names(mods)<-gams

# store model info
mod_info<-vector("list",length=length(mods))
names(mod_info)<-names(mods)
for(j in seq(mods)){
  info1<-setNames(data.frame(coefs=coef(mods[[j]]), confint(mods[[j]])), c("coefs","loCI","upCI"))  %>%
    mutate(pred = rownames(.),se = (upCI-coefs)/1.96, model = names(mods[j])) %>%
    filter(!is.na(coefs))
  
  mod_info[[j]]<-info1
}
mod_info

# create function to identify where the influence of mum behavior on juvenile behavior is significant and moderated by sex
sig_int<-function(m){
  m1<-m[6,]
  sig<- m1$loCI > 0 & m1$upCI > 0 | m1$loCI < 0 & m1$upCI < 0
  return(sig)
}

# select models
select <- unlist(lapply(mod_info, sig_int))
sig_mods<-mods[select]
names(sig_mods)

# adjust confidence intervals for false discovery rate
adj_mod_info <- vector("list",length=length(sig_mods))
adj_mod_names <- names(sig_mods)
names(adj_mod_info) <- adj_mod_names

for(j in seq(sig_mods)){
  info <- setNames(data.frame(coef(sig_mods[[j]]), 
                              confint(sig_mods[[j]]),
                              confint(sig_mods[[j]], level = 1 - 0.05/6*4)), c("coefs","loCI","upCI","consloCI","consupCI")) %>%
    mutate(pred = rownames(.), se = (upCI - coefs)/1.96, model = names(sig_mods)[j]) %>%
    filter(!is.na(coefs))
  
  adj_mod_info[[j]]<-info
}

adj_mod_info


# (3) Avg lifetime models (Mother's behavior in subjects' development) -----

#Juvenile grooming with any groupmate, predicted by grooming with mothers over lifetime
lavg.mj.gg.all<-gamlss(j.p.gg ~ z.(mj.gg_avg_lifetime) + sex + z.(age)  + z.(m.rank.sp) +  z.(mj.gg_avg_lifetime)*sex + random(mum) , family="BE",data=mjdf_ajp)

#Juvenile grooming with any peer, predicted by grooming with mothers over lifetime
lavg.mj.gg.prs <- gamlss(p.gg.j ~ z.(mj.gg_avg_lifetime) + sex + z.(age)  + z.(m.rank.sp) +  z.(mj.gg_avg_lifetime)*sex+ random(mum) , family="BE",data=mjdf_ajp, method=RS(50))

#Juvenile playing with any peer, predicted by grooming with mothers over lifetime
lavg.pl.prsg<-gamlss(p.pl.j ~ z.(p.gg.af_avg_lifetime) + sex + z.(age)  + z.(m.rank.sp) + z.(p.gg.af_avg_lifetime)*sex+ random(mum) , family="BEZI",data=mjdf_ajp, method=RS(50))

#lavg.rc.all<-gamlss(j.p.rc ~ z.(m.p.rc_not_juv_avg_lifetime) + sex + z.(age)  + z.(m.rank.sp) +  z.(m.p.rc_not_juv_avg_lifetime)*sex + random(mum) , family="BE",data=mjdf_ajp)

summary(lavg.mj.gg.all)
summary(lavg.mj.gg.prs)
summary(lavg.pl.prsg)


# (4) Mum juvenile association over time ----

# Mother juvenile time grooming
time.gg<-gamlss(mj.gg ~ cs(mum_juv_year) + sex +  cs(mum_juv_year)*sex + random(mum), data= mjdf_year_ajp, family="BEZI")
summary(time.gg)

# Mother juvenile time near
time.rc <- gamlss(mj.rc ~ cs(mum_juv_year) + sex + cs(mum_juv_year)*sex + random(mum), data= mjdf_year_ajp, family="BEZI", method = RS(50))
summary(time.rc)


# (5) Visualize mum behavior in first year predict juvenile behavior ----

#colorblind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot1 <- mjdf_ajp %>%
  ggplot(.) + geom_point(aes(x = 100*fy.mj.gg, y = 100*j.p.gg, color = sex, shape = sex), show.legend = F) +
  geom_smooth(aes(x = 100*fy.mj.gg, y = 100*j.p.gg, group = sex, color = sex, linetype = sex), method = "lm", show.legend = F) +
  scale_colour_manual(values=cbPalette) +
  labs(y = "", x = "% Time mom groom w juvenile during infancy") +
  ggtitle("a) % Time juvenile groom with groupmates *") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

plot2 <- mjdf_ajp %>%
  ggplot(.) + geom_point(aes(x = 100*fy.mj.gg, y = 100*p.gg.j, color = sex, shape = sex)) +
  geom_smooth(aes(x = 100*fy.mj.gg, y = 100*p.gg.j, group = sex, color = sex, linetype = sex), method = "lm") +
  scale_colour_manual(values=cbPalette) +
  labs(y = "", x = "% Time mom groom w juvenile during infancy") +
  ggtitle("b) % Time juvenile groom with peers *") +
  theme(plot.title = element_text(hjust = 0.5, size = 12), legend.position = c(0.95,0.9))

grid.arrange(plot1, plot2, nrow = 1, ncol = 2)

plot3 <- mjdf_ajp %>%
  ggplot(.) + geom_point(aes(x = 100*fy.p.gg.af, y = 100*p.pl.j, color = sex, shape = sex), show.legend = F) +
  geom_smooth(aes(x = 100*fy.p.gg.af, y = 100*p.pl.j, group = sex, color = sex, linetype = sex), method = "lm", show.legend = F) +
  scale_colour_manual(values=cbPalette) +
  labs(y = "", x = "% Time mom groom w peers during juvenile infancy") +
  ggtitle("c) % Time juvenile play with peers *") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

plot4 <- mjdf_ajp %>%
  ggplot(.) + geom_point(aes(x = 100*fy.m.p.rc_not_juv, y = 100*j.p.rc, color = sex, shape = sex), show.legend = F) +
  geom_smooth(aes(x = 100*fy.m.p.rc_not_juv, y = 100*j.p.rc, group = sex, color = sex, linetype = sex), method = "lm", show.legend = F) +
  scale_colour_manual(values=cbPalette) +
  labs(y = "", x = "% Time mom near groupmates during \n juvenile infancy") +
  ggtitle("d) % Time juvenile near groupmates") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))
grid.arrange(plot3, plot4, nrow = 1, ncol = 2)

# (6) Visualize mum juvenile association over lifetime ----

#colorblind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

mjdf_year_ajp %>%
  mutate(`a) grooming offspring` = mj.gm, `b) groomed by offspring` = mj.gg - mj.gm,
         `c) near offspring` =  mj.rc) %>%
  gather(key = aff, value = percent_time, `a) grooming offspring`, `b) groomed by offspring`, `c) near offspring`) %>%
  ggplot(.) +
  geom_point(aes(y = 100*percent_time, x = mum_juv_year, color = sex, shape = sex), position = "jitter") +
  geom_smooth(aes(y = 100*percent_time, x = mum_juv_year, color = sex, linetype = sex),method = "loess") +
  scale_colour_manual(values=cbPalette) +
  labs(y = "", x = "Year of offspring life", color = "offspring sex", linetype = "offspring sex", shape = "offspring sex") +
  scale_x_continuous(breaks = seq(min(dat$mum_juv_year), max(dat$mum_juv_year), by = 1)) +
  ggtitle("% Time mother spent") +
  theme(strip.background =element_rect(fill="white"), plot.title = element_text(hjust = 0.5)) +
  facet_wrap(facets = vars(aff), scales = "free_y")
