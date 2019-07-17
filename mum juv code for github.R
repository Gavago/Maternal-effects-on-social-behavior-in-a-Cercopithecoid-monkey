#install.packages("RCurl")
#install.packages("tidyverse")
#install.packages("gamlss")

library(RCurl)
library(tidyverse)
library(gamlss)
library(lubridate)
z.<-function(x) scale(x)
select <- dplyr::select



#######
# Load and prep data
#######

# See data dictionaries under README.md file - https://github.com/Gavago/Maternal-effects-on-social-behavior-in-a-Cercopithecoid-monkey/blob/master/README.md

# Data for first year and overall development models
x <- getURL("https://raw.githubusercontent.com/Gavago/Maternal-effects-on-social-behavior-in-a-Cercopithecoid-monkey/master/maternal%20effects%20blue%20monkeys%20social.csv")
mjdf <- read_csv(file = x) %>%   # "mum juv data frame"
  mutate(juv = as.factor(juv), mum = as.factor(mum))
dim(mjdf) #41 21

# Data for mum juv association over each year of juv's life
y <- getURL("https://raw.githubusercontent.com/Gavago/Maternal-effects-on-social-behavior-in-a-Cercopithecoid-monkey/master/maternal-effects-blue-monkeys-social-groom-and-near-data-by-year.csv")
mum_juv_over_time <- read_csv(file = y)
dat <- mum_juv_over_time %>%
  mutate(juv = as.factor(juv), mum = as.factor(mum)) %>%
  filter(!end %in% c(ymd("2015-08-01"), ymd("2010-10-30"), ymd("2012-04-09")))
dim(dat) # 142 9
# filter out incomplete final years - 
# ending 08-01 - start of juv observation period
# ending 10-30 - last year fletcher (mum) alive
# ending 04-09 - last year artemann (mum) alive

########
# First-year models: relationship bt juv sociliaty and mother's sociality in juv's first year of life
########

# response = Juvenile giving and receiving grooming ========

#mother grooming any partner in juv first year influence on juv groom any partner
fy.gg.all<-gamlss(j.p.gg ~ z.(fy.m.p.gg_not_juv) + sex + z.(age)  + z.(m.rank.sp) + z.(fy.m.p.gg_not_juv)*sex+ random(mum) , family="BE",data=mjdf)
summary(fy.gg.all)

#mother grooming w juv in juv first year influence on juv gm w any partner
fy.mj.gg.all<-gamlss(j.p.gg ~ z.(fy.mj.gg) + sex + z.(age)  + z.(m.rank.sp) + z.(fy.mj.gg)*sex + random(mum) , family="BE",data=mjdf)
summary(fy.mj.gg.all)

#mother grooming w peers in juv first year influence on juv groom w peers
fy.gg.prs<-gamlss(p.gg.j ~ z.(fy.p.gg.af) + sex + z.(age)  + z.(m.rank.sp) +  z.(fy.p.gg.af)*sex+ random(mum) , family="BE",data=mjdf)
summary(fy.gg.prs)

#mother grooming with juv in juv first year influence on juv groom w peers
fy.mj.gg.prs <- gamlss(p.gg.j ~ z.(fy.mj.gg) + sex + z.(age)  + z.(m.rank.sp) +  z.(fy.mj.gg)*sex+ random(mum) , family="BE",data=mjdf, method=RS(50))
summary(fy.mj.gg.prs)

# response = Juvenile resting in 1 m and sitting in contact ===========

#Mother time resting or sitting in contact (rc) with any partner in juv first year influence on juvenile rc with any partner
fy.rc.all<-gamlss(j.p.rc ~ z.(fy.m.p.rc_not_juv) + sex + z.(age)  + z.(m.rank.sp) +  z.(fy.m.p.rc_not_juv)*sex + random(mum) , family="BE",data=mjdf)
summary(fy.rc.all)

#Mother time resting or sitting in contact (rc) with juv in juv first year influence on juvenile rc with any partner
fy.mj.rc.all<-gamlss(j.p.rc ~ z.(fy.mj.rc) + sex + z.(age)  + z.(m.rank.sp) +  z.(fy.mj.rc)*sex+ random(mum) , family="BE",data=mjdf)
summary(fy.mj.rc.all)

#Mother time resting or sitting in contact (rc) with adult female peers in juv first year influence on juvenile rc with juv peers
fy.rc.prs<-gamlss(p.rc.j ~ z.(fy.p.rc.af) + sex + z.(age) + z.(m.rank.sp) + z.(fy.p.rc.af)*sex+ random(mum) , family="BE",data=mjdf, method=RS(50))
summary(fy.rc.prs)

#Mother time resting or sitting in contact (rc) with juv in juv first year influence on juvenile rc with juv peers
fy.mj.rc.prs<-gamlss(p.rc.j ~ z.(fy.mj.rc) + sex + z.(age)  + z.(m.rank.sp) +  z.(fy.mj.rc)*sex+ random(mum) , family="BE",data=mjdf, method=RS(50))
summary(fy.mj.rc.prs)

# response = Juvenile playing =========

#mother grooming w any partner in juv first year influence on juv play w any partner
fy.pl.allg<-gamlss(j.p.pl ~ z.(fy.m.p.gg_not_juv) + sex + z.(age)  + z.(m.rank.sp) +  z.(fy.m.p.gg_not_juv)*sex+ random(mum) , family="BE",data=mjdf, method = RS(50))
summary(fy.pl.allg)

#mother grooming w juv in juv first year influence on juv play w any partner
fy.mj.pl.allg<-gamlss(j.p.pl ~ z.(fy.mj.gg) + sex + z.(age)  + z.(m.rank.sp) +  z.(fy.mj.gg)*sex+ random(mum) , family="BE",data=mjdf, method=RS(50))
summary(fy.mj.pl.allg)

#mother grooming w adult female peers in juv first year influence on juv play w juv peers
fy.pl.prsg<-gamlss(p.pl.j ~ z.(fy.p.gg.af) + sex + z.(age)  + z.(m.rank.sp) + z.(fy.p.gg.af)*sex+ random(mum) , family="BEZI",data=mjdf, method=RS(50))
summary(fy.pl.prsg)

#mother grooming w juv in juv first year influence on juv play w juv peers
fy.mj.pl.prsg<-gamlss(p.pl.j ~ z.(fy.mj.gg) + sex + z.(age)  + z.(m.rank.sp) + z.(fy.mj.gg)*sex+ random(mum) , family="BEZI",data=mjdf, method=RS(50))
summary(fy.mj.pl.prsg)

# gather first-year GAM outputs =======

#get gam names
gams <- Filter( function(x) 'gamlss' %in% class(get(x)), ls())
length(gams) # should be 12
gams #names

#save outputs in list
mods<-list()
for(i in seq(gams)){
  mods[[i]]<-get(gams[i])  
}
names(mods)<-gams

#extract info - non sig model info presented in supplemental table 1

mod.info<-vector("list",length=length(mods))
names(mod.info)<-names(mods)
for(j in seq(mods)){
  m<-mods[[j]]
  info1<-setNames(data.frame(coefs=coef(m), confint(m)), c("coefs","loCI","upCI"))  %>%
    mutate(pred = rownames(.),se = (upCI-coefs)/1.96, model = names(mods[j])) %>%
    filter(!is.na(coefs))
  
  mod.info[[j]]<-info1
}

mod.info

# find which models have significant interactions ====

sig.int<-function(m){ #in which models is influence of mum beh on juv beh significantly by sex
  m1<-m[6,] #row 6 is interaction
  sig<- (m1$loCI>0 & m1$upCI>0) | (m1$loCI<0 & m1$upCI<0)
  return(sig)
}

select<-unlist(lapply(mod.info, sig.int))
sig.mods<-mods[select] #sig interactions
names(sig.mods)


# correct for False Coverage Rate in confidence intervals - sig mod info presented in Table 3 ======

sig.mod.info<-vector("list",length=length(sig.mods))
sig.mod.names <- names(sig.mods)
names(sig.mod.info)<- sig.mod.names

for(j in seq(sig.mods)){
  m<-sig.mods[[j]]
  info <- setNames(data.frame(coef(m), 
                              confint(m),
                              confint(m, level = 1 - 0.05/6*4)), c("coefs","loCI","upCI","consloCI","consupCI")) %>%
    mutate(pred = rownames(.), se = (upCI - coefs)/1.96, model = names(sig.mods)[j]) %>%
    filter(!is.na(coefs))
  
  sig.mod.info[[j]]<-info
}
sig.mod.info


# trade-off models - grooming w mothers vs. others *as a juv* ======

#juv time grooming w mom trade off w grooming non mothers
gg_trade_off_mom_others <- gamlss(j.p.gg ~ z.(jm.gg) + sex + z.(age) + z.(m.rank.sp) + random(mum) , family="BE",data=mjdf, method=RS(50))
summary(gg_trade_off_mom_others)
confint(gg_trade_off_mom_others)

gg_trade_off_mom_others_int <- gamlss(j.p.gg ~ z.(jm.gg) + sex + z.(age) + z.(m.rank.sp) + z.(jm.gg)*sex + random(mum) , family="BE",data=mjdf, method=RS(50))
summary(gg_trade_off_mom_others_int)
confint(gg_trade_off_mom_others_int)

# juv time grooming w mom trade off w peers
gg_trade_off_mom_peers <- gamlss(p.gg.j ~ z.(jm.gg) + sex + z.(age) + z.(m.rank.sp) + random(mum) , family="BE",data=mjdf, method=RS(50))
summary(gg_trade_off_mom_peers)

gg_trade_off_mom_peers_int <- gamlss(p.gg.j ~ z.(jm.gg) + sex + z.(age) + z.(m.rank.sp) + z.(jm.gg)*sex + random(mum) , family="BE",data=mjdf, method=RS(50))
summary(gg_trade_off_mom_peers_int)


#######
# Overall development (or lifetime) models: relationsip bt juv sociality and avg mother's sociality throughout juv's development
#######

#avg mother grooming w juv over each year of juv life influence on juv gm w any partner
lavg.mj.gg.all<-gamlss(j.p.gg ~ z.(mj.gg_avg_lifetime) + sex + z.(age)  + z.(m.rank.sp) +  z.(mj.gg_avg_lifetime)*sex + random(mum) , family="BE",data=mjdf)
#avg mother grooming w juv over each year of juv life influence on juv gm w peers
lavg.mj.gg.prs <- gamlss(p.gg.j ~ z.(mj.gg_avg_lifetime) + sex + z.(age)  + z.(m.rank.sp) +  z.(mj.gg_avg_lifetime)*sex+ random(mum) , family="BE",data=mjdf, method=RS(50))
#avg mother grooming w adult female peers over each year of juv life influence on juv play w peers
lavg.pl.prsg<-gamlss(p.pl.j ~ z.(p.gg.af_avg_lifetime) + sex + z.(age)  + z.(m.rank.sp) + z.(p.gg.af_avg_lifetime)*sex+ random(mum) , family="BEZI",data=mjdf, method=RS(50))
#avg mother resting in 1 m or in sitting in contact with any partner influence on juv rc w any partner
lavg.rc.all<-gamlss(j.p.rc ~ z.(m.p.rc_not_juv_avg_lifetime) + sex + z.(age)  + z.(m.rank.sp) +  z.(m.p.rc_not_juv_avg_lifetime)*sex + random(mum) , family="BE",data=mjdf)

summary(lavg.mj.gg.all)
summary(lavg.mj.gg.prs)
summary(lavg.pl.prsg)
summary(lavg.rc.all)

lifetime_mods<- list(lavg.mj.gg.all, lavg.mj.gg.prs, lavg.pl.prsg, lavg.rc.all)
lifetime.mod.names <- c("lavg.mj.gg.all", "lavg.mj.gg.prs", "lavg.pl.prsg", "lavg.rc.all")
names(lifetime_mods)<- lifetime.mod.names


lifetime.mod.info<-vector("list",length=length(lifetime_mods))
names(lifetime.mod.info)<- lifetime.mod.names


for(j in seq(lifetime_mods)){
  m<-lifetime_mods[[j]]
  info <- setNames(data.frame(coef(m), 
                              confint(m)),
                              c("coefs","loCI","upCI")) %>%
    mutate(pred = rownames(.), se = (upCI - coefs)/1.96, model = lifetime.mod.names[j]) %>%
    filter(!is.na(coefs))
  
  lifetime.mod.info[[j]]<-info
}

lifetime.mod.info # info presented in table 4


#######
# Mother juv association over time
#######

time.rc<-gamlss(mj.rc ~ cs(mum_juv_year) + sex + cs(mum_juv_year)*sex + random(mum), data= dat, family="BEZI", method = RS(50))
summary(time.rc)
confint(time.rc)

time.gg<-gamlss(mj.gg ~ cs(mum_juv_year) + sex +  cs(mum_juv_year)*sex + random(mum), data= dat, family="BEZI")
summary(time.gg)
confint(time.gg)


######
# Visualizations
######

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #colorblind palette

# Figure 2 ------
plot1 <- mjdf %>%
  ggplot(.) + geom_point(aes(x = 100*fy.mj.gg, y = 100*j.p.gg, color = sex, shape = sex), show.legend = F) +
  geom_smooth(aes(x = 100*fy.mj.gg, y = 100*j.p.gg, group = sex, color = sex, linetype = sex), method = "lm", show.legend = F) +
  scale_colour_manual(values=cbPalette) +
  labs(y = "", x = "% Time mom groom w subject during infancy") +
  ggtitle("a) % Time juvenile subject groom w any partner *") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

plot2 <- mjdf %>%
  ggplot(.) + geom_point(aes(x = 100*fy.mj.gg, y = 100*p.gg.j, color = sex, shape = sex)) +
  geom_smooth(aes(x = 100*fy.mj.gg, y = 100*p.gg.j, group = sex, color = sex, linetype = sex), method = "lm") +
  scale_colour_manual(values=cbPalette) +
  labs(y = "", x = "% Time mom groom w subject during infancy") +
  ggtitle("b) % Time juvenile subject groom with peers *") +
  theme(plot.title = element_text(hjust = 0.5, size = 12), legend.position = c(0.95,0.9))

plot3 <- mjdf %>%
  ggplot(.) + geom_point(aes(x = 100*fy.p.gg.af, y = 100*p.pl.j, color = sex, shape = sex), show.legend = F) +
  geom_smooth(aes(x = 100*fy.p.gg.af, y = 100*p.pl.j, group = sex, color = sex, linetype = sex), method = "lm", show.legend = F) +
  scale_colour_manual(values=cbPalette) +
  labs(y = "", x = "% Time mom groom w peers during subject infancy") +
  ggtitle("c) % Time juvenile subject play with peers *") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

plot1
plot2
plot3

# Figure 3 ------
plot4 <- dat %>%
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
plot4
