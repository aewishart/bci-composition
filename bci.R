---
title: "BCI-3species - Wishart et al"
output: pdf_document
---
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

#This is the code for the manuscript studying the relationship between body condition indices and body composition in three squirrel species. bioRxiv preprint DOI: 10.1101/2023.01.31.524791 Feb 2023

#Load libraries
```{r}
library(tidyverse)
library(lubridate)
library(lme4)
library(lmerTest)
library(MCMCglmm)
library(ggplot2)
library(visreg)
library(plotrix)
library(viridis)
library(multcomp)
library(AICcmodavg)
library(smatr)

select = dplyr::select #necessary as MASS also has a select function
filter = dplyr::filter

ctrl <- lmerControl(optCtrl=list(xtol_rel=1e-6)) # Warning otherwise
```
#Define functions 
```{r}
std <- function(x) sd(x)/sqrt(length(x))

CV <- function(x){
        (sd(x)/mean(x))*100
}

```


#Get data 
```{r}
rsq <- read.csv("BCI-rsq.csv")

btpd <- read.csv("BCI-btpd.csv")

cgs <- read.csv("BCI-cgs.csv") #CGS data for both seasons

cgsspr <- cgs %>%
  filter(season == "spring")


cgsfall <- cgs %>%
  filter(season == "fall")

```
#Select species to anlyze here by uncommenting data <- [species of interest]

#NOTE - this is to ensure each species is run with identical code. Figure and table text referencing species and/or season MUST be adjusted accordingly. e.g., if running RSQ code, must manually change figure titles to Red squirrel - otherwise they may read the wrong species/season. 
```{r}

#Red squirrel
#data <- rsq  #uncomment if running code for this species 

#Black-tailed prairie dog
#data <- btpd  #comment out if running code for a different species 

#CGS
#data <- cgsfall

#uncomment for median daysdiff in autumn for CGS
#which(rank(data$days.diffprehib) == length(data$days.diffprehib) / 2 + 1)

data <- cgsspr

```


#Split data by sex, scale, and rejoin
```{r}
#Females
data.f <- data %>%
  filter(sex == "female")

data.f <- data.f %>%
  mutate(logwgt = log(wgt +1 ), 
    logrhf = log(rhf.avg +1 ), 
    logzyg = log(zyg.avg +1 ), 
    loglean = log(lean +1 ), 
    logfat = log(fat + 1), 
    logleanp = log(leanp +1 ), 
    logfatp = log(fatp + 1))

#Scale within sex
data.f <- data.f %>% 
   mutate(across(c("logwgt", "logrhf", "logzyg", "loglean", "logfat", "logleanp", "logfatp"), ~(c(scale(.)))))


data.f <- data.f %>%
 mutate(
    wgt.scl = logwgt, 
    rhf.scl = logrhf, 
    zyg.scl = logzyg, 
    lean.scl = loglean, 
    fat.scl = logfat, 
    leanp.scl = logleanp,
    fatp.scl = logfatp)

summary(data.f)

#Relationship between zyg and rhf in females?
cor.test(data.f$rhf.scl, data.f$zyg.scl) 



#Males 
data.m <- data%>%
  filter(sex == "male")

data.m <- data.m %>%
  mutate(logwgt = log(wgt +1 ), 
    logrhf = log(rhf.avg +1 ), 
    logzyg = log(zyg.avg +1 ), 
    loglean = log(lean +1 ), 
    logfat = log(fat + 1), 
    logleanp = log(leanp +1 ), 
    logfatp = log(fatp + 1))

#Scale within sex
data.m <- data.m %>% 
   mutate(across(c("logwgt", "logrhf", "logzyg", "loglean", "logfat", "logleanp", "logfatp"), ~(c(scale(.)))))


data.m <- data.m %>%
 mutate(
    wgt.scl = logwgt, 
    rhf.scl = logrhf, 
    zyg.scl = logzyg, 
    lean.scl = loglean, 
    fat.scl = logfat, 
    leanp.scl = logleanp,
    fatp.scl = logfatp)

summary(data.m)

#Relationship between zyg and rhf in males?
cor.test(data.m$rhf.scl, data.m$zyg.scl) 


#Rejoin male and afemale data 
data <- rbind(data.m, data.f)
```

#Relationship between RHF and mass
```{r, fig.width=5,fig.height=3}
#Linear model of RHF with mass and sex
rhf.m1 = lm(rhf.avg ~ wgt, 
            data = data)
summary(rhf.m1) 

{plot((wgt.scl) ~ rhf.avg, data= data, xlab= "Right hind foot length (scaled)", ylab = " Body Mass (scaled)", main = "Regression of scaled body mass on scaled RHF")}

figa <- ggplot(data = data, 
              aes(x=wgt, y=rhf.avg, fill=sex, shape = sex)) +
 geom_point(aes(shape=sex,color=sex),size=2,stroke = 1)+
  scale_shape_manual(values = c(16:17))+ 
  scale_fill_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue'))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
   ylab("Right hind foot length (mm)") + xlab("Body mass (g)") +
   ggtitle('C) Ground squirrels, pre-winter') +                                                            
   stat_smooth(method = "lm", se = TRUE, alpha = 0.5, aes(color = sex))+                                
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
      theme(#legend.position = "top",
   # legend.title = element_blank(),
  #  legend.key = element_blank(),
   # legend.text = element_text(size = 10),
    axis.ticks = element_blank(),
    axis.title=element_text(size=12), 
    axis.text = element_text(size = 12),
    strip.background = element_rect(fill="white", color = "black", size = 1),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill=NA, size = 1)) #
 # theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) #+
#  facet_wrap(~grid, ncol = 1)
figa
```
#Relationship between Zyg and mass
```{r, fig.width=5,fig.height=3}
#ZYG by mass and sex 
zyg.m1 = lm(zyg.avg ~ wgt + sex, 
            data = data)
summary(zyg.m1) 

{plot((wgt.scl) ~ zyg.avg, data= data, xlab= "Zygomatic width (scaled)", ylab = " Body Mass (scaled)", main = "Regression of scaled body mass on scaled zyg")}


figa <- ggplot(data = data, 
              aes(x=wgt, y=zyg.avg, fill=sex, shape = sex)) +
 geom_point(aes(shape=sex,color=sex),size=2,stroke = 1)+
  scale_shape_manual(values = c(16:17))+ 
  scale_fill_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue'))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
   ylab("Zygomatic width (mm)") + xlab("Body mass (g)") +
   ggtitle('B) Prairie dogs, pre-winter') +
   stat_smooth(method = "lm", se = TRUE, alpha = 0.5, aes(color = sex))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
      theme(#legend.position = "top",
   # legend.title = element_blank(),
  #  legend.key = element_blank(),
   # legend.text = element_text(size = 10),
    axis.ticks = element_blank(),
    axis.title=element_text(size=12), 
    axis.text = element_text(size = 12),
    strip.background = element_rect(fill="white", color = "black", size = 1),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill=NA, size = 1)) #+
 # theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) #+
#  facet_wrap(~grid, ncol = 1)
figa
```
#Relationship between zyg and RHF
```{r, fig.width=5,fig.height=3}

rzm1= lm(rhf.avg ~ zyg.avg + sex, 
            data = data)
summary(rzm1) 


figS1 <- ggplot(data = data, 
              aes(x=rhf.avg, y=zyg.avg, fill=sex, shape = sex)) +
 geom_point(aes(shape=sex,color=sex),size=2,stroke = 1)+
  scale_shape_manual(values = c(16:17))+ 
  scale_fill_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue'))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','#0072b2','turquoise4','green','blue')) +
   ylab("Zygomatic width (mm)") + xlab("Right hind foot length (mm)") +
   ggtitle('A) Red squirrels, pre-winter') +                                                        #<-------Change title to match data
   stat_smooth(method = "lm", se = TRUE, alpha = 0.5, aes(color = sex))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','#0072b2','turquoise4','green','blue')) +
      theme(#legend.position = "top",
   # legend.title = element_blank(),
  #  legend.key = element_blank(),
   # legend.text = element_text(size = 10),
    axis.ticks = element_blank(),
    axis.title=element_text(size=12), 
    axis.text = element_text(size = 12),
    strip.background = element_rect(fill="white", color = "black", size = 1),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill=NA, size = 1)) +
  annotate("text", x=Inf, y=-Inf, hjust=1.2,vjust=-0.3,label = bquote("R[adj.]^2 < 0.01"), parse = TRUE) #<-------Change R2 value to match data
 # theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) #+
#  facet_wrap(~grid, ncol = 1)
figS1




#Relationship between zyg and rhf ?
#Females
cor.test(data.f$rhf.avg, data.f$zyg.avg) 


#Relationship between zyg and rhf ?
#Males
cor.test(data.m$rhf.avg, data.m$zyg.avg)

filter()
```


#Summary data for tables
```{r}

data %>% count(sex) 
#summarize by sex 
summbysex <- data %>%                                       
  group_by(sex) %>%                        
  summarize(wgt.mn = mean(wgt), 
            wgt.sem = std(wgt), 
            rhf.mn = mean(rhf.avg), 
            rhf.sem = std(rhf.avg), 
            zyg.mn = mean(zyg.avg), 
            zyg.sem = std(zyg.avg),
            fat.mn = mean(fat), 
            fat.sem = std(fat), 
            lean.mn = mean(lean), 
            lean.sem = std(lean), 
            fatp.mn = mean(fatp), 
            leanp.mn = mean(leanp))  


#Coefficient of variation 
CoV <- data %>%                                       
  group_by(sex) %>%                        
  summarize(cvwgt = CV(wgt), 
            cvrhf = CV(rhf.avg), 
            cvzyg = CV(zyg.avg), 
            cvfat = CV(fat), 
            cvlean = CV(lean), 
            cvfatp = CV(fatp), 
            cvleanp = CV(leanp)) 

```


#Sex differences  + boxplots to visualize
```{r}
#Mass
massbox1 <- ggplot(data, aes(sex, wgt))
massbox1 + geom_boxplot() +
               xlab("Sex") + ylab("Body Mass (g)")+
              theme(axis.text.x = element_text(size=20), axis.title = element_text(size = 20), axis.text.y = element_text(size=20))

wgt1 <- t.test(wgt ~ sex, data = data)
wgt1



#RHF
rhfbox1 <- ggplot(data, aes(sex, rhf.avg))
rhfbox1 + geom_boxplot() +
               xlab("Sex") + ylab("RHF (mm)") + ggtitle("Fall") + 
              theme(axis.text.x = element_text(size=15), axis.title = element_text(size = 18), axis.text.y =
                      element_text(size=15))

rhf.avg1 <- t.test(rhf.avg ~ sex, data = data)
rhf.avg1



#Zyg

zygbox3<- ggplot(data, aes(sex, zyg.avg))
zygbox3 + geom_boxplot() +
               xlab("Sex") + ylab("Zyg (mm)")+ ggtitle("Fall") + 
              theme(axis.text.x = element_text(size=15), axis.title = element_text(size = 18), axis.text.y = element_text(size=15))

zyg.avg1 <- t.test(zyg.avg ~ sex, data = data)
zyg.avg1



#Fat
fatbox4<- ggplot(data, aes(sex, fat))
fatbox4 + geom_boxplot() +
               xlab("Sex") + ylab("Fat %")+  ggtitle("Spring") + 
              theme(axis.text.x = element_text(size=15), axis.title = element_text(size = 18), axis.text.y = element_text(size=15), aspect.ratio = 1)

fat1 <- t.test(fat ~ sex, data = data)
fat1
 


#Lean
leanbox5<- ggplot(data, aes(sex, lean))
leanbox5 + geom_boxplot() +
               xlab("Sex") + ylab("Lean %")+ ggtitle("Fall") + 
              theme(axis.text.x = element_text(size=15), axis.title = element_text(size = 18), axis.text.y = element_text(size=15), aspect.ratio = 1)

lean1 <- t.test(lean ~ sex, data = data)
lean1

```

#Calculate BCI from singular skeletal measurements - RHF or ZYG on BM resids, Females only (Male code follows)
```{r}
####
#RHF
####
rhfm1.f <- lm(wgt.scl ~ rhf.scl, data=data.f)
attributes(rhfm1.f)
summary(rhfm1.f)

plot(rhfm1.f, which=1)



#Check normality 
plot(rhfm1.f)
shapiro.test(rhfm1.f$residuals)

#Residuals from this regression of log body mass ~ log RHF will serve as the body condition index RHF (BCI.RHF).

data.f$BCI.RHF <- rhfm1.f$residuals 





####
#ZYG
####
zygm1.f <- lm(wgt.scl ~ zyg.scl, data=data.f)
attributes(zygm1.f)
summary(zygm1.f)

plot(zygm1.f, which=1)


#Check normality 
plot(zygm1.f)
shapiro.test(zygm1.f$residuals)

#Residuals from this regression of log body mass ~ log RHF will serve as the body condition index RHF (BCI.RHF).

data.f$BCI.ZYG <- zygm1.f$residuals 

```

#Calculate BCI from singular skeletal measurements - RHF or ZYG on BM resids, Males only
```{r}

####
#RHF
####
rhfm1.m <- lm(wgt.scl ~ rhf.scl, data=data.m)
attributes(rhfm1.m)
summary(rhfm1.m)

plot(rhfm1.m, which=1)

#Check normality 
plot(rhfm1.m)
shapiro.test(rhfm1.m$residuals)

#Residuals from this regression of log body mass ~ log RHF will serve as the body condition index RHF (BCI.RHF).

data.m$BCI.RHF <- rhfm1.m$residuals 


####
#ZYG
####
zygm1.m <- lm(wgt.scl ~ zyg.scl, data=data.m)
attributes(zygm1.m)
summary(zygm1.m)

plot(zygm1.m, which=1)


#Check normality 
plot(zygm1.m)
shapiro.test(zygm1.m$residuals)

#Residuals from this regression of log body mass ~ log RHF will serve as the body condition index RHF (BCI.RHF).

data.m$BCI.ZYG <- zygm1.m$residuals 

```





#NEW - Scaled Mass Index (SMI)

#Females 

```{r}

#RHF SMI
head(data.f)

plot(data.f$rhf.avg, data.f$wgt)
plot(log(data.f$rhf.avg), log(data.f$wgt))

#Look at data 
loglm <- lm(log(wgt) ~ log(rhf.avg), data = data.f)
plot(loglm)

#Equation 1 to get slope of the standardized major axis (SMA) regression on ln-transformed data
bSMA <- coef(sma(log(wgt) ~ log(rhf.avg), data = data.f))[2]

bSMA

#L0 - arithmetic mean of x variable aka linear measure 
L0 = mean(data.f$rhf.avg)
L0

#Equation to calculate SMI and put into dataframe
data.f$SMI.rhf <- (data.f$wgt * (L0/data.f$rhf.avg)^bSMA)



#ZYG SMI

plot(data.f$zyg.avg, data.f$wgt)
plot(log(data.f$zyg.avg), log(data.f$wgt))

#Look at data 
loglm <- lm(log(wgt) ~ log(zyg.avg), data = data.f)
plot(loglm)

#Equation 1 to get slope of the standardized major axis (SMA) regression on ln-transformed data
bSMA <- coef(sma(log(wgt) ~ log(zyg.avg), data = data.f))[2]

bSMA

#L0 - arithmetic mean of x variable aka linear measure 
L0 = mean(data.f$zyg.avg)
L0

#Calculate SMI and put into dataframe
data.f$SMI.zyg <- (data.f$wgt * (L0/data.f$zyg.avg)^bSMA)

```

#Males SMI  

```{r}

#RHF SMI
head(data.m)

plot(data.m$rhf.avg, data.m$wgt)
plot(log(data.m$rhf.avg), log(data.m$wgt))

#Look at data 
loglm <- lm(log(wgt) ~ log(rhf.avg), data = data.m)
plot(loglm)

#Equation 1 to get slope of the standardized major axis (SMA) regression on ln-transformed data
bSMA <- coef(sma(log(wgt) ~ log(rhf.avg), data = data.m))[2]

bSMA

#L0 - arithmetic mean of x variable aka linear measure 
L0 = mean(data.m$rhf.avg)
L0

#Equation 2 to calculate SMI 
data.m$SMI.rhf <- (data.m$wgt * (L0/data.m$rhf.avg)^bSMA)



#ZYG SMI

plot(data.m$zyg.avg, data.m$wgt)
plot(log(data.m$zyg.avg), log(data.m$wgt))

#Look at data 
loglm <- lm(log(wgt) ~ log(zyg.avg), data = data.m)
plot(loglm)

#Equation 1 to get slope of the standardized major axis (SMA) regression on ln-transformed data
bSMA <- coef(sma(log(wgt) ~ log(zyg.avg), data = data.m))[2]

bSMA

#L0 - arithmetic mean of x variable aka linear measure 
L0 = mean(data.m$zyg.avg)
L0

#Equation 2 to calculate SMI 
data.m$SMI.zyg <- (data.m$wgt * (L0/data.m$zyg.avg)^bSMA)

```


#Join data
```{r} 
#Join the split-by-sex data together
data.bysex <- rbind(data.m, data.f)
head(data.bysex)

data.bysex <- data.bysex %>%
  dplyr::select(squirrel_id, 
                BCI.RHF, 
                BCI.ZYG, 
                SMI.rhf, 
                SMI.zyg)

data <- inner_join(data, data.bysex, by = c("squirrel_id"))

```

#Write data to CSV - this gets read in later for 3-species comparison once this code has been run for each species/season. 
```{r}

#Uncomment the species to which data is set, reiterate for each species.







#Red squirrel
# rsq <- data %>%
#   dplyr::select(sex,
# squirrel_id,
# wgt,
# rhf.avg,
# zyg.avg,
# fat,
# fat.sd,
# lean,
# lean.sd,
# leanp,
# fatp,
# logwgt,
# logrhf,
# logzyg,
# loglean,
# logfat,
# wgt.scl,
# rhf.scl,
# zyg.scl,
# lean.scl,
# fat.scl,
# BCI.ZYG,
# SMI.zyg)%>%
#   collect()%>%
#   mutate(BCI = BCI.ZYG,
#          SMI = SMI.zyg)%>%
#  collect()%>%
# dplyr::select(-BCI.ZYG, -SMI.zyg)
# 
# write.csv(rsq, "rsq-composition.csv")


# #Prairie dog
# btpd <- data %>%
#   mutate(season = "fall")%>%
#   dplyr::select(sex,
#                 squirrel_id,
#                 wgt,
#                 rhf.avg,
#                 zyg.avg,
#                 fat,
#                 fat.sd,
#                 lean,
#                 lean.sd,
#                 season,
#                 leanp,
#                 fatp,
#                 logwgt,
#                 logrhf,
#                 logzyg,
#                 loglean,
#                 logfat,
#                 wgt.scl,
#                 rhf.scl,
#                 zyg.scl,
#                 lean.scl,
#                 fat.scl,
# BCI.ZYG,
# SMI.zyg)%>%
#   collect()%>%
#   mutate(BCI = BCI.ZYG,
#          SMI = SMI.zyg)%>%
#  collect()%>%
# dplyr::select(-BCI.ZYG, -SMI.zyg)
# 
# write.csv(btpd, "btpd-composition.csv")

# #Ground squirrel fall
# cgsfall <- data %>%
#   mutate(fatsd = fat.sd,
#          leansd = lean.sd)%>%
#   collect()%>%
#   dplyr::select(sex,
# squirrel_id,
# wgt,
# rhf.avg,
# zyg.avg,
# fat,
# fat.sd,
# lean,
# lean.sd,
# season,
# leanp,
# fatp,
# logwgt,
# logrhf,
# logzyg,
# loglean,
# logfat,
# wgt.scl,
# rhf.scl,
# zyg.scl,
# lean.scl,
# fat.scl,
# BCI.RHF,
# SMI.rhf)%>%
#   collect()%>%
#   mutate(BCI = BCI.RHF,
#          SMI = SMI.rhf)%>%
#  collect()%>%
# dplyr::select(-BCI.RHF, -SMI.rhf)
# 
# write.csv(cgsfall, "cgsfall-composition.csv")
# # 

#Ground squirrel spring
# cgsspr <- data %>%
#   mutate(fatsd = fat.sd,
#          leansd = lean.sd)%>%
#   collect()%>%
#   dplyr::select(sex,
# squirrel_id,
# wgt,
# rhf.avg,
# zyg.avg,
# fat,
# fat.sd,
# lean,
# lean.sd,
# season,
# leanp,
# fatp,
# logwgt,
# logrhf,
# logzyg,
# loglean,
# logfat,
# wgt.scl,
# rhf.scl,
# zyg.scl,
# lean.scl,
# fat.scl,
# BCI.RHF,
# SMI.rhf)%>%
#    collect()%>%
#    mutate(BCI = BCI.RHF,
#           SMI = SMI.rhf)%>%
#   collect()%>%
#  dplyr::select(-BCI.RHF, -SMI.rhf)
# 
#  write.csv(cgsspr, "cgsspr-composition.csv")

```


#MODELS

#Lean/fat by BCI/wgt  
#Analyzing BCIs against fat, lean. All data
```{r}
#BCI ZYG 

#Fat-BCIZYG - RSQ, BTPD <<--- note of which species are to be calculated with which BCI. Here, red squirrels and prairie dogs are noted for the ZYG BCI
fatzygm1 <- lm(fat ~
                 BCI.ZYG*sex, 
               data = data)
fatzygm1.resid1 = resid(fatzygm1)
summary(fatzygm1)
AICc(fatzygm1)


#Lean-BCIZYG - RSQ, BTPD 
leanzygm1 <- lm(lean ~
                 BCI.ZYG*sex, 
               data = data)
leanzygm1.resid1 = resid(leanzygm1)
summary(leanzygm1)
AICc(leanzygm1)

visreg(fatzygm1, "BCI.ZYG", by = "sex")
visreg(leanzygm1, "BCI.ZYG", by = "sex")



#BCI RHF 

#Fat-BCIRHF- CGS <<--- Columbian ground squirrels should use RHF BCI for both seasons
fatrhfm1 <- lm(fat ~
                 BCI.RHF*sex, 
               data = data)
fatrhfm1.resid1 = resid(fatrhfm1)
summary(fatrhfm1)
AICc(fatrhfm1)


#Lean-BCIRHF- CGS
leanrhfm1 <- lm(lean ~
                 BCI.RHF*sex, 
               data = data)
leanrhfm1.resid1 = resid(leanrhfm1)
summary(leanrhfm1)
AICc(leanrhfm1)

visreg(fatrhfm1, "BCI.RHF", by = "sex")
visreg(leanrhfm1, "BCI.RHF", by = "sex")



#WGT

#Fat- RSQ, BTPD, CGS
fatwgtm1 <- lm(fat ~
                 wgt*sex, 
               data = data)
fatwgtm1.resid1 = resid(fatwgtm1)
summary(fatwgtm1)
AICc(fatwgtm1)

#Lean-SQ, BTPD, CGS
leanwgtm1 <- lm(lean ~
                 wgt*sex, 
               data = data)
leanwgtm1.resid1 = resid(leanwgtm1)
summary(leanwgtm1)
AICc(leanwgtm1)


visreg(fatwgtm1, "wgt", by = "sex")
visreg(leanwgtm1, "wgt", by = "sex")

```


#Lean/Fat by BCI/mass
#Figures - BCI against lean/fat, mass against lean/fat
###NOTE: MUST change titles depending on what species you are running code for!
```{r, fig.width=5,fig.height=3}

#FAT -by ZYG - RSQ, BTPD
fig2azyg<-ggplot(data = data, 
              aes(x=BCI.ZYG, y=fat, fill = sex, shape = sex)) +
 geom_point(aes(shape=sex,color=sex),size=2,stroke = 1)+
  scale_shape_manual(values = c(16:17))+ 
  scale_fill_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue'))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
   ylab("Fat (g) ") + xlab("Body condition index (ZW index)") +
  #theme_void() +
  ggtitle('A) Red squirrel fat, pre-winter') +                                                    #<-------Change title to match data
   stat_smooth(method = "lm", se = TRUE, alpha = 0.5, aes(color = sex))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
      theme(#legend.position = "top",
   # legend.title = element_blank(),
  #  legend.key = element_blank(),
   # legend.text = element_text(size = 10),
    axis.ticks = element_blank(),
    axis.title=element_text(size=12), 
    axis.text = element_text(size = 12),
    strip.background = element_rect(fill="white", color = "black", size = 1),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill=NA, size = 1))  +
  annotate("text", x=Inf, y=-Inf, hjust=1.5,vjust=-0.5,label = bquote("R[adj.]^2 == 0.22"), parse = TRUE)    #<-------Change R2 value to match data
 # theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) #+
#  facet_wrap(~grid, ncol = 1)
fig2azyg

#FAT -by RHF - CGS
fig2arhf<-ggplot(data = data, 
              aes(x=BCI.RHF, y=fat, fill = sex, shape = sex)) +
 geom_point(aes(shape=sex,color=sex),size=2,stroke = 1)+
  scale_shape_manual(values = c(16:17))+ 
  scale_fill_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue'))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
   ylab("Fat (g) ") + xlab("Body condition index (RHF index)") +
  #theme_void() +
  ggtitle('A) Ground squirrel fat, pre-winter') +                                        #<-------Change title to match data
   stat_smooth(method = "lm", se = TRUE, alpha = 0.5, aes(color = sex))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
      theme(#legend.position = "top",
   # legend.title = element_blank(),
  #  legend.key = element_blank(),
   # legend.text = element_text(size = 10),
    axis.ticks = element_blank(),
    axis.title=element_text(size=12), 
    axis.text = element_text(size = 12),
    strip.background = element_rect(fill="white", color = "black", size = 1),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill=NA, size = 1)) +
  annotate("text", x=Inf, y=-Inf, hjust=1.5,vjust=-0.5,label = bquote("R[adj.]^2 == 0.80"), parse = TRUE)   #<-------Change R2 value to match data
 # theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) #+
#  facet_wrap(~grid, ncol = 1)
fig2arhf


#LEAN -by BCI  ZYG - RSQ, BTPD
fig2bzyg<-ggplot(data = data, 
              aes(x=BCI.ZYG, y=lean, fill = sex, shape = sex)) +
 geom_point(aes(shape=sex,color=sex),size=2,stroke = 1)+
  scale_shape_manual(values = c(16:17))+ 
  scale_fill_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue'))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
   ylab("Lean (g) ") + xlab("Body condition index (ZW index)") +
  #theme_void() +
  ggtitle('B) Red squirrel lean, pre-winter') +
   stat_smooth(method = "lm", se = TRUE, alpha = 0.5, aes(color = sex))+                            #<-------Change title to match data
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
      theme(#legend.position = "top",
   # legend.title = element_blank(),
  #  legend.key = element_blank(),
   # legend.text = element_text(size = 10),
    axis.ticks = element_blank(),
    axis.title=element_text(size=12), 
    axis.text = element_text(size = 12),
    strip.background = element_rect(fill="white", color = "black", size = 1),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill=NA, size = 1))  +
  annotate("text", x=Inf, y=-Inf, hjust=1.5,vjust=-0.5,label = bquote("R[adj.]^2 == 0.64"), parse = TRUE)   #<-------Change R2 value to match data 
 # theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) #+
#  facet_wrap(~grid, ncol = 1)
fig2bzyg



#LEAN -by BCI  RHF - CGS
fig2brhf<-ggplot(data = data, 
              aes(x=BCI.RHF, y=lean, fill = sex, shape = sex)) +
 geom_point(aes(shape=sex,color=sex),size=2,stroke = 1)+
  scale_shape_manual(values = c(16:17))+ 
  scale_fill_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue'))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
   ylab("Lean (g) ") + xlab("Body condition index (RHF index)") +
  #theme_void() +
  ggtitle('B) Ground squirrel lean, pre-winter') +                                                   #<-------Change title to match data
   stat_smooth(method = "lm", se = TRUE, alpha = 0.5, aes(color = sex))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
      theme(#legend.position = "top",
   # legend.title = element_blank(),
  #  legend.key = element_blank(),
   # legend.text = element_text(size = 10),
    axis.ticks = element_blank(),
    axis.title=element_text(size=12), 
    axis.text = element_text(size = 12),
    strip.background = element_rect(fill="white", color = "black", size = 1),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill=NA, size = 1))+
  annotate("text", x=Inf, y=-Inf, hjust=1.5,vjust=-0.5,label = bquote("R[adj.]^2 == 0.83"), parse = TRUE) #+ #<-------Change R2 value to match data
 # theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) #+
#  facet_wrap(~grid, ncol = 1)
fig2brhf




#FAT by WGT - RSQ, BTPD, CGS
fig2c <-ggplot(data = data, 
              aes(x=wgt, y = fat, fill = sex, shape = sex)) +
 geom_point(aes(shape=sex,color=sex),size=2,stroke = 1)+
  scale_shape_manual(values = c(16:17))+ 
  scale_fill_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue'))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
   ylab("Fat (g)") + xlab("Body mass (g)") +
  #theme_void() +
  ggtitle('C) Red squirrel fat, pre-winter') +                                                      #<-------Change title to match data
   stat_smooth(method = "lm", se = TRUE, alpha = 0.5, aes(color = sex))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
      theme(#legend.position = "top",
   # legend.title = element_blank(),
  #  legend.key = element_blank(),
   # legend.text = element_text(size = 10),
    axis.ticks = element_blank(),
    axis.title=element_text(size=12), 
    axis.text = element_text(size = 12),
    strip.background = element_rect(fill="white", color = "black", size = 1),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill=NA, size = 1)) +
  annotate("text", x=Inf, y=-Inf, hjust=1.5,vjust=-0.5,label = bquote("R[adj.]^2 == 0.22 "), parse = TRUE) #<-------Change R2 value to match data
 # theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) #+
#  facet_wrap(~grid, ncol = 1)
fig2c



#LEAN by WGT - RSQ, BTPD, CGS
fig2d <-ggplot(data = data, 
              aes(x=wgt, y = lean, fill = sex, shape = sex)) +
 geom_point(aes(shape=sex,color=sex),size=2,stroke = 1)+
  scale_shape_manual(values = c(16:17))+ 
  scale_fill_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue'))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
   ylab("Lean (g)") + xlab("Body mass (g)") +
  #theme_void() +
  ggtitle('D) Red squirrel lean, prewinter') +                                              #<-------Change title to match data
   stat_smooth(method = "lm", se = TRUE, alpha = 0.5, aes(color = sex))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
      theme(#legend.position = "top",
   # legend.title = element_blank(),
  #  legend.key = element_blank(),
   # legend.text = element_text(size = 10),
    axis.ticks = element_blank(),
    axis.title=element_text(size=12), 
    axis.text = element_text(size = 12),
    strip.background = element_rect(fill="white", color = "black", size = 1),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill=NA, size = 1)) +
  annotate("text", x=Inf, y=-Inf, hjust=1.5,vjust=-0.5,label = bquote("R[adj.]^2 == 0.78 "), parse = TRUE) #<-------Change R2 value to match data
 # theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) #+
#  facet_wrap(~grid, ncol = 1)
fig2d


```


#All data - read in and rbind 
```{r}
#This reads in the data written to CSV above - you MUST have run the code for each species (RSQ, BTPD, and CGS in both seasons) and have the data written to your working directory to proceed. 

rsq <- read.csv("rsq-composition.csv")
rsq$species <- "rsq"

btpd <- read.csv("btpd-composition.csv")
btpd$species <- "btpd"

cgsspr <- read.csv("cgsspr-composition.csv")
cgsspr$species <- "cgsspr"

cgsfall <- read.csv("cgsfall-composition.csv")
cgsfall$species <- "cgsfall"


#Red squirrel doesn't have season (all pre-winter) - add in column here. 

rsq$season <- "fall"
rsq <- rsq %>% dplyr::relocate(season, .after = lean.sd) #relocate season column


#Remove the unused BCI for each critter
rsq$BCI.RHF <- NULL
btpd$BCI.RHF <- NULL
cgsspr$BCI.ZYG<- NULL
cgsfall$BCI.ZYG <- NULL



#Rowbind all data
squirrels <- rbind(rsq, btpd, cgsspr, cgsfall)


squirrels <- squirrels %>%
    mutate(speciessex = paste(species, sex, sep = ''), 
           species = as.factor(species)) #Create species-sex variable


#Species variable as factor 

squirrels$species <- factor(x = squirrels$species, levels = c("rsq", "btpd", "cgsfall", "cgsspr"))


#Models 

fat.3sp <- lm (fat ~
                  BCI*species, 
                data = squirrels)
summary(fat.3sp)

plot(fat.3sp)
fatfit <- visreg(fat.3sp, "BCI", by = "species")

visreg(fat.3sp, "BCI", by = "species", gg=TRUE) + theme_classic()

lean.3sp <- lm (lean ~
                  BCI*species, 
                data = squirrels)
summary(lean.3sp)
leanfit <- visreg(lean.3sp, "BCI", by = "species")

visreg(lean.3sp, "BCI", by = "species", gg=TRUE) + theme_classic()


AICc(fat.3sp)
AICc(lean.3sp)


# New facet label names for species 
species.labs <- c("Red squirrels (pre-winter)", "Prairie dogs (pre-winter)", "Ground squirrels (pre-winter)", "Ground squirrels (spring)" )
names(species.labs) <- c("rsq", "btpd", "cgsfall", "cgsspr") 


#All three species

fig6 <-ggplot(filter(fatfit$fit), aes(BCI, visregFit))+
    geom_line(colour='black', size=1)+
    geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), 
                alpha=.3)+
    geom_point(data=filter(fatfit$res), 
               aes(BCI, visregRes), 
               size=1, alpha=.3, 
               position = "jitter") + 
    xlab('Body condition index (species-specfiic)')+
    ylab('Effect on fat (g)') + 
  theme_classic() + 
   facet_wrap(~species, labeller = labeller(species = species.labs),  ncol = 2, scales = "free")+ #scales = "free" allows axis ranges to vary in each panel 
  theme(text=element_text(size=15))
fig6

fig7 <-ggplot(filter(leanfit$fit), aes(BCI, visregFit))+
    geom_line(colour='black', size=1)+
    geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), 
                alpha=.3)+
    geom_point(data=filter(leanfit$res), 
               aes(BCI, visregRes), 
               size=1, alpha=.3, 
               position = "jitter") + 
    xlab('Body condition index (species-specfiic)')+
    ylab('Effect on lean (g)') + 

  theme_classic() + 
   facet_wrap(~species, labeller = labeller(species = species.labs),  ncol = 2, , scales = "free") + 
  theme(text=element_text(size=15))
fig7
```

#3 species but with body mass 
```{r}

squirrels %>%
  ggplot(aes(x = species, y = wgt.scl)) + 
  geom_boxplot() 


fat.3sp <- lm (fat ~
                  wgt.scl*species, 
                data = squirrels)
summary(fat.3sp)
fatfit <- visreg(fat.3sp, "wgt.scl", by = "species")


lean.3sp <- lm (lean ~
                  wgt.scl*species, 
                data = squirrels)
summary(lean.3sp)
leanfit <- visreg(lean.3sp, "wgt.scl", by = "species")

AICc(fat.3sp)
AICc(lean.3sp)

#All three species

fig6mass <-ggplot(filter(fatfit$fit), aes(wgt.scl, visregFit))+
    geom_line(colour='black', size=1)+
    geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), 
                alpha=.3)+
    geom_point(data=filter(fatfit$res), 
               aes(wgt.scl, visregRes), 
               size=1, alpha=.3, 
               position = "jitter") + 
    xlab('Body mass (scaled)')+
    ylab('Effect on fat (g)') + 
  theme_classic() + 
   facet_wrap(~species, labeller = labeller(species = species.labs),  ncol = 2, scales = "free")+ 
  theme(text=element_text(size=15))
fig6mass

fig7mass <-ggplot(filter(leanfit$fit), aes(wgt.scl, visregFit))+
    geom_line(colour='black', size=1)+
    geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), 
                alpha=.3)+
    geom_point(data=filter(leanfit$res), 
               aes(wgt.scl, visregRes), 
               size=1, alpha=.3, 
               position = "jitter") + 
    xlab('Body mass (scaled)')+
    ylab('Effect on lean (g)') + 
  theme_classic() + 
   facet_wrap(~species, labeller = labeller(species = species.labs),  ncol = 2, scales = "free")+ 
  theme(text=element_text(size=15))
fig7mass
```




#####




#Analyzing SMIs against fat, lean. All data
```{r}
#SMI ZYG 

#Fat-SMI zyg - RSQ, BTPD <<--- note of which species are to be calculated with which SMI Here, red squirrels and prairie dogs are noted for the ZYG SMI
fatzygm1 <- lm(fat ~
                 SMI.zyg*sex, 
               data = data)
fatzygm1.resid1 = resid(fatzygm1)
summary(fatzygm1)
AICc(fatzygm1)


#Lean-BCIZYG - RSQ, BTPD 
leanzygm1 <- lm(lean ~
                 SMI.zyg*sex, 
               data = data)
leanzygm1.resid1 = resid(leanzygm1)
summary(leanzygm1)
AICc(leanzygm1)

visreg(fatzygm1, "SMI.zyg", by = "sex")
visreg(leanzygm1, "SMI.zyg", by = "sex")



#SMI RHF 

#Fat-SMI RHF- CGS <<--- Columbian ground squirrels should use RHF SMI for both seasons
fatrhfm1 <- lm(fat ~
                 SMI.rhf*sex, 
               data = data)
fatrhfm1.resid1 = resid(fatrhfm1)
summary(fatrhfm1)
AICc(fatrhfm1)


#Lean-BCIRHF- CGS
leanrhfm1 <- lm(lean ~
                 SMI.rhf*sex, 
               data = data)
leanrhfm1.resid1 = resid(leanrhfm1)
summary(leanrhfm1)
AICc(leanrhfm1)

visreg(fatrhfm1, "SMI.rhf", by = "sex")
visreg(leanrhfm1, "SMI.rhf", by = "sex")



#WGT

#Fat- RSQ, BTPD, CGS
fatwgtm1 <- lm(fat ~
                 wgt*sex, 
               data = data)
fatwgtm1.resid1 = resid(fatwgtm1)
summary(fatwgtm1)
AICc(fatwgtm1)

#Lean-SQ, BTPD, CGS
leanwgtm1 <- lm(lean ~
                 wgt*sex, 
               data = data)
leanwgtm1.resid1 = resid(leanwgtm1)
summary(leanwgtm1)
AICc(leanwgtm1)


visreg(fatwgtm1, "wgt", by = "sex")
visreg(leanwgtm1, "wgt", by = "sex")

```


#Lean/Fat by SMI/mass
#Figures - SMI against lean/fat, mass against lean/fat
###NOTE: must change titles depending on what species you are running code for!
```{r, fig.width=5,fig.height=3}

#FAT -by ZYG - RSQ, BTPD
fig2azyg<-ggplot(data = data, 
              aes(x=SMI.zyg, y=fat, fill = sex, shape = sex)) +
 geom_point(aes(shape=sex,color=sex),size=2,stroke = 1)+
  scale_shape_manual(values = c(16:17))+ 
  scale_fill_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue'))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
   ylab("Fat (g) ") + xlab("Scaled mass index (ZW index)") +
  #theme_void() +
  ggtitle('A) Prairie dog fat, pre-winter') +
   stat_smooth(method = "lm", se = TRUE, alpha = 0.5, aes(color = sex))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
      theme(#legend.position = "top",
   # legend.title = element_blank(),
  #  legend.key = element_blank(),
   # legend.text = element_text(size = 10),
    axis.ticks = element_blank(),
    axis.title=element_text(size=12), 
    axis.text = element_text(size = 12),
    strip.background = element_rect(fill="white", color = "black", size = 1),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill=NA, size = 1)) #+
 # theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) #+
#  facet_wrap(~grid, ncol = 1)
fig2azyg

#FAT -by RHF - CGS
fig2arhf<-ggplot(data = data, 
              aes(x=SMI.rhf, y=fat, fill = sex, shape = sex)) +
 geom_point(aes(shape=sex,color=sex),size=2,stroke = 1)+
  scale_shape_manual(values = c(16:17))+ 
  scale_fill_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue'))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
   ylab("Fat (g) ") + xlab("Scaled massindex (RHF index)") +
  #theme_void() +
  ggtitle('A) Ground squirrel fat, spring') +
   stat_smooth(method = "lm", se = TRUE, alpha = 0.5, aes(color = sex))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
      theme(#legend.position = "top",
   # legend.title = element_blank(),
  #  legend.key = element_blank(),
   # legend.text = element_text(size = 10),
    axis.ticks = element_blank(),
    axis.title=element_text(size=12), 
    axis.text = element_text(size = 12),
    strip.background = element_rect(fill="white", color = "black", size = 1),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill=NA, size = 1)) #+
 # theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) #+
#  facet_wrap(~grid, ncol = 1)
fig2arhf


#LEAN -by SMI  ZYG - RSQ, BTPD
fig2bzyg<-ggplot(data = data, 
              aes(x=SMI.zyg, y=lean, fill = sex, shape = sex)) +
 geom_point(aes(shape=sex,color=sex),size=2,stroke = 1)+
  scale_shape_manual(values = c(16:17))+ 
  scale_fill_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue'))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
   ylab("Lean (g) ") + xlab("Scaled mass index (ZW index)") +
  #theme_void() +
  ggtitle('B) Prairie dog lean, pre-winter') +
   stat_smooth(method = "lm", se = TRUE, alpha = 0.5, aes(color = sex))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
      theme(#legend.position = "top",
   # legend.title = element_blank(),
  #  legend.key = element_blank(),
   # legend.text = element_text(size = 10),
    axis.ticks = element_blank(),
    axis.title=element_text(size=12), 
    axis.text = element_text(size = 12),
    strip.background = element_rect(fill="white", color = "black", size = 1),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill=NA, size = 1)) #+
 # theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) #+
#  facet_wrap(~grid, ncol = 1)
fig2bzyg

#LEAN -by SMI  RHF - CGS
fig2brhf<-ggplot(data = data, 
              aes(x=SMI.rhf, y=lean, fill = sex, shape = sex)) +
 geom_point(aes(shape=sex,color=sex),size=2,stroke = 1)+
  scale_shape_manual(values = c(16:17))+ 
  scale_fill_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue'))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
   ylab("Lean (g) ") + xlab("Scaled mass index (RHF index)") +
  #theme_void() +
  ggtitle('B) Ground squirrel lean, spring') +
   stat_smooth(method = "lm", se = TRUE, alpha = 0.5, aes(color = sex))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
      theme(#legend.position = "top",
   # legend.title = element_blank(),
  #  legend.key = element_blank(),
   # legend.text = element_text(size = 10),
    axis.ticks = element_blank(),
    axis.title=element_text(size=12), 
    axis.text = element_text(size = 12),
    strip.background = element_rect(fill="white", color = "black", size = 1),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill=NA, size = 1)) #+
 # theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) #+
#  facet_wrap(~grid, ncol = 1)
fig2brhf




#FAT by WGT - RSQ, BTPD, CGS
fig2c <-ggplot(data = data, 
              aes(x=wgt, y = fat, fill = sex, shape = sex)) +
 geom_point(aes(shape=sex,color=sex),size=2,stroke = 1)+
  scale_shape_manual(values = c(16:17))+ 
  scale_fill_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue'))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
   ylab("Fat (g)") + xlab("Body mass (g)") +
  #theme_void() +
  ggtitle('C) Ground squirrel fat, spring') +
   stat_smooth(method = "lm", se = TRUE, alpha = 0.5, aes(color = sex))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
      theme(#legend.position = "top",
   # legend.title = element_blank(),
  #  legend.key = element_blank(),
   # legend.text = element_text(size = 10),
    axis.ticks = element_blank(),
    axis.title=element_text(size=12), 
    axis.text = element_text(size = 12),
    strip.background = element_rect(fill="white", color = "black", size = 1),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill=NA, size = 1)) #+
 # theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) #+
#  facet_wrap(~grid, ncol = 1)
fig2c



#LEAN by WGT - RSQ, BTPD, CGS
fig2d <-ggplot(data = data, 
              aes(x=wgt, y = lean, fill = sex, shape = sex)) +
 geom_point(aes(shape=sex,color=sex),size=2,stroke = 1)+
  scale_shape_manual(values = c(16:17))+ 
  scale_fill_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue'))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
   ylab("Lean (g)") + xlab("Body mass (g)") +
  #theme_void() +
  ggtitle('D) Ground squirrel lean, spring') +
   stat_smooth(method = "lm", se = TRUE, alpha = 0.5, aes(color = sex))+
  scale_color_manual(values=c('#E69F00','#0072b2','darkorchid4','turquoise4','green','blue')) +
      theme(#legend.position = "top",
   # legend.title = element_blank(),
  #  legend.key = element_blank(),
   # legend.text = element_text(size = 10),
    axis.ticks = element_blank(),
    axis.title=element_text(size=12), 
    axis.text = element_text(size = 12),
    strip.background = element_rect(fill="white", color = "black", size = 1),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill=NA, size = 1)) #+
 # theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) #+
#  facet_wrap(~grid, ncol = 1)
fig2d



```


#All data - read in and rbind 
```{r}
#This reads in the data written to CSV above - you MUST have run the code for each species (RSQ, BTPD, and CGS in both seasons) and have the data written to your working directory to proceed. 

rsq <- read.csv("rsq-composition.csv")
rsq$species <- "rsq"

btpd <- read.csv("btpd-composition.csv")
btpd$species <- "btpd"

cgsspr <- read.csv("cgsspr-composition.csv")
cgsspr$species <- "cgsspr"

cgsfall <- read.csv("cgsfall-composition.csv")
cgsfall$species <- "cgsfall"


#Red squirrel doesn't have season (all pre-winter) - add in column here. 

rsq$season <- "fall"
rsq <- rsq %>% dplyr::relocate(season, .after = lean.sd) #relocate season column


#Remove the unused BCI for each critter
rsq$SMI.rhf <- NULL
btpd$SMI.rhf <- NULL
cgsspr$SMI.zyg<- NULL
cgsfall$SMI.zyg <- NULL



#Rowbind all data
squirrels <- rbind(rsq, btpd, cgsspr, cgsfall)


squirrels <- squirrels %>%
    mutate(speciessex = paste(species, sex, sep = ''), 
           species = as.factor(species)) #Create species-sex variable

#Scale within species
squirrels <- squirrels %>%
  mutate(SMI.scl = scale(SMI, center = TRUE, scale = TRUE))


#Species variable as factor 

squirrels$species <- factor(x = squirrels$species, levels = c("rsq", "btpd", "cgsfall", "cgsspr"))



# New facet label names for species 
species.labs <- c("Red squirrels (pre-winter)", "Prairie dogs (pre-winter)", "Ground squirrels (pre-winter)", "Ground squirrels (spring)" )
names(species.labs) <- c("rsq", "btpd", "cgsfall", "cgsspr") 


#All three species

fig6 <-ggplot(filter(fatfit$fit), aes(SMI, visregFit))+
    geom_line(colour='black', size=1)+
    geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), 
                alpha=.3)+
    geom_point(data=filter(fatfit$res), 
               aes(SMI, visregRes), 
               size=1, alpha=.3, 
               position = "jitter") + 
    xlab('Body condition index (species-specfiic)')+
    ylab('Effect on fat (g)') + 
  theme_classic() + 
   facet_wrap(~species, labeller = labeller(species = species.labs),  ncol = 2, scales = "free")+ #scales = "free" allows axis ranges to vary in each panel 
  theme(text=element_text(size=15))
fig6

fig7 <-ggplot(filter(leanfit$fit), aes(SMI, visregFit))+
    geom_line(colour='black', size=1)+
    geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), 
                alpha=.3)+
    geom_point(data=filter(leanfit$res), 
               aes(SMI, visregRes), 
               size=1, alpha=.3, 
               position = "jitter") + 
    xlab('Body condition index (species-specfiic)')+
    ylab('Effect on lean (g)') + 

  theme_classic() + 
   facet_wrap(~species, labeller = labeller(species = species.labs),  ncol = 2) + #, , scales = "free") + 
  theme(text=element_text(size=15))
fig7
```

#3 species but with body mass 
```{r}

squirrels %>%
  ggplot(aes(x = species, y = wgt.scl)) + 
  geom_boxplot() 


fat.3sp <- lm (fat ~
                  wgt.scl*species, 
                data = squirrels)
summary(fat.3sp)
fatfit <- visreg(fat.3sp, "wgt.scl", by = "species")


lean.3sp <- lm (lean ~
                  wgt.scl*species, 
                data = squirrels)
summary(lean.3sp)
leanfit <- visreg(lean.3sp, "wgt.scl", by = "species")

AICc(fat.3sp)
AICc(lean.3sp)

#All three species

fig6mass <-ggplot(filter(fatfit$fit), aes(wgt.scl, visregFit))+
    geom_line(colour='black', size=1)+
    geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), 
                alpha=.3)+
    geom_point(data=filter(fatfit$res), 
               aes(wgt.scl, visregRes), 
               size=1, alpha=.3, 
               position = "jitter") + 
    xlab('Body mass (scaled)')+
    ylab('Effect on fat (g)') + 
  theme_classic() + 
   facet_wrap(~species, labeller = labeller(species = species.labs),  ncol = 2, scales = "free")+ 
  theme(text=element_text(size=15))
fig6mass

fig7mass <-ggplot(filter(leanfit$fit), aes(wgt.scl, visregFit))+
    geom_line(colour='black', size=1)+
    geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr), 
                alpha=.3)+
    geom_point(data=filter(leanfit$res), 
               aes(wgt.scl, visregRes), 
               size=1, alpha=.3, 
               position = "jitter") + 
    xlab('Body mass (scaled)')+
    ylab('Effect on lean (g)') + 
  theme_classic() + 
   facet_wrap(~species, labeller = labeller(species = species.labs),  ncol = 2, scales = "free")+ 
  theme(text=element_text(size=15))
fig7mass
```

