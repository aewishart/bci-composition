---
title: "BCI-3species - Wishart et al"
output: pdf_document
---
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

#This code
#Load libraries
```{r}
library(tidyverse)
library(lubridate)
library(lme4)
library(lmerTest)
library(MCMCglmm)
library(ggplot2)
library(visreg)

select = dplyr::select #necessary as MASS also has a select function
filter = dplyr::filter
```


#Get Data 

#CGS
```{r}
cgs <- read.csv("BCIdata_cgs.csv") #Read in data
head(cgs)

#Remove individuals for whom we only have
cgs <- cgs %>%
  filter(!is.na(rhf.avg), !is.na(zyg.avg), !is.na(fat), !is.na(lean), age.class %in% "Adult")%>%
   mutate(sex = recode(sex, "M" = "male", "F" = "female")) %>%
  mutate(leanp = (lean/wgt)*100, fatp = (fat/wgt)*100)

#Split on season 

cgsspr <- cgs %>% #Spring only 
  filter(season %in% "spring") 
nrow(cgsspr)

cgsfall <- cgs %>% #Fall only 
  filter(season %in% "fall") 
nrow(cgsfall)


cgsfall %>% count(sex)
cgsspr %>% count(sex)

#Summary stats
cgsfall.summ.m <- cgsfall %>% #males
          filter(sex == "male")%>%
            summarize(meanwgt = mean(wgt),
            wgtsd = sd(wgt),
            meanrhf = mean(rhf.avg),
            sdrhf = sd(rhf.avg), 
           meanzyg = mean(zyg.avg), 
           sdzyg = sd(zyg.avg), 
           meanfat = mean(fat), 
           sdfat = sd(fat), 
            meanlean = mean(lean), 
            sdlean = sd(lean),
           meanfatp = mean(fatp), 
           sdfatp = sd(fatp), 
            meanleanp = mean(leanp), 
            sdleanp = sd(leanp))

cgsfall.summ.f <- cgsfall %>% #females
          filter(sex == "female")%>%
            summarize(meanwgt = mean(wgt),
            wgtsd = sd(wgt),
            meanrhf = mean(rhf.avg),
            sdrhf = sd(rhf.avg), 
           meanzyg = mean(zyg.avg), 
           sdzyg = sd(zyg.avg), 
        meanfat = mean(fat), 
           sdfat = sd(fat), 
            meanlean = mean(lean), 
            sdlean = sd(lean),
           meanfatp = mean(fatp), 
           sdfatp = sd(fatp), 
            meanleanp = mean(leanp), 
            sdleanp = sd(leanp))


cgsspr.summ.m <- cgsspr %>% #males
          filter(sex == "male")%>%
            summarize(meanwgt = mean(wgt),
            wgtsd = sd(wgt),
            meanrhf = mean(rhf.avg),
            sdrhf = sd(rhf.avg), 
           meanzyg = mean(zyg.avg), 
           sdzyg = sd(zyg.avg), 
           meanfat = mean(fat), 
           sdfat = sd(fat), 
            meanlean = mean(lean), 
            sdlean = sd(lean),
           meanfatp = mean(fatp), 
           sdfatp = sd(fatp), 
            meanleanp = mean(leanp), 
            sdleanp = sd(leanp))

cgsspr.summ.f <- cgsspr %>% #females
          filter(sex == "female")%>%
            summarize(meanwgt = mean(wgt),
            wgtsd = sd(wgt),
            meanrhf = mean(rhf.avg),
            sdrhf = sd(rhf.avg), 
           meanzyg = mean(zyg.avg), 
           sdzyg = sd(zyg.avg), 
            meanfat = mean(fat), 
           sdfat = sd(fat), 
            meanlean = mean(lean), 
            sdlean = sd(lean),
           meanfatp = mean(fatp), 
           sdfatp = sd(fatp), 
            meanleanp = mean(leanp), 
            sdleanp = sd(leanp))


cgsall <- rbind(cgsspr, cgsfall)
```

#Red Squirrel
```{r}
#Limit red squirrel data to FALL ONLY. Get means and SD
rsq <- read.csv("BCIdata_rsq.csv") #Read in data


rsq <- rsq %>% #Fall only 
  filter(season %in% "fall") %>%
    dplyr::select(-X)
head(rsq)
nrow(rsq)
#71 rows


#Remove individuals for whom we only have
rsq <- rsq %>%
  filter(!is.na(rhf.avg), !is.na(zyg.avg), !is.na(fat), !is.na(lean))%>%
   mutate(sex = recode(sex, "M" = "male", "F" = "female")) %>%
  mutate(leanp = (lean/wgt)*100, fatp = (fat/wgt)*100)

nrow(rsq)
n_distinct(rsq$squirrel.id)


#Summary stats
rsq %>% count(sex)

rsq.summ.m <- rsq %>% #males
          filter(sex == "male")%>%
            summarize(meanwgt = mean(wgt),
            wgtsd = sd(wgt),
            meanrhf = mean(rhf.avg),
            sdrhf = sd(rhf.avg), 
           meanzyg = mean(zyg.avg), 
           sdzyg = sd(zyg.avg), 
           meanfat = mean(fatp), 
           sdfat = sd(fatp), 
            meanlean = mean(lean), 
            sdlean = sd(lean),
           meanleanp = mean(leanp), 
            sdleanp = sd(leanp))
rsq.summ.m 

rsq.summ.f <- rsq  %>% #females
          filter(sex == "female")%>%
            summarize(meanwgt = mean(wgt),
            wgtsd = sd(wgt),
            meanrhf = mean(rhf.avg),
            sdrhf = sd(rhf.avg), 
           meanzyg = mean(zyg.avg), 
           sdzyg = sd(zyg.avg), 
           meanfat = mean(fatp), 
           sdfat = sd(fatp), 
            meanlean = mean(lean), 
            sdlean = sd(lean),
            meanleanp = mean(leanp), 
            sdleanp = sd(leanp))
rsq.summ.f

```



#Black-Tailed Prairie Dog
```{r}
btpd <- read.csv("BCIdata_btpd.csv") #Read in data

#btpd <- btpd %>% #Fall only 
 # filter(season %in% "fall") %>%
   # dplyr::select(-X)
head(btpd)
nrow(btpd)
#38 rows


#Remove individuals for whom we only have
btpd <- btpd %>%
  filter(!is.na(rhf.avg), !is.na(zyg.avg), !is.na(fat), !is.na(lean))%>%
   mutate(sex = recode(sex, "M" = "male", "F" = "female")) %>%
  mutate(leanp = (lean/wgt)*100, fatp = (fat/wgt)*100)

nrow(btpd)
n_distinct(btpd$squirrel.id)

#Summary stats
btpd %>% count(sex)

btpd.summ.m <- btpd %>% #males
          filter(sex == "male")%>%
            summarize(meanwgt = mean(wgt),
            wgtsd = sd(wgt),
            meanrhf = mean(rhf.avg),
            sdrhf = sd(rhf.avg), 
           meanzyg = mean(zyg.avg), 
           sdzyg = sd(zyg.avg), 
           meanfatp = mean(fatp), 
           sdfatp = sd(fatp), 
            meanleanp = mean(leanp), 
            sdleanp = sd(leanp),
           meanfat = mean(fat), 
           sdfat = sd(fat), 
            meanlean = mean(lean), 
            sdlean = sd(lean))
btpd.summ.m

btpd.summ.f <- btpd %>% #females
          filter(sex == "female")%>%
            summarize(meanwgt = mean(wgt),
            wgtsd = sd(wgt),
            meanrhf = mean(rhf.avg),
            sdrhf = sd(rhf.avg), 
           meanzyg = mean(zyg.avg), 
           sdzyg = sd(zyg.avg), 
           meanfat = mean(fatp), 
           sdfat = sd(fatp), 
            meanlean = mean(leanp), 
            sdleanp = sd(leanp),
           meanfat = mean(fat), 
           sdfat = sd(fat), 
            meanlean = mean(lean), 
            sdlean = sd(lean), 
            sdlean = sd(leanp))

btpd.summ.f

```

#Select data here
```{r}


#Red squirrel
#data <- rsq

#Black-tailed prairie dog
data <- btpd

#CGS
#data <- cgsfall

#data <- cgsspr

#data <- cgsall

```

#Look for outliers
```{r}
#Look for any outliers in morphometric data

#Body mass
par(mfrow=c(1, 2))
plot(data$wgt, data$rhf.avg, main="With Outliers", xlab="Body Mass", ylab="RHF Avg (mm)", pch="*", col="red", cex=2) 
abline(lm(wgt ~ rhf.avg, data=data), col="blue", lwd=3, lty=2)

outlier_values <- boxplot.stats(data$wgt)$out  # outlier values.
boxplot(data$wgt, main="Mass (g)", boxwex=0.1)
mtext(paste("Outliers: ", paste(outlier_values, collapse=", ")), cex=0.6)

print(outlier_values)#RHF Outliers to remove - look at each of these cases individually and give explanation to keep/remove

#Removing Outliers - example code from red squirrel data. This can be done for any variable that needs outliers removed. 

#With red squirrel data, we remove the smallest individual, squirrel ID = 24827; a male caught only once (fall QMR scan) with no other records. Can't be sure whether it was a juvenile or adult. 

data <- data %>%
  filter(!squirrel.id %in% (24827)) #We remove this record based on the squirrel ID (squirrel.id = 24827) as there is only one record for that individual. You may have to change what variable you are filtering on (may be better to remove based off a value)


#RHF Avg
par(mfrow=c(1, 2))
plot(data$wgt, data$rhf.avg, main="With Outliers", xlab="Body Mass", ylab="RHF Avg (mm)", pch="*", col="red", cex=2) 
abline(lm(wgt ~ rhf.avg, data=data), col="blue", lwd=3, lty=2)

outlier_values <- boxplot.stats(data$rhf.avg)$out  # outlier values.
boxplot(data$rhf.avg, main="RHF Avg (mm)", boxwex=0.1)
mtext(paste("Outliers: ", paste(outlier_values, collapse=", ")), cex=0.6)

print(outlier_values)#RHF Outliers to remove 


#Zyg Avg
par(mfrow=c(1, 2))
plot(data$wgt, data$zyg.avg, main="With Outliers", xlab="Body Mass", ylab="Zyg Avg (mm)", pch="*", col="red", cex=2) 
abline(lm(wgt~ rhf.avg, data=data), col="blue", lwd=3, lty=2)

outlier_values <- boxplot.stats(data$zyg.avg)$out  # outlier values.
boxplot(data$zyg.avg, main="Zyg Avg (mm)", boxwex=0.1)
mtext(paste("Outliers: ", paste(outlier_values, collapse=", ")), cex=0.6)

print(outlier_values) #Zyg Outliers to potentially remove 


#Cooks Distance for RHF and mass
mod <- lm(rhf.avg ~ wgt, data=data)
cooksd <- cooks.distance(mod)

par(mfrow=c(1, 1))
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels


#Histogram of rhf and zyg to help visualize these data. 

hist(data$rhf.avg)
mean(data$rhf.avg, na.rm = TRUE)
sd(data$rhf.avg, na.rm = TRUE)

hist(data$zyg.avg)
mean(data$zyg.avg, na.rm = TRUE)
sd(data$zyg.avg, na.rm = TRUE)

hist(data$wgt)
mean(data$wgt, na.rm = TRUE)
sd(data$wgt, na.rm = TRUE)

#Relationship between zyg and mass
qplot(wgt, zyg.avg, 
      color = sex,
      xlab = "Body Mass (g)", 
      ylab = "Zyg width (mm)",
      data = data)+ 
     stat_smooth(method="lm", alpha=0.2) + 
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()

zm1 <- ggplot(data, aes(wgt, zyg.avg, colour=sex)) +
 geom_point() +  stat_smooth(method="lm", alpha=0.2) +
  scale_color_manual(values=c("orange2", "dodgerblue4", "red")) +
  xlab("Body Mass (g)") + ylab("Zygomatic width (mm)") +
 # ggtitle("Columbian Ground Squirrel") +
   theme_bw()
ggExtra::ggMarginal(zm1, type = "density", groupColour =TRUE, groupFill = TRUE)

zyg.m1 = lm(zyg.avg ~ wgt + sex, 
            data = data)
summary(zyg.m1) 

#Relationship between rhf and mass
qplot(wgt, rhf.avg, 
      color = sex,
      xlab = "Body Mass (g)", 
      ylab = "RHF length (mm)",
      data = data)+ 
     stat_smooth(method="lm", alpha=0.2) + 
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()

rm1 <- ggplot(data, aes(wgt, rhf.avg, colour=sex)) +
 geom_point() +  stat_smooth(method="lm", alpha=0.2) +
  scale_color_manual(values=c("orange2", "dodgerblue4", "red")) +
  xlab("Body Mass (g)") + ylab("Right hind foot length(mm)") +
 # ggtitle("Columbian Ground Squirrel") +
   theme_bw()
ggExtra::ggMarginal(rm1, type = "density", groupColour =TRUE, groupFill = TRUE)

rhf.m1 = lm(rhf.avg ~ wgt + sex, 
            data = data)
summary(rhf.m1) 

#Relationship between zyg and RHF
rzm1= lm(rhf.avg ~ zyg.avg + sex, 
            data = data)
summary(rhf.m1) 

rz1 <- ggplot(data, aes(rhf.avg, zyg.avg, colour=sex)) +
 geom_point() +  stat_smooth(method="lm", alpha=0.2) +
  scale_color_manual(values=c("orange2", "dodgerblue4", "red")) +
  xlab("Right hind foot length (mm)") + ylab("Zygomatic width (mm)") +
#  ggtitle("Red Squirrel") +
   theme_bw()
ggExtra::ggMarginal(rz1, type = "density", groupColour =TRUE, groupFill = TRUE)



```
#Summary data 
```{r}
#Mass

massbox1 <- ggplot(data, aes(sex, wgt))
massbox1 + geom_boxplot() +
               xlab("Sex") + ylab("Body Mass (g)")+
              theme(axis.text.x = element_text(size=15), axis.title = element_text(size = 18), axis.text.y = element_text(size=15))

wgt1 <- aov(wgt ~ sex, data = data)
summary(wgt1)



#RHF
mean(data$rhf.avg)
sd(data$rhf.avg)

mean(data$rhf.avg)
sd(data$rhf.avg)

rhfbox1 <- ggplot(data, aes(sex, rhf.avg))
rhfbox1 + geom_boxplot() +
               xlab("Sex") + ylab("RHF (mm)") + ggtitle("Fall") + 
              theme(axis.text.x = element_text(size=15), axis.title = element_text(size = 18), axis.text.y =
                      element_text(size=15))

rhf1 <- aov(rhf.avg ~ sex, data = data)
summary(rhf1)

qplot(wgt, rhf.avg,
      color = sex,
      xlab = "wgt (g)",
      ylab = "rhf.avg (mm)",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()


#Zyg
mean(data$zyg.avg)
sd(data$zyg.avg)



zygbox3<- ggplot(data, aes(sex, zyg.avg))
zygbox3 + geom_boxplot() +
               xlab("Sex") + ylab("Zyg (mm)")+ ggtitle("Fall") + 
              theme(axis.text.x = element_text(size=15), axis.title = element_text(size = 18), axis.text.y = element_text(size=15))

zyg1 <- aov(zyg.avg ~ sex, data = data)
summary(zyg1)
 
qplot(wgt, zyg.avg,
      color = sex,
      xlab = "wgt (g)",
      ylab = "zyg.avg (mm)",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()


#Fat
mean(data$fat)
sd(data$fat)

mean(data$fatp)
sd(data$fatp)


fatbox4<- ggplot(data, aes(sex, fat))
fatbox4 + geom_boxplot() +
               xlab("Sex") + ylab("Fat %")+  ggtitle("Fall") + 
              theme(axis.text.x = element_text(size=15), axis.title = element_text(size = 18), axis.text.y = element_text(size=15))

fat1 <- aov(fat ~ sex, data = data)
summary(fat1)

qplot(wgt, fat,
      color = sex,
      xlab = "wgt (g)",
      ylab = "fat (g)",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()

qplot(wgt, fatp,
      color = sex,
      xlab = "wgt (g)",
      ylab = "fat (%)",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()


#Lean
mean(data$lean)
sd(data$lean)

mean(data$leanp)
sd(data$leanp)

leanbox5<- ggplot(data, aes(sex, lean))
leanbox5 + geom_boxplot() +
               xlab("Sex") + ylab("Lean %")+ ggtitle("Fall") + 
              theme(axis.text.x = element_text(size=15), axis.title = element_text(size = 18), axis.text.y = element_text(size=15))

lean1 <- aov(lean ~ sex, data = data)
summary(lean1)

qplot(wgt, lean,
      color = sex,
      xlab = "wgt (g)",
      ylab = "lean (g)",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()


qplot(wgt, leanp,
      color = sex,
      xlab = "wgt (g)",
      ylab = "lean (%)",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()

#lookhere
#Relationship between lean and rhf
qplot(zyg.avg, fatp,
  #    color = sex,
      xlab = "rhf.avg",
      ylab = "lean %",
      main = "Spring",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()
```



#Split data by sex (and rejoin)

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

#Relationship between zyg and rhf ?
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

data <- rbind(data.m, data.f)

#Relationship between zyg and rhf ?
cor.test(data.m$rhf.scl, data.m$zyg.scl) #RSQ no 

```

#Zygomatic and RHF first look 
```{r}

#Zygo and mass

zm2 <- ggplot(data, aes(wgt.scl, zyg.scl, colour=sex)) +
 geom_point() +  stat_smooth(method="lm", alpha=0.2) +
  scale_color_manual(values=c("orange2", "dodgerblue4", "red")) +
  xlab("Body Mass (scaled") + ylab("Zygomatic width (scaled)") +
 # ggtitle("Columbian Ground Squirrel") +
   theme_bw()
ggExtra::ggMarginal(zm1, type = "density", groupColour =TRUE, groupFill = TRUE)

#Relationship between zyg and mass?
cor.test(data$wgt.scl, data$zyg.scl) #all
cor.test(data.m$wgt.scl, data.m$zyg.scl) #males
cor.test(data.f$wgt.scl, data.f$zyg.scl) #females


#RHF and mass

rz2 <- ggplot(data, aes(zyg.scl, rhf.scl, colour=sex)) +
 geom_point() +  stat_smooth(method="lm", alpha=0.2) +
  scale_color_manual(values=c("orange2", "dodgerblue4", "red")) +
  xlab("Right hind foot (scaled)") + ylab("Zygomatic width (scaled)") +
 # ggtitle("Columbian Ground Squirrel") +
   theme_bw()
ggExtra::ggMarginal(rz2, type = "density", groupColour =TRUE, groupFill = TRUE)

#Relationship between rhf and mass?
cor.test(data$wgt.scl, data$rhf.scl) #all
cor.test(data.m$wgt.scl, data.m$rhf.scl) #males
cor.test(data.f$wgt.scl, data.f$rhf.scl) #females

```

#Calculate BCI from skeletal measurements - PCA based. All data. 
```{r}

#Log (RHF)
qplot(rhf.scl, wgt.scl,
       colour = sex,
      xlab = "Right hind foot",
      ylab = "Body Mass",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()


#Log(Zyg)
qplot(zyg.scl, wgt.scl,
       colour = sex,
      xlab = "Zygomatic width",
      ylab = "Body Mass",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()


#Principal components analysis to collapse RHF and Zyg into an estimate of skeletal size

skel <- as.data.frame(cbind(data$rhf.scl, data$zyg.scl))
head(skel)

skel.pca <- prcomp(skel, center=TRUE, scale.=TRUE)
skel.pca #Look at standard deviations, rotations for each of the newly generated componenents 


summary(skel.pca)

plot(prcomp(skel, center = TRUE, scale = TRUE))

pc.score <- skel.pca$x #PC1 and PC2 put into an object pc.score that we then cbind to the dataset. 

data <- cbind(data, pc.score)

head(data)
```

#Look at how skeletel measurements correlate with PC1 - all data 
```{r}
#Relationship between logzyg and PCs

qplot(PC1, zyg.scl,
      color = sex,
      xlab = "PC1",
      ylab = "Zygomatic width (scaled)",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()


#Relationship between logrhf and PCs

qplot(PC1, rhf.scl,
      color = sex,
      xlab = "PC1",
      ylab = "Right hind foot (scaled)",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()



#Get residuals of body mass regressed on PCs. The residuals are what will function as the BCI.

lm1.PC1 <- lm(wgt.scl ~ PC1, data=data)
attributes(lm1.PC1)
summary(lm1.PC1)
visreg(lm1.PC1)

{plot((wgt.scl) ~ PC1, data= data, xlab= "PC1", ylab = "Body Mass (scaled)", main = "Regression of scaled body mass on PC1")
abline(lm(wgt.scl ~ PC1, data=data), lty=2)}




cor.test(data$PC1, data$wgt.scl)


# Check normality 

plot(lm1.PC1)
shapiro.test(lm1.PC1$residuals)


#Residuals from this regression of log body mass ~ PC1 will serve as the body condition index (BCI). We keep in PC2 

data$BCI.1 <- lm1.PC1$residuals 

#Check to see how body condition index changes with body mass, body size and various skeletal measurements. 

plot(BCI.1 ~ PC1, data=data, xlab="PC1", ylab="Regression residuals (BCI)")

plot(BCI.1 ~ logwgt, data=data, xlab="Body mass (g)", ylab="Regression residuals (BCI)")

```

#Calculate BCI from skeletal measurements - PCA based. Females
```{r}

#Log (RHF)
qplot(rhf.scl, wgt.scl,
       colour = sex,
      xlab = "Right hind foot (scaled)",
      ylab = "Body Mass (Scaled)",
      data = data.f)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()


#Log(Zyg)
qplot(zyg.scl, wgt.scl,
       colour = sex,
      xlab = "zyg.scl",
      ylab = "wgt.scl",
      data = data.f)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()


#Principal components analysis to collapse RHF and Zyg into an estimate of skeletal size

skel.f <- as.data.frame(cbind(data.f$rhf.scl, data.f$zyg.scl))
head(skel.f)

skel.pca.f <- prcomp(skel.f, center=TRUE, scale.=TRUE)
skel.pca.f #Look at standard deviations, rotations for each of the newly generated componenents 


summary(skel.pca.f)

plot(prcomp(skel.f, center = TRUE, scale = TRUE))

pc.score.f <- skel.pca.f$x #PC1 and PC2 put into an object pc.score that we then cbind to the dataset. 

data.f <- cbind(data.f, pc.score.f)

head(data.f)
```

#Look at how skeletel measurements correlate with PC1 - Females only
```{r}
#Relationship between logzyg and PCs

qplot(PC1, zyg.scl,
      color = sex,
      xlab = "PC1",
      ylab = "Zygomatic width (scaled)",
      data = data.f)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()


#Relationship between logrhf and PCs

qplot(PC1, rhf.scl,
      color = sex,
      xlab = "PC1",
      ylab = "Right hind foot length (scaled)",
      data = data.f)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()


#Get residuals of body mass regressed on PCs. The residuals are what will function as the BCI.

lm1.PC1 <- lm(wgt.scl ~ PC1, data=data.f)
attributes(lm1.PC1)
summary(lm1.PC1)
visreg(lm1.PC1)

{plot((wgt.scl) ~ PC1, data= data.f, xlab= "PC1", ylab = "Body Mass (scaled)", main = "Regression of scaled body mass on PC1")
abline(lm(logwgt ~ PC1, data=data.f), lty=2)}



cor.test(data.f$PC1, data.f$wgt.scl)

# Check normality 

plot(lm1.PC1)
shapiro.test(lm1.PC1$residuals)

data.f$BCI.1 <- lm1.PC1$residuals 


#Check to see how body condition index changes with body mass, body size and various skeletal measurements.

plot(BCI.1 ~ PC1, data=data.f, xlab="PC1", ylab="Regression residuals (BCI)")

plot(BCI.1 ~ wgt.scl, data=data.f, xlab="Body mass (scaled)", ylab="Regression residuals (BCI)")
plot(BCI.1 ~ rhf.scl, data=data.f, xlab="RHF scaled", ylab="Regression residuals (BCI)")
plot(BCI.1 ~ zyg.scl, data=data.f, xlab="Zyg scaled)", ylab="Regression residuals (BCI)")

```
#Calculate BCI from skeletal measurements - PCA based, Males
```{r}

#Log (RHF)
qplot(rhf.scl, wgt.scl,
       colour = sex,
      xlab = "Rigt hind foot (scaled)",
      ylab = "Body mass (scaled)",
      data = data.m)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()


#Log(Zyg)
qplot(zyg.scl, wgt.scl,
       colour = sex,
      xlab = "Zygomatic width (scaled)",
      ylab = "Body Mass (scaled)",
      data = data.m)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()


#Principal components analysis to collapse RHF and Zyg into an estimate of skeletal size

skel.m <- as.data.frame(cbind(data.m$rhf.scl, data.m$zyg.scl))
head(skel.m)

skel.pca.m <- prcomp(skel.m, center=TRUE, scale.=TRUE)
skel.pca.m #Look at standard deviations, rotations for each of the newly generated componenents 


summary(skel.pca.m)

plot(prcomp(skel.m, center = TRUE, scale = TRUE))

pc.score.m <- skel.pca.m$x #PC1 and PC2 put into an object pc.score that we then cbind to the dataset. 

data.m <- cbind(data.m, pc.score.m)

head(data.m)
```


#Look at how skeletel measurements correlate with PC1 - Males only
```{r}
#Relationship between logzyg and PCs

qplot(PC1, zyg.scl,
      color = sex,
      xlab = "PC1",
      ylab = "Zygomatic width (scaled)",
      data = data.m)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()


#Relationship between logrhf and PCs

qplot(PC1, rhf.scl,
      color = sex,
      xlab = "PC1",
      ylab = "Right hind foot (scaled)",
      data = data.m)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()



#Get residuals of body mass regressed on PCs. The residuals are what will function as the BCI.

lm1.PC1 <- lm(wgt.scl ~ PC1, data=data.m)
attributes(lm1.PC1)
summary(lm1.PC1)
visreg(lm1.PC1)

{plot((logwgt) ~ PC1, data= data.m, xlab= "PC1", ylab = "Body Mass (scaled)", main = "Regression of scaled body mass on PC1")
abline(lm(logwgt ~ PC1, data=data.m), lty=2)}



cor.test(data.m$PC1, data.m$wgt.scl)

# Check normality 

plot(lm1.PC1)
shapiro.test(lm1.PC1$residuals)

data.m$BCI.1 <- lm1.PC1$residuals 


#Check to see how body condition index changes with body mass, body size and various skeletal measurements.

plot(BCI.1 ~ PC1, data=data.m, xlab="PC1", ylab="Regression residuals (BCI)")

plot(BCI.1 ~ wgt.scl, data=data.m, xlab="Body mass (scaled)", ylab="Regression residuals (BCI)")
plot(BCI.1 ~ rhf.scl, data=data.m, xlab="Right hind foot (scaled)", ylab="Regression residuals (BCI)")
plot(BCI.1 ~ zyg.scl, data=data.m, xlab="Zygomatic width (scaled)", ylab="Regression residuals (BCI)")
```



#Calculate BCI from singular skeletal measurements - RHF or ZYG on BM resids, all data
```{r}
####
#RHF
####
rhfm1 <- lm(wgt.scl ~ rhf.scl, data=data)
attributes(rhfm1)
summary(rhfm1)
visreg(rhfm1)

{plot((wgt.scl) ~ rhf.scl, data= data, xlab= "Right hind foot (scaled)", ylab = "Body Mass (scaled)", main = "Regression of scaled body mass on scaled RHF")
abline(lm(wgt.scl ~ rhf.scl, data = data), lty=2)}


#Check normality 
plot(rhfm1)
shapiro.test(rhfm1$residuals)

#Residuals from this regression of log body mass ~ log RHF will serve as the body condition index RHF (BCI.RHF).

data$BCI.RHF <- rhfm1$residuals 


####
#ZYG
####
zygm1 <- lm(wgt.scl ~ zyg.scl, data=data)
attributes(zygm1)
summary(zygm1)
visreg(zygm1)

{plot((wgt.scl) ~ zyg.scl, data= data, xlab= "Zygomatic width (scaled)", ylab = "Body Mass (scaled)", main = "Regression of scaled body mass on scaled zyg")
abline(lm(wgt.scl ~ zyg.scl, data=data), lty=2)}


#Check normality 
plot(zygm1)
shapiro.test(zygm1$residuals)


#Residuals from this regression of log body mass ~ log RHF will serve as the body condition index RHF (BCI.RHF).

data$BCI.ZYG <- zygm1$residuals 

```

#Calculate BCI from singular skeletal measurements - RHF or ZYG on BM resids, Females only
```{r}

####
#RHF
####
rhfm1.f <- lm(wgt.scl ~ rhf.scl, data=data.f)
attributes(rhfm1.f)
summary(rhfm1.f)

{plot((wgt.scl) ~ rhf.scl, data= data.f, xlab= "Right hind foot (scaled)", ylab = " Body Mass (scaled)", main = "Regression of scaled body mass on scaled RHF - Females only")
abline(lm(logwgt ~ logrhf, data=data.f), lty=2)}


#Check normality 
plot(rhfm1.f)
shapiro.test(rhfm1.f$residuals)

#Residuals from this regression of log body mass ~ log RHF will serve as the body condition index RHF (BCI.RHF).

data.f$BCI.RHFbysex <- rhfm1.f$residuals 


####
#ZYG
####
zygm1.f <- lm(wgt.scl ~ zyg.scl, data=data.f)
attributes(zygm1.f)
summary(zygm1.f)

{plot((wgt.scl) ~ zyg.scl, data= data.f, xlab= "Zygomatic width (scaled)", ylab = "Body Mass (scaled)", main = "Regression of scaled body mass on scaled zygomatic width - Females only")
abline(lm(wgt.scl ~ zyg.scl, data=data.f), lty=2)}


#Check normality 
plot(zygm1.f)
shapiro.test(zygm1.f$residuals)

#Residuals from this regression of log body mass ~ log RHF will serve as the body condition index RHF (BCI.RHF).

data.f$BCI.ZYGbysex <- zygm1.f$residuals 

```

#Calculate BCI from singular skeletal measurements - RHF or ZYG on BM resids, Males only
```{r}

####
#RHF
####
rhfm1.m <- lm(wgt.scl ~ rhf.scl, data=data.m)
attributes(rhfm1.m)
summary(rhfm1.m)

{plot((wgt.scl) ~ rhf.scl, data= data.m, xlab= "Right hind food (scaled)", ylab = "Body Mass (scaled)", main = "Regression of scaled body mass on scaled RHF - Males only")
abline(lm(wgt.scl ~ rhf.scl, data=data.m), lty=2)}


#Check normality 
plot(rhfm1.m)
shapiro.test(rhfm1.m$residuals)

#Residuals from this regression of log body mass ~ log RHF will serve as the body condition index RHF (BCI.RHF).

data.m$BCI.RHFbysex <- rhfm1.m$residuals 


####
#ZYG
####
zygm1.m <- lm(wgt.scl ~ zyg.scl, data=data.m)
attributes(zygm1.m)
summary(zygm1.m)

{plot((wgt.scl) ~ zyg.scl, data= data.m, xlab= "Zygomatic width (scaled)", ylab = "Body mass (scaled)", main = "Regression of scaled body mass on scaled zygomatic width - Males only")
abline(lm(wgt.scl ~ zyg.scl, data=data.m), lty=2)}


#Check normality 
plot(zygm1.m)
shapiro.test(zygm1.m$residuals)

#Residuals from this regression of log body mass ~ log RHF will serve as the body condition index RHF (BCI.RHF).

data.m$BCI.ZYGbysex <- zygm1.m$residuals 

```



#Join data
```{r} 
#Join the split-by-sex data together
data.bysex <- rbind(data.m, data.f)
head(data.bysex)

data.bysex <- data.bysex %>%
  mutate(BCI.1.bysex = BCI.1) %>%
  select(squirrel.id, BCI.1.bysex, BCI.RHFbysex, BCI.ZYGbysex)

data <- inner_join(data, data.bysex, by = c("squirrel.id"))



```



#Analyzing BCI against fat, lean. All data

```{r}
#FAT

qplot(BCI.1, fat.scl,
      color = sex,
      xlab = "BCI.1",
      ylab = "Fat (scaled)",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()

qplot(BCI.1.bysex, fat.scl,
     color = sex,
      xlab = "BCI.1.bysex",
      ylab = "Fat (scaled)",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()


qplot(BCI.1, fatp.scl,
      color = sex,
      xlab = "BCI.1",
      ylab = "Fat (% scaled)",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()

qplot(BCI.1.bysex, fatp.scl,
     color = sex,
      xlab = "BCI.1.bysex",
      ylab = "Fat (% scaled)",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()


fatm1 <- lm(fat.scl ~ BCI.1 + sex, data=data) #Regress residuals from body mass/PC1 on absolute fat
fat.resid1 = resid(fatm1)
summary(fatm1)


fatm2 <- lm(fat.scl ~ BCI.1.bysex + sex, data=data) #Regress residuals from body mass/PC1 on absolute fat
fat.resid2 = resid(fatm2)
summary(fatm2)

#Model comparisons

summary(fatm1) #BCI calculated for total
summary(fatm2) #BCI calculated by sex


cor.test(data$fat.scl, data$BCI.1) 
cor.test(data$fat.scl, data$BCI.1.bysex) 

AIC(fatm1) # by pop
AIC(fatm2) #by sex
AIC(fatm1) - AIC(fatm2) 



## LEAN

leanm1 <- lm(lean.scl ~ BCI.1, data=data) #Regress residuals from body mass/PC1 on absolute lean
lean.resid1 = resid(leanm1)
summary(leanm1)

leanm2 <- lm(lean.scl ~ BCI.1.bysex, data=data) #Regress residuals from body mass/PC1 on absolute lean
lean.resid2 = resid(leanm2)
summary(leanm2)
 
qplot(BCI.1, lean.scl,
       color = sex,
      xlab = "BCI.1",
      ylab = "Lean (scaled)",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()

qplot(BCI.1.bysex, lean.scl,
      color = sex,
      xlab = "BCI.1.bysex",
      ylab = "Lean (scaled)",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()

qplot(BCI.1, leanp.scl,
       color = sex,
      xlab = "BCI.1",
      ylab = "Lean (% scaled)",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()

qplot(BCI.1.bysex, leanp.scl,
      color = sex,
      xlab = "BCI.1.bysex",
      ylab = "Lean (% scaled)",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()

#Model comparisons

summary(leanm1) #BCI by pop
summary(leanm2) #BCI by sex

cor.test(data$loglean, data$BCI.1) 
cor.test(data$loglean, data$BCI.1.bysex) 

AIC(leanm1)#BCI by pop
AIC(leanm2)#BCI by sex

AIC(leanm1) - AIC(leanm2)



###RHF BCI
#FAT

qplot(BCI.RHF, fat.scl,
      color = sex,
      xlab = "BCI-RHF",
      ylab = "Fat (g scaled)",
      main = "Population-derived BCI - RHF only",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()

qplot(BCI.RHFbysex, fat.scl,
      color = sex,
      xlab = "BCI-RHF",
      ylab = "Fat (g scaled)",
      main = "By-Sex BCI - RHF only",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()

qplot(BCI.RHF, fatp.scl,
      color = sex,
      xlab = "BCI-RHF",
      ylab = "Fat (% scaled)",
      main = "Population-derived BCI - RHF only",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()

qplot(BCI.RHFbysex, fatp.scl,
      color = sex,
      xlab = "BCI-RHF",
      ylab = "Fat (g scaled)",
      main = "By-Sex BCI - RHF only",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()


fat.BCIRm1 <- lm(fat.scl ~ BCI.RHF + sex, data=data) #Regress residuals from body mass/PC1 on absolute fat
fat.resid1 = resid(fat.BCIRm1)
summary(fat.BCIRm1)

fat.BCIRm2 <- lm(logfat ~ BCI.RHFbysex + sex, data=data) #Regress residuals from body mass/PC1 on absolute fat
fat.resid2 = resid(fat.BCIRm2)
summary(fat.BCIRm2)

#Model comparisons

summary(fat.BCIRm1) #BCI calculated for total
summary(fat.BCIRm2) #BCI calculated by sex

cor.test(data$logfat, data$BCI.RHF) 
cor.test(data$logfat, data$BCI.RHFbysex) 

AIC(fat.BCIRm1) #by pop
AIC(fat.BCIRm2) #by sex

AIC(fat.BCIRm1) - AIC(fat.BCIRm2) 

#LEAN

qplot(BCI.RHF, lean.scl,
      color = sex,
      xlab = "BCI-RHF",
      ylab = "Lean (scaled)",
      main = "Population-derived BCI - RHF only",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()

qplot(BCI.RHFbysex, lean.scl,
      color = sex,
      xlab = "BCI-RHF",
      ylab = "Lean (scaled)",
      main = "By-Sex BCI - RHF only",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()


qplot(BCI.RHF, leanp.scl,
      color = sex,
      xlab = "BCI-RHF",
      ylab = "Lean (scaled)",
      main = "Population-derived BCI - RHF only",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()

qplot(BCI.RHFbysex, leanp.scl,
      color = sex,
      xlab = "BCI-RHF",
      ylab = "Lean (% scaled)",
      main = "By-Sex BCI - RHF only",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()


leanBCIRm1 <- lm(lean.scl ~ BCI.RHF, data=data) #Regress residuals from body mass/PC1 on absolute lean
lean.resid1 = resid(leanBCIRm1)
summary(leanBCIRm1)

lean.BCIRm2 <- lm(lean.scl ~ BCI.RHFbysex, data=data) #Regress residuals from body mass/PC1 on absolute lean
lean.resid2 = resid(lean.BCIRm2)
summary(lean.BCIRm2)

#Model comparisons

summary(leanBCIRm1) #BCI by pop
summary(lean.BCIRm2) #BCI by sex

cor.test(data$lean.scl, data$BCI.RHF) 
cor.test(data$lean.scl, data$BCI.RHFbysex) 

AIC(leanBCIRm1) #by pop
AIC(lean.BCIRm2) # by sex
AIC(leanBCIRm1) - AIC(lean.BCIRm2) 



### Zyg BCI
#FAT

qplot(BCI.ZYG, fat.scl,
      color = sex,
      xlab = "BCI-Zyg",
      ylab = "Fat (g scaled)",
      main = "Population-derived BCI - Zyg only",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()

qplot(BCI.ZYGbysex, fat.scl,
      color = sex,
      xlab = "BCI-Zyg",
      ylab = "Fat (g scaled)",
      main = "By-Sex BCI - Zyg only",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()

qplot(BCI.ZYG, fatp.scl,
      color = sex,
      xlab = "BCI-Zyg",
      ylab = "Fat (% scaled)",
      main = "Population-derived BCI - Zyg only",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()

qplot(BCI.ZYGbysex, fatp.scl,
      color = sex,
      xlab = "BCI-Zyg",
      ylab = "Fat (% scaled)",
      main = "By-Sex BCI - Zyg only",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()


fat.BCIZm1 <- lm(fat.scl ~ BCI.ZYG + sex, data=data) #Regress residuals from body mass/PC1 on absolute fat
fat.resid1 = resid(fat.BCIZm1)
summary(fat.BCIZm1)

fat.BCIZm2 <- lm(fat.scl ~ BCI.ZYGbysex + sex, data=data) #Regress residuals from body mass/PC1 on absolute fat
fat.resid2 = resid(fat.BCIZm2)
summary(fat.BCIZm2)

#Model comparisons

summary(fat.BCIZm1) #BCI calculated for total
summary(fat.BCIZm2) #BCI calculated by sex

cor.test(data$fat.scl, data$BCI.ZYG) 
cor.test(data$fat.scl, data$BCI.ZYGbysex) 

AIC(fat.BCIZm1) #by pop
AIC(fat.BCIZm2) #by sex

AIC(fat.BCIZm1) - AIC(fat.BCIZm2) 


#LEAN

qplot(BCI.ZYG, lean.scl,
      color = sex,
      xlab = "BCI-Zyg",
      ylab = "Lean (g scaled)",
      main = "Population-derived BCI - Zyg only",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()

qplot(BCI.ZYGbysex, lean.scl,
      color = sex,
      xlab = "BCI-Zyg",
      ylab = "Lean (g scaled)",
      main = "By-Sex BCI - ZYG only",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()

qplot(BCI.ZYG, leanp.scl,
      color = sex,
      xlab = "BCI-Zyg",
      ylab = "Lean (% scaled)",
      main = "Population-derived BCI - Zyg only",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()

qplot(BCI.ZYGbysex, leanp.scl,
      color = sex,
      xlab = "BCI-Zyg",
      ylab = "Lean (% scaled)",
      main = "By-Sex BCI - ZYG only",
      data = data)+
     stat_smooth(method="lm", alpha=0.2) +
    scale_color_manual(values=c("orange2", "dodgerblue4", "red"))+
     theme_bw()

lean.BCIZm1 <- lm(lean.scl ~ BCI.ZYG + sex, data=data) #Regress residuals from body mass/PC1 on absolute lean
lean.resid1 = resid(lean.BCIZm1)
summary(lean.BCIZm1)


lean.BCIZm2 <- lm(lean.scl ~ BCI.ZYGbysex + sex, data=data) #Regress residuals from body mass/PC1 on absolute lean
lean.resid2 = resid(lean.BCIZm2)
summary(lean.BCIZm2)

lean.BCIZm2.2 <- lm(lean.scl ~ BCI.ZYGbysex, data=data) 


#Model comparisons

summary(lean.BCIZm1) #BCI calculated for total
summary(lean.BCIZm2) #BCI calculated by sex

cor.test(data$lean.scl, data$BCI.ZYG) 
cor.test(data$lean.scl, data$BCI.ZYGbysex) 

AIC(lean.BCIZm1) # by pop
AIC(lean.BCIZm2) #by sex
AIC(lean.BCIZm1) - AIC(lean.BCIZm2)


```
