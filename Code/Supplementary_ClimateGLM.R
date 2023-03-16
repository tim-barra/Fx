library(rTPC)
library(nls.multstart)
library(broom)
library(readr)
library(raster)
library(geodata)
library(caret)
library(rsq)
hosts<-c("Arabica","Robusta")
cols<-c('#BB5566','#004488')
#################################################################
#### 1.Data preparation (same as in Methods_ClimateNLS.R)
## bioclimatic variables of WorldClim data (1970-2000)
bio.data<-raster::getData('worldclim', var='bio', res=10, lon=5, lat=5)
## bioclimatic variables of CMIP6 data (2060-2080)
future.bio.data<-stack("Rdata/data/CMIP6/wc2.1_10m_bioc_HadGEM3-GC31-LL_ssp245_2061-2080.tif")
crs(future.bio.data) <- "+proj=utm +zone=1 +datum=WGS84"
## severity records of Oduor et al (2003)
severe_data <- read_csv("Rdata/data/Severity.csv")

# use coordinates of villages to extract bioclimatic values
fx_coords = cbind(severe_data$Longitude,severe_data$Latitude)
fx_bio=read.csv("fx_bio.csv")
# rename columns of bioclimatic dataset
colnames(fx_bio)<-c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")

# combine severity data with bioclimatic values
fx_serv=cbind(severe_data,fx_bio)
# convert severity data from percentage to decimal form
fx_serv$Severityp=fx_serv$Severity/100

#################################################################
#### 2.Binomial GLM modelling
## model parameters have been selected according to VIF and p-value
# GLM for arabica strains
Climtest_ara2 <- glm(Severityp ~ bio1+bio12+bio14,  family = binomial, data = na.omit(subset(fx_serv,Host=="Arabica")))
# inspect model parameters
summary(Climtest_ara2)
car::Anova(Climtest_ara2)
rsq(Climtest_ara2, adj=TRUE)

# GLM for robusta strains
Climtest_rob2 <- glm(Severityp ~ bio1+bio12+bio14,  family = binomial, data = na.omit(subset(fx_serv,Host=="Robusta")))
# inspect model parameters
summary(Climtest_rob2)
rsq(Climtest_rob2, adj=TRUE)
car::Anova(Climtest_rob2)

#################################################################
#### 3.Predicting past and future severity
# past (1970-2000) predictions
pg_ara_mix <- predict(bio.data2, Climtest_ara2, ext=e,type="response")
pg_rob_mix <- predict(bio.data2, Climtest_rob2, ext=e,type="response")
# future (2060-2080) predictions
ft_ara_mix <- predict(future.bio.data, Climtest_ara2, ext=e,type="response")
ft_rob_mix <- predict(future.bio.data, Climtest_rob2, ext=e,type="response")

cellStats(pg_ara_mix, stat='mean', na.rm=TRUE)
cellStats(pg_rob_mix, stat='mean', na.rm=TRUE)
cellStats(ft_ara_mix, stat='mean', na.rm=TRUE)
cellStats(ft_rob_mix, stat='mean', na.rm=TRUE)
cellStats(pg_ara_mix, stat='sd', na.rm=TRUE)
cellStats(pg_rob_mix, stat='sd', na.rm=TRUE)

# compare significance of change
cellStats(pg_ara_mix, stat='mean', na.rm=TRUE)-cellStats(ft_ara_mix, stat='mean', na.rm=TRUE)
cellStats(pg_rob_mix, stat='mean', na.rm=TRUE)-cellStats(ft_rob_mix, stat='mean', na.rm=TRUE)
t.test(ft_ara_mix,pg_ara_mix)
t.test(ft_rob_mix,pg_rob_mix)

# plot arabica severity predictions
plot(pg_ara_mix, main='Arabica severity under 1970-2000 climate',ext=e,zlim=c(0,0.85))
plot(wrld_simpl, add=TRUE, border='dark grey') 
plot(ft_ara_mix-pg_ara_mix, main='Arabica severity under 2060-2080 climate',ext=e,col=bpy.colors(255),zlim=c(-0.4,0.4))
plot(wrld_simpl, add=TRUE, border='dark grey') 
# plot robusta severity predictions
plot(pg_rob_mix, main='Robusta severity under 1970-2000 climate',ext=e,zlim=c(0,1))
plot(wrld_simpl, add=TRUE, border='dark grey') 
plot(ft_rob_mix-pg_rob_mix, main='Robusta severity under 2060-2080 climate',ext=e,col=bpy.colors(255),zlim=c(-0.6,0.6))
plot(wrld_simpl, add=TRUE, border='dark grey') 

