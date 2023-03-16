###This contains code and data files to implement climate analyses in 

##Temperature contributes to host specialization of coffee wilt disease (Fusarium xylarioides) on arabica and robusta coffee crops

#Xiuhan Zhang1, Lily D. Peck1, Julie Flood2, Matthew J. Ryan2, Timothy G. Barraclough1,3

#1Department of Life Sciences, Imperial College London, Silwood Park Campus, Ascot, Berkshire SL5 7PY, UK
#2CABI, Bakeham Lane, Egham, Surrey TW20 9TY UK
#3Department of Biology, University of Oxford, 11a Mansfield Rd, Oxford OX1 3SZ, UK.

##Authors: Xiuhan Zhang and Timothy G. Barraclough (tim.barraclough@biology.ox.ac.uk)

library(rTPC)
library(nls.multstart)
library(broom)
library(readr)
library(raster)
library(geodata)
library(nlraa)
hosts<-c("Arabica","Robusta")
cols<-c('#BB5566','#004488')

setwd("Data/Climatic modelling")

#################################################################
#### 1.Data preparation
## bioclimatic variables of WorldClim data (1970-2000)
bio.data<-raster::getData('worldclim', var='bio', res=10, lon=5, lat=5)
## bioclimatic variables of CMIP6 data (2060-2080)
future.bio.data<-cmip6_world(model="HadGEM3-GC31-LL",var="bio",ssp="245",res=10,time="2061-2080",path=getwd())
future.bio.data<-stack(future.bio.data)
names(future.bio.data)<-names(bio.data)

## severity records of Oduor et al (2003)
severe_data <- read_csv("Severity.csv")

# use coordinates of villages to extract bioclimatic values
fx_coords = cbind(severe_data$Longitude,severe_data$Latitude)

##read in 2002 average bioclim variables, year of sampling of severity
##these were extracted from monthly worldclim historic data and calculated for 2002 averages
fx_bio<-read.csv("fx_bio.csv")

# combine severity data with bioclimatic values
fx_serv=cbind(severe_data,fx_bio)

#################################################################
#### 2.NLS modelling
###  2.1.Min monthly temperatures 
# plot recorded severity points in Oduor et al (2003)
plot(fx_serv$bio6,fx_serv$Severity,type="n",xlab="Minimum Temperature (°C)",ylab="Severity (%)")
points(fx_serv$bio6[fx_serv$Host=="Arabica"],fx_serv$Severity[fx_serv$Host=="Arabica"],col='#BB5566')
points(fx_serv$bio6[fx_serv$Host=="Robusta"],fx_serv$Severity[fx_serv$Host=="Robusta"],col='#004488')
legend("topleft", legend=hosts, col=cols,bty="n",lty=1,cex=1.2)

params<-NULL
summ.CI<- NULL
y_predict<- NULL
##   2.1.1.ARABICA
# set up starting values and limits for NLS for arabica population
start_vals = get_start_vals(fx_serv$bio6[fx_serv$Host=="Arabica"], fx_serv$Severity[fx_serv$Host=="Arabica"], model_name = 'sharpeschoolhigh_1981')
low_lims = get_lower_lims(fx_serv$bio6[fx_serv$Host=="Arabica"], fx_serv$Severity[fx_serv$Host=="Arabica"], model_name = 'sharpeschoolhigh_1981')
upper_lims = get_upper_lims(fx_serv$bio6[fx_serv$Host=="Arabica"], fx_serv$Severity[fx_serv$Host=="Arabica"], model_name = 'sharpeschoolhigh_1981')
# run NLS for arabica population
fit_ara_minT <- nls_multstart(Severity~sharpeschoolhigh_1981(temp = bio6, r_tref,e,eh,th, tref = 15),
                     data = data.frame(fx_serv[fx_serv$Host=="Arabica",]),
                     iter = 500,
                     start_lower = start_vals - 10,
                     start_upper = start_vals + 10,
                     lower = low_lims,
                     upper = upper_lims,
                     supp_errors = 'Y')
# predict fitted values for smooth range of temperatures between chosen range
new_data <- data.frame(bio6 = seq(0, 20, 0.25))
y_predict=rbind(y_predict,predict(fit_ara_minT, new_data))
params<-rbind(params,calc_params(fit_ara_minT))
simulated.values<-simulate_nls(fit_ara_minT,100,3,data=data.frame(fx_serv[fx_serv$Host=="Arabica",]))
results<-NULL
predictions<-NULL
for (j in (1:100)) {
  e<-fx_serv[fx_serv$Host=="Arabica",]
  e$Severity<-simulated.values[,j]
  
  start_vals <- get_start_vals(e$bio6, e$Severity, model_name = 'sharpeschoolhigh_1981')
  low_lims <- get_lower_lims(e$bio6, e$Severity, model_name = 'sharpeschoolhigh_1981')
  upper_lims <- get_upper_lims(e$bio6, e$Severity, model_name = 'sharpeschoolhigh_1981')
  upper_lims[4]<-35
  
  fit <- nls_multstart(Severity~sharpeschoolhigh_1981(temp = bio6, r_tref,e,eh,th, tref = 15),
                       data = e,
                       iter = 500,
                       start_lower = start_vals - 10,
                       start_upper = start_vals + 10,
                       lower = low_lims,
                       upper = upper_lims,
                       supp_errors = 'Y')
  results<-rbind(results,calc_params(fit))
  if(!is.null(fit)){
    predictions<-rbind(predictions,predict(fit, new_data))}
}
## Plot confidence intervals
pre.95<-apply(predictions,2,function(x) quantile(x,c(0.025,0.975),na.rm=T))
polygon(x = c(new_data[,1], rev(new_data[,1])),
        y = c(pre.95[1,],rev(pre.95[2,])),
        col =  adjustcolor(cols[1], alpha.f = 0.10), border = NA)
CI.95<-apply(results,2,function(x) quantile(x,c(0.025,0.975),na.rm=T))
summ.CI<-rbind(summ.CI,apply(CI.95,2,function(x) paste(round(x,2),collapse="-")))

##   2.1.2.ROBUSTA
# set up starting values and limits for NLS for robusta population
start_vals <- get_start_vals(fx_serv$bio6[fx_serv$Host=="Robusta"], fx_serv$Severity[fx_serv$Host=="Robusta"], model_name = 'sharpeschoolhigh_1981')
low_lims <- get_lower_lims(fx_serv$bio6[fx_serv$Host=="Robusta"], fx_serv$Severity[fx_serv$Host=="Robusta"], model_name = 'sharpeschoolhigh_1981')
upper_lims <- get_upper_lims(fx_serv$bio6[fx_serv$Host=="Robusta"], fx_serv$Severity[fx_serv$Host=="Robusta"], model_name = 'sharpeschoolhigh_1981')
# run NLS for robusta population
fit_rob_minT <- nls_multstart(Severity~sharpeschoolhigh_1981(temp = bio6, r_tref,e,eh,th, tref = 15),
                     data = data.frame(fx_serv[fx_serv$Host=="Robusta",]),
                     iter = 500,
                     start_lower = start_vals - 10,
                     start_upper = start_vals + 10,
                     lower = low_lims,
                     upper = upper_lims,
                     supp_errors = 'Y')
# predict fitted values for smooth range of temperatures between chosen range
new_data <- data.frame(bio6 = seq(0, 20, 0.25))
y_predict=rbind(y_predict,predict(fit_rob_minT, new_data))
params<-rbind(params,calc_params(fit_rob_minT))
simulated.values<-simulate_nls(fit_rob_minT,100,3,data=data.frame(fx_serv[fx_serv$Host=="Robusta",]))
results<-NULL
predictions<-NULL
for (j in (1:100)) {
  e<-fx_serv[fx_serv$Host=="Robusta",]
  e$Severity<-simulated.values[,j]
  
  start_vals <- get_start_vals(e$bio6, e$Severity, model_name = 'sharpeschoolhigh_1981')
  low_lims <- get_lower_lims(e$bio6, e$Severity, model_name = 'sharpeschoolhigh_1981')
  upper_lims <- get_upper_lims(e$bio6, e$Severity, model_name = 'sharpeschoolhigh_1981')
  upper_lims[4]<-35
  
  fit <- nls_multstart(Severity~sharpeschoolhigh_1981(temp = bio6, r_tref,e,eh,th, tref = 15),
                       data = e,
                       iter = 500,
                       start_lower = start_vals - 10,
                       start_upper = start_vals + 10,
                       lower = low_lims,
                       upper = upper_lims,
                       supp_errors = 'Y')
  results<-rbind(results,calc_params(fit))
  if(!is.null(fit)){
    predictions<-rbind(predictions,predict(fit, new_data))}
}
## Plot confidence intervals
pre.95<-apply(predictions,2,function(x) quantile(x,c(0.025,0.975),na.rm=T))
polygon(x = c(new_data[,1], rev(new_data[,1])),
        y = c(pre.95[1,],rev(pre.95[2,])),
        col =  adjustcolor(cols[2], alpha.f = 0.10), border = NA)
CI.95<-apply(results,2,function(x) quantile(x,c(0.025,0.975),na.rm=T))
summ.CI<-rbind(summ.CI,apply(CI.95,2,function(x) paste(round(x,2),collapse="-")))
# calculate thermal performance values
lines(new_data[,1] ,y_predict[1,],col='#BB5566')
lines(new_data[,1] ,y_predict[2,],col='#004488')
params.min<-params
summ.CI.min<-summ.CI

######
###  2.2.Annual mean temperatures 
# plot recorded severity points in Oduor et al (2003)
plot(fx_serv$bio1,fx_serv$Severity,type="n",xlab="Annual mean Temperature (°C)",ylab="Severity (%)")
points(fx_serv$bio1[fx_serv$Host=="Arabica"],fx_serv$Severity[fx_serv$Host=="Arabica"],col='#BB5566')
points(fx_serv$bio1[fx_serv$Host=="Robusta"],fx_serv$Severity[fx_serv$Host=="Robusta"],col='#004488')
legend("topleft", legend=hosts, col=cols,lty=1,bty="n",cex=1.2)
params<-NULL
summ.CI<- NULL
y_predict<- NULL
##   2.2.1.ARABICA
# set up starting values and limits for NLS for arabica population
start_vals = get_start_vals(fx_serv$bio1[fx_serv$Host=="Arabica"], fx_serv$Severity[fx_serv$Host=="Arabica"], model_name = 'sharpeschoolhigh_1981')
low_lims = get_lower_lims(fx_serv$bio1[fx_serv$Host=="Arabica"], fx_serv$Severity[fx_serv$Host=="Arabica"], model_name = 'sharpeschoolhigh_1981')
upper_lims = get_upper_lims(fx_serv$bio1[fx_serv$Host=="Arabica"], fx_serv$Severity[fx_serv$Host=="Arabica"], model_name = 'sharpeschoolhigh_1981')
# run NLS for arabica population
fit_ara_meanT <- nls_multstart(Severity~sharpeschoolhigh_1981(temp = bio1, r_tref,e,eh,th, tref = 15),
                              data = data.frame(fx_serv[fx_serv$Host=="Arabica",]),
                              iter = 500,
                              start_lower = start_vals - 10,
                              start_upper = start_vals + 10,
                              lower = low_lims,
                              upper = upper_lims,
                              supp_errors = 'Y')
# predict fitted values for smooth range of temperatures between chosen range
new_data <- data.frame(bio1 = seq(0, 30, 0.25))
y_predict=rbind(y_predict,predict(fit_ara_meanT, new_data))
params<-rbind(params,calc_params(fit_ara_meanT))
simulated.values<-simulate_nls(fit_ara_meanT,100,3,data=data.frame(fx_serv[fx_serv$Host=="Arabica",]))
results<-NULL
predictions<-NULL
for (j in (1:100)) {
  e<-fx_serv[fx_serv$Host=="Arabica",]
  e$Severity<-simulated.values[,j]
  
  start_vals <- get_start_vals(e$bio1, e$Severity, model_name = 'sharpeschoolhigh_1981')
  low_lims <- get_lower_lims(e$bio1, e$Severity, model_name = 'sharpeschoolhigh_1981')
  upper_lims <- get_upper_lims(e$bio1, e$Severity, model_name = 'sharpeschoolhigh_1981')
  upper_lims[4]<-35
  
  fit <- nls_multstart(Severity~sharpeschoolhigh_1981(temp = bio1, r_tref,e,eh,th, tref = 15),
                       data = e,
                       iter = 500,
                       start_lower = start_vals - 10,
                       start_upper = start_vals + 10,
                       lower = low_lims,
                       upper = upper_lims,
                       supp_errors = 'Y')
  results<-rbind(results,calc_params(fit))
  if(!is.null(fit)){
    predictions<-rbind(predictions,predict(fit, new_data))}
}
## Plot confidence intervals
pre.95<-apply(predictions,2,function(x) quantile(x,c(0.025,0.975),na.rm=T))
polygon(x = c(new_data[,1], rev(new_data[,1])),
        y = c(pre.95[1,],rev(pre.95[2,])),
        col =  adjustcolor(cols[1], alpha.f = 0.10), border = NA)
CI.95<-apply(results,2,function(x) quantile(x,c(0.025,0.975),na.rm=T))
summ.CI<-rbind(summ.CI,apply(CI.95,2,function(x) paste(round(x,2),collapse="-")))

##   2.2.2.ROBUSTA
# set up starting values and limits for NLS for robusta population
start_vals <- get_start_vals(fx_serv$bio1[fx_serv$Host=="Robusta"], fx_serv$Severity[fx_serv$Host=="Robusta"], model_name = 'sharpeschoolhigh_1981')
low_lims <- get_lower_lims(fx_serv$bio1[fx_serv$Host=="Robusta"], fx_serv$Severity[fx_serv$Host=="Robusta"], model_name = 'sharpeschoolhigh_1981')
upper_lims <- get_upper_lims(fx_serv$bio1[fx_serv$Host=="Robusta"], fx_serv$Severity[fx_serv$Host=="Robusta"], model_name = 'sharpeschoolhigh_1981')
# run NLS for robusta population
fit_rob_meanT <- nls_multstart(Severity~sharpeschoolhigh_1981(temp = bio1, r_tref,e,eh,th, tref = 15),
                              data = data.frame(fx_serv[fx_serv$Host=="Robusta",]),
                              iter = 500,
                              start_lower = start_vals - 10,
                              start_upper = start_vals + 10,
                              lower = low_lims,
                              upper = upper_lims,
                              supp_errors = 'Y')
# predict fitted values for smooth range of temperatures between chosen range
new_data <- data.frame(bio1 = seq(0, 30, 0.25))
y_predict=rbind(y_predict,predict(fit_rob_meanT, new_data))
params<-rbind(params,calc_params(fit_rob_meanT))
simulated.values<-simulate_nls(fit_rob_meanT,100,3,data=data.frame(fx_serv[fx_serv$Host=="Robusta",]))
results<-NULL
predictions<-NULL
for (j in (1:100)) {
  e<-fx_serv[fx_serv$Host=="Robusta",]
  e$Severity<-simulated.values[,j]
  
  start_vals <- get_start_vals(e$bio1, e$Severity, model_name = 'sharpeschoolhigh_1981')
  low_lims <- get_lower_lims(e$bio1, e$Severity, model_name = 'sharpeschoolhigh_1981')
  upper_lims <- get_upper_lims(e$bio1, e$Severity, model_name = 'sharpeschoolhigh_1981')
  upper_lims[4]<-35
  
  fit <- nls_multstart(Severity~sharpeschoolhigh_1981(temp = bio1, r_tref,e,eh,th, tref = 15),
                       data = e,
                       iter = 500,
                       start_lower = start_vals - 10,
                       start_upper = start_vals + 10,
                       lower = low_lims,
                       upper = upper_lims,
                       supp_errors = 'Y')
  results<-rbind(results,calc_params(fit))
  if(!is.null(fit)){
    predictions<-rbind(predictions,predict(fit, new_data))}
}
## Plot confidence intervals
pre.95<-apply(predictions,2,function(x) quantile(x,c(0.025,0.975),na.rm=T))
polygon(x = c(new_data[,1], rev(new_data[,1])),
        y = c(pre.95[1,],rev(pre.95[2,])),
        col =  adjustcolor(cols[2], alpha.f = 0.10), border = NA)
CI.95<-apply(results,2,function(x) quantile(x,c(0.025,0.975),na.rm=T))
summ.CI<-rbind(summ.CI,apply(CI.95,2,function(x) paste(round(x,2),collapse="-")))
# calculate thermal performance values
lines(new_data[,1] ,y_predict[1,],col='#BB5566')
lines(new_data[,1] ,y_predict[2,],col='#004488')
params.mean<-params
summ.CI.mean<-summ.CI

######
###  2.3.Max monthly temperatures 
# plot recorded severity points in Oduor et al (2003)
plot(fx_serv$bio5,fx_serv$Severity,type="n",xlab="Maximum Temperature (°C)",ylab="Severity (%)")
points(fx_serv$bio5[fx_serv$Host=="Arabica"],fx_serv$Severity[fx_serv$Host=="Arabica"],col='#BB5566')
points(fx_serv$bio5[fx_serv$Host=="Robusta"],fx_serv$Severity[fx_serv$Host=="Robusta"],col='#004488')
legend("topleft", legend=hosts, col=cols,lty=1,bty="n",cex=1.2)
params<-NULL
summ.CI<- NULL
y_predict<- NULL
##   2.3.1.ARABICA
# set up starting values and limits for NLS for arabica population
start_vals = get_start_vals(fx_serv$bio5[fx_serv$Host=="Arabica"], fx_serv$Severity[fx_serv$Host=="Arabica"], model_name = 'sharpeschoolhigh_1981')
low_lims = get_lower_lims(fx_serv$bio5[fx_serv$Host=="Arabica"], fx_serv$Severity[fx_serv$Host=="Arabica"], model_name = 'sharpeschoolhigh_1981')
upper_lims = get_upper_lims(fx_serv$bio5[fx_serv$Host=="Arabica"], fx_serv$Severity[fx_serv$Host=="Arabica"], model_name = 'sharpeschoolhigh_1981')
# run NLS for arabica population
fit_ara_maxT <- nls_multstart(Severity~sharpeschoolhigh_1981(temp = bio5, r_tref,e,eh,th, tref = 15),
                               data = data.frame(fx_serv[fx_serv$Host=="Arabica",]),
                               iter = 500,
                               start_lower = start_vals - 10,
                               start_upper = start_vals + 10,
                               lower = low_lims,
                               upper = upper_lims,
                               supp_errors = 'Y')
# predict fitted values for smooth range of temperatures between chosen range
new_data <- data.frame(bio5 = seq(0, 40, 0.25))
y_predict=rbind(y_predict,predict(fit_ara_maxT, new_data))
params<-rbind(params,calc_params(fit_ara_maxT))
simulated.values<-simulate_nls(fit_ara_maxT,100,3,data=data.frame(fx_serv[fx_serv$Host=="Arabica",]))
results<-NULL
predictions<-NULL
for (j in (1:100)) {
  e<-fx_serv[fx_serv$Host=="Arabica",]
  e$Severity<-simulated.values[,j]
  
  start_vals <- get_start_vals(e$bio5, e$Severity, model_name = 'sharpeschoolhigh_1981')
  low_lims <- get_lower_lims(e$bio5, e$Severity, model_name = 'sharpeschoolhigh_1981')
  upper_lims <- get_upper_lims(e$bio5, e$Severity, model_name = 'sharpeschoolhigh_1981')
  upper_lims[4]<-35
  
  fit <- nls_multstart(Severity~sharpeschoolhigh_1981(temp = bio5, r_tref,e,eh,th, tref = 15),
                       data = e,
                       iter = 500,
                       start_lower = start_vals - 10,
                       start_upper = start_vals + 10,
                       lower = low_lims,
                       upper = upper_lims,
                       supp_errors = 'Y')
  results<-rbind(results,calc_params(fit))
  if(!is.null(fit)){
    predictions<-rbind(predictions,predict(fit, new_data))}
}
## Plot confidence intervals
pre.95<-apply(predictions,2,function(x) quantile(x,c(0.025,0.975),na.rm=T))
polygon(x = c(new_data[,1], rev(new_data[,1])),
        y = c(pre.95[1,],rev(pre.95[2,])),
        col =  adjustcolor(cols[1], alpha.f = 0.10), border = NA)
CI.95<-apply(results,2,function(x) quantile(x,c(0.025,0.975),na.rm=T))
summ.CI<-rbind(summ.CI,apply(CI.95,2,function(x) paste(round(x,2),collapse="-")))

##   2.3.2.ROBUSTA
# set up starting values and limits for NLS for robusta population
start_vals <- get_start_vals(fx_serv$bio5[fx_serv$Host=="Robusta"], fx_serv$Severity[fx_serv$Host=="Robusta"], model_name = 'sharpeschoolhigh_1981')
low_lims <- get_lower_lims(fx_serv$bio5[fx_serv$Host=="Robusta"], fx_serv$Severity[fx_serv$Host=="Robusta"], model_name = 'sharpeschoolhigh_1981')
upper_lims <- get_upper_lims(fx_serv$bio5[fx_serv$Host=="Robusta"], fx_serv$Severity[fx_serv$Host=="Robusta"], model_name = 'sharpeschoolhigh_1981')
# run NLS for robusta population
fit_rob_maxT <- nls_multstart(Severity~sharpeschoolhigh_1981(temp = bio5, r_tref,e,eh,th, tref = 15),
                               data = data.frame(fx_serv[fx_serv$Host=="Robusta",]),
                               iter = 500,
                               start_lower = start_vals - 10,
                               start_upper = start_vals + 10,
                               lower = low_lims,
                               upper = upper_lims,
                               supp_errors = 'Y')
# predict fitted values for smooth range of temperatures between chosen range
new_data <- data.frame(bio5 = seq(0, 40, 0.25))
y_predict=rbind(y_predict,predict(fit_rob_maxT, new_data))
params<-rbind(params,calc_params(fit_rob_maxT))
simulated.values<-simulate_nls(fit_rob_maxT,100,3,data=data.frame(fx_serv[fx_serv$Host=="Robusta",]))
results<-NULL
predictions<-NULL
for (j in (1:100)) {
  e<-fx_serv[fx_serv$Host=="Robusta",]
  e$Severity<-simulated.values[,j]
  
  start_vals <- get_start_vals(e$bio5, e$Severity, model_name = 'sharpeschoolhigh_1981')
  low_lims <- get_lower_lims(e$bio5, e$Severity, model_name = 'sharpeschoolhigh_1981')
  upper_lims <- get_upper_lims(e$bio5, e$Severity, model_name = 'sharpeschoolhigh_1981')
  upper_lims[4]<-35
  
  fit <- nls_multstart(Severity~sharpeschoolhigh_1981(temp = bio5, r_tref,e,eh,th, tref = 15),
                       data = e,
                       iter = 500,
                       start_lower = start_vals - 10,
                       start_upper = start_vals + 10,
                       lower = low_lims,
                       upper = upper_lims,
                       supp_errors = 'Y')
  results<-rbind(results,calc_params(fit))
  if(!is.null(fit)){
    predictions<-rbind(predictions,predict(fit, new_data))}
}
## Plot confidence intervals
pre.95<-apply(predictions,2,function(x) quantile(x,c(0.025,0.975),na.rm=T))
polygon(x = c(new_data[,1], rev(new_data[,1])),
        y = c(pre.95[1,],rev(pre.95[2,])),
        col =  adjustcolor(cols[2], alpha.f = 0.10), border = NA)
CI.95<-apply(results,2,function(x) quantile(x,c(0.025,0.975),na.rm=T))
summ.CI<-rbind(summ.CI,apply(CI.95,2,function(x) paste(round(x,2),collapse="-")))
# calculate thermal performance values
lines(new_data[,1] ,y_predict[1,],col='#BB5566')
lines(new_data[,1] ,y_predict[2,],col='#004488')
params.max<-params
summ.CI.max<-summ.CI

######
###  2.4.Shortcut section to thermal performance values
calc_params(fit_ara_minT)
calc_params(fit_rob_minT)
calc_params(fit_ara_meanT)
calc_params(fit_rob_meanT)
calc_params(fit_ara_maxT)
calc_params(fit_rob_maxT)

#################################################################
#### 3.Spatial predictions
###  3.1.Predict severity using corresponding temperature records from 1970-2000
# temperatures in bio.data are x10, so they need to be divided first
## Past (1970-2000) predictions
pg_ara_minT <- predict(bio.data$bio6/10, fit_ara_minT, ext=e,type="response")
pg_rob_minT <- predict(bio.data$bio6/10, fit_rob_minT, ext=e,type="response")
pg_ara_meanT <- predict(bio.data$bio1/10, fit_ara_meanT, ext=e,type="response")
pg_rob_meanT <- predict(bio.data$bio1/10, fit_rob_meanT, ext=e,type="response")
pg_ara_maxT <- predict(bio.data$bio5/10, fit_ara_maxT, ext=e,type="response")
pg_rob_maxT <- predict(bio.data$bio5/10, fit_rob_maxT, ext=e,type="response")
cellStats(pg_ara_meanT/100, stat='mean', na.rm=TRUE)
cellStats(pg_ara_maxT/100, stat='mean', na.rm=TRUE)
cellStats(pg_ara_minT/100, stat='mean', na.rm=TRUE)
cellStats(pg_rob_meanT/100, stat='mean', na.rm=TRUE)
cellStats(pg_rob_maxT/100, stat='mean', na.rm=TRUE)
cellStats(pg_rob_minT/100, stat='mean', na.rm=TRUE)
cellStats(pg_ara_meanT/100, stat='sd', na.rm=TRUE)
cellStats(pg_ara_maxT/100, stat='sd', na.rm=TRUE)
cellStats(pg_ara_minT/100, stat='sd', na.rm=TRUE)
cellStats(pg_rob_meanT/100, stat='sd', na.rm=TRUE)
cellStats(pg_rob_maxT/100, stat='sd', na.rm=TRUE)
cellStats(pg_rob_minT/100, stat='sd', na.rm=TRUE)

## Future (2060-2080) predictions
ft_ara_minT <- predict(future.bio.data$bio6, fit_ara_minT, ext=e,type="response")
ft_rob_minT <- predict(future.bio.data$bio6, fit_rob_minT, ext=e,type="response")
ft_ara_meanT <- predict(future.bio.data$bio1, fit_ara_meanT, ext=e,type="response")
ft_rob_meanT <- predict(future.bio.data$bio1, fit_rob_meanT, ext=e,type="response")
ft_ara_maxT <- predict(future.bio.data$bio5, fit_ara_maxT, ext=e,type="response")
ft_rob_maxT <- predict(future.bio.data$bio5, fit_rob_maxT, ext=e,type="response")
# test for significance of change
cellStats((ft_ara_meanT-pg_ara_meanT)/100, stat='mean', na.rm=TRUE)
t.test(ft_ara_meanT,pg_ara_meanT)
cellStats((ft_ara_maxT-pg_ara_maxT)/100, stat='mean', na.rm=TRUE)
t.test(ft_ara_maxT,pg_ara_maxT)
cellStats((ft_ara_minT-pg_ara_minT)/100, stat='mean', na.rm=TRUE)
t.test(ft_ara_minT,pg_ara_minT)
cellStats((ft_rob_meanT-pg_rob_meanT)/100, stat='mean', na.rm=TRUE)
t.test(ft_rob_meanT,pg_rob_meanT)
cellStats((ft_rob_maxT-pg_rob_maxT)/100, stat='mean', na.rm=TRUE)
t.test(ft_rob_maxT,pg_rob_maxT)
cellStats((ft_rob_minT-pg_rob_minT)/100, stat='mean', na.rm=TRUE)
t.test(ft_rob_minT,pg_rob_minT)
###  3.2.Plot spatial predictions
par(mfrow=c(1,1))
# arabica severity ~ min temperature
plot(pg_ara_minT/100, main='Arabica severity under 1970-2000 min temperature',ext=e,zlim=c(0,0.06))
plot(wrld_simpl, add=TRUE, border='dark grey') 
plot(ft_ara_minT/100-pg_ara_minT/100, main='Arabica severity under 2060-2080 min temperature',ext=e,col=bpy.colors(255),zlim=c(-0.05,0.05))
plot(wrld_simpl, add=TRUE, border='dark grey') 
# robusta severity ~ min temperature
plot(pg_rob_minT/100, main='Robusta severity under 1970-2000 min temperature',ext=e,zlim=c(0,0.5))
plot(wrld_simpl, add=TRUE, border='dark grey') 
plot(ft_rob_minT/100-pg_rob_minT/100, main='Robusta severity under 2060-2080 min temperature',ext=e,col=bpy.colors(255),zlim=c(-0.5,0.5))
plot(wrld_simpl, add=TRUE, border='dark grey') 
# arabica severity ~ mean temperature
plot(pg_ara_meanT/100, main='Arabica severity under 1970-2000 mean temperature',ext=e,zlim=c(0,0.06))
plot(wrld_simpl, add=TRUE, border='dark grey') 
plot(ft_ara_meanT/100-pg_ara_meanT/100, main='Arabica severity under 2060-2080 mean temperature',ext=e,col=bpy.colors(255),zlim=c(-0.05,0.05))
plot(wrld_simpl, add=TRUE, border='dark grey') 
# robusta severity ~ mean temperature
plot(pg_rob_meanT/100, main='Robusta severity under 1970-2000 mean temperature',ext=e,zlim=c(0,0.5))
plot(wrld_simpl, add=TRUE, border='dark grey') 
plot(ft_rob_meanT/100-pg_rob_meanT/100, main='Robusta severity under 2060-2080 mean temperature',ext=e,col=bpy.colors(255),zlim=c(-0.5,0.5))
plot(wrld_simpl, add=TRUE, border='dark grey') 
# arabica severity ~ max temperature
plot(pg_ara_maxT/100, main='Arabica severity under 1970-2000 max temperature',ext=e,zlim=c(0,0.06))
plot(wrld_simpl, add=TRUE, border='dark grey') 
plot(ft_ara_maxT/100-pg_ara_maxT/100, main='Arabica severity under 2060-2080 max temperature',ext=e,col=bpy.colors(255),zlim=c(-0.05,0.05))
plot(wrld_simpl, add=TRUE, border='dark grey') 
# robusta severity ~ max temperature
plot(pg_rob_maxT/100, main='Robusta severity under 1970-2000 max temperature',ext=e,zlim=c(0,0.5))
plot(wrld_simpl, add=TRUE, border='dark grey') 
plot(ft_rob_maxT/100-pg_rob_maxT/100, main='Robusta severity under 2060-2080 max temperature',ext=e,col=bpy.colors(255),zlim=c(-0.5,0.5))
plot(wrld_simpl, add=TRUE, border='dark grey') 
points(robu.coords,col = '#004488', pch = 4, cex = 0.5)

