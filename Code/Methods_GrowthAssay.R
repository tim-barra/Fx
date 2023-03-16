###This contains code and data files to implement assay analyses in 

##Temperature contributes to host specialization of coffee wilt disease (Fusarium xylarioides) on arabica and robusta coffee crops

#Xiuhan Zhang1, Lily D. Peck1, Julie Flood2, Matthew J. Ryan2, Timothy G. Barraclough1,3

#1Department of Life Sciences, Imperial College London, Silwood Park Campus, Ascot, Berkshire SL5 7PY, UK
#2CABI, Bakeham Lane, Egham, Surrey TW20 9TY UK
#3Department of Biology, University of Oxford, 11a Mansfield Rd, Oxford OX1 3SZ, UK.

##Authors: Xiuhan Zhang and Timothy G. Barraclough (tim.barraclough@biology.ox.ac.uk)

library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(rsq)
library(rTPC)
library(nls.multstart)
library(broom)
library(nlraa)

###########################################
#####
# Specify the Fx hosts and colours for them
hosts<-c("Arabica","Robusta")
cols<-c('#BB5566','#004488')

##put 1 plot on one graph
par(mfrow=c(1,1))

###########################################
##1.Mean growth rates between 4-8 days
#####
##1.1.Data processing
# Read data for mean growth rates 
Growth_mean <- read.csv("Growth_mean.csv")
Growth_mean.stat<-Growth_mean
# Specify temperature and strain as factor for LMM analysis
Growth_mean.stat$Temperature=as.factor(Growth_mean.stat$Temperature)
#Growth_mean.stat$Temperature2[Growth_mean.stat$Temperature2==35]<-40
Growth_mean.stat$Strain=as.factor(Growth_mean.stat$Strain)

#####
##1.2.Run LMM for growth rates
# Perform LMM for mean growth rates vs host and temperature,
# round and strain included as random effects
lmem <- lmer(Mean_Growth ~ Host * Temperature + (1|Round)+(1|Strain), data = Growth_mean.stat) 
# Calculate mean and SE of estimates and adjusted R-squared of LMM
summary(lmem)
emmeans(lmem,list(pairwise~Host*Temperature),adjust="tukey")
rsq(lmem,adj=TRUE)
# Produce data frame of LMM results
LMMGrowth <- as.data.frame(emmeans(lmem,list(pairwise~Host*Temperature),adjust="tukey")[[1]])

#####
##1.3.Run NLS for growth rates 
# Plot the background for NLS curves
plot(c(10,40),range(Growth_mean$Mean_Growth),type="n",
     xlab="Temperature (°C)",ylab="Growth rate (mm/day)",
     cex.lab=1.2,las=1)
# Add legend of different population
legend("topleft", legend=hosts, col=cols,lty=1,bty="n",cex=1.2)

# Clear the dataframe for parameter outputs 
params<-NULL
# Loop through populations
for (i in (1:length(hosts))) {
  # Select population
  host<-hosts[i]
  # Extract data for just this population
  subs<-(Growth_mean$Host==host)
  
  # Extract sensible starting guesses for the parameters using Sharpes-Schoolfield model
  start_vals <- get_start_vals(Growth_mean$Temperature[subs], Growth_mean$Mean_Growth[subs], model_name = 'sharpeschoolhigh_1981')
  low_lims <- get_lower_lims(Growth_mean$Temperature[subs], Growth_mean$Mean_Growth[subs], model_name = 'sharpeschoolhigh_1981')
  upper_lims <- get_upper_lims(Growth_mean$Temperature[subs], Growth_mean$Mean_Growth[subs], model_name = 'sharpeschoolhigh_1981')
  upper_lims[4]<-35
 
  # Run nls fitting for equation of choice - this one fits the Sharpe-Schoolfield equation
  fit <- nls_multstart(Mean_Growth~sharpeschoolhigh_1981(temp = Temperature, r_tref,e,eh,th, tref = 15),
                       data = Growth_mean[subs,],
                       iter = 500,
                       start_lower = start_vals - 10,
                       start_upper = start_vals + 10,
                       lower = low_lims,
                       upper = upper_lims,
                       supp_errors = 'Y')
  
  # Predict fitted values for smooth range of temperatures between chosen range
  new_data <- data.frame(Temperature = seq(10, 40, 0.25))
  # Plot lines of predicted values into the empty plot
  lines(new_data[,1] ,predict(fit, new_data),col=cols[i])
  # Plot the points and error bars for the LMM predictions
  subs1<-(LMMGrowth$Host==host)
  points(as.numeric(as.character(LMMGrowth$Temperature[subs1])),LMMGrowth$emmean[subs1],col=cols[i],pch=16)
  arrows(as.numeric(as.character(LMMGrowth$Temperature[subs1])),LMMGrowth$emmean[subs1]+LMMGrowth$SE[subs1], as.numeric(as.character(LMMGrowth$Temperature[subs1])),LMMGrowth$emmean[subs1]-LMMGrowth$SE[subs1], code=3, angle=90, length=0.1,col=cols[i])
  # Record output parameter estimate values into object
  params<-rbind(params,calc_params(fit))
  
  #show data points
  #points(jitter(as.numeric(Growth_mean$Temperature[Growth_mean$Host==hosts[i]]),0.25),Growth_mean$Mean_Growth[Growth_mean$Host==hosts[i]],col=cols[i],pch=20,cex=0.5)

}
# Report parameter estimates for each population
params.growth<-params

##FIND CONFIDENCE LIMITS FOR PARAMETERS
##work out confidence intervals
#simulated.values<-simulate(lmem,100)
plot(c(10,40),range(Growth_mean$Mean_Growth),type="n",
     xlab="Temperature (°C)",ylab="Growth rate (mm/day)",
     cex.lab=1.2,las=1)
# Add legend of different population
legend("topright", legend=hosts, col=cols,lty=1,bty="n",cex=1.2)

# Clear the dataframe for parameter outputs 
params<-NULL
summ.CI<- NULL
y_predict<- NULL
# Loop through populations
for (i in (1:length(hosts))) {
  # Select population
  host<-hosts[i]
  # Extract data for just this population
  subs<-(Growth_mean$Host==host)
  
  # Extract sensible starting guesses for the parameters using Sharpes-Schoolfield model
  start_vals <- get_start_vals(Growth_mean$Temperature[subs], Growth_mean$Mean_Growth[subs], model_name = 'sharpeschoolhigh_1981')
  low_lims <- get_lower_lims(Growth_mean$Temperature[subs], Growth_mean$Mean_Growth[subs], model_name = 'sharpeschoolhigh_1981')
  upper_lims <- get_upper_lims(Growth_mean$Temperature[subs], Growth_mean$Mean_Growth[subs], model_name = 'sharpeschoolhigh_1981')
  upper_lims[4]<-35

 # Run nls fitting for equation of choice - this one fits the Sharpe-Schoolfield equation
  fit1 <- nls_multstart(Mean_Growth~sharpeschoolhigh_1981(temp = Temperature, r_tref,e,eh,th, tref = 15),
                       data = Growth_mean[subs,],
                       iter = 500,
                       start_lower = start_vals - 10,
                       start_upper = start_vals + 10,
                       lower = low_lims,
                       upper = upper_lims,
                       supp_errors = 'Y')
  new_data <- data.frame(Temperature = seq(10, 40, 0.25))
  # Plot lines of predicted values into the empty plot
  #lines(new_data[,1] ,predict(fit1, new_data),col=cols[i])
  y_predict=rbind(y_predict,predict(fit1, new_data))
  # Plot the points and error bars for the LMM predictions

  # Record output parameter estimate values into object
  params<-rbind(params,calc_params(fit1))
  simulated.values<-simulate_nls(fit1,100,3,data=Growth_mean[subs,])

  results<-NULL
  predictions<-NULL
  for (j in (1:100)) {
  	e<-Growth_mean[subs,]
  	e$Mean_Growth<-simulated.values[,j]
  
    start_vals <- get_start_vals(e$Temperature, e$Mean_Growth, model_name = 'sharpeschoolhigh_1981')
    low_lims <- get_lower_lims(e$Temperature, e$Mean_Growth, model_name = 'sharpeschoolhigh_1981')
    upper_lims <- get_upper_lims(e$Temperature, e$Mean_Growth, model_name = 'sharpeschoolhigh_1981')
    upper_lims[4]<-35
  
  	fit <- nls_multstart(Mean_Growth~sharpeschoolhigh_1981(temp = Temperature, r_tref,e,eh,th, tref = 15),
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
  pre.95<-apply(predictions,2,function(x) quantile(x,c(0.025,0.975),na.rm=T))
  polygon(x = c(new_data[,1], rev(new_data[,1])),
          y = c(pre.95[1,],rev(pre.95[2,])),
          col =  adjustcolor(cols[i], alpha.f = 0.10), border = NA)
	CI.95<-apply(results,2,function(x) quantile(x,c(0.025,0.975),na.rm=T))
	summ.CI<-rbind(summ.CI,apply(CI.95,2,function(x) paste(round(x,2),collapse="-")))
}
##Plot lines of model predictions
lines(new_data[,1] ,y_predict[1,],col='#BB5566')
lines(new_data[,1] ,y_predict[2,],col='#004488')
##Plot means and error bars
subs_ara<-(LMMGrowth$Host=="Arabica")
subs_rob<-(LMMGrowth$Host=="Robusta")
points(as.numeric(as.character(LMMGrowth$Temperature[subs_ara])),LMMGrowth$emmean[subs_ara],col='#BB5566',pch=16)
arrows(as.numeric(as.character(LMMGrowth$Temperature[subs_ara])),LMMGrowth$emmean[subs_ara]+LMMGrowth$SE[subs_ara], as.numeric(as.character(LMMGrowth$Temperature[subs_ara])),LMMGrowth$emmean[subs_ara]-LMMGrowth$SE[subs_ara], code=3, angle=90, length=0.1,col='#BB5566')
points(as.numeric(as.character(LMMGrowth$Temperature[subs_rob])),LMMGrowth$emmean[subs_rob],col='#004488',pch=16)
arrows(as.numeric(as.character(LMMGrowth$Temperature[subs_rob])),LMMGrowth$emmean[subs_rob]+LMMGrowth$SE[subs_rob], as.numeric(as.character(LMMGrowth$Temperature[subs_rob])),LMMGrowth$emmean[subs_rob]-LMMGrowth$SE[subs_rob], code=3, angle=90, length=0.1,col='#004488')
##95% CI
params.growth<-params
summ.CI.growth<-summ.CI

###########################################
##2.Spore density
#####
##2.1.Data processing
# Read data for unit spore concentration (calculated data,
# raw data available in Sporecount(raw_data).csv)
Sporesden <- read.csv("Sporesden.csv")
##convert density to 1000 count per mm2 by multiplying up to 10mL then divide by 1000
Sporesden$density<-10*Sporesden$density/1000
Sporesden.stat<-Sporesde
# Specify temperature and strain as factor for LMM analysis
Sporesden.stat$Temperature=as.factor(Sporesden$Temperature)
Sporesden.stat$Strain=as.factor(Sporesden$Strain)

#####
##2.2.Run LMM for unit spore concentration
# Perform LMM for unit spore conc vs host and temperature,
# round and strain included as random effects
lmem_spore=lmer(density ~ Host * Temperature + (1|Round)+(1|Strain), data = Sporesden.stat) 
# Calculate mean and SE of estimates and adjusted R-squared of LMM
summary(lmem_spore)
emmeans(lmem_spore,list(pairwise~Host*Temperature),adjust="tukey")
rsq(lmem_spore,adj=TRUE)
##anova to see if host X temperature interaction is significant overall
anova(lmem_spore)

##check result for just robusta, reported in results
lmem_spore_robusta=lmer(density ~ Temperature + (1|Round)+(1|Strain), data = Sporesden.stat, subset=Host=="Robusta") 
anova(lmem_spore_robusta)
# Produce data frame of LMM results
LMMSpore<-as.data.frame(emmeans(lmem_spore,list(pairwise~Host*Temperature),adjust="tukey",type="response")[[1]])

#####
##2.3.Run NLS for unit spore concentration
# Plot the background for NLS curves
plot(c(10,40),range(Sporesden$density),type="n",
     xlab="Temperature (°C)",ylab="Sporulation rate (1000/mm2)",
     cex.lab=1.2,las=1)
# Add legend of different population
legend("topright", legend=hosts, col=cols,lty=1,bty="n",cex=1.2)

# Clear the dataframe for parameter outputs 
params<-NULL
summ.CI<- NULL
y_predict<- NULL
# Loop through populations
for (i in (1:length(hosts))) {
   # Select population
  host<-hosts[i]
  # Extract data for just this population
  subs<-(Sporesden$Host==host)

# Extract sensible starting guesses for the parameters using Sharpes-Schoolfield model
  start_vals <- get_start_vals(Sporesden$Temperature[subs], Sporesden$density[subs], model_name = 'sharpeschoolhigh_1981')
  low_lims <- get_lower_lims(Sporesden$Temperature[subs], Sporesden$density[subs], model_name = 'sharpeschoolhigh_1981')
  upper_lims <- get_upper_lims(Sporesden$Temperature[subs], Sporesden$density[subs], model_name = 'sharpeschoolhigh_1981')
    upper_lims[4]<-35

  # Run nls fitting for equation of choice - this one fits the Sharpe-Schoolfield equation
  fit0 <- nls_multstart(density~sharpeschoolhigh_1981(temp = Temperature, r_tref,e,eh,th, tref = 15),
                       data = Sporesden[subs,],
                       iter = 500,
                       start_lower = start_vals - 10,
                       start_upper = start_vals + 10,
                       lower = low_lims,
                       upper = upper_lims,
                       supp_errors = 'Y')
  new_data <- data.frame(Temperature = seq(10, 40, 0.25))

  # Record output parameter estimate values into object
  params<-rbind(params,calc_params(fit0))
  simulated.values<-simulate_nls(fit0,100,3,data=Sporesden[subs,])
  y_predict=rbind(y_predict,predict(fit0, new_data))
  ##array to save simulated results
  results<-NULL
  predictions<-NULL
  for (j in (1:100)) {
    e<-Sporesden[subs,]
    e$density<-simulated.values[,j]
    
    start_vals <- get_start_vals(e$Temperature, e$density, model_name = 'sharpeschoolhigh_1981')
    low_lims <- get_lower_lims(e$Temperature, e$density, model_name = 'sharpeschoolhigh_1981')
    upper_lims <- get_upper_lims(e$Temperature, e$density, model_name = 'sharpeschoolhigh_1981')
    upper_lims[4]<-35
    fit <- nls_multstart(density~sharpeschoolhigh_1981(temp = Temperature, r_tref,e,eh,th, tref = 15),
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
    
  pre.95<-apply(predictions,2,function(x) quantile(x,c(0.025,0.975),na.rm=T))
  polygon(x = c(new_data[,1], rev(new_data[,1])),
          y = c(pre.95[1,],rev(pre.95[2,])),
          col =  adjustcolor(cols[i], alpha.f = 0.10), border = NA)
  #lines(new_data[,1] ,pre.95[1,],col=cols[i],lty=2)
  #lines(new_data[,1] ,pre.95[2,],col=cols[i],lty=2)
  CI.95<-apply(results,2,function(x) quantile(x,c(0.025,0.975),na.rm=T))
  summ.CI<-rbind(summ.CI,apply(CI.95,2,function(x) paste(round(x,2),collapse="-")))
}
##Plot lines of model predictions
lines(new_data[,1] ,y_predict[1,],col='#BB5566')
lines(new_data[,1] ,y_predict[2,],col='#004488')
##Plot means and error bars
subs_ara<-(LMMSpore$Host=="Arabica")
subs_rob<-(LMMSpore$Host=="Robusta")
points(as.numeric(as.character(LMMSpore$Temperature[subs_ara])),LMMSpore$emmean[subs_ara],col='#BB5566',pch=16)
arrows(as.numeric(as.character(LMMSpore$Temperature[subs_ara])),LMMSpore$emmean[subs_ara]+LMMSpore$SE[subs_ara], as.numeric(as.character(LMMSpore$Temperature[subs_ara])),LMMSpore$emmean[subs_ara]-LMMSpore$SE[subs_ara], code=3, angle=90, length=0.1,col='#BB5566')
points(as.numeric(as.character(LMMSpore$Temperature[subs_rob])),LMMSpore$emmean[subs_rob],col='#004488',pch=16)
arrows(as.numeric(as.character(LMMSpore$Temperature[subs_rob])),LMMSpore$emmean[subs_rob]+LMMSpore$SE[subs_rob], as.numeric(as.character(LMMSpore$Temperature[subs_rob])),LMMSpore$emmean[subs_rob]-LMMSpore$SE[subs_rob], code=3, angle=90, length=0.1,col='#004488')
##95% CI
summ.CI.spore<-summ.CI
params.spore<-params

###########################################
##3.Spore germination
#####
##3.1.Data processing
# Read data for germination rate
Germrate <- read_csv("Germrate.csv")
# Specify temperature and strain as factor for LMM analysis
Germrate.stat<-Germrate
Germrate.stat$Temperature=as.factor(Germrate$Temperature)
Germrate.stat$Strain=as.factor(Germrate$Strain)

#####
##3.2.Run LMM for germination rate
# Perform LMM for unit spore conc vs host and temperature,
# round and strain included as random effects
lmem_germ=lmer(Rate ~ Host * Temperature + (1|Round)+(1|Strain), data = Germrate.stat) 
summary(lmem_germ)
emmeans(lmem_germ,list(pairwise~Host*Temperature),adjust="tukey")
rsq(lmem_germ,adj=TRUE)
anova(lmem_germ)
# Produce data frame of LMM results
LMMgerm<-as.data.frame(emmeans(lmem_germ,list(pairwise~Host*Temperature),adjust="tukey",type="response")[[1]])

#####
##3.3.Run NLS for unit spore concentration
# Plot the background for NLS curves
plot(c(10,40),range(Germrate$Rate),type="n",xlab="Temperature (°C)",
     ylab="Proportion germinated spores",cex.lab=1.2,las=1)
# Add legend of different population
legend("bottom", legend=hosts, col=cols,bty="n",lty=1,cex=1.2)

# Clear the dataframe for parameter outputs 
params<-NULL
summ.CI<- NULL
y_predict<- NULL
# Loop through populations
for (i in (1:length(hosts))) {
   # Select population
  host<-hosts[i]
  # Extract data for just this population
  subs<-(Germrate$Host==host)

  # Extract sensible starting guesses for the parameters using Sharpes-Schoolfield model
  start_vals <- get_start_vals(Germrate$Temperature[subs], Germrate$Rate[subs], model_name = 'sharpeschoolhigh_1981')
  low_lims <- get_lower_lims(Germrate$Temperature[subs], Germrate$Rate[subs], model_name = 'sharpeschoolhigh_1981')
  upper_lims <- get_upper_lims(Germrate$Temperature[subs], Germrate$Rate[subs], model_name = 'sharpeschoolhigh_1981')
  
  ##restrict upper limit for th to within data range
  upper_lims[4]<-35
  
  # Run nls fitting for equation of choice - this one fits the Sharpe-Schoolfield equation
  fit2 <- nls_multstart(Rate~sharpeschoolhigh_1981(temp = Temperature, r_tref,e,eh,th, tref = 15),
                       data = Germrate[subs,],
                       iter = 500,
                       start_lower = start_vals - 10,
                       start_upper = start_vals + 10,
                       lower = low_lims,
                       upper = upper_lims,
                       supp_errors = 'Y')
  new_data <- data.frame(Temperature = seq(10, 40, 0.25))
  # Plot lines of predicted values into the empty plot
  #lines(new_data[,1] ,predict(fit1, new_data),col=cols[i])
  y_predict=rbind(y_predict,predict(fit2, new_data))
  # Record output parameter estimate values into object
  params<-rbind(params,calc_params(fit2))
  simulated.values<-simulate_nls(fit2,100,3,data=Germrate[subs,])

##array to save simulated results
  results<-NULL
  predictions<-NULL
  for (j in (1:100)) {
  	e<-Germrate[subs,]
  	e$Rate<-simulated.values[,j]
  	
  	  # Extract sensible starting guesses for the parameters using Sharpes-Schoolfield model
    start_vals <- get_start_vals(e$Temperature, e$Rate, model_name = 'sharpeschoolhigh_1981')
    low_lims <- get_lower_lims(e$Temperature, e$Rate, model_name = 'sharpeschoolhigh_1981')
    upper_lims <- get_upper_lims(e$Temperature, e$Rate, model_name = 'sharpeschoolhigh_1981')
    upper_lims[4]<-35
  
   fit <- nls_multstart(Rate~sharpeschoolhigh_1981(temp = Temperature, r_tref,e,eh,th, tref = 15),
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
  pre.95<-apply(predictions,2,function(x) quantile(x,c(0.025,0.975),na.rm=T))
  polygon(x = c(new_data[,1], rev(new_data[,1])),
          y = c(pre.95[1,],rev(pre.95[2,])),
          col =  adjustcolor(cols[i], alpha.f = 0.10), border = NA)
  #lines(new_data[,1] ,pre.95[1,],col=cols[i],lty=2)
  #lines(new_data[,1] ,pre.95[2,],col=cols[i],lty=2)
  CI.95<-apply(results,2,function(x) quantile(x,c(0.025,0.975),na.rm=T))
  summ.CI<-rbind(summ.CI,apply(CI.95,2,function(x) paste(round(x,2),collapse="-")))
}
##Plot lines of model predictions
lines(new_data[,1] ,y_predict[1,],col='#BB5566')
lines(new_data[,1] ,y_predict[2,],col='#004488')
##Plot means and error bars
subs_ara<-(LMMgerm$Host=="Arabica")
subs_rob<-(LMMgerm$Host=="Robusta")
points(as.numeric(as.character(LMMgerm$Temperature[subs_ara])),LMMgerm$emmean[subs_ara],col='#BB5566',pch=16)
arrows(as.numeric(as.character(LMMgerm$Temperature[subs_ara])),LMMgerm$emmean[subs_ara]+LMMgerm$SE[subs_ara], as.numeric(as.character(LMMgerm$Temperature[subs_ara])),LMMgerm$emmean[subs_ara]-LMMgerm$SE[subs_ara], code=3, angle=90, length=0.1,col='#BB5566')
points(as.numeric(as.character(LMMgerm$Temperature[subs_rob])),LMMgerm$emmean[subs_rob],col='#004488',pch=16)
arrows(as.numeric(as.character(LMMgerm$Temperature[subs_rob])),LMMgerm$emmean[subs_rob]+LMMgerm$SE[subs_rob], as.numeric(as.character(LMMgerm$Temperature[subs_rob])),LMMgerm$emmean[subs_rob]-LMMgerm$SE[subs_rob], code=3, angle=90, length=0.1,col='#004488')
##95% CI
summ.CI.Germ<-summ.CI
params.Germ<-params

##Pull together parameters and CI
param.table<-data.frame(rbind(params.growth,params.spore,params.germ))
param.table<-cbind(Host=rep(hosts,3),Trait=rep(c("Growth","Sporul","Germin"),each=2),param.table)
write.csv(param.table,file="param.table.csv")

CL.table<-data.frame(rbind(summ.CI.growth,summ.CI.spore,summ.CI.germ))
CL.table<-cbind(Host=rep(hosts,3),Trait=rep(c("Growth","Sporul","Germin"),each=2),CL.table)
write.csv(CL.table,file="CL.table.csv")

