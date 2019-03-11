# Description: Mini-project 7: modelling water availability taking slope into account
# Author: Group 19
# Date: 05-03-2019

library(raster)
library(RColorBrewer)

source("scripts/empirical.R")

# Load data (Created in preprocessing.R)
ET_stack <- stack("data/ET_idw_re")
Sret_r <- raster("data/Sret_re") # Raster that assigns S values based on "curve number" approach
precip_stack <- stack("data/precip_stack") # precip data for 8 day period, interpolated spatially from 3 points
FC_r <- raster("data/AdaptedFieldCapacity1.tif") 
FC_r_re <- resample(FC_r, precip_stack) # Field capacity raster based on soil type (in cm)
soil_depth <- raster("data/BDRICM_M_250m.tif")
smabs_stack <- stack("data/smabs_stack_re") # Remotely sensed soil moisture multiplied by soil depth
slope_re <- raster("data/slopes_re")


#########
# Model #
#########
# Initialize soil water availability at .5 of FC
S <- FC_r_re * 10 # FC in cm to FC in mm

# Instead, initialize siol water availability using data
S <- subset(smabs_stack, 1)

# for loop with length = number of timesteps in the rasters 
S_c <- c()   # Vector to store soil moisture snapshots in
deltaS_c <- c() # Vector to store soil moisture differences in
ET_stack10 <- ET_stack * 10 # Reverse incorrect scaling in preprocessing step
for (i in 1:nlayers(ET_stack)){
  # Select the right precip, ET and precip layer
  precip_i <- subset(precip_stack, i)
  ET_i <- subset(ET_stack10, i)
  
  # Calculate run off from Precip raster + Sret raster
  runoff_r <- overlay(Sret_r, precip_i, fun = calc_runoff2)  
  
  # Calculate Leaching (L) from precipitation (P)
  leaching_r <- calc(precip_i, fun = calc_leaching)
  
  # Calculate dS = P + run_on - ET - run_off - L 
  deltaS <- precip_i - ET_i - runoff_r - leaching_r
  
  # Add dS to the raster containing current S and store intermediate raster with datetime in output folder
  S <- S + deltaS
  S_c <- c(S_c, S)
  deltaS_c <- c(deltaS_c, deltaS)
}
Sm_stack <- stack(S_c)
Sm_stack <- trim(Sm_stack)
Sm_stack <- resample(Sm_stack, smabs_stack)

# writeRaster(Sm_stack, "data/sm_stack_calculated", format = "raster")

###########################################
## Visualization and statistical analysis #
###########################################
Sm_stack <- stack("data/sm_stack_calculated")

JDs <- seq(1,366,8) # Julian days
datelst <- as.character(as.Date(JDs, origin=as.Date("2016-01-01")))

# calculating mean directly out of the rasterstack and add it to a list
mean_sm <- c()
for (i in 1:nlayers(Sm_stack)){
  val <- values(Sm_stack[[i]])
  mean <- mean(val, na.rm = T)
  mean_sm[[i]] <- mean
}

# Calcualte area wide mean for remote sensing data
mean_sm_data <- c() 
for (i in 1:nlayers(smabs_stack)){
  val <- values(smabs_stack[[i]])
  mean <- mean(val, na.rm = T)
  mean_sm_data[[i]] <- mean
}
names(mean_sm) <- datelst
names(mean_sm_data) <- datelst

# Get seasonal RMSE
precip_spring <- subset(Sm_stack, 1:11)
precip_summer <- subset(Sm_stack, 12:23)
precip_autumn <- subset(Sm_stack, 24:35)
precip_winter <- subset(Sm_stack, 36:46)

mean_spring <- trim(mean(precip_spring))
mean_summer <- trim(mean(precip_summer))
mean_autumn <- trim(mean(precip_autumn))
mean_winter <- trim(mean(precip_winter))

precip_springD <- subset(smabs_stack, 1:11)
precip_summerD <- subset(smabs_stack, 12:23)
precip_autumnD <- subset(smabs_stack, 24:35)
precip_winterD <- subset(smabs_stack, 36:46)

mean_springD <- trim(mean(precip_springD))
mean_summerD <- trim(mean(precip_summerD))
mean_autumnD <- trim(mean(precip_autumnD))
mean_winterD <- trim(mean(precip_winterD))

#comparing the created dataset with the soilmoisture dataset RMSE
rmsetrySum<-sqrt(mean(values(mean_summer - mean_summerD)^2,na.rm=TRUE))
maxVsu <- max(maxValue(mean_summer), maxValue(mean_summerD))
rmseNSum <- rmsetrySum / maxVsu

rmsetrySpr<-sqrt(mean(values(mean_spring - mean_springD)^2,na.rm=TRUE))
maxVs <- max(maxValue(mean_spring), maxValue(mean_springD))
rmseNspring <- rmsetrySpr / maxVs

rmsetryAut<-sqrt(mean(values(mean_autumn - mean_autumnD)^2,na.rm=TRUE))
maxVa <- max(maxValue(mean_autumn), maxValue(mean_autumnD))
rmseNaut <- rmsetryAut / maxVa

rmsetryWin<-sqrt(mean(values(mean_winter - mean_winterD)^2,na.rm=TRUE))
maxV <- max(maxValue(mean_winter), maxValue(mean_winterD))
rmseNWin <- rmsetryWin / maxV

abs_min <- min(minValue(mean_springD), minValue(mean_autumnD), minValue(mean_winterD), minValue(mean_summerD),
               minValue(mean_spring), minValue(mean_autumn), minValue(mean_winter), minValue(mean_summer))
abs_max <- max(maxVsu, maxVs, maxVa, maxV)

par(mfrow = (c(1,2)))
plot(mean_spring)
cols <- brewer.pal(99, name="YlOrBr")
plot(mean_spring, col=cols, breaks=seq(abs_min, abs_max, 50))

par(mfrow = (c(1,2)))
plot(mean_springD, col=cols, zlim=c(abs_min,abs_max), legend=F, main="Remotely sensed")
plot(mean_spring, col=cols, zlim=c(abs_min, abs_max), legend.args=list(text='soil moisture (mm)', side=4, font=2, line=2.5, cex=0.8), main="Modelled")
mtext(paste0("Mean spring soil moisture model vs. data. RMSE: ", round(rmsetrySpr), " mm"), side = 1, line = -21, outer = TRUE)

par(mfrow = (c(1,2)))
plot(mean_summerD, col=cols, zlim=c(abs_min,abs_max), legend=F, main="Remotely sensed")
plot(mean_summer, col=cols, zlim=c(abs_min, abs_max), legend.args=list(text='soil moisture (mm)', side=4, font=2, line=2.5, cex=0.8), main="Modelled")
mtext(paste0("Mean summer soil moisture model vs. data. RMSE: ", round(rmsetrySum), " mm"), side = 1, line = -21, outer = TRUE)

par(mfrow = (c(1,2)))
plot(mean_autumnD, col=cols, zlim=c(abs_min,abs_max), legend=F, main="Remotely sensed")
plot(mean_autumn, col=cols, zlim=c(abs_min, abs_max), legend.args=list(text='soil moisture (mm)', side=4, font=2, line=2.5, cex=0.8), main="Modelled")
mtext(paste0("Mean autumn soil moisture model vs. data. RMSE: ", round(rmsetryAut), " mm"), side = 1, line = -21, outer = TRUE)

par(mfrow = (c(1,2)))
plot(mean_winterD, col=cols, zlim=c(abs_min,abs_max), legend=F, main="Remotely sensed")
plot(mean_winter, col=cols, zlim=c(abs_min, abs_max), legend.args=list(text='soil moisture (mm)', side=4, font=2, line=2.5, cex=0.8), main="Modelled")
mtext(paste0("Mean winter soil moisture model vs. data. RMSE: ", round(rmsetryWin), " mm"), side = 1, line = -21, outer = TRUE)

# Check correlation with physical parameters

# Correlation with soil depth
soil_depth <- resample(soil_depth, mean_winter)
cor_stackw <- stack(soil_depth, mean_winterD)
cor_stacks <- stack(soil_depth, mean_springD)
cor_stacksu <- stack(soil_depth, mean_summerD)
cor_stackau <- stack(soil_depth, mean_autumnD)

cor()

print("winter")
layerStats(cor_stackw, stat='pearson')
print(mean(values(cor_stackw)))
print("spring")
layerStats(cor_stacks, stat='pearson')
print(mean(values(cor_stacks)))
print("summer")
layerStats(cor_stacksu, stat='pearson')
print(mean(values(cor_stacksu)))
print("autumn")
layerStats(cor_stackau, stat='pearson')
print(mean(values(cor_stackau)))

# Plot correlation with slope
slope_cor <- function(season_mean){
  print(layerStats(stack(season_mean, slope_re), stat="pearson"))
}
plot(values(mean_winter), values(slope_re))
plot(values(mean_winterD), values(slope_re))

plot(values(mean_summer), values(slope_re))
plot(values(mean_summerD), values(slope_re))

par(mfrow = c(1,2))
plot(values(mean_summer), values(soil_depth))
plot(values(mean_winter), values(soil_depth))



## Plotting mean over area
png("images/model_sensed_SM.png")
labs <- names(mean_sm)[seq(1, 48, 12)]
plot(mean_sm, type = "l", col="red", ylim=c(0, 300), ylab = "mean soil moisture (mm)",
     xaxt='n', xlab = "date", main="Modelled vs. remotely sensed soil moisture for study area")
axis(1,at=seq(1, 48, 12),labels=labs)
lines(mean_sm_data, col="blue")
legend(1, 95, legend=c("Modelled", "Remotely sensed"),
       col=c("red", "blue"), lty=1, cex=0.8)
dev.off()