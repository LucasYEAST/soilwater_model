# Description: Mini-project 7: modelling water availability taking slope into account
# Author: Group 19
# Date: 06-03-2019
library(MODIS)
library(gdalUtils)
library(raster)
library(gstat)

##########################
# EVAPOTRANSPIRATION     #
##########################

# I'm only interested in the first subdataset and I can use gdal_translate to convert it to a .tif
dat_files <- list.files("data")
for (fn in dat_files){
  if (grepl('MOD16A2', fn) & grepl('2016', fn)){
    sds <- get_subdatasets(paste0("data/", fn))
    # tiff_fn <- paste0("data/ET2_", substr(fn,10,16))
    tiff_fn <- paste0(fn, ".tif")
    print(tiff_fn)
    gdal_translate(sds[1], dst_dataset = tiff_fn)
  }
}

indices <- grep(".hdf.tif", list.files())
hdftifs <- list.files()[indices]

unique_dates <- c()
for (fn in hdftifs){
  unique_dates <- c(unique_dates, substr(fn,1,16))
}

unique_dates <- unique(unique_dates)
mosaiced_rasters <- list()
empty_raster <- raster("MOD16A2.A2016001.h08v05.006.2017117173504.hdf.tif")
empty_raster[] <- NA
i <- 1
for (uniqdat in unique_dates){
  mosaic_r <- empty_raster
  for (fn in hdftifs){
    if (grepl(uniqdat, fn)){
      r <- raster(fn)
      r[r >= 2999] <- NA
      mosaic_r <- mosaic(mosaic_r, r, fun=mean)
    }
  }
  mosaiced_rasters[[i]] <- mosaic_r
  print(length(mosaiced_rasters))
  print(i)
  i <- i + 1
  
}

# Create mosaic
ET_stack_mosaiced <- stack(mosaiced_rasters)

# Project to WGS84
proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
ET_stack_projected <- projectRaster(ET_stack_mosaiced, crs = proj)
ET_stack_projected <- stack("data/ET_stack_projected")

# Crop to right extent
Sret <- raster("data/Curvenr100.tif")
ET_stack_cropped <- crop(ET_stack_projected, Sret)
ET_stack_scaled <- ET_stack_cropped / 10

# Aggregate to save time
ET_stack_ag <- aggregate(ET_stack_scaled, 5)

# Interpolate voids using Kriging
interpolated_c <- c()
for (i in 1:nlayers(ET_stack_scaled)){ #
  layer <- subset(ET_stack_ag,i)
  # ET_ag <- aggregate(layer, 2)
  ET_ag_df <- as.data.frame(layer, xy = TRUE, na.rm = TRUE)
  mg <- gstat(id = paste0("layer.",i), formula = ET_ag_df[,3]~1, locations = ~x+y,
              data = ET_ag_df, nmax=10, set=list(idp = 3))
  z <- interpolate(layer, mg)
  interpolated_c <- c(interpolated_c, z)
}
ET_stack_idw <- stack(interpolated_c)
writeRaster(ET_stack_idw, "data/ET_stack_idw", format="raster")


#######################
# Precipitation Data  #
#######################

# 3 point interpolated data
precip_files <- list.files("data/precipitation_rasters/", pattern = ".tif$", full.names = T)
# file no. 41 is corrupted, instead repeat file no. 40
precip_files[[41]] <- precip_files[[40]]

precip_stack <- stack(precip_files)
# Extent was flipped around the wrong way
extent(precip_stack) <- extent(extent(precip_stack)[c(3,4,1,2)])

#################
# Soil moisture #
#################
JDs <- seq(0,371,1)
datelst <- as.character(as.Date(JDs, origin=as.Date("2016-01-01")))
datestrs <- gsub("-", "", datelst)
rs <- c()
for (i in seq(1,length(datestrs),8)){
  days_dates8 <- datestrs[c(i:(i+7))]
  stack_files <- c()
  for (date in days_dates8){
    fns <- list.files("data/soil_moisture/", pattern = date, recursive = T, full.names = T)
    stack_files <- c(stack_files, fns)
  }
  stack <- stack(stack_files)
  stack_crop <- crop(stack, precip_stack)
  stack_crop[stack_crop <= -9990] <- NA
  days8_r <- mean(stack_crop, na.rm = T)
  rs <- c(rs, days8_r)
}

# get absolute soil moisture
soil_depth <- raster("data/BDRICM_M_250m.tif")
soil_depth_cropped <- crop(soil_depth, precip_stack)

# Multiply the percentage sm by the soil depth in mm
sm_stack <- stack(rs)
sm_stack_redepth <- resample(sm_stack, soil_depth_cropped) 
smabs_stack <- sm_stack_redepth * soil_depth * 10
smabs_stack_re <- resample(smabs_stack, precip_stack)


# Precipitation had smaller extent, so resample everything to precipitation data
ET_idw_re <- resample(ET_stack_idw, precip_stack)
Sret <- raster("data/Curvenr100.tif")
Sret_re <- resample(Sret, precip_stack)

# Save data
writeRaster(ET_idw_re, "data/ET_idw_re", format="raster")
writeRaster(Sret_re, "data/Sret_re", format="raster")
writeRaster(precip_stack, "data/precip_stack", format="raster")
writeRaster(sm_stack, "data/sm_stack", format = "raster")
writeRaster(smabs_stack_re, "data/smabs_stack_re", format = "raster", overwrite=T)


