## SST by region avergaes
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sp)
library(rgdal)
library(RNetCDF)
library(raster)

GOMstrat <-c(01240, 01260:01300, 01360:01400,03570:03900,01410:01600,01310:01351)
SNEstrat<-c(01010:01120, 03010:03140,03450,03470:03510,03910,03990,03540,03530)
MABstrat <- c(01610:01760,03150:03440)
GBstrat <-c(01130:01230, 01250,03460,03520,03550,03560,02160,02170,02210,01990)
allstrat <- c(GOMstrat, SNEstrat, GBstrat, MABstrat)

prj_crs <- CRS("+proj=aea +lat_1=27.33333333333333 +lat_2=40.66666666666666 +lat_0=34 +lon_0=-78 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

prj_wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

setwd("~/Dropbox/NMFS data/shapefiles/")
strat <- st_read(".", "BTS_strata")
st_crs(strat) <- st_crs(prj_wgs)

allstrat <- c(GOMstrat, SNEstrat, MABstrat, GBstrat)
strat_all_sf <- strat[strat$STRATA %in% allstrat,]
strat_all_sp <- as(strat_all_sf, "Spatial")
st_crs(strat_all_sf) <- st_crs(prj_wgs)

strat_all_sf$region <- ifelse(strat_all_sf$STRATA %in% GOMstrat, "GOMstrat", NA)
strat_all_sf$region <- ifelse(strat_all_sf$STRATA %in% SNEstrat, "SNEstrat", strat_all_sf$region)
strat_all_sf$region <- ifelse(strat_all_sf$STRATA %in% GBstrat, "GBstrat", strat_all_sf$region)
strat_all_sf$region <- ifelse(strat_all_sf$STRATA %in% MABstrat, "MABstrat", strat_all_sf$region)

sf::sf_use_s2(FALSE)
sf_entire_region <- strat_all_sf %>% group_by() %>% summarise()

setwd("~/Dropbox/My Mac (MacBook-Pro)/Documents/PhD/SST_download/REYNOLDS_NCDC_L4_MONTHLY_V5_1963_2020")
list.sst <- list.files(".", pattern = ".tif$", recursive = TRUE)
rast.sst <- lapply(list.sst, raster)
rast.sst <- rast.sst[rep(seq(1,12), 58) >= 9]
rast.sst <- rast.sst[57:232]

rast.sst <- stack(rast.sst)

rast.sst.yr <- stackApply(rast.sst, rep(1:44, each=4), fun = mean)
names(rast.sst.yr) <- paste("sst_yr_", 1977:2020, sep="")

rast.avg_entire <- raster::extract(rast.sst.yr, sf_entire_region, fun=mean)
rast.avg_entire <- data.frame(t(rast.avg_entire))
rast.avg_entire$YEAR <- 1977:2020

rast.max_entire <- raster::extract(rast.sst.yr, sf_entire_region, fun=max)
rast.max_entire <- data.frame(t(rast.max_entire))
rast.max_entire$YEAR <- 1963:2020

### Pull in data from bio_weighted_all.R

#saveRDS(rast.avg_entire, "rastavg_entire.rds")
#saveRDS(rast.sst.yr, "yearlySST_REYNOLDS.rds")

rast.avg_entire1 <- readRDS("rastavg_entire.rds")
rast.avg_entire <- read.csv("rastavg_entire_1977_2020_9_12.csv")


#To find the correct url, go to the opendap link that goes to the folder structure. Navigate to the first day you want, click on "html" and it'll bring you to the Opendap form where you can constrain by lat/lon and figure out what dimensions you need to pull
dim_yr <- rep(1981:2020, each=365)
dim_yday <- rep(1:365, 40)
dim_unique <- paste(dim_yr, dim_yday, sep="_")

dimm <- data.frame(cbind(dim_yr, dim_yday, dim_unique))
dimm <- dimm[244:nrow(dimm),]

#sst_stack_newer <- stack() DONE THROUGH 12/31/1986
for(i in unique(dimm$dim_unique)){
  r <- dimm[dimm$dim_unique == i,]
  yr <- r$dim_yr
  dayy <- (sprintf("%03s", r$dim_yday))
  fakedate <- ymd(paste(yr, "0101", sep="")) + as.numeric(dayy) - 1
  ymdate <- paste(yr, sprintf("%02d",month(fakedate)), sprintf("%02d", day(fakedate)), sep="")
  t <- open.nc(paste("https://podaac-opendap.jpl.nasa.gov/opendap/allData/ghrsst/data/GDS2/L4/GLOB/NCEI/AVHRR_OI/v2/", yr, "/", dayy, "/", ymdate, "120000-NCEI-L4_GHRSST-SSTblend-AVHRR_OI-GLOB-v02.0-fv02.0.nc?analysed_sst[0:1:0][440:1:560][380:1:500],lat[440:1:560],lon[380:1:500]", sep=""))
  t.t <- var.get.nc(t, "analysed_sst")
  t.lat <- var.get.nc(t, "lat")
  t.lon <- var.get.nc(t, "lon")
  
  t.rast <- (flip(raster(t(t.t), xmn = min(t.lon), xmx = max(t.lon), ymn = min(t.lat), ymx = max(t.lat))))
  t.rast <- t.rast - 273.15
  proj4string(t.rast) <- CRS("+init=epsg:4326")
  t.rast <- stack(t.rast)
  
  names(t.rast) <- paste("sst", ymdate, sep="")
  sst_stack_newer <- stack(sst_stack_newer, t.rast)
  print(i)
}


### Older data
dim_yr1 <- rep(1963:2020, each=12)
dim_ymon1 <- rep(1:12, 58)
dim_unique1 <- paste(dim_yr1, dim_ymon1, sep="_")

dimm1 <- data.frame(cbind(dim_yr1, dim_ymon1, dim_unique1))
#dimm1 <- dimm1[244:nrow(dimm1),]

sst_stack_older <- stack()
for(i in unique(dimm1$dim_unique1)){
  r <- dimm1[dimm1$dim_unique1 == i,]
  yr <- r$dim_yr1
  monn <- (sprintf("%02s", r$dim_ymon1))
  #fakedate <- ymd(paste(yr, "0101", sep="")) + as.numeric(dayy) - 1
  ymdate <- paste(yr, monn, sep="")
  t <- open.nc(paste("https://podaac-opendap.jpl.nasa.gov/opendap/allData/ersst/L4/ncei/v5/monthly/netcdf/", yr, "/", monn, "/ersst.v5.", ymdate, ".nc?lat[54:1:69],lon[137:1:153],sst[0:1:0][0:1:0][54:1:69][137:1:153]", sep=""))
  t.t <- var.get.nc(t, "sst")
  t.lat <- var.get.nc(t, "lat")
  t.lon <- var.get.nc(t, "lon")
  t.lon <- t.lon - 360
  
  t.rast <- (flip(raster(t(t.t), xmn = min(t.lon), xmx = max(t.lon), ymn = min(t.lat), ymx = max(t.lat))))
  #t.rast <- t.rast - 273.15
  proj4string(t.rast) <- CRS("+init=epsg:4326")
  t.rast <- stack(t.rast)
  
  names(t.rast) <- paste("sst", yr, "_", monn, sep="")
  sst_stack_older <- stack(sst_stack_older, t.rast)
  print(i)
}

#setwd("~/Dropbox/My Mac (MacBook-Pro)/Documents/PhD/SST_download/REYNOLDS_NCDC_L4_MONTHLY_V5_1963_2020")
#writeRaster(sst_stack_older, paste(names(sst_stack_older), ".tif", sep=""), bylayer = TRUE)

