#Prep biomass weighted metrics
#poleward shifts 
# This version is compiled from bio_analysis_all_20220216.R and Bio_analysis_by_group_20220208.R

library(spatialEco)
library(sp)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)
library(rnaturalearthdata)
library(raster)
library(tidyr)
library(dplyr)
library(beepr)
library(ggplot2)
library(reshape2)
library(lubridate)
library(viridis)
library(scales)

#Define the projections & download land. Projection based on custom projection in Roberts et al., 2016. 
prj_crs <- CRS("+proj=aea +lat_1=27.33333333333333 +lat_2=40.66666666666666 +lat_0=34 +lon_0=-78 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
prj_wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
land <- ne_countries(scale=10, country=c("United States of America", "Canada"), returnclass = "sf")

############### BOTTOM TRAWL ANALYSIS AND SETUP ---------------
#Download most recent BT data from the Nye Lab Dropbox. For access email Janet Nye. 
setwd("~/Dropbox/NMFS data/")
load("NEFSC_BTS_2021_all_seasons.RData")

#Define the strata used based on continuously sampled strata and define subregions by code from Caro.
GOMstrat <-c(01240, 01260:01300, 01360:01400,03570:03900,01410:01600,01310:01351)
SNEstrat<-c(01010:01120, 03010:03140,03450,03470:03510,03910,03990,03540,03530)
MABstrat <- c(01610:01760,03150:03440)
GBstrat <-c(01130:01230, 01250,03460,03520,03550,03560,02160,02170,02210,01990)
allstrat <- c(GOMstrat, SNEstrat, GBstrat, MABstrat)

#Select survdat dataset
survdat <- survey[[1]]
#Create table of combos of the SVSPP numbers and the species common names
SVSPP_table <- unique(data.frame(survdat$SVSPP, survdat$COMNAME))
SVSPP_table2 <- unique(data.frame(survdat$SVSPP, survdat$COMNAME))

#The BT data are structured where there are multiple rows per tow per species. They are based on length classes and more. Important to consider this depending on your desired data goal. (e.g. presence-absence requires creation of pseudoabsences in tows where your modeled species was not recorded).
#Select the columns that are relevant for my analyses and DO NOT have row-specific variation (i.e. length). 
survdat_base <- survdat %>% 
  dplyr::select(-c("CATCHSEX",  "LENGTH", "NUMLEN"))
#Calc unique rows from the subset. This gets rid of repeated rows per spp per trawl
survdat_base <- unique(survdat_base)

#Select fall and filter for all the strata we want to assess
survdat_fall <- survdat_base[survdat_base$SEASON == "FALL",]
survdat_fall <- dplyr::filter(survdat_fall, STRATUM %in% allstrat)

#Bin in 5 year categories - may or may not use this in final
survdat_fall$yr5 <- ifelse(survdat_fall$YEAR >= 1963 & survdat_fall$YEAR <= 1964, "1963-1965", NA)
survdat_fall$yr5 <- ifelse(survdat_fall$YEAR >= 1965 & survdat_fall$YEAR <= 1969, "1966-1970", survdat_fall$yr5)
survdat_fall$yr5 <- ifelse(survdat_fall$YEAR >= 1970 & survdat_fall$YEAR <= 1974, "1971-1975", survdat_fall$yr5)
survdat_fall$yr5 <- ifelse(survdat_fall$YEAR >= 1975 & survdat_fall$YEAR <= 1979, "1976-1980", survdat_fall$yr5)
survdat_fall$yr5 <- ifelse(survdat_fall$YEAR >= 1980 & survdat_fall$YEAR <= 1984, "1981-1985", survdat_fall$yr5)
survdat_fall$yr5 <- ifelse(survdat_fall$YEAR >= 1985 & survdat_fall$YEAR <= 1989, "1986-1990", survdat_fall$yr5)
survdat_fall$yr5 <- ifelse(survdat_fall$YEAR >= 1990 & survdat_fall$YEAR <= 1994, "1991-1995", survdat_fall$yr5)
survdat_fall$yr5 <- ifelse(survdat_fall$YEAR >= 1995 & survdat_fall$YEAR <= 1999, "1996-2000", survdat_fall$yr5)
survdat_fall$yr5 <- ifelse(survdat_fall$YEAR >= 2000 & survdat_fall$YEAR <= 2004, "2001-2005", survdat_fall$yr5)
survdat_fall$yr5 <- ifelse(survdat_fall$YEAR >= 2005 & survdat_fall$YEAR <= 2009, "2006-2010", survdat_fall$yr5)
survdat_fall$yr5 <- ifelse(survdat_fall$YEAR >= 2010 & survdat_fall$YEAR <= 2014, "2011-2015", survdat_fall$yr5)
survdat_fall$yr5 <- ifelse(survdat_fall$YEAR >= 2015 & survdat_fall$YEAR <= 2019, "2016-2019", survdat_fall$yr5)

survdat_fall$yr10 <- ifelse(survdat_fall$YEAR >= 1980 & survdat_fall$YEAR <= 1989, "1980-1989", NA)
survdat_fall$yr10 <- ifelse(survdat_fall$YEAR >= 1990 & survdat_fall$YEAR <= 1999, "1990-1999", survdat_fall$yr10)
survdat_fall$yr10 <- ifelse(survdat_fall$YEAR >= 2000 & survdat_fall$YEAR <= 2009, "2000-2009", survdat_fall$yr10)
survdat_fall$yr10 <- ifelse(survdat_fall$YEAR >= 2010 & survdat_fall$YEAR <= 2019, "2010-2019", survdat_fall$yr10)

#Define the strata subregions as a column
survdat_fall$stratcat <- ifelse(survdat_fall$STRATUM %in% GOMstrat, "GOMstrat", NA)
survdat_fall$stratcat <- ifelse(survdat_fall$STRATUM %in% SNEstrat, "SNEstrat", survdat_fall$stratcat)
survdat_fall$stratcat <- ifelse(survdat_fall$STRATUM %in% GBstrat, "GBstrat", survdat_fall$stratcat)
survdat_fall$stratcat <- ifelse(survdat_fall$STRATUM %in% MABstrat, "MABstrat", survdat_fall$stratcat)

#Now we begin the stratified mean latitude calculations
#From the Nye database, load the stratum area shp. This has taken a number of forms over the years depending on who has changed filenames in the group folder
load("stratum.area.RData")
setwd("~/Dropbox/NMFS data/")
load("stratum.area.Rdata")
setwd("~/Dropbox/NMFS data/shapefiles/")
strat <- st_read(".", "BTS_strata")

#Merge stratum sf with the survdat data prepared above
survdat_fall <- merge(survdat_fall, stratum.area, by.x="STRATUM", by.y="STRATUM", all.x=T)

#Define location of Cape Hatteras (based on Nye et al 2009)
ch_loc <- data.frame(t(c("ch", -75,35)))
colnames(ch_loc) <- c("loc", "lon", "lat")

#Turn into sf object, define WGS projection, then reproject to Roberts et al 2016 custom projection
ch_loc <- ch_loc %>%
  st_as_sf(coords = c("lon", "lat") ,crs=4326) %>%
  st_transform(crs= crs(prj_crs))

#Newer version - based on Lesley's analysis. Created a raster in ArcPro using a land and offshore mask to calculate an along-shelf euclidean distance. Mask was custom built by L. Thorne. Euclidean distance was calculated from a horizontal line at 35N. Email J. Stepanuk for access.
setwd("~/Documents/Github/cet_prey/calcs_envr/")
euc_dist_along_shelf <- raster("euc_dist_along_shelf.tif")

#Turn the BT data to an sf object and project
survdat_fall_sf <- survdat_fall %>%
  st_as_sf(coords=c("LON", "LAT"), crs = 4326, remove = FALSE) %>%
  st_transform(crs= crs(prj_crs))

#For initial method, calculate straight line euclidean distance with the st_distance function
survdat_fall_sf$dCH <- t(st_distance(ch_loc, survdat_fall_sf))
#For second method, calculate along shelf euclidean distance by extracting those values from the ArcPro raster
survdat_fall_sf$dCH_AS <- raster::extract(euc_dist_along_shelf, survdat_fall_sf)

#Calculate the YEARLY biomass weighted mean along-shelf distance from the Cape Hatteras region 
#This value is housed in "wt_mean"
survdat_fall_min <- st_drop_geometry(cbind(survdat_fall_sf, st_coordinates(survdat_fall_sf)))

### This calculates the biomass weighted by GROUP

survdat_fall_min$group <- ifelse(survdat_fall_min$SVSPP %in% c(32, 851, 205, 859, 427, 430, 30, 31, 44),  "Clupeid", NA)
survdat_fall_min$group <- ifelse(survdat_fall_min$SVSPP %in% c(121, 578, 572, 209, 124, 898, 570, 702, 571, 582, 744, 860, 208, 212, 211, 874, 822, 745, 577, 877), "Scombrid", survdat_fall_min$group)
survdat_fall_min$group <- ifelse(survdat_fall_min$SVSPP %in% c(181, 734, 38), "Sandlance", survdat_fall_min$group)
survdat_fall_min$group <- ifelse((survdat_fall_min$SVSPP >= 501 & survdat_fall_min$SVSPP <= 505) | (survdat_fall_min$SVSPP >= 510 & survdat_fall_min$SVSPP <= 513), "Squid_BT", survdat_fall_min$group)
survdat_fall_min$group <- ifelse(survdat_fall_min$SVSPP == 297 | survdat_fall_min$SVSPP == 298 | survdat_fall_min$SVSPP == 306 | survdat_fall_min$SVSPP == 307, "Shrimp_BT", survdat_fall_min$group)
survdat_fall_min$group <- ifelse(survdat_fall_min$SVSPP == 306, "Northern_shrimp", survdat_fall_min$group)
survdat_fall_min$group <- ifelse(survdat_fall_min$SVSPP %in% c(73:75), "Large_Gadid", survdat_fall_min$group)
survdat_fall_min$group <- ifelse(survdat_fall_min$SVSPP %in% c(189, 1, 262, 192, 131, 750, 176, 168, 249, 193, 155, 164), "Misc_Fish", survdat_fall_min$group) 
survdat_fall_min$group <- ifelse(survdat_fall_min$SVSPP %in% c(81, 454, 80, 86, 79, 69, 77, 72, 455, 78, 76), "Small_Gadid", survdat_fall_min$group)
# survdat_fall_min$group <- ifelse(survdat_fall_min$SVSPP %in% c(790, 792, 110, 793, 777, 100, 104, 785, 786, 788, 109, 795, 775, 776, 773, 791, 787, 117, 789, 783, 103, 853, 774, 873, 106, 107, 105), "Flatfish", survdat_fall_min$group)
survdat_fall_min$group <- ifelse(survdat_fall_min$SVSPP %in% c(51, 54, 56, 55, 229) , "Meso_pelagic", survdat_fall_min$group)

#survdat_fall_min$group <- ifelse(survdat_fall_min$group == "Shrimp_BT" | survdat_fall_min$group == "Northern_shrimp", "Shrimp_BT", survdat_fall_min$group)


#Calculate the YEARLY biomass weighted mean along-shelf distance from the Cape Hatteras region 
#This value is housed in "wt_mean"
group_str_bio <-survdat_fall_min %>%
  group_by(group, YEAR, SEASON) %>%
  summarise(wt_mean=weighted.mean(dCH_AS,BIOMASS, na.rm=TRUE))

group_str_bio$wt_mean <- as.numeric(group_str_bio$wt_mean)
group_str_comb <- group_str_bio #This is so we can move to combining ECOMON data later

## Now do the above by region
group_str_bio_reg <-survdat_fall_min %>%
  group_by(group, YEAR, SEASON, stratcat) %>%
  summarise(wt_mean=weighted.mean(dCH_AS,BIOMASS, na.rm=TRUE))

group_str_bio_reg$wt_mean <- as.numeric(group_str_bio_reg$wt_mean)
group_str_comb_reg <- group_str_bio_reg #This is so we can move to combining ECOMON data later

### EXTRACT SST
#This was calculated in "download_and_average_SST.R"
setwd("~/Dropbox/My Mac (MacBook-Pro)/Documents/PhD/SST_download/REYNOLDS_NCDC_L4_MONTHLY_V5_1963_2020")
rast.avg_entire <- read.csv("rastavg_entire_1977_2020_9_12.csv") #entire site
group_str_comb_sst <- merge(group_str_comb, rast.avg_entire, by.x="YEAR", by.y="YEAR")
group_str_comb_sst <- merge(group_str_comb_sst, bt_all, by.x="YEAR", by.y="YEAR")


setwd("~/Dropbox/My Mac (MacBook-Pro)/Documents/PhD/SST_download/REYNOLDS_NCDC_L4_MONTHLY_V5_1963_2020")
rast.avg_entire_reg <- read.csv("avgSST_reynolds_byregion.csv", header=T) #by region
rast.avg_entire_reg_melt <- melt(rast.avg_entire_reg[,-1], id.vars = "YEAR")
colnames(rast.avg_entire_reg_melt)[3] <- "SURFTEMP"
group_str_comb_reg_sst <- merge(group_str_comb_reg, rast.avg_entire_reg_melt, by.x=c("YEAR", "stratcat"), by.y=c("YEAR", "variable"), all.x=T)

group_str_comb_reg_sst <- merge(group_str_comb_reg_sst, bt_region, by.x=c("YEAR", "stratcat"), by.y=c("YEAR", "stratcat"))




####################### Ecomon
#Downloaded from https://www.bco-dmo.org/dataset/3327
#It is requested that people use the following standard acknowledgement when publishing or presenting results that utilize NMFS data extracted from BCO-DMO: "This study utilized data collected by the NOAA's Northeast Fisheries Science Center as part of an ongoing mission to monitor and assess the Northeast Continental Shelf ecosystem." Source: NOAA NMFS NEFSC Oceanography Branch"
setwd("~/Dropbox/Ecomon_v3_7_donotdistribute/")
ZPD  <- read.csv("EcoMon_Plankton_Data_v3_7_do not distribute.csv")

#SPP codes has been hand altered to group species as zooplankton, shrimp, or neither based on Smith et al. 2015. 
spp_codes <- read.csv("spp_codes.csv", header=T)
spp_codes$vol <- ifelse(grepl("m3", spp_codes$COLUMN.DESCRIPTION), 1, 0)
spp_codes <- spp_codes[spp_codes$vol == 1,]
spp_codes$tax_num <- as.numeric(spp_codes$TAXA.CODES.IN.DATABASE)
spp_codes$taxa_edited <- spp_codes$TAXA.CODES.IN.DATABASE
spp_codes$taxa_edited <- gsub("and", ",", spp_codes$taxa_edited)
spp_codes$taxa_edited <- gsub("-", ":", spp_codes$taxa_edited)
spp_codes$combod <- ifelse(grepl(",", spp_codes$taxa_edited) & !grepl(":", spp_codes$taxa_edited), 1, 0)
grptax_seq <- spp_codes$taxa_edited[spp_codes$combod==1]
spp_codes$combod2 <- ifelse(grepl(":", spp_codes$taxa_edited), 1, 0)
grptax_rng <- spp_codes$taxa_edited[spp_codes$combod2==1]

grptax_seq <- gsub("\\s+", "", grptax_seq, perl=T)
grptax_seq <- as.numeric(unlist(strsplit(grptax_seq, ",")))

expand.dash <- function(dashed) {
  limits <- as.numeric(as.character(unlist(strsplit(dashed, ':'))))
  seq(limits[1], limits[2])
}
grptax_rng <- unlist(strsplit(grptax_rng, ","))
grptax_rng_ext <- grptax_rng[!grepl(":", grptax_rng)]
zzz <- grptax_rng[grepl(":", grptax_rng)]
for(i in 1:length(zzz)){
  qqq <- zzz[i]
  ppp <- expand.dash(qqq)
  grptax_rng_ext <- c(grptax_rng_ext, ppp)
}

grptax_tot <- c(grptax_seq, grptax_rng_ext)

spp_codes$notingroup <- as.numeric(spp_codes$TAXA.CODES.IN.DATABASE) %in% grptax_tot
spp_codes$includeinfinal <- ifelse(spp_codes$notingroup == "FALSE" & !is.na(spp_codes$tax_num), 1, 0)
spp_codes$includeinfinal <- ifelse(spp_codes$combod ==1 | spp_codes$combod2 == 1, 1, spp_codes$includeinfinal)

spp_codes_sub <- spp_codes[spp_codes$includeinfinal == 1,]
#Melt the dataset by most columns to move from species specific columns
ZPD_melt <- melt(ZPD, id=c("cruise_name", "station", "zoo_gear", "ich_gear", "lat", "lon", "date", "time", "depth", "sfc_temp", "sfc_salt", "btm_temp", "btm_salt", "volume_1m2"))

ZPD_melt <- merge(ZPD_melt, spp_codes_sub, by.x="variable", by.y="COLUMN.NAME", all.x=T)
ZPD_melt$value <- as.numeric(as.character(ZPD_melt$value))

#Get MDY
ZPD_melt$date <- dmy(ZPD_melt$date)
ZPD_melt$month <- month(ZPD_melt$date)
ZPD_melt$year <- year(ZPD_melt$date)

#Subset by fall only
ZPD_fall <- ZPD_melt[ZPD_melt$month >= 9 & ZPD_melt$month <= 11,]
ZPD_fall <- ZPD_fall[!is.na(ZPD_fall$group_multi),]

#Read in BT Strata - remember the ECOMON strata are the same as the BT ones
strat_zp <- st_read(".", "EcomonStrata_v4b")
st_crs(strat_zp) <- st_crs(prj_wgs)

setwd("~/Dropbox/NMFS data/shapefiles/")
strat <- st_read(".", "BTS_strata")
st_crs(strat) <- st_crs(prj_wgs)

#Choose strata we want, and assign subregions
allstrat <- c(GOMstrat, SNEstrat, MABstrat, GBstrat)
strat_reg <- strat[strat$STRATA %in% allstrat,]
strat_reg_sp <- as(strat_reg, "Spatial")
st_crs(strat_reg) <- st_crs(strat_reg)
strat_reg$stratcat <- ifelse(strat_reg$STRATA %in% GOMstrat, "GOMstrat", NA)
strat_reg$stratcat <- ifelse(strat_reg$STRATA %in% SNEstrat, "SNEstrat", strat_reg$stratcat)
strat_reg$stratcat <- ifelse(strat_reg$STRATA %in% GBstrat, "GBstrat", strat_reg$stratcat)
strat_reg$stratcat <- ifelse(strat_reg$STRATA %in% MABstrat, "MABstrat", strat_reg$stratcat)

#Convert to sf 
ZPD_fall_sf <- st_as_sf(ZPD_fall, coords=c("lon", "lat"))
st_crs(ZPD_fall_sf) <- st_crs(prj_wgs)

#Something happening with the closure of these polygons. Get around this with sf_use_s2
#strat_reg[2] is the strata by number
sf_use_s2(FALSE)
ZPD_fall_sf <- st_intersection(ZPD_fall_sf, strat_zp)

#Aggregate abundance by zoop/shrimp and year
ZPD_fall_sf_year <- aggregate(value ~ group_multi_fish + year, data=ZPD_fall_sf, FUN=sum)

#Normalize the ZPD data, we don't actually use this again
ZPD_temp <- NULL
for (i in unique(ZPD_fall_sf_year$group_multi)){
  p <- ZPD_fall_sf_year[ZPD_fall_sf_year$group_multi == i,]
  a <- (p$value - min(p$value))/(max(p$value) - min(p$value))
  w <- data.frame(p, a)
  ZPD_temp <- rbind(ZPD_temp, w)
}

#Assign 5 year categories as before with BT data
ZPD_fall_sf$yr5 <- ifelse(ZPD_fall_sf$year >= 1963 & ZPD_fall_sf$year <= 1964, "1963-1965", NA)
ZPD_fall_sf$yr5 <- ifelse(ZPD_fall_sf$year >= 1965 & ZPD_fall_sf$year <= 1970, "1966-1970", ZPD_fall_sf$yr5)
ZPD_fall_sf$yr5 <- ifelse(ZPD_fall_sf$year >= 1970 & ZPD_fall_sf$year <= 1974, "1971-1975", ZPD_fall_sf$yr5)
ZPD_fall_sf$yr5 <- ifelse(ZPD_fall_sf$year >= 1975 & ZPD_fall_sf$year <= 1979, "1976-1980", ZPD_fall_sf$yr5)
ZPD_fall_sf$yr5 <- ifelse(ZPD_fall_sf$year >= 1980 & ZPD_fall_sf$year <= 1984, "1981-1985", ZPD_fall_sf$yr5)
ZPD_fall_sf$yr5 <- ifelse(ZPD_fall_sf$year >= 1985 & ZPD_fall_sf$year <= 1989, "1986-1990", ZPD_fall_sf$yr5)
ZPD_fall_sf$yr5 <- ifelse(ZPD_fall_sf$year >= 1990 & ZPD_fall_sf$year <= 1994, "1991-1995", ZPD_fall_sf$yr5)
ZPD_fall_sf$yr5 <- ifelse(ZPD_fall_sf$year >= 1995 & ZPD_fall_sf$year <= 1999, "1996-2000", ZPD_fall_sf$yr5)
ZPD_fall_sf$yr5 <- ifelse(ZPD_fall_sf$year >= 2000 & ZPD_fall_sf$year <= 2004, "2001-2005", ZPD_fall_sf$yr5)
ZPD_fall_sf$yr5 <- ifelse(ZPD_fall_sf$year >= 2005 & ZPD_fall_sf$year <= 2009, "2006-2010", ZPD_fall_sf$yr5)
ZPD_fall_sf$yr5 <- ifelse(ZPD_fall_sf$year >= 2010 & ZPD_fall_sf$year <= 2014, "2011-2015", ZPD_fall_sf$yr5)
ZPD_fall_sf$yr5 <- ifelse(ZPD_fall_sf$year >= 2015 & ZPD_fall_sf$year <= 2019, "2016-2019", ZPD_fall_sf$yr5)

ZPD_fall_sf$yr10 <- ifelse(ZPD_fall_sf$year >= 1977 & ZPD_fall_sf$year <= 1989, "1977-1989", NA)
ZPD_fall_sf$yr10 <- ifelse(ZPD_fall_sf$year >= 1990 & ZPD_fall_sf$year <= 1999, "1990-1999", ZPD_fall_sf$yr10)
ZPD_fall_sf$yr10 <- ifelse(ZPD_fall_sf$year >= 2000 & ZPD_fall_sf$year <= 2009, "2000-2009", ZPD_fall_sf$yr10)
ZPD_fall_sf$yr10 <- ifelse(ZPD_fall_sf$year >= 2010 & ZPD_fall_sf$year <= 2019, "2010-2019", ZPD_fall_sf$yr10)

#Assign strata categories as before with BT data
GOMstrat_zp <- unique(strat_zp$Name[strat_zp$Region == "GOM"])
SNEstrat_zp <- unique(strat_zp$Name[strat_zp$Region == "SNE"])
GBstrat_zp <- unique(strat_zp$Name[strat_zp$Region == "GB"])
MABstrat_zp <- unique(strat_zp$Name[strat_zp$Region == "MAB"])

ZPD_fall_sf$stratcat <- ifelse(ZPD_fall_sf$Name %in% GOMstrat_zp, "GOMstrat", NA)
ZPD_fall_sf$stratcat <- ifelse(ZPD_fall_sf$Name %in% SNEstrat_zp, "SNEstrat", ZPD_fall_sf$stratcat)
ZPD_fall_sf$stratcat <- ifelse(ZPD_fall_sf$Name %in% GBstrat_zp, "GBstrat", ZPD_fall_sf$stratcat)
ZPD_fall_sf$stratcat <- ifelse(ZPD_fall_sf$Name %in% MABstrat_zp, "MABstrat", ZPD_fall_sf$stratcat)

#Project so everything matches Roberts et al 2016 custom projection
ZPD_fall_sf_prj <- st_transform(ZPD_fall_sf, crs(prj_crs))
#Calculate distance from CH using first method: straight line euclidean distance from 35,-75
ZPD_fall_sf_prj$dCH <- t(st_distance(ch_loc, ZPD_fall_sf_prj))
#Calculate distance from CH using second method: along shelf dist from 35N
ZPD_fall_sf_prj$dCH_AS <- raster::extract(euc_dist_along_shelf, ZPD_fall_sf_prj)

#Because this file is so massive, drop the sf designation for the ZPD data but keep the projected coordinates as a column in the df
zpd_df <- st_drop_geometry(cbind(ZPD_fall_sf_prj, st_coordinates(ZPD_fall_sf_prj)))
zpd_df$group <- ifelse(zpd_df$TAXA.NAME %in% preyspp_zoop, "zoop", NA)
zpd_df$group <- ifelse(zpd_df$TAXA.NAME %in% preyspp_shrp, "shrimp", zpd_df$group)

#zpd_df$group_multi <- ifelse(zpd_df$group_multi == "euphausiid", "shrimp", zpd_df$group_multi)
#Calculate the annual along-shelf weighted mean distance from Cape Hatteras using along-shelf Euclidean Distance
group_zpd_str_bio <-zpd_df %>%
  group_by(group_multi, year) %>%
  summarise(wt_mean= with(na.omit(data.frame(dCH_AS, value)), weighted.mean(dCH_AS, value, na.rm=T)))

group_zpd_str_bio$wt_mean <- as.numeric(group_zpd_str_bio$wt_mean)
group_zpd_str_comb <- data.frame(group_zpd_str_bio) #Do this for matching up later without losing all the work we've done so far

#Again, extract SST annual to the rows
group_zpd_str_merge <- merge(group_zpd_str_comb, rast.avg_entire, by.x="year", by.y="YEAR")
group_zpd_str_merge <- merge(group_zpd_str_merge, bt_all, by.x="year", by.y="YEAR")

#do same as above but with regions
group_zpd_str_bio_REG <-zpd_df %>%
  group_by(group_multi, year, stratcat) %>%
  summarise(wt_mean=weighted.mean(dCH_AS,value, na.rm=TRUE))

group_zpd_str_bio_REG$wt_mean <- as.numeric(group_zpd_str_bio_REG$wt_mean)
group_zpd_str_comb_REG <- data.frame(group_zpd_str_bio_REG)
#extract regional SST
group_zpd_str_merge_region <- merge(group_zpd_str_comb_REG, rast.avg_entire_reg_melt, by.x=c("year", "stratcat"), by.y=c("YEAR", "variable"), all.x=T)

#extract regional BT
group_zpd_str_merge_region <- merge(group_zpd_str_merge_region, bt_region, by.x=c("year", "stratcat"), by.y=c("YEAR", "stratcat"))

###### Combine BT and ECOMON data ----------------
#Select similar columns for BT and Ecomon, assign correct columnames, merge
group_fish_dch <- group_str_comb_sst %>% dplyr::select(YEAR, group, wt_mean, t.rast.avg_entire.)
group_zoop_dch <- group_zpd_str_merge %>% dplyr::select(year, group_multi, wt_mean, t.rast.avg_entire.)
colnames(group_fish_dch)[1:2] <- colnames(group_zoop_dch)[1:2]

group_allprey_dch <- data.frame(rbind(group_fish_dch, group_zoop_dch))

### Do same for regions
group_fish_dch_reg <- group_str_comb_reg_sst %>% dplyr::select(YEAR, group, stratcat, wt_mean, SURFTEMP)
group_zoop_dch_reg <- group_zpd_str_merge_region %>% dplyr::select(year, group_multi, stratcat, wt_mean, SURFTEMP)
colnames(group_fish_dch_reg)[1:2] <- colnames(group_zoop_dch_reg)[1:2]

group_allprey_reg_dch <- data.frame(rbind(group_fish_dch_reg, group_zoop_dch_reg))

##### END RESULTS #####
#survdat_fall - all spp all replicates all years, but for fall only
#survdat_fall_min - all spp one replicate per trawl each, fall only

#ZPD fall df - zoops that are already filtered for relevant groups, fall only (9-11)

#group_allprey_reg_dch - the annual distance from CH including annual surftemp and bottemp for each species AND zoop groups