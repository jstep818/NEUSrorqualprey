#calculate stratified biomass
#This comes directly after S01 - the biomass weighted calcs

####################3 STRATIFIED BIOMASS CALCS -------
#remotes::install_github("noaa-edab/ecodata",build_vignettes=TRUE)
#remotes::install_github("NOAA-EDAB/survdat",build_vignettes = TRUE)
library(ecodata)
library(survdat)
library(data.table)

survdat_fall$group <- ifelse(survdat_fall$SVSPP %in% c(32, 851, 205, 859, 427, 430, 30, 31, 44),  "Clupeid", NA)
survdat_fall$group <- ifelse(survdat_fall$SVSPP %in% c(121, 578, 572, 209, 124, 898, 570, 702, 571, 582, 744, 860, 208, 212, 211, 874, 822, 745, 577, 877), "Scombrid", survdat_fall$group)
survdat_fall$group <- ifelse(survdat_fall$SVSPP %in% c(181, 734, 38), "Sandlance", survdat_fall$group)
survdat_fall$group <- ifelse((survdat_fall$SVSPP >= 501 & survdat_fall$SVSPP <= 505) | (survdat_fall$SVSPP >= 510 & survdat_fall$SVSPP <= 513), "Squid_BT", survdat_fall$group)
survdat_fall$group <- ifelse(survdat_fall$SVSPP == 297 | survdat_fall$SVSPP == 298 | survdat_fall$SVSPP == 306 | survdat_fall$SVSPP == 307, "Shrimp_BT", survdat_fall$group)
survdat_fall$group <- ifelse(survdat_fall$SVSPP == 306, "Northern_shrimp", survdat_fall$group)
survdat_fall$group <- ifelse(survdat_fall$SVSPP %in% c(73:75), "Large_Gadid", survdat_fall$group)
survdat_fall$group <- ifelse(survdat_fall$SVSPP %in% c(189, 1, 262, 192, 131, 750, 176, 168, 249, 193, 155, 164), "Misc_Fish", survdat_fall$group)
survdat_fall$group <- ifelse(survdat_fall$SVSPP %in% c(81, 454, 80, 86, 79, 69, 77, 72, 455, 78, 76), "Small_Gadid", survdat_fall$group)
# survdat_fall$group <- ifelse(survdat_fall$SVSPP %in% c(790, 792, 110, 793, 777, 100, 104, 785, 786, 788, 109, 795, 775, 776, 773, 791, 787, 117, 789, 783, 103, 853, 774, 873, 106, 107, 105), "Flatfish", survdat_fall$group)
survdat_fall$group <- ifelse(survdat_fall$SVSPP %in% c(51, 54, 56, 55, 229) , "Meso_pelagic", survdat_fall$group)
survdat_fall$CATCHSEX <- 0

#Calculate stratified mean for all regions combined, but then also do it by strata regions
t <- calc_stratified_mean(data.table(survdat_fall), areaDescription = "STRATA", groupDescription = "group", filterBySeason = "FALL")
t_gom <- calc_stratified_mean(data.table(survdat_fall[survdat_fall$stratcat == "GOMstrat",-c(25)]), areaDescription = "STRATA", groupDescription = "group", filterBySeason = "FALL")
t_gb <- calc_stratified_mean(data.table(survdat_fall[survdat_fall$stratcat == "GBstrat",-c(25)]), areaDescription = "STRATA", groupDescription = "group", filterBySeason = "FALL")
t_sne <- calc_stratified_mean(data.table(survdat_fall[survdat_fall$stratcat == "SNEstrat",-c(25)]), areaDescription = "STRATA", groupDescription = "group", filterBySeason = "FALL")
t_mab <- calc_stratified_mean(data.table(survdat_fall[survdat_fall$stratcat == "MABstrat",-c(25)]), areaDescription = "STRATA", groupDescription = "group", filterBySeason = "FALL")

#Assign category for each region
t_gom$stratcat <- "GOMstrat"
t_gb$stratcat <- "GBstrat"
t_sne$stratcat <- "SNEstrat"
t_mab$stratcat <- "MABstrat"

#Combine regions 
t_region <- rbind(t_gom, t_gb)
t_region <- rbind(t_region, t_sne)
t_region <- rbind(t_region, t_mab)

#Optional - read in and plot the strata
#areaPolygon <- sf::st_read(dsn = system.file("extdata","EPU.shp",package="survdat"), quiet=T)
#plot_shapefile(areaPolygon)

##### ALL BIOMASS PER REGION
#Calculate stratified mean for all regions combined, but then also do it by strata regions
#t <- calc_stratified_mean(data.table(survdat_fall), areaDescription = "STRATA", groupDescription = "group", filterBySeason = "FALL")
big_survdat_fall <- survdat_fall
big_survdat_fall$isprey <- ifelse(!is.na(big_survdat_fall$group), "prey", "not_prey")

bigt_gom <- calc_stratified_mean(data.table(big_survdat_fall[big_survdat_fall$stratcat == "GOMstrat",-c(25)]), areaDescription = "STRATA", groupDescription = "isprey", filterBySeason = "FALL")
bigt_gb <- calc_stratified_mean(data.table(big_survdat_fall[big_survdat_fall$stratcat == "GBstrat",-c(25)]), areaDescription = "STRATA", groupDescription = "isprey",  filterBySeason = "FALL")
bigt_sne <- calc_stratified_mean(data.table(big_survdat_fall[big_survdat_fall$stratcat == "SNEstrat",-c(25)]), areaDescription = "STRATA", groupDescription = "isprey", filterBySeason = "FALL")
bigt_mab <- calc_stratified_mean(data.table(big_survdat_fall[big_survdat_fall$stratcat == "MABstrat",-c(25)]), areaDescription = "STRATA", groupDescription = "isprey", filterBySeason = "FALL")

#Assign category for each region
t_gom$stratcat <- "GOMstrat"
t_gb$stratcat <- "GBstrat"
t_sne$stratcat <- "SNEstrat"
t_mab$stratcat <- "MABstrat"

#Combine regions 
t_region <- rbind(t_gom, t_gb)
t_region <- rbind(t_region, t_sne)
t_region <- rbind(t_region, t_mab)

#Assign groupings based on the Smith et al paper
t$group <- ifelse(t$SVSPP %in% c(32, 851, 205, 859, 427, 430, 30, 31, 44),  "Clupeid", NA)
t$group <- ifelse(t$SVSPP %in% c(121, 578, 572, 209, 124, 898, 570, 702, 571, 582, 744, 860, 208, 212, 211, 874, 822, 745, 577, 877), "Scombrid", t$group)
t$group <- ifelse(t$SVSPP %in% c(181, 734, 38), "Sandlance", t$group)
t$group <- ifelse((t$SVSPP >= 501 & t$SVSPP <= 505) | (t$SVSPP >= 510 & t$SVSPP <= 513), "Squid_BT", t$group)
t$group <- ifelse(t$SVSPP == 297 | t$SVSPP == 298 | t$SVSPP == 306 | t$SVSPP == 307, "Shrimp_BT", t$group)
t$group <- ifelse(t$SVSPP == 306, "Northern_shrimp", t$group)
t$group <- ifelse(t$SVSPP %in% c(73:75), "Large_Gadid", t$group)
t$group <- ifelse(t$SVSPP %in% c(189, 1, 262, 192, 131, 750, 176, 168, 249, 193, 155, 164), "Misc_Fish", t$group) 
t$group <- ifelse(t$SVSPP %in% c(81, 454, 80, 86, 79, 69, 77, 72, 455, 78, 76), "Small_Gadid", t$group)
# t$group <- ifelse(t$SVSPP %in% c(790, 792, 110, 793, 777, 100, 104, 785, 786, 788, 109, 795, 775, 776, 773, 791, 787, 117, 789, 783, 103, 853, 774, 873, 106, 107, 105), "Flatfish", t$group)
t$group <- ifelse(t$SVSPP %in% c(51, 54, 56, 55, 229) , "Meso_pelagic", t$group)

#Assign which prey groupings belong to which predators - top 3 or 4
t$bp_prey <- ifelse(t$group == "Shrimp_BT" | t$group == "Clupeid" |t$group == "Sandlance", 1, 0) # also zoop
t$mn_prey <- ifelse(t$group == "Shrimp_BT" | t$group == "Misc_Fish" |t$group == "Sandlance", 1, 0)
t$ba_prey <- ifelse(t$group == "Shrimp_BT" | t$group == "Clupeid" |t$group == "Sandlance" | t$group == "Large_Gadid", 1, 0)
t$bb_prey <- ifelse(t$group == "Shrimp_BT" | t$group == "Clupeid", 1, 0) #also zoop

t <- merge(t, SVSPP_table, by.x = "SVSPP", by.y = "survdat.SVSPP", all.x=T)

#Do the same as above but for the by_region dataset
t_region$group <- ifelse(t_region$SVSPP %in% c(32, 851, 205, 859, 427, 430, 30, 31, 44),  "Clupeid", NA)
t_region$group <- ifelse(t_region$SVSPP %in% c(121, 578, 572, 209, 124, 898, 570, 702, 571, 582, 744, 860, 208, 212, 211, 874, 822, 745, 577, 877), "Scombrid", t_region$group)
t_region$group <- ifelse(t_region$SVSPP %in% c(181, 734, 38), "Sandlance", t_region$group)
t_region$group <- ifelse((t_region$SVSPP >= 501 & t_region$SVSPP <= 505) | (t_region$SVSPP >= 510 & t_region$SVSPP <= 513), "Squid_BT", t_region$group)
t_region$group <- ifelse(t_region$SVSPP == 297 | t_region$SVSPP == 298 | t_region$SVSPP == 306 | t_region$SVSPP == 307, "Shrimp_BT", t_region$group)
t_region$group <- ifelse(t_region$SVSPP %in% c(73:75), "Large_Gadid", t_region$group)
t_region$group <- ifelse(t_region$SVSPP %in% c(189, 1, 262, 192, 131, 750, 176, 168, 249, 193, 155, 164), "Misc_Fish", t_region$group) 
t_region$group <- ifelse(t_region$SVSPP %in% c(81, 454, 80, 86, 79, 69, 77, 72, 455, 78, 76), "Small_Gadid", t_region$group)
t_region$group <- ifelse(t_region$SVSPP %in% c(790, 792, 110, 793, 777, 100, 104, 785, 786, 788, 109, 795, 775, 776, 773, 791, 787, 117, 789, 783, 103, 853, 774, 873, 106, 107, 105), "Flatfish", t_region$group)
t_region$group <- ifelse(t_region$SVSPP %in% c(51, 54, 56, 55, 229) , "Meso_pelagic", t_region$group)

t_region$bp_prey <- ifelse(t_region$group == "Shrimp_BT" | t_region$group == "Clupeid" |t_region$group == "Sandlance" | t_region$group == "Scombrid", 1, 0) # also zoop
t_region$mn_prey <- ifelse(t_region$group == "Shrimp_BT" | t_region$group == "Misc_Fish"| t_region$group == "Scombrid" | t_region$group == "Clupeid" |t_region$group == "Sandlance", 1, 0)
t_region$ba_prey <- ifelse(t_region$group == "Shrimp_BT" | t_region$group == "Clupeid" |t_region$group == "Sandlance" | t_region$group == "Large_Gadid", 1, 0)
t_region$bb_prey <- ifelse(t_region$group == "Shrimp_BT" | t_region$group == "Clupeid" | t_region$group == "Squid_BT", 1, 0) #also zoop

t_region$sig <- ifelse(t_region$group == "Flatfish" | t_region$group == "Large_Gadid" | t_region$group == "Misc_Fish" | t_region$group == "shrimp" | t_region$group == "Small_Gadid", 1, 0)

#Merge back with the table of SVSPP numbers and common names to plot by common name
t_region <- merge(t_region, SVSPP_table, by.x = "SVSPP", by.y = "survdat.SVSPP", all.x=T)

# add sst
t_region <- merge(t_region, rast.avg_entire_reg_melt, by.x=c("YEAR", "stratcat"), by.y=c("YEAR", "variable"), all.x=T)
t <- merge(t, rast.avg_entire, by.x=c("YEAR"), by.y=c("YEAR"), all.x=T)

#extract regional BT
t_region <- merge(t_region, bt_region, by.x=c("YEAR", "stratcat"), by.y=c("YEAR", "stratcat"))


##### Boop over to the scaled biomass calculation to scale the biomass by proportion of diet for each species

########## ECOMON ###########
#### BRING IN ZOOPS ----------------------------
#First have to munge all this together to match the bottom trawl data to apply the survdat functions to it
ZPD_fall_sf$STRATUM <- ZPD_fall_sf$Name
ZPD_fall_sf$SEASON <- "FALL"
ZPD_fall_sf$STATION <- ZPD_fall_sf$station
ZPD_fall_sf$CRUISE6 <- ZPD_fall_sf$cruise_name
ZPD_fall_sf$YEAR <- ZPD_fall_sf$year
ZPD_fall_sf$BIOMASS <- ZPD_fall_sf$value
ZPD_fall_sf$ABUNDANCE <- ZPD_fall_sf$value
ZPD_fall_sf$CATCHSEX <- NA

#Calculate stratified mean for the whole region, then by each sub region
zpd_sel <- ZPD_fall_sf %>% dplyr::select(group_multi, STRATUM, SEASON, STATION, CRUISE6, YEAR, BIOMASS, ABUNDANCE, CATCHSEX, stratcat)
zpd_sel <- st_drop_geometry(cbind(zpd_sel, st_coordinates(zpd_sel)))
colnames(zpd_sel)[11:12] <- c("LON", "LAT")

zpd_strat_mean <- calc_stratified_mean(data.table(zpd_sel), areaDescription = "Name", groupDescription = "group_multi", filterBySeason = "FALL", areaPolygon = strat_zp)

zpd_strat_mean_gom <- calc_stratified_mean(data.table(zpd_sel[zpd_sel$stratcat == "GOMstrat",]), areaDescription = "Name", groupDescription = "group_multi", filterBySeason = "FALL", areaPolygon = strat_zp)

zpd_strat_mean_gb <- calc_stratified_mean(data.table(zpd_sel[zpd_sel$stratcat == "GBstrat",]), areaDescription = "Name", groupDescription = "group_multi", filterBySeason = "FALL", areaPolygon = strat_zp)

zpd_strat_mean_sne <- calc_stratified_mean(data.table(zpd_sel[zpd_sel$stratcat == "SNEstrat",]), areaDescription = "Name", groupDescription = "group_multi", filterBySeason = "FALL", areaPolygon = strat_zp)

zpd_strat_mean_mab <- calc_stratified_mean(data.table(zpd_sel[zpd_sel$stratcat == "MABstrat",]), areaDescription = "Name", groupDescription = "group_multi", filterBySeason = "FALL", areaPolygon = strat_zp)


#### BIG ECO
bigzpd_sel <- zpd_sel
bigzpd_sel$isgroup <- ifelse(!is.na(bigzpd_sel$group_multi), "group", "not")
bigzpd_strat_mean_gom <- calc_stratified_mean(data.table(bigzpd_sel[bigzpd_sel$stratcat == "GOMstrat",]), areaDescription = "Name", groupDescription = "isgroup", filterBySeason = "FALL", areaPolygon = strat_zp)

bigzpd_strat_mean_gb <- calc_stratified_mean(data.table(bigzpd_sel[bigzpd_sel$stratcat == "GBstrat",]), areaDescription = "Name", groupDescription = "isgroup", filterBySeason = "FALL", areaPolygon = strat_zp)

bigzpd_strat_mean_sne <- calc_stratified_mean(data.table(bigzpd_sel[bigzpd_sel$stratcat == "SNEstrat",]), areaDescription = "Name", groupDescription = "isgroup", filterBySeason = "FALL", areaPolygon = strat_zp)

bigzpd_strat_mean_mab <- calc_stratified_mean(data.table(bigzpd_sel[bigzpd_sel$stratcat == "MABstrat",]), areaDescription = "Name", groupDescription = "isgroup", filterBySeason = "FALL", areaPolygon = strat_zp)

#Assign subregion names
zpd_strat_mean_gom$stratcat <- "GOMstrat"
zpd_strat_mean_gb$stratcat <- "GBstrat"
zpd_strat_mean_sne$stratcat <- "SNEstrat"
zpd_strat_mean_mab$stratcat <- "MABstrat"

#Bind all the subregions together
zpd_region <- rbind(zpd_strat_mean_gom, zpd_strat_mean_gb)
zpd_region <- rbind(zpd_region, zpd_strat_mean_sne)
zpd_region <- rbind(zpd_region, zpd_strat_mean_mab)

# add sst
zpd_region <- merge(zpd_region, rast.avg_entire_reg_melt, by.x=c("YEAR", "stratcat"), by.y=c("YEAR", "variable"), all.x=T)
zpd_strat_mean <- merge(zpd_strat_mean, rast.avg_entire, by.x="YEAR", by.y="YEAR", all.x=T)

#extract regional BT
zpd_region <- merge(zpd_region, bt_region, by.x=c("YEAR", "stratcat"), by.y=c("YEAR", "stratcat"))

##### COMBINE #####
## Combine for all regions
zpd_strat_mean_sub <- zpd_strat_mean %>% 
  filter(group_multi != "NA") %>% 
  dplyr::select(YEAR, strat.biomass, group_multi, t.rast.avg_entire.)
t_sub <- t %>% 
  filter(group != "NA") %>% 
  dplyr::select(YEAR, strat.biomass, group, t.rast.avg_entire.)

colnames(zpd_strat_mean_sub)[3] <- "group"

stratbio_allregion <- data.frame(rbind(zpd_strat_mean_sub, t_sub))

zpd_region_comb <- zpd_region %>% 
  filter(group_multi != "NA") %>% 
  dplyr::select(YEAR, strat.biomass, group_multi, stratcat, SURFTEMP)
t_region_comb <- t_region %>% 
  filter(group != "NA") %>% 
  dplyr::select(YEAR, strat.biomass, group, stratcat, SURFTEMP)

colnames(zpd_region_comb)[3] <- "group"

stratbio_comb <- data.frame(rbind(zpd_region_comb, t_region_comb))

## standardize biomass
## for 1990 - 2019, calculate mean and sd. then do (value - mean) / sd
std_ano_bio <- NULL
for(i in unique(stratbio_comb$group)){
  ww <- stratbio_comb[stratbio_comb$group == i & stratbio_comb$YEAR >= 1980,]
  for(j in unique(ww$stratcat)){
    www <- ww[ww$stratcat == j,]
    mm <- mean(www$strat.biomass)
    maxx <- max(www$strat.biomass)
    ss <- sd(www$strat.biomass)
    www$std_ano <- (www$strat.biomass - mm) / ss
    www$prop <- www$strat.biomass / maxx
    std_ano_bio <- rbind(std_ano_bio, www)
  }
}

std_ano_bio$biggrp <- ifelse(std_ano_bio$group == "Clupeid" | std_ano_bio$group == "Flatfish" | std_ano_bio$group == "Meso_pelagic" | std_ano_bio$group == "Misc_Fish" | std_ano_bio$group == "Sandlance" | std_ano_bio$group == "Scombrid" | std_ano_bio$group == "Small_Gadid", "BT Fish", NA)
std_ano_bio$biggrp <- ifelse(std_ano_bio$group == "Northern_shrimp" | std_ano_bio$group == "Shrimp_BT" | std_ano_bio$group == "Squid_BT", "BT Squid and Shrimp", std_ano_bio$biggrp)
std_ano_bio$biggrp <- ifelse(std_ano_bio$group == "euphausiid" | std_ano_bio$group == "fish" | std_ano_bio$group == "lg_zoop" | std_ano_bio$group == "shrimp" | std_ano_bio$group == "sm_zoop", "Zoops and Shrimps", std_ano_bio$biggrp)

std_ano_bio %>% 
  #filter(group == "Zoops and Shrimps") %>%
  ggplot( aes(x=YEAR, y=group, height=prop, group=group, color=group, fill=group)) + 
  #geom_point(size=0.5, alpha = 0.7) + 
  #geom_line() +
  #geom_ma(ma_fun = SMA, n=5, size=1) +
  geom_density_ridges(stat="identity", scale=1) +
  geom_smooth(method = "loess", se=F, span=0.8) + 
  facet_wrap(.~stratcat, scales = "free_y")

summary(lm(std_ano ~ YEAR*group, std_ano_bio[std_ano_bio$stratcat == "GOMstrat",]))





### group for whales
survdat_fall$Bp_food <- ifelse(survdat_fall$group == "Shrimp_BT" | survdat_fall$group == "Northern_shrimp" |survdat_fall$group == "euphausiid" | survdat_fall$group == "lg_zoop" | survdat_fall$group == "sm_zoop" | survdat_fall$group == "shrimp" | survdat_fall$group == "fish" | survdat_fall$group == "Clupeid" | survdat_fall$group == "Sandlance" , "Bp_food", NA)

survdat_fall$Mn_food <- ifelse(survdat_fall$group == "Shrimp_BT" | survdat_fall$group == "Northern_shrimp" | survdat_fall$group == "shrimp" | survdat_fall$group == "Misc_Fish" | survdat_fall$group == "Clupeid" | survdat_fall$group == "Sandlance" , "Mn_food", survdat_fall$whgrp)

survdat_fall$Bb_food <- ifelse(survdat_fall$group == "Shrimp_BT" | survdat_fall$group == "Northern_shrimp" |survdat_fall$group == "euphausiid" | survdat_fall$group == "lg_zoop" | survdat_fall$group == "sm_zoop" | survdat_fall$group == "shrimp" | survdat_fall$group == "fish", "Bb_food", survdat_fall$whgrp)

survdat_fall$Ba_food <- ifelse(survdat_fall$group == "Shrimp_BT" | survdat_fall$group == "Northern_shrimp" |survdat_fall$group == "Sandlance" | survdat_fall$group == "Large_Gadid"| survdat_fall$group == "shrimp" | survdat_fall$group == "Clupeid", "Ba_food", survdat_fall$whgrp)

##### eco

survdat_fall$Bp_food <- ifelse(survdat_fall$group == "euphausiid" | survdat_fall$group == "lg_zoop" | survdat_fall$group == "sm_zoop" | survdat_fall$group == "shrimp" | survdat_fall$group == "fish" , "Bp_food", NA)

survdat_fall$Mn_food <- ifelse( survdat_fall$group == "shrimp" , "Mn_food", survdat_fall$whgrp)

survdat_fall$Bb_food <- ifelse(survdat_fall$group == "euphausiid" | survdat_fall$group == "lg_zoop" | survdat_fall$group == "sm_zoop" | survdat_fall$group == "shrimp" | survdat_fall$group == "fish", "Bb_food", survdat_fall$whgrp)

survdat_fall$Ba_food <- ifelse( survdat_fall$group == "shrimp", "Ba_food", survdat_fall$whgrp)

survdat_fall_wh <- survdat_fall[,-23]
survdat_fall_wh_bp <- survdat_fall[,-c(23,25, 27, 28, 29)]
survdat_fall_wh_bp$whales <- "Bp"
colnames(survdat_fall_wh_bp)[24] <- "food"
survdat_fall_wh_mn <- survdat_fall[,-c(23,25, 26, 28, 29)]
survdat_fall_wh_mn$whales <- "Mn"
colnames(survdat_fall_wh_mn)[24] <- "food"
survdat_fall_wh_bb <- survdat_fall[,-c(23,25, 26, 27, 29)]
survdat_fall_wh_bb$whales <- "Bb"
colnames(survdat_fall_wh_bb)[24] <- "food"
survdat_fall_wh_ba <- survdat_fall[,-c(23,25, 26, 27, 28)]
survdat_fall_wh_ba$whales <- "Ba"
colnames(survdat_fall_wh_ba)[24] <- "food"

survdat_fall_wh <- rbind(survdat_fall_wh_bp, survdat_fall_wh_ba)
survdat_fall_wh <- rbind(survdat_fall_wh, survdat_fall_wh_bb)
survdat_fall_wh <- rbind(survdat_fall_wh, survdat_fall_wh_mn)

t_wh <- calc_stratified_mean(data.table(survdat_fall_wh), areaDescription = "STRATA", groupDescription = "food", filterBySeason = "FALL")
t_wh_gom <- calc_stratified_mean(data.table(survdat_fall_wh[survdat_fall_wh$stratcat == "GOMstrat",-c(25)]), areaDescription = "STRATA", groupDescription = "food", filterBySeason = "FALL")
t_wh_gb <- calc_stratified_mean(data.table(survdat_fall_wh[survdat_fall_wh$stratcat == "GBstrat",-c(25)]), areaDescription = "STRATA", groupDescription = "food", filterBySeason = "FALL")
t_wh_sne <- calc_stratified_mean(data.table(survdat_fall_wh[survdat_fall_wh$stratcat == "SNEstrat",-c(25)]), areaDescription = "STRATA", groupDescription = "food", filterBySeason = "FALL")
t_wh_mab <- calc_stratified_mean(data.table(survdat_fall_wh[survdat_fall_wh$stratcat == "MABstrat",-c(25)]), areaDescription = "STRATA", groupDescription = "food", filterBySeason = "FALL")

#Assign category for each region
t_wh_gom$stratcat <- "GOMstrat"
t_wh_gb$stratcat <- "GBstrat"
t_wh_sne$stratcat <- "SNEstrat"
t_wh_mab$stratcat <- "MABstrat"

#Combine regions 
t_wh_region <- rbind(t_wh_gom, t_wh_gb)
t_wh_region <- rbind(t_wh_region, t_wh_sne)
t_wh_region <- rbind(t_wh_region, t_wh_mab)

ggplot(t_wh_region[!is.na(t_wh_region$food) & t_wh_region$YEAR >= 1980,], aes(x=YEAR, y=strat.biomass, color=food)) + geom_point() + geom_smooth() + facet_wrap(.~stratcat, scales="free_y")

## standardize biomass for whales
## for 1990 - 2019, calculate mean and sd. then do (value - mean) / sd
std_ano_bio_wh <- NULL
for(i in unique(t_wh_region$food)){
  ww <- t_wh_region[t_wh_region$food == i & t_wh_region$YEAR >= 1980,]
  for(j in unique(ww$stratcat)){
    www <- ww[ww$stratcat == j,]
    mm <- mean(www$strat.biomass)
    ss <- sd(www$strat.biomass)
    www$std_ano <- (www$strat.biomass - mm) / ss
    std_ano_bio_wh <- rbind(std_ano_bio_wh, www)
  }
}

ggplot(std_ano_bio_wh, aes(x=YEAR, y=std_ano, group=food, color=food)) + geom_point() + geom_smooth(se=T) + facet_wrap(.~stratcat, scales = "free_y")



###create table of groups, spp numbers, common names
fish_names <- unique(t_region[,c(3,18,12)])
zp_names <- unique(ZPD_melt[,c(27,18,21)])

colnames(zp_names) <- colnames(fish_names)
critter_names <- data.frame(rbind(fish_names, zp_names))


zpd_region1 <- merge(zpd_region, critter_names[,2:3], by.x=c("TAXA.NAME"), by.y="survdat.COMNAME")
zpd_region1 <- aggregate(strat.biomass ~ stratcat + YEAR + group, zpd_region1, FUN=sum)
#Now to calculate the lm results for each of these groupings
corrtest_whale_bio <- NULL
for(i in unique(dual_comb_anom$stratcat)){
  w <- dual_comb_anom[dual_comb_anom$stratcat == i & dual_comb_anom$YEAR >= 1977 & dual_comb_anom$YEAR <= 2019,]
  w$new <- ifelse(!is.na(w$zoop_abundance), w$zoop_abundance, w$strat.biomass)
  for(j in unique(w$group)){
    o <- w[w$group == j,]
    if(nrow(o) > 0){
      r <- lm(o$stdano ~ o$YEAR)
      p <- c(i, j, summary(r)$coefficients[8], r$coefficients[2], summary(r)$r.squared)
      corrtest_whale_bio <- rbind(corrtest_whale_bio, p)
    }
  }
}
corrtest_whale_bio <- data.frame(corrtest_whale_bio)
colnames(corrtest_whale_bio) <-  c("reg", "group", "pval", "slope", "r2")
corrtest_whale_bio$pval <- as.numeric(corrtest_whale_bio$pval)
corrtest_whale_bio$slope <- as.numeric(corrtest_whale_bio$slope)
corrtest_whale_bio$r2 <- as.numeric(corrtest_whale_bio$r2)
corrtest_whale_bio$sig <- ifelse(corrtest_whale_bio$pval <= 0.05, 1, 0) #is it sig @ 0.05?

