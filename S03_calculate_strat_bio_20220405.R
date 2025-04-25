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
t_gom <- calc_stratified_mean(data.table(survdat_fall[survdat_fall$stratcat == "GOMstrat",]), areaDescription = "STRATA", groupDescription = "group", filterBySeason = "FALL")
t_gb <- calc_stratified_mean(data.table(survdat_fall[survdat_fall$stratcat == "GBstrat",]), areaDescription = "STRATA", groupDescription = "group", filterBySeason = "FALL")
t_sne <- calc_stratified_mean(data.table(survdat_fall[survdat_fall$stratcat == "SNEstrat",]), areaDescription = "STRATA", groupDescription = "group", filterBySeason = "FALL")
t_mab <- calc_stratified_mean(data.table(survdat_fall[survdat_fall$stratcat == "MABstrat",]), areaDescription = "STRATA", groupDescription = "group", filterBySeason = "FALL")

#Assign category for each region
t_gom$stratcat <- "GOMstrat"
t_gb$stratcat <- "GBstrat"
t_sne$stratcat <- "SNEstrat"
t_mab$stratcat <- "MABstrat"

#Combine regions
t_region <- rbind(t_gom, t_gb)
t_region <- rbind(t_region, t_sne)
t_region <- rbind(t_region, t_mab)

##### ALL BIOMASS PER REGION
#Calculate stratified mean for all regions combined, but then also do it by strata regions

big_survdat_fall <- survdat_fall
big_survdat_fall$isprey <- ifelse(!is.na(big_survdat_fall$group), "prey", "not_prey")
big_survdat_fall <- big_survdat_fall[,!("group")]

bigt_gom <- calc_stratified_mean(data.table(big_survdat_fall[big_survdat_fall$stratcat == "GOMstrat",]), areaDescription = "STRATA", groupDescription = "isprey", filterBySeason = "FALL")
bigt_gb <- calc_stratified_mean(data.table(big_survdat_fall[big_survdat_fall$stratcat == "GBstrat",]), areaDescription = "STRATA", groupDescription = "isprey",  filterBySeason = "FALL")
bigt_sne <- calc_stratified_mean(data.table(big_survdat_fall[big_survdat_fall$stratcat == "SNEstrat",]), areaDescription = "STRATA", groupDescription = "isprey", filterBySeason = "FALL")
bigt_mab <- calc_stratified_mean(data.table(big_survdat_fall[big_survdat_fall$stratcat == "MABstrat",]), areaDescription = "STRATA", groupDescription = "isprey", filterBySeason = "FALL")

#Assign category for each region
bigt_gom$stratcat <- "GOMstrat"
bigt_gb$stratcat <- "GBstrat"
bigt_sne$stratcat <- "SNEstrat"
bigt_mab$stratcat <- "MABstrat"

#Combine regions
bigt_region <- rbind(bigt_gom, bigt_gb)
bigt_region <- rbind(bigt_region, bigt_sne)
bigt_region <- rbind(bigt_region, bigt_mab)

### TABLE 3 ###
#Now to calculate the lm results for each of these groupings
corrtest_fishgroup_bio <- NULL
for(i in unique(bigt_region$stratcat)){
  w <- bigt_region[bigt_region$stratcat == i & bigt_region$YEAR >= 1980 & bigt_region$YEAR <= 2019,]
  # w$new <- ifelse(!is.na(w$zoop_abundance), w$zoop_abundance, w$strat.biomass)
  for(j in unique(w$isprey)){
    o <- w[w$isprey == j,]
    if(nrow(o) > 0){
      r <- lm(o$strat.biomass ~ o$YEAR)
      p <- c(i, j, summary(r)$coefficients[8], r$coefficients[2], summary(r)$r.squared)
      corrtest_fishgroup_bio <- rbind(corrtest_fishgroup_bio, p)
    }
  }
}
corrtest_fishgroup_bio <- data.frame(corrtest_fishgroup_bio)
colnames(corrtest_fishgroup_bio) <-  c("reg", "group", "pval", "slope", "r2")
corrtest_fishgroup_bio$pval <- as.numeric(corrtest_fishgroup_bio$pval)
corrtest_fishgroup_bio$slope <- as.numeric(corrtest_fishgroup_bio$slope)
corrtest_fishgroup_bio$r2 <- as.numeric(corrtest_fishgroup_bio$r2)
corrtest_fishgroup_bio$sig <- ifelse(corrtest_fishgroup_bio$pval <= 0.0253, 1, 0) #is it sig @ 0.05?


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

#Assign subregion names
zpd_strat_mean_gom$stratcat <- "GOMstrat"
zpd_strat_mean_gb$stratcat <- "GBstrat"
zpd_strat_mean_sne$stratcat <- "SNEstrat"
zpd_strat_mean_mab$stratcat <- "MABstrat"

#Bind all the subregions together
zpd_region <- rbind(zpd_strat_mean_gom, zpd_strat_mean_gb)
zpd_region <- rbind(zpd_region, zpd_strat_mean_sne)
zpd_region <- rbind(zpd_region, zpd_strat_mean_mab)


#### BIG ECO
bigzpd_sel <- zpd_sel
bigzpd_sel$isgroup <- ifelse(!is.na(bigzpd_sel$group_multi), "group", "not")
bigzpd_strat_mean_gom <- calc_stratified_mean(data.table(bigzpd_sel[bigzpd_sel$stratcat == "GOMstrat",]), areaDescription = "Name", groupDescription = "isgroup", filterBySeason = "FALL", areaPolygon = strat_zp)

bigzpd_strat_mean_gb <- calc_stratified_mean(data.table(bigzpd_sel[bigzpd_sel$stratcat == "GBstrat",]), areaDescription = "Name", groupDescription = "isgroup", filterBySeason = "FALL", areaPolygon = strat_zp)

bigzpd_strat_mean_sne <- calc_stratified_mean(data.table(bigzpd_sel[bigzpd_sel$stratcat == "SNEstrat",]), areaDescription = "Name", groupDescription = "isgroup", filterBySeason = "FALL", areaPolygon = strat_zp)

bigzpd_strat_mean_mab <- calc_stratified_mean(data.table(bigzpd_sel[bigzpd_sel$stratcat == "MABstrat",]), areaDescription = "Name", groupDescription = "isgroup", filterBySeason = "FALL", areaPolygon = strat_zp)

#Assign subregion names
bigzpd_strat_mean_gom$stratcat <- "GOMstrat"
bigzpd_strat_mean_gb$stratcat <- "GBstrat"
bigzpd_strat_mean_sne$stratcat <- "SNEstrat"
bigzpd_strat_mean_mab$stratcat <- "MABstrat"

#Bind all the subregions together
bigzpd_region <- rbind(bigzpd_strat_mean_gom, bigzpd_strat_mean_gb)
bigzpd_region <- rbind(bigzpd_region, bigzpd_strat_mean_sne)
bigzpd_region <- rbind(bigzpd_region, bigzpd_strat_mean_mab)

### TABLE 3 ###
#Now to calculate the lm results for each of these groupings
corrtest_zoopgroup_bio <- NULL
for(i in unique(bigzpd_region$stratcat)){
  w <- bigzpd_region[bigzpd_region$stratcat == i & bigzpd_region$YEAR >= 1980 & bigzpd_region$YEAR <= 2019,]
  # w$new <- ifelse(!is.na(w$zoop_abundance), w$zoop_abundance, w$strat.biomass)
  for(j in unique(w$isgroup)){
    o <- w[w$isgroup == j,]
    if(nrow(o) > 0){
      r <- lm(o$strat.biomass ~ o$YEAR)
      p <- c(i, j, summary(r)$coefficients[8], r$coefficients[2], summary(r)$r.squared)
      corrtest_zoopgroup_bio <- rbind(corrtest_zoopgroup_bio, p)
    }
  }
}
corrtest_zoopgroup_bio <- data.frame(corrtest_zoopgroup_bio)
colnames(corrtest_zoopgroup_bio) <-  c("reg", "group", "pval", "slope", "r2")
corrtest_zoopgroup_bio$pval <- as.numeric(corrtest_zoopgroup_bio$pval)
corrtest_zoopgroup_bio$slope <- as.numeric(corrtest_zoopgroup_bio$slope)
corrtest_zoopgroup_bio$r2 <- as.numeric(corrtest_zoopgroup_bio$r2)
corrtest_zoopgroup_bio$sig <- ifelse(corrtest_zoopgroup_bio$pval <= 0.0253, 1, 0) #is it sig @ 0.05?


##### COMBINE #####
## Combine for all regions
zpd_strat_mean_sub <- zpd_strat_mean %>%
  filter(group_multi != "NA") %>%
  dplyr::select(YEAR, strat.biomass, group_multi)
t_sub <- t %>%
  filter(group != "NA") %>%
  dplyr::select(YEAR, strat.biomass, group)

colnames(zpd_strat_mean_sub)[3] <- "group"

stratbio_allregion <- data.frame(rbind(zpd_strat_mean_sub, t_sub))

zpd_region_comb <- zpd_region %>%
  filter(group_multi != "NA") %>%
  dplyr::select(YEAR, strat.biomass, group_multi, stratcat)
t_region_comb <- t_region %>%
  filter(group != "NA") %>%
  dplyr::select(YEAR, strat.biomass, group, stratcat)

colnames(zpd_region_comb)[3] <- "group"

stratbio_comb <- data.frame(rbind(zpd_region_comb, t_region_comb))

### TABLE 2 ###
#Now to calculate the lm results for each of these groupings
corrtest_biomass <- NULL
for(i in unique(stratbio_comb$stratcat)){
  w <- stratbio_comb[stratbio_comb$stratcat == i & stratbio_comb$YEAR >= 1980 & stratbio_comb$YEAR <= 2019,]
  # w$new <- ifelse(!is.na(w$zoop_abundance), w$zoop_abundance, w$strat.biomass)
  for(j in unique(w$group)){
    o <- w[w$group == j,]
    if(nrow(o) > 0){
      r <- lm(o$strat.biomass ~ o$YEAR)
      p <- c(i, j, summary(r)$coefficients[8], r$coefficients[2], summary(r)$r.squared)
      corrtest_biomass <- rbind(corrtest_biomass, p)
    }
  }
}
corrtest_biomass <- data.frame(corrtest_biomass)
colnames(corrtest_biomass) <-  c("reg", "group", "pval", "slope", "r2")
corrtest_biomass$pval <- as.numeric(corrtest_biomass$pval)
corrtest_biomass$slope <- as.numeric(corrtest_biomass$slope)
corrtest_biomass$r2 <- as.numeric(corrtest_biomass$r2)
corrtest_biomass$sig <- ifelse(corrtest_biomass$pval <= 0.0037, 1, 0) #is it sig @ 0.05?

##### GROUP BY PREFERRED LARGE WHALE PREY ITEMS #########
### group for whales
survdat_fall$Bp_food <- ifelse(survdat_fall$group == "Shrimp_BT" | survdat_fall$group == "Northern_shrimp" |survdat_fall$group == "euphausiid" | survdat_fall$group == "lg_zoop" | survdat_fall$group == "sm_zoop" | survdat_fall$group == "shrimp" | survdat_fall$group == "fish" | survdat_fall$group == "Clupeid" | survdat_fall$group == "Sandlance" , "Bp_food", NA)

survdat_fall$Mn_food <- ifelse(survdat_fall$group == "Shrimp_BT" | survdat_fall$group == "Northern_shrimp" | survdat_fall$group == "shrimp" | survdat_fall$group == "Misc_Fish" | survdat_fall$group == "Clupeid" | survdat_fall$group == "Sandlance" , "Mn_food", NA)

survdat_fall$Bb_food <- ifelse(survdat_fall$group == "Shrimp_BT" | survdat_fall$group == "Northern_shrimp" |survdat_fall$group == "euphausiid" | survdat_fall$group == "lg_zoop" | survdat_fall$group == "sm_zoop" | survdat_fall$group == "shrimp" | survdat_fall$group == "fish", "Bb_food", NA)

survdat_fall$Ba_food <- ifelse(survdat_fall$group == "Shrimp_BT" | survdat_fall$group == "Northern_shrimp" |survdat_fall$group == "Sandlance" | survdat_fall$group == "Large_Gadid"| survdat_fall$group == "shrimp" | survdat_fall$group == "Clupeid", "Ba_food", NA)

##### eco

survdat_fall$Bp_food <- ifelse(survdat_fall$group == "euphausiid" | survdat_fall$group == "lg_zoop" | survdat_fall$group == "sm_zoop" | survdat_fall$group == "shrimp" | survdat_fall$group == "fish" , "Bp_food", survdat_fall$Bp_food)

survdat_fall$Mn_food <- ifelse( survdat_fall$group == "shrimp" , "Mn_food", survdat_fall$Mn_food)

survdat_fall$Bb_food <- ifelse(survdat_fall$group == "euphausiid" | survdat_fall$group == "lg_zoop" | survdat_fall$group == "sm_zoop" | survdat_fall$group == "shrimp" | survdat_fall$group == "fish", "Bb_food", survdat_fall$Bb_food)

survdat_fall$Ba_food <- ifelse( survdat_fall$group == "shrimp", "Ba_food", survdat_fall$Ba_food)

survdat_fall_wh <- survdat_fall
survdat_fall_wh_bp <- survdat_fall %>% 
  dplyr::select(-c("COMNAME", "ID", "Mn_food", "Bb_food", "Ba_food", "group")) #[,-c(24, 27, 28, 29)]
survdat_fall_wh_bp$whales <- "Bp"
colnames(survdat_fall_wh_bp)[colnames(survdat_fall_wh_bp) == "Bp_food"] <- "food"
survdat_fall_wh_mn <- survdat_fall %>% 
  dplyr::select(-c("COMNAME", "ID", "Bp_food", "Bb_food", "Ba_food", "group")) #[,-c(24, 26, 28, 29)]
survdat_fall_wh_mn$whales <- "Mn"
colnames(survdat_fall_wh_mn)[colnames(survdat_fall_wh_mn) == "Mn_food"] <- "food"
survdat_fall_wh_bb <- survdat_fall %>% 
  dplyr::select(-c("COMNAME", "ID", "Mn_food", "Bp_food", "Ba_food", "group"))  #[,-c(24, 26, 27, 29)]
survdat_fall_wh_bb$whales <- "Bb"
colnames(survdat_fall_wh_bb)[colnames(survdat_fall_wh_bb) == "Bb_food"] <- "food"
survdat_fall_wh_ba <- survdat_fall %>% 
  dplyr::select(-c("COMNAME", "ID", "Mn_food", "Bb_food", "Bp_food", "group")) #[,-c(24, 26, 27, 28)]
survdat_fall_wh_ba$whales <- "Ba"
colnames(survdat_fall_wh_ba)[colnames(survdat_fall_wh_ba) == "Ba_food"] <- "food"

survdat_fall_wh <- rbind(survdat_fall_wh_bp, survdat_fall_wh_ba)
survdat_fall_wh <- rbind(survdat_fall_wh, survdat_fall_wh_bb)
survdat_fall_wh <- rbind(survdat_fall_wh, survdat_fall_wh_mn)

t_wh <- calc_stratified_mean(data.table(survdat_fall_wh), areaDescription = "STRATA", groupDescription = "food", filterBySeason = "FALL")
t_wh_gom <- calc_stratified_mean(data.table(survdat_fall_wh[survdat_fall_wh$stratcat == "GOMstrat",]), areaDescription = "STRATA", groupDescription = "food", filterBySeason = "FALL")
t_wh_gb <- calc_stratified_mean(data.table(survdat_fall_wh[survdat_fall_wh$stratcat == "GBstrat",]), areaDescription = "STRATA", groupDescription = "food", filterBySeason = "FALL")
t_wh_sne <- calc_stratified_mean(data.table(survdat_fall_wh[survdat_fall_wh$stratcat == "SNEstrat",]), areaDescription = "STRATA", groupDescription = "food", filterBySeason = "FALL")
t_wh_mab <- calc_stratified_mean(data.table(survdat_fall_wh[survdat_fall_wh$stratcat == "MABstrat",]), areaDescription = "STRATA", groupDescription = "food", filterBySeason = "FALL")

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



##create table of groups, spp numbers, common names
fish_names <- unique(survdat_fall %>% dplyr::select(SVSPP, COMNAME, group))
zp_names <- unique(ZPD_melt %>% dplyr::select(taxa_edited, TAXA.NAME, smith_group))

colnames(zp_names) <- colnames(fish_names)
critter_names <- data.frame(rbind(fish_names, zp_names))

### TABLE 4 ###
#Now to calculate the lm results for each of these groupings
corrtest_whale_bio <- NULL
for(i in unique(t_wh_region$stratcat)){
  w <- t_wh_region[t_wh_region$stratcat == i & t_wh_region$YEAR >= 1980 & t_wh_region$YEAR <= 2019,]
  # w$new <- ifelse(!is.na(w$zoop_abundance), w$zoop_abundance, w$strat.biomass)
  for(j in unique(w$food)){
    o <- w[w$food == j,]
    if(nrow(o) > 0){
      r <- lm(o$strat.biomass ~ o$YEAR)
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
corrtest_whale_bio$sig <- ifelse(corrtest_whale_bio$pval <= 0.013, 1, 0) #is it sig @ 0.05?

### Do this in 5-year increments for Figure 5 ###
survdat_fall_clusteryear <- survdat_fall[survdat_fall$YEAR >= 1980,]
survdat_fall_clusteryear$CATCHSEX <- 0
survdat_fall_clusteryear$YEAR <- survdat_fall_clusteryear$yr5
#Calculate stratified mean for all regions combined, but then also do it by strata regions
t_5 <- calc_stratified_mean(data.table(survdat_fall_clusteryear), areaDescription = "STRATA", groupDescription = "group", filterBySeason = "FALL")
t_gom_5 <- calc_stratified_mean(data.table(survdat_fall_clusteryear[survdat_fall_clusteryear$stratcat == "GOMstrat",]), areaDescription = "STRATA", groupDescription = "group", filterBySeason = "FALL")
t_gb_5 <- calc_stratified_mean(data.table(survdat_fall_clusteryear[survdat_fall_clusteryear$stratcat == "GBstrat",]), areaDescription = "STRATA", groupDescription = "group", filterBySeason = "FALL")
t_sne_5 <- calc_stratified_mean(data.table(survdat_fall_clusteryear[survdat_fall_clusteryear$stratcat == "SNEstrat",]), areaDescription = "STRATA", groupDescription = "group", filterBySeason = "FALL")
t_mab_5 <- calc_stratified_mean(data.table(survdat_fall_clusteryear[survdat_fall_clusteryear$stratcat == "MABstrat",]), areaDescription = "STRATA", groupDescription = "group", filterBySeason = "FALL")

#Assign category for each region
t_gom_5$stratcat <- "GOMstrat"
t_gb_5$stratcat <- "GBstrat"
t_sne_5$stratcat <- "SNEstrat"
t_mab_5$stratcat <- "MABstrat"

#Combine regions 
t_region_5 <- rbind(t_gom_5, t_gb_5)
t_region_5 <- rbind(t_region_5, t_sne_5)
t_region_5 <- rbind(t_region_5, t_mab_5)

t_region_5$bp_prey <- ifelse(t_region_5$group == "Shrimp_BT" | t_region_5$group == "Clupeid" |t_region_5$group == "Sandlance" , "Preferred Prey", "Other Prey") # also zoop
t_region_5$mn_prey <- ifelse(t_region_5$group == "Shrimp_BT" | t_region_5$group == "Misc_Fish"| t_region_5$group == "Clupeid" |t_region_5$group == "Sandlance", "Preferred Prey", "Other Prey")
t_region_5$ba_prey <- ifelse(t_region_5$group == "Shrimp_BT" | t_region_5$group == "Clupeid" |t_region_5$group == "Sandlance" | t_region_5$group == "Large_Gadid", "Preferred Prey", "Other Prey")
t_region_5$bb_prey <- ifelse(t_region_5$group == "Shrimp_BT" | t_region_5$group == "Clupeid" | t_region_5$group == "Squid_BT", "Preferred Prey", "Other Prey") #also zoop

t_region_5$group2 <- t_region_5$group
t_region_5$group2 <- ifelse(t_region_5$group2 == "Large_Gadid", "Large Gadid", t_region_5$group2)
t_region_5$group2 <- ifelse(t_region_5$group2 == "Meso_pelagic", "Mesopelagic Fish", t_region_5$group2)
t_region_5$group2 <- ifelse(t_region_5$group2 == "Misc_Fish", "Miscellaneous Fish", t_region_5$group2)
t_region_5$group2 <- ifelse(t_region_5$group2 == "Northern_shrimp", "Northern Shrimp", t_region_5$group2)
t_region_5$group2 <- ifelse(t_region_5$group2 == "Shrimp_BT", "Shrimp", t_region_5$group2)
t_region_5$group2 <- ifelse(t_region_5$group2 == "Small_Gadid", "Small Gadid", t_region_5$group2)
t_region_5$group2 <- ifelse(t_region_5$group2 == "Squid_BT", "Squid", t_region_5$group2)


t_region2 <- t_region_5
t_region2$fin_prey <- ifelse(t_region2$bp_prey == "Preferred Prey", t_region2$strat.biomass, NA)
t_region2$humpback_prey <- ifelse(t_region2$mn_prey == "Preferred Prey", t_region2$strat.biomass, NA)
t_region2$minke_prey <- ifelse(t_region2$ba_prey == "Preferred Prey", t_region2$strat.biomass, NA)
t_region2$sei_prey <- ifelse(t_region2$bb_prey == "Preferred Prey", t_region2$strat.biomass, NA)


p <- aggregate(fin_prey ~ YEAR + stratcat, data=t_region2, FUN=sum)
s <- aggregate(minke_prey ~ YEAR + stratcat, data=t_region2, FUN=sum)
pp <- aggregate(humpback_prey ~ YEAR + stratcat, data=t_region2, FUN=sum)
ss <- aggregate(sei_prey ~ YEAR + stratcat, data=t_region2, FUN=sum)

t_region2 <- merge(p, s, by.x=c("YEAR", "stratcat"), by.y=c("YEAR", "stratcat"))
t_region2 <- merge(t_region2, pp, by.x=c("YEAR", "stratcat"), by.y=c("YEAR", "stratcat"))
t_region2 <- merge(t_region2, ss, by.x=c("YEAR", "stratcat"), by.y=c("YEAR", "stratcat"))
t_region2 <- melt(t_region2, id.vars=c("YEAR", "stratcat"))
