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

survdat_fall_clusteryear <- survdat_fall[survdat_fall$YEAR >= 1980,]

survdat_fall_clusteryear$group <- ifelse(survdat_fall_clusteryear$SVSPP %in% c(51, 54, 56, 55, 229) , "Meso_pelagic", survdat_fall_clusteryear$group)
survdat_fall_clusteryear$CATCHSEX <- 0
survdat_fall_clusteryear$YEAR <- survdat_fall_clusteryear$yr5
#Calculate stratified mean for all regions combined, but then also do it by strata regions
t <- calc_stratified_mean(data.table(survdat_fall_clusteryear), areaDescription = "STRATA", groupDescription = "group", filterBySeason = "FALL")
t_gom <- calc_stratified_mean(data.table(survdat_fall_clusteryear[survdat_fall_clusteryear$stratcat == "GOMstrat",-c(19)]), areaDescription = "STRATA", groupDescription = "group", filterBySeason = "FALL")
t_gb <- calc_stratified_mean(data.table(survdat_fall_clusteryear[survdat_fall_clusteryear$stratcat == "GBstrat",-c(19)]), areaDescription = "STRATA", groupDescription = "group", filterBySeason = "FALL")
t_sne <- calc_stratified_mean(data.table(survdat_fall_clusteryear[survdat_fall_clusteryear$stratcat == "SNEstrat",-c(19)]), areaDescription = "STRATA", groupDescription = "group", filterBySeason = "FALL")
t_mab <- calc_stratified_mean(data.table(survdat_fall_clusteryear[survdat_fall_clusteryear$stratcat == "MABstrat",-c(19)]), areaDescription = "STRATA", groupDescription = "group", filterBySeason = "FALL")

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

## standardize biomass
## for 1990 - 2019, calculate mean and sd. then do (value - mean) / sd
std_ano_bio_5yr <- NULL
for(i in unique(t_region$group)){
  ww <- t_region[t_region$group == i,]
  for(j in unique(ww$stratcat)){
    www <- ww[ww$stratcat == j,]
    mm <- mean(www$strat.biomass)
    ss <- sd(www$strat.biomass)
    maxx <- max(www$strat.biomass)
    www$std_ano <- (www$strat.biomass - mm) / ss
    www$prop <- www$strat.biomass / maxx
    std_ano_bio_5yr <- rbind(std_ano_bio_5yr, www)
  }
}

std_ano_bio_5yr$biggrp <- ifelse(std_ano_bio_5yr$group == "Clupeid" | std_ano_bio_5yr$group == "Flatfish" | std_ano_bio_5yr$group == "Meso_pelagic" | std_ano_bio_5yr$group == "Misc_Fish" | std_ano_bio_5yr$group == "Sandlance" | std_ano_bio_5yr$group == "Scombrid" | std_ano_bio_5yr$group == "Small_Gadid", "BT Fish", NA)
std_ano_bio_5yr$biggrp <- ifelse(std_ano_bio_5yr$group == "Northern_shrimp" | std_ano_bio_5yr$group == "Shrimp_BT" | std_ano_bio_5yr$group == "Squid_BT", "BT Squid and Shrimp", std_ano_bio_5yr$biggrp)
std_ano_bio_5yr$biggrp <- ifelse(std_ano_bio_5yr$group == "euphausiid" | std_ano_bio_5yr$group == "fish" | std_ano_bio_5yr$group == "lg_zoop" | std_ano_bio_5yr$group == "shrimp" | std_ano_bio_5yr$group == "sm_zoop", "Zoops and Shrimps", std_ano_bio_5yr$biggrp)

std_ano_bio_5yr %>% 
  #filter(group == "Zoops and Shrimps") %>%
  ggplot( aes(x=YEAR, y=group, group=group, height = prop,  fill=YEAR)) + 
  #geom_point() + 
  #geom_line() +
  #geom_ma(ma_fun = SMA, n=5) +
  #geom_smooth(method = "loess", se=F, span = 1) + 
  geom_ridgeline_gradient(stat="identity", scale=1) +
  scale_fill_viridis_d(direction = -1, guide = "none") +
  facet_wrap(.~factor(stratcat, levels=c("GOMstrat", "GBstrat", "MABstrat", "SNEstrat")), scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        text = element_text(size=16))

#Assign which prey groupings belong to which predators - top 3 or 4
t$bp_prey <- ifelse(t$group == "Shrimp_BT" | t$group == "Clupeid" |t$group == "Sandlance", 1, 0) # also zoop
t$mn_prey <- ifelse(t$group == "Shrimp_BT" | t$group == "Misc_Fish" |t$group == "Sandlance"| t$group == "Clupeid", 1, 0)
t$ba_prey <- ifelse(t$group == "Shrimp_BT" | t$group == "Clupeid" |t$group == "Sandlance" | t$group == "Large_Gadid", 1, 0)
t$bb_prey <- ifelse(t$group == "Shrimp_BT" | t$group == "Clupeid"| t_region$group == "Squid_BT", 1, 0) #also zoop

t <- merge(t, SVSPP_table, by.x = "SVSPP", by.y = "survdat.SVSPP", all.x=T)

#Do the same as above but for the by_region dataset

t_region$bp_prey <- ifelse(t_region$group == "Shrimp_BT" | t_region$group == "Clupeid" |t_region$group == "Sandlance" , 1, 0) # also zoop
t_region$mn_prey <- ifelse(t_region$group == "Shrimp_BT" | t_region$group == "Misc_Fish"| t_region$group == "Clupeid" |t_region$group == "Sandlance", 1, 0)
t_region$ba_prey <- ifelse(t_region$group == "Shrimp_BT" | t_region$group == "Clupeid" |t_region$group == "Sandlance" | t_region$group == "Large_Gadid", 1, 0)
t_region$bb_prey <- ifelse(t_region$group == "Shrimp_BT" | t_region$group == "Clupeid" | t_region$group == "Squid_BT", 1, 0) #also zoop

t_region$bp_prey <- ifelse(t_region$group == "Shrimp_BT" | t_region$group == "Clupeid" |t_region$group == "Sandlance" , "Preferred Prey", "Other Prey") # also zoop
t_region$mn_prey <- ifelse(t_region$group == "Shrimp_BT" | t_region$group == "Misc_Fish"| t_region$group == "Clupeid" |t_region$group == "Sandlance", "Preferred Prey", "Other Prey")
t_region$ba_prey <- ifelse(t_region$group == "Shrimp_BT" | t_region$group == "Clupeid" |t_region$group == "Sandlance" | t_region$group == "Large_Gadid", "Preferred Prey", "Other Prey")
t_region$bb_prey <- ifelse(t_region$group == "Shrimp_BT" | t_region$group == "Clupeid" | t_region$group == "Squid_BT", "Preferred Prey", "Other Prey") #also zoop

t_region$sig <- ifelse(t_region$group == "Flatfish" | t_region$group == "Large_Gadid" | t_region$group == "Misc_Fish" | t_region$group == "shrimp" | t_region$group == "Small_Gadid", 1, 0)

#Merge back with the table of SVSPP numbers and common names to plot by common name
t_region <- merge(t_region, SVSPP_table, by.x = "SVSPP", by.y = "survdat.SVSPP", all.x=T)

# add sst

t <- merge(t, rast.avg_entire, by.x=c("YEAR"), by.y=c("YEAR"), all.x=T)

##### Boop over to the scaled biomass calculation to scale the biomass by proportion of diet for each species

t$bp_prey <- ifelse(t$group == "Shrimp_BT" | t$group == "Clupeid" |t$group == "Sandlance" , 1, 0) # also zoop
t$mn_prey <- ifelse(t$group == "Shrimp_BT" | t$group == "Misc_Fish"| t$group == "Clupeid" |t$group == "Sandlance", 1, 0)
t$ba_prey <- ifelse(t$group == "Shrimp_BT" | t$group == "Clupeid" |t$group == "Sandlance" | t$group == "Large_Gadid", 1, 0)
t$bb_prey <- ifelse(t$group == "Shrimp_BT" | t$group == "Clupeid" | t$group == "Squid_BT", 1, 0) #also zoop

t$bp_prey <- ifelse(t$group == "Shrimp_BT" | t$group == "Clupeid" |t$group == "Sandlance", "Preferred Prey", "Other Prey") # also zoop
t$mn_prey <- ifelse(t$group == "Shrimp_BT" | t$group == "Misc_Fish"| t$group == "Clupeid" |t$group == "Sandlance", 1, 0)
t$ba_prey <- ifelse(t$group == "Shrimp_BT" | t$group == "Clupeid" |t$group == "Sandlance" | t$group == "Large_Gadid", 1, 0)
t$bb_prey <- ifelse(t$group == "Shrimp_BT" | t$group == "Clupeid" | t$group == "Squid_BT", 1, 0) #also zoop

t_region$group2 <- t_region$group
t_region$group2 <- ifelse(t_region$group2 == "Large_Gadid", "Large Gadid", t_region$group2)
t_region$group2 <- ifelse(t_region$group2 == "Meso_pelagic", "Mesopelagic Fish", t_region$group2)
t_region$group2 <- ifelse(t_region$group2 == "Misc_Fish", "Miscellaneous Fish", t_region$group2)
t_region$group2 <- ifelse(t_region$group2 == "Northern_shrimp", "Northern Shrimp", t_region$group2)
t_region$group2 <- ifelse(t_region$group2 == "Shrimp_BT", "Shrimp", t_region$group2)
t_region$group2 <- ifelse(t_region$group2 == "Small_Gadid", "Small Gadid", t_region$group2)
t_region$group2 <- ifelse(t_region$group2 == "Squid_BT", "Squid", t_region$group2)


t_region2 <- t_region
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

#  ggplot for stacked bars
coll <- c("#e3657b", "#ff8b63", "#ffb165", "#eed948", "#74a76a", "#3a613c", "#6bd3aa", "#4b8497", "#bdcaf6", "#9a68af")

t %>% filter(!is.na(group)) %>%
  group_by(YEAR) %>%
  mutate(freq = strat.biomass / sum(strat.biomass)) %>% 
  ggplot(aes(fill=factor(group), y=strat.biomass, x=YEAR)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  #scale_fill_manual(values=coll)
  ggtitle("Bb")

setwd("/Users/jstepanuk/Dropbox/my_papers/InPrep_Stepanuketal_GOMprey")
pdf("bts_biomass_groups_scaledtotal_20220723.pdf", width=8, height=6)
t_region %>% filter(!is.na(group)) %>%
  group_by(stratcat,YEAR) %>%
  mutate(freq = strat.biomass / sum(strat.biomass)) %>% 
  ggplot(aes(fill=factor(group2), y=strat.biomass, x=YEAR)) + 
  geom_bar(position="stack", stat="identity") +
  facet_wrap(.~factor(stratcat, levels=c("GOMstrat", "GBstrat", "SNEstrat", "MABstrat"), labels = c("Gulf of Maine", "George's Bank", "Southern New England", "Mid-Atlantic Bight"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        text=element_text(size=16)) +
  scale_fill_manual(values=coll) +
  #scale_fill_manual(values=c("goldenrod", "darkgreen")) +
  ggtitle("Prey Group Relative Proportion - BTS") +
  labs(x="Year", y="Biomass", fill = "")
dev.off()

setwd("/Users/jstepanuk/Dropbox/my_papers/InPrep_Stepanuketal_GOMprey")
pdf("bts_prey_grouped_20220723.pdf", width=8, height=6)
ggplot(t_region2, aes(x=YEAR, y=value, group=variable, fill=variable)) + 
  geom_bar(stat="identity", position="dodge", color="black", size=0.3) + 
  facet_wrap(.~factor(stratcat, levels=c("GOMstrat", "GBstrat", "SNEstrat", "MABstrat"), labels = c("Gulf of Maine", "George's Bank", "Southern New England", "Mid-Atlantic Bight"))) +
  labs(x="Year", y="Biomass", fill = "") +
  scale_fill_manual(values=c("#E69F00", "#009E73", "#0072B2", "#F0E442")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        text=element_text(size=16)) 
dev.off()

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
zpd_sel <- ZPD_fall_sf %>% dplyr::select(group_multi, STRATUM, SEASON, STATION, CRUISE6, YEAR, BIOMASS, ABUNDANCE, CATCHSEX, stratcat, yr5)
zpd_sel <- st_drop_geometry(cbind(zpd_sel, st_coordinates(zpd_sel)))
colnames(zpd_sel)[12:13] <- c("LON", "LAT")

zpd_clusteryear <- zpd_sel[zpd_sel$YEAR > 1980,]
zpd_clusteryear$YEAR <- zpd_clusteryear$yr5

zpd_strat_mean <- calc_stratified_mean(data.table(zpd_clusteryear), areaDescription = "Name", groupDescription = "group_multi", filterBySeason = "FALL", areaPolygon = strat_zp)

zpd_strat_mean_gom <- calc_stratified_mean(data.table(zpd_clusteryear[zpd_clusteryear$stratcat == "GOMstrat",]), areaDescription = "Name", groupDescription = "group_multi", filterBySeason = "FALL", areaPolygon = strat_zp)

zpd_strat_mean_gb <- calc_stratified_mean(data.table(zpd_clusteryear[zpd_clusteryear$stratcat == "GBstrat",]), areaDescription = "Name", groupDescription = "group_multi", filterBySeason = "FALL", areaPolygon = strat_zp)

zpd_strat_mean_sne <- calc_stratified_mean(data.table(zpd_clusteryear[zpd_clusteryear$stratcat == "SNEstrat",]), areaDescription = "Name", groupDescription = "group_multi", filterBySeason = "FALL", areaPolygon = strat_zp)

zpd_strat_mean_mab <- calc_stratified_mean(data.table(zpd_clusteryear[zpd_clusteryear$stratcat == "MABstrat",]), areaDescription = "Name", groupDescription = "group_multi", filterBySeason = "FALL", areaPolygon = strat_zp)

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

zpd_strat_mean$mnba_prey <- ifelse(zpd_strat_mean$group_multi == "shrimp" | zpd_strat_mean$group_multi == "euphausiid",  "Preferred Prey", "Other Prey") # also zoop
zpd_strat_mean$bbbp_prey <- ifelse(zpd_strat_mean$group_multi == "shrimp" | zpd_strat_mean$group_multi == "euphausiid" | zpd_strat_mean$group_multi == "sm_zoop"| zpd_strat_mean$group_multi == "lg_zoop" , "Preferred Prey", "Other Prey")

zpd_region$mnba_prey <- ifelse(zpd_region$group_multi == "shrimp" | zpd_region$group_multi == "euphausiid",  "Preferred Prey", "Other Prey") # also zoop
zpd_region$bbbp_prey <- ifelse(zpd_region$group_multi == "shrimp" | zpd_region$group_multi == "euphausiid" | zpd_region$group_multi == "sm_zoop"| zpd_region$group_multi == "lg_zoop" , "Preferred Prey", "Other Prey")

zpd_region$group_multi2 <- zpd_region$group_multi
zpd_region$group_multi2 <- ifelse(zpd_region$group_multi2 == "euphausiid", "Euphausiid", zpd_region$group_multi2)
zpd_region$group_multi2 <- ifelse(zpd_region$group_multi2 == "lg_zoop", "Mysids & Amphipods", zpd_region$group_multi2)
zpd_region$group_multi2 <- ifelse(zpd_region$group_multi2 == "sm_zoop", "Copepods", zpd_region$group_multi2)
zpd_region$group_multi2 <- ifelse(zpd_region$group_multi2 == "shrimp", "Small Shrimp", zpd_region$group_multi2)

zpd_region2 <- zpd_region
zpd_region2$fin_sei_prey <- ifelse(zpd_region2$bbbp_prey == "Preferred Prey", zpd_region2$strat.biomass, NA)
zpd_region2$humpback_minke_prey <- ifelse(zpd_region2$mnba_prey == "Preferred Prey", zpd_region2$strat.biomass, NA)
p <- aggregate(fin_sei_prey ~ YEAR + stratcat, data=zpd_region2, FUN=sum)
s <- aggregate(humpback_minke_prey ~ YEAR + stratcat, data=zpd_region2, FUN=sum)

zpd_region2 <- merge(p, s, by.x=c("YEAR", "stratcat"), by.y=c("YEAR", "stratcat"))
zpd_region2 <- melt(zpd_region2, id.vars=c("YEAR", "stratcat"))


zpd_strat_mean %>% filter(group_multi!=c("")) %>%
  group_by(YEAR) %>%
  mutate(freq = strat.biomass / sum(strat.biomass)) %>% 
  ggplot(aes(fill=factor(group_multi), y=freq, x=YEAR)) + 
  geom_bar(position="stack", stat="identity")

setwd("/Users/jstepanuk/Dropbox/my_papers/InPrep_Stepanuketal_GOMprey")
pdf("zoop_biomass_groups_freq_20220723.pdf", width=8, height=6)
zpd_region %>% filter(group_multi!=c("")) %>%
  group_by(stratcat, YEAR) %>%
  mutate(freq = strat.biomass / sum(strat.biomass)) %>% 
  ggplot(aes(fill=factor(group_multi2), y=freq, x=YEAR)) + 
  facet_wrap(.~factor(stratcat, levels=c("GOMstrat", "GBstrat", "SNEstrat", "MABstrat"), labels = c("Gulf of Maine", "George's Bank", "Southern New England", "Mid-Atlantic Bight"))) +
  geom_bar(position="stack", stat="identity") +
  #scale_fill_manual(values=c("goldenrod", "darkgreen")) +
  labs(x="Year", y="Abundance", fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        text=element_text(size=16)) 
dev.off()

setwd("/Users/jstepanuk/Dropbox/my_papers/InPrep_Stepanuketal_GOMprey")
pdf("zoop_prey_grouped_20220727_mn.pdf", width=7, height=6)
ggplot(zpd_region2[zpd_region2$variable == "humpback_minke_prey",], aes(x=YEAR, y=value, group=variable, fill=variable)) + 
  geom_bar(stat="identity", position="dodge", color="black", size=0.3) + 
  facet_wrap(.~factor(stratcat, levels=c("GOMstrat", "GBstrat", "SNEstrat", "MABstrat"), labels = c("Gulf of Maine", "George's Bank", "Southern New England", "Mid-Atlantic Bight"))) +
  labs(x="Year", y="Abundance", fill = "") +
  scale_y_continuous(limits = c(0, 4e06), labels = function(x) format(x, scientific = TRUE)) +
  theme_minimal() +
  scale_fill_manual(values=c( "#008499")) + #"#ECC828",
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        text=element_text(size=16),
        legend.position = "bottom") 
dev.off()

sss <- 0.0005 #relational of axes of biomass and abundance
# ggplot(data=group_allprey_dch[group_allprey_dch$year >= 1977 & group_allprey_dch$year <= 2019 ,], aes(x=year, y=wt_mean/1000, color=t.rast.avg_entire., group=as.character(group))) + 
#   geom_point(data=group_allprey_dch, aes(y =wt_mean / 1000), color="darkgreen", size=0.8) +
#   geom_smooth(data=group_allprey_dch, aes(y=wt_mean/1000), method="lm", se=T, color="darkgreen", fill="darkgreen") +
#   
#   #Below plots any category with an abundance - aka ecomon groups
#   #geom_point(data=dual_comb, aes(y=zoop_abundance/sss), alpha=0.4, pch=16) +
#   geom_smooth(data=group_allprey_dch, aes(y=t.rast.avg_entire./sss), method="lm", se=T, color="blue", fill= "blue") +
#   
ggplot(zpd_region2, aes(x=YEAR, group=variable, fill=variable)) + 
  geom_bar(data=zpd_region2[zpd_region2$variable == "fin_sei_prey",], aes(y= value), stat="identity", position="dodge", color="black", size=0.3) + 
  geom_bar(data=zpd_region2[zpd_region2$variable == "humpback_minke_prey",], aes(y= value/sss), stat="identity", position="dodge", color="black", size=0.3) + 
  facet_wrap(.~factor(stratcat, levels=c("GOMstrat", "GBstrat", "SNEstrat", "MABstrat"), labels = c("Gulf of Maine", "George's Bank", "Southern New England", "Mid-Atlantic Bight"))) +
  labs(x="Year",fill = "") +
  theme_minimal() +
  scale_fill_manual(values=c("#ECC828", "#008499")) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        text=element_text(size=16)) +
  scale_y_continuous(
    # Features of the first axis
    name = "Fin and Sei whale prey abundance",
    # Add a second axis and specify its features
    sec.axis = sec_axis( trans=~.*sss, name="humpback and minke whale prey abundance")
  ) 


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
    ss <- sd(www$strat.biomass)
    www$std_ano <- (www$strat.biomass - mm) / ss
    std_ano_bio <- rbind(std_ano_bio, www)
  }
}

std_ano_bio$biggrp <- ifelse(std_ano_bio$group == "Clupeid" | std_ano_bio$group == "Flatfish" | std_ano_bio$group == "Meso_pelagic" | std_ano_bio$group == "Misc_Fish" | std_ano_bio$group == "Sandlance" | std_ano_bio$group == "Scombrid" | std_ano_bio$group == "Small_Gadid", "BT Fish", NA)
std_ano_bio$biggrp <- ifelse(std_ano_bio$group == "Northern_shrimp" | std_ano_bio$group == "Shrimp_BT" | std_ano_bio$group == "Squid_BT", "BT Squid and Shrimp", std_ano_bio$biggrp)
std_ano_bio$biggrp <- ifelse(std_ano_bio$group == "euphausiid" | std_ano_bio$group == "fish" | std_ano_bio$group == "lg_zoop" | std_ano_bio$group == "shrimp" | std_ano_bio$group == "sm_zoop", "Zoops and Shrimps", std_ano_bio$biggrp)

std_ano_bio %>% 
  #filter(group == "Zoops and Shrimps") %>%
  ggplot( aes(x=YEAR, y=std_ano, group=group, color=group)) + geom_point() + geom_smooth(method = "lm", se=F) + facet_wrap(.~stratcat, scales = "free_y")

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

