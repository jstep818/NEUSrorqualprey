library(rgdal)
library(sp)
library(graticule)
library(rgeos)
library(raster)
library(ggplot2)
library(reshape2)
library(lubridate)
library(RColorBrewer)
library(sf)
library(dplyr)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(lwgeom)
library(mapview)
library(tmap)

setwd("~/Documents/Github/cet_prey/calcs_envr/")
euc_dist_along_shelf <- raster("euc_dist_along_shelf.tif")

euc_dist_along_shelf_sub <- reclassify(euc_dist_along_shelf, cbind(2000000, 7000000, NA), right = TRUE)
euc_dist_along_shelf_sub[is.infinite(euc_dist_along_shelf_sub)]  <- NA

natl_coastline <- ne_countries(scale=10, returnclass = "sf")
natl_coastline <- st_crop(natl_coastline, c(xmin=-81, xmax=-55, ymin=30, ymax=50))

natl_coastline <- natl_coastline %>%
  st_transform(crs = "+proj=aea +lat_0=34 +lon_0=-78 +lat_1=27.3333333333333 +lat_2=40.6666666666667 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>%
  group_by( featurecla) %>%  
  summarise(Project = unique(featurecla), do_union = TRUE)

euc_dist_along_shelf_sub_con1 <- st_centroid(st_as_sf(rasterToContour(euc_dist_along_shelf_sub, levels= group_allprey_dch$wt_mean[group_allprey_dch$group == "Meso_pelagic" & group_allprey_dch$year >= 1980])))
euc_dist_along_shelf_sub_con1$yr <- 1980:2019

euc_dist_along_shelf_sub_con2 <- (st_as_sf(rasterToContour(euc_dist_along_shelf_sub, levels= group_allprey_dch$wt_mean[group_allprey_dch$group == "LONGFIN SQUID"])))
euc_dist_along_shelf_sub_con2$yr <- 1977:2018

strat_all_simple <- strat_all_sf %>%
  group_by(region) %>%
  summarize(geometry = st_union(geometry)) 

setwd( "/Users/jstepanuk/Dropbox/my_papers/InPrep_Stepanuketal_GOMprey")
pdf("latshift_map_20220727.pdf", width=6, height = 6)
ggplot(data=natl_coastline) +
  #layer_spatial(euc_dist_along_shelf_sub / 1000) +
  #
  #scale_fill_viridis_c(na.value=NA) +
  geom_sf(data=strat_all_simple, fill = "gray90", color="black", size=1) +
  geom_sf(data=st_as_sf(raster::rasterToContour(euc_dist_along_shelf_sub / 1000, levels = seq(0, 1300, 100))), aes(color=factor(level, levels=rev(seq(0,1300,100)))), size=2, lineend="round")+
  scale_color_viridis_d(option = "C") +
  geom_sf(data=strat_all_simple, fill = "NA", color="black", size=1) +
  geom_sf(data=natl_coastline, color=NA, fill="gray20") +
  # geom_sf(data=euc_dist_along_shelf_sub_con1[euc_dist_along_shelf_sub_con1$yr >= 1990,], aes(color=as.numeric(yr))) +
  # geom_sf(data=euc_dist_along_shelf_sub_con2[euc_dist_along_shelf_sub_con2$yr >= 1990,], aes(color=as.numeric(yr))) +
  #geom_sf(data=extentt_orth, fill=NA, color="black", size=1.5)+
  #annotate(geom="rect", xmin=-78, xmax=-63, ymin=34, ymax=46, color="black", fill="green", alpha=0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=0)) +
  labs(color="Distance from Cape Hatteras (km)") +
  coord_sf(crs="+proj=aea +lat_0=34 +lon_0=-78 +lat_1=27.3333333333333 +lat_2=40.6666666666667 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs", xlim=c(-33343.4,1036657), ylim=c(130000,1354041))
dev.off()

group_allprey_dch_sub <- merge(group_allprey_dch, unique(data.frame(t_region$group)), by.x="group_multi", by.y="t_region.group", all.x=T, all.y=F)

group_allprey_dch_sub_yr <- group_allprey_dch_sub[(group_allprey_dch_sub$year == 1981) & !is.na(group_allprey_dch_sub$group_multi),]

#zzmelt_slope
group_trends <- merge(group_allprey_dch_sub_yr, corrtest_biolat, by.x=c("group_multi"), by.y=c("group"))
group_trends$wt_mean[12] <- 251585.73
group_trends$wt_meankm <- group_trends$wt_mean / 1000
group_trends$slopescaled <- group_trends$slope * 39 
group_trends$endloc <- group_trends$wt_meankm + group_trends$slopescaled

#group_trends <- merge(group_trends, group_allprey_dch_sub_yr, by.x=c("X1"), by.y=c("group_multi"))
#group_trends <- group_trends[!is.na(group_trends$wt_mean),]

ggplot(group_trends, aes(x=group_multi, y=wt_meankm)) +
  geom_point(data=group_trends, aes(x=group_multi, y=wt_meankm), color="black", pch=95, size=5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90),
        strip.background = element_rect(fill=NA,colour="grey50"),
        panel.border = element_rect(fill=NA, colour="grey50"),
        strip.text.x = element_text(angle=90)
  ) +
  #facet_grid(~t_region.group, space="free_x", scales="free_x",) +
  geom_segment(data=group_trends[group_trends$sig == 1 & group_trends$year == 1980,], #& group_trends$slopegrp == "pole_year_slope"
               aes(y = wt_meankm, yend = endloc, 
                   x=group_multi, xend=group_multi, color=slope),
               arrow = arrow(length = unit(0.2, "cm")), size=2) +
  scale_color_viridis_c()

setwd( "/Users/jstepanuk/Dropbox/my_papers/InPrep_Stepanuketal_GOMprey")
pdf("latshift_plot_20220727.pdf", width=6, height = 5)
ggplot(group_allprey_dch, aes(x=group_multi, y=wt_mean / 1000)) +
  geom_point(data=group_allprey_dch, aes(x=group_multi, y=wt_mean / 1000, color=year), pch=95, size=8) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90),
        strip.background = element_rect(fill=NA,colour="grey50"),
        panel.border = element_rect(fill=NA, colour="grey50"),
        strip.text.x = element_text(angle=90)
  ) +
  #facet_grid(~t_region.group, space="free_x", scales="free_x",) +
   geom_segment(data=group_trends[group_trends$sig == 1 ,], 
                aes(y = wt_meankm, yend = endloc, 
                    x=group_multi, xend=group_multi), color="black",
                arrow = arrow(length = unit(0.2, "cm")), size=2) +
  scale_color_viridis_c() +
  labs(color="Year", x="", y="Distance from Cape Hatteras (km)") +
  scale_y_continuous(breaks=seq(0,1300,100))
dev.off()
