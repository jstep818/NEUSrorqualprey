library(vegan)
library(ggnewscale)
#use t_region from bio_analysis_all_20211213.R




#### Old stsuff
# 
# t_region_species_scores$group <- ifelse(t_region_species_scores$SVSPP %in% c(32, 851, 205, 859, 427, 430, 30, 31, 44),  "Clupeid", NA)
# t_region_species_scores$group <- ifelse(t_region_species_scores$SVSPP %in% c(121, 578, 572, 209, 124, 898, 570, 702, 571, 582, 744, 860, 208, 212, 211, 874, 822, 745, 577, 877), "Scombrid", t_region_species_scores$group)
# t_region_species_scores$group <- ifelse(t_region_species_scores$SVSPP %in% c(181, 734, 38), "Sandlance", t_region_species_scores$group)
# t_region_species_scores$group <- ifelse((t_region_species_scores$SVSPP >= 501 & t_region_species_scores$SVSPP <= 505) | (t_region_species_scores$SVSPP >= 510 & t_region_species_scores$SVSPP <= 513), "Squid_BT", t_region_species_scores$group)
# t_region_species_scores$group <- ifelse(t_region_species_scores$SVSPP == 297 | t_region_species_scores$SVSPP == 298 | t_region_species_scores$SVSPP == 306 | t_region_species_scores$SVSPP == 307, "Shrimp_BT", t_region_species_scores$group)
# t_region_species_scores$group <- ifelse(t_region_species_scores$SVSPP %in% c(73:75), "Large_Gadid", t_region_species_scores$group)
# t_region_species_scores$group <- ifelse(t_region_species_scores$SVSPP %in% c(189, 1, 262, 192, 131, 750, 176, 168, 249, 193, 155, 164), "Misc_Fish", t_region_species_scores$group) 
# t_region_species_scores$group <- ifelse(t_region_species_scores$SVSPP %in% c(81, 454, 80, 86, 79, 69, 77, 72, 455, 78, 76), "Small_Gadid", t_region_species_scores$group)
# t_region_species_scores$group <- ifelse(t_region_species_scores$SVSPP %in% c(790, 792, 110, 793, 777, 100, 104, 785, 786, 788, 109, 795, 775, 776, 773, 791, 787, 117, 789, 783, 103, 853, 774, 873, 106, 107, 105), "Flatfish", t_region_species_scores$group)
# t_region_species_scores$group <- ifelse(t_region_species_scores$SVSPP %in% c(51, 54, 56, 55, 229) , "Meso_pelagic", t_region_species_scores$group)
# 
# t_region_species_scores$bp_prey <- ifelse(t_region_species_scores$group == "Shrimp_BT" | t_region_species_scores$group == "Clupeid" |t_region_species_scores$group == "Sandlance", 1, 0) # also zoop
# t_region_species_scores$mn_prey <- ifelse(t_region_species_scores$group == "Shrimp_BT" | t_region_species_scores$group == "Misc_Fish" |t_region_species_scores$group == "Sandlance", 1, 0)
# t_region_species_scores$ba_prey <- ifelse(t_region_species_scores$group == "Shrimp_BT" | t_region_species_scores$group == "Clupeid" |t_region_species_scores$group == "Sandlance" | t_region_species_scores$group == "Large_Gadid", 1, 0)
# t_region_species_scores$bb_prey <- ifelse(t_region_species_scores$group == "Shrimp_BT" | t_region_species_scores$group == "Clupeid", 1, 0) #also zoop
# 
# #Hull data - https://chrischizinski.github.io/rstats/vegan-ggplot2/
# t_region_species_scores <- t_region_species_scores[t_region_species_scores$spp != "NA ",]
# 
# group.misc <- t_region_species_scores[t_region_species_scores$group == "Misc_Fish", ][chull(t_region_species_scores[t_region_species_scores$group == "Misc_Fish", c("NMDS1", "NMDS2")]), ]  
# 
# group.shr <- t_region_species_scores[t_region_species_scores$group == "Shrimp_BT", ][chull(t_region_species_scores[t_region_species_scores$group == "Shrimp_BT", c("NMDS1", "NMDS2")]), ]  
# 
# group.scom <- t_region_species_scores[t_region_species_scores$group == "Scombrid" & t_region_species_scores$NMDS1 != "NaN", ][chull(t_region_species_scores[t_region_species_scores$group == "Scombrid"& t_region_species_scores$NMDS1 != "NaN", c("NMDS1", "NMDS2")]), ]  
# 
# group.sl <- t_region_species_scores[t_region_species_scores$group == "Sandlance" & t_region_species_scores$NMDS1 != "NaN", ][chull(t_region_species_scores[t_region_species_scores$group == "Sandlance" & t_region_species_scores$NMDS1 != "NaN", c("NMDS1", "NMDS2")]), ]  
# 
# group.clup <- t_region_species_scores[t_region_species_scores$group == "Clupeid", ][chull(t_region_species_scores[t_region_species_scores$group == "Clupeid", c("NMDS1", "NMDS2")]), ]  
# 
# group.squid <- t_region_species_scores[t_region_species_scores$group == "Squid_BT", ][chull(t_region_species_scores[t_region_species_scores$group == "Squid_BT", c("NMDS1", "NMDS2")]), ]  
# 
# group.lggad <- t_region_species_scores[t_region_species_scores$group == "Large_Gadid", ][chull(t_region_species_scores[t_region_species_scores$group == "Large_Gadid", c("NMDS1", "NMDS2")]), ]  
# 
# hull.data <- rbind(group.misc, group.shr)  #combine group.a and group.b
# hull.data <- rbind(hull.data, group.scom)
# hull.data <- rbind(hull.data, group.sl)
# hull.data <- rbind(hull.data, group.clup)
# hull.data <- rbind(hull.data, group.squid)
# hull.data <- rbind(hull.data, group.lggad)



### SINGLE NMDS
#t_region_wide <- reshape(t_region[!is.na(t_region$group),], idvar = "YEAR", timevar = "survdat.COMNAME", v.names="strat.biomass", direction = "wide")
t_region_wide <- reshape(stratbio_allregion[!is.na(stratbio_allregion$group) & (stratbio_allregion$YEAR >= 1980 & stratbio_allregion$YEAR <= 2019),], idvar = "YEAR", timevar = "group", v.names="strat.biomass", direction = "wide")
colnames(t_region_wide) <- gsub("strat.biomass", "", colnames(t_region_wide))
rownames(t_region_wide) <- paste0("y",as.character(t_region_wide$YEAR))

t_region_wide_sub <- t_region_wide %>%
  #dplyr::select( 17:72)
  #dplyr::select( 16:527)
  #dplyr::select( 16:22)
  dplyr::select(4:17)
#dplyr::select( 16:52)
t_region_wide_sub[is.na(t_region_wide_sub)] <- 0
t_region_wide_sub <- t_region_wide_sub %>%
  filter_all(any_vars(. != 0))
t_region_wide_sub <- data.frame(t_region_wide_sub, row.names = as.character(t_region_wide$YEAR))

t_region_mds <- metaMDS(t_region_wide_sub, distance = "bray", k=2, trymax = 200)

t_region_scores <-as.data.frame(scores(t_region_mds, display = "sites"))
t_region_scores$year <- rownames(t_region_scores)
t_region_species_scores <- as.data.frame(scores(t_region_mds, "species"))
t_region_species_scores$spp <- gsub("\\.", " ", rownames(t_region_species_scores))
t_region_species_scores$spp <- gsub("SHRIMP  PINK  BROWN   WHITE ", "SHRIMP (PINK -BROWN - WHITE)", t_region_species_scores$spp)
t_region_species_scores <- merge(t_region_species_scores, SVSPP_table, by.x="spp", by.y = "survdat.COMNAME", all.x=T)
colnames(t_region_species_scores)[4] <- "SVSPP"


#### Add in SST
yrs_match <- data.frame(t_region_wide$YEAR, rep(1, nrow(t_region_wide)))
colnames(yrs_match) <- c("YEAR", "reps")
yrs_match <- merge(rast.avg_entire, yrs_match, by.x="YEAR", by.y="YEAR", all.x=T)
yrs_match$rnms <- rownames(yrs_match)
yrs_id <- yrs_match$rnms[!is.na(yrs_match$rnms[yrs_match$reps == 1])]
col_id <- match(nn, colnames(rast.avg_entire_reg))

#env.fit <- envfit(t_region_mds, rast.avg_entire_reg[yrs_id, c(col_id, 6)], perm = 999)
env.fit <- envfit(t_region_mds, rast.avg_entire[yrs_id, c(3,2)], perm = 999)
#env.fit 
spp.scrs <- as.data.frame(scores(env.fit, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))

a <- ggplot() +
  geom_point(data = t_region_scores, aes(x = NMDS1, y = NMDS2, color = as.numeric(year)), size = 3) +
  scale_color_viridis_c(name = "Year") +
  geom_text(data = t_region_species_scores, aes(x = NMDS1, y = NMDS2, label = spp),
            alpha = 0.8, size = 4) +
  annotate(geom = "label", x = 0, y = 1.25, size = 6,
           label = paste(nn, "- Stress: ", round(t_region_mds$stress, digits = 3))) +
  theme_minimal() +
  theme(legend.position = "right",
        text = element_text(size = 14)) +
  geom_abline(intercept=0, slope=0) +
  geom_vline(xintercept = 0) +
  ggnewscale::new_scale("color") +
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2, color=Species),
               arrow = arrow(length = unit(0.5, "cm")), size=1.2) +
  scale_color_manual(name="Env Variables", labels = c("SST", "Year"), values=c("Orange", "Purple"))
#geom_segment(data=t_region_species_scores, aes(xend=NMDS1, yend=NMDS2, group=group,
#x=0, y=0), color=NA)
#geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,group=group),alpha=0, color=NA) 
nmds_plots[[i]] <- a
nmds_env[[i]] <- env.fit
nmds_nmds[[i]] <- t_region_mds
}
cowplot::plot_grid(nmds_plots[[1]], nmds_plots[[2]], nmds_plots[[3]], nmds_plots[[4]], ncol=2)
