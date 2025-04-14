#calculate correlations


############################ RUNNING MODELS -----
#Calculate dch per 1) year, 2) sst, 3) BT for each species
corrtest_biolat <- NULL
for(j in unique(group_allprey_dch$group_multi)[!is.na(unique(group_allprey_dch$group_multi))]){
  o <- group_allprey_dch[group_allprey_dch$group_multi == j & group_allprey_dch$year >= 1980 & group_allprey_dch$year <=2019,]
  if(nrow(o) > 2 & length(o$wt_mean[!is.na(o$wt_mean)]) > 2){
    r1 <- lm(o$wt_mean/1000 ~ o$year)
    r2 <- lm(o$wt_mean/1000 ~ o$t.rast.avg_entire.)
    #r3 <- lm(o$wt_mean/1000 ~ o$bot_temp)
    p1 <- c("YEAR", j, summary(r1)$coefficients[1], 
            summary(r1)$coefficients[2], 
            summary(r1)$coefficients[8], 
            summary(r1)$coefficients[4],  
            summary(r1)$coefficients[6], 
            summary(r1)$adj.r.squared)
    # p2 <- c("SST", j, summary(r2)$coefficients[1], 
    #         summary(r2)$coefficients[2], 
    #         summary(r2)$coefficients[8], 
    #         summary(r2)$coefficients[4], 
    #         summary(r2)$coefficients[6], 
    #         summary(r2)$adj.r.squared)
    # p3 <- c("BT", j, summary(r3)$coefficients[1], 
    #         summary(r3)$coefficients[2], 
    #         summary(r3)$coefficients[8], 
    #         summary(r3)$coefficients[4], 
    #         summary(r3)$coefficients[6], 
    #         summary(r3)$adj.r.squared)

    p4 <- rbind(p1, p2)
    #p4 <- rbind(p4, p3)
    corrtest_biolat <- rbind(corrtest_biolat, p1)
  }
}
#}
corrtest_biolat <- data.frame(corrtest_biolat)
colnames(corrtest_biolat) <- c("type", "group", "est", "slope",  "pval", "std_err_slope", "tval", "adj_r2")
corrtest_biolat$pval <- as.numeric(corrtest_biolat$pval)
corrtest_biolat$slope <- as.numeric(corrtest_biolat$slope)
corrtest_biolat$adj_r2 <- as.numeric(corrtest_biolat$adj_r2)
corrtest_biolat$sig <- ifelse(corrtest_biolat$pval <= 0.003, 1, 0) # & corrtest_biolat$slope >= 0
corrtest_biolat$sigpos <- ifelse(corrtest_biolat$pval <= 0.003 & corrtest_biolat$slope >= 0, 1, 0) 
corrtest_biolat$signeg <- ifelse(corrtest_biolat$pval <= 0.003 & corrtest_biolat$slope < 0, 1, 0) 


summary(lm(wt_mean ~ year*group_multi, data=group_allprey_dch))
emtrends(lm(wt_mean ~ year*group_multi, data=group_allprey_dch), specs = "group_multi", var="year", )

library(lme4)
summary(lmer(wt_mean ~ year + (1+year|group_multi), data=group_allprey_dch))
car::Anova(lmer(wt_mean ~ year + (1+year|group_multi), data=group_allprey_dch))
plot(lmer(wt_mean ~ year + (1+year|group_multi), data=group_allprey_dch))


anova(lmer(wt_mean ~ year + (year|group_multi), data=group_allprey_dch, REML=F), lmer(wt_mean ~ year + (1|group_multi), data=group_allprey_dch, REML=F))

## plot it out
group_allprey_dch$biggrp <- ifelse(group_allprey_dch$group_multi == "Clupeid"  | group_allprey_dch$group_multi == "Meso_pelagic" | group_allprey_dch$group_multi == "Misc_Fish" | group_allprey_dch$group_multi == "Sandlance" | group_allprey_dch$group_multi == "Scombrid" | group_allprey_dch$group_multi == "Small_Gadid", "BT Fish", NA)
group_allprey_dch$biggrp <- ifelse(group_allprey_dch$group_multi == "Northern_shrimp" | group_allprey_dch$group_multi == "Shrimp_BT" | group_allprey_dch$group_multi == "Squid_BT", "BT Squid and Shrimp", group_allprey_dch$biggrp)
group_allprey_dch$biggrp <- ifelse(group_allprey_dch$group_multi == "euphausiid" | group_allprey_dch$group_multi == "fish" | group_allprey_dch$group_multi == "lg_zoop" | group_allprey_dch$group_multi == "shrimp" | group_allprey_dch$group_multi == "sm_zoop", "Zoops and Shrimps", group_allprey_dch$biggrp)

ggplot(group_allprey_dch[group_allprey_dch$year >= 1980,], aes(y=wt_mean / 1000, x= year, group = group_multi, color=group_multi, fill=group_multi)) + geom_smooth(method="lm", se=T) + facet_wrap(.~biggrp) + geom_point()

# Calculate biomass per year per region vs. 1) year, 2) sst, 3) bt
std_ano_bio$stratcat2 <- ifelse(std_ano_bio$stratcat == "GOMstrat" | std_ano_bio$stratcat == "GBstrat", "North", "South")

corrtest_biomass <- NULL
for(i in unique(std_ano_bio$stratcat)){
  ww <- std_ano_bio[std_ano_bio$stratcat == i,]
  for(j in unique(ww$group)){
    o <- ww[ww$group == j & ww$YEAR >= 1980 & ww$YEAR <=2019,]
    if(nrow(o) > 2 & length(o$strat.biomass[!is.na(o$strat.biomass)]) > 2){
      r1 <- lm(o$strat.biomass ~ o$YEAR)
      #r2 <- lm(o$strat.biomass ~ o$SURFTEMP)
      #r3 <- lm(o$strat.biomass ~ o$bot_temp)
      p1 <- c("YEAR",i, j, summary(r1)$coefficients[1], 
              summary(r1)$coefficients[2], 
              summary(r1)$coefficients[8], 
              summary(r1)$coefficients[4],  
              summary(r1)$coefficients[6], 
              summary(r1)$adj.r.squared)
      # p2 <- c("SST",i, j, summary(r2)$coefficients[1], 
      #         summary(r2)$coefficients[2], 
      #         summary(r2)$coefficients[8], 
      #         summary(r2)$coefficients[4], 
      #         summary(r2)$coefficients[6], 
      #         summary(r2)$adj.r.squared)
      # p3 <- c("BT",i, j, summary(r3)$coefficients[1], 
      #         summary(r3)$coefficients[2], 
      #         summary(r3)$coefficients[8], 
      #         summary(r3)$coefficients[4], 
      #         summary(r3)$coefficients[6], 
      #         summary(r3)$adj.r.squared)
      # 
      # p4 <- rbind(p1, p2)
      # p4 <- rbind(p4, p3)
      corrtest_biomass <- rbind(corrtest_biomass, p1)
    }
  }
}
corrtest_biomass <- data.frame(corrtest_biomass)
colnames(corrtest_biomass) <- c("type","stratcat", "group", "est", "slope",  "pval", "std_err_slope", "tval", "adj_r2")
corrtest_biomass$pval <- as.numeric(corrtest_biomass$pval)
corrtest_biomass$slope <- as.numeric(corrtest_biomass$slope)
corrtest_biomass$adj_r2 <- as.numeric(corrtest_biomass$adj_r2)
corrtest_biomass$sig <- ifelse(corrtest_biomass$pval <= 0.05 , 1, 0) #& corrtest_biomass$slope >= 0
corrtest_biomass$sigpos <- ifelse(corrtest_biomass$pval <= 0.05 & corrtest_biomass$slope >= 0, 1, 0) 
corrtest_biomass$signeg <- ifelse(corrtest_biomass$pval <= 0.05 & corrtest_biomass$slope < 0, 1, 0)


#All biomass
corrtest_biomass2 <- NULL
  for(j in unique(ww$group)){
    o <- ww[ww$group == j & ww$YEAR >= 1980 & ww$YEAR <=2019,]
    if(nrow(o) > 2 & length(o$strat.biomass[!is.na(o$strat.biomass)]) > 2){
      r1 <- lm(o$strat.biomass ~ o$YEAR)
      #r2 <- lm(o$strat.biomass ~ o$SURFTEMP)
      #r3 <- lm(o$strat.biomass ~ o$bot_temp)
      p1 <- c("YEAR",i, j, summary(r1)$coefficients[1], 
              summary(r1)$coefficients[2], 
              summary(r1)$coefficients[8], 
              summary(r1)$coefficients[4],  
              summary(r1)$coefficients[6], 
              summary(r1)$adj.r.squared)
      # p2 <- c("SST",i, j, summary(r2)$coefficients[1], 
      #         summary(r2)$coefficients[2], 
      #         summary(r2)$coefficients[8], 
      #         summary(r2)$coefficients[4], 
      #         summary(r2)$coefficients[6], 
      #         summary(r2)$adj.r.squared)
      # p3 <- c("BT",i, j, summary(r3)$coefficients[1], 
      #         summary(r3)$coefficients[2], 
      #         summary(r3)$coefficients[8], 
      #         summary(r3)$coefficients[4], 
      #         summary(r3)$coefficients[6], 
      #         summary(r3)$adj.r.squared)
      # 
      # p4 <- rbind(p1, p2)
      # p4 <- rbind(p4, p3)
      corrtest_biomass2 <- rbind(corrtest_biomass2, p1)
    }
  }
###########

### calculate correlations


zzall <- merge(t_region[,-17], corrtest_biomass_sst_gom, by.x="survdat.COMNAME", by.y="group", all.x=T, all.y=T)
colnames(zzall)[26:29] <- c("bio_sst_gom_pval", "bio_sst_gom_slope", "bio_sst_gom_r2", "bio_sst_gom_sig")
zzall <- merge(zzall, corrtest_biomass_sst_gb, by.x="survdat.COMNAME", by.y="group", all.x=T, all.y=T)
colnames(zzall)[30:33] <- c("bio_sst_gb_pval", "bio_sst_gb_slope", "bio_sst_gb_r2", "bio_sst_gb_sig")
zzall <- merge(zzall, corrtest_biomass_sst_sne, by.x="survdat.COMNAME", by.y="group", all.x=T, all.y=T)
colnames(zzall)[34:37] <- c("bio_sst_sne_pval", "bio_sst_sne_slope", "bio_sst_sne_r2", "bio_sst_sne_sig")
zzall <- merge(zzall, corrtest_biomass_sst_MAB, by.x="survdat.COMNAME", by.y="group", all.x=T, all.y=T)
colnames(zzall)[38:41] <- c("bio_sst_mab_pval", "bio_sst_mab_slope", "bio_sst_mab_r2", "bio_sst_mab_sig")
zzall <- merge(zzall, corrtest_biolat_sst, by.x="survdat.COMNAME", by.y="group", all.x=T, all.y=T)
colnames(zzall)[42:45] <- c("pole_sst_pval", "pole_sst_slope", "pole_sst_r2", "pole_sst_sig")
zzall <- merge(zzall, corrtest_biolat_yr, by.x="survdat.COMNAME", by.y="group", all.x=T, all.y=T)
colnames(zzall)[46:49] <- c("pole_year_pval", "pole_year_slope", "pole_year_r2", "pole_year_sig")

zzall$bio_sst_slopenorm <- ((2*(zzall$bio_sst_slope - min(zzall$bio_sst_slope, na.rm=T))) / (max(zzall$bio_sst_slope, na.rm=T) - min(zzall$bio_sst_slope, na.rm=T))) - 1
zzall$pole_sst_slopenorm <- ((2*(zzall$pole_sst_slope - min(zzall$pole_sst_slope, na.rm=T))) / (max(zzall$pole_sst_slope, na.rm=T) - min(zzall$pole_sst_slope, na.rm=T))) - 1
zzall$pole_year_slopenorm <- ((2*(zzall$pole_year_slope - min(zzall$pole_year_slope, na.rm=T))) / (max(zzall$pole_year_slope, na.rm=T) - min(zzall$pole_year_slope, na.rm=T))) - 1

plotfor_regions <- c("bio_sst_gom_slope", "bio_sst_gb_slope", "bio_sst_sne_slope", "bio_sst_mab_slope")
plotfor_regions_lab <- c("bio_sst_gom_sig", "bio_sst_gb_sig", "bio_sst_sne_sig", "bio_sst_mab_sig")

plot_for_dch <- c("pole_sst_slope", "pole_year_slope")
plot_for_dch_lab <- c("pole_sst_sig", "pole_year_sig")

zzmelt <- data.frame(cbind(corrtest_biolat$group, corrtest_biolat$slope))
zzmelt$pval <- corrtest_biolat$pval
zzmelt$src <- "biolat"

zzz <- data.frame(cbind(corrtest_biomass$group[corrtest_biomass$stratcat == "GOMstrat"], corrtest_biomass$slope[corrtest_biomass$stratcat == "GOMstrat"]))
zzz$pval <- corrtest_biomass$pval[corrtest_biomass$stratcat == "GOMstrat"]
zzz$src <- "GOMstrat"

zzs <- data.frame(cbind(corrtest_biomass$group[corrtest_biomass$stratcat == "GBstrat"], corrtest_biomass$slope[corrtest_biomass$stratcat == "GBstrat"]))
zzs$pval <- corrtest_biomass$pval[corrtest_biomass$stratcat == "GBstrat"]
zzs$src <- "GBstrat"

zzd <- data.frame(cbind(corrtest_biomass$group[corrtest_biomass$stratcat == "SNEstrat"], corrtest_biomass$slope[corrtest_biomass$stratcat == "SNEstrat"]))
zzd$pval <- corrtest_biomass$pval[corrtest_biomass$stratcat == "SNEstrat"]
zzd$src <- "SNEstrat"

zzf <- data.frame(cbind(corrtest_biomass$group[corrtest_biomass$stratcat == "MABstrat"], corrtest_biomass$slope[corrtest_biomass$stratcat == "MABstrat"]))
zzf$pval <- corrtest_biomass$pval[corrtest_biomass$stratcat == "MABstrat"]
zzf$src <- "MABstrat"

zzmelt <- rbind(zzmelt, zzz)
zzmelt <- rbind(zzmelt, zzs)
zzmelt <- rbind(zzmelt, zzd)
zzmelt <- rbind(zzmelt, zzf)

zzmelt$sigstar <- ifelse(zzmelt$pval <= 0.003333333333, "*", NA)
zzmelt$posneg <- ifelse(zzmelt$X2 > 0 & zzmelt$sigstar == "*", "+", NA)
zzmelt$posneg <- ifelse(zzmelt$X2 < 0 & zzmelt$sigstar == "*", "-", zzmelt$posneg)
zzmelt$X2 <- ifelse(zzmelt$sigstar == "*", as.numeric(zzmelt$X2), NA)

# zzmelt_slope <- unique(melt(zzall, id.vars = c("survdat.COMNAME", "group"), measure.vars = plot_for_dch))
# zzmelt_slope <- zzmelt_slope[!is.na(zzmelt_slope$group),]
# zzmelt_slope$nm <- substr(zzmelt_slope$variable, 9, 12)
# zzmelt_slope$nm <- substr(zzmelt_slope$variable, 6, 9)
# 
# zzmelt_lab <- unique(melt(zzall, id.vars = c("survdat.COMNAME", "group"), measure.vars = plot_for_dch_lab))
# zzmelt_lab <- zzmelt_lab[!is.na(zzmelt_lab$group),]
# zzmelt_lab$nm <- substr(zzmelt_lab$variable, 9, 12)
# zzmelt_lab$nm <- substr(zzmelt_lab$variable, 6, 9)
# 
# zzmelt_slope <- merge(zzmelt_slope, zzmelt_lab, by.x=c("survdat.COMNAME", "nm", "group"), by.y=c("survdat.COMNAME", "nm", "group"))
# colnames(zzmelt_slope)[4:7] <- c("slopegrp", "slopeval", "siggrp", "sigval")
# zzmelt_slope$sigstar <- as.character(ifelse(zzmelt_slope$sigval == 1, "*", NA))
# zzmelt_slope$sigslope <- ifelse(zzmelt_slope$sigval == 1, zzmelt_slope$slopeval, NA)
# 


# biocorr <- merge(corrtest_biomass, critter_names[,c(2,3)], by.x="group", by.y="survdat.COMNAME", all.x=T)
#  #, all.x=T
# biocorr$sigslop <- ifelse(biocorr$sig == 1, biocorr$slope, NA)
# biocorr$sigstar <- ifelse(biocorr$sigslop > 0, "+", "-")
# dchcorr <- merge(corrtest_biolat, critter_names[,c(2,3)], by.x="group", by.y="survdat.COMNAME", all.x=T) #, all.x=T
# dchcorr$sigslop <- ifelse(dchcorr$sig == 1, dchcorr$slope, NA)
# dchcorr$sigstar <- ifelse(dchcorr$sigslop > 0, "+", "-")
# dchcorr <- dchcorr[dchcorr$group.y == "zoop" | dchcorr$group.y == "shrimp" | !is.na(dchcorr$group.y),]

cols <- colorRampPalette(c("#59919f", "#92b9d0", "#f3de96", "#fcdf75","#f98f8d", "#e51d50", "#820338"))
my.at <- c(-330, seq(-40,-5, 5), seq(-4.75, 4.75, 0.25), seq(5, 40, 5), seq(50, 180, 20), 450)


### SUPPLEMENTARY (BIG FIGS)
setwd("~/Dropbox/my_papers/InPrep_Stepanuketal_GOMprey/")
# ggplot(biocorr[biocorr$stratcat == "MABstrat" & !is.na(biocorr$group.y) & biocorr$group.y != "zoop" & biocorr$group.y != "shrimp",] %>% arrange(type), aes(x=reorder(type, group.y), y=stringr::str_to_sentence(group), fill=sigslop)) +
pdf("spp_bt_supplementary_20220421.pdf", width=8, height=8)
ggplot(biocorr[biocorr$type == "BT" & !is.na(biocorr$group.y) & biocorr$group.y != "zoop" & biocorr$group.y != "shrimp",] %>% arrange(stratcat), aes(x=factor(stratcat, levels=c("GOMstrat", "GBstrat", "SNEstrat", "MABstrat")), y=stringr::str_to_sentence(group), fill=sigslop)) +
  geom_tile(color="black") +
  #geom_tile(aes(color=stratcat), color="black", na.value=NA) + #facet_wrap(.~group) +
  geom_text(aes(label=sigstar), color="black") +
  scale_x_discrete(NULL, expand = c(0, 0), labels = c("Gulf of Maine", "George's Bank", "Southern New England", "Mid-Atlantic Bight"), position="top") +
  scale_y_discrete(NULL, expand = c(0, 0)) +
  scale_fill_gradientn(colours = rev(viridis(10, option="F")[4:10]),
                       #values = rescale(x = c(-2, 0, 2, 5, 10, 20, 30, 40), from = c(-2,40)),
                       #oob = squish,
                       #limits = c(-2, 40),
                       name="Slope",
                       na.value = "NA") +
  #scale_fill_viridis(option = "C") +
  #colorspace::scale_fill_continuous_divergingx(palette = 'Tropic', mid = 0.01, l3 = 0, p3 = .2, p4 = .6) + 
  theme_bw() +
  #scale_fill_gradient(trans = "exp") +
  #coord_cartesian(clip = "off") +
  theme(axis.text.x = element_text(angle=60,vjust = 0.5, hjust = 0),
        strip.placement = "inside",
        strip.background = element_rect(fill=NA,colour="grey50"),
        strip.text.y.left = element_text(angle=0),
        panel.spacing=unit(0,"cm"),
        panel.border = element_rect(color="black")) +
  facet_grid(group.y~., space="free_y", scales="free_y", switch="y") +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5), color="black")
dev.off()
  

ggplot(dchcorr %>% arrange(type), aes(x=reorder(type, group.y), y=stringr::str_to_sentence(group), fill=sigslop)) +
  geom_tile(color="black") + #facet_wrap(.~group) +
  #geom_text(aes(label=sigstar), color="black") +
  scale_x_discrete(NULL, expand = c(0, 0), position="top") +
  scale_y_discrete(NULL, expand = c(0, 0)) +
#this is for the poleward shift slopes
scale_fill_gradientn(colours = rev(viridis(10, option="F")[4:10]),
                     values = rescale(x = c(-100, -10, 0, 2, 5, 10, 50, 100, 200, 300, 400, 450), from = c(-100,450)),
                     oob = squish,
                     limits = c(-100, 450),
                     name="Slope") +
  theme_bw() +
  #scale_fill_gradient(trans = "exp") +
  #coord_cartesian(clip = "off") +
  theme(axis.text.x = element_text(angle=60,vjust = 0.5, hjust = 0),
        strip.placement = "inside",
        strip.background = element_rect(fill=NA,colour="grey50"),
        strip.text.y.left = element_text(angle=0),
        panel.spacing=unit(0,"cm")) +
  facet_grid(group.y~., space="free_y", scales="free_y", switch="y")


### Main text (proportion figs)
# biocorr_agg <- biocorr
# biocorr_agg <- aggregate(signeg ~ group.y + type + stratcat, data=biocorr, FUN=sum)
# biocorr_tab <- aggregate(signeg ~ group.y + type + stratcat, data=biocorr, FUN=length)
# colnames(biocorr_tab)[4] <- "length"
# biocorr_agg <- merge(biocorr_agg, biocorr_tab, by.x=c("group.y", "type", "stratcat"), by.y=c("group.y", "type", "stratcat"))
# biocorr_agg$prop <- biocorr_agg$sig / biocorr_agg$length
# biocorr_agg$src <- ifelse(biocorr_agg$group.y == "zoop" | biocorr_agg$group.y == "shrimp", "ECOMON", "BTS")
# 
# dchcorr_agg <- dchcorr
# dchcorr_agg <- aggregate(sigpos ~ group.y + type , data=dchcorr, FUN=sum)
# dchcorr_tab <- aggregate(sigpos ~ group.y + type , data=dchcorr, FUN=length)
# colnames(dchcorr_tab)[3] <- "length"
# dchcorr_agg <- merge(dchcorr_agg, dchcorr_tab, by.x=c("group.y", "type"), by.y=c("group.y", "type"))
# dchcorr_agg$prop <- dchcorr_agg$sig / dchcorr_agg$length
# dchcorr_agg$src <- ifelse(dchcorr_agg$group.y == "zoop" | dchcorr_agg$group.y == "shrimp", "ECOMON", "BTS")

zzmelt$src2 <- zzmelt$src
zzmelt$src2 <- ifelse(zzmelt$src2 == "biolat", "Poleward Shift", zzmelt$src2)
zzmelt$src2 <- ifelse(zzmelt$src2 == "GOMstrat", "Gulf of Maine", zzmelt$src2)
zzmelt$src2 <- ifelse(zzmelt$src2 == "GBstrat", "George's Bank", zzmelt$src2)
zzmelt$src2 <- ifelse(zzmelt$src2 == "SNEstrat", "Southern New England", zzmelt$src2)
zzmelt$src2 <- ifelse(zzmelt$src2 == "MABstrat", "Mid-Atlantic Bight", zzmelt$src2)

zzmelt$src2 <- factor(zzmelt$src2, levels=c("Poleward Shift", "Gulf of Maine", "George's Bank", "Southern New England", "Mid-Atlantic Bight"))

zzmelt$group2 <- zzmelt$X1
zzmelt$group2 <- ifelse(zzmelt$group2 == "Large_Gadid", "Large Gadid", zzmelt$group2)
zzmelt$group2 <- ifelse(zzmelt$group2 == "Meso_pelagic", "Mesopelagic Fish", zzmelt$group2)
zzmelt$group2 <- ifelse(zzmelt$group2 == "Misc_Fish", "Miscellaneous Fish", zzmelt$group2)
zzmelt$group2 <- ifelse(zzmelt$group2 == "Northern_shrimp", "Northern Shrimp", zzmelt$group2)
zzmelt$group2 <- ifelse(zzmelt$group2 == "Shrimp_BT", "Large Shrimp", zzmelt$group2)
zzmelt$group2 <- ifelse(zzmelt$group2 == "Small_Gadid", "Small Gadid", zzmelt$group2)
zzmelt$group2 <- ifelse(zzmelt$group2 == "Squid_BT", "Squid", zzmelt$group2)
zzmelt$group2 <- ifelse(zzmelt$group2 == "sm_zoop", "Copepods", zzmelt$group2)
zzmelt$group2 <- ifelse(zzmelt$group2 == "lg_zoop", "Mysids & Amphipods", zzmelt$group2)
zzmelt$group2 <- ifelse(zzmelt$group2 == "shrimp", "Small Shrimp", zzmelt$group2)
zzmelt$group2 <- ifelse(zzmelt$group2 == "euphausiid", "Euphausiid", zzmelt$group2)

zzmelt$group2 <- factor(zzmelt$group2, levels=c("Small Gadid", "Large Gadid", "Mesopelagic Fish", "Miscellaneous Fish", "Scombrid", "Sandlance", "Flatfish", "Clupeid", "Large Shrimp", "Northern Shrimp", "Squid", "Euphausiid", "Copepods", "Mysids & Amphipods", "Small Shrimp", ""))

setwd("~/Dropbox/my_papers/InPrep_Stepanuketal_GOMprey/")
pdf("biomass_trends_20220723.pdf", width=8, height=6)
ggplot(zzmelt[zzmelt$group2 != "",] %>% arrange(group2), aes(x=src2, y=group2, fill=as.numeric(X2))) +
  geom_tile(color="black") + #facet_wrap(.~group) +
  geom_text(aes(label=posneg), color="white") +
  scale_x_discrete(NULL, expand = c(0, 0), position="top") +
  scale_y_discrete(NULL, expand = c(0, 0)) +
  theme_classic() +
  labs(title = "Prey group changes - 1980-2019", fill="Slope") +
  scale_fill_viridis(option = "C",  na.value="NA") +
  #scale_fill_gradient(trans = "exp") +
  #coord_cartesian(clip = "off") +
  theme(axis.text.x = element_text(angle=60,vjust = 0.5, hjust = 0),
        strip.placement = "inside",
        strip.background = element_rect(fill=NA,colour="grey50"),
        strip.text.y.left = element_text(angle=0),
        panel.spacing=unit(0.2,"cm"),
        panel.grid.minor = element_line(color="black", size=3),
        panel.grid.major = element_blank(),
        panel.border = element_rect(color="black", fill=NA),
        text = element_text(size=16)
        ) #+
  # facet_grid(src~., space="free_y", scales="free_y", switch="y") +
  # facet_wrap(.~factor(stratcat, levels=c("GOMstrat", "GBstrat", "SNEstrat", "MABstrat"),labels = c("Gulf of Maine", "George's Bank", "Southern New England", "Mid-Atlantic Bight")))
dev.off()



pdf("poleward_proportion_20220422.pdf", width=8, height=6)
ggplot(dchcorr_agg %>% arrange(src), aes(x=reorder(type, group.y), y=stringr::str_to_sentence(group.y), fill=prop)) +
  geom_tile(color="black") + #facet_wrap(.~group) +
  geom_text(aes(label=round(prop, 2)), color="black") +
  scale_x_discrete(NULL, expand = c(0, 0), position="top") +
  scale_y_discrete(NULL, expand = c(0, 0)) +
  theme_bw() +
  labs(title = "Poleward shift vs. environmental covariates") +
  scale_fill_viridis(option="C", limits=c(0.0001,0.5)) + #direction = -1
  #scale_fill_gradient(trans = "exp") +
  #coord_cartesian(clip = "off") +
  theme(axis.text.x = element_text(angle=60,vjust = 0.5, hjust = 0),
        strip.placement = "inside",
        strip.background = element_rect(fill=NA,colour="grey50"),
        strip.text.y.left = element_text(angle=0),
        panel.spacing=unit(0,"cm")) +
  facet_grid(src~., space="free_y", scales="free_y", switch="y")
dev.off()


##### breakdown by species by location

###########
biocorr_sppsig <- aggregate(sig ~ group + group.y +stratcat, biocorr[biocorr$type != "YEAR",], FUN=sum)
dchcorr_sppsig <- aggregate(sig ~ group + group.y, dchcorr[dchcorr$type != "YEAR",], FUN=sum)
biocorr_sppsig$sigyn <- ifelse(biocorr_sppsig$sig > 0 , 1 , 0)
dchcorr_sppsig$sigyn <- ifelse(dchcorr_sppsig$sig > 0 , 1 , 0)


biocorr_sppsig <- merge(t_region, biocorr_sppsig, by.x=c("stratcat", "survdat.COMNAME"), by.y=c("stratcat", "group"))
dchcorr_sppsig <- merge(t_region, dchcorr_sppsig, by.x=c( "survdat.COMNAME"), by.y=c( "group"))

#zzsignif <- aggregate(Ba_prey_scaled ~ YEAR + stratcat + pole_vs_sst, zzt[!is.na(zzt$group),], FUN=sum)
zzsignif <- aggregate(Ba_prey_scaled ~ YEAR + stratcat + sigyn, biocorr_sppsig, FUN=sum)
zzsignif[,4][zzsignif[,4] == 0] <- NA
zza <- ggplot(zzsignif[zzsignif$YEAR >= 1990,], aes(x=YEAR, y=Ba_prey_scaled, color=factor(as.character(sigyn)))) +
  geom_point(size=0.6) +
  geom_smooth(method="loess", size=0.6) +
  #geom_smooth(aes(color=factor(as.character(pole_vs_year))), method="loess", linetype = "dashed") +
  facet_wrap(.~factor(stratcat, levels = c("GOMstrat", "GBstrat", "SNEstrat", "MABstrat")), ncol=1) +
  labs(color = "Significant", y="", x="", title="Minke") +
  ylim(c(0,30))+
  theme_minimal() +
  theme(legend.position = "none",
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        text = element_text(size=14)) +
  scale_color_manual(values=c("blue", "orange"))

zzleg <- ggpubr::get_legend(ggplot(zzsignif, aes(x=YEAR, y=Ba_prey_scaled, color=factor(as.character(sigyn)))) +
                              geom_point(size=0.6) +
                              labs(color = "Significant", y="", x="", title="Minke") +
                              theme_minimal() +
                              ylim(c(0,30))+
                              theme(legend.position = "bottom",
                                    plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
                                    text = element_text(size=14)) +
                              scale_color_manual(values=c("blue", "orange")))


#zzsignif <- aggregate(Bb_prey_scaled ~ YEAR + stratcat + pole_vs_sst, zzt[!is.na(zzt$group),], FUN=sum)
zzsignif <- aggregate(Bb_prey_scaled ~ YEAR + stratcat + sigyn, biocorr_sppsig, FUN=sum)
zzsignif[,4][zzsignif[,4] == 0] <- NA
zzb <- ggplot(zzsignif[zzsignif$YEAR >= 1990,], aes(x=YEAR, y=Bb_prey_scaled, color=factor(as.character(sigyn)))) +
  geom_point(size=0.6) +
  geom_smooth(method="loess", size=0.6) +
  #geom_smooth(aes(color=factor(as.character(pole_vs_year))), method="loess", linetype = "dashed") +
  facet_wrap(.~factor(stratcat, levels = c("GOMstrat", "GBstrat", "SNEstrat", "MABstrat")), ncol=1) +
  labs(color = "Significant", y="", x="", title="Sei") +
  ylim(c(0,30))+
  theme_minimal() +
  theme(legend.position = "none",
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        text = element_text(size=14)) +
  scale_color_manual(values=c("blue", "orange"))

#zzsignif <- aggregate(Bp_prey_scaled ~ YEAR + stratcat + pole_vs_sst, zzt[!is.na(zzt$group),], FUN=sum)
zzsignif <- aggregate(Bp_prey_scaled ~ YEAR + stratcat + sigyn, biocorr_sppsig, FUN=sum)
zzsignif[,4][zzsignif[,4] == 0] <- NA
zzc <- ggplot(zzsignif[zzsignif$YEAR >= 1990,], aes(x=YEAR, y=Bp_prey_scaled, color=factor(as.character(sigyn)))) +
  geom_point(size=0.6) +
  geom_smooth( method="loess", size=0.6) +
  #geom_smooth(aes(color=factor(as.character(pole_vs_year))), method="loess", linetype = "dashed") +
  facet_wrap(.~factor(stratcat, levels = c("GOMstrat", "GBstrat", "SNEstrat", "MABstrat")), ncol=1) +
  labs(color = "Significant", y="", x="", title="Fin") +
  ylim(c(0,30))+
  theme_minimal() +
  theme(legend.position = "none",
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        text = element_text(size=14)) +
  scale_color_manual(values=c("blue", "orange"))

#zzsignif <- aggregate(Mn_prey_scaled ~ YEAR + stratcat + pole_vs_sst, zzt[!is.na(zzt$group),], FUN=sum)
zzsignif <- aggregate(Mn_prey_scaled ~ YEAR + stratcat + sigyn, biocorr_sppsig, FUN=sum)
zzsignif[,4][zzsignif[,4] == 0] <- NA
zzd <- ggplot(zzsignif[zzsignif$YEAR >= 1990,], aes(x=YEAR, y=Mn_prey_scaled, color=factor(as.character(sigyn)))) +
  geom_point(size=0.6) +
  geom_smooth( method="loess", size=0.6) +
  #geom_smooth(aes(color=factor(as.character(pole_vs_year))), method="loess", linetype = "dashed") +
  facet_wrap(.~factor(stratcat, levels = c("GOMstrat", "GBstrat", "SNEstrat", "MABstrat")), ncol=1) +
  labs(color = "Significant", y="", x="", title="Humpback") +
  ylim(c(0,30))+
  theme_minimal() +
  theme(legend.position = "none",
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        text = element_text(size=14)) +
  scale_color_manual(values=c("blue", "orange"))

zze <- gridExtra::grid.arrange(zza, zzb, zzc, zzd, ncol=4, left = "Prey Availability", bottom="Year")
pdf("sig_groups_biomass_20220422.pdf", width=8, height=6)
gridExtra::grid.arrange(zze, zzleg, ncol=1, heights = c(15,1))
dev.off()

###### now do for poleward distance
#zzsignif <- aggregate(Ba_prey_scaled ~ YEAR + stratcat + pole_vs_sst, zzt[!is.na(zzt$group),], FUN=sum)
zzsignif <- aggregate(Ba_prey_scaled ~ YEAR  + sigyn, dchcorr_sppsig, FUN=sum)
zzsignif[,3][zzsignif[,3] == 0] <- NA
zza <- ggplot(zzsignif[zzsignif$YEAR >= 1990,], aes(x=YEAR, y=Ba_prey_scaled, color=factor(as.character(sigyn)))) +
  geom_point(size=0.6) +
  geom_smooth(method="loess", size=0.6) +
  #geom_smooth(aes(color=factor(as.character(pole_vs_year))), method="loess", linetype = "dashed") +
  labs(color = "Significant", y="", x="", title="Minke") +
  ylim(c(0,30))+
  theme_minimal() +
  theme(legend.position = "none",
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        text = element_text(size=14)) +
  scale_color_manual(values=c("blue", "orange"))

zzleg <- ggpubr::get_legend(ggplot(zzsignif, aes(x=YEAR, y=Ba_prey_scaled, color=factor(as.character(sigyn)))) +
                              geom_point(size=0.6) +
                              labs(color = "Significant", y="", x="", title="Minke") +
                              theme_minimal() +
                              theme(legend.position = "bottom",
                                    plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
                                    text = element_text(size=14)) +
                              scale_color_manual(values=c("blue", "orange")))


#zzsignif <- aggregate(Bb_prey_scaled ~ YEAR + stratcat + pole_vs_sst, zzt[!is.na(zzt$group),], FUN=sum)
zzsignif <- aggregate(Bb_prey_scaled ~ YEAR + sigyn, dchcorr_sppsig, FUN=sum)
zzsignif[,3][zzsignif[,3] == 0] <- NA
zzb <- ggplot(zzsignif[zzsignif$YEAR >= 1990,], aes(x=YEAR, y=Bb_prey_scaled, color=factor(as.character(sigyn)))) +
  geom_point(size=0.6) +
  geom_smooth(method="loess", size=0.6) +
  #geom_smooth(aes(color=factor(as.character(pole_vs_year))), method="loess", linetype = "dashed") + +
  labs(color = "Significant", y="", x="", title="Sei") +
  ylim(c(0,30))+
  theme_minimal() +
  theme(legend.position = "none",
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        text = element_text(size=14)) +
  scale_color_manual(values=c("blue", "orange"))

#zzsignif <- aggregate(Bp_prey_scaled ~ YEAR + stratcat + pole_vs_sst, zzt[!is.na(zzt$group),], FUN=sum)
zzsignif <- aggregate(Bp_prey_scaled ~ YEAR + sigyn, dchcorr_sppsig, FUN=sum)
zzsignif[,3][zzsignif[,3] == 0] <- NA
zzc <- ggplot(zzsignif[zzsignif$YEAR >= 1990,], aes(x=YEAR, y=Bp_prey_scaled, color=factor(as.character(sigyn)))) +
  geom_point(size=0.6) +
  geom_smooth( method="loess", size=0.6) +
  #geom_smooth(aes(color=factor(as.character(pole_vs_year))), method="loess", linetype = "dashed") +
  labs(color = "Significant", y="", x="", title="Fin") +
  ylim(c(0,30))+
  theme_minimal() +
  theme(legend.position = "none",
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        text = element_text(size=14)) +
  scale_color_manual(values=c("blue", "orange"))

#zzsignif <- aggregate(Mn_prey_scaled ~ YEAR + stratcat + pole_vs_sst, zzt[!is.na(zzt$group),], FUN=sum)
zzsignif <- aggregate(Mn_prey_scaled ~ YEAR  + sigyn, dchcorr_sppsig, FUN=sum)
zzsignif[,3][zzsignif[,3] == 0] <- NA
zzd <- ggplot(zzsignif[zzsignif$YEAR >= 1990,], aes(x=YEAR, y=Mn_prey_scaled, color=factor(as.character(sigyn)))) +
  geom_point(size=0.6) +
  geom_smooth( method="loess", size=0.6) +
  #geom_smooth(aes(color=factor(as.character(pole_vs_year))), method="loess", linetype = "dashed") +
  labs(color = "Significant", y="", x="", title="Humpback") +
  ylim(c(0,30))+
  theme_minimal() +
  theme(legend.position = "none",
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        text = element_text(size=14)) +
  scale_color_manual(values=c("blue", "orange"))

zze <- gridExtra::grid.arrange(zza, zzb, zzc, zzd, ncol=2, left = "Prey Availability", bottom="Year")

pdf("sig_groups_dch_20220422.pdf", width=8, height=6)
gridExtra::grid.arrange(zze, zzleg, ncol=1, heights = c(15,1))
dev.off()
