## Aggregate and plot out scaled prey for each predator species

t_agg_all_s <- aggregate(cbind(Bp_prey_scaled, Ba_prey_scaled, Mn_prey_scaled, Bb_prey_scaled) ~ YEAR + stratcat, t_region, function(x){sum(x, na.rm=T)})

t_agg_all_s_melt <- melt(t_agg_all_s, id.vars=c("YEAR", "stratcat"))

t_agg_all_s_melt$spp2 <- ifelse(t_agg_all_s_melt$variable == "Bp_prey_scaled", "Fin Whale", NA)
t_agg_all_s_melt$spp2 <- ifelse(t_agg_all_s_melt$variable == "Ba_prey_scaled", "Minke Whale", t_agg_all_s_melt$spp2)
t_agg_all_s_melt$spp2 <- ifelse(t_agg_all_s_melt$variable == "Bb_prey_scaled", "Sei Whale", t_agg_all_s_melt$spp2)
t_agg_all_s_melt$spp2 <- ifelse(t_agg_all_s_melt$variable == "Mn_prey_scaled", "Humpback Whale", t_agg_all_s_melt$spp2)


wh_pt <- ggplot(data=t_agg_all_s_melt[t_agg_all_s_melt$YEAR >=1990,], aes(x=YEAR, y=value, color=spp2, fill=spp2)) + 
  #Below plots any category with a biomass - aka BT groups
  geom_point(alpha=0.8, pch=16) + 
  geom_smooth( method="loess", span=0.8, se=T) +
  
  #Below plots any category with an abundance - aka ecomon groups
  #geom_point(data=dual_comb, aes(y=zoop_abundance/sss), alpha=0.4, pch=16) +
  #geom_smooth(data=dual_comb, aes(y=zoop_abundance/sss), method="lm", se=T) +
  
  scale_y_continuous(
    # Features of the first axis
    #    name = "Bottom Trawl Biomass",
    # Add a second axis and specify its features
    #sec.axis = sec_axis( trans=~.*sss, name="Ecomon Zoop Abundance"),
    #trans="log10"
  ) +
  scale_color_manual(values = c("Minke Whale" = "#5BC0EB",
                                "Sei Whale"="#FDE008",
                                "Fin Whale"="#9BC53D",
                                "Humpback Whale" = "#E55934")) +
  scale_fill_manual(values = c("Minke Whale" = "#5BC0EB",
                               "Sei Whale"="#FDE008",
                               "Fin Whale"="#9BC53D",
                               "Humpback Whale" = "#E55934")) +
  
  facet_wrap(.~factor(stratcat, levels=c("GOMstrat", "GBstrat", "SNEstrat", "MABstrat"), labels = c("Gulf of Maine", "George's Bank", "Southern New England", "Mid-Atlantic Bight"))) +
  theme_minimal() +
  labs(x= "Year", y="Prey Availability", color = "Species", fill="Species") +
  theme(text = element_text(size=16),
        axis.text.x = element_text(angle=45, vjust = 1, hjust=1),
        legend.position = "right",
        plot.margin = unit(c(0, 0, 0, 0.5), "cm"))

setwd("~/Dropbox/my_papers/InPrep_Stepanuketal_GOMprey/")
pdf("whales_20220422.pdf", width=8, height=6)
wh_pt
dev.off()


####### prey by year - no scale
t_agg_all_gr <- aggregate(strat.biomass ~ YEAR + stratcat + group, t_region, function(x){sum(x, na.rm=T)})

t_agg_all_gr_melt <- melt(t_agg_all_gr, id.vars=c("YEAR", "stratcat", "group"))

t_agg_all_gr_melt$grp <- ifelse(t_agg_all_gr_melt$group == "Clupeid", "Clupeids", NA)
t_agg_all_gr_melt$grp <- ifelse(t_agg_all_gr_melt$group == "Flatfish", "Flatfish", t_agg_all_gr_melt$grp)
t_agg_all_gr_melt$grp <- ifelse(t_agg_all_gr_melt$group == "Large_Gadid", "Large Gadids", t_agg_all_gr_melt$grp)
t_agg_all_gr_melt$grp <- ifelse(t_agg_all_gr_melt$group == "Meso_pelagic", "Mesopelagic Fish", t_agg_all_gr_melt$grp)
t_agg_all_gr_melt$grp <- ifelse(t_agg_all_gr_melt$group == "Misc_Fish", "Miscellaneous Fish", t_agg_all_gr_melt$grp)
t_agg_all_gr_melt$grp <- ifelse(t_agg_all_gr_melt$group == "Sandlance", "Sandlance", t_agg_all_gr_melt$grp)
t_agg_all_gr_melt$grp <- ifelse(t_agg_all_gr_melt$group == "Scombrid", "Scombrids", t_agg_all_gr_melt$grp)
t_agg_all_gr_melt$grp <- ifelse(t_agg_all_gr_melt$group == "Shrimp_BT", "Shrimp (BTS)", t_agg_all_gr_melt$grp)
t_agg_all_gr_melt$grp <- ifelse(t_agg_all_gr_melt$group == "Small_Gadid", "Small Gadids", t_agg_all_gr_melt$grp)
t_agg_all_gr_melt$grp <- ifelse(t_agg_all_gr_melt$group == "Squid_BT", "Squid", t_agg_all_gr_melt$grp)


fish_pt <- ggplot(data=t_agg_all_gr_melt[t_agg_all_gr_melt$YEAR >=1990,], aes(x=YEAR, y=value, color=grp, fill=grp)) + 
  #Below plots any category with a biomass - aka BT groups
  geom_point(alpha=0.8, size=1) + 
  geom_smooth( method="loess", span=0.8, se=T, alpha=0.2) +
  scale_y_continuous()+
  facet_wrap(.~factor(stratcat, levels=c("GOMstrat", "GBstrat", "SNEstrat", "MABstrat"), labels = c("Gulf of Maine", "George's Bank", "Southern New England", "Mid-Atlantic Bight"))) +
  theme_minimal() +
  labs(x= "Year", y="Stratified biomass", color = "Prey Groups", fill="Prey Groups") +
  scale_color_manual(values=c("#D00000", "#FFBA08", "#73237A", "#98FF1A", "#8FE388", "#1B998B", "#3185FC", "#5D2E8C",  "#FF7B9C", "#FF9B85"), guide= guide_legend(ncol=3)) +
  scale_fill_manual(values=c("#D00000", "#FFBA08", "#73237A", "#98FF1A", "#8FE388", "#1B998B", "#3185FC", "#5D2E8C",  "#FF7B9C", "#FF9B85"), guide= guide_legend(ncol=3)) +
  theme(text = element_text(size=18),
        axis.text.x = element_text(angle=45, vjust = 1, hjust=1),
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0.5), "cm")) 
#fish_leg <- ggpubr::get_legend(fish_pt)
fish_pt

zoop_pt <- ggplot(data=zpd_region1[zpd_region1$YEAR >= 1990,], aes(x=YEAR, y=strat.biomass, color=group, fill=group)) + 
  #Below plots any category with a biomass - aka BT groups
  #geom_point(alpha=0.8, size=1) + 
  geom_smooth( method="loess", span=0.8, se=T, alpha=0.2) +
  
  scale_y_continuous(
    # #Features of the first axis
    #    name = "Bottom Trawl Biomass",
    # #Add a second axis and specify its features
    # sec.axis = sec_axis( trans=~.*100, name="Ecomon Zoop Abundance")
  ) +
  # scale_color_manual(values = c("Minke Whale" = "#5BC0EB",
  #                               "Sei Whale"="#FDE008",
  #                               "Fin Whale"="#9BC53D",
  #                               "Humpback Whale" = "#E55934")) +
  # scale_fill_manual(values = c("Minke Whale" = "#5BC0EB",
  #                              "Sei Whale"="#FDE008",
  #                              "Fin Whale"="#9BC53D",
  #                              "Humpback Whale" = "#E55934")) +
  
  facet_wrap(.~factor(stratcat, levels=c("GOMstrat", "GBstrat", "SNEstrat", "MABstrat"), labels = c("Gulf of Maine", "George's Bank", "Southern New England", "Mid-Atlantic Bight"))) +
  scale_color_manual(labels=c("Shrimp (EcoMon)", "Zooplankton"), values=c("#94C595", "#372772")) +
  scale_fill_manual(labels=c("Shrimp (EcoMon)", "Zooplankton"), values=c("#94C595", "#372772")) +
  theme_minimal() +
  labs(x= "Year", y="Stratified abundance", color = "", fill="") +
  theme(text = element_text(size=18),
        axis.text.x = element_text(angle=45, vjust = 1, hjust=1),
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0.5), "cm"))
#zoop_leg <- ggpubr::get_legend(zoop_pt)
zoop_pt



setwd("~/Dropbox/my_papers/InPrep_Stepanuketal_GOMprey/")
pdf("prey_20220422.pdf", width=12, height=8)
cowplot::plot_grid(fish_pt, zoop_pt,fish_leg, zoop_leg, ncol=2, rel_heights = c(0.8, 0.2))
dev.off()
