#stratified bottom temp
bt_survdat <- survdat_fall
bt_survdat$BIOMASS <- bt_survdat$BOTTEMP
bt_survdat$SVSPP <- 1
bt_survdat$COMNAME <- "b"
bt_survdat$ABUNDANCE <- 1
bt_survdat <- unique(bt_survdat)

bt_test_gom <- strat_prep(data.table(bt_survdat[bt_survdat$stratcat == "GOMstrat",-c(25:26)]))
bt_test_gom <- unique(bt_test_gom)

bt_test_gb <- strat_prep(data.table(bt_survdat[bt_survdat$stratcat == "GBstrat",-c(25:26)]))
bt_test_gb <- unique(bt_test_gb)

bt_test_sne <- strat_prep(data.table(bt_survdat[bt_survdat$stratcat == "SNEstrat",-c(25:26)]))
bt_test_sne <- unique(bt_test_sne)

bt_test_mab <- strat_prep(data.table(bt_survdat[bt_survdat$stratcat == "MABstrat",-c(25:26)]))
bt_test_mab <- unique(bt_test_mab)

bt_gom <- bt_test_gom %>%
  group_by(YEAR) %>%
  summarise(wt_mean=weighted.mean(BOTTEMP,W.h, na.rm=TRUE))
bt_gb <- bt_test_gb %>%
  group_by(YEAR) %>%
  summarise(wt_mean=weighted.mean(BOTTEMP,W.h, na.rm=TRUE))
bt_sne <- bt_test_sne %>%
  group_by(YEAR) %>%
  summarise(wt_mean=weighted.mean(BOTTEMP,W.h, na.rm=TRUE))
bt_mab <- bt_test_mab %>%
  group_by(YEAR) %>%
  summarise(wt_mean=weighted.mean(BOTTEMP,W.h, na.rm=TRUE))

#Assign category for each region
bt_gom$stratcat <- "GOMstrat"
bt_gb$stratcat <- "GBstrat"
bt_sne$stratcat <- "SNEstrat"
bt_mab$stratcat <- "MABstrat"

#Combine regions 
bt_region <- rbind(bt_gom, bt_gb)
bt_region <- rbind(bt_region, bt_sne)
bt_region <- rbind(bt_region, bt_mab)

colnames(bt_region)[2] <- "bot_temp"
ggplot(bt_region, aes(x=YEAR, y=bot_temp, color=stratcat)) + geom_point() + geom_smooth()
