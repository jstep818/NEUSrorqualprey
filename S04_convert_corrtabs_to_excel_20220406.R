#Calculate lm outputs and put into table

# dch by envr vars
#disaggregate
table_corrtest_dch <-  reshape(corrtest_biolat, idvar = "group", timevar = "type",  direction = "wide")
table_corrtest_dch <- merge(table_corrtest_dch, critter_names[,2:3], by.x = "group", by.y="survdat.COMNAME")

table_corrtest_dch <- table_corrtest_dch %>%
  rename(TAXA.NAME= group) %>%
  rename(group = group.y) %>%
  arrange(group) %>%
  mutate_at(2:22, funs(as.numeric(as.character(.)))) %>%
  mutate_at(c(3,4,6), funs(round(., 0))) %>%
  mutate_at(c(5,7,8), funs(round(., 2))) %>%
  mutate_at("TAXA.NAME", stringr::str_to_sentence)

table_corrtest_bio <-  reshape(corrtest_biomass, idvar = c("group", "stratcat"), timevar = "type",  direction = "wide")
table_corrtest_bio <- merge(table_corrtest_bio, critter_names[,2:3], by.x = "group", by.y="survdat.COMNAME")

table_corrtest_bio <- table_corrtest_bio %>%
  rename(TAXA.NAME= group) %>%
  rename(group = group.y) %>%
  arrange(stratcat, group) %>%
  mutate_at(3:23, funs(as.numeric(as.character(.)))) %>%
  mutate_at(c(3,6, 10, 13, 17, 20), funs(round(., 0))) %>%
  mutate_at(c(4,5,7,8, 11, 12, 14, 15, 18, 19, 21, 22), funs(round(., 2))) %>%
  mutate_at("TAXA.NAME", stringr::str_to_sentence)

write.csv(table_corrtest_dch, "~/Dropbox/my_papers/InPrep_Stepanuketal_GOMprey/table_corrtest_dch.csv")
write.csv(table_corrtest_bio, "~/Dropbox/my_papers/InPrep_Stepanuketal_GOMprey/table_corrtest_bio.csv")


zpd_spp_table <- ZPD_fall
zpd_spp_table <- aggregate(value ~ TAXA.NAME + group_multi, data=zpd_spp_table, FUN=sum)
zpd_spp_table <- zpd_spp_table[zpd_spp_table$group_multi != "",]

setwd("~/Dropbox/my_papers/InPrep_Stepanuketal_GOMprey/")
write.csv(zpd_spp_table, "zoop_taxa.csv")
