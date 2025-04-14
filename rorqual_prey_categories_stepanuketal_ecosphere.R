data$group <- ifelse(data$SVSPP %in% c(32, 851, 205, 859, 427, 430, 30, 31, 44),  "Clupeid", NA)
data$group <- ifelse(data$SVSPP %in% c(121, 578, 572, 209, 124, 898, 570, 702, 571, 582, 744, 860, 208, 212, 211, 874, 822, 745, 577, 877), "Scombrid", data$group)
data$group <- ifelse(data$SVSPP %in% c(181, 734, 38), "Sandlance", data$group)
data$group <- ifelse((data$SVSPP >= 501 & data$SVSPP <= 505) | (data$SVSPP >= 510 & data$SVSPP <= 513), "Squid_BT", data$group)
data$group <- ifelse(data$SVSPP == 297 | data$SVSPP == 298 | data$SVSPP == 306 | data$SVSPP == 307, "Shrimp_BT", data$group)
data$group <- ifelse(data$SVSPP == 306, "Northern_shrimp", data$group)
data$group <- ifelse(data$SVSPP %in% c(73:75), "Large_Gadid", data$group)
data$group <- ifelse(data$SVSPP %in% c(189, 1, 262, 192, 131, 750, 176, 168, 249, 193, 155, 164), "Misc_Fish", data$group)
data$group <- ifelse(data$SVSPP %in% c(81, 454, 80, 86, 79, 69, 77, 72, 455, 78, 76), "Small_Gadid", data$group)
# data$group <- ifelse(data$SVSPP %in% c(790, 792, 110, 793, 777, 100, 104, 785, 786, 788, 109, 795, 775, 776, 773, 791, 787, 117, 789, 783, 103, 853, 774, 873, 106, 107, 105), "Flatfish", data$group)

data$bp_prey <- ifelse(data$group == "Shrimp_BT" | data$group == "Clupeid" |data$group == "Sandlance" , 1, 0) # also zoop
data$mn_prey <- ifelse(data$group == "Shrimp_BT" | data$group == "Misc_Fish"| data$group == "Clupeid" |data$group == "Sandlance", 1, 0)
data$ba_prey <- ifelse(data$group == "Shrimp_BT" | data$group == "Clupeid" |data$group == "Sandlance" | data$group == "Large_Gadid", 1, 0)
data$bb_prey <- ifelse(data$group == "Shrimp_BT" | data$group == "Clupeid" | data$group == "Squid_BT", 1, 0) #also zoop