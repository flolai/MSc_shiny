ngs <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/shiny_reactive_MSc/final_use_csv/full_data/DO_NOT_USE/ngs_count.csv')
meta <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/shiny_reactive_MSc/final_use_csv/full_data/DO_NOT_USE/meta_area.csv')
samp <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/shiny_reactive_MSc/final_use_csv/full_data/samplelist.csv')

all <- merge(samp, meta, by.x = 1, by.y =1)
all <- merge (all, ngs,by.x = 1, by.y =1)
all<- all[order(all$Sample_type),]
row.names(all)<- all$Sample_name
all <- all[,-c(1:6)]
all <- t(all)

#wilcox.test(all[1, 1:21], all[1, 22:42])

wilcoxtest <- function(df, grp1, grp2) {
  x <- df[grp1]
  y <- df[grp2]
  x <- as.numeric(x)
  y <- as.numeric(y)  
  pvalues = wilcox.test(x, y)
  pvalues$p.value
}

# P6 vs PAK
allpvalue <- as.data.frame(apply(all, 1, wilcoxtest, grp1 = c(1:21), grp2 = c(22:42)))
log2_all <- log2(all)
p6_FC <- apply(log2_all[,1:21], 1, median)
pak_FC <- apply(log2_all[,22:42], 1, median)

log2_fc <- as.data.frame(p6_FC - pak_FC)
pvalues_fc <- merge(allpvalue, log2_fc, by.x = 0, by.y = 0)

p6 <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/shiny_reactive_MSc/shiny_processed_data/P6_vs_PAK/P6_all_time.csv')
pak <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/shiny_reactive_MSc/shiny_processed_data/P6_vs_PAK/PAK_all_time.csv')
p6_FC <- merge(p6, pvalues_fc, by.x = 1, by.y = 1)
pak_FC <- merge(pak, pvalues_fc, by.x = 1, by.y = 1)

write.csv(p6_FC, '/home/florence/Documents/Bioinformatics_MSc/Project/shiny_reactive_MSc/shiny_processed_data/P6_vs_PAK/P6_FC.csv' )
write.csv(pak_FC, '/home/florence/Documents/Bioinformatics_MSc/Project/shiny_reactive_MSc/shiny_processed_data/P6_vs_PAK/pak_FC.csv')


#treated and control

control <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/shiny_reactive_MSc/shiny_processed_data/treated_control_P6_PAK/control.csv')
treated <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/shiny_reactive_MSc/shiny_processed_data/treated_control_P6_PAK/treated.csv')
no_24_samp_list <- all[-(grep('0|24', all$Time_point)),]
no_24_samp_list<- no_24_samp_list[order(no_24_samp_list$Treated_Control),]
row.names(no_24_samp_list)<- no_24_samp_list$Sample_name
no_24_samp_list <- no_24_samp_list[,-c(1:6)]
no_24_samp_list <- t(no_24_samp_list)
t_c_pvalue <- as.data.frame(apply(no_24_samp_list, 1, wilcoxtest, grp1 = c(1:12), grp2 = c(13:24)))
log2_t_c <- log2(no_24_samp_list)

control_FC <- apply(log2_t_c[,1:12], 1, median)
treated_FC <- apply(log2_t_c[,12:24], 1, median)

log2_fc_t_c <- as.data.frame(control_FC - treated_FC)
pvalues_fc_t_C <- merge(t_c_pvalue, log2_fc_t_c, by.x = 0, by.y = 0)

control_FC <- merge(control, pvalues_fc_t_C, by.x = 1, by.y = 1)
treated_FC <- merge(treated, pvalues_fc_t_C, by.x = 1, by.y = 1)

write.csv(control_FC, '/home/florence/Documents/Bioinformatics_MSc/Project/shiny_reactive_MSc/shiny_processed_data/treated_control_P6_PAK/control_FC.csv')
write.csv(treated_FC, '/home/florence/Documents/Bioinformatics_MSc/Project/shiny_reactive_MSc/shiny_processed_data/treated_control_P6_PAK/treated_FC.csv')

#PAK treated control

PAK_treated <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/shiny_reactive_MSc/shiny_processed_data/PAK_control_treated/PAK_treated.csv')
PAK_control <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/shiny_reactive_MSc/shiny_processed_data/PAK_control_treated/PAK_control.csv')
PAK_samp_list <- all[(grep('PAK', all$Sample_type)),]
PAK_samp_list <- PAK_samp_list[-(grep('0|24', PAK_samp_list$Time_point)),]
PAK_samp_list<- PAK_samp_list[order(PAK_samp_list$Treated_Control),]
row.names(PAK_samp_list)<- PAK_samp_list$Sample_name
PAK_samp_list <- PAK_samp_list[,-c(1:6)]
PAK_samp_list <- t(PAK_samp_list)
#p-value calculation
pak_pvalue <- as.data.frame(apply(PAK_samp_list, 1, wilcoxtest, grp1 = c(1:6), grp2 = c(7:12)))
log2_pak <- log2(PAK_samp_list)

pak_control_FC <- apply(log2_pak[,1:6], 1, median)
pak_treated_FC <- apply(log2_pak[,7:12], 1, median)

log2_pak_t_c <- as.data.frame(pak_control_FC - pak_treated_FC)
pvalues_log2_pak <- merge(pak_pvalue, log2_pak_t_c, by.x = 0, by.y = 0)

pak_control_FC <- merge(PAK_control, pvalues_log2_pak, by.x = 1, by.y = 1)
pak_treated_FC <- merge(PAK_treated, pvalues_log2_pak, by.x = 1, by.y = 1)

write.csv(pak_treated_FC, '/home/florence/Documents/Bioinformatics_MSc/Project/shiny_reactive_MSc/shiny_processed_data/PAK_control_treated/pak_treated_FC.csv')
write.csv(pak_control_FC , '/home/florence/Documents/Bioinformatics_MSc/Project/shiny_reactive_MSc/shiny_processed_data/PAK_control_treated/pak_control_FC.csv')