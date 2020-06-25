ngs <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/full_data/ngs_count.csv')
meta <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/full_data/meta_area.csv')
samp <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/full_data/samplelist.csv')

all <- merge(samp, meta, by.x = 1, by.y =1)
all <- merge (all, ngs,by.x = 1, by.y =1)
all<- all[order(all$Sample_type),]
row.names(all)<- all$Sample_name
all <- all[,-c(1:6)]
all <- t(all)

wilcox.test(all[1, 1:21], all[1, 22:42])

wilcoxtest <- function(df, grp1, grp2) {
  x <- df[grp1]
  y <- df[grp2]
  x <- as.numeric(x)
  y <- as.numeric(y)  
  pvalues = wilcox.test(x, y)
  pvalues$p.value
}

allpvalue <- as.data.frame(apply(all, 1, wilcoxtest, grp1 = c(1:21), grp2 = c(22:42)))
log2_all <- log2(all)
p6_FC <- apply(log2_all[,1:21], 1, median)
pak_FC <- apply(log2_all[,22:42], 1, median)

log2_fc <- as.data.frame(p6_FC - pak_FC)
pvalues_fc <- merge(allpvalue, log2_fc, by.x = 0, by.y = 0)

p6 <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/shiny_processed_data/P6_PAK/p6_oplsda_all_time.csv', sep = '\t')
pak <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/shiny_processed_data/P6_PAK/pak_oplsda_all_time.csv', sep = '\t')
p6_FC <- merge(p6, pvalues_fc, by.x = 1, by.y = 1)
pak_FC <- merge(pak, pvalues_fc, by.x = 1, by.y = 1)

write.csv(p6_FC, '/home/florence/Documents/Bioinformatics_MSc/Project/shiny_processed_data/P6_PAK/p6_FC.csv' )
write.csv(pak_FC, '/home/florence/Documents/Bioinformatics_MSc/Project/shiny_processed_data/P6_PAK/pak_FC.csv')
