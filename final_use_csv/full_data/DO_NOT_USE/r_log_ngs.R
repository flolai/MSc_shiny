

count <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/full_data/ngs_count.csv', row.names = 1)

count <- t(count)
count[count==0.001]<- 0

library("DESeq2")

samp_info <- as.matrix(read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/full_data/samplelist.csv', row.names = 1))

all(rownames(samp_info) %in% colnames(count))         
count <- count[, rownames(samp_info)]

dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = samp_info, 
                              design = ~ Sample_type)

rlog_count <- rlog(dds, blind = TRUE)

rlog_ngs_count <- as.data.frame(assay(rlog_count))

rlog_ngs_count <-t(rlog_ngs_count)


meta <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/full_data/meta_area.csv', row.names = 1)

meta_log2 <- log2(meta)

write.csv(rlog_ngs_count,'/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/full_data/rlog_ngs.csv' )
write.csv(meta_log2, '/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/full_data/meta_log2.csv')
