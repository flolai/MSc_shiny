
count <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/full_data/ngs_count.csv', row.names = 1)

count <- t(count)
count[count==0.001]<- 0

library("DESeq2")

samp_info <- as.matrix(read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/shiny_reactive_MSc/final_use_csv/full_data/samplelist.csv', row.names = 1))

all(rownames(samp_info) %in% colnames(count))         
count <- count[, rownames(samp_info)]

dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = samp_info, 
                              design = ~ Sample_type)

rlog_count <- rlog(dds, blind = TRUE)

rlog_ngs_count <- as.data.frame(assay(rlog_count))

rlog_ngs_count <-t(rlog_ngs_count)


meta <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/shiny_reactive_MSc/final_use_csv/full_data/meta_area.csv', row.names = 1)

meta_log <- log10(meta)
ngs <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/shiny_reactive_MSc/final_use_csv/full_data/rlog_ngs.csv',  row.names = 1)

meta_no_24 <- meta_log[-(grep('24', row.names(meta_log))),]
meta_no_24 <- meta_no_24[-(grep('QC', row.names(meta_no_24))),]
meta_no_24 <- meta_no_24[-(grep('0', row.names(meta_no_24))),]

ngs_no_24 <- ngs[-(grep('24', row.names(ngs))),]
ngs_no_24 <- ngs_no_24[-(grep('0', row.names(ngs_no_24))),]

P6_meta <- meta_no_24[(grep('P6', row.names(meta_no_24))),]
P6_ngs <- ngs_no_24[(grep('P6', row.names(ngs_no_24))),]

PAK_meta <- meta_no_24[(grep('PAK', row.names(meta_no_24))),]

PAK_ngs <- ngs_no_24[(grep('PAK', row.names(ngs_no_24))),]



treated_meta <- meta_no_24[-(grep('C', row.names(meta_no_24))),]
treated_meta <- treated_meta[-(grep('0', row.names(treated_meta))),]
treated_ngs <- ngs_no_24[-(grep('C', row.names(ngs_no_24))),]
treated_ngs <- treated_ngs[-(grep('0', row.names(treated_ngs))),]
