library(rtracklayer)
library(dplyr)

ngs<- read.csv("/home/florence/Documents/Bioinformatics_MSc/Project/raw_ngs/counts.txt", sep = "", row.names = 1)
ngs <- ngs[ , -which(names(ngs) %in% c("Length"))]
names(ngs) <- gsub(pattern = "X.d.projects.u.lc001.aligned.", replacement = "", x = names(ngs))
names(ngs) <- gsub(pattern = ".fastq.gz_trimmed.fastq.gz.subread.BAM", replacement = "", x = names(ngs))
ngs_sample<- read.csv("/home/florence/Documents/Bioinformatics_MSc/Project/raw_ngs/RNA_seq_sample_name.txt", sep = "\t")
Sample_name<- as.data.frame(gsub("Illumina HiSeq 1500 sequencing; RNA-Seq of P. aeruginosa: ", "\\1", ngs_sample$experiment_title))
colnames(Sample_name)<- c("Sample_name")
Sample_name<- as.data.frame(gsub("AKpmrB", "\\1", Sample_name$Sample_name))
colnames(Sample_name)<- c("Sample_name")
Sample_name<- as.data.frame(gsub("h", "\\1", Sample_name$Sample_name))
colnames(Sample_name)<- c("Sample_name")

name<- Sample_name %>% separate(Sample_name, c("a", "b","c", "d"))
Sample_name<- unite(name, "samp", c("a", "d", "b", "c"), sep = ".")
Sample_name<- as.data.frame(gsub("NA.", "\\1", Sample_name$samp))
colnames(Sample_name)<- c("Sample_name")
ngs_sample <-cbind(ngs_sample, Sample_name)

#export ngs_sample for minor manual fix of sample name inconsistancy
write.csv(ngs_sample, '/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/ngs_sample_temp.csv', row.names = FALSE)
ngs_sample <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/ngs_sample_temp.csv')
ngs<- t(ngs)
ngs<- merge(ngs_sample, ngs, by.x = "run_accession", by.y = 0)
ngs <- ngs[ , -which(names(ngs) %in% c("run_accession"))]
ngs[ngs == 0] <- 0.001


write.csv(ngs, '/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/ngs_count.csv', row.names = FALSE)

gff<-readGFF("/home/florence/Documents/Bioinformatics_MSc/Project/raw_ngs/PAK_annotation.gff3", version=3)
ngs_ID<-as.data.frame(gff)
ngs_ID <- ngs_ID[ , which(names(ngs_ID) %in% c("ID", "Name", "gene"))]

write.csv(ngs_ID, '/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/ngs_gene_name.csv', row.names = FALSE)
