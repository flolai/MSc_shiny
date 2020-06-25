test_p6_ngs <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/P6_treated_control/p6_ngs.csv')
test_p6_meta <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/P6_treated_control/p6_meta.csv')
test_p6_samp <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/P6_treated_control/p6_samp.csv')


p6_no_t0_ngs <- test_p6_ngs[- grep("0", test_p6_ngs$Sample_name),]
p6_no_t0_meta <- test_p6_meta[- grep("0", test_p6_meta$Sample_name),]
p6_no_t0_samp <- test_p6_samp[- grep("0", test_p6_samp$Sample_name),]

write.csv(p6_no_t0_ngs, '/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/P6_treated_control/p6_no_t0_ngs.csv', row.names = F)
write.csv(p6_no_t0_meta, '/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/P6_treated_control/p6_no_t0_meta.csv', row.names = F)
write.csv(p6_no_t0_samp, '/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/P6_treated_control/p6_no_t0_samp.csv', row.names = F)

class <- as.character(test_p6_samp$Treated_Control)
p6_meta <- test_p6_meta[- grep("5", test_p6_meta$Sample_name),]
p6_meta <- p6_meta[- grep("QC", p6_meta$Sample_name),]
p6_meta <- p6_meta[- grep("4.4", p6_meta$Sample_name),]
p6_meta <- p6_meta[- grep("4.1", p6_meta$Sample_name),]
p6_meta <- p6_meta[- grep("4.0", p6_meta$Sample_name),]
row.names(p6_meta)<- p6_meta[,c(1)]

p6_meta <- log2(p6_meta[,-c(1)])
p6_meta <- scale(p6_meta, center = TRUE, scale = TRUE)
p6_meta[is.na(p6_meta)] <- 0.001


test_model <- opls(p6_meta , scaleC= "none")

test_oplsda <- opls(p6_meta , class, scaleC= "none", predI = 1, orthoI = NA ,permI = 100)
