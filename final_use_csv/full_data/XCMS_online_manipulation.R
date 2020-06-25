library(tidyr)

neg_comp_name<- read.delim('/home/florence/Documents/Bioinformatics_MSc/Project/raw_file_meta/XCMS_online_neg/results/report_1395910.tsv', header = TRUE, sep = '\t')
PeakID <-sprintf("neg_%d",neg_comp_name$featureidx)
neg_comp_name<-cbind(PeakID,neg_comp_name)
neg_area<- read.delim('/home/florence/Documents/Bioinformatics_MSc/Project/raw_file_meta/XCMS_online_neg/results/XCMS.annotated.Report_1395910.tsv', header = TRUE, sep = '\t')
neg_area<- cbind(PeakID, neg_area)
neg_area<- neg_area[-c(2:12)]
neg_area <- neg_area[ , -which(names(neg_area) %in% c("isotopes","adduct", "pcgroup"))]
rownames(neg_area)<- neg_area$PeakID
neg_area<- t(neg_area)
neg_comp_name <-separate(neg_comp_name, col = METLIN, , into = c("ppm", "Putative.metabolite"), sep = "::")
neg_comp_name <-neg_comp_name[ , which(names(neg_comp_name) %in% c("PeakID","mzmed","rtmed","Putative.metabolite"))]



pos_comp_name<- read.delim('/home/florence/Documents/Bioinformatics_MSc/Project/raw_file_meta/XCMS_online_pos/report_1395904.tsv', header = TRUE, sep = '\t')
PeakID<-sprintf("pos_%d",pos_comp_name$featureidx)
pos_comp_name<-cbind(PeakID,pos_comp_name)
pos_area<- read.delim('/home/florence/Documents/Bioinformatics_MSc/Project/raw_file_meta/XCMS_online_pos/XCMS.annotated.Report_1395904.tsv', header = TRUE, sep = '\t')
pos_area<- cbind(PeakID, pos_area)
pos_area<- pos_area[-c(2:12)]
pos_area <- pos_area[ , -which(names(pos_area) %in% c("isotopes","adduct", "pcgroup"))]
rownames(pos_area)<- pos_area$PeakID
pos_area<- t(pos_area)
pos_comp_name <-separate(pos_comp_name, col = METLIN, , into = c("ppm", "Putative.metabolite"), sep = "::")
pos_comp_name <-pos_comp_name[ , which(names(pos_comp_name) %in% c("PeakID","mzmed","rtmed","Putative.metabolite"))]

meta_area <- cbind(neg_area, pos_area)
meta_comp_name <- rbind(neg_comp_name, pos_comp_name)
names(meta_comp_name)[names(meta_comp_name) == "mzmed"] <- "Mass"
names(meta_comp_name)[names(meta_comp_name) == "rtmed"] <- "RT"


lip<- read.csv("/home/florence/Documents/Bioinformatics_MSc/Project/raw_file_meta/authors_lipid_data.csv", header = TRUE)
rownames(lip)<- lip$PeakID

lip_name <- lip[ , which(names(lip) %in% c("PeakID","Mass", "RT", "Putative.metabolite"))]

lip_area<- t(lip[,-c(1:8)])
tol_area<- merge(meta_area, lip_area, by.x = 0, by.y = 0)
tol_area[tol_area == 0] <- 0.001
names(tol_area)[names(tol_area) == "Row.names"] <- "Sample_name"
comp_name <- rbind(meta_comp_name, lip_name)


write.csv(tol_area, '/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv//meta_area.csv', row.names = FALSE)
write.csv(comp_name, '/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv//meta_comp_name.csv', row.names = FALSE)
