full_ngs <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/full_data/ngs_count.csv')
full_meta<- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/full_data/meta_area.csv')

no_24_ngs <- full_ngs[- grep("24", full_ngs$Sample_name),]
no_24_meta <- full_meta[- grep("24", full_meta$Sample_name),]


full_samp <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/full_data/samplelist.csv')
no_24_samp <- full_samp[- grep("24", full_samp$Sample_name),]

P6_ngs <- no_24_ngs[- grep("PAK", no_24_ngs$Sample_name),]
P6_meta <- no_24_meta[- grep("PAK", no_24_meta$Sample_name),]
P6_samp <- no_24_samp[- grep("PAK", no_24_samp$Sample_name),]

write.csv(P6_ngs, '/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/P6_treated_control/p6_ngs.csv', row.names = F)
write.csv(P6_meta, '/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/P6_treated_control/p6_meta.csv', row.names = F)
write.csv(P6_samp, '/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/P6_treated_control/p6_samp.csv', row.names = F)


PAK_ngs <- no_24_ngs[- grep("P6", no_24_ngs$Sample_name),]
PAK_meta <- no_24_meta[- grep("P6", no_24_meta$Sample_name),]
PAK_samp <- no_24_samp[- grep("P6", no_24_samp$Sample_name),]

write.csv(PAK_ngs, '/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/PAK_treated_control/pak_ngs.csv', row.names = F)
write.csv(PAK_meta, '/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/PAK_treated_control/pak_meta.csv', row.names = F)
write.csv(PAK_samp, '/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/PAK_treated_control/pak_samp.csv', row.names = F)

no_pak_treat_samp <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/SUS_plots/samplelist_no_PAK_P.csv')


no_pak_p_ngs <- merge(full_ngs, no_pak_treat_samp, by.x = 1, by.y = 1)
no_pak_p_meta <- merge(full_meta, no_pak_treat_samp, by.x = 1, by.y = 1)

write.csv(no_pak_p_meta, '/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/SUS_plots/no_pak_p_meta.csv', row.names = F)

write.csv(no_pak_p_ngs, '/home/florence/Documents/Bioinformatics_MSc/Project/final_use_csv/SUS_plots/no_pak_p_ngs.csv', row.names = F)


