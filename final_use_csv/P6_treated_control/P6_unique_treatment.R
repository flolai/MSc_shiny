
#Reading files for genes and metabolites that are for both P6 and PAK
all_treated <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/shiny_processed_data/treated.csv', sep='\t')
all_control <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/shiny_processed_data/control.csv', sep='\t')

#Reading PAK treated and control data files
PAK_treated <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/shiny_processed_data/PAK_treated.csv')
PAK_control <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/shiny_processed_data/PAK_control.csv', sep='\t')

#Reading genes and metabolites that are different from PAK
P6 <- read.csv('/home/florence/Documents/Bioinformatics_MSc/Project/shiny_processed_data/p6_oplsda_all_time.csv', sep='\t')


library(dplyr)

names(PAK_treated)[1] <- paste("ID")
names(PAK_control)[1] <- paste("ID")
names(all_treated)[1] <- paste("ID")
names(all_control)[1] <- paste("ID")

#finding genes and metabolites that changes in treated and control samples but not unique to PAK
not_in_PAK_treated <- anti_join(all_treated, PAK_treated, by = c("ID" = "ID"))
not_in_PAK_control <- anti_join(all_control, PAK_control, by = c("ID" = "ID"))

names(P6)[1] <- paste("ID")


#finding genes and metabolite that changes in treated and control samples that are unique to P6
p6_treated <- inner_join(P6, not_in_PAK_treated, by = c("ID" = "ID"))
p6_control <- inner_join(P6, not_in_PAK_control, by = c("ID" = "ID"))


write.csv(share_control, '/home/florence/Documents/Bioinformatics_MSc/Project/shiny_processed_data/shared_P6_PAK/share_control.csv', row.names = F)
write.csv(share_treated, '/home/florence/Documents/Bioinformatics_MSc/Project/shiny_processed_data/shared_P6_PAK/share_treated.csv', row.names = F)

