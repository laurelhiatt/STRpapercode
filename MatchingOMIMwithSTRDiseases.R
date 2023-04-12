setwd("/Users/quinlan/Documents/Quinlan-PhD/UDN+STRdb") ###or wherever you want to be
list.files()

data <- read.csv("KnownSTRdOMIMwithPhenotype.txt", stringsAsFactors = FALSE, header = FALSE, sep='\t')
data2 <- read.csv("STR disease loci - Disease loci - 2022.csv", stringsAsFactors = TRUE, header = TRUE)

names(data)[1] <- 'gene_id'

names(data)[9] <- 'OMIM.ID'

merged <- merge(data, data2[, c("OMIM.ID", "disease", "disease_id")], by="OMIM.ID")

write.table(merged, 'mergedSTRdiseaseswithphenotypes.tsv', row.names = FALSE, quote=FALSE, sep='\t')