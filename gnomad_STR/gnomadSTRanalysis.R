library('dplyr')
library('stringr')
library('ggplot2')
library('cowplot')
theme_set(theme_cowplot())
options(stringsAsFactors = FALSE)

# Get data from here: https://gnomad.broadinstitute.org/downloads#v3-short-tandem-repeats
gnomADSTRcalls = read.csv('/Users/quinlan/Documents/Quinlan-PhD/UDN+STRdb/gnomad_STR/gnomAD_STR_genotypes__2022_01_20.tsv.gz', sep = '\t')

STR_table <- read.csv('/Users/quinlan/Documents/Quinlan-PhD/UDN+STRdb/gnomad_STR/STR_table_04072023.csv')


#cleaning up table to match gnomad syntax
STR_table$gene[STR_table$gene == "AFF2/FMR2"] <- "AFF2"
STR_table$gene[STR_table$gene == "MARCH6"] <- "MARCHF6"
STR_table$gene[STR_table$gene == "C9orf72"] <- "C9ORF72"
STR_table$gene[STR_table$gene == "ATXN8OS/ATXN8"] <- "ATXN8OS"
STR_table$gene[STR_table$gene == "NUTM2B-AS1/NUTM2B-AS1"] <- "NUTM2B-AS1"
STR_table$gene[STR_table$gene == "ARX" & STR_table$stop_hg38 == 25013697] <- "ARX_1"
STR_table$gene[STR_table$gene == "ARX" & STR_table$stop_hg38 == 25013565] <- "ARX_2"
STR_table$gene[STR_table$gene == "HOXA13" & STR_table$stop_hg38 == 27199966] <- "HOXA13_1"
STR_table$gene[STR_table$gene == "HOXA13" & STR_table$stop_hg38 == 27199861] <- "HOXA13_2"
STR_table$gene[STR_table$gene == "HOXA13" & STR_table$stop_hg38 == 27199732] <- "HOXA13_3"

unique(gnomADSTRcalls$Id) -> genelist
unique(STR_table$gene) -> tablegenelist
sort(tablegenelist) -> tablegenelist


# what are the differences between the table and the gnomad gene list
tablegenelist[!(tablegenelist %in% genelist)]

genelist[!(genelist %in% tablegenelist)]

# what are the same genesbetween the table and the gnomad gene list
tablegenelist[(tablegenelist %in% genelist)]

STR_table_clean <-subset(STR_table, select=c("disease_id", "gene",
                                             "repeatunit_ref", "Inheritance",
                                             "type", "normal_min", "normal_max",
                                             "intermediate_min", "intermediate_max",
                                             "pathogenic_min", "pathogenic_max",
                                             "repeat_unit_lenth", "age_onset_min",
                                             "age_onset_max", "novel", "repeatunit_path_norm"))

STR_table_clean <- STR_table_clean %>%
  rename("Id" = "gene")


gnomADSTRcalls$AgeMax = as.numeric(str_sub(gnomADSTRcalls$Age,-2,-1)) # get last two char in Age. Assumes all <100

total <- merge(STR_table_clean,gnomADSTRcalls,by="Id")


total <- sample_n(total, 10000)

# Subset to a single locus to start with (should loop through all loci later)
#gnomADSTRcalls = subset(gnomADSTRcalls, Id == "ATXN8OS")

# Make a numerical age to simplify plotting


# Check Allele 2 is always the largest
table(gnomADSTRcalls$Allele2 >= gnomADSTRcalls$Allele1)

ggplot(gnomADSTRcalls) +
  geom_histogram(aes(x = Allele1, fill = 'Allele1'), alpha = 0.5) +
  geom_histogram(aes(x = Allele2, fill = 'Allele2'), alpha = 0.5)

ggplot(subset(gnomADSTRcalls, Age != 'age_not_available')) +
  geom_histogram(aes(x = Allele2, fill = Age), alpha = 0.5) + facet_wrap(~AgeMax > 50)

ggplot(subset(gnomADSTRcalls, Population == 'nfe')) +
  geom_histogram(aes(x = Allele2, fill = Age), alpha = 0.5) + facet_wrap(~AgeMax > 50, scales = 'free_y') +
  scale_y_log10()

ggplot(gnomADSTRcalls) +
  geom_histogram(aes(x = Allele2, fill = PcrProtocol), alpha = 0.5) + facet_wrap(~AgeMax > 50, scales = 'free_y')

wilcox.test(subset(gnomADSTRcalls, Age != 'age_not_available' & Population == 'nfe')$Allele2,

            subset(gnomADSTRcalls, Age == 'age_not_available' & Population == 'nfe')$Allele2)

summary(glm(AgeMax ~ Allele2, data = subset(gnomADSTRcalls, Age != 'age_not_available')))

ggplot(subset(gnomADSTRcalls, Age != 'age_not_available'), aes(x = AgeMax, y = Allele2)) + geom_point()

table(subset(gnomADSTRcalls, Age != 'age_not_available')$Population)

ggplot(subset(total, Age != 'age_not_available'), aes(x = type, y = Allele2, col=Id)) + geom_point()

ggplot(total) + geom_point(aes(x = Id, y = Allele2UsingOfftargetRegions, col=type)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(~type, scales = "free_x", space = "free_x")

ggplot(total) + geom_point(aes(x = Id, y = Allele2UsingOfftargetRegions, col=type)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_grid(~type, scales = "free_x") + ylim(0,4000)

ggplot(subset(total, Age != 'age_not_available'), aes(x = type, y = Allele1)) + geom_point()


ggplot(subset(total, Age != 'age_not_available'), aes(x = AgeMax, y = age_onset_min, col=Id)) + geom_point()
ggplot(subset(total), aes(x = AgeMax, y = age_onset_min, col=Id)) + geom_point()


ggplot(total) + geom_point(aes(x = age_onset_max, y = Allele2UsingOfftargetRegions, col=Id))

ggplot(total, aes(x = age_onset_max, y = Allele2, col=Id)) + geom_point()

