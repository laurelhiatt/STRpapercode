


### clean up ucsc
# awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) print i, $i }' ucscrepeatmaskerhg38.bed
# awk -F'\t' 'BEGIN {OFS=FS} {print $6, $7, $8, $0}' ucscrepeatmaskerhg38 > ucscrepeatmaskerhg38.bed
#awk -F'\t' -v values="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY" 'BEGIN { split(values, arr, " ") } $1 ~ "^("arr[1]"|"arr[2]"|"arr[3]"|"arr[4]"|"arr[5]"|"arr[6]"|"arr[7]"|"arr[8]"|"arr[9]"|"arr[10]"|"arr[11]"|"arr[12]"|"arr[13]"|"arr[14]"|"arr[15]"|"arr[16]"|"arr[17]"|"arr[18]"|"arr[19]"|"arr[20]"|"arr[21]"|"arr[22]"|"arr[23]"|"arr[24]")$" { print }' ucscsimplerepeatshg38.bed > ucscsimplerepeatshg38_filtered.bed
#awk -F'\t' 'BEGIN {OFS=FS} {print $2, $3, $4, $0}' ucscsimplerepeatshg38 > ucscsimplerepeatshg38.bed


### clean up gnomadSTR
awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) print i, $i }' gnomAD_STR_genotypes__2022_01_20.tsv
awk -F'\t' 'BEGIN {OFS=FS} {print $4, $5, $6, $0}' gnomAD_STR_genotypes__2022_01_20.tsv > gnomAD_STR_genotypes__2022_01_20.bed
tail -n +2 gnomAD_STR_genotypes__2022_01_20.bed > gnomAD_STR_genotypes__2022_01_20_noheader.bed



bedtools intersect -a /Users/quinlan/Documents/Quinlan-PhD/UDN+STRdb/ucscrepeatmaskerhg38_filtered.bed -b /Users/quinlan/Documents/Quinlan-PhD/UDN+STRdb/gnomAD_STR_genotypes__2022_01_20.bed -header -wb > intersectionofgnomADandrepeats.txt

# uniq intersectionofgnomADandrepeats.txt > intersectionofgnomADandrepeats_uniq.txt


sort -k1,1 -k2,2n ucscrepeatmaskerhg38_filtered.bed > ucscrepeatmaskerhg38_filtered_sorted.bed
sort -k1,1 -k2,2n gnomAD_STR_genotypes__2022_01_20_noheader.bed > gnomAD_STR_genotypes__2022_01_20_noheader_sorted.bed
cut -f 1-6 gnomAD_STR_genotypes__2022_01_20_noheader_sorted.bed > gnomAD_STR_genotypes__2022_01_20_noheader_sorted_small.bed

bedtools closest -a /Users/quinlan/Documents/Quinlan-PhD/UDN+STRdb/gnomAD_STR_genotypes__2022_01_20_noheader_sorted_extrasmall.bed -b /Users/quinlan/Documents/Quinlan-PhD/UDN+STRdb/ucscrepeatmaskerhg38_filtered_sorted.bed -d > closestgnomADandrepeats.txt

#> windowofgnomADandrepeats.txt



awk '{ $2 = $2 - 500; $3 = $3 + 500; print}' /Users/quinlan/Documents/Quinlan-PhD/UDN+STRdb/gnomAD_STR_genotypes__2022_01_20_noheader_sorted_extrasmall.bed > /Users/quinlan/Documents/Quinlan-PhD/UDN+STRdb/gnomAD_STR_genotypes__2022_01_20_noheader_sorted_extrasmall_wiggle.bed
awk '{ $2 = $2 - 500; $3 = $3 + 500; print}' /Users/quinlan/Documents/Quinlan-PhD/UDN+STRdb/ucscrepeatmaskerhg38_filtered_sorted.bed > /Users/quinlan/Documents/Quinlan-PhD/UDN+STRdb/ucscrepeatmaskerhg38_filtered_sorted_wiggle.bed
awk '{ $2 = $2 - 500; $3 = $3 + 500; print}' /Users/quinlan/Documents/Quinlan-PhD/UDN+STRdb/ucscsimplerepeatshg38_filtered_sorted.bed > /Users/quinlan/Documents/Quinlan-PhD/UDN+STRdb/ucscsimplerepeatshg38_filtered_sorted_wiggle.bed


bedtools intersect -a /Users/quinlan/Documents/Quinlan-PhD/UDN+STRdb/gnomAD_STR_genotypes__2022_01_20_noheader_sorted_extrasmall_wiggle.bed -b ucscrepeatmaskerhg38_filtered_sorted_wiggle.bed -header -wb > intersectionofgnomADandrepeats_wiggle.txt

bedtools intersect -a /Users/quinlan/Documents/Quinlan-PhD/UDN+STRdb/gnomAD_STR_genotypes__2022_01_20_noheader_sorted_extrasmall_wiggle.bed -b ucscsimplerepeatshg38_filtered_sorted_wiggle.bed -header -wb > intersectionofgnomADandSIMPLErepeats_wiggle.txt

# perl -p -i -e 's/ /\t/g' ucscrepeatmaskerhg38_filtered_sorted_wiggle.bed
# perl -p -i -e 's/ /\t/g' gnomAD_STR_genotypes__2022_01_20_noheader_sorted_extrasmall_wiggle.bed
perl -p -i -e 's/ /\t/g' ucscsimplerepeatshg38_filtered_sorted_wiggle.bed