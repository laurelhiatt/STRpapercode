write.table(result, file='/Users/quinlan/Documents/Quinlan-PhD/UDN+STRdb/STRpapercode/gnomad_STR/pathogenicwithmotifandage_results.tsv', quote=FALSE, sep='\t', row.names = NA)


pathogenic = read.csv('/Users/quinlan/Documents/Quinlan-PhD/UDN+STRdb/STRpapercode/gnomad_STR/pathogenic_results.tsv', sep = '\t', stringsAsFactors = FALSE)

pathogenic_motif = read.csv('/Users/quinlan/Documents/Quinlan-PhD/UDN+STRdb/STRpapercode/gnomad_STR/pathogenicwithmotif_results.tsv', sep = '\t', stringsAsFactors = FALSE)

pathogenic_motif_age = read.csv('/Users/quinlan/Documents/Quinlan-PhD/UDN+STRdb/STRpapercode/gnomad_STR/pathogenicwithmotifandage_results.tsv', sep = '\t', stringsAsFactors = FALSE)

mismatched_ids <- vector()

# Loop through each row of the pathogenic dataframe
for (i in 1:nrow(pathogenic)) {
  # Get the corresponding row in the pathogenic_motif dataframe
  matching_row <- pathogenic_motif[pathogenic_motif$Id == pathogenic$Id[i], ]

  # Compare the percent values
  if (pathogenic$percent[i] != matching_row$percent) {
    # If they don't match, add the ID to the mismatched_ids vector
    mismatched_ids <- c(mismatched_ids, pathogenic$Id[i])
  }
}

# Print the mismatched IDs
print(mismatched_ids)


for (i in 1:nrow(pathogenic_motif)) {
  Id_val <- pathogenic_motif$Id[i]
  percent_val <- pathogenic_motif$percent[i]
  percent_val_motif <- pathogenic_motif_age$percent[pathogenic_motif_age$Id == Id_val]

  if (length(percent_val_motif) > 0) {
    percent_val_motif <- percent_val_motif[1]

    if (percent_val != percent_val_motif) {
      cat("Id:", Id_val, "has a percent difference of", percent_val_motif - percent_val, "\n")
    }
  } else {
    cat("Id:", Id_val, "not found in pathogenic_motif data frame\n")
  }
}
