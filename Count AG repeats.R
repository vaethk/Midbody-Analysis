library("Biostrings")
library(dplyr)
library(stringr)

fastaFile <- readDNAStringSet("~/3'UTR/UTR.fa")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence)


all_counts <- c()
biggest <- c()
count <- c()

for(i in seq_len(nrow(df))) {
  row <- df[i, 2]
  for (j in seq_len(str_length(row))){
    section <- substr(row, j, (j + 99))
    count <- c(count, countPattern("AG", section))
    biggest <- max(count)
  }
  all_counts <- c(all_counts, biggest)
  biggest <- c()
  count <- c()
}

seperate()

seq_name <- as.data.frame(seq_name)
Names <- tidyr::separate(seq_name, col = seq_name, sep = '_', into = c('ensembl_gene_id_version', 'other'))

Names <- as.data.frame(Names[1])

Names <- tidyr::separate(Names, col = ensembl_gene_id_version, sep = '[.]', into = c('ensembl_gene_id', 'other'))
Names <- as.data.frame(Names[1])

numAG <- data.frame(Names, sequence, all_counts)

write.csv(numAG,"~/3'UTR/AG_repeat_counts.csv", row.names = TRUE)

for(i in seq_len(nrow(numAG))) {
  name <- numAG[i, 1]
  MB <- dat[i,1]
  
}


