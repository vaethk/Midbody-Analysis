---
title: "Counting AG in 3'UTR"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(biomaRt)
library(Biostrings)
library(stringi)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggpubr)

c("biomaRt", "Biostrings", "stringi", "stringr", "ggplot2", "ggpubr")
```

```{r}
mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host='www.ensembl.org')
t2g <- biomaRt::getBM(attributes = c('ensembl_transcript_id', 'ensembl_gene_id', 'external_gene_name', 'refseq_mrna'), mart = mart)
t2g <- dplyr::rename(t2g, target_id= ensembl_transcript_id, ext_gene = external_gene_name)
ens2gene <- t2g[,c(2,3)]
colnames(ens2gene)[2] <- 'Gene'
ens2gene <- unique(ens2gene)
 here()
transcripts <- read.csv("transcript_id.csv")
write.csv(t2g$target_id, "transcript_id.csv")
```

```{r, Get 3'UTRs}
#use this if you need to get the 3'UTR sequences 

Localized_genes <- data.frame(Gene = c("AFAP1L1", "OGT", "PAXBP1","FNBP4" , "TRAK2", "RAB13", "CPLX2", "NET1", "CDC42BPG", "GDF11","AKAP12", "TRP53INP2", "KIF5B")) %>% left_join(., ens2gene, by = "Gene") #my sequences of interest

THREE_UTR <- getBM(attributes = c('ensembl_gene_id',
                                'ensembl_transcript_id',
                                '3utr'),
                 filters = 'ensembl_gene_id',
                 values = list(Localized_genes$ensembl_gene_id), #list of genes, all genes, or one gene
                 mart = mart)

THREE_UTR <- THREE_UTR %>%
  as.data.frame() %>% 
  filter(THREE_UTR$`3utr` != "Sequence unavailable") #only want isoforms that have a 3'UTR sequence included

THREE_UTR <- THREE_UTR %>% 
  mutate(length = nchar(THREE_UTR$`3utr`)) %>% 
  group_by(., ensembl_gene_id) %>% 
  filter(length == max(length)) #keep only the longest isoform 
```

```{r, Already have the 3'UTR sequences}
#use this if you have the 3'UTR sequences-- make a df! 

fastaFile <- readDNAStringSet("~/3'UTR/UTR.fa")
seq_name = names(fastaFile)
  seq_name <- as.data.frame(seq_name)
  seq_name <- tidyr::separate(seq_name, col = seq_name, sep = '[.]', into = c('Gene', 'other')) 
  seq_name <- seq_name %>% select(Gene)
sequence = paste(fastaFile)

df <- data.frame(seq_name, sequence) %>% 
  mutate(seq_length = apply(df, 2, nchar)[,2]) #the names were annoying in my file 

```

```{r, count A/G - if your df have >1 3'UTR sequence of interest}
#Quantifying max A/G count in 100bp chunks - this ONLY saves the max AG
#df has all 3'UTRs so it take a little to run but not forever 

all_counts <- c()
biggest <- c()

for(i in seq_len(nrow(df))) {
  row <- df[i, 2]
  for (j in seq_len(str_length(row))){
    section <- substr(row, j, (j + 99))
    A_G <- stringi::stri_count_fixed(section, c("A", "G"))
    biggest <- max (c(biggest, sum(A_G)))
  }
  all_counts <- c(all_counts, biggest)
  biggest <- c()
  count <- c()
}

df$AG_content <- all_counts

df <- df %>% 
  mutate(AG_rich = ifelse(AG_content > 60, "AG_rich", "Not"))

df_noseq <- df %>% select( -sequence, -'ifelse(AG_content > 70, "AG_rich", "Not")')
```

```{r, this is code that stores the A/G count in the windows across the 3'UTR }
#THREE_UTR is a smaller df that I used to look at  a few genes/make sure my code worked. it might take forever to run on all the genes 
for (j in 1:nrow(THREE_UTR)) {
  count <- c()
  biggest <- c()
  for (i in 1:str_length(THREE_UTR$`3utr`[j])) {
    section <- substr(THREE_UTR$`3utr`[j], i, (i + 99)) 
    A_G <- stringi::stri_count_fixed(section, c("A", "G"))
    count <- c(count, sum(A_G))
    biggest <- max(c(biggest, sum(A_G)))
  }
  gene <- as.character(Localized_genes$Gene[j])
  count <- as.data.frame(count)
  count$start_window <- 1:nrow(count)
  assign(gene, count)
}

#the annoyng thing about this code is that it stores each individal gene as its own df :/ 
```


```{r}
Localized_genes <- res_df_sig_MB %>% left_join(., df, by = "Gene")
Not_Localized_genes <- res_df_sig_soma %>% left_join(., df, by = "Gene")
colnames(df)[1] <- "ensembl_gene_id"

AG_plot_data <- df %>% 
  left_join(., res_df, by = "ensembl_gene_id") %>% 
  filter(., padj < 0.05) %>% 
  mutate(Localized = ifelse(log2FoldChange > 0, "Midbody", "Not Localized")) %>% 
  mutate(window_cutoff = ifelse(AG_content > 70, "Contains High A/G Window", "Does Not Contain High A/G Window"))

ggplot(AG_plot_data, aes(Localized, AG_content, color = Localized)) + geom_violin() + 
  theme_classic(16) + 
  xlab(element_blank()) + 
  ylab("Max A/G Content in 100nt Window") + 
  scale_color_brewer(palette = "Paired") + 
  theme(legend.position = "none")

ggplot(AG_plot_data, aes(window_cutoff, log2FoldChange)) + 
  geom_boxplot(aes(color = window_cutoff)) + 
  theme_classic(16) + 
  xlab(element_blank()) + 
  ylab("Midbody / Whole Cell, Log2") + 
  scale_color_brewer(palette = "Paired") + 
  theme(legend.position = "none")
```

```{r}
FIVE_UTR <- getBM(attributes = c('ensembl_gene_id',
                                'ensembl_transcript_id',
                                '5utr'),
                 filters = 'ensembl_gene_id',
                 values = list(transcripts$x),
                 mart = ensembl)

FIVE_UTR <- FIVE_UTR %>%
  as.data.frame() %>% 
  filter(FIVE_UTR$`5utr` != "Sequence unavailable") #only want isoforms that have a 3'UTR sequence included

FIVE_UTR <- FIVE_UTR %>% 
  mutate(length = nchar(FIVE_UTR$`5utr`)) %>% 
  group_by(., ensembl_gene_id) %>% 
  filter(length == max(length)) #keep only the longest isoform 

all_counts_CT <- c()
 
for(i in 1:nrow(FIVE_UTR)) {
  row <- FIVE_UTR[i, 2]
  for (j in 1:str_length(row)){
    section <- substr(row, j, (j + 14))
    C_T <- stringi::stri_count_fixed(section, c("C", "T"))
    biggest <- max (c(biggest, sum(C_T)))
  }
  all_counts_CT <- c(all_counts_CT, biggest)
  biggest <- c()
} 
  
  
all_counts_CT <- data.frame(all_counts_CT) 
all_counts_CT$ensembl_transcript_id <- FIVE_UTR$ensembl_transcript_id
all_counts_CT <- all_counts_CT %>% mutate(threshold = ifelse(all_counts_CT > 12, "High", "Low"))


df_CT <- all_counts_CT %>% 
  left_join(., FIVE_UTR, by = "ensembl_transcript_id") %>% 
  left_join(., res_df, by = "ensembl_gene_id") %>%
  filter(FIVE_UTR$`5utr` != "Sequence unavailable") %>% 
  mutate(Localized = ifelse(log2FoldChange > 0, "Midbody", "Not Localized")) %>% 
  mutate(window_cutoff = ifelse(all_counts_CT > 12, "Contains C/T rich window", "Does Not Contain C/G Rich Window"))

ggplot(df_CT, aes(Localized, all_counts_CT, color = Localized)) + geom_violin() + 
  theme_classic(16) + 
  xlab(element_blank()) + 
  ylab("Max C/T Content in 25nt Window") + 
  scale_color_brewer(palette = "Paired") + 
  theme(legend.position = "none") + 
  stat_compare_means( method = "t.test", ref.group = "Not Localized")

ggplot(df_CT, aes(window_cutoff, log2FoldChange)) + 
  geom_boxplot(aes(color = window_cutoff)) + 
  theme_classic(16) + 
  xlab(element_blank()) + 
  ylab("Midbody / Whole Cell, Log2") + 
  scale_color_brewer(palette = "Paired") + 
  theme(legend.position = "none") + 
  stat_compare_means( method = "t.test", ref.group = "Does Not Contain C/G Rich Window")



```


```{r}
PRTE <- "^(C|T)(T|C)(T|C)(C|T)(C|T|A)T(T|C)(T|C)(C|T)"

all_counts_PRTE <- c()
 
for(i in 1:nrow(FIVE_UTR)) {
  row <- FIVE_UTR[i, 2]
  count <- str_count(row, PRTE)
  all_counts_PRTE <- c(all_counts_PRTE, count)
} 


all_counts_PRTE <- data.frame(all_counts_PRTE) 
all_counts_PRTE$ensembl_gene_id <- FIVE_UTR$ensembl_gene_id

PRTE_plot <- all_counts_PRTE %>% left_join(., res_df, by = "ensembl_gene_id") %>% 
  mutate(PRTE_yes = ifelse(all_counts_PRTE > 0 , "Contains PRTE", "Does Not Contain PRTE")) %>% 
  mutate(localized = ifelse(log2FoldChange > 0, "MB", "WC"))


ggplot(PRTE_plot, aes(PRTE_yes, log2FoldChange)) + 
  geom_boxplot(aes(color = PRTE_yes)) + 
  theme_classic(16) + 
  xlab(element_blank()) + 
  ylab("Midbody / Whole Cell, Log2") + 
  scale_color_brewer(palette = "Paired") + 
  theme(legend.position = "none") + 
  stat_compare_means( method = "t.test", ref.group = "Does Not Contain PRTE")

ggplot(PRTE_plot, aes(localized, all_counts_RTE)) + 
  geom_violin(aes(color = localized)) + 
  theme_classic(16) + 
  xlab(element_blank()) + 
  ylab("Midbody / Whole Cell, Log2") + 
  scale_color_brewer(palette = "Paired") + 
  theme(legend.position = "none") + 
  stat_compare_means( method = "t.test", ref.group = "WC")

```
```{r}
AB_PRTE <- all_counts_PRTE %>% left_join(., AB_data, by = "ensembl_gene_id") %>% 
  mutate(PRTE_yes = ifelse(all_counts_PRTE > 0 , "Contains PRTE", "Does Not Contain PRTE")) %>% 
  mutate(localized = ifelse(Log2FoldChange_AB > 0, "Apical", "Basal"))

ggplot(AB_PRTE, aes(PRTE_yes, Log2FoldChange_AB)) + 
  geom_boxplot(aes(color = PRTE_yes)) + 
  theme_classic(16) + 
  xlab(element_blank())  + 
  scale_color_brewer(palette = "Paired") + 
  theme(legend.position = "none") + 
  stat_compare_means( method = "t.test", ref.group = "Does Not Contain PRTE")

ggplot(AB_PRTE, aes(localized, all_counts_PRTE)) + 
  geom_violin(aes(color = localized)) + 
  theme_classic(16) + 
  xlab(element_blank()) + 
  ylab("Midbody / Whole Cell, Log2") + 
  scale_color_brewer(palette = "Paired") + 
  theme(legend.position = "none") + 
  stat_compare_means( method = "t.test", ref.group = "WC")
```
