library(tidyverse)
library(dplyr)
library(DESeq2)
library(EdgeR)

kallisto_files <- dir("/datastore/Chidimma/Trainome_data/Training_data/Kallisto_output/")
kallisto_files

#create function to read kallisto output files
read_kallisto_output <- function(file){
  df <- read.table(file, header = TRUE)
  df <- separate(df, target_id, c("transcript_ID", "gene_ID", "Havana_gene_ID",
                                  "Havana_transcript_ID", "transcript_name",
                                  "gene_name", "sequence_length", "transcript_biotype"), sep = "\\|")
  df <- df%>%select(transcript_ID, gene_ID, transcript_name, gene_name,
                    transcript_biotype, length, est_counts, tpm)
  #extract all those with est_counts above 1
  df <- filter(df, est_counts > 1)
  return(df)
}

one <- read_kallisto_output("/datastore/Chidimma/Trainome_data/Training_data/Kallisto_output/1-FP12w0R.tsv")

#extracting only long_non-coding RNAS based on EBI's definition

lncRNAs_one <- one %>% filter(transcript_biotype %in% c( "processed_transcript", "lncRNA",
                                                        "lincRNA", "3prime_overlapping_ncrna", 
                                                        "antisense", "non_coding", 
                                                        "sense_intronic", "sense_overlapping",
                                                        "TEC", "known_ncrna", "bidirectional_promoter_lncrna",
                                                       "macro_lncRNA" ), est_counts > 1)
#function to extract all the lncRNAs according to EBI's definition
read_ebi_lncRNAs <- function(file){
  df <- read_kallisto_output(file)
  df <- df %>% filter(transcript_biotype %in% c( "processed_transcript", "lncRNA",
                                                "lincRNA", "3prime_overlapping_ncrna", 
                                                "antisense", "non_coding", 
                                                "sense_intronic", "sense_overlapping",
                                                "TEC", "known_ncrna", "bidirectional_promoter_lncrna",
                                                "macro_lncRNA" ), est_counts >1)
  return(df)
}

lncRNAS_STRICT <- ebi_lncRNAs("/datastore/Chidimma/Trainome_data/Training_data/Kallisto_output/1-FP12w0R.tsv")

#function to extract strictly those termed lncRNAs

read_lncRNAs <- function(file){
  df <-read_kallisto_output(file)
  df <- filter(df, transcript_biotype == "lncRNA")
  return(df)
}
lncRNAS_STRICT <- lncRNAs("/datastore/Chidimma/Trainome_data/Training_data/Kallisto_output/1-FP12w0R.tsv")    


#function to extract protein_coding genes
read_protein_coding_genes <- function(file){
  df <- read_kallisto_output(file)
  df <- filter(df, transcript_biotype == "protein_coding")
  return(df)
}

one_protein_coding <- read_protein_coding_genes("/datastore/Chidimma/Trainome_data/Training_data/Kallisto_output/1-FP12w0R.tsv")

#read those described as "processed_transcripts"
read_processed_transcript <- function(file){
  df <- read_kallisto_output(file)
  df <- filter(df, transcript_biotype == "processed_transcript")
  return(df)
}

one_processed <- read_processed_transcript("/datastore/Chidimma/Trainome_data/Training_data/Kallisto_output/1-FP12w0R.tsv")
#Loop and read the entire file.
imports <- list()
for (i in 1:length(kallisto_files)){
  temp <- read_kallisto_output(kallisto_files[i])
  temp$file_name <- kallisto_files[[i]]
  imports[[i]] <- temp
}

df <- bind_rows(imports)

#Loop and read the entire lncRNA EBI def

import2 <- list()
for (i in 1:length(kallisto_files)){
  lnc <- read_ebi_lncRNAs(kallisto_files[i])
  lnc$file_name <- kallisto_files[[i]]
  import2[[i]] <- lnc
  
}
lncRNAs_all <- bind_rows(import2)


#loop through and read those strictly lncRNA
import3 <- list()
for (i in 1:length(kallisto_files)){
  specific_lncs <- read_lncRNAs(kallisto_files[i])
  specific_lncs$file_name <- kallisto_files[[i]]
  import3[[i]] <- specific_lncs
}
lncRNA_ALL_STRICT <- bind_rows(import3)


import4 <- list()
for (i in 1:length(kallisto_files)){
  processed_transcripts <- read_processed_transcript(kallisto_files[i])
  processed_transcripts$file_name <- kallisto_files[[i]]
  import4[[i]] <- processed_transcripts
}

processed_transcripts_all <- bind_rows(import4)
