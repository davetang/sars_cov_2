#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(BiocParallel))

my_fasta <- readDNAStringSet("../raw/genbank_efetch.fa")
my_seq_len <- sapply(X = my_fasta, FUN = length)
my_seq <- my_fasta[my_seq_len > 29000]

my_param <- MulticoreParam(workers = 16)

my_comb <- combn(x = 1:length(my_seq[1:3]), m = 2)

my_comb_list <- lapply(apply(X = my_comb, MARGIN = 2, FUN = list), unlist)

align_wrapper <- function(x){
   pairwiseAlignment(pattern = my_seq[x[1]], subject = my_seq[x[2]], type = "overlap")
}

system.time(
   my_alignments <- bplapply(my_comb_list, align_wrapper, BPPARAM = my_param)
)

saveRDS(my_alignments, "../result/alignments.rds")

