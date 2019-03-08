library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Hs.eg.db)

# Load in RNA seq data - counts per million
cmd = paste0("cut -d ' ' -f 1,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37 ", base_path,
             "/input/rnaseq/edgeR_all_description.txt")
cpm <- read.table(pipe(cmd), header = TRUE, stringsAsFactors = FALSE)
cpm <- unique(cpm)

# Load in RNAseq data - raw counts
cmd = paste0("cut -d ' ' -f 1,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21 ", base_path,
             "/input/rnaseq/edgeR_all_description.txt")
seqdata = read.table(pipe(cmd), header = TRUE, stringsAsFactors = FALSE)
seqdata <- unique(seqdata)

## Store gene names as rownames
row.names(seqdata) = seqdata$Gene_symbol
row.names(cpm) = cpm$Gene_symbol