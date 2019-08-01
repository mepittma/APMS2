base_path = "/Users/student/Documents/PollardLab/APMS2"
out_path = paste0(base_path, "/output/pathogenicity_scoring")

source(paste0(base_path, "/code/load_data/load_variants.R"))

####################### Kihyun ######################## 

mut_table = DNV_cases[,cols]
gene_data = mut_table
analysis_type = "saintq_n"
ints = read.csv(paste0(base_path, "/intermediate/interactome_lists/",
                       analysis_type,"/combined_APMS_interactome_G001_T05_N1.csv"),
                stringsAsFactors = FALSE, header = TRUE)
int_genes = unique(ints$Prey_genename)
int_genes = unique(unlist(strsplit(int_genes, "; ")))

# mutperkb
exon_coords <- read.table(paste0(base_path, "/input/databases/gene_start_stop.txt"), 
                          sep = "\t", header = TRUE)
exon_coords$CDS.Length <- NULL
exon_coords$geneLength <- exon_coords$Gene.end..bp. - exon_coords$Gene.start..bp.
exon_coords <- unique(exon_coords)

freq_table <- as.data.frame(table(gene_data$Gene))
freq_table$geneLength <- exon_coords$geneLength[match(freq_table$Var1, exon_coords$Gene.name)]
freq_table$mutperkb <- (freq_table$Freq / freq_table$geneLength) * 1000

gene_data$gene_length <- freq_table$geneLength[match(gene_data$Gene, freq_table$Var1)]
gene_data$n_mutations <- freq_table$Freq[match(gene_data$Gene, freq_table$Var1)]
gene_data$mutperkb <- freq_table$mutperkb[match(gene_data$Gene, freq_table$Var1)]
# Normalize the mutperkb score to be between 0.5 and 1 - we don't want to penalize
# rarely-mutated genes too harshly
range51 <- function(x){(0.5*(x-min(x))/(max(x)-min(x)) + 0.5)}
gene_data$mutperkb[is.na(gene_data$mutperkb)] = 0.000001
gene_data$mutperkb = gene_data$mutperkb + 0.000001
gene_data$norm_mutperkb = range51(log(gene_data$mutperkb))

# Is the gene in Barbara's interactome
gene_data$Interactome = "No"
gene_data$Interactome[which(gene_data$Gene %in% int_genes)] <- "Yes"
gd <- gene_data[,c("Gene", "gene_length", "n_mutations","mutperkb", "pLI.Score","Interactome")]
gd <- unique(gd)

gd$pLI.rank = rank(gd$pLI.Score)
gd$mutperkb.rank = rank(gd$mutperkb)
gd$AddedRank = rowSums(gd[,c("pLI.rank", "mutperkb.rank")])
final_rank <- gd[order(-gd$AddedRank),]
ofile = paste0(base_path, "/output/pathogenicity_scoring/saintq_n/Kihyun_genes_ranked_NApLIhigh.csv")

write.csv(final_rank, ofile, row.names = FALSE)


# Coerce NAs to 0s
gd$pLI.Score[which(gd$pLI.Score %in% c("NA","N/A"))] <- 0
gd$pLI.Score[is.na(gd$pLI.Score)] <- 0
gd$pLI.rank = rank(gd$pLI.Score)
gd$AddedRank = rowSums(gd[,c("pLI.rank", "mutperkb.rank")])
final_rank <- gd[order(-gd$AddedRank),]
ofile = paste0(base_path, "/output/pathogenicity_scoring/saintq_n/Kihyun_genes_ranked_NApLI_low.csv")

write.csv(final_rank, ofile, row.names = FALSE)



# Combine variants into a single object
cols <- c("Blinded.ID", "Cardiac.Category", "EM", "NDD", "CHR", "POS", 
          "REF", "ALT", "Gene","pLI.Score")

mut_table = DNV_cases[,cols]
mut_table$Cardiac.Category[which(mut_table$Cardiac.Category == "CTD ")] <- "CTD"
mut_table$Cardiac.Category[which(mut_table$Cardiac.Category == " CTD (TGA)")] <- "CTD (TGA)"
mut_table$Cardiac.Category[which(mut_table$Cardiac.Category == "other")] <- "OTHER"
mut_table$Cardiac.Category[which(mut_table$Cardiac.Category == "Other")] <- "OTHER"
mut_table$Cardiac.Category[which(mut_table$Cardiac.Category == "Other (AVC)")] <- "other (AVC)"

diag <- read.csv(paste0(base_path, "/input/variants/diagnoses.csv"), stringsAsFactors = FALSE)
mut_table$Diagnosis <- diag$Cardiac.Diagnoses[match(mut_table$Blinded.ID, diag$Blinded.ID)]

analysis_type = "saintq_n"
ints = read.csv(paste0(base_path, "/intermediate/interactome_lists/",
                       analysis_type,"/combined_APMS_interactome_G001_T05_N1.csv"),
                stringsAsFactors = FALSE, header = TRUE)
int_genes = unique(ints$Prey_genename)
int_genes = unique(unlist(strsplit(int_genes, "; ")))

# mutperkb
gene_data <- mut_table
exon_coords <- read.table(paste0(base_path, "/input/databases/gene_start_stop.txt"), 
                          sep = "\t", header = TRUE)
exon_coords$CDS.Length <- NULL
exon_coords$geneLength <- exon_coords$Gene.end..bp. - exon_coords$Gene.start..bp.
exon_coords <- unique(exon_coords)

freq_table <- as.data.frame(table(gene_data$Gene))
freq_table$geneLength <- exon_coords$geneLength[match(freq_table$Var1, exon_coords$Gene.name)]
freq_table$mutperkb <- (freq_table$Freq / freq_table$geneLength) * 1000

gene_data$gene_length <- freq_table$geneLength[match(gene_data$Gene, freq_table$Var1)]
gene_data$n_mutations <- freq_table$Freq[match(gene_data$Gene, freq_table$Var1)]
gene_data$mutperkb <- freq_table$mutperkb[match(gene_data$Gene, freq_table$Var1)]
# Normalize the mutperkb score to be between 0.5 and 1 - we don't want to penalize
# rarely-mutated genes too harshly
range51 <- function(x){(0.5*(x-min(x))/(max(x)-min(x)) + 0.5)}
gene_data$mutperkb[is.na(gene_data$mutperkb)] = 0.000001
gene_data$mutperkb = gene_data$mutperkb + 0.000001
gene_data$norm_mutperkb = range51(log(gene_data$mutperkb))



genelist = c()
n_total_patients_with_mutation = c()
n_NDD_positive = c()
n_NDD_unknown = c()
n_CTD = c()
n_HTX = c()
n_LVO = c()
n_Other_category = c()
n_ToF = c()
n_Truncus_Arteriosus = c()

for (gene in c(unique(gene_data$Gene))){
  
  genelist = c(genelist, gene)
  
  subs = mut_table[mut_table$Gene == gene,]
  n_indi = length(unique(subs$Blinded.ID))
  n_total_patients_with_mutation = c(n_total_patients_with_mutation, n_indi)
  
  sub_NDD_yes = subs[which(subs$NDD == "Yes"),]
  n_NDD = length(unique(sub_NDD_yes$Blinded.ID))
  n_NDD_positive = c(n_NDD_positive, n_NDD)
  
  sub_NDD_unk = subs[which(subs$NDD == "Unknown"),]
  n_unk = length(unique(sub_NDD_unk$Blinded.ID))
  n_NDD_unknown = c(n_NDD_unknown, n_unk)
  
  sub_CTD = subs[which(subs$Cardiac.Category %in% c("CTD", "CTD (TGA)")),]
  sub_HTX = subs[which(subs$Cardiac.Category =="HTX"),]
  sub_LVO = subs[which(subs$Cardiac.Category =="LVO"),]
  sub_other = subs[which(subs$Cardiac.Category %in% c("OTHER", "other (AVC)")),]
  
  n_CTD = c(n_CTD, nrow(sub_CTD))
  n_HTX = c(n_HTX, nrow(sub_HTX))
  n_LVO = c(n_LVO, nrow(sub_LVO))
  n_Other_category = c(n_Other_category, nrow(sub_other))
  
  sub_ToF = subs[grepl("ToF", subs$Diagnosis) | grepl("TETRALOGY OF FALLOT", subs$Diagnosis),]
  n_tf = length(unique(sub_ToF$Blinded.ID))
  n_ToF = c(n_ToF, n_tf)
  
  sub_TA = subs[grepl("TRUNCUS ARTERIOSUS", subs$Diagnosis),]
  n_ta = length(unique(sub_TA$Blinded.ID))
  n_Truncus_Arteriosus = c(n_Truncus_Arteriosus, n_ta)
  
}

# Merge gene_df together
diag_df = as.data.frame(cbind(genelist, n_total_patients_with_mutation, n_NDD_positive, n_NDD_unknown,
                n_CTD, n_HTX, n_LVO, n_Other_category, n_ToF, n_Truncus_Arteriosus))

# Is the gene in Barbara's interactome
gene_data$Interactome = "No"
gene_data$Interactome[which(gene_data$Gene %in% int_genes)] <- "Yes"
gd <- gene_data[,c("Gene", "gene_length", "n_mutations","mutperkb", "pLI.Score","Interactome")]
gd <- unique(gd)

# Merge the two
gene_data1 <- merge(gd, diag_df, by.x = "Gene", by.y = "genelist", all.x = T, all.y = F)
gd = gene_data1

gd$pLI.rank = rank(gd$pLI.Score)
gd$mutperkb.rank = rank(gd$mutperkb)
gd$ToF.rank = rank(gd$n_ToF)
gd$TA.rank = rank(gd$n_Truncus_Arteriosus)
gd$AddedRank = rowSums(gd[,c("pLI.rank", "mutperkb.rank", "ToF.rank", "TA.rank")])
final_rank <- gd[order(-gd$AddedRank),]
ofile = paste0(base_path, "/output/pathogenicity_scoring/saintq_n/genes_ranked_NApLIhigh_diagnoses.csv")

write.csv(final_rank, ofile, row.names = FALSE)


# Coerce NAs to 0s
gd$pLI.Score[which(gd$pLI.Score %in% c("NA","N/A"))] <- 0
gd$pLI.Score[is.na(gd$pLI.Score)] <- 0
gd$pLI.rank = rank(gd$pLI.Score)
gd$ToF.rank = rank(gd$n_ToF)
gd$TA.rank = rank(gd$n_Truncus_Arteriosus)
gd$AddedRank = rowSums(gd[,c("pLI.rank", "mutperkb.rank", "ToF.rank", "TA.rank")])
final_rank <- gd[order(-gd$AddedRank),]
ofile = paste0(base_path, "/output/pathogenicity_scoring/saintq_n/genes_ranked_NApLI_low_diagnoses.csv")

write.csv(final_rank, ofile, row.names = FALSE)
