# Deliverable specs:
#   Number of variants that gene has in the dataset, normalized by gene length
#   pLI score of the gene
#   PhyloP consesrvation
#   CADD score of the variant
#   Scoring based on other mutations present in the proband.
#   Known protein domain (based on Barbara’s annotation)
#   Connectivity degree with variants CHD
#   Connectivity degree with interactome
#   Belongs to a pathway enriched in CHD known genes.


# Function to perform fisher exact test given input list and reactome ID
get_path_fisher <- function(gList, react_id, reactome, total_prots, gene2UP){
  prot_tab = reactome[which(reactome$PathID == react_id),]
  prot_tab$GeneSymbol = gene2UP$GeneSymbol[match(prot_tab$UniProtID, gene2UP$UniProt)]
  
  n_common = length(prot_tab$GeneSymbol[!is.na(prot_tab$GeneSymbol)])
  n_genes_not_in_path = length(gList) - n_common
  n_path_not_in_genes = length(prot_tab$GeneSymbol[is.na(prot_tab$GeneSymbol)])
  n_neither = length(total_prots) - (n_common + n_genes_not_in_path + n_path_not_in_genes)
  
  fisherMat <- matrix(c(n_common, n_genes_not_in_path, n_path_not_in_genes, n_neither), nrow = 2,
                      dimnames = list(Complex = c("Yes","No"),
                                      Interactome = c("Yes", "No")))
  
  # Fisher test
  fish = fisher.test(fisherMat, alternative = 'greater')
  return(fish[1])
  
}

base_path = "/Users/student/Documents/PollardLab/APMS2"
out_path = paste0(base_path, "/output/pathogenicity_scoring")

source(paste0(base_path, "/code/load_data/load_variants.R"))

# Combine variants into a single object
cols <- c("Blinded.ID", "Cardiac.Category", "EM", "NDD", "CHR", "POS", "REF", "ALT", "Gene","Variant.Class")
#mut_table <- rbind(DNV_cases[,cols], LoF_cases[,cols])
mut_table = DNV_cases[,cols]

quant_unimap = read.table(paste0(base_path, "/input/aliases/combined_unimap.txt"),
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)
quant_unimap <- unique(quant_unimap)

######################## Step 1: Select genes in interactomes with DNVs ######################## 

# Load in protein scoring
tabG = read.table(file = paste0(base_path, "/intermediate/interaction_scoring/saintq/GATA4_protein_norm_true.txt"), 
                  sep = "\t", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, comment.char = "")
tabT = read.table(file = paste0(base_path,"/intermediate/interaction_scoring/saintq/TBX5_peptide_norm_true.txt"),
                  sep = "\t", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, comment.char = "")

# Choose significant genes
tab_sigG = subset(tabG, tabG$BFDR < 0.001)
tab_sigT = subset(tabT, tabT$BFDR < 0.05)
cols = c("Bait", "Prey", "BFDR")
tab_sig = rbind(tab_sigG[,cols], tab_sigT[,cols])

# Remove non-nuclear genes from interactome
names(tab_sig)[which(names(tab_sig) == "Prey_proteinname")] <- "Prey"
non_nuclear = read.table(paste0(base_path, "/input/databases/localization_non_nuclear_BINGO_.txt"),
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE)
non_nuclear = non_nuclear[which(!non_nuclear$subcellular.location == "Nucleus"),]
tab_sig = tab_sig[which(!tab_sig$Prey %in% non_nuclear$UniprotID),]

############# NOTE: Remove blacklist genes, if this winds up being applicable #############

# Get gene names for interactome
interactome_prots = tab_sig$Prey
interactome_genes = quant_unimap$GeneSymbol[match(interactome_prots, quant_unimap$UniProt)]
interactome_genes = unlist(strsplit(interactome_genes, "; "))

# Reduce variant list to those that fall in interactome genes
mut_table <- mut_table[which(mut_table$Gene %in% interactome_genes),]
  
######################## Step 2: Number of variants, normalized by length ######################## 
exon_coords <- read.table(paste0(base_path, "/input/databases/gene_start_stop.txt"), 
                            sep = "\t", header = TRUE)
exon_coords$CDS.Length <- NULL
exon_coords$geneLength <- exon_coords$Gene.end..bp. - exon_coords$Gene.start..bp.
exon_coords <- unique(exon_coords)
  
freq_table <- as.data.frame(table(mut_table$Gene))
freq_table$geneLength <- exon_coords$geneLength[match(freq_table$Var1, exon_coords$Gene.name)]
freq_table$mutperkb <- (freq_table$Freq / freq_table$geneLength) * 1000
  
mut_table$mutperkb <- freq_table$mutperkb[match(mut_table$Gene, freq_table$Var1)]
  
######################## Step 3: PhyloP conservation score ######################## 
  
# Load in PhyloP data
phylop <- read.table(paste0(base_path, "/input/databases/phyloP.txt"), 
                     sep = "\t", header = FALSE, stringsAsFactors = FALSE)
names(phylop) = c("CHR", "POS", "STOP", "phyloP")
phylop = phylop[which(phylop$POS %in% mut_table$POS),]
  
# Create stop column in gene data
mut_table$STOP = 0
for (i in c(1:nrow(mut_table))){
  mut_table$STOP[i] = max(c(nchar(mut_table$REF[i]),
                            nchar(mut_table$ALT[i]))) + as.numeric(mut_table$POS[i]) - 1
}
  
# Merge by Chr/Start/Stop data - why different lengths?
mg = merge(mut_table, phylop, by = c("CHR", "POS","STOP"), all.x=TRUE)
mg[which(duplicated(mg[,c('CHR','POS','STOP')])),]

gene_data <- mg
  
######################## Step 4: Cell-specific expression score ######################## 

broad_expr = read.csv(paste0(base_path, "/input/rnaseq/BroadPopulations_AverageExpression.csv"))
myo_expr = read.csv(paste0(base_path, "/input/rnaseq/Myocardium_AverageExpression.csv"))

broad_expr$Myocardium <- NULL
row.names(myo_expr) = myo_expr[,1]
myo_expr[,1] <- NULL
myo_expr$max <- NULL
myo_expr$max = apply(myo_expr, 1, function(x) max(x))

broad_expr$myocardium = myo_expr$max[match(broad_expr$X, row.names(myo_expr))]
broad_expr$myo_spec = broad_expr$myocardium / (rowSums(broad_expr[,c(2:7)]))
broad_expr$myo_spec[is.nan(broad_expr$myo_spec)] = 0
broad_expr$myo_spec[is.infinite(broad_expr$myo_spec)] = 200

# Homology between mouse and human gene names
homologene = read.table(paste0(base_path, "/input/databases/homologene.data"), 
                        sep="\t", header = TRUE, stringsAsFactors = FALSE, fill=T)
hom_human = homologene[which(homologene$taxon == "9606"),]
hom_mouse = homologene[which(homologene$taxon == "10090"),]
homologene = merge(hom_human, hom_mouse, by = "homologene")
broad_expr$Gene = homologene$symbol.x[match(broad_expr$X, homologene$symbol.y)]

### Note: this is not ideal; is there a better homology conversion method?
library(dplyr)
broad_expr$Gene[is.na(broad_expr$Gene)] = toupper(broad_expr$X[is.na(broad_expr$Gene)])

gene_data$specificity_score = broad_expr$myo_spec[match(gene_data$Gene, broad_expr$Gene)]

  
######################## Step 5: CADD scores ########################

# Save out vcf file to run through CADD
#write.table(gene_data[,c("CHR", "POS", "Blinded.ID", "REF", "ALT")], 
#            paste0(base_path, "/intermediate/CADD/vcf_CADD_saintq_n_input.vcf"),
#            row.names = FALSE, quote = FALSE, sep = "\t")

#### MANUAL: UPLOAD FILE TO https://cadd.gs.washington.edu/score ####

CADD = read.table(paste0(base_path, "/intermediate/CADD/saintq_n_interactors_CADD.tsv"),
           sep = "\t", header = FALSE, stringsAsFactors = FALSE)
names(CADD) = c("CHR", "POS", "REF", "ALT", "Raw.score", "CADD")
gd = merge(gene_data, CADD[,c('CHR','POS','REF','ALT','CADD')], by = c("CHR","POS","REF","ALT"), all.x=T)
#gd$scaled.CADD = range01(as.numeric(gd$Raw.score))
#gd$scaled.CADD[which(gd$scaled.CADD == 0)] <- min(gd$scaled.CADD[which(gd$scaled.CADD != 0)])/2 # remove 0 values
#gene_data = gd[,c("Blinded.ID", "Cardiac.Category", "NDD", "CHR", "POS", "REF", "ALT", "Gene", 
#                  "pLI.Score", "norm_mutperkb", "phyloP", "specificity_score", "scaled.CADD")]


######################## Step 6: Binary - other mutations, known protein domain, CHD pathway ########################

# Determine if that person has another mutation - 0 for yes, 1 for no
#### NOTE: This includes whether the variant itself occurs in a known gene
known = read.table(paste0(base_path, "/input/databases/known_genes.txt"), stringsAsFactors = F)
known_genes = known$V1
cols <- c("Blinded.ID", "Cardiac.Category", "EM", "NDD", "CHR", "POS", "REF", "ALT", "Gene","pLI.Score")
all_muts <- rbind(DNV_cases[,cols], LoF_cases[,cols])

gd$no_known_muts = 1
for (i in c(1:nrow(gd))){
  muts = all_muts[which(all_muts$Blinded.ID == gd$Blinded.ID[i]),]
  if (any(muts$Gene %in% known_genes)){
    gd$no_known_muts[i] = 0
  }
}

# Determine if this occurs in a known protein domain
###### OR is a lof variant - we expect these proteins to not do their job #######
domains = read.csv(paste0(base_path, "/intermediate/domain_annotations.csv"), stringsAsFactors = F)
gene_data = merge(gd, domains[,c("CHR","POS","REF","ALT","Known.interaction","NEARBY.VARIANTS","BAIT",
                                 "Variant.protein.domain","Protein.domains.affected")], 
                  by = c("CHR","POS","REF","ALT"), all.x = TRUE)
gene_data$protein_domain_or_lof = 1
neg_str = c("no specific domain","non specific domain","not known domain","")
gene_data$protein_domain_or_lof[which(gene_data$Protein.domains.affected %in% neg_str)] <- 0.5
gene_data$protein_domain_or_lof[which(gene_data$Variant.Class %in% c("non","frameshift","startloss","splice",
                                                                     "stoploss"))] <- 1

# Determine if this occurs in a pathway enriched for known CHD genes

# Edit reactome file if necessary
# NOTE: Line 113879 was somehow corrupted; it contains no human pathways, so I did not correct it.
react = read.table(paste0(base_path, "/input/databases/reactome.txt"), 
                   sep = "\t", stringsAsFactors = F, comment.char = "")
react <- react[which(react$V6 == "Homo sapiens"), c('V1','V2','V4','V6')]
names(react) <- c("UniProtID","PathID","PathName","Species")

# Gene to UP accession file
conv = read.table(paste0(base_path, "/input/aliases/reactome_unimap.txt"), 
                  sep = "\t", header = T, stringsAsFactors = F)
s <- strsplit(conv$UniProt.Accession, split = "; ")
gene2UP = data.frame(GeneSymbol = rep(conv$Gene.Symbol, sapply(s, length)), UniProt = unlist(s))

# For each reactome pathway, calculate significance of enrichment for known genes
react$enrichment_p = ">0.05"
total_prots = unique(react$UniProtID)
sig_paths = c()
n_paths = length(unique(react$PathID))
for (react_ID in unique(react$PathID)){
  fish = get_path_fisher(known_genes, react_ID, react, total_prots, gene2UP)
  if (fish <= 0.05){
    react$enrichment_p[which(react$PathID == react_ID)] <- unlist(fish)
    sig_paths = c(sig_paths, react_ID)
  }
}

# Examine significantly-enriched pathways; if necessary, save out uniprot IDs in those pathways
check <- react[which(react$enrichment_p != ">0.05"),c("UniProtID","PathID","PathName","enrichment_p")]
#write.table(unique(check$UniProtID), file = paste0(base_path, "/input/aliases/sig_reactome_prots.txt"),
#            row.names = F, col.names = F, quote = F)
check$Gene = quant_unimap$GeneSymbol[match(check$UniProtID, quant_unimap$UniProt)]

# Write out info to examine
write.csv(check, file = paste0(base_path, "/intermediate/enriched_pathway_info.csv"), row.names = F, quote = F)

# If the gene is in an enriched pathway, pathway = 1; else, pathway = 0.5
### NOTE: None of the DNVs are in enriched pathways....alternative to consider all pathways, not just enriched??
gene_data$enriched_path = 0.5
gene_data$enriched_path[which(gene_data$Gene %in% check$Gene)] <- 1


######################## Step 7: Connectivity degrees - interactome, known CHD genes ########################
PPI = read.table(paste0(base_path, "/input/databases/iRefIndex.txt"), 
                 sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Connectivity to interactome genes
first_interactome = PPI[which(PPI$aliasA %in% interactome_genes | PPI$aliasB %in% interactome_genes),]
gene_data$interactome_node_degree = 0
for(i in c(1:nrow(gene_data))){
  lines = unique(first_interactome[which(first_interactome$aliasA == gene_data$Gene[i] | 
                                    first_interactome$aliasB == gene_data$Gene[i]), c('aliasA','aliasB')])
  
  # Remove self-loops
  lines = lines[which(lines$aliasA != lines$aliasB),]
  
  # Remove transitive interactions
  swap_lines = lines[,c('aliasB','aliasA')]
  names(swap_lines) = c('aliasA','aliasB')
  both = rbind(lines, swap_lines)
  both = unique(both)
  degree = nrow(both)/2
  gene_data$interactome_node_degree[i] = degree
}

# Connectivity to known CHD genes
first_known_chd = PPI[which(PPI$aliasA %in% known_genes | PPI$aliasB %in% known_genes),]
gene_data$known_chd_node_degree = 0
for(i in c(1:nrow(gene_data))){
  lines = unique(first_known_chd[which(first_known_chd$aliasA == gene_data$Gene[i] | 
                                           first_known_chd$aliasB == gene_data$Gene[i]), c('aliasA','aliasB')])
  
  # Remove self-loops
  lines = lines[which(lines$aliasA != lines$aliasB),]
  
  # Remove transitive interactions
  swap_lines = lines[,c('aliasB','aliasA')]
  names(swap_lines) = c('aliasA','aliasB')
  both = rbind(lines, swap_lines)
  both = unique(both)
  degree = nrow(both)/2
  gene_data$known_chd_node_degree[i] = degree
}

######################## Step 8: Impute missing, rank-normalize, and sum ########################

# Update pLI score with most recent gnomAD release
gnomAD = read.table(paste0(base_path, "/input/databases/gnomad_pli.txt"), sep = "\t", header=T, stringsAsFactors = F)
gnomAD <- gnomAD[,c("gene","pLI","oe_lof","oe_mis")]
names(gnomAD) <- c("Gene","pLI_gnomAD","oe_lof","oe_mis")
gene_data <- merge(gene_data, gnomAD, by = "Gene", all.x = TRUE)

# Save relevant oe
gene_data$oe = NA
for(i in c(1:nrow(gene_data))){
  if (gene_data$Variant.Class[i] %in% c("mis","misD")){
    gene_data$oe[i] = gene_data$oe_mis[i]
  } else if(gene_data$Variant.Class[i] %in% c("non","frameshift","startloss","splice","stoploss")){
    gene_data$oe[i] = gene_data$oe_lof[i]
  }
}

# Median-permute missing data
medify <- function(dataframe, col_name){
  med = median(as.numeric(dataframe[,col_name][!is.na(dataframe[,col_name])]))
  dataframe[,col_name][which(is.na(dataframe[,col_name]))] <- med
  return(dataframe)
}

gene_data = medify(gene_data, "phyloP")
gene_data = medify(gene_data, "specificity_score")

# Create rank_cols for applicable columns: mutperkb, phyloP, specificity_score, CADD, pLI_gnomAD, oe
rank_col <- function(dataframe, col_name_list){
  for(col_name in col_name_list){
    rank_col_name = paste0(col_name, "_rank")
    if (col_name != "oe"){
      dataframe[,rank_col_name] <- rank(dataframe[,col_name])
    } else{
      dataframe[,rank_col_name] <- rank(-dataframe[,col_name])
    }
  }
  return(dataframe)
}

rank_data <- rank_col(gene_data, c("mutperkb","specificity_score","CADD","phyloP","interactome_node_degree",
                                   "known_chd_node_degree","pLI_gnomAD","oe"))
rank_data$rank_sum =rowSums(rank_data[,c("mutperkb_rank","specificity_score_rank","CADD_rank","phyloP_rank","oe_rank",
                                          "interactome_node_degree_rank","known_chd_node_degree_rank",
                                          "pLI_gnomAD_rank")])

# Save out rank order before removing known genes
gd <- rank_data[,c('CHR','POS','REF','ALT','Gene','Blinded.ID','Cardiac.Category','NDD',"CADD","phyloP","pLI_gnomAD",
                   "oe","mutperkb","specificity_score", "interactome_node_degree","known_chd_node_degree",
                   "no_known_muts","protein_domain_or_lof", "enriched_path","rank_sum")]
f1 <- gd[order(-gd$rank_sum),]

rank_data$binary_multiplied = rank_data$rank_sum * 
  rank_data$no_known_muts * rank_data$protein_domain_or_lof * rank_data$enriched_path

# Save out desired columns
gd <- rank_data[,c('CHR','POS','REF','ALT','Gene','Blinded.ID','Cardiac.Category','NDD',"CADD","phyloP","pLI_gnomAD",
                   "oe","mutperkb","specificity_score", "interactome_node_degree","known_chd_node_degree",
                   "no_known_muts","protein_domain_or_lof", "enriched_path","rank_sum","binary_multiplied")]
f2 <- gd[order(-gd$binary_multiplied),]


# Save out
write.csv(f1, paste0(base_path, "/output/pathogenicity_scoring/saintq_n/GATA4-TBX5_interactors_simple_rank_sum.csv"),
          quote=F, row.names = F)
write.csv(f2, paste0(base_path, "/output/pathogenicity_scoring/saintq_n/GATA4-TBX5_interactors_binary_multiplied.csv"),
          quote=F, row.names = F)


# # # # # # OLD # # # # # # 

  ######################## Step 6: Combine scores and sorts ########################
  
  # Replace missing pLI with median)
  median_pLI = median(as.numeric(gene_data$pLI.Score[!is.na(gene_data$pLI.Score)]))
  gene_data$pLI.Score[is.na(gene_data$pLI.Score)] <- median_pLI
  
  # Add pseudocount to pLI
  gene_data$pLI.Score[which(as.numeric(gene_data$pLI.Score) == 0)] <- min(
    as.numeric(gene_data$pLI.Score[which(as.numeric(gene_data$pLI.Score) != 0)]))/2
  
  # Reduce dataframe
  numerical_data = as.data.frame(
    gene_data[,c("pLI.Score","norm_mutperkb","specificity_score","phyloP","scaled.CADD")])
  df2 <- mutate_all(numerical_data, function(x) as.numeric(as.character(x)))
  df2$mult_score = apply(df2, 1, prod)
  df2$add_score = rowSums(df2[,c(1:5)])
  
  final <- cbind(gene_data[,c("Gene","Blinded.ID","Cardiac.Category","NDD", "CHR","POS", "REF", "ALT")], df2)
  final_mult <- final[order(-final$mult_score),]
  final_add <- final[order(-final$add_score),]
  
  ofile = paste0(base_path, "/output/pathogenicity_scoring/", analysis_type, 
                 "/all_interactomes_median_permuted_scores_G001_T05_N1.csv")
  
  
  write.csv(final_mult, ofile,
            quote = FALSE, row.names = FALSE)
  
  # Create a file that instead uses ranks
  final$pLI.rank = rank(final$pLI.Score)
  final$mutperkb.rank = rank(final$norm_mutperkb)
  final$specificity.rank = rank(final$specificity_score)
  final$phyloP.rank = rank(final$phyloP)
  final$CADD.rank = rank(final$scaled.CADD)
  
  # lower rank = lower number
  # highest numbers are most interesting
  final$AddedRank = rowSums(final[,c("pLI.rank", "mutperkb.rank", "specificity.rank", "phyloP.rank", "CADD.rank")])
  final_rank <- final[order(-final$AddedRank),c("Gene","Blinded.ID","Cardiac.Category","NDD", "CHR","POS", 
                                                "REF", "ALT",
                                                "pLI.rank", "mutperkb.rank", "specificity.rank", "phyloP.rank", 
                                                "CADD.rank",
                                                "AddedRank")]
  ofile = paste0(base_path, "/output/pathogenicity_scoring/", analysis_type, 
                 "/all_interactomes_ranked_G001_T05_N1.csv")
  
  write.csv(final, ofile, quote = FALSE, row.names = FALSE)
  
  ######################## GATA4/TBX5 only ########################
  
  spec_genes = unique(ints$Prey_genename[which(ints$Bait == "GATA4" | ints$Bait == "TBX5")])
  GATAs = final_mult[which(final_mult$Gene %in% spec_genes),]
  
  ofile = paste0(base_path, "/output/pathogenicity_scoring/",analysis_type, 
                 "/GATA4-TBX5_interactomes_median_permuted_scores_G001_T05.csv")
  write.csv(GATAs, ofile, quote = FALSE, row.names = FALSE)
  
  # For ranked scores
  GATAs = final[which(final_mult$Gene %in% spec_genes),]
  ofile = paste0(base_path, "/output/pathogenicity_scoring/",analysis_type, 
                 "/GATA4-TBX5_interactomes_ranked_G001_T05.csv")
  write.csv(GATAs, ofile,quote = FALSE, row.names = FALSE)
  
  ######################## Nearby variants? ########################
  res = c()
  for(chr in c(1:22, "X")){
    var_list = final_mult[which(final_mult$CHR == chr),]
    if(nrow(var_list) > 1){
      
      # Sort numerically
      coords = sort(as.numeric(var_list$POS))
      genes = var_list[match(coords, var_list$POS),]
      
      # Check whether any mutations next to each other are within 800bp in either direction
      #### This is janky and assumes that each variant has a position not shared on any other chromosome
      for (i in c(1:(length(coords)-1))){
        if(coords[i] - coords[i+1] <= 800 && coords[i] - coords[i+1] >= -800){
          success= paste0("Nearby variants! chr", chr, ":", coords[i], ",", coords[i+1],
                          " in gene ", genes[i,"Gene"])
          print(success)
          res = c(res, coords[i], coords[i+1])
        }
      }
    }
  }
  
  domain = final_mult[which(final_mult$POS %in% res),]
  domain = domain[order(domain$Gene, domain$mult_score),]
  ofile = paste0(base_path, "/output/pathogenicity_scoring/",
                 analysis_type,"/nearby_variants_G001_T05_N1.csv")
  write.csv(domain, file = ofile,
            row.names = FALSE, quote = FALSE)
  
}



### Step 1.5: save out a file that includes interactome data for all variants
cols <- c("Gene","pLI.Score", "type")
DNV_cases$type = "DNV"
LoF_cases$type = "LoF"
rec_cases$type = "recessive"


for(analysis_type in c("saintq_n", "original")){
  ints = read.table(paste0(base_path, "/intermediate/interactome_lists/",
                           analysis_type,"/combined_APMS_interactome_G001_T05_N1.csv"),
                         sep = ",", header = TRUE, stringsAsFactors = FALSE)
  ints = ints[,c("Bait","Prey","Prey_genename","BFDR")]
  
  mut_data <- rbind(DNV_cases[,cols], LoF_cases[,cols], rec_cases[,cols])
  mut_data = mut_data[which(mut_data$Gene %in% ints$Prey_genename),]
  counts = as.matrix(table(mut_data$Gene, mut_data$type))
  cdf = as.data.frame(cbind(row.names(counts),counts[,c(1:3)]))
  cdf$pLI.score = mut_data$pLI.Score[match(cdf$V1, mut_data$Gene)]
  
  # Connect to baits in interactome scoring
  counts = cdf
  counts$GATA4 = 0
  counts$TBX5 = 0
  counts$NKX25 = 0
  
  GATA4_genes = unlist(strsplit(ints$Prey_genename[which(ints$Bait == "GATA4")], "; "))
  TBX5_genes = unlist(strsplit(ints$Prey_genename[which(ints$Bait == "TBX5")], "; "))
  NKX25_genes = unlist(strsplit(ints$Prey_genename[which(ints$Bait == "NKX25")], "; "))
  
  for (i in c(1:nrow(counts))){
    if(counts$V1[i] %in% GATA4_genes){
      counts$GATA4[i] = 1
    }
    if(counts$V1[i] %in% TBX5_genes){
      counts$TBX5[i] = 1
    }
    if(counts$V1[i] %in% NKX25_genes){
      counts$NKX25[i] = 1
    }
  }
  
  names(counts)[1] = "Gene"
  write.csv(counts, file = paste0(base_path, "/output/pathogenicity_scoring/", analysis_type, 
                                  "/variant_count_G001_T05_n1.csv"))
}

# Filter to only include DNVs
for (analysis_type in c("original", "saintq_n")){
  outpath = paste0(base_path, "/output/pathogenicity_scoring/", analysis_type)
  f1 = "/all_interactomes_median_permuted_scores_G001_T05_N1.csv"
  f2 = "/GATA4-TBX5_interactomes_median_permuted_scores_G001_T05.csv"
  f3 = "/nearby_variants_G001_T05_N1.csv"
  f4 = "/variant_count_G001_T05_n1.csv"
  f5 = "/all_interactomes_ranked_G001_T05_N1.csv"
  f6 = "/GATA4-TBX5_interactomes_ranked_G001_T05.csv"
  
  for(infile in c(f1, f2, f3, f4, f5, f6)){
    tab = read.csv(paste0(outpath, infile))
    DNV_tab = tab[which(tab$POS %in% DNV_cases$POS),]
    write.csv(DNV_tab, file = paste0(base_path, "/output/pathogenicity_scoring/", analysis_type,
                                     "/DNV_only/", infile))
  }
  
}


