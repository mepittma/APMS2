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
#### NOTE: This DOES NOT include whether the variant itself occurs in a known gene
known = read.table(paste0(base_path, "/input/databases/known_genes.txt"), stringsAsFactors = F)
known_genes = known$V1
cols <- c("Blinded.ID", "Cardiac.Category", "EM", "NDD", "CHR", "POS", "REF", "ALT", "Gene","pLI.Score")
all_muts <- rbind(DNV_cases[,cols], LoF_cases[,cols])

gd$no_known_muts = 1
for (i in c(1:nrow(gd))){
  muts = all_muts[which(all_muts$Blinded.ID == gd$Blinded.ID[i]),]
  muts = muts[which(!(muts$CHR == gd$CHR[i] & muts$POS == gd$POS[i])),] # don't include the variant itself
  if (any(muts$Gene %in% known_genes)){
    gd$no_known_muts[i] = 0
  }
}

# Column to state whether the gene the variant is in is known or not
gd$known_gene = "unknown"
gd$known_gene[which(gd$Gene %in% known_genes)] <- "known"

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
enriched_paths = read.csv(paste0(base_path, "/intermediate/InnateDB_sig_paths.csv"), stringsAsFactors = F)
enriched_paths = enriched_paths[which(enriched_paths$SourceDatabase %in% c("REACTOME", "PID BIOCARTA")),]
path2gene = read.csv(paste0(base_path, "/intermediate/path2gene_table.csv"), stringsAsFactors = F)
path_genes = path2gene[which(path2gene$Pathway.Name %in% enriched_paths$PathwayName),]
gene_data$enriched_path = 0.5
gene_data$enriched_path[which(gene_data$Gene %in% path_genes$Query.Xref)] <- 1


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
gd <- rank_data[,c('CHR','POS','REF','ALT','Gene','Variant.Class','Blinded.ID','Cardiac.Category','NDD',"CADD","phyloP","pLI_gnomAD",
                   "oe","mutperkb","specificity_score", "interactome_node_degree","known_chd_node_degree",
                   "no_known_muts","known_gene","protein_domain_or_lof", "enriched_path","rank_sum")]
f1 <- gd[order(-gd$rank_sum),]

rank_data$binary_multiplied = rank_data$rank_sum * 
  rank_data$no_known_muts * rank_data$protein_domain_or_lof * rank_data$enriched_path

# Save out desired columns
gd <- rank_data[,c('CHR','POS','REF','ALT','Gene','Variant.Class','Blinded.ID','Cardiac.Category','NDD',"CADD","phyloP","pLI_gnomAD",
                   "oe","mutperkb","specificity_score", "interactome_node_degree","known_chd_node_degree",
                   "no_known_muts","known_gene","protein_domain_or_lof", "enriched_path","rank_sum","binary_multiplied")]
f2 <- gd[order(-gd$binary_multiplied),]


# Save out
write.csv(f1, paste0(base_path, "/output/pathogenicity_scoring/saintq_n/GATA4-TBX5_interactors_simple_rank_sum.csv"),
          quote=F, row.names = F)
write.csv(f2, paste0(base_path, "/output/pathogenicity_scoring/saintq_n/GATA4-TBX5_interactors_binary_multiplied.csv"),
          quote=F, row.names = F)

