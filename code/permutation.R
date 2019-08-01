######################## Setup ######################## 
base_path = "/Users/student/Documents/PollardLab/APMS2"
#out_path = paste0(base_path, "/output/permutations/saintq_n")
#out_path = paste0(base_path, "/output/permutations/original")
#out_path = paste0(base_path, "/output/permutations/saintq_n_cutoff_G02_T02_N1")
#out_path = paste0(base_path, "/output/permutations/saintq_n_cutoff_G001_T05_N1")
#out_path = paste0(base_path, "/output/permutations/saintq_n_enrichedComp")
#out_path = paste0(base_path, "/output/permutations/original_enrichedComp")
#out_path = paste0(base_path, "/output/permutations/saintq_n_unexpanded")
#out_path = paste0(base_path, "/output/permutations/saintq_n_unexpanded_G001_T05_N1")
#out_path = paste0(base_path, "/output/permutations/conference_figs")
#out_path = paste0(base_path, "/output/permutations/g296s")
out_path = paste0(base_path, "/output/july22_test")

library(ggplot2)
source(paste0(base_path, "/code/functions/permutation_functions.R"))
#get_fisher, get_OR, get_prots, get_sig_compl, perm_viz, permute_status
source(paste0(base_path, "/code/functions/score_interactomes.R"))
#get_comps, get_prots, get_fisher, iref_response_list

# Read in corum data
corum <- read.table(paste0(base_path, "/input/databases/corum_human_CytoscapeFormatted.txt"), 
                    sep = "\t", quote = "", fill = TRUE, header = TRUE, stringsAsFactors = FALSE)
corum <- as.data.frame(corum[,c(1:8)])

# Read in variant data
source(paste0(base_path, "/code/load_data/load_variants.R"))

# Read in alias conversion file
unimap = read.table(paste0(base_path, "/input/aliases/UniMap.txt"),header = TRUE, sep = "\t", stringsAsFactors = FALSE)
quant_unimap = read.table(paste0(base_path, "/input/aliases/combined_unimap.txt"),header = TRUE, sep = "\t", stringsAsFactors = FALSE)
quant_unimap <- unique(quant_unimap)


######################## Permutation test: Corum-expanded ########################

n_perm = 1000

for(int_type in c("GATA4", "NKX25", "TBX5")){
  
  #### Saintq_norm ####
  # Load in GATA4 scores (saintq protein scoring)
  if (int_type == "GATA4"){
    tab = read.table(file = paste0(base_path, "/intermediate/interaction_scoring/saintq/GATA4_protein_norm_true.txt"), 
                         sep = "\t", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, comment.char = "")
  } else{
    tab = read.table(file = paste0(base_path,"/intermediate/interaction_scoring/saintq/", int_type, "_peptide_norm_true.txt"),
                         sep = "\t", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, comment.char = "")
  }

  #### original interactomes ####
  #ints = read.table(paste0(base_path, "/input/precomp_interactomes/interactomes.tsv"), 
  #                   sep = "\t",fill = TRUE, header = TRUE, stringsAsFactors = FALSE)
  #tab = ints[which(ints$Bait == int_type),]
  
  # Choose significant genes
  #tab_sig = subset(tab, tab$BFDR < 0.05)
  if(int_type == "NKX25"){
    tab_sig = subset(tab, tab$BFDR < 0.1)
  } else if (int_type == "TBX5"){
    tab_sig = subset(tab, tab$BFDR < 0.05)
  } else if (int_type == "GATA4") {
    tab_sig = subset(tab, tab$BFDR < 0.001)
  }
  
  # Remove non-nuclear genes from interactome
  names(tab_sig)[which(names(tab_sig) == "Prey_proteinname")] <- "Prey"
  non_nuclear = read.table(paste0(base_path, "/input/databases/localization_non_nuclear_BINGO_.txt"),
                           sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  non_nuclear = non_nuclear[which(!non_nuclear$subcellular.location == "Nucleus"),]
  tab_sig = tab_sig[which(!tab_sig$Prey %in% non_nuclear$UniprotID),]
  
  #### Corum expansion ####
  #complexes = get_comps(tab_sig$Prey,corum, enriched=FALSE)
  #complexes = get_comps(tab_sig$Prey, corum, enriched=TRUE)
  complexes = get_comps("NA", corum, enriched=TRUE)
  interactome_prots = unique(c(get_prots(complexes, corum), tab_sig$Prey))
  
  # Get gene names for interactome
  interactome_genes = quant_unimap$GeneSymbol[match(interactome_prots, quant_unimap$UniProt)]
  interactome_genes = unlist(strsplit(interactome_genes, "; "))
  
  # Remove blacklist genes from interactome ####NOTE: Should this be before corum expansion?
  blacklist = read.table(paste0(base_path, "/intermediate/rnaseq/", int_type, "_negativeFC_blacklist.txt"))
  geneList <- interactome_genes[which(!interactome_genes %in% blacklist)]
  
  # Create plots and permutation records
  for(mut_type in c("DNV", "LoF", "syn-DNV")){
    
    if(mut_type == "DNV"){cases = DNV_cases; ctrls = DNV_ctrls
    } else if(mut_type == "LoF"){cases = LoF_cases; ctrls = LoF_ctrls
    } else if(mut_type == "syn-DNV"){cases = DNV_case_syn; ctrls = DNV_ctrl_syn
    }
    
    true_OR = get_OR(cases, ctrls, geneList, fisher = TRUE, int_type=int_type, mut_type)
    permList <- permute_status(cases, ctrls, geneList, n_perm)
    img = perm_viz(permList, true_OR, mut_type, int_type, n_tests=3)
    ggsave(file= paste0(out_path, "/", int_type, "_", mut_type, ".svg"), plot = img)
    
    # Save out record of permutation data
    true_row = as.data.frame(true_OR)
    true_row$type = "true"
    perm_rows = as.data.frame(permList)
    perm_rows$type = "permuted"
    names(true_row) = c("OddsRatio", "type")
    names(perm_rows) = c("OddsRatio", "type")
    perm_record = rbind(true_row, perm_rows)
    
    write.table(perm_record, 
                file = paste0(out_path, "/perm_score_lists/", int_type, "_",mut_type,".txt"),
                row.names = FALSE, quote = FALSE)
    
    
  }
}

# # # # # # GATA4-TBX5 combined interactome # # # # # # # 
int_type = "GATA4-TBX5 combined"
n_perm = 1000
  
#### Saintq_norm ####
# Load in GATA4 scores (saintq protein scoring)
tabG = read.table(file = paste0(base_path, "/intermediate/interaction_scoring/saintq/GATA4_protein_norm_true.txt"), 
                 sep = "\t", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, comment.char = "")
tabT = read.table(file = paste0(base_path,"/intermediate/interaction_scoring/saintq/TBX5_peptide_norm_true.txt"),
                  sep = "\t", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, comment.char = "")

#### original interactomes ####
#ints = read.table(paste0(base_path, "/input/precomp_interactomes/interactomes.tsv"), 
#                   sep = "\t",fill = TRUE, header = TRUE, stringsAsFactors = FALSE)
#tab = ints[which(ints$Bait == "GATA4" | ints$Bait == "TBX5"),]

# Choose significant genes
#tab_sig = subset(tab, tab$BFDR < 0.05)
tab_sigG = subset(tabG, tabG$BFDR < 0.001)
#tab_sigG = subset(tabG, tabG$BFDR < 0.05)
tab_sigT = subset(tabT, tabT$BFDR < 0.05)
cols = c("Bait", "Prey", "BFDR")
tab_sig = rbind(tab_sigG[,cols], tab_sigT[,cols])

# Remove non-nuclear genes from interactome
names(tab_sig)[which(names(tab_sig) == "Prey_proteinname")] <- "Prey"
non_nuclear = read.table(paste0(base_path, "/input/databases/localization_non_nuclear_BINGO_.txt"),
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE)
non_nuclear = non_nuclear[which(!non_nuclear$subcellular.location == "Nucleus"),]
tab_sig = tab_sig[which(!tab_sig$Prey %in% non_nuclear$UniprotID),]

#### Corum expansion ####
#complexes = get_comps(tab_sig$Prey,corum, enriched=FALSE)
#complexes = get_comps(tab_sig$Prey, corum, enriched=TRUE)
complexes = get_comps("NA", corum, enriched=TRUE)
interactome_prots = unique(c(get_prots(complexes, corum), tab_sig$Prey))

# Get gene names for interactome
interactome_genes = quant_unimap$GeneSymbol[match(interactome_prots, quant_unimap$UniProt)]
interactome_genes = unlist(strsplit(interactome_genes, "; "))

# Create plots and permutation records
for(mut_type in c("DNV", "LoF", "syn-DNV")){
  
  if(mut_type == "DNV"){cases = DNV_cases; ctrls = DNV_ctrls
  } else if(mut_type == "LoF"){cases = LoF_cases; ctrls = LoF_ctrls
  } else if(mut_type == "syn-DNV"){cases = DNV_case_syn; ctrls = DNV_ctrl_syn
  }
  
  true_OR = get_OR(cases, ctrls, interactome_genes, fisher = TRUE, int_type=int_type, mut_type)
  permList <- permute_status(cases, ctrls, interactome_genes, n_perm)
  img = perm_viz(permList, true_OR, mut_type, int_type, n_tests=3)
  ggsave(file= paste0(out_path, "/", int_type, "_", mut_type, ".svg"), plot = img)
  
  # Save out record of permutation data
  true_row = as.data.frame(true_OR)
  true_row$type = "true"
  perm_rows = as.data.frame(permList)
  perm_rows$type = "permuted"
  names(true_row) = c("OddsRatio", "type")
  names(perm_rows) = c("OddsRatio", "type")
  perm_record = rbind(true_row, perm_rows)
  
  write.table(perm_record, 
              file = paste0(out_path, "/perm_score_lists/", int_type, "_",mut_type,".txt"),
              row.names = FALSE, quote = FALSE)
  
  
}

# # # # # # new GATA4 WT runs # # # # # # # 
int_type = "wt_gata4"
  
#### Saintq_norm ####
tab = read.table(paste0(base_path, "/intermediate/interaction_scoring/saintq/", 
                    "wt_gata4_protein_norm_true.txt"),
             sep = "\t", stringsAsFactors=F, comment.char="", header=T)
tab_sig = subset(tab, tab$BFDR < 0.05)

# Remove non-nuclear genes from interactome
names(tab_sig)[which(names(tab_sig) == "Prey_proteinname")] <- "Prey"
non_nuclear = read.table(paste0(base_path, "/input/databases/localization_non_nuclear_BINGO_.txt"),
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE)
non_nuclear = non_nuclear[which(!non_nuclear$subcellular.location == "Nucleus"),]
tab_sig = tab_sig[which(!tab_sig$Prey %in% non_nuclear$UniprotID),]

#### Corum expansion ####
#complexes = get_comps(tab_sig$Prey,corum, enriched=FALSE)
#complexes = get_comps(tab_sig$Prey, corum, enriched=TRUE)
complexes = get_comps("NA", corum, enriched=TRUE)
interactome_prots = unique(c(get_prots(complexes, corum), tab_sig$Prey))

# Get gene names for interactome
interactome_genes = quant_unimap$GeneSymbol[match(interactome_prots, quant_unimap$UniProt)]
interactome_genes = unlist(strsplit(interactome_genes, "; "))
write.table(interactome_genes[!is.na(interactome_genes)],
            file = paste0(base_path, "/intermediate/interactome_lists/saintq_n/new_gata4_wt.txt"),
            quote = F, row.names=F, col.names = F)

# Remove blacklist genes from interactome ####NOTE: Should this be before corum expansion?
#blacklist = read.table(paste0(base_path, "/intermediate/rnaseq/", int_type, "_negativeFC_blacklist.txt"))
#geneList <- interactome_genes[which(!interactome_genes %in% blacklist)]
geneList <- interactome_genes

# Create plots and permutation records
for(mut_type in c("DNV", "LoF", "syn-DNV")){
  
  if(mut_type == "DNV"){cases = DNV_cases; ctrls = DNV_ctrls
  } else if(mut_type == "LoF"){cases = LoF_cases; ctrls = LoF_ctrls
  } else if(mut_type == "syn-DNV"){cases = DNV_case_syn; ctrls = DNV_ctrl_syn
  }
  
  true_OR = get_OR(cases, ctrls, geneList, fisher = TRUE, int_type=int_type, mut_type)
  permList <- permute_status(cases, ctrls, geneList, n_perm)
  img = perm_viz(permList, true_OR, mut_type, int_type, n_tests=3)
  ggsave(file= paste0(out_path, "/", int_type, "_", mut_type, ".svg"), plot = img)
  
  # Save out record of permutation data
  true_row = as.data.frame(true_OR)
  true_row$type = "true"
  perm_rows = as.data.frame(permList)
  perm_rows$type = "permuted"
  names(true_row) = c("OddsRatio", "type")
  names(perm_rows) = c("OddsRatio", "type")
  perm_record = rbind(true_row, perm_rows)
  
  write.table(perm_record, 
              file = paste0(out_path, "/perm_score_lists/", int_type, "_",mut_type,".txt"),
              row.names = FALSE, quote = FALSE)
  
  
}

# # # # # # combining GATA4 WT runs # # # # # # # 

int_type = "wt_gata4_combined"

#### Saintq_norm ####
tab = read.table(paste0(base_path, "/intermediate/interaction_scoring/saintq/",
                        "wt_gata4_combined_protein_norm_true_wtcombined.txt"),
                 sep = "\t", stringsAsFactors=F, comment.char="", header=T)
tab_sig = subset(tab, tab$BFDR < 0.05)

# Remove non-nuclear genes from interactome
names(tab_sig)[which(names(tab_sig) == "Prey_proteinname")] <- "Prey"
non_nuclear = read.table(paste0(base_path, "/input/databases/localization_non_nuclear_BINGO_.txt"),
                         sep = "\t", header = TRUE, stringsAsFactors = FALSE)
non_nuclear = non_nuclear[which(!non_nuclear$subcellular.location == "Nucleus"),]
tab_sig = tab_sig[which(!tab_sig$Prey %in% non_nuclear$UniprotID),]

#### Corum expansion ####
#complexes = get_comps(tab_sig$Prey,corum, enriched=FALSE)
#complexes = get_comps(tab_sig$Prey, corum, enriched=TRUE)
complexes = get_comps("NA", corum, enriched=TRUE)
interactome_prots = unique(c(get_prots(complexes, corum), tab_sig$Prey))

# Get gene names for interactome
interactome_genes = quant_unimap$GeneSymbol[match(interactome_prots, quant_unimap$UniProt)]
interactome_genes = unlist(strsplit(interactome_genes, "; "))
write.table(interactome_genes[!is.na(interactome_genes)],
            file = paste0(base_path, "/intermediate/interactome_lists/saintq_n/gata4_wtcombined.txt"),
            quote = F, row.names=F, col.names = F)

# Remove blacklist genes from interactome ####NOTE: Should this be before corum expansion?
#blacklist = read.table(paste0(base_path, "/intermediate/rnaseq/", int_type, "_negativeFC_blacklist.txt"))
#geneList <- interactome_genes[which(!interactome_genes %in% blacklist)]
geneList <- interactome_genes

# Create plots and permutation records
for(mut_type in c("DNV", "LoF", "syn-DNV")){
  
  if(mut_type == "DNV"){cases = DNV_cases; ctrls = DNV_ctrls
  } else if(mut_type == "LoF"){cases = LoF_cases; ctrls = LoF_ctrls
  } else if(mut_type == "syn-DNV"){cases = DNV_case_syn; ctrls = DNV_ctrl_syn
  }
  
  true_OR = get_OR(cases, ctrls, geneList, fisher = TRUE, int_type=int_type, mut_type)
  permList <- permute_status(cases, ctrls, geneList, n_perm)
  img = perm_viz(permList, true_OR, mut_type, int_type, n_tests=3)
  ggsave(file= paste0(out_path, "/", int_type, "_", mut_type, ".svg"), plot = img)
  
  # Save out record of permutation data
  true_row = as.data.frame(true_OR)
  true_row$type = "true"
  perm_rows = as.data.frame(permList)
  perm_rows$type = "permuted"
  names(true_row) = c("OddsRatio", "type")
  names(perm_rows) = c("OddsRatio", "type")
  perm_record = rbind(true_row, perm_rows)
  
  write.table(perm_record, 
              file = paste0(out_path, "/perm_score_lists/", int_type, "_",mut_type,".txt"),
              row.names = FALSE, quote = FALSE)
  
  
}
