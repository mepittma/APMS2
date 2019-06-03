base_path = "/Users/student/Documents/PollardLab/APMS2"
out_path = paste0(base_path, "/output/permutations/geneSet_saintq_n_unexpanded_G001_T05_N1")

library(ggplot2)
source(paste0(base_path, "/code/functions/permutation_functions.R"))
source(paste0(base_path, "/code/functions/score_interactomes.R"))

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

######################## Establish background ########################
background = read.table(paste0(base_path, "/input/databases/HGNC.txt"), sep = "\t",
                              header = TRUE, stringsAsFactors = FALSE)
background_genes = background$Approved.symbol

######################## Permutation test ########################

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
  
  # Get list of all possible prey names
  tab$Prey_geneName = quant_unimap$GeneSymbol[match(tab$Prey, quant_unimap$UniProt)]
  
  # Choose significant genes
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
  
  # Get gene names for interactome
  interactome_genes = tab_sig$Prey_geneName
  interactome_genes = unlist(strsplit(interactome_genes, "; "))
  
  # Remove blacklist genes from interactome
  blacklist = read.table(paste0(base_path, "/intermediate/rnaseq/", int_type, "_negativeFC_blacklist.txt"))
  geneList <- interactome_genes[which(!interactome_genes %in% blacklist)]
  background_genes <- background_genes[which(!background_genes %in% blacklist)]
  
  # Create plots and permutation records
  for(mut_type in c("DNV", "LoF", "syn-DNV")){
    
    if(mut_type == "DNV"){cases = DNV_cases; ctrls = DNV_ctrls
    } else if(mut_type == "LoF"){cases = LoF_cases; ctrls = LoF_ctrls
    } else if(mut_type == "syn-DNV"){cases = DNV_case_syn; ctrls = DNV_ctrl_syn
    }
    
    true_OR = get_OR(cases, ctrls, geneList, fisher = TRUE, int_type=int_type, mut_type)
    permList <- permute_geneList(cases, ctrls, geneList, background_genes, n_perm)
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