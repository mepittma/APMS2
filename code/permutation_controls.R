######################## Setup ######################## 
base_path = "/Users/student/Documents/PollardLab/APMS2"
#out_path = paste0(base_path, "/output/permutations/controls_expanded")
out_path = paste0(base_path, "/output/permutations/controls_unexpanded")

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


# Load in interactome data.
ints <- read.csv(paste0(base_path, "/input/precomp_interactomes/control_interactomes.csv"), 
                 fill = TRUE, header = TRUE, stringsAsFactors = FALSE)

# Save out protein IDs of any interacting corum members
## 02/11/19 - This was translated using an online web service 
# (bioDBnet) into GeneSymbols, then appended to the unimap file used below
#write.table(unique(c(int_comp$Interactor1, int_comp$Interactor2)), 
#            file = paste0(base_path, "/input/interactomes/control_toMap.txt"),
#            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Get gene names of corum complex members
unimap <- read.table(paste0(base_path, "/input/aliases/UniMap.txt"),header = TRUE, sep = "\t", fill = TRUE)
unimap <- unique(unimap)

# Permutation analysis for corum-expanded interactomes
for(int_type in c("FOXA2", "GATA1", "GATA2", "GATA3")){
  
  tab_sig = ints[which(ints$Bait == int_type),]

  #int_conv <- read.table(paste0(base_path, "/input/databases/control_UPaccessions.txt"), 
  #                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  #int_prots <- unlist(strsplit(
  #  int_conv$UniProt.Accession[which(int_conv$Gene.Symbol.and.Synonyms %in% tab_sig$Prey_genename)],
  #  "; "))
  
  #complexes = get_comps(int_prots,corum, enriched=FALSE)
  #interactome_prots = unique(c(get_prots(complexes, corum), tab_sig$Prey))
  
  # Get gene names for interactome
  #interactome_genes = quant_unimap$GeneSymbol[match(interactome_prots, quant_unimap$UniProt)]
  interactome_genes = tab_sig$Prey_genename
  geneList = unlist(strsplit(interactome_genes, "; "))
  
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

######### HEK cell control interactomes #########

# Load in interactome data.
ints <- read.csv(paste0(base_path, "/input/precomp_interactomes/PPG_TF_APMS.csv"), 
                 fill = TRUE, header = TRUE, stringsAsFactors = FALSE)
names(ints) = c("Bait", "Prey_proteinname", "MIST", "SAINT","etc","etc2","Prey_genename")

#int_prots = unique(c(ints$Prey_proteinname))
# Find corum complexes with members in the interactomes.
#int_complexes <- unique(c(as.vector(corum$Complex.id[which(
#  corum$Interactor1 %in% as.vector(int_prots)| corum$Interactor2 %in% as.vector(int_prots))])))
#int_comp <- corum[which(corum$Complex.id %in% int_complexes),]

# Save out protein IDs of any interacting corum members
## 02/21/19 - This was translated using an online web service 
# (bioDBnet) into GeneSymbols, then appended to the unimap file used below
#write.table(unique(c(int_comp$Interactor1, int_comp$Interactor2, int_prots)), 
#            file = paste0(base_path, "/input/interactomes/PPG_toMap.txt"),
#            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Get gene names of corum complex members
#unimap <- read.table(paste0(base_path, "/input/UniMap.txt"),header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
#unimap <- unique(unimap)
#int_comp$geneName1 <- unimap$GeneSymbol[match(int_comp$Interactor1, as.vector(unimap$UniProt))]
#int_comp$geneName2 <- unimap$GeneSymbol[match(int_comp$Interactor2, as.vector(unimap$UniProt))]
#ints$Prey_genename <- unimap$GeneSymbol[match(ints$Prey_proteinname, as.vector(unimap$UniProt))]

#int_names <- unique(c(int_comp$geneName1, int_comp$geneName2, ints$Prey_genename))
#mut_names <- unique(c(LoF_cases$Gene, LoF_ctrls$Gene, DNV_cases$Gene, DNV_ctrls$Gene,
#                      rec_cases$Gene, rec_ctrls$Gene))
#int_names[!(int_names %in% mut_names)]

#int_types = unique(ints$Bait)
int_types = c("GATA4", "TBX5")

# Run permutation
n_perm = 1000
for (int_type in int_types){
  
  #compl <- unique(int_comp$Complex.id[which(int_comp$Interactor1 %in% ints$Prey_proteinname[which(ints$Bait == int_type)]  |
  #                                            int_comp$Interactor2 %in% ints$Prey_proteinname[which(ints$Bait == int_type)])])
  #geneList = unique(c(ints$Prey_genename[which(ints$Bait == int_type)], int_type, 
  #                    as.vector(int_comp$geneName1[which(int_comp$Complex.id %in% compl)]),
  #                    as.vector(int_comp$geneName2[which(int_comp$Complex.id %in% compl)])))
  sig_tab = ints[which(ints$Bait == int_type),]
  geneList = sig_tab$Prey_genename
  
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

# Write out the genes that have mutations 
for (int_type in int_types){
#  compl <- unique(int_comp$Complex.id[which(int_comp$Interactor1 %in% ints$Prey_proteinname[which(ints$Bait == int_type)]  |
#                                              int_comp$Interactor2 %in% ints$Prey_proteinname[which(ints$Bait == int_type)])])
#  geneList = unique(c(ints$Prey_genename[which(ints$Bait == int_type)], int_type, 
#                      as.vector(int_comp$geneName1[which(int_comp$Complex.id %in% compl)]),
#                      as.vector(int_comp$geneName2[which(int_comp$Complex.id %in% compl)])))
  
  sig_tab = ints[which(ints$Bait == int_type),]
  geneList = sig_tab$Prey_genename
  
  mutated = c()
  
  for(i in 1:nrow(DNV_cases)){
    if(DNV_cases[i,'Gene'] %in% geneList){
      
      mutated = c(mutated, DNV_cases[i, 'Gene'])
      
    }
  }
  mut_tab = table(mutated)
  write.csv(mut_tab, paste0(out_path, "/mutation_lists/", int_type, "_mutationCounts.csv"), row.names = FALSE)
  
}

#########################
n_perm = 1000
out_path = paste0(base_path, "/output/permutations/bgt_controls")
# Load in interactome data.
ints <- read.table(paste0(base_path, "/input/precomp_interactomes/bgt_TBX5_GATA4.tsv"), 
                 fill = TRUE, header = TRUE, stringsAsFactors = FALSE, sep = "\t", comment.char = "")
for(int_type in c("gata4", "tbx5")){
  if(int_type=="tbx5"){cutoff=0.01
  } else if (int_type == "gata4"){cutoff = 0.05}
  tab_sig = ints[which(ints$Bait == int_type & ints$BFDR <= cutoff),]
  
  interactome_genes = quant_unimap$GeneSymbol[match(tab_sig$Prey, quant_unimap$UniProt)]
  geneList = unlist(strsplit(interactome_genes, "; "))
  
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
  
# Write out the genes that have mutations 
for (int_type in c("gata4", "tbx5")){
  #  compl <- unique(int_comp$Complex.id[which(int_comp$Interactor1 %in% ints$Prey_proteinname[which(ints$Bait == int_type)]  |
  #                                              int_comp$Interactor2 %in% ints$Prey_proteinname[which(ints$Bait == int_type)])])
  #  geneList = unique(c(ints$Prey_genename[which(ints$Bait == int_type)], int_type, 
  #                      as.vector(int_comp$geneName1[which(int_comp$Complex.id %in% compl)]),
  #                      as.vector(int_comp$geneName2[which(int_comp$Complex.id %in% compl)])))
  
  if(int_type=="tbx5"){cutoff=0.01
  } else if (int_type == "gata4"){cutoff = 0.05}
  tab_sig = ints[which(ints$Bait == int_type & ints$BFDR <= cutoff),]
  interactome_genes = quant_unimap$GeneSymbol[match(tab_sig$Prey, quant_unimap$UniProt)]
  geneList = unlist(strsplit(interactome_genes, "; "))
  
  mutated = c()
  
  for(i in 1:nrow(DNV_cases)){
    if(DNV_cases[i,'Gene'] %in% geneList){
      
      mutated = c(mutated, DNV_cases[i, 'Gene'])
      
    }
  }
  mut_tab = table(mutated)
  write.csv(mut_tab, paste0(out_path, "/mutation_lists/", int_type, "_mutationCounts.csv"), row.names = FALSE)
  
}
