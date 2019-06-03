######################## Setup ######################## 
base_path = "/Users/student/Documents/PollardLab/APMS2"
#out_path = paste0(base_path, "/output/permutations/saintq_n")
#out_path = paste0(base_path, "/output/permutations/original")
#out_path = paste0(base_path, "/output/permutations/saintq_n_cutoff_G02_T02_N1")
out_path = paste0(base_path, "/output/permutations/saintq_n_cutoff_G001_T05_N1")
#out_path = paste0(base_path, "/output/permutations/saintq_n_enrichedComp")
#out_path = paste0(base_path, "/output/permutations/original_enrichedComp")

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

################
n_perm = 1000
int_type = "GATA4-TBX5-combined"

for(analysis_type in c("original", "saintq_n")){
  out_path = paste0(base_path, "/output/permutations/", analysis_type, "_cutoff_G001_T05_N1")
  
  tab_sig = read.csv(file = paste0(base_path, "/intermediate/interactome_lists/",
                                 analysis_type, "/combined_APMS_interactome_G001_T05_N1.csv"),
                 header = TRUE, stringsAsFactors = FALSE)
  
  #### Corum expansion ####
  complexes = get_comps(tab_sig$Prey,corum, enriched=FALSE)
  #complexes = get_comps(tab_sig$Prey, corum, enriched=TRUE)
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
}



