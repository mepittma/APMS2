# This file creates permutation visualizations for lab meeting July24.
base_path = "/Users/student/Documents/PollardLab/APMS2"
out_path = paste0(base_path, "/output/july22_test")

known_genes = read.csv("/Users/student/Documents/PollardLab/PCGC/input/known_genes.txt", stringsAsFactors = F,
                       header = F)
known_genes = known_genes$V1

# Load in the variants being tested in this file
source(paste0(base_path, "/code/load_data/load_variants.R"))

# Overwrite LoF cases/ctrls; load rare inherited synonymous
LoF_cases = read.csv(paste0(out_path, "/variants/WES_LoF_cases.csv"), stringsAsFactors = F)
LoF_ctrls = read.csv(paste0(out_path, "/variants/WES_LoF_ctrls.csv"), stringsAsFactors = F)

isyn_cases = read.csv(paste0(out_path, "/variants/WES_inhsyn_cases.csv"), stringsAsFactors = F)
isyn_ctrls = read.csv(paste0(out_path, "/variants/WES_inhsyn_ctrls.csv"), stringsAsFactors = F)

# Select a random subset of inherited synonymous variants (100,000 is a lot to work with for prelim exploration)
isyn_cases = isyn_cases[sample(nrow(isyn_cases), nrow(LoF_cases)),]
isyn_ctrls = isyn_ctrls[sample(nrow(isyn_ctrls), nrow(LoF_ctrls)),]

# Read in alias conversion file
unimap = read.table(paste0(base_path, "/input/aliases/UniMap.txt"),header = TRUE, sep = "\t", stringsAsFactors = FALSE)
quant_unimap = read.table(paste0(base_path, "/input/aliases/combined_unimap.txt"),header = TRUE, sep = "\t", stringsAsFactors = FALSE)
quant_unimap <- unique(quant_unimap)


######################## Permutation test ########################
n_perm = 1000
int_type = "GATA4-TBX5"

# Load in GATA4 scores (saintq protein scoring)
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
interactome_prots = unique(tab_sig$Prey)

# Get gene names for interactome
interactome_genes = quant_unimap$GeneSymbol[match(interactome_prots, quant_unimap$UniProt)]
interactome_genes = unlist(strsplit(interactome_genes, "; "))

source(paste0(base_path, "/code/functions/permutation_functions.R"))
#get_fisher, get_OR, get_prots, get_sig_compl, perm_viz, permute_status

# Create plots and permutation records
n_perm=1000
for(int_type in c("all_genes", "unknown_genes")){
  
  if (int_type == "unknown_genes"){
    interactome_genes = interactome_genes[which(!interactome_genes %in% known_genes)]
  }
  
  #for(mut_type in c("LoF","syn","DNV","synDNV")){
  for(mut_type in c("DNV","synDNV")){
    if(mut_type == "DNV"){cases = DNV_cases; ctrls = DNV_ctrls
    } else if(mut_type == "LoF"){cases = LoF_cases; ctrls = LoF_ctrls
    } else if(mut_type == "syn"){cases = isyn_cases; ctrls = isyn_ctrls
    } else if(mut_type == "synDNV"){cases = DNV_case_syn; ctrls = DNV_ctrl_syn}
    
    true_OR = get_OR(cases, ctrls, interactome_genes, fisher = TRUE, int_type=int_type, mut_type)
    permList <- permute_status(cases, ctrls, gene_list=interactome_genes, n_perm=n_perm)
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

# For July 24 lab meeting, we want to see permutations for GATA4-TBX5 interactomes combined, both with and
# without known CHD-associated genes. I suck at life so we're re-creating them from perm_score_lists files.
gnum = "all"
gnum = "unknown"
mut_type = "LoF"
mut_type = "syn"
int_type = "GATA4-TBX5"

perms = read.table(paste0(base_path, "/output/july22_test/perm_score_lists/", gnum, "_genes_", mut_type, ".txt"),
                       stringsAsFactors = F, header = T)
perm_list = perms$OddsRatio[which(perms$type == "permuted")]
true_OR = perm_list[1,1]

  
nGreater = length(perm_list[which(perm_list >= true_OR)])
pval = (nGreater/length(perm_list))
if(nGreater == 0 ){
  pval = "< 0.001"
}
  
  
par(lwd=2)
h = hist(perm_list)
#pdf(paste0(out_path, "/", int_type, "_", mut_type, "_",gnum,"_genes.pdf"))
plot(h, main=paste0(mut_type, 
                    " mutations in ", int_type, " interactome"), 
     sub = paste0("p: ", pval),
     xlab = "Odds Ratio", lwd=2)
abline(v = 0.75, lty="dotted", lwd=5, col="red")
dev.off()
  


