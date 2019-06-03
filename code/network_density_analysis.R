base_path = "/Users/student/Documents/PollardLab/APMS2"
out_path = paste0(base_path, "/output/network_density")
library(ggplot2)

# Read in variant data
source(paste0(base_path, "/code/load_data/load_variants.R"))

# Read in PPI table
PPI = read.table(paste0(base_path, "/input/databases/iRefIndex.txt"), 
                 sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# # # # # # # # # #
# Functions

# Calculate network density
get_density <- function(sub_net){
  n_nodes = length(unique(c(sub_net$aliasA, sub_net$aliasB)))
  pc = (n_nodes*(n_nodes-1))/2
  ac = nrow(sub_net)
  return(ac/pc)
}

# Permute subnetwork to create a distribution of random subnetwork densities
permute_subnet <- function(ppi, n_nodes, n_perm){
  
  perm_list = c()
  unique_nodes = unique(c(ppi$aliasA, ppi$aliasB))
  
  for (i in 1:n_perm){
    
    # Choose random nodes to create a network of the same size
    #### ISSUE: WHY TF ARE RANDOM_SUB_NODES THE SAME EVERY TIME???
    random_sub_nodes = sample(unique_nodes, n_nodes)
    random_sub = ppi[which(ppi$aliasA %in% random_sub_nodes | ppi$aliasB %in% random_sub_nodes),]
    nrow(random_sub)
    
    # Find density of the random network
    perm_list = c(perm_list, get_density(random_sub))
    
  }
  return(perm_list)
}

# Function to visualize a list of permuted odds ratios compared to the true odds ratio
perm_viz <- function(perm_list, true_OR, mut_type, int_type, n_tests){
  
  nGreater = length(perm_list[which(perm_list >= true_OR)])
  pval = (nGreater/length(perm_list))
  if(nGreater == 0 ){
    pval = "< 0.001"
  }
  
  pdf(paste0(out_path, "/", int_type, "_", mut_type, ".pdf"))
  hist(perm_list, main=paste0("Subnetwork permutation of ",mut_type, 
                              " mutations in ", int_type, "s"), 
       sub = paste0("p: ", pval),
       xlab = "Network Density")
  abline(v = true_OR, lty="dotted", lwd="5", col = "red")
  dev.off()
  
  df = as.data.frame(perm_list)
  img <- ggplot(data=df, aes(df$perm_list)) + geom_histogram() + 
    geom_vline(xintercept = true_OR, col = "red", linetype="dashed") +
    labs(title=paste0("Subnetwork permutation of ",mut_type," mutations in ", int_type, "s"),
         subtitle = paste0("p = ", pval),
         x="Network Density", y = "Frequency")
  
  return(img)
}

# # # # # # # # # # ARE VARIANTS SIGNIFICANTLY CONNECTED? # # # # # # # # # # 

# Check to see if any genes are not represented in the PPI table
all_genes = unique(c(DNV_case_syn$Gene, DNV_ctrl_syn$Gene, DNV_cases$Gene, DNV_ctrls$Gene,
                     LoF_cases$Gene, LoF_ctrls$Gene, rec_cases$Gene, rec_ctrls$Gene))
length(all_genes[which(!all_genes %in% unique(c(PPI$aliasA, PPI$aliasB)))])
#2219 mutated genes do not have entries in iRefIndex

# Number of total connections
n_conn_total = nrow(PPI)

# For each type of mutation, test if a random network of n nodes is as dense as the true
for(mut_type in c("DNV", "LoF", "syn-DNV")){
  
  if(mut_type == "DNV"){cases = DNV_cases; ctrls = DNV_ctrls
  } else if(mut_type == "LoF"){cases = LoF_cases; ctrls = LoF_ctrls
  } else if(mut_type == "syn-DNV"){cases = DNV_case_syn; ctrls = DNV_ctrl_syn
  }
  
  for (status in c("case","ctrl")){
    if(status == "case"){cases = cases
    } else if(status == "ctrl"){cases = ctrls
    }
    
    # Compare true variant subnetwork to randomly-permuted subnetworks
    sub_net = PPI[which(PPI$aliasA %in% cases$Gene | PPI$aliasB %in% cases$Gene),]
    true_density = get_density(sub_net)
    
    n_perm = 1000
    n_nodes = length(unique(cases$Gene))
    permList = permute_subnet(PPI, n_nodes, n_perm)
    img = perm_viz(permList, true_density, mut_type, status, n_tests=6)
    ggsave(file= paste0(out_path, "/", status, "_", mut_type, ".svg"), plot = img)
  }
  

}

# ODD RESULTS...
# Why would synonymous DNVs occur in nodes that are more connected in cases over controls?

# # # # # # # # # # ARE GATA4 AND TBX5 SIGNIFICANTLY CONNECTED? # # # # # # # # # # 

# GATA4 and TBX5 alone
quant_unimap = read.table(paste0(base_path, "/input/aliases/combined_unimap.txt"),header = TRUE, sep = "\t", stringsAsFactors = FALSE)
quant_unimap <- unique(quant_unimap)
for(int_type in c("GATA4","TBX5")){
  
  if (int_type == "GATA4"){
    tab = read.table(file = paste0(base_path, "/intermediate/interaction_scoring/saintq/GATA4_protein_norm_true.txt"), 
                     sep = "\t", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, comment.char = "")
  } else{
    tab = read.table(file = paste0(base_path,"/intermediate/interaction_scoring/saintq/", int_type, "_peptide_norm_true.txt"),
                     sep = "\t", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, comment.char = "")
  }
  
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
  
  # Get gene names for interactome
  interactome_genes = quant_unimap$GeneSymbol[match(tab_sig$Prey, quant_unimap$UniProt)]
  interactome_genes = unlist(strsplit(interactome_genes, "; "))
  
  # Remove blacklist genes from interactome ####NOTE: Should this be before corum expansion?
  blacklist = read.table(paste0(base_path, "/intermediate/rnaseq/", int_type, "_negativeFC_blacklist.txt"))
  geneList <- interactome_genes[which(!interactome_genes %in% blacklist)]
  
  # Permutation: are these interactome genes more connected than other random subnetworks?
  sub_net = PPI[which(PPI$aliasA %in% geneList | PPI$aliasB %in% geneList),]
  true_density = get_density(sub_net)
  
  n_perm = 100
  n_nodes = length(unique(geneList))
  permList = permute_subnet(PPI, n_nodes, n_perm)
  img = perm_viz(permList, true_density, int_type, "network_density", n_tests=6)
  ggsave(file= paste0(out_path, "/network_density_", int_type, ".svg"), plot = img)

}

# GATA4 and TBX5 together
int_type = "GATA4-TBX5"
tabG = read.table(file = paste0(base_path, "/intermediate/interaction_scoring/saintq/GATA4_protein_norm_true.txt"), 
                  sep = "\t", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, comment.char = "")
tabT = read.table(file = paste0(base_path,"/intermediate/interaction_scoring/saintq/TBX5_peptide_norm_true.txt"),
                  sep = "\t", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, comment.char = "")

# Choose significant genes
tab_sigG = subset(tabG, tabG$BFDR < 0.001)
tab_sigT = subset(tabT, tabT$BFDR < 0.05)
cols = c("Bait", "Prey", "BFDR")
tab_sig = rbind(tabG[,cols], tabT[,cols])
tab_sig = tab_sig[which(!tab_sig$Prey %in% non_nuclear$UniprotID),]

# Get gene names for interactome
interactome_genes = quant_unimap$GeneSymbol[match(tab_sig$Prey, quant_unimap$UniProt)]
interactome_genes = unlist(strsplit(interactome_genes, "; "))

# Remove blacklist genes from interactome
geneList <- interactome_genes[which(!interactome_genes %in% blacklist)]

# Permutation: are these interactome genes more connected than other random subnetworks?
sub_net = PPI[which(PPI$aliasA %in% geneList | PPI$aliasB %in% geneList),]
true_density = get_density(sub_net)

n_perm = 100
n_nodes = length(unique(geneList))
permList = permute_subnet(PPI, n_nodes, n_perm)
img = perm_viz(permList, true_density, int_type, "network_density", n_tests=6)
ggsave(file= paste0(out_path, "/network_density_", int_type, ".svg"), plot = img)

