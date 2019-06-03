base_path = "/Users/student/Documents/PollardLab/APMS2"
out_path = paste0(base_path, "/output/network_density")
library(plyr)
library(dplyr)
library(ggplot2)

# Read in variant data
source(paste0(base_path, "/code/load_data/load_variants.R"))

# Read in PPI table
PPI = read.table(paste0(base_path, "/input/databases/iRefIndex.txt"), 
                 sep = "\t", header = TRUE, stringsAsFactors = FALSE)


# Check to see if any genes are not represented in the PPI table
all_genes = unique(c(DNV_case_syn$Gene, DNV_ctrl_syn$Gene, DNV_cases$Gene, DNV_ctrls$Gene,
                     LoF_cases$Gene, LoF_ctrls$Gene, rec_cases$Gene, rec_ctrls$Gene))
length(all_genes[which(!all_genes %in% unique(c(PPI$aliasA, PPI$aliasB)))])
#2219 mutated genes do not have entries in iRefIndex

# # # # # # # # # # CREATE EDGE FILE # # # # # # # # # # 
case_LoF_genes = unique(c(LoF_cases$Gene, DNV_cases$Gene))
ctrl_LoF_genes = unique(c(LoF_ctrls$Gene, DNV_cases$Gene))

iRef_cases = unique(PPI[which(PPI$aliasA %in% case_LoF_genes & PPI$aliasB %in% case_LoF_genes),])
iRef_ctrls = unique(PPI[which(PPI$aliasA %in% ctrl_LoF_genes & PPI$aliasB %in% ctrl_LoF_genes),])

iRef_cases = iRef_cases[which(!duplicated(iRef_cases[,c("aliasA", "aliasB")])),]
iRef_ctrls = iRef_ctrls[which(!duplicated(iRef_ctrls[,c("aliasA", "aliasB")])),]

write.csv(iRef_cases, paste0(out_path, "/Cytoscape_LoF_network.csv"), 
          row.names = FALSE, quote = FALSE)

# # # # # # # # # # CALCULATE NODE DEGREE # # # # # # # # # #
edge_1 = as.data.frame(table(PPI$aliasA))
edge_2 = as.data.frame(table(PPI$aliasB))
names(edge_1) = c("Var1", "Freq1")
both = merge(edge_1, edge_2, all=TRUE)
both[is.na(both)] <- 0
both$degree = rowSums(both[,c("Freq1","Freq")])

# Identify nodes that cannot be permuted based on node degree
degree_freq = as.data.frame(table(both$degree))
singletons = as.vector(degree_freq$Var1[which(degree_freq$Freq == 1)])
uni_nodes = as.vector(both$Var1[which(both$degree %in% singletons)])

# # # # # # # # # # NETWORK PERMUTATION # # # # # # # # # #
net_perm = 1000

# Save disease-degree of each node
degree_frame = data.frame(matrix(NA, nrow = net_perm, ncol = length(case_LoF_genes)))
names(degree_frame) = case_LoF_genes

# Save overall density of each permuted network
perm_dens = c()

# These degree values can be swapped around
swappable_degree = as.vector(degree_freq$Var1[which(degree_freq$Freq > 1)])

for(i in c(1:net_perm)){
  print(i)
  both$new = both$Var1
  node_lookup = both[,c("Var1", "new")]
  names(node_lookup) = c("old","new")
  
  # Node switching step
  for (k in swappable_degree){
    Ki = as.vector(both$Var1[which(both$degree == k)])
    swap = as.vector(sample(Ki))
    node_lookup$new[which(node_lookup$old %in% Ki)] = 
      swap[match(node_lookup$old[which(node_lookup$old %in% Ki)], Ki)]
  }
  
  G0 = PPI[,c("aliasA", "aliasB", "interactionType")]
  G0 = G0[which(!duplicated(G0[,c("aliasA", "aliasB")])),]
  G0$newA = node_lookup$new[match(G0$aliasA, node_lookup$old)]
  G0$newB = node_lookup$new[match(G0$aliasB, node_lookup$old)]
  G0 = G0[,c("newA","newB")]
  
  # Edge switching step
  G_unique = G0[which(G0$newA %in% uni_nodes | G0$newB %in% uni_nodes),]
  edge_perm = 2
  for(g in c(1:edge_perm)){
    seen = c()
    while(length(seen) < length(uni_nodes)){
      swap = as.data.frame(G_unique[sample(nrow(G_unique), 2),])
      AD = setNames(data.frame(as.character(swap[1,"newA"]), as.character(swap[2,"newB"])), 
                    c("newA", "newB"))
      BC = setNames(data.frame(as.character(swap[1,"newB"]), as.character(swap[2,"newA"])),
                    c("newA","newB"))
      
      if(nrow(match_df(G0,AD)) + nrow(match_df(G0,BC)) < 1){
        # Remove originals
        G0 <- G0[!(G0$newA == swap[1,1] & G0$newB == swap[1,2]),]
        G0 <- G0[!(G0$newA == swap[2,1] & G0$newB == swap[2,2]),]
        # Replace with new
        G0 <- rbind(G0, AD)
        G0 <- rbind(G0, BC)
      }
      seen = unique(c(seen, as.character(swap[1,"newA"]), 
                      as.character(swap[1,"newB"]), as.character(swap[2,"newA"]), 
                      as.character(swap[2,"newB"])))
    }
  }
  
  # Get disease subnetwork density
  disease_net = G0[which(G0$newA %in% case_LoF_genes | G0$newB %in% case_LoF_genes),]
  direct_disease_net = G0[which(G0$newA %in% case_LoF_genes & G0$newB %in% case_LoF_genes),]
  
  perm_dens = c(perm_dens, nrow(direct_disease_net))
  
  # Get number of connections to disease proteins
  for (j in c(1:ncol(degree_frame))){
    degree_frame[i,j] = nrow(direct_disease_net[which(
      direct_disease_net$newA == names(degree_frame)[j] | 
        direct_disease_net$newB == names(degree_frame)[j]),])
  }
}

# Save out permutation data to use later
today <- format(Sys.Date(), format = "%Y%b%d")
write.csv(degree_frame[c(1:7),], 
          file = paste0(base_path, "/intermediate/network_permutations/LoFDNV_", today, ".csv"),
          quote = FALSE, row.names = FALSE)

# # # # # # # # # # CALCULATE P-VALUES # # # # # # # # # #
true_net_dens = nrow(PPI[which(PPI$aliasA %in% case_LoF_genes & PPI$aliasB %in% case_LoF_genes),])
gene_dens <- vector(mode="integer", length=ncol(degree_frame))
for (j in c(1:ncol(degree_frame))){
  gene_PPI = PPI[which(PPI$aliasA == names(degree_frame[j]) |
                         PPI$aliasB == names(degree_frame[j])),]
  gene_dens[j] = nrow(gene_PPI[which(gene_PPI$aliasA %in% case_LoF_genes &
                                       gene_PPI$aliasB %in% case_LoF_genes),])
}

# Network permutation p-value
nGreater = length(perm_dens[which(perm_dens >= true_net_dens)])
pval = (nGreater/length(perm_dens))
write.table(paste0("Network connectivity p-value: ", pval), file = paste0(out_path, "/network_pval_LoFDNV.txt"))

# Node p-values & node file
nodes = data.frame(matrix(NA, nrow = length(case_LoF_genes)), stringsAsFactors = FALSE)
nodes$Gene = names(degree_frame)
nodes$connections = gene_dens
nodes$nGreater = NA
row.names(nodes) <- nodes$Gene
for(g in nodes$Gene){
  nodes[g,"nGreater"] = length(degree_frame[,g][which(degree_frame[,g] >= nodes[g,"connections"])])
}
nodes$pval = nodes$nGreater/nrow(degree_frame)

# BH correction
nodes$BH = p.adjust(nodes$pval, method = "BH")
nodes_adj = nodes[,c("Gene", "BH")]

write.csv(nodes_adj, file = paste0(out_path, "/Cytoscape_LoFDNV_nodes.csv"), 
          row.names = FALSE, quote = FALSE)

# # # # # # SCRATCH # # # # # #

# Read in permutation results
#nodes_adj = read.csv(file = paste0(out_path, "/Cytoscape_LoFDNV_nodes.csv"), stringsAsFactors = FALSE)
#rank1 = ((nrow(nodes_adj[which(nodes_adj$BH == 0),])+1)/2)
#nodes_adj$BH[which(nodes_adj$BH == 0)] <- (rank1/nrow(nodes_adj)) * 0.05
write.csv(nodes_adj, file = paste0(out_path, "/Cytoscape_LoFDNV_nodes.csv"), 
          row.names = FALSE, quote = FALSE)
