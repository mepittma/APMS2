# This script saves out files to observe the differences between the original saintq interaction scoring and the
# highest-performing scores discovered in interactome_selection.R

# These files were created by running the saintq algorithm at the protein level for the GATA4 interactome
# and the peptide level for TBX5 and NKX2-5. Normalization was enabled. A BFDR cutoff of 0.05 was selected.
# A blacklist of genes with significant negative differential expression between WT and KO cell lines was filtered out.
# This was determined by genes with negative fold change and an FDR cutoff of 0.05 in edgeR differential expression.

base_path = "/Users/student/Documents/PollardLab/APMS2"
unimap = read.table(paste0(base_path, "/input/aliases/combined_unimap.txt"),header = TRUE, sep = "\t", stringsAsFactors = FALSE)
unimap <- unique(unimap)

# Load in the normalized data

GATA_tab = read.table(file = paste0(base_path, "/intermediate/interaction_scoring/saintq/GATA4_protein_norm_true.txt"), 
                      sep = "\t", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, comment.char = "")
TBX5_tab = read.table(file = paste0(base_path,"/intermediate/interaction_scoring/saintq/TBX5_peptide_norm_true.txt"),
                      sep = "\t", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, comment.char = "")
NKX25_tab = read.table(file = paste0(base_path,"/intermediate/interaction_scoring/saintq/NKX25_peptide_norm_true.txt"),
                       sep = "\t", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, comment.char = "")
comb = rbind(GATA_tab, NKX25_tab[,c("Bait","Prey","X.Rep","AvgP","BFDR")])
norm_ints = rbind(comb, TBX5_tab[,c("Bait","Prey","X.Rep","AvgP","BFDR")])
norm_ints = norm_ints[which(norm_ints$BFDR <= 0.05),]
norm_ints$Prey_geneName = unimap$GeneSymbol[match(norm_ints$Prey, unimap$UniProt)]

write.table(norm_ints, 
            file = paste0(base_path, "/input/precomp_interactomes/saintq_n_interactomes.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Load in the originally-used data
int_file = paste0(base_path, "/input/precomp_interactomes/interactomes.csv")
ints <- read.table(int_file, 
                   sep = "\t",fill = TRUE, header = TRUE, stringsAsFactors = FALSE)

# Compare members of norm_ints and ints
GATA4_norm_prots = norm_ints$Prey[which(norm_ints$Bait == "GATA4")]
TBX5_norm_prots = norm_ints$Prey[which(norm_ints$Bait == "TBX5")]
NKX25_norm_prots = norm_ints$Prey[which(norm_ints$Bait == "NKX25")]

GATA4_f_prots = ints$Prey_proteinname[which(ints$Bait == "GATA4")]
TBX5_f_prots = ints$Prey_proteinname[which(ints$Bait == "TBX5")]
NKX25_f_prots = ints$Prey_proteinname[which(ints$Bait == "NKX25")]

GATA4_f_prots[which(!GATA4_f_prots %in% GATA4_norm_prots)] # 2 of the original didn't make it
TBX5_f_prots[which(!TBX5_f_prots %in% TBX5_norm_prots)] # 18 of the original didn't make it
NKX25_f_prots[which(!NKX25_f_prots %in% NKX25_norm_prots)] # 24 of the original didn't make it

GATA4_norm_prots[which(!GATA4_norm_prots %in% GATA4_f_prots)] # 302 new interactions
TBX5_norm_prots[which(!TBX5_norm_prots %in% TBX5_f_prots)] # 58 new interactions
NKX25_norm_prots[which(!NKX25_norm_prots %in% NKX25_f_prots)] # 1 new interaction

# Does blacklisting do anything?
GATA4_norm_genes = norm_ints$Prey_geneName[which(norm_ints$Bait == "GATA4")]
TBX5_norm_genes = norm_ints$Prey_geneName[which(norm_ints$Bait == "TBX5")]
NKX25_norm_genes = norm_ints$Prey_geneName[which(norm_ints$Bait == "NKX25")]

GATA4_f_genes = ints$Prey_genename[which(ints$Bait == "GATA4")]
TBX5_f_genes = ints$Prey_genename[which(ints$Bait == "TBX5")]
NKX25_f_genes = ints$Prey_genename[which(ints$Bait == "NKX25")]

TBX5_blacklist = read.table(paste0(base_path, "/intermediate/rnaseq/TBX5_negativeFC_blacklist.txt"))
length(TBX5_norm_genes[which (TBX5_norm_genes %in% TBX5_blacklist)])
length(TBX5_f_genes[which (TBX5_f_genes %in% TBX5_blacklist)])

NKX25_blacklist = read.table(paste0(base_path, "/intermediate/rnaseq/NKX25_negativeFC_blacklist.txt"))
length(NKX25_norm_genes[which (NKX25_norm_genes %in% NKX25_blacklist)])
length(NKX25_f_genes[which (NKX25_f_genes %in% NKX25_blacklist)])

GATA4_blacklist = read.table(paste0(base_path, "/intermediate/rnaseq/GATA4_negativeFC_blacklist.txt"))
length(GATA4_norm_genes[which (GATA4_norm_genes %in% GATA4_blacklist)])
length(GATA4_f_genes[which (GATA4_f_genes %in% GATA4_blacklist)])
