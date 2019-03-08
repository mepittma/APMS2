base_path = "/Users/student/Documents/PollardLab/APMS2"

####################################################
# Read in the data
####################################################

# Load in interactome data.
int_file = paste0(base_path, "/input/precomp_interactomes/interactomes.csv")
ints <- read.table(int_file, 
                   sep = "\t",fill = TRUE, header = TRUE, stringsAsFactors = FALSE)

# Read in corum data
corum <- read.table(paste0(base_path, "/input/databases/corum_human_CytoscapeFormatted.txt"), 
                    sep = "\t", fill = TRUE, header = TRUE, stringsAsFactors = FALSE)
corum <- as.data.frame(corum[,c(1:8)])

# Read in variant data
source(paste0(base_path, "/code/load_data/load_variants.R"))

# Read in alias conversion file
unimap <- read.table(paste0(base_path, "/input/aliases/UniMap.txt"),header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#######################################
# Save out expanded corum interactomes
######################################

for(int_type in c("TBX5", "NKX25", "GATA4")){
  
  # Select complexes with interactome members 
  int_prots = ints$Prey_proteinname[which(ints$Bait == int_type)]
  complexes = unique(c(as.vector(corum$Complex.id[which(
    corum$Interactor1 %in% int_prots| corum$Interactor2 %in% int_prots)])))
  int_comp <- corum[which(corum$Complex.id %in% complexes),]
  
  # Save out names of complexes
  complex_names = unique(corum[which(corum$Complex.id %in% complexes), c("Complex.id", "Complex.name")])
  write.table(complex_names, 
              file = paste0(base_path, "/output/enrichment_analysis/expandedInteractomes/", int_type, "_complexNames.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  # Save out complex members
  int_comp$geneName1 <- unimap$GeneSymbol[match(int_comp$Interactor1, as.vector(unimap$UniProt))]
  int_comp$geneName2 <- unimap$GeneSymbol[match(int_comp$Interactor2, as.vector(unimap$UniProt))]
  write.table(int_comp[,c("Interactor1", "Interactor2", "geneName1", "geneName2", "Complex.name", "Complex.id")],
              file = paste0(base_path, "/output/enrichment_analysis/expandedInteractomes/cytoscapeFormat_", int_type, ".txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE )
  
  complex_members = unique(c(int_comp$geneName1, int_comp$geneName2))
  complex_members = unlist(strsplit(complex_members, "; "))
  write.table(complex_members, 
              file = paste0(base_path, "/output/enrichment_analysis/expandedInteractomes/", int_type, "_complexMembers.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
}

#######################################
# Save out complexes in common
######################################

# Find any prey proteins in common
common_prey = ints$Prey_proteinname[which(duplicated(ints$Prey_proteinname))]
common_table = ints[which(ints$Prey_proteinname %in% common_prey),]
write.table(common_table, 
            file = paste0(base_path, "/output/enrichment_analysis/commonComplexes/original_interactome_overlaps.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Find any complexes in common
GATA_compl = unique(corum$Complex.id[which(corum$Interactor1 %in% ints$Prey_proteinname[which(ints$Bait == "GATA4")]  |
                                            corum$Interactor2 %in% ints$Prey_proteinname[which(ints$Bait == "GATA4")])])
TBX5_compl = unique(corum$Complex.id[which(corum$Interactor1 %in% ints$Prey_proteinname[which(ints$Bait == "TBX5")]  |
                                            corum$Interactor2 %in% ints$Prey_proteinname[which(ints$Bait == "TBX5")])])
NKX25_compl = unique(corum$Complex.id[which(corum$Interactor1 %in% ints$Prey_proteinname[which(ints$Bait == "NKX25")]  |
                                             corum$Interactor2 %in% ints$Prey_proteinname[which(ints$Bait == "NKX25")])])

combined = c(GATA_compl, TBX5_compl, NKX25_compl)
common = combined[which(duplicated(combined))]
common_frame = unique(corum[which(corum$Complex.id %in% common), c("Interactor1", "Interactor2","Complex.name", "Complex.id")])
common_frame$geneName1 <- unimap$GeneSymbol[match(common_frame$Interactor1, as.vector(unimap$UniProt))]
common_frame$geneName2 <- unimap$GeneSymbol[match(common_frame$Interactor2, as.vector(unimap$UniProt))]

common_frame$GATA4 = 0
common_frame$TBX5 = 0
common_frame$NKX25 = 0

for (i in c(1:nrow(common_frame))){
  if (common_frame[i,4] %in% GATA_compl){common_frame[i,7] = 1} 
  if (common_frame[i,4] %in% TBX5_compl){common_frame[i,8] = 1} 
  if (common_frame[i,4] %in% NKX25_compl){common_frame[i,9] = 1} 
}

write.table(common_frame,
            file = paste0(base_path, "/output/enrichment_analysis/commonComplexes/commonComplexes_cytoscapeFormat.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE )
write.table(unique(common_frame[,c("Complex.name", "Complex.id")]),
            file = paste0(base_path, "/output/enrichment_analysis/commonComplexes/commonComplexes_list.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE )

# Find any complex members in common between the three
GATA_prots = unique(c(corum$Interactor1[which(corum$Complex.id %in% GATA_compl)], 
                      corum$Interactor2[which(corum$Complex.id %in% GATA_compl)]))
TBX5_prots = unique(c(corum$Interactor1[which(corum$Complex.id %in% TBX5_compl)], 
                      corum$Interactor2[which(corum$Complex.id %in% TBX5_compl)]))
NKX25_prots = unique(c(corum$Interactor1[which(corum$Complex.id %in% NKX25_compl)], 
                      corum$Interactor2[which(corum$Complex.id %in% NKX25_compl)]))

combined = c(GATA_prots, TBX5_prots, NKX25_prots)
common = combined[which(duplicated(combined))]
protein_connectors = as.data.frame(cbind(common, unimap$GeneSymbol[match(common, as.vector(unimap$UniProt))]))

protein_connectors$GATA4_expanded = 0
protein_connectors$TBX5_expanded = 0
protein_connectors$NKX25_expanded = 0

for (i in c(1:nrow(protein_connectors))){
  if (protein_connectors[i,1] %in% GATA_prots){protein_connectors[i,3] = 1} 
  if (protein_connectors[i,1] %in% TBX5_prots){protein_connectors[i,4] = 1} 
  if (protein_connectors[i,1] %in% NKX25_prots){protein_connectors[i,5] = 1} 
}

names(protein_connectors)[c(1,2)] = c("UniProt", "GeneSymbol")

write.table(protein_connectors,
            file = paste0(base_path, "/output/enrichment_analysis/commonComplexes/expanded_interactome_overlaps.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE )

#######################################################
# Save out corum complexes enriched in each interactome
#######################################################

# Overrepresentation test - Given a list of genes, analyzes whether the supplied list is significantly associated 
# with a particular complex (instead of randomly scattered throughout the whole set of possible complexes).
# Background set: all interactors in the CORUM file

# Function to isolate proteins given a complex ID
get_prots <- function(compl_id, corum){
  return(unique(c(corum$Interactor1[which(corum$Complex.id == compl_id)],
                  corum$Interactor2[which(corum$Complex.id == compl_id)])))
}

# Function to perform fisher exact test given input list and complex id
get_fisher <- function(gList, compl_id, corum, total_prots){
  compl_prots = get_prots(compl_id, corum)
  
  int_comp = length(gList[gList %in% compl_prots])
  int_non_comp = length(gList[!(gList %in% compl_prots)])
  non_int_comp = length(compl_prots[!(compl_prots %in% gList)])
  non_int_non_comp = length(total_prots) - (int_comp + int_non_comp + non_int_comp)
  
  fisherMat <- matrix(c(int_comp, int_non_comp, non_int_comp, non_int_non_comp), nrow = 2,
                      dimnames = list(Complex = c("Yes","No"),
                                      Interactome = c("Yes", "No")))
  
  # Fisher test
  fish = fisher.test(fisherMat, alternative = 'greater')
  return(fish[1])
  
}

# Get list of all background proteins
total_prots = unique(c(corum$Interactor1, corum$Interactor2))
compl_list = unique(corum$Complex.id)

for(int_type in c("TBX5", "NKX25", "GATA4")){
  
  int_prots = ints$Prey_proteinname[which(ints$Bait == int_type)]

  sig_complexes = c()
  for (compl in compl_list){
    fish = get_fisher(int_prots, compl, corum, total_prots)
    
    if (fish <= 0.05){
      sig_complexes = c(sig_complexes, compl)
    }
  }
  
  sig_frame = unique(corum[which(corum$Complex.id %in% sig_complexes),c("Complex.id", "Complex.name", "subunits..UniProt.IDs.")])
  write.table(sig_frame, 
              file = paste0(base_path, "/output/enrichment_analysis/enriched4interactors/", int_type, "_complexes.txt"),
              quote = FALSE, row.names = FALSE, col.names = TRUE)
  
}

#######################################################
# Save out complexes enriched for PCGC variants
#######################################################

# Function to calculate the p-value of a group of genes
get_poisson <- function(case_table, ctrl_table, gene_list){
  
  case_ct = 0
  ctrl_ct = 0
  case_non_ct = 0
  ctrl_non_ct = 0
  
  # Count the number of times a case individual had a damaging mutation in an interactome gene
  # and the number of times a case individual had a damaging mutation in a non-interactome gene
  for(i in 1:nrow(case_table)){
    if(case_table[i,'Gene'] %in% gene_list){
      case_ct = case_ct + 1
    } else{
      case_non_ct = case_non_ct + 1
    }
  }
  
  # Count the number of times a control individual had a damaging mutation in an interactome gene
  # and the number of times a control individual had a damaging mutation in a non-interactome gene
  for(i in 1:nrow(ctrl_table)){
    if(ctrl_table[i,'Gene'] %in% gene_list){
      ctrl_ct = ctrl_ct + 1
    } else{
      ctrl_non_ct = ctrl_non_ct + 1
    }
  }
  
  # Add pseudocount if any of the values are equal to zero
  if (0 %in% c(case_ct, ctrl_non_ct, ctrl_ct, case_non_ct)){
    case_ct = case_ct + 0.5
    ctrl_non_ct = ctrl_non_ct + 0.5
    ctrl_ct = ctrl_ct + 0.5
    case_non_ct = case_non_ct + 0.5
  }
  
  # Create data table
  my_data <- data.frame(Status = c(0,1,0,1),
                        Inter = c(1,1,0,0),
                        Freq=c(ctrl_ct,case_ct,ctrl_non_ct,case_non_ct),
                        stringsAsFactors=FALSE) 
  
  odds = (case_ct * ctrl_non_ct) / (ctrl_ct * case_non_ct)
  fit <- glm(Status ~ Inter, weights = Freq, data = my_data,family = poisson())
  return(coef(summary(fit))[2,4])
}


# Find corum complexes with members in the interactomes.
int_prots <- unique(c(as.vector(ints$Prey_proteinname), "P43694", "Q99593", "P52952"))
int_genes <- unique(c(as.vector(ints$Prey_genename), "GATA4", "TBX5", "NKX25"))
int_complexes <- unique(c(as.vector(corum$Complex.id[which(
  corum$Interactor1 %in% int_prots| corum$Interactor2 %in% int_prots)])))
int_comp <- corum[which(corum$Complex.id %in% int_complexes),]

# Get gene names of corum complex members
int_comp$geneName1 <- unimap$GeneSymbol[match(int_comp$Interactor1, as.vector(unimap$UniProt))]
int_comp$geneName2 <- unimap$GeneSymbol[match(int_comp$Interactor2, as.vector(unimap$UniProt))]

# Create a lookup table of corum complex <-> p-value
complex_lookup <- data.frame(matrix(NA, nrow = length(int_complexes), ncol = 10), stringsAsFactors = FALSE)
names(complex_lookup) <- c("ComplexID", "ComplexName", "DNV_pval", "LoF_pval", "rec_pval", 
                           "combined_pval", "GATA4", "TBX5", "NKX25", "Genes")
complex_lookup[,1:2] <- unique(int_comp[,3:4])
for(i in c(1:length(int_complexes))){
  gene_list <- unique(c(as.vector(int_comp$geneName1[which(int_comp$Complex.id == complex_lookup$ComplexID[i])]), 
                        as.vector(int_comp$geneName2[which(int_comp$Complex.id == complex_lookup$ComplexID[i])])))
  gene_list <- unlist(strsplit(gene_list, split=";"))
  
  # Get poisson test p-value for each mutation type
  complex_lookup[i, 'DNV_pval'] <- get_poisson(DNV_cases, DNV_ctrls, gene_list)
  complex_lookup[i, 'LoF_pval'] <- get_poisson(LoF_cases, LoF_ctrls, gene_list)
  complex_lookup[i, 'rec_pval'] <- get_poisson(rec_cases, rec_ctrls, gene_list)
  complex_lookup[i, 'combined_pval'] <- get_poisson(combined_cases, combined_ctrls, gene_list)
  
  # Mark membership
  complex_lookup[i, 7:9] <- 0
  if (any(gene_list %in% ints$Prey_genename[which(ints$Bait == "GATA4")])){
    complex_lookup[i, 'GATA4'] <- 1
  }
  if (any(gene_list %in% ints$Prey_genename[which(ints$Bait == "TBX5")])){
    complex_lookup[i, 'TBX5'] <- 1
  }
  if (any(gene_list %in% ints$Prey_genename[which(ints$Bait == "NKX25")])){
    complex_lookup[i, 'NKX25'] <- 1
  }
  
  # Append gene list
  complex_lookup[i, 'Genes'] <- paste(gene_list, collapse = ";")
  
}

# Order by combined pvalue and save out
write.csv(complex_lookup[order(complex_lookup$combined_pval),], 
          file = paste0(base_path, "/output/enrichment_analysis/enriched4PCGC/interactome_complex_enrichment.csv"), 
          row.names = FALSE) 


