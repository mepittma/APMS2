base_path = "/Users/student/Documents/PollardLab/APMS2"

####################################################
# Read in the data
####################################################

# Read in alias conversion file
unimap = read.table(paste0(base_path, "/input/aliases/UniMap.txt"),header = TRUE, sep = "\t", stringsAsFactors = FALSE)
quant_unimap = read.table(paste0(base_path, "/input/aliases/combined_unimap.txt"),header = TRUE, sep = "\t", stringsAsFactors = FALSE)
quant_unimap <- unique(quant_unimap)

####################################################
# Read in the functions
####################################################
source(paste0(base_path, "/code/functions/evidence2inputs.R"))
#create_artMSfiles, run_artMS

source(paste0(base_path, "/code/functions/score_interactions.R"))
#run_FCA, run_FCB, run_saintx, run_saintq, write_params

source(paste0(base_path, "/code/functions/score_interactomes.R"))
#get_comps, get_prots, get_fisher, iref_response_list

source(paste0(base_path, "/code/functions/permutation_functions.R"))
#get_fisher, get_OR, get_prots, get_sig_compl, perm_viz, permute_status

####################################################
# Break apart evidence file
####################################################

ev = read.table(paste0(base_path, "/input/evidence/expTFs/bg15_evidence.txt"),
                sep = "\t", header = T, stringsAsFactors = FALSE)
ctrl = ev[which(grepl("gko",ev$Experiment)),]
old = ev[which(grepl("old",ev$Experiment)),]
wt = ev[which(grepl("wt_g",ev$Experiment)),]
g296 = ev[which(grepl("296",ev$Experiment)),]

write.table(rbind(ctrl,old), 
            file = paste0(base_path, "/input/evidence/expTFs/old_gata4_evidence.txt"),
            sep = "\t", row.names = F, quote = F)
write.table(rbind(ctrl,wt), 
            file = paste0(base_path, "/input/evidence/expTFs/wt_gata4_evidence.txt"),
            sep = "\t", row.names = F, quote = F)
write.table(rbind(ctrl,g296), 
            file = paste0(base_path, "/input/evidence/expTFs/g296s_gata4_evidence.txt"),
            sep = "\t", row.names = F, quote = F)

####################################################
# Run quantification
####################################################
# APMS data for the three interactomes
for (int_type in c("old_gata4", "wt_gata4", "g296s_gata4")){
  
  ev_file = paste0(base_path, "/input/evidence/expTFs/", int_type, "_evidence.txt")
  key_file= paste0(base_path, "/input/evidence/expTFs/keys_", int_type,".txt")
  
  # Creates files necessary to run FC, saintq, and saint analysis - uncomment if these files do not exist
  #create_artMSfiles(int_type, key_file)
  #setwd(paste0(base_path, "/intermediate/MS_QC/", int_type, "_msspc"))
  #run_artMS(int_type, ev_file, key_file, "msspc")
  #setwd(paste0(base_path, "/intermediate/MS_QC/", int_type, "_msint"))
  #run_artMS(int_type, ev_file, key_file, "msint")
  
  # Run saintq with all possible parameter files 
  for(level in c("peptide", "protein")){
    
    setwd(paste0(base_path, "/input/evidence/interdependence/saintq_inputs/", int_type))
    in_file = paste0("saintq_input_", level, "s.txt")
    
    normYN = "true"
    param_out = paste0( level, "_norm", normYN, ".txt")
    write_params(normYN, in_file, level, param_out)
    out_file = paste0(base_path, "/intermediate/interaction_scoring/saintq/", int_type, "_", level, "_norm_", normYN, ".txt")
    run_saintq(paste0(base_path, "/input/evidence/interdependence/saintq_inputs/", int_type), param_out, in_file, out_file)
    
  }
  
}

####################################################
# Interdependence analysis 
####################################################
# Which proteins are different across runs?
# Using protein level, normalization
level = "protein"
old_gata4 = read.table(paste0(base_path, 
                              "/intermediate/interaction_scoring/saintq/", 
                              "old_gata4_", level, "_norm_true.txt"),
                       sep = "\t", stringsAsFactors=F, comment.char="", header=T)
old_gata4_interactors = old_gata4$Prey[which(old_gata4$BFDR < 0.001)]

# Confirm that these are the same as the old experiments
confirm = read.table(file = paste0(base_path, "/intermediate/interaction_scoring/saintq/GATA4_protein_norm_true.txt"), 
                 sep = "\t", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, comment.char = "")
confirm_interactors = confirm$Prey[which(confirm$BFDR < 0.001)]
different = c(confirm_interactors[!confirm_interactors %in% old_gata4_interactors],
              old_gata4_interactors[!old_gata4_interactors %in% confirm_interactors])

######## ADDRESS: why are there 18 proteins with different BFDR cutoffs? 
# Answer: saintq has different Bayesian FDRs depending on which proteins are looked at...
# But what's the difference between the two source files:
#   /input/evidence/interdependence/saintq_inputs/old_gata4/saintq_input_peptides.txt
#   /input/evidence/inout/evidence/expTFs/saintq_inputs/GATA4/saintq_input_peptides.txt
# It appears there are 8 additional lines in old_gata4 file... (e.g. O75352)

# Check what got lost/gained in the new WT analysis
wt_gata4 = read.table(paste0(base_path, 
                             "/intermediate/interaction_scoring/saintq/", 
                             "wt_gata4_", level, "_norm_true.txt"),
                      sep = "\t", stringsAsFactors=F, comment.char="", header=T)
wt_gata4_interactors = wt_gata4$Prey[which(wt_gata4$BFDR < 0.001)]
lost = old_gata4_interactors[!old_gata4_interactors %in% wt_gata4_interactors] #208 proteins!
new = wt_gata4_interactors[!wt_gata4_interactors %in% old_gata4_interactors] #40 proteins!
lost = confirm_interactors[!confirm_interactors %in% wt_gata4_interactors] #204 proteins!
new = wt_gata4_interactors[!wt_gata4_interactors %in% confirm_interactors] #40 proteins!
common = intersect(wt_gata4_interactors, old_gata4_interactors) # only 30.......

# # # # # # # # UNDER CONSTRUCTION # # # # # # # # 
# Read in results file and compare to non-KO interactome
if(grepl("GATA4", int_type)){
  WT_str = "GATA4"
  level = "protein"
  comp_file = paste0(base_path, "/intermediate/interaction_scoring/saintq/GATA4_", level, "_norm_true.txt")
} else if(grepl("NKX25", int_type)){
  WT_str = "NKX25"
  level = "peptide"
  comp_file = paste0(base_path, "/intermediate/interaction_scoring/saintq/NKX25_", level, "_norm_true.txt")
}

KO = read.table(paste0(base_path, "/intermediate/interaction_scoring/saintq/", int_type, "_", level, "_norm_true.txt"),
                sep = "\t", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, comment.char = "")
WT = read.table(comp_file, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, comment.char = "")

sig_KO = KO[which(KO$BFDR <= 0.05),]
sig_WT = WT[which(WT$BFDR <= 0.05),]

sig_KO$Prey_genename = unimap$GeneSymbol[match(sig_KO$Prey, unimap$UniProt)]
sig_WT$Prey_genename = unimap$GeneSymbol[match(sig_WT$Prey, unimap$UniProt)]

common_genes = sig_KO$Prey_genename[which(sig_KO$Prey_genename %in% sig_WT$Prey_genename)]

write.csv(sig_KO, file = paste0(base_path, "/intermediate/interdependence/", int_type, "_significant_interactions.csv"),
          row.names = FALSE, quote = FALSE)
write.csv(sig_WT, file = paste0(base_path, "/intermediate/interdependence/", WT_str, "_significant_interactions.csv"),
          row.names = FALSE, quote = FALSE)
write.table(common_genes, file = paste0(base_path, "/intermediate/interdependence/", WT_str, "_", int_type, "_common_genes.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
