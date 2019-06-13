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
# Break apart/examine evidence file
####################################################

ev = read.table(paste0(base_path, "/input/evidence/expTFs/bg15_evidence.txt"),
                sep = "\t", header = T, stringsAsFactors = FALSE)

# PCA clustering of experiments
library(ggfortify)
ev_subset <- ev[,c("Proteins","MS.MS.count","Experiment")]
ev_subset$GeneName <- quant_unimap$GeneSymbol[match(ev_subset$Proteins,quant_unimap$UniProt)]
#Aggregate by protein
agg_sub = data.frame(matrix(NA, nrow=length(unique(ev_subset$Proteins)), ncol = length(unique(ev_subset$Experiment))))
names(agg_sub) <- unique(ev_subset$Experiment)
row.names(agg_sub) <- unique(ev_subset$Proteins)
for (i in c(1:nrow(ev_subset))){
  agg_sub[ev_subset$Proteins[i],ev_subset$Experiment[i]] <- ev_subset$MS.MS.count[i]
}
agg_sub[is.na(agg_sub)] <- 0
agg <- t(agg_sub)
exp_types = c("original_wt", "original_wt", "gko", "original_wt", "gko", "gko",
              "new_wt", "new_wt", "new_wt", "g296s", "g296s", "g296s", "gko")
agg_data <- as.data.frame(agg)
agg_data$exp <- as.vector(exp_types)
PC <- prcomp(agg)
PCi <- data.frame(PC$x, Experiment = exp_types)
ggplot(PCi, aes(x=PC1, y = PC2, col=Experiment))+
  geom_point(size=5, alpha=0.5)+
  theme_classic()

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
level = "protein"
wt_gata4 = read.table(paste0(base_path, 
                             "/intermediate/interaction_scoring/saintq/", 
                             "wt_gata4_", level, "_norm_true.txt"),
                      sep = "\t", stringsAsFactors=F, comment.char="", header=T)
wt_gata4_interactors = wt_gata4$Prey[which(wt_gata4$BFDR < 0.05)]

lost = old_gata4_interactors[!old_gata4_interactors %in% wt_gata4_interactors] #208 proteins!
new = wt_gata4_interactors[!wt_gata4_interactors %in% old_gata4_interactors] #40 proteins!
lost = confirm_interactors[!confirm_interactors %in% wt_gata4_interactors] #204 proteins!
new = wt_gata4_interactors[!wt_gata4_interactors %in% confirm_interactors] #40 proteins!
common = intersect(wt_gata4_interactors, old_gata4_interactors) # only 30...39 when threshold is relaxed to 0.05
# Should I test the permutation using the overlapping variants?

# Read in the input_proteins files to understand what is different between these experiments
new_input = read.table(paste0(base_path, 
                      "/input/evidence/interdependence/saintq_inputs/wt_gata4/saintq_input_proteins.txt"),
                       sep = "\t", stringsAsFactors=F, comment.char="", skip=2, header = T)
old_input = read.table(paste0(base_path, 
                              "/input/evidence/interdependence/saintq_inputs/old_gata4/saintq_input_proteins.txt"),
                       sep = "\t", stringsAsFactors=F, comment.char="", skip=2, header = T)
names(old_input)[1:5] <- c("Proteins","Sequence","wt.1","wt.2","wt.3")

# rbind based on Proteins
all_input <- merge(x = new_input, y = old_input, by = "Proteins", all = TRUE)
all_input$Sequence.y <- NULL
autoplot(prcomp(all_input[,c(3:ncol(all_input))]))

####################################################
# Quantify combined 
####################################################

# Write out separate files

ev = read.table(paste0(base_path, "/input/evidence/expTFs/bg15_evidence_wtcombined.txt"),
                sep = "\t", header = T, stringsAsFactors = FALSE)
ctrl = ev[which(grepl("gko",ev$Experiment)),]
wt = ev[which(grepl("wt_g",ev$Experiment)),]
g296 = ev[which(grepl("296",ev$Experiment)),]

write.table(rbind(ctrl,wt), 
            file = paste0(base_path, "/input/evidence/expTFs/wt_gata4_evidence_wtcombined.txt"),
            sep = "\t", row.names = F, quote = F)
write.table(rbind(ctrl,g296), 
            file = paste0(base_path, "/input/evidence/expTFs/g296s_gata4_evidence_wtcombined.txt"),
            sep = "\t", row.names = F, quote = F)

# APMS data for the three interactomes
int_type = "wt_gata4_combined"
  
ev_file = paste0(base_path, "/input/evidence/expTFs/", int_type, "_evidence.txt")
key_file= paste0(base_path, "/input/evidence/expTFs/keys_", int_type,".txt")

# Creates files necessary to run FC, saintq, and saint analysis - uncomment if these files do not exist
create_artMSfiles(int_type, key_file)
setwd(paste0(base_path, "/intermediate/MS_QC/", int_type, "_msspc"))
run_artMS(int_type, ev_file, key_file, "msspc")
setwd(paste0(base_path, "/intermediate/MS_QC/", int_type, "_msint"))
run_artMS(int_type, ev_file, key_file, "msint")

# Run saintq with all possible parameter files 
for(level in c("peptide", "protein")){
  
  setwd(paste0(base_path, "/input/evidence/interdependence/saintq_inputs/", int_type))
  in_file = paste0("saintq_input_", level, "s.txt")
  
  normYN = "true"
  param_out = paste0( level, "_norm", normYN, ".txt")
  write_params(normYN, in_file, level, param_out)
  out_file = paste0(base_path, "/intermediate/interaction_scoring/saintq/", int_type, "_", level, "_norm_", normYN, "_wtcombined.txt")
  run_saintq(paste0(base_path, "/input/evidence/interdependence/saintq_inputs/", int_type), param_out, in_file, out_file)
  
}


# Check what got lost/gained in the new WT analysis
combined_wt_gata4 = read.table(paste0(base_path, 
                             "/intermediate/interaction_scoring/saintq/", 
                             "wt_gata4_combined_protein_norm_true_wtcombined.txt"),
                      sep = "\t", stringsAsFactors=F, comment.char="", header=T)
wt_gata4_combined_interactors = combined_wt_gata4$Prey[which(wt_gata4$BFDR < 0.001)]

wt_gata4_combined_interactors %in% wt_gata4_interactors
wt_gata4_combined_interactors %in% old_gata4_interactors

# Check what got lost/gained in the mutant (compared to wt_gata4_new)
g296s = read.table(paste0(base_path,"/intermediate/interaction_scoring/saintq/", 
                          "g296s_gata4_protein_norm_true.txt"),
                   sep = "\t", stringsAsFactors=F, comment.char="", header=T)
g296s_interactors = g296s$Prey[which(g296s$BFDR < 0.05)]
interactome_genes = quant_unimap$GeneSymbol[match(g296s_interactors, quant_unimap$UniProt)]
interactome_genes = unlist(strsplit(interactome_genes, "; "))
write.table(interactome_genes[!is.na(interactome_genes)],
            file = paste0(base_path, "/intermediate/interactome_lists/saintq_n/g296s.txt"),
            quote = F, row.names=F, col.names = F)
lost = wt_gata4_interactors[!(wt_gata4_interactors %in% g296s_interactors)]
gained = g296s_interactors[!(g296s_interactors %in% wt_gata4_interactors)]
