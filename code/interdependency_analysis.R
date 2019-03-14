base_path = "/Users/student/Documents/PollardLab/APMS2"

source(paste0(base_path, "/code/functions/evidence2inputs.R"))
#create_artMSfiles, run_artMS

source(paste0(base_path, "/code/functions/score_interactions.R"))
#run_FCA, run_FCB, run_saintx, run_saintq, write_params

unimap = read.table(paste0(base_path, "/input/aliases/combined_unimap.txt"),header = TRUE, sep = "\t", stringsAsFactors = FALSE)
unimap <- unique(unimap)

for (int_type in c("NKX25_GKO", "GATA4_TKO_NKO")){
  
  ev_file = paste0(base_path, "/input/evidence/interdependence/", int_type, "_evidence.txt")
  key_file= paste0(base_path, "/input/evidence/interdependence/keys_", int_type, ".txt")
  
  # Creates files necessary to run FC, saintq, and saint analysis - uncomment if these files do not exist
  create_artMSfiles(int_type, key_file)
  setwd(paste0(base_path, "/intermediate/MS_QC/", int_type, "_msspc"))
  run_artMS(int_type, ev_file, key_file, "msspc")
  setwd(paste0(base_path, "/intermediate/MS_QC/", int_type, "_msint"))
  run_artMS(int_type, ev_file, key_file, "msint")
  
  # Run saintq with desired parameters
  for(level in c("peptide", "protein")){
    
    setwd(paste0(base_path, "/input/evidence/interdependence/saintq_inputs/", int_type))
    in_file = paste0("saintq_input_", level, "s.txt")
    
    normYN = "true"
    param_out = paste0( level, "_norm", normYN, ".txt")
    write_params(normYN, in_file, level, param_out)
    out_file = paste0(base_path, "/intermediate/interaction_scoring/saintq/", int_type, "_", level, "_norm_", normYN, ".txt")
    run_saintq(paste0(base_path, "/input/evidence/interdependence/saintq_inputs/", int_type),param_out, in_file, out_file)
    
  }
  
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
  
}
