library(artMS)


####################################################
# create_artMSfiles: writes out key file
####################################################
create_artMSfiles <- function(int_type, key_file){
  
  # Create artMS input files
  #keys = read.table(paste0(base_path, "/input/evidence/expTFs/annotation_", int_type, ".txt"), 
  #                  na.strings=c("", "NA"), sep = "\t", header = TRUE)
  keys = read.table(paste0(base_path, "/input/evidence/interdependence/annotation_", int_type, ".txt"), 
                    na.strings=c("", "NA"), sep = "\t", header = TRUE)
  keys = keys[which(!is.na(keys$Condition)),]
  keys = keys[,c('Raw.file', 'IsotopeLabelType', 'Condition', 'BioReplicaSaint', 'Run', "SAINT")]
  names(keys) = c("RawFile", "IsotopeLabelType", "Condition", "BioReplicate", "Run", "SAINT")
  keys$BioReplicate = gsub("_","-",keys$BioReplicate)
  
  write.table(keys, file =  key_file,
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  #contrast = write.table(paste0("control-",int_type), 
  #                       file = paste0(base_path, "/input/evidence/expTFs/contrast_", int_type, ".txt"),
  #                       quote = FALSE, row.names = FALSE, col.names = FALSE)
}

####################################################
# run_artMS: creates saintq and saintExpress inputs
####################################################
run_artMS <- function(int_type, ev_file, key_file, quant_var){
  
  if(int_type != "NKX25"){
    artmsQuantification(yaml_config_file = paste0(base_path, "/input/evidence/interdependence/", int_type, "_config.yaml"))
  }
  
  artmsEvidenceToSaintExpress(evidence_file = ev_file, 
                              keys_file = key_file, 
                              ref_proteome_file = paste0(base_path, "/input/databases/UP000005640_9606.fasta"),
                              quant_variable = quant_var,
                              output_file = paste0(base_path, "/input/evidence/interdependence/saintx_inputs/",int_type, "_", quant_var, ".txt"))
  
  artmsEvidenceToSAINTq(evidence_file = ev_file, 
                        keys_file = key_file, 
                        output_dir = paste0(base_path, "/input/evidence/interdependence/saintq_inputs/", int_type))
}