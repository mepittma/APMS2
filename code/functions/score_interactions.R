########################################
# Functions for FC-A and FC-B
########################################
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

get_FCA <- function(quant_table, int_type){
  
  quant_table$FCA = "NA"
  
  # Calculate normalization factors for each experiment
  norm_factors = c(NA, NA)
  
  exp_cols = which(grepl(int_type,names(quant_table)))
  ctrl_cols = which(grepl("control",names(quant_table)))
  all_cols = c(exp_cols, ctrl_cols)
  
  for (j in all_cols){
    norm_factor = sum(as.numeric(unlist(quant_table[,j])))
    norm_factors = c(norm_factors, norm_factor)
    aveN = mean(norm_factors[c(ctrl_cols)])
  }
  
  # Calculate fold change for each protein
  for (i in c(1:nrow(quant_table))){
    ctrl_norms = c()
    for(k in ctrl_cols){
      norm = quant_table[i,k]/norm_factors[k]
      ctrl_norms = c(ctrl_norms, norm)
    }
    control_norm = sum(ctrl_norms)/length(ctrl_norms)
    
    FCs = c()
    for (k in exp_cols){
      fc = (quant_table[i,k]/norm_factors[k] + 1/aveN) / (control_norm + 1/aveN)
      FCs = c(FCs, fc)
    }
    
    final = mean(unlist(FCs))
    quant_table[i,"FCA"] = log(final)
  }
  
  # Normalize all the fold changes to be between 0 and 1
  scaled = range01(as.numeric(quant_table$FCA))
  quant_table$scaled_FC = scaled
  
  write.table(quant_table, file = paste0(base_path, "/intermediate/interaction_scoring/FCA/", int_type, "_FCA.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}


get_FCB <- function(quant_table, int_type){
  
  quant_table$FCB = "NA"
  
  # Calculate normalization factors for each experiment
  norm_factors = c(NA, NA)
  
  exp_cols = which(grepl(int_type,names(quant_table)))
  ctrl_cols = which(grepl("control",names(quant_table)))
  all_cols = c(exp_cols, ctrl_cols)
  
  
  for (j in all_cols){
    norm_factor = sum(as.numeric(unlist(quant_table[,j])))
    norm_factors = c(norm_factors, norm_factor)
    aveN = mean(norm_factors[c(ctrl_cols)])
  }
  
  # Find the three highest spectral counts across control experiments
  highest = c()
  for (j in ctrl_cols){
    norm_cts = as.vector(quant_table[,j]) / norm_factors[j]
    highest = c(highest, tail(sort(norm_cts), n=3))
  }
  
  top3 = tail(sort(norm_cts), n=3)
  control_norm = mean(top3)
  
  # Calculate fold change for each protein
  for (i in c(1:nrow(quant_table))){
    
    FCs = c()
    for (k in exp_cols){
      fc = (quant_table[i,k]/norm_factors[k] + 1/aveN) / (control_norm + 1/aveN)
      FCs = c(FCs, fc)
    }
    
    final = mean(unlist(FCs))
    quant_table[i,"FCB"] = log(final)
  }
  
  # Normalize all the fold changes to be between 0 and 1
  scaled = range01(as.numeric(quant_table$FCB))
  quant_table$scaled_FC = scaled
  
  write.table(quant_table, file = paste0(base_path, "/intermediate/interaction_scoring/FCB/", int_type, "_FCB.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

########################################
# Run saintExpress
########################################
run_saintx <- function(int_file, prey_file, bait_file, out_file){
  path2saint = paste0(base_path, "/code/SAINTexpress_v3.6.1__2015-05-03/bin/SAINTexpress-spc")
  
  setwd(paste0(base_path, "/intermediate/interaction_scoring/saintExpress/"))

  cmd = paste(path2saint, int_file, prey_file, bait_file)
  
  system(cmd)
  
  cmd = paste0("mv list.txt ", out_file)
  system(cmd)
}

########################################
# Run saintq
########################################

write_params <- function(normyn, input_file, input_level, param_out){
  if(input_level == "protein"){
    content = paste0("normalize_control=",normyn,"\ninput_filename=",input_file,"\ninput_level=",input_level,
                     "\nprotein_colname=Proteins\ncompress_n_ctrl=100\ncompress_n_rep=100")
  }else if (input_level == "peptide"){
    content = paste0("normalize_control=",normyn,"\ninput_filename=",input_file,"\ninput_level=",input_level,
                     "\nprotein_colname=Proteins\npep_colname=Sequence\ncompress_n_ctrl=100\ncompress_n_rep=100\n",
                     "min_n_pep=3\nbest_prop_pep=0.5")
  }
  
  write.table(content, file = param_out, quote = FALSE, col.names = FALSE, row.names = FALSE)
}

run_saintq <- function(wd, paramfileName, infileName, newName){
  
  setwd(wd)
  path2saintq = paste0(base_path, "/code/saintq/bin/saintq")
  cmd = paste(path2saintq, paramfileName)
  system(cmd)
  
  out_file = paste0("scores_list__", infileName,"__.tsv")
  cmd = paste0("mv ", out_file, " ", newName)
  system(cmd)
  
}
