########################################
# Functions for FC-A and FC-B
########################################
get_FCA <- function(quant_table){
  
  quant_table$FCA = "NA"
  
  # Calculate normalization factors for each experiment
  norm_factors = c(NA, NA)
  for (j in c(3:8)){
    norm_factor = sum(as.numeric(unlist(quant_table[,j])))
    norm_factors = c(norm_factors, norm_factor)
    aveN = mean(norm_factors[c(6:8)])
  }
  
  # Calculate fold change for each protein
  for (i in nrow(quant_table)){
    control_norm = (quant_table[i,6]/norm_factors[6] + quant_table[i,7]/norm_factors[7] + quant_table[i,8]/norm_factors[8])/3
    
    FCs = c()
    for (k in c(3:5)){
      fc = (quant_table[i,k]/norm_factors[k] + 1/aveN) / (control_norm + 1/aveN)
      FCs = c(FCs, fc)
    }
    
    final = mean(unlist(FCs))
    quant_table[i,9] = final
  }
  
  return(quant_table)
}


get_FCB <- function(quant_table){
  
  quant_table$FCB = "NA"
  
  # Calculate normalization factors for each experiment
  norm_factors = c(NA, NA)
  for (j in c(3:8)){
    norm_factor = sum(as.numeric(unlist(quant_table[,j])))
    norm_factors = c(norm_factors, norm_factor)
    aveN = mean(norm_factors[c(6:8)])
  }
  
  # Find the three highest spectral counts across control experiments
  highest = c()
  for (j in c(6:8)){
    norm_cts = as.vector(quant_table[,j]) / norm_factors[j]
    highest = c(highest, tail(sort(norm_cts), n=3))
  }
  
  top3 = tail(sort(norm_cts), n=3)
  control_norm = mean(top3)
  
  # Calculate fold change for each protein
  for (i in c(1:nrow(quant_table))){
    
    FCs = c()
    for (k in c(3:5)){
      fc = (quant_table[i,k]/norm_factors[k] + 1/aveN) / (control_norm + 1/aveN)
      FCs = c(FCs, fc)
    }
    
    final = mean(unlist(FCs))
    quant_table[i,9] = final
  }
  
  return(quant_table)
}

########################################
# Run saintExpress
########################################
run_saintx <- function(int_file, prey_file, bait_file, out_file, topoFile){
  topoPath = paste0(base_path, "/input/topAvg.txt")
  path2saint = paste0(base_path, "/code/SAINTexpress_v3.6.1__2015-05-03/bin/SAINTexpress-spc")
  
  setwd(paste0(base_path, "/intermediate/saintExpress/"))
  
  if(topoFile == "Y"){
    cmd = paste(path2saint, int_file, prey_file, bait_file, topoPath)
  } else if(topoFile == "N"){
    cmd = paste(path2saint, int_file, prey_file, bait_file)
  }
  
  system(cmd)
  
  cmd = paste0("mv list.txt ", out_file)
  system(cmd)
}

########################################
# Run saintq
########################################
run_saintq <- function(wd, paramfileName, infileName, newName){
  
  setwd(wd)
  path2saintq = paste0(base_path, "/code/saintq/bin/saintq")
  cmd = paste(path2saintq, paramfileName)
  system(cmd)
  
  out_file = paste0("scores_list__", infileName,"__.tsv")
  cmd = paste0("mv ", out_file, " ", newName)
  system(cmd)
  
}

########################################
# get_interactomes
########################################

# Saves out a text file with a list of genes in the interactome based on cutoff parameters
get_interactomes <- function(score_file, FDR_cutoff, )