base_path = "/Users/student/Documents/PollardLab/APMS2"
library(pROC)

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

# Read in RNAseq data
source(paste0(base_path, "/code/load_data/load_rnaseq.R"))

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
# Run quantification/pre-processing
####################################################
# APMS data for the three interactomes
for (int_type in c("TBX5", "NKX25", "GATA4")){
  
  ev_file = paste0(base_path, "/input/evidence/expTFs/", int_type, "_evidence.txt")
  key_file= paste0(base_path, "/input/evidence/expTFs/keys_", int_type, ".txt")
  
  # Creates files necessary to run FC, saintq, and saint analysis - uncomment if these files do not exist
  #create_artMSfiles(int_type, key_file)
  #setwd(paste0(base_path, "/intermediate/MS_QC/", int_type, "_msspc"))
  #run_artMS(int_type, ev_file, key_file, "msspc")
  #setwd(paste0(base_path, "/intermediate/MS_QC/", int_type, "_msint"))
  #run_artMS(int_type, ev_file, key_file, "msint")
  
  # Load in APMS quantification file - SAINTq msspc file
  quant = read.table(paste0(base_path, "/input/evidence/expTFs/saintq_inputs/", int_type, "/saintq_input_proteins.txt"), 
                     sep = "\t", skip = 2 , header = TRUE)
  
  # Run FCA and FCB to save those intermediate files
  get_FCA(quant, int_type)
  get_FCB(quant, int_type)
  
  # Run saintExpress for peptide and protein
  for (spec_type in c("msspc", "msint")){
    int_file = paste0(base_path, "/input/evidence/expTFs/saintx_inputs/", int_type, "_", spec_type, "-saint-interactions.txt")
    prey_file = paste0(base_path, "/input/evidence/expTFs/saintx_inputs/",int_type, "_", spec_type, "-saint-preys.txt")
    bait_file = paste0(base_path, "/input/evidence/expTFs/saintx_inputs/",int_type, "_", spec_type, "-saint-baits.txt")
    
    out_file = paste0(base_path, "/intermediate/interaction_scoring/saintExpress/", int_type, "_", spec_type, ".txt")
    run_saintx(int_file, prey_file, bait_file, out_file)
  }
  
  
  # Run saintq with all possible parameter files 
  for(level in c("peptide", "protein")){
    
    setwd(paste0(base_path, "/input/evidence/expTFs/saintq_inputs/", int_type))
    in_file = paste0("saintq_input_", level, "s.txt")
    
    for (normYN in c("true", "false")){
      param_out = paste0( level, "_norm", normYN, ".txt")
      write_params(normYN, in_file, level, param_out)
      out_file = paste0(base_path, "/intermediate/interaction_scoring/saintq/", int_type, "_", level, "_norm_", normYN, ".txt")
      run_saintq(paste0(base_path, "/input/evidence/expTFs/saintq_inputs/", int_type),param_out, in_file, out_file)
    }
    
  }
  
}

####################################################
# Create blacklists for proteins with low expression in KO
####################################################

# ONLY TBX5 and NKX25 as of 3/7/2019

## Analysis for each interactome type
for (int_type in c("TBX5", "NKX25")){
  
  # Subset - cpm data
  ctrl_cpm <- cpm[,c(2:6)]
  indx <- sapply(cpm[,c(2:ncol(cpm))], is.factor)
  cpm[indx] <- lapply(cpm[indx], function(x) as.numeric(as.character(x)))
  cpm_cols = grepl(paste0(substr(int_type,1,1), "KO"), names(cpm))
  exp_cpm = cpm[,which(cpm_cols)]
  
  # Subset - raw data
  ctrl_seq <- seqdata[,c(2:6)]
  exp_seq = seqdata[,which(cpm_cols)]
  countdata = cbind(ctrl_seq, exp_seq)
  
  # Identify genes with at least 0.5 cpm in at least 2 samples
  temp = cbind(ctrl_cpm, exp_cpm)
  thresh <- temp > 0.5
  keep <- rowSums(thresh) >= 2
  counts.keep <- countdata[keep,]
  
  # Convert to an edgeR object
  dgeObj <- DGEList(counts.keep)
  dgeObj <- calcNormFactors(dgeObj) #TMM normalization
  group = c("WT","WT","WT","WT","WT","KO","KO","KO","KO","KO")
  design = model.matrix(~group)
  
  # MDS
  #plotMDS(dgeObj, labels=group, cex=0.75, xlim=c(-4, 5))
  
  # Estimating dispersion
  dgeObj <- estimateCommonDisp(dgeObj)
  dgeObj <- estimateGLMTrendedDisp(dgeObj)
  dgeObj <- estimateTagwiseDisp(dgeObj)
  #plotBCV(dgeObj)
  
  # Testing for differential expression
  fit <- glmFit(dgeObj, design)
  lrt.BvsL <- glmLRT(fit, coef=2)
  tab = as.data.frame(lrt.BvsL$table)
  tab$FDR = p.adjust(tab$PValue, method = "fdr", n = length(tab$PValue))
  sigmaybe <- tab[which(tab$FDR <= 0.05 & tab$logFC <= -0.5),]
  
  # Save out blacklist
  blacklist = row.names(sigmaybe)
  write.table(blacklist, file = paste0(base_path, "/intermediate/rnaseq/", int_type, "_negativeFC_blacklist.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  
}



####################################################
# Score interactome lists for all parameter choices
####################################################

for(int_type in c("GATA4", "TBX5", "NKX25")){
  
  # Create a 34x3 matrix to score results
  method_vec = c()
  DNV_vec = c()
  LoF_vec = c()
  
  for(FDR_cutoff in c(0.05, 0.1)){
    
    for(KO_status in c("Y", "N")){
      
      ### FC ###
      
      for(FC in c("FCA", "FCB")){
        
        # Read in data, get gene names
        tab = read.table(file = paste0(base_path, "/intermediate/interaction_scoring/",FC, "/", int_type, "_", FC,".txt"),
                             sep = "\t", header = TRUE, stringsAsFactors = FALSE)
        tab$Genes = quant_unimap$GeneSymbol[match(tab$Proteins, quant_unimap$UniProt)]
        
        # Choose significant genes
        tab_sig = subset(tab, tab[,FC] > quantile(tab[,FC], prob = 1 - FDR_cutoff))
        
        # Corum expansion
        complexes = get_comps(tab_sig$Proteins,corum)
        interactome_prots = unique(c(get_prots(complexes, corum), tab_sig$Proteins))
        
        # Get gene names for interactome
        interactome_genes = quant_unimap$GeneSymbol[match(interactome_prots, quant_unimap$UniProt)]
        interactome_genes = unlist(strsplit(interactome_genes, "; "))
        
        # Remove blacklist genes from interactome and coerce FC score to 0
        if (KO_status == "Y"){
          blacklist = read.table(paste0(base_path, "/intermediate/rnaseq/", int_type, "_negativeFC_blacklist.txt"))
          tab$scaled_FC[which(tab$Genes %in% blacklist)] <- 0
          interactome_genes <- interactome_genes[which(!interactome_genes %in% blacklist)]
        }
        
        # Get true (0 = no interaction, 1 = interaction)
        int_prot = quant_unimap$UniProt[which(quant_unimap$GeneSymbol == int_type)]
        #response = iref_response_list(int_type, tab$Proteins, iref)
        
        #iRefAUC = auc(response = response, predictor=FCA_tab$scaled_FCA)
        DNV_fish = get_fisher(DNV_cases, DNV_ctrls, interactome_genes)
        LoF_fish = get_fisher(LoF_cases, LoF_ctrls, interactome_genes)
        
        # Save into results table
        method = paste0(FC, "_FDR", FDR_cutoff, "_KOfilt", KO_status)
        method_vec = c(method_vec, method)
        DNV_vec = c(DNV_vec, DNV_fish[1])
        LoF_vec = c(LoF_vec, LoF_fish[1])
      }
      
      ######################
      ### Original lists ###
      tab = read.csv(paste0(base_path, "/input/precomp_interactomes/interactomes.csv"), sep = "\t", 
                     stringsAsFactors = FALSE)
      interactome_genes = tab$Prey_genename[which(tab$Bait == int_type)]
      
      # Remove blacklist genes from interactome and coerce FC score to 0
      if (KO_status == "Y"){
        blacklist = read.table(paste0(base_path, "/intermediate/rnaseq/", int_type, "_negativeFC_blacklist.txt"))
        interactome_genes <- interactome_genes[which(!interactome_genes %in% blacklist)]
      }
      
      # Fisher p
      DNV_fish = get_fisher(DNV_cases, DNV_ctrls, interactome_genes)
      LoF_fish = get_fisher(LoF_cases, LoF_ctrls, interactome_genes)
      
      # Save into results table
      method = paste0("original_KOfilt", KO_status)
      method_vec = c(method_vec, method)
      DNV_vec = c(DNV_vec, DNV_fish[1])
      LoF_vec = c(LoF_vec, LoF_fish[1])
      
      ############################### saintx
      for(level_type in c("protein", "peptide")){
        
        # Read in data, get gene names
        if(level_type == "protein"){
          tab = read.table(file = paste0(base_path, "/intermediate/interaction_scoring/saintExpress/", int_type, "_msspc.txt"),
                             sep = "\t", header = TRUE, stringsAsFactors = FALSE)
        } else if (level_type == "peptide"){
          tab = read.table(file = paste0(base_path, "/intermediate/interaction_scoring/saintExpress/", int_type, "_msint.txt"),
                           sep = "\t", header = TRUE, stringsAsFactors = FALSE)
        }
        tab$Genes = quant_unimap$GeneSymbol[match(tab$Prey, quant_unimap$UniProt)]
          
        # Choose significant genes
        tab_sig = subset(tab, tab$BFDR < FDR_cutoff)
          
        # Corum expansion
        complexes = get_comps(tab_sig$Prey,corum)
        interactome_prots = unique(c(get_prots(complexes, corum), tab_sig$Prey))
          
        # Get gene names for interactome
        interactome_genes = quant_unimap$GeneSymbol[match(interactome_prots, quant_unimap$UniProt)]
        interactome_genes = unlist(strsplit(interactome_genes, "; "))
          
        # Remove blacklist genes from interactome and coerce BFDR score to 1
        if (KO_status == "Y"){
          blacklist = read.table(paste0(base_path, "/intermediate/rnaseq/", int_type, "_negativeFC_blacklist.txt"))
          tab$BFDR[which(tab$Genes %in% blacklist)] <- 1
          interactome_genes <- interactome_genes[which(!interactome_genes %in% blacklist)]
        }
          
        # Get true (0 = no interaction, 1 = interaction)
        int_prot = quant_unimap$UniProt[which(quant_unimap$GeneSymbol == int_type)]
        #response = iref_response_list(int_prot, tab$Prey, iref)
          
        #iRefAUC = auc(response = response, predictor=1 - tab$BFDR)
        DNV_fish = get_fisher(DNV_cases, DNV_ctrls, interactome_genes)
        LoF_fish = get_fisher(LoF_cases, LoF_ctrls, interactome_genes)
          
        # Save into results table
        method = paste0("saintx_", level_type, "_FDR", FDR_cutoff, "_KOfilt", KO_status)
        method_vec = c(method_vec, method)
        DNV_vec = c(DNV_vec, DNV_fish[1])
        LoF_vec = c(LoF_vec, LoF_fish[1])
        
        ############################ saintq
        for(norm_status in c("Y", "N")){
          
          if(norm_status == "Y"){n = "true"
          } else{n = "false"
          }
          
          tab = read.table(file = paste0(base_path, "/intermediate/interaction_scoring/saintq/", 
                                         int_type, "_", level_type, "_norm_", n,".txt"),
                           sep = "\t", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, comment.char = "")
          
          # Choose significant genes
          tab_sig = subset(tab, tab$BFDR < FDR_cutoff)
          
          # Corum expansion
          complexes = get_comps(tab_sig$Prey,corum)
          interactome_prots = unique(c(get_prots(complexes, corum), tab_sig$Prey))
          
          # Get gene names for interactome
          interactome_genes = quant_unimap$GeneSymbol[match(interactome_prots, quant_unimap$UniProt)]
          interactome_genes = unlist(strsplit(interactome_genes, "; "))
          
          # Remove blacklist genes from interactome and coerce BFDR score to 1
          if (KO_status == "Y"){
            blacklist = read.table(paste0(base_path, "/intermediate/rnaseq/", int_type, "_negativeFC_blacklist.txt"))
            tab$BFDR[which(tab$Genes %in% blacklist)] <- 1
            interactome_genes <- interactome_genes[which(!interactome_genes %in% blacklist)]
          }
          
          # Get true (0 = no interaction, 1 = interaction)
          int_prot = quant_unimap$UniProt[which(quant_unimap$GeneSymbol == int_type)]
          #response = iref_response_list(int_prot, tab$Prey, iref)
          
          #iRefAUC = auc(response = response, predictor=1 - tab$BFDR)
          DNV_fish = get_fisher(DNV_cases, DNV_ctrls, interactome_genes)
          LoF_fish = get_fisher(LoF_cases, LoF_ctrls, interactome_genes)
          
          # Save into results table
          method = paste0("saintq_", level_type, "_FDR", FDR_cutoff, "_KOfilt", KO_status,"_norm",norm_status)
          method_vec = c(method_vec, method)
          DNV_vec = c(DNV_vec, DNV_fish[1])
          LoF_vec = c(LoF_vec, LoF_fish[1])
          
        }
        
      }
    }
    
  }
  res = cbind(method_vec, DNV_vec, LoF_vec)
  names(res) = c("Method", "DNV_pvalue", "LoF_pvalue")
  write.table(res, file = paste0(base_path, "/output/interactome_selection/", int_type, "score_table.tsv"),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}



####################################################
# Methodical subsets of parameter choices
####################################################


####################################################
# Hand-picked parameter set comparison
####################################################
