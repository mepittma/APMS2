########################################
# Permutation of interactome enrichment
########################################


# Function to isolate proteins given a complex ID
get_prots <- function(compl_id, corum){
  return(unique(c(corum$Interactor1[which(corum$Complex.id == compl_id)],
                  corum$Interactor2[which(corum$Complex.id == compl_id)])))
}

# Function to perform fisher exact test given input list and complex id
get_complex_fisher <- function(gList, compl_id, corum, total_prots){
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

# Function to get significant corum complex IDs given a list of proteins in an interactome
get_sig_compl <- function(prot_list, corum){
  
  total_prots = unique(c(corum$Interactor1, corum$Interactor2))
  compl_list = unique(corum$Complex.id)
  sig_complexes = c()
  
  for ( compl in compl_list ){
    
    fish = get_fisher(prot_list, compl, corum, total_prots)
    
    if (fish <= 0.05){
      sig_complexes = c(sig_complexes, compl)
    }
  }
  return(sig_complexes)
}

# Function to complete n_perm number of permutations on a given input mutation table
permute_status <- function(case_table, ctrl_table, gene_list, n_perm){
  
  perm_list = c()
  case_ids = case_table$Blinded.ID
  ctrl_ids = ctrl_table$Blinded.ID
  mut_table <- rbind(case_table[,c('Blinded.ID', "Gene")], ctrl_table[,c('Blinded.ID', "Gene")])
  
  for (i in 1:n_perm){
    
    # Create a scrambled lookup table for the ids 
    scram = mut_table$Blinded.ID[sample(length(mut_table$Blinded.ID))]
    lookup = as.data.frame(cbind(mut_table$Blinded.ID, scram))
    
    # Replace mut_table IDs with the scrambled IDs
    perm_data <- mut_table
    perm_data$Blinded.ID <- lookup$scram[match(mut_table$Blinded.ID, lookup$V1)]
    
    # Separate data by case/control ID
    perm_cases = perm_data[which(perm_data$Blinded.ID %in% case_ids),]
    perm_ctrls = perm_data[which(perm_data$Blinded.ID %in% ctrl_ids),]
    perm_list = c(perm_list, get_OR(perm_cases, perm_ctrls, gene_list))
  }
  
  return(perm_list)
}


# Function to calculate the odds ratio of a group of genes
get_OR <- function(case_table, ctrl_table, gene_list, fisher=FALSE, int_type = NULL, mut_type = NULL){
  
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
  
  if(fisher==TRUE){
    
    # Write out contingency table
    fname = paste0(out_path, "/contingency_tables/",int_type, "_", mut_type, "muts_contTable.csv")
    fisherMat <- matrix(c(case_ct, ctrl_ct, case_non_ct, ctrl_non_ct), nrow = 2,
                        dimnames = list(CHD = c("Yes","No"),
                                        Interactome = c("Yes", "No")))
    
    # Fisher test
    fish = fisher.test(fisherMat, alternative = 'greater')
    
    # Write out pvalue
    #cat(paste("Fisher exact test p-value: ", fish[1]),file=fname,append=TRUE)
  }
  
  odds = (case_ct * ctrl_non_ct) / (ctrl_ct * case_non_ct)
  return(odds)
}

############################
# Permutation visualization
############################


# Function to visualize a list of permuted odds ratios compared to the true odds ratio
perm_viz <- function(perm_list, true_OR, mut_type, int_type, n_tests){
  
  nGreater = length(perm_list[which(perm_list >= true_OR)])
  pval = (nGreater/length(perm_list))
  if(nGreater == 0 ){
    pval = "< 0.001"
  }
  
  pdf(paste0(out_path, "/", int_type, "_", mut_type, ".pdf"))
  hist(perm_list, main=paste0("Permutation of ",mut_type, 
                              " mutations in ", int_type, " interactome"), 
       sub = paste0("p: ", pval),
       xlab = "Odds Ratio")
  abline(v = true_OR, lty="dotted", lwd="5", col = "red")
  dev.off()
  
  df = as.data.frame(perm_list)
  img <- ggplot(data=df, aes(df$perm_list)) + geom_histogram() + 
    geom_vline(xintercept = true_OR, col = "red", linetype="dashed") +
    labs(title=paste0("Permutation of ",mut_type," mutations in ", int_type, " interactome"),
         subtitle = paste0("p = ", pval),
         x="Odds Ratio", y = "Frequency")

  return(img)
}

