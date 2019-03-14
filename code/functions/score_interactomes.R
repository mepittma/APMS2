
#################################################
# Corum expansion functions
#################################################

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
    
    fish = get_complex_fisher(prot_list, compl, corum, total_prots)
    
    if (fish <= 0.05){
      sig_complexes = c(sig_complexes, compl)
    }
  }
  return(sig_complexes)
}

# Function to get a list of complex IDs given a list of interactors
get_comps <- function(int_list, corum, enriched = FALSE){
  
  if(enriched == TRUE){
    return(get_sig_compl(int_list, corum))
  } else {
    return(unique(c(corum$Complex.id[which(corum$Interactor1 %in% int_list)],
                    corum$Complex.id[which(corum$Interactor2 %in% int_list)])))
  }

}

# Function to isolate proteins given a complex ID
get_prots <- function(complex_list, corum){
  prots = c()
  for (compl_id in complex_list){
    prots = (unique(c(prots, corum$Interactor1[which(corum$Complex.id == compl_id)],
                    corum$Interactor2[which(corum$Complex.id == compl_id)])))
  }
  return(prots)

}

#################################################
# PCGC variant scoring
#################################################

# Get fisher p-value
get_fisher <- function(case_table, ctrl_table, gene_list){
  
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
  
  fisherMat <- matrix(c(case_ct, ctrl_ct, case_non_ct, ctrl_non_ct), nrow = 2,
                      dimnames = list(CHD = c("Yes","No"),
                                      Interactome = c("Yes", "No")))
    
  # Fisher test
  fish = fisher.test(fisherMat, alternative = 'greater')
  return(fish)
}


#################################################
# iRefIndex scoring
#################################################

iref = read.table(paste0(base_path, "/input/databases/iRefIndex.txt"), 
                  sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Checks whether a bait-prey pair exists in the database
iref_check <- function(bait, prey, db){
  # Returns TRUE if the interaction is in iRefDB, FALSE if not
  n_match = sum(db$uidA==bait & db$uidB==prey) + sum(db$uidA==prey & db$uidB==bait)
  if(n_match > 0){
    return(1)
  } else {
    return(0)
  }
}

# Creates a list of values indicating the "ground truth" of an interaction
iref_response_list <- function(bait, int_list, db){
  
  # Returns the true membership of each potential interaction
  responses = c()
  for (interactor in int_list){
    responses = c(responses, iref_check(bait, interactor, db))
  }
  
  return(responses)
}
