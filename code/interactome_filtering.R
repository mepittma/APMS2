base_path = "/Users/student/Documents/PollardLab/APMS2"

source(paste0(base_path, "/code/functions/permutation_functions.R"))
#get_fisher, get_OR, get_prots, get_sig_compl, perm_viz, permute_status
source(paste0(base_path, "/code/functions/score_interactomes.R"))
#get_comps, get_prots, get_fisher, iref_response_list

# Read in corum data
corum <- read.table(paste0(base_path, "/input/databases/corum_human_CytoscapeFormatted.txt"), 
                    sep = "\t", quote = "", fill = TRUE, header = TRUE, stringsAsFactors = FALSE)
corum <- as.data.frame(corum[,c(1:8)])

# Read in variant data
source(paste0(base_path, "/code/load_data/load_variants.R"))

# Read in alias conversion file
unimap = read.table(paste0(base_path, "/input/aliases/combined_unimap.txt"),header = TRUE, sep = "\t", stringsAsFactors = FALSE)
unimap <- unique(unimap)

for(int_type in c("GATA4", "NKX25", "TBX5")){
  
  for(analysis_type in c("original", "saintq_n")){
    
    ## Read in interactome data
    if(analysis_type == "original"){
      ints = read.table(paste0(base_path, "/input/precomp_interactomes/interactomes.tsv"), 
                         sep = "\t",fill = TRUE, header = TRUE, stringsAsFactors = FALSE)
      tab = ints[which(ints$Bait == int_type),]
    } else{
      if (int_type == "GATA4"){
        tab = read.table(file = paste0(base_path, "/intermediate/interaction_scoring/saintq/GATA4_protein_norm_true.txt"), 
                         sep = "\t", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, comment.char = "")
        BFDR_cutoff = 0.001
      } else if (int_type == "TBX5") {
        tab = read.table(file = paste0(base_path,"/intermediate/interaction_scoring/saintq/", int_type, "_peptide_norm_true.txt"),
                         sep = "\t", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, comment.char = "")
        BFDR_cutoff = 0.05
      } else if (int_type == "NKX25") {
        tab = read.table(file = paste0(base_path,"/intermediate/interaction_scoring/saintq/", int_type, "_peptide_norm_true.txt"),
                         sep = "\t", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, comment.char = "")
        BFDR_cutoff = 0.1
      }
    }
    
    ## Add genename column
    names(tab)[which(names(tab) == "Prey_proteinname")] <- "Prey"
    tab$Prey_genename = unimap$GeneSymbol[match(tab$Prey, unimap$UniProt)]
    
    ## Filter generously
    tab_sig = subset(tab, tab$BFDR < 0.1)
    
    ## Remove any genes in blacklists
    non_nuclear = read.table(paste0(base_path, "/input/databases/localization_non_nuclear_BINGO_.txt"),
                             sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    non_nuclear = non_nuclear[which(!non_nuclear$subcellular.location == "Nucleus"),]
    if(length(tab_sig$Prey[tab_sig$Prey %in% non_nuclear$UniprotID]) > 0){
      print(paste0("non-nuclear gene found: ", int_type, "interactome analyzed by ", 
                   analysis_type, ", genes: ", tab_sig$Prey[tab_sig$Prey %in% non_nuclear$UniprotID]))
    }
    tab_sig = tab_sig[which(!tab_sig$Prey %in% non_nuclear$UniprotID),]

    blacklist = read.table(paste0(base_path, "/intermediate/rnaseq/", int_type, "_negativeFC_blacklist.txt"))
    if(length(tab_sig$Prey_genename[tab_sig$Prey_genename %in% blacklist]) > 0){
      print(paste0("blacklisted gene found: ", int_type, "interactome analyzed by ", 
                   analysis_type, ", genes: ", tab_sig$Prey_genename[tab_sig$Prey_genename %in% blacklist]))
    }
    tab_sig = tab_sig[which(!tab_sig$Prey_genename %in% blacklist),]
    
    # Save out APMS interactome
    write.csv(tab_sig, 
              paste0(base_path, "/intermediate/interactome_lists/", analysis_type, "/", int_type, "_APMS_interactome_cutoff_1.csv"),
              row.names = FALSE, quote= FALSE)

    # Save out more specific filters
    if(int_type %in% c("GATA4", "TBX5")){
      tab_sig = subset(tab, tab$BFDR < 0.05)
      write.csv(tab_sig, 
                paste0(base_path, "/intermediate/interactome_lists/", analysis_type, "/", int_type, "_APMS_interactome_cutoff_05.csv"),
                row.names = FALSE, quote= FALSE)
    }
    
    if(int_type == "GATA4"){
      tab_sig = subset(tab, tab$BFDR < 0.001)
      write.csv(tab_sig, 
                paste0(base_path, "/intermediate/interactome_lists/", analysis_type, "/", int_type, "_APMS_interactome_cutoff_001.csv"),
                row.names = FALSE, quote= FALSE)
    }
  } 
}
  


# Create single file with more specific filtering
for(analysis_type in c("original", "saintq_n")){
  
  Gtab = read.csv(paste0(base_path, "/intermediate/interactome_lists/", analysis_type, "/GATA4_APMS_interactome_cutoff_1.csv"))
  Ttab = read.csv(paste0(base_path, "/intermediate/interactome_lists/", analysis_type, "/TBX5_APMS_interactome_cutoff_1.csv"))
  Ntab = read.csv(paste0(base_path, "/intermediate/interactome_lists/", analysis_type, "/NKX25_APMS_interactome_cutoff_1.csv"))
  
  Ttab$X.Pep = NULL
  Ntab$X.Pep = NULL
  
  Gtab = Gtab[which(Gtab$BFDR <= 0.001),]
  Ttab = Ttab[which(Ttab$BFDR <= 0.05),]
  Ntab = Ntab[which(Ntab$BFDR <= 0.1),]
  
  comb_tab = rbind(Gtab, Ttab)
  comb_tab = rbind(comb_tab, Ntab)
  
  write.csv(comb_tab, 
            paste0(base_path, "/intermediate/interactome_lists/", analysis_type, "/combined_APMS_interactome_G001_T05_N1.csv"),
            row.names = FALSE, quote= FALSE)
}

# Corum expansion for more specific filtering - one file listing all protein/gene names, one file with all interactions
for(analysis_type in c("original", "saintq_n")){
  
  for (int_type in c("GATA4","TBX5","NKX25")){
    if(int_type == "GATA4"){
      tab = read.csv(paste0(base_path, "/intermediate/interactome_lists/", analysis_type, "/GATA4_APMS_interactome_cutoff_001.csv"))
      cutoff_str = "G001"
    } else if(int_type == "TBX5") {
      tab = read.csv(paste0(base_path, "/intermediate/interactome_lists/", analysis_type, "/TBX5_APMS_interactome_cutoff_05.csv"))
      cutoff_str = "T05"
    } else if (int_type == "NKX25"){
      tab = read.csv(paste0(base_path, "/intermediate/interactome_lists/", analysis_type, "/NKX25_APMS_interactome_cutoff_1.csv"))
      cutoff_str = "N1"
    }
    
    # Select corum interactors from tab list
    complexes = get_comps(tab$Prey,corum, enriched=FALSE)
    corum_subset = corum[which(corum$Complex.id %in% complexes),]
    
    # Combine into a single dataset with column for source
    tab$source = "APMS"
    tab$X.Pep = NULL
    names(tab) = c("GeneName1", "Interactor2", "x.Rep", "AvgP", "BFDR", "GeneName2", "source")
    tab$Interactor1 = unimap$UniProt[match(tab$GeneName1, unimap$GeneSymbol)]
    
    if(nrow(corum_subset) > 0){
      corum_subset$source = "CORUM"
      corum_subset$GeneName1 = unimap$GeneSymbol[match(corum_subset$Interactor1, unimap$UniProt)]
      corum_subset$GeneName2 = unimap$GeneSymbol[match(corum_subset$Interactor2, unimap$UniProt)]
      expanded = rbind(corum_subset[,c("Interactor1", "Interactor2", "GeneName1", "GeneName2", "source")],
                       tab[,c("Interactor1", "Interactor2", "GeneName1", "GeneName2", "source")])
    } else {
      expanded = tab[,c("Interactor1", "Interactor2", "GeneName1", "GeneName2", "source")]
    }
    
    
    write.csv(expanded, 
              paste0(base_path, "/intermediate/interactome_lists/", analysis_type, "/", 
                     int_type, "_corum_expanded_interactome_" ,cutoff_str, ".csv"),
              row.names = FALSE, quote= FALSE)
  }

}

# Write out combined corum-expanded interactomes
for(analysis_type in c("original", "saintq_n")){
  GATA = read.csv(paste0(base_path, "/intermediate/interactome_lists/", analysis_type, "/GATA4", 
                         "_corum_expanded_interactome_G001.csv"), stringsAsFactors = FALSE)
  TBX5 = read.csv(paste0(base_path, "/intermediate/interactome_lists/", analysis_type, "/TBX5", 
                         "_corum_expanded_interactome_T05.csv"), stringsAsFactors = FALSE)
  NKX25 = read.csv(paste0(base_path, "/intermediate/interactome_lists/", analysis_type, "/NKX25", 
                          "_corum_expanded_interactome_N1.csv"), stringsAsFactors = FALSE)
  
  combined_corum = rbind(GATA,TBX5)
  combined_corum = rbind(combined_corum, NKX25)
  
  write.csv(combined_corum, 
            paste0(base_path, "/intermediate/interactome_lists/", analysis_type, 
                   "/combined_corum_expanded_interactome_G001T05N1.csv"),
            row.names = FALSE, quote= FALSE)
}

