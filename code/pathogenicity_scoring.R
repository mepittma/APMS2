base_path = "/Users/student/Documents/PollardLab/APMS2"
out_path = paste0(base_path, "/output/pathogenicity_scoring")

source(paste0(base_path, "/code/load_data/load_variants.R"))

# Combine variants into a single object
cols <- c("Blinded.ID", "Cardiac.Category", "EM", "NDD", "CHR", "POS", "REF", "ALT", "Gene","pLI.Score")
mut_table <- rbind(DNV_cases[,cols], LoF_cases[,cols])

######################## Step 1: Select genes in interactomes with LoF and DNV ######################## 

for (analysis_type in c("original", "sainq_n")){
  ints = read.csv(paste0(base_path, "/intermediate/interactome_lists/",
                           analysis_type,"/combined_APMS_interactome_G001_T05_N1.csv"),
                    stringsAsFactors = FALSE, header = TRUE)
  int_genes = unique(ints$Prey_genename)
  int_genes = unique(unlist(strsplit(int_genes, "; ")))
  
  candidates = int_genes[which(int_genes %in% mut_table$Gene)]
  gene_data = mut_table[which(mut_table$Gene %in% candidates),]
  
  ######################## Step 2: Number of variants, normalized by length ######################## 
  exon_coords <- read.table(paste0(base_path, "/input/databases/gene_start_stop.txt"), 
                            sep = "\t", header = TRUE)
  exon_coords$CDS.Length <- NULL
  exon_coords$geneLength <- exon_coords$Gene.end..bp. - exon_coords$Gene.start..bp.
  exon_coords <- unique(exon_coords)
  
  freq_table <- as.data.frame(table(gene_data$Gene))
  freq_table$geneLength <- exon_coords$geneLength[match(freq_table$Var1, exon_coords$Gene.name)]
  freq_table$mutperkb <- (freq_table$Freq / freq_table$geneLength) * 1000
  
  gene_data$mutperkb <- freq_table$mutperkb[match(gene_data$Gene, freq_table$Var1)]
  
  # Normalize the mutperkb score to be between 0.5 and 1 - we don't want to penalize
  # rarely-mutated genes too harshly
  range51 <- function(x){(0.5*(x-min(x))/(max(x)-min(x)) + 0.5)}
  gene_data$mutperkb[is.na(gene_data$mutperkb)] = 0.000001
  gene_data$mutperkb = gene_data$mutperkb + 0.000001
  gene_data$norm_mutperkb = range51(log(gene_data$mutperkb))
  
  ######################## Step 3: PhyloP conservation score ######################## 
  
  ### Note : is there a better way to handle indels?
  
  # Load in PhyloP data
  phylop <- read.table(paste0(base_path, "/input/databases/phyloP.txt"), 
                       sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  names(phylop) = c("CHR", "POS", "STOP", "phyloP")
  phylop = phylop[which(phylop$POS %in% gene_data$POS),]
  
  # Create stop column in gene data
  gene_data$STOP = 0
  for (i in c(1:nrow(gene_data))){
    gene_data$STOP[i] = max(c(nchar(gene_data$REF[i]),
                              nchar(gene_data$ALT[i]))) + as.numeric(gene_data$POS[i]) - 1
  }
  
  # Merge by Chr/Start/Stop data 
  mg = merge(gene_data, phylop, by = c("CHR", "POS"), all.x=TRUE)
  median_phylop = median(as.numeric(mg$phyloP[!is.na(mg$phyloP)]))
  mg$phyloP[is.na(mg$phyloP)] <- median_phylop
  
  # Normalize between 0 and 1
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  mg$phyloP = range01(as.numeric(mg$phyloP))
  mg$phyloP[which(mg$phyloP == 0)] <- min(mg$phyloP[which(mg$phyloP != 0)])/2 # remove 0 values
  
  #Add to gene_data
  gene_data$phyloP <- mg$phyloP
  
  ######################## Step 4: Cell-specific expression score ######################## 
  
  broad_expr = read.csv(paste0(base_path, "/input/rnaseq/BroadPopulations_AverageExpression.csv"))
  myo_expr = read.csv(paste0(base_path, "/input/rnaseq/Myocardium_AverageExpression.csv"))
  
  broad_expr$Myocardium <- NULL
  row.names(myo_expr) = myo_expr[,1]
  myo_expr[,1] <- NULL
  myo_expr$max <- NULL
  myo_expr$max = apply(myo_expr, 1, function(x) max(x))
  
  broad_expr$myocardium = myo_expr$max[match(broad_expr$X, row.names(myo_expr))]
  broad_expr$myo_spec = broad_expr$myocardium / (rowSums(broad_expr[,c(2:7)]))
  broad_expr$myo_spec[is.nan(broad_expr$myo_spec)] = 0
  broad_expr$myo_spec[is.infinite(broad_expr$myo_spec)] = 200
  
  # Homology between mouse and human gene names
  homologene = read.table(paste0(base_path, "/input/databases/homologene.txt"), 
                          sep="\t", header = TRUE, stringsAsFactors = FALSE)
  hom_human = homologene[which(homologene$NCBI.Taxon.ID == "9606"),]
  hom_mouse = homologene[which(homologene$NCBI.Taxon.ID == "10090"),]
  homologene = merge(hom_human, hom_mouse, by = "HomoloGene.ID")
  broad_expr$Gene = homologene$Symbol.x[match(broad_expr$X, homologene$Symbol.y)]
  
  ### Note: this is not ideal; is there a better homology conversion method?
  library(dplyr)
  broad_expr$Gene[is.na(broad_expr$Gene)] = toupper(broad_expr$X[is.na(broad_expr$Gene)])
  
  gene_data$specificity_score = broad_expr$myo_spec[match(gene_data$Gene, broad_expr$Gene)]
  median_spec = median(gene_data$specificity_score[!is.na(gene_data$specificity_score)])
  gene_data$specificity_score[is.na(gene_data$specificity_score)] <- median_spec
  
  # Normalize between 0 and 1
  gene_data$specificity_score = range01(log(as.numeric(gene_data$specificity_score)))
  gene_data$specificity_score[which(gene_data$specificity_score == 0)] <- min(gene_data$specificity_score[which(gene_data$specificity_score != 0)])/2 # remove 0 values
  
  ######################## Step 5: CADD scores ########################
  
  # Save out vcf file to run through CADD
  #write.table(gene_data[,c("CHR", "POS", "Blinded.ID", "REF", "ALT")], 
  #            paste0(base_path, "/intermediate/CADD/vcf_CADD_", analysis_type, "_input.vcf"),
  #            row.names = FALSE, quote = FALSE, sep = "\t")
  
  #### MANUAL: UPLOAD FILE TO https://cadd.gs.washington.edu/score ####
  
  CADD = read.table(paste0(base_path, "/intermediate/CADD/vcf_CADD_", analysis_type, "_output.tsv"),
             sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  names(CADD) = c("CHR", "POS", "REF", "ALT", "Raw.score", "CADD")
  gd = merge(gene_data, CADD)
  gd$scaled.CADD = range01(as.numeric(gd$Raw.score))
  gd$scaled.CADD[which(gd$scaled.CADD == 0)] <- min(gd$scaled.CADD[which(gd$scaled.CADD != 0)])/2 # remove 0 values
  
  gene_data = gd[,c("Blinded.ID", "Cardiac.Category", "NDD", "CHR", "POS", "REF", "ALT", "Gene", 
                    "pLI.Score", "norm_mutperkb", "phyloP", "specificity_score", "scaled.CADD")]
  
  ######################## Step 6: Combine scores and sorts ########################
  
  # Replace missing pLI with median)
  median_pLI = median(as.numeric(gene_data$pLI.Score[!is.na(gene_data$pLI.Score)]))
  gene_data$pLI.Score[is.na(gene_data$pLI.Score)] <- median_pLI
  
  # Add pseudocount to pLI
  gene_data$pLI.Score[which(as.numeric(gene_data$pLI.Score) == 0)] <- min(as.numeric(gene_data$pLI.Score[which(as.numeric(gene_data$pLI.Score) != 0)]))/2
  
  # Reduce dataframe
  numerical_data = as.data.frame(
    gene_data[,c("pLI.Score","norm_mutperkb","specificity_score","phyloP","scaled.CADD")])
  df2 <- mutate_all(numerical_data, function(x) as.numeric(as.character(x)))
  df2$mult_score = apply(df2, 1, prod)
  df2$add_score = rowSums(df2[,c(1:5)])
  
  final <- cbind(gene_data[,c("Gene","Blinded.ID","Cardiac.Category","NDD", "CHR","POS", "REF", "ALT")], df2)
  final_mult <- final[order(-final$mult_score),]
  final_add <- final[order(-final$add_score),]
  
  ofile = paste0(base_path, "/output/pathogenicity_scoring/", analysis_type, 
                 "/all_interactomes_median_permuted_scores_G001_T05_N1.csv")
  
  
  write.csv(final_mult, ofile,
            quote = FALSE, row.names = FALSE)
  
  # Create a file that instead uses ranks
  final$pLI.rank = rank(final$pLI.Score)
  final$mutperkb.rank = rank(final$norm_mutperkb)
  final$specificity.rank = rank(final$specificity_score)
  final$phyloP.rank = rank(final$phyloP)
  final$CADD.rank = rank(final$scaled.CADD)
  
  # lower rank = lower number
  # highest numbers are most interesting
  final$AddedRank = rowSums(final[,c("pLI.rank", "mutperkb.rank", "specificity.rank", "phyloP.rank", "CADD.rank")])
  final_rank <- final[order(-final$AddedRank),c("Gene","Blinded.ID","Cardiac.Category","NDD", "CHR","POS", "REF", "ALT",
                                                "pLI.rank", "mutperkb.rank", "specificity.rank", "phyloP.rank", "CADD.rank",
                                                "AddedRank")]
  ofile = paste0(base_path, "/output/pathogenicity_scoring/", analysis_type, 
                 "/all_interactomes_ranked_G001_T05_N1.csv")
  
  write.csv(final, ofile, quote = FALSE, row.names = FALSE)
  
  ######################## GATA4/TBX5 only ########################
  
  spec_genes = unique(ints$Prey_genename[which(ints$Bait == "GATA4" | ints$Bait == "TBX5")])
  GATAs = final_mult[which(final_mult$Gene %in% spec_genes),]
  
  ofile = paste0(base_path, "/output/pathogenicity_scoring/",analysis_type, 
                 "/GATA4-TBX5_interactomes_median_permuted_scores_G001_T05.csv")
  write.csv(GATAs, ofile, quote = FALSE, row.names = FALSE)
  
  # For ranked scores
  GATAs = final[which(final_mult$Gene %in% spec_genes),]
  ofile = paste0(base_path, "/output/pathogenicity_scoring/",analysis_type, 
                 "/GATA4-TBX5_interactomes_ranked_G001_T05.csv")
  write.csv(GATAs, ofile,quote = FALSE, row.names = FALSE)
  
  ######################## Nearby variants? ########################
  res = c()
  for(chr in c(1:22, "X")){
    var_list = final_mult[which(final_mult$CHR == chr),]
    if(nrow(var_list) > 1){
      
      # Sort numerically
      coords = sort(as.numeric(var_list$POS))
      genes = var_list[match(coords, var_list$POS),]
      
      # Check whether any mutations next to each other are within 800bp in either direction
      #### This is janky and assumes that each variant has a position not shared on any other chromosome
      for (i in c(1:(length(coords)-1))){
        if(coords[i] - coords[i+1] <= 800 && coords[i] - coords[i+1] >= -800){
          success= paste0("Nearby variants! chr", chr, ":", coords[i], ",", coords[i+1],
                          " in gene ", genes[i,"Gene"])
          print(success)
          res = c(res, coords[i], coords[i+1])
        }
      }
    }
  }
  
  domain = final_mult[which(final_mult$POS %in% res),]
  domain = domain[order(domain$Gene, domain$mult_score),]
  ofile = paste0(base_path, "/output/pathogenicity_scoring/",
                 analysis_type,"/nearby_variants_G001_T05_N1.csv")
  write.csv(domain, file = ofile,
            row.names = FALSE, quote = FALSE)
  
}



### Step 1.5: save out a file that includes interactome data for all variants
cols <- c("Gene","pLI.Score", "type")
DNV_cases$type = "DNV"
LoF_cases$type = "LoF"
rec_cases$type = "recessive"


for(analysis_type in c("saintq_n", "original")){
  ints = read.table(paste0(base_path, "/intermediate/interactome_lists/",
                           analysis_type,"/combined_APMS_interactome_G001_T05_N1.csv"),
                         sep = ",", header = TRUE, stringsAsFactors = FALSE)
  ints = ints[,c("Bait","Prey","Prey_genename","BFDR")]
  
  mut_data <- rbind(DNV_cases[,cols], LoF_cases[,cols], rec_cases[,cols])
  mut_data = mut_data[which(mut_data$Gene %in% ints$Prey_genename),]
  counts = as.matrix(table(mut_data$Gene, mut_data$type))
  cdf = as.data.frame(cbind(row.names(counts),counts[,c(1:3)]))
  cdf$pLI.score = mut_data$pLI.Score[match(cdf$V1, mut_data$Gene)]
  
  # Connect to baits in interactome scoring
  counts = cdf
  counts$GATA4 = 0
  counts$TBX5 = 0
  counts$NKX25 = 0
  
  GATA4_genes = unlist(strsplit(ints$Prey_genename[which(ints$Bait == "GATA4")], "; "))
  TBX5_genes = unlist(strsplit(ints$Prey_genename[which(ints$Bait == "TBX5")], "; "))
  NKX25_genes = unlist(strsplit(ints$Prey_genename[which(ints$Bait == "NKX25")], "; "))
  
  for (i in c(1:nrow(counts))){
    if(counts$V1[i] %in% GATA4_genes){
      counts$GATA4[i] = 1
    }
    if(counts$V1[i] %in% TBX5_genes){
      counts$TBX5[i] = 1
    }
    if(counts$V1[i] %in% NKX25_genes){
      counts$NKX25[i] = 1
    }
  }
  
  names(counts)[1] = "Gene"
  write.csv(counts, file = paste0(base_path, "/output/pathogenicity_scoring/", analysis_type, 
                                  "/variant_count_G001_T05_n1.csv"))
}

# Filter to only include DNVs
for (analysis_type in c("original", "saintq_n")){
  outpath = paste0(base_path, "/output/pathogenicity_scoring/", analysis_type)
  f1 = "/all_interactomes_median_permuted_scores_G001_T05_N1.csv"
  f2 = "/GATA4-TBX5_interactomes_median_permuted_scores_G001_T05.csv"
  f3 = "/nearby_variants_G001_T05_N1.csv"
  f4 = "/variant_count_G001_T05_n1.csv"
  f5 = "/all_interactomes_ranked_G001_T05_N1.csv"
  f6 = "/GATA4-TBX5_interactomes_ranked_G001_T05.csv"
  
  for(infile in c(f1, f2, f3, f4, f5, f6)){
    tab = read.csv(paste0(outpath, infile))
    DNV_tab = tab[which(tab$POS %in% DNV_cases$POS),]
    write.csv(DNV_tab, file = paste0(base_path, "/output/pathogenicity_scoring/", analysis_type,
                                     "/DNV_only/", infile))
  }
  
}


