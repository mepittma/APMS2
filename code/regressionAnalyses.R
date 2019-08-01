####MSSTATS DATA FORMATTING####
library(MSstats)
base_path = "/Users/student/Documents/PollardLab/APMS2"
setwd(paste0(base_path, "/intermediate/interaction_scoring/regression"))
source(paste0(base_path, "/code/functions/MSstats_workaround.R"))
source(paste0(base_path, "/code/functions/regression_functions.R"))

###STEP 0 - prepare data
#Convert evidence file to MSstats
for (int_type in c("GATA4", "TBX5", "NKX25")){
  ev_file = read.table(paste0(base_path, "/input/evidence/expTFs/",int_type,"_evidence.txt"),
                       sep = "\t", header=T, stringsAsFactors = F)
  an_file = read.table(paste0(base_path, "/input/evidence/expTFs/annotation_",int_type,".txt"),sep = "\t", 
                       header=T, stringsAsFactors = F)
  pr_file = read.table(paste0(base_path, "/input/evidence/expTFs/", int_type, "_proteinGroups.txt"),sep = "\t", 
                       header=T, stringsAsFactors = F)
  
  msStatsInputFiltered <- MaxQtoMSstatsFormat(ev_file, an_file, pr_file)
  #This does: format the data. remove iRT. remove truncated. convert bad qvalues into 0 (if requested). 
  #remove non unique peptides. remove with very few intensities across runs.
  
  ###STEP 1 - dataProcess
  #Run DataProcess. All hits. No normalization (the model design does a good job and keeps variability to 
  #the minimum. also difficult to find non-BAG3 interacting proteins to use as controls.)
  processedFilteredNoNorm <- dataProcess(msStatsInputFiltered, censoredInt = '0', normalization = F)
  #, featureSubset = "top3")
  
  ###STEP 3 - quantification
  #This will produce a table with summarized intensity values per protein and replicate. Note its logvalues.
  quantificationFiltered <- quantification(processedFilteredNoNorm, type="Sample")
  
  #Generates 'quant' for dowsntream analyses (model analyses and others)
  quant <- quantificationFiltered
  row.names(quant) <- quant[,1]
  quant <- quant[-1]
  colnames(quant) <- gsub("\\-", "_", colnames(quant))
  colnames(quant) <- gsub("\\+", "", colnames(quant))
  quantBackup <- quant
  
  ####MODEL ANALYSIS FOR INTERACTORS (vs controls). Captures replicate variability much better!####
  #FIRST OF ALL, imputing and/or cleaning data?
  quant <- quantBackup
  #quant <- remove_low_replicate(quant, minReplicates = 2, missing=NA, replaceWith = 0) 
  #Replace by 0 those values for proteins that are present in less than minReplicates. Function code at bottom
  quant <- impute_data(quant, missingVal=0, downshift = 3, byColumn = TRUE, width = 0.2) 
  #IMPUTE NA DATA - REMOVE IF NECESSARY. impute those that are NA or 0. 
  boxplot(quant)
  
  samples <- colnames(quant)
  variant <- gsub("_.*$", "", samples)
  variant <- as.factor(variant)
  variant <- relevel(variant, ref = "control") #set control
  
  replicate = c()
  idx = 1
  for(i in range(1,length(unique(variant)))){
    n = sum(variant==variant[idx])
    replicate = c(replicate, rep(1:n))
    idx = idx + n
  }
  
  replicate <- as.factor(replicate)
  phenoData <- data.frame(cbind(Samples=samples, Variant=variant, Replicate=replicate))
  # head(PhenoData)
  
  #Limma helps with the model functions.
  # source("https://bioconductor.org/biocLite.R")
  # biocLite("limma")
  require(limma)
  
  # The following model controls for replicate variability. 
  # It also identifies all hits significantly different from controls.
  design <- model.matrix(~0+variant+replicate)
  fit <- lmFit(quant, design=design) #Note this is with log values!
  
  #Here  we get significant factors for each one 'individually'.
  #Iterate through baits.
  combinedTable <- data.frame()
  for(condition in colnames(design)[grepl("variant", colnames(design))]){
    #if(condition=="variantcontrol"){next}
    print(condition)
    #contrastM <- makeContrasts(paste0(condition, "-variantControl"), levels=design)
    #contrastM <- makeContrasts(paste0(condition, "+variantBAG3_WT-2*variantControl"), levels=design)
    contrastM <- makeContrasts(condition, levels=design)
    contrastFit <- contrasts.fit(fit, contrastM)
    contrastFit <- eBayes(contrastFit, proportion=0.3, trend=TRUE)
    t <- topTable(contrastFit, number=200)
    t$gene <- rownames(t)
    t$condition <- condition
    combinedTable <- rbind(combinedTable, t)
  }
  
  interactorTable <- combinedTable[combinedTable$adj.P.Val<0.05 & combinedTable$logFC>=0,] 
  print(nrow(interactorTable))
  #(also filter out things where theres more on control than bait, though it seems none  of the significant 
  # ones were like that)
  
  # Convert UPIDs to genes
  quant_unimap = read.table(paste0(base_path, "/input/aliases/combined_unimap.txt"),
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  quant_unimap <- unique(quant_unimap)
  interactorTable$GeneName = quant_unimap$GeneSymbol[match(interactorTable$gene, quant_unimap$UniProt)]
  combinedTable$GeneName = quant_unimap$GeneSymbol[match(combinedTable$gene, quant_unimap$UniProt)]
  write.table(interactorTable, file = paste0(base_path, "/intermediate/interaction_scoring/regression/", 
                                             int_type, "_regression.txt"), sep = "\t", row.names = F)
}


####MODEL ANALYSIS FOR WT VS VARIANTS####

# Load in variant line data
ev_file = read.table(paste0(base_path, "/input/evidence/expTFs/g296s_vs_wt_evidence.txt"), sep = "\t", 
                     header=T, stringsAsFactors = F)
an_file = read.table(paste0(base_path, "/input/evidence/expTFs/g296s_annotation_bgt15_2.txt"),sep = "\t", 
                     header=T, stringsAsFactors = F)
pr_file = read.table(paste0(base_path, "/input/evidence/expTFs/g296s_proteinGroups.txt"),sep = "\t", 
                     header=T, stringsAsFactors = F)

msStatsInputFiltered <- MaxQtoMSstatsFormat(ev_file, an_file, pr_file)
processedFilteredNoNorm <- dataProcess(msStatsInputFiltered, censoredInt = '0', normalization = F)
#, featureSubset = "top3")
quantificationFiltered <- quantification(processedFilteredNoNorm, type="Sample")

#Generates 'quant' for dowsntream analyses (model analyses and others)
quant <- quantificationFiltered
row.names(quant) <- quant[,1]
quant <- quant[-1]
colnames(quant) <- gsub("\\-", "_", colnames(quant))
colnames(quant) <- gsub("\\+", "", colnames(quant))
quantBackup <- quant


for (run_type in c("gained", "lost")){
  
  int_type = paste0("g296s_vs_wt_",run_type)
  if(run_type == "gained"){
      ref = "wt"
    } else{
      ref = "g296s"
  }
  
  samples <- colnames(quant)
  variant <- gsub("_.*$", "", samples)
  variant <- as.factor(variant)
  variant <- relevel(variant, ref = ref)
  replicate = c()
  idx = 1
  for(i in range(1,length(unique(variant)))){
    n = sum(variant==variant[idx])
    replicate = c(replicate, rep(1:n))
    idx = idx + n
  }
  
  replicate <- as.factor(replicate)
  phenoData <- data.frame(cbind(Samples=samples, Variant=variant, Replicate=replicate))
  # head(PhenoData)
  boxplot(quant)
  
  #Limma helps with the model functions.
  # source("https://bioconductor.org/biocLite.R")
  # biocLite("limma")
  require(limma)
  
  # The following model controls for replicate variability. 
  # It also identifies all hits significantly different from controls.
  design <- model.matrix(~0+variant+replicate)
  fit <- lmFit(quant, design=design) #Note this is with log values!
  
  #Here  we get significant factors for each one 'individually'.
  #Iterate through baits.
  combinedTable <- data.frame()
  for(condition in colnames(design)[grepl("variant", colnames(design))]){
    #if(condition=="variantcontrol"){next}
    print(condition)
    #contrastM <- makeContrasts(paste0(condition, "-variantControl"), levels=design)
    #contrastM <- makeContrasts(paste0(condition, "+variantBAG3_WT-2*variantControl"), levels=design)
    contrastM <- makeContrasts(condition, levels=design)
    contrastFit <- contrasts.fit(fit, contrastM)
    contrastFit <- eBayes(contrastFit, proportion=0.3, trend=TRUE)
    t <- topTable(contrastFit, number=200)
    t$gene <- rownames(t)
    t$condition <- condition
    combinedTable <- rbind(combinedTable, t)
  }
  
  interactorTable <- combinedTable[combinedTable$adj.P.Val<0.05 & combinedTable$logFC>=0,] 
  print(nrow(interactorTable))
  #(also filter out things where theres more on control than bait, though it seems none  of the significant 
  # ones were like that)
  
  # Convert UPIDs to genes
  quant_unimap = read.table(paste0(base_path, "/input/aliases/combined_unimap.txt"),
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  quant_unimap <- unique(quant_unimap)
  interactorTable$GeneName = quant_unimap$GeneSymbol[match(interactorTable$gene, quant_unimap$UniProt)]
  combinedTable$GeneName = quant_unimap$GeneSymbol[match(combinedTable$gene, quant_unimap$UniProt)]
  write.table(interactorTable, file = paste0(base_path, "/intermediate/interaction_scoring/regression/", 
                                             int_type, "_regression.txt"), sep = "\t", row.names = F)
}



#All genes that were significant, one by one bait. 
#I did this independently for each pairwise combination just because I had multiple conditions with different controls
#in my case I had one of mutants vs wild type and one of drugs vs DMSO, but you could see this being used for 
#different WT baits vs mutant baits.
combinedTable2 <- data.frame()
for(condition in unique(gsub("_.$", "", variant))){
  if(condition %in% c("Control", "BAG3_WT")){next}
  if(condition=="BAG3_Bort"){
    variant <- relevel(variant, ref = "BAG3_DMSO")
    quantShort <- quant[,grepl(paste0(condition, "|DMSO"), colnames(quant))]
    quantShort <- quantShort[genes,]
    variantShort <- droplevels(variant[grepl(paste0(condition, "|DMSO"), colnames(quant))])
  }else{
    variant <- relevel(variant, ref = "BAG3_WT")
    quantShort <- quant[,grepl(paste0(condition, "|WT"), colnames(quant))]
    quantShort <- quantShort[genes,]
    variantShort <- droplevels(variant[grepl(paste0(condition, "|WT"), colnames(quant))])
  }
  replicateShort <- droplevels(replicate[1:8])
  print(replicateShort)
  print(head(quantShort))
  # print(variantShort)
  #vector of BAG3 levels to use for quantitation
  BAG3 <- quantShort[row.names(quantShort)=="BAG3",]
  BAG3 <- as.matrix(BAG3, ncol=1)
  BAG3 <- as.vector(BAG3)
  design2 <- model.matrix(~variantShort+replicateShort+BAG3) #Normalize by BAG3 levels
  fit2 <- lmFit(quantShort, design=design2)
  fit2 <- eBayes(fit2, proportion=0.3)
  t <- topTable(fit2, coef=2, number=200)
  t$gene <- rownames(t)
  t$condition <- condition
  #print(head(t))
  combinedTable2 <- rbind(combinedTable2, t)
}
combinedTable2$P.Value <- combinedTable2$P.Value
combinedTable2$sig <- combinedTable2$P.Value<=0.05# & abs(combinedTable2$logFC)>=log2(2) 
#I will later use this for a volcano plot or heatmap.
sigTable2 <- combinedTable2[combinedTable2$sig,]  #Final table with all the significantly different hits. 



