base_path = "/Users/student/Documents/PollardLab/APMS2"
out_path = paste0(base_path, "/output/pathogenicity_scoring")

source(paste0(base_path, "/code/load_data/load_variants.R"))

# Combine variants into a single object
cols <- c("Blinded.ID", "Cardiac.Category", "EM", "NDD", "CHR", "POS", "REF", "ALT", "Gene","pLI.Score")
mut_table <- rbind(DNV_cases[,cols], LoF_cases[,cols])

######################## Step 1: Select genes in interactomes with LoF and DNV ######################## 

norm_ints = read.table(paste0(base_path, "/input/precomp_interactomes/saintq_n_interactomes.csv"),
            sep = ",", header = TRUE, stringsAsFactors = FALSE)
og_ints = read.table(paste0(base_path, "/input/precomp_interactomes/interactomes.tsv"), 
                   sep = "\t",fill = TRUE, header = TRUE, stringsAsFactors = FALSE)

int_genes = unique(c(norm_ints$Prey_geneName, og_ints$Prey_genename))
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
  gene_data$STOP[i] = max(c(nchar(gene_data$REF[i]), nchar(gene_data$ALT[i]))) + as.numeric(gene_data$POS[i]) - 1
}

# Merge by Chr/Start/Stop data 
mg = merge(gene_data, phylop, by = c("CHR", "POS"), all.x=TRUE)
median_phylop = median(mg$phyloP[!is.na(mg$phyloP)])
mg$phyloP[is.na(mg$phyloP)] <- median_phylop

# Normalize between 0 and 1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
mg$phyloP = range01(as.numeric(mg$phyloP))

#Add to gene_data
gene_data$phyloP <- mg$phyloP

######################## Step 4: Downweight individuals with NDD ######################## 

# 0.5 penalty for individuals with NDD
gene_data$NDD_score = 1
gene_data$NDD_score[which(gene_data$NDD == "Yes")] = 0.5

######################## Step 5: Cell-specific expression score ######################## 

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

######################## Step 6: Combine scores and sorts ########################

# Replace missing pLI with median)
median_pLI = median(as.numeric(gene_data$pLI.Score[!is.na(gene_data$pLI.Score)]))
gene_data$pLI.Score[is.na(gene_data$pLI.Score)] <- median_pLI

# Add pseudocount to pLI
gene_data$pLI.Score[which(as.numeric(gene_data$pLI.Score) == 0)] <- 0.001

# Reduce dataframe
numerical_data = as.data.frame(
  gene_data[,c("pLI.Score","norm_mutperkb","NDD_score","specificity_score","phyloP")])
df2 <- mutate_all(numerical_data, function(x) as.numeric(as.character(x)))
df2$mult_score = apply(df2, 1, prod)
df2$add_score = rowSums(df2[,c(1:5)])

final <- cbind(gene_data[,c("Gene","Blinded.ID","Cardiac.Category","CHR","POS", "REF", "ALT")],
               df2)
final_mult <- final[order(-final$mult_score),]
final_add <- final[order(-final$add_score),]

write.csv(final_mult, paste0(base_path, "/output/pathogenicity_scoring/median_permuted_scores.csv"),
          quote = FALSE, row.names = FALSE)

######################## GATA4 only ########################

# All of these mutations are in the GATA4 interactome...  
GATA_genes = unique(c(norm_ints$Prey_genename[which(norm_ints$Bait == "GATA4")],
                      og_ints$Prey_genename[which(og_ints$Bait == "GATA4")]))
GATAs = final_mult[which(final_mult$Gene %in% GATA_genes),]

write.csv(GATAs, paste0(base_path, "/output/pathogenicity_scoring/GATA4_median_permuted_scores.csv"),
          quote = FALSE, row.names = FALSE)

######################## Nearby variants? ########################
res = c()
for(chr in c(1:22, "X")){
  var_list = final_mult[which(final_mult$CHR == chr),]
  if(nrow(var_list) > 1){
    
    # Sort numerically
    coords = sort(as.numeric(var_list$POS))
    genes = var_list[match(coords, var_list$POS),]
    
    # Check whether any mutations next to each other are within 800bp in either direction
    for (i in c(1:(length(coords)-1))){
      if(coords[i] - coords[i+1] <= 800 && coords[i] - coords[i+1] >= -800){
        success= paste0("Nearby variants! chr", chr, ":", coords[i], ",", coords[i+1], " in gene ", genes[i,"Gene"])
        print(success)
        res = c(res, genes[i,"Gene"])
      }
    }
  }
}

domain = final_mult[which(final_mult$Gene %in% res),]
domain = domain[order(domain$Gene, domain$mult_score),]
write.csv(domain, file = paste0(base_path, "/output/pathogenicity_scoring/nearby_variants.csv"),
            row.names = FALSE, quote = FALSE)
