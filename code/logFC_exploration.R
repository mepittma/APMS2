# This file attempts to replicate the various normalizations performed by saint_q on protein/peptide intensity data.

# Missing data imputation:
  # control purifications only
  # if all control prey protein/peptide values are missing, set intensity to 90% of the lowest-observed across 
#all proteins in each control run
  # if there is at least one prey protein/peptide intensity, the missing intensities are replaced by 90% of the 
#lowest-observed in that protein
  # average total intensity across all bait purifications is the same as the average total intensity in the controls

base_path = "/Users/student/Documents/PollardLab/APMS2"

quant_unimap = read.table(paste0(base_path, "/input/aliases/combined_unimap.txt"),
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)
quant_unimap <- unique(quant_unimap)

#### TEST CASE: creating the GATA4 blacklist - comparing intensity and RNA-seq fold-changes
library("tidyr")
library("tibble")

rnaseq = read.csv(paste0(base_path,"/input/rnaseq/DifferentialExpressionResults_koStudy_Sigma.csv"))
int_type = "GATA4"
level = "protein"
  
# Step 1: Deconvolute protein complexes
ev_file = paste0(base_path, "/input/evidence/expTFs/",int_type,"_evidence.txt")
  
ev = read.table(ev_file, sep = "\t", header = T, stringsAsFactors = FALSE)
s <- strsplit(ev$Proteins, split=";")
df_ <- data.frame(Proteins = unlist(s), Intensity = rep(ev$Intensity, sapply(s, length)), 
                  Experiment = rep(ev$Experiment, sapply(s, length)), stringsAsFactors = F)
mut_df = rowid_to_column(df_)
mdf = mut_df %>% spread(Experiment,Intensity)
mdf = aggregate(mdf[,c(3:ncol(mdf))], by=list(Proteins=mdf$Proteins), FUN=sum, na.rm=T)
mdf[mdf == 0] <- NA

mko = mdf[,c("Proteins",names(mdf)[which(grepl("KO|Cont", names(mdf)))])]
mwt = mdf[,c(names(mdf)[which(!grepl("KO|Cont", names(mdf)))])]

gene_lookup <- unique(ev[,c("Proteins","Gene.names")])



# Step 2: If there are rows with all NA, replace with 90% min observed in each control run
min_list = c()
for (i in c(2:ncol(mko))){
  min_list = c(min_list, min(mko[,i][!is.na(mko[,i])]))
}
min_list = 0.9 * min_list

# This is the hackiest shit I've ever written
for (i in c(1:nrow(mko))){
  if (is.na(mko[i,2])){
    if (is.na(mko[i,3])){
      if (is.na(mko[i,4])){
        if (is.na(mko[i, ncol(mko)])){
          mko[i,] = c(mko[i,1],min_list)
        }
      }
    }
  }
}
mko[which(mko$Proteins == "P43694"),]

# Step 3: if there is at least one prey protein/peptide intensity, replace with 90% min-observed in that protein
for (i in c(1:nrow(mko))){
  min_val = 0.9 * as.numeric(min(mko[i,2:ncol(mko)][!is.na(mko[i,2:ncol(mko)])]))
  mko[i,][is.na(mko[i,])] <- min_val
}

# Step 4: normalize so that the average total intensity across all bait purifications is the same as 
# the average total intensity in the controls
mko$MeanIntensity = 0
for (i in c(1:nrow(mko))){
  mko$MeanIntensity[i] = mean(as.numeric(mko[i,c(2:(ncol(mko)-1))]))
}
ctrl_mean = mean(mko$MeanIntensity)

mwt$MeanIntensity = 0
for (i in c(1:nrow(mwt))){
  mwt$MeanIntensity[i] = mean(as.numeric(mwt[i,c(2:(ncol(mwt)-1))]), na.rm = T)
}
bait_mean = mean(mwt$MeanIntensity, na.rm=T)

# Coerce averages to be the same across all bait/all ctrl
factr = ctrl_mean/bait_mean #ctrl measurements are about half of the wt bait measurements
mwt$wt_NormalizedIntensity = mwt$MeanIntensity * factr
mwt$wt_NormalizedIntensity[is.na(mwt$wt_NormalizedIntensity)] <- 1 #add a pseudocount of 1

# Compare logFC between the two
gene_lookup <- unique(ev[,c("Proteins","Gene.names")])

norm_intensities = as.data.frame(cbind(mko$Proteins, mko$MeanIntensity, mwt$wt_NormalizedIntensity[match(mko$Proteins, mwt$Proteins)]), stringsAsFactors = F)
names(norm_intensities) = c("Proteins","KO_Intensity","WT_Intensity")
norm_intensities$GeneName = gene_lookup$Gene.names[match(norm_intensities$Proteins, gene_lookup$Proteins)]
norm_intensities$logFC_intensity = log2(as.numeric(as.character(norm_intensities$KO_Intensity))/as.numeric(as.character(norm_intensities$WT_Intensity)))

norm_intensities$logFC_rnaseq = rnaseq$logFC.TreatmentGKO[match(norm_intensities$GeneName, rnaseq$hgnc_symbol)]
norm_intensities$FDR_rnaseq = rnaseq$FDR[match(norm_intensities$GeneName, rnaseq$hgnc_symbol)]

# # # # # # # 
# Compare for blacklist items

g_blacklist = rnaseq[which(rnaseq$logFC.TreatmentGKO < -1 & rnaseq$FDR < 0.05),]
n_blacklist = rnaseq[which(rnaseq$logFC.TreatmentNKO < -1 & rnaseq$FDR < 0.05),]
t_blacklist = rnaseq[which(rnaseq$logFC.TreatmentTKO < -1 & rnaseq$FDR < 0.05),]

ints = read.csv(paste0(base_path, "/intermediate/interactome_lists/saintq_n/combined_APMS_interactome_G001_T05_N1.csv"))

f_gata = ints[which(ints$Bait == "GATA4" & ints$Prey_genename %in% g_blacklist$hgnc_symbol),]
f_black_g = g_blacklist[which(g_blacklist$hgnc_symbol %in% f_gata$Prey_genename),]

norm_intensities[which(norm_intensities$GeneName %in% f_black_g$hgnc_symbol),]

# Load in the results from saintq analysis
gata4 = norm_intensities
out_file = paste0(base_path, "/intermediate/interaction_scoring/saintq/", 
                             int_type, "_", level, "_norm_true.txt")
saintq = read.table(out_file, stringsAsFactors=F, comment.char="", header=T)
gata4$saintq_FDR = saintq$BFDR[match(gata4$Proteins, saintq$Prey)]
sus = gata4[which(gata4$saintq_FDR < 0.001 & gata4$logFC_intensity>-1),]

#### REPEAT: creating the TBX5 blacklist - comparing intensity and RNA-seq fold-changes

int_type = "TBX5"
level = "peptide"

# Step 1: Deconvolute protein complexes
ev_file = paste0(base_path, "/input/evidence/expTFs/",int_type,"_evidence.txt")

ev = read.table(ev_file, sep = "\t", header = T, stringsAsFactors = FALSE)
s <- strsplit(ev$Proteins, split=";")
df_ <- data.frame(Proteins = unlist(s), Intensity = rep(ev$Intensity, sapply(s, length)), 
                  Experiment = rep(ev$Experiment, sapply(s, length)), stringsAsFactors = F)
mut_df = rowid_to_column(df_)
mdf = mut_df %>% spread(Experiment,Intensity)
mdf = aggregate(mdf[,c(3:ncol(mdf))], by=list(Proteins=mdf$Proteins), FUN=sum, na.rm=T)
mdf[mdf == 0] <- NA

mko = mdf[,c("Proteins",names(mdf)[which(grepl("KO|Cont", names(mdf)))])]
mwt = mdf[,c(names(mdf)[which(!grepl("KO|Cont", names(mdf)))])]

gene_lookup <- unique(ev[,c("Proteins","Gene.names")])



# Step 2: If there are rows with all NA, replace with 90% min observed in each control run
min_list = c()
for (i in c(2:ncol(mko))){
  min_list = c(min_list, min(mko[,i][!is.na(mko[,i])]))
}
min_list = 0.9 * min_list

# This is the hackiest shit I've ever written
for (i in c(1:nrow(mko))){
  if (is.na(mko[i,2])){
    if (is.na(mko[i,3])){
      if (is.na(mko[i,4])){
        if (is.na(mko[i, ncol(mko)])){
          mko[i,] = c(mko[i,1],min_list)
        }
      }
    }
  }
}
mko[which(mko$Proteins == "P43694"),]

# Step 3: if there is at least one prey protein/peptide intensity, replace with 90% min-observed in that protein
for (i in c(1:nrow(mko))){
  min_val = 0.9 * as.numeric(min(mko[i,2:ncol(mko)][!is.na(mko[i,2:ncol(mko)])]))
  mko[i,][is.na(mko[i,])] <- min_val
}

# Step 4: normalize so that the average total intensity across all bait purifications is the same as 
# the average total intensity in the controls
mko$MeanIntensity = 0
for (i in c(1:nrow(mko))){
  mko$MeanIntensity[i] = mean(as.numeric(mko[i,c(2:(ncol(mko)-1))]))
}
ctrl_mean = mean(mko$MeanIntensity)

mwt$MeanIntensity = 0
for (i in c(1:nrow(mwt))){
  mwt$MeanIntensity[i] = mean(as.numeric(mwt[i,c(2:(ncol(mwt)-1))]), na.rm = T)
}
bait_mean = mean(mwt$MeanIntensity, na.rm=T)

# Coerce averages to be the same across all bait/all ctrl
factr = ctrl_mean/bait_mean #ctrl measurements are about half of the wt bait measurements
mwt$wt_NormalizedIntensity = mwt$MeanIntensity * factr
mwt$wt_NormalizedIntensity[is.na(mwt$wt_NormalizedIntensity)] <- 1 #add a pseudocount of 1

# Compare logFC between the two
gene_lookup <- unique(ev[,c("Proteins","Gene.names")])

norm_intensities = as.data.frame(cbind(mko$Proteins, mko$MeanIntensity, mwt$wt_NormalizedIntensity[match(mko$Proteins, mwt$Proteins)]), stringsAsFactors = F)
names(norm_intensities) = c("Proteins","KO_Intensity","WT_Intensity")
norm_intensities$GeneName = gene_lookup$Gene.names[match(norm_intensities$Proteins, gene_lookup$Proteins)]
norm_intensities$logFC_intensity = log2(as.numeric(as.character(norm_intensities$KO_Intensity))/as.numeric(as.character(norm_intensities$WT_Intensity)))

norm_intensities$logFC_rnaseq = rnaseq$logFC.TreatmentGKO[match(norm_intensities$GeneName, rnaseq$hgnc_symbol)]
norm_intensities$FDR_rnaseq = rnaseq$FDR[match(norm_intensities$GeneName, rnaseq$hgnc_symbol)]

# Load in the results from saintq analysis
tbx5 = norm_intensities
out_file = paste0(base_path, "/intermediate/interaction_scoring/saintq/", 
                  int_type, "_", level, "_norm_true.txt")
saintq = read.table(out_file, stringsAsFactors=F, comment.char="", header=T)
ints = saintq[which(saintq$Bait == "TBX5" & saintq$BFDR <= 0.05),]
tbx5$saintq_FDR = saintq$BFDR[match(tbx5$Proteins, saintq$Prey)]

# Write out members that are for some reason saintq-significant, but have a positive logFC
sus = tbx5[which(tbx5$saintq_FDR < 0.05 & tbx5$logFC_intensity>-1),]


#### Creating volcano plots #### 

# Make a basic volcano plot
vulcan <- function(res){
  
  # Pseudo-count to prevent infinity p-values
  res$saintq_FDR[which(res$saintq_FDR == 0)] <- 0.0001
  
  with(res, plot(logFC_intensity, -log10(saintq_FDR), pch=20, main="LogFC (KO/WT) Intensity Volcano Plot",xlim=c(-13,6)))
  
  # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
  with(subset(res, saintq_FDR<0.05 ), points(logFC_intensity, -log10(saintq_FDR), pch=20, col="red"))
  with(subset(res, abs(logFC_intensity)>1), points(logFC_intensity, -log10(saintq_FDR), pch=20, col="orange"))
  with(subset(res, saintq_FDR<0.05 & abs(logFC_intensity)>1), points(logFC_intensity, -log10(saintq_FDR), pch=20, col="green"))
  
  # Label points with the textxy function from the calibrate plot
  library(calibrate)
  with(subset(res, saintq_FDR<.05 & abs(logFC_intensity)<1), 
       textxy(logFC_intensity, -log10(saintq_FDR), labs=GeneName, cex=.8))
}




