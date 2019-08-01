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

rnaseq = read.csv(paste0(base_path,"/input/rnaseq/DifferentialExpressionResults_koStudy_Sigma.csv"))
int_type = "GATA4"
level = "protein"
  
ev_file = paste0(base_path, "/input/evidence/expTFs/",int_type,"_evidence.txt")
out_file = paste0(base_path, "/intermediate/interaction_scoring/saintq/", int_type, "_", level, "_norm_true.txt")
  
ev = read.table(ev_file, sep = "\t", header = T, stringsAsFactors = FALSE)
ko = ev[which(grepl("KO|Cont",ev$Experiment)),c("Proteins","Intensity","Experiment")]
wt = ev[which(!grepl("KO|Cont",ev$Experiment)),c("Proteins","Intensity","Experiment")]

gene_lookup <- unique(ev[,c("Proteins","Gene.names")])

# Normalize the control dataframe
library("tidyr")
library("tibble")
ko = rowid_to_column(ko)
mko = ko %>% spread(Experiment,Intensity)
mko = aggregate(mko[,c(3:ncol(mko))], by=list(Proteins=mko$Proteins), FUN=sum, na.rm=T)
mko[mko == 0] <- NA

# Step 1: If there are rows with all NA, replace with 90% min observed in each control run
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

# Step 2: if there is at least one prey protein/peptide intensity, replace with 90% min-observed in that protein
for (i in c(1:nrow(mko))){
  min_val = 0.9 * as.numeric(min(mko[i,2:ncol(mko)][!is.na(mko[i,2:ncol(mko)])]))
  mko[i,][is.na(mko[i,])] <- min_val
}

# Step 3: normalize so that the average total intensity across all bait purifications is the same as 
# the average total intensity in the controls
mko$MeanIntensity = 0
for (i in c(1:nrow(mko))){
  mko$MeanIntensity[i] = mean(as.numeric(mko[i,c(2:(ncol(mko)-1))]))
}
ctrl_mean = mean(mko$MeanIntensity)

wt = rowid_to_column(wt)
mwt = wt %>% spread(Experiment,Intensity)
mwt = aggregate(mwt[,c(3:ncol(mwt))], by=list(Proteins=mwt$Proteins), FUN=sum, na.rm=T)
mwt[mwt == 0] <- NA

mwt$MeanIntensity = 0
for (i in c(1:nrow(mwt))){
  mwt$MeanIntensity[i] = mean(as.numeric(mwt[i,c(2:(ncol(mwt)-1))]), na.rm = T)
}
bait_mean = mean(mwt$MeanIntensity, na.rm=T)

# Coerce averages to be the same across all bait/all ctrl
factr = ctrl_mean/bait_mean
mwt$wt_NormalizedIntensity = mwt$MeanIntensity * factr
mwt$wt_NormalizedIntensity[is.na(mwt$NormalizedIntensity)] <- 0

# Compare logFC between the two
gene_lookup <- unique(ev[,c("Proteins","Gene.names")])

norm_intensities = as.data.frame(cbind(mko$Proteins, mko$MeanIntensity, mwt$NormalizedIntensity[match(mko$Proteins, mwt$Proteins)]))
names(norm_intensities) = c("Proteins","KO_Intensity","WT_Intensity")
norm_intensities$GeneName = gene_lookup$Gene.names[match(norm_intensities$Proteins, gene_lookup$Proteins)]

