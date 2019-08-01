base_path = "/Users/student/Documents/PollardLab/APMS2"

quant_unimap = read.table(paste0(base_path, "/input/aliases/combined_unimap.txt"),
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)
quant_unimap <- unique(quant_unimap)

rnaseq = read.csv(paste0(base_path,"/input/rnaseq/DifferentialExpressionResults_koStudy_Sigma.csv"))
g_blacklist = rnaseq[which(rnaseq$logFC.TreatmentGKO < -1 & rnaseq$FDR < 0.05),]
n_blacklist = rnaseq[which(rnaseq$logFC.TreatmentNKO < -1 & rnaseq$FDR < 0.05),]
t_blacklist = rnaseq[which(rnaseq$logFC.TreatmentTKO < -1 & rnaseq$FDR < 0.05),]

ints = read.csv(paste0(base_path, "/intermediate/interactome_lists/saintq_n/combined_APMS_interactome_G001_T05_N1.csv"))

f_gata = ints[which(ints$Bait == "GATA4" & ints$Prey_genename %in% g_blacklist$hgnc_symbol),]
f_black_g = g_blacklist[which(g_blacklist$hgnc_symbol %in% f_gata$Prey_genename),]

f_tbx5 = ints[which(ints$Bait == "TBX5" & ints$Prey_genename %in% t_blacklist$hgnc_symbol),]
f_black_t = t_blacklist[which(t_blacklist$hgnc_symbol %in% f_tbx5$Prey_genename),]

f_nkx2 = ints[which(ints$Bait == "NKX25" & ints$Prey_genename %in% n_blacklist$hgnc_symbol),]
f_black_n = n_blacklist[which(n_blacklist$hgnc_symbol %in% f_nkx2$Prey_genename),]

# Check to see whether these logFCs are drastically different from the intensity data
g_intensities = read.table(paste0(base_path, "/intermediate/interaction_scoring/saintq/GATA4_mergedIntensities.txt"), sep = "\t",
                           header = T, stringsAsFactors = F)
g_intensities$logFC.intensity = log2(g_intensities$Intensity_KO/g_intensities$Intensity_WT)
t_intensities = read.table(paste0(base_path, "/intermediate/interaction_scoring/saintq/TBX5_mergedIntensities.txt"), sep = "\t",
                           header = T, stringsAsFactors = F)
t_intensities$logFC.intensity = log2(t_intensities$Intensity_KO/t_intensities$Intensity_WT)
n_intensities = read.table(paste0(base_path, "/intermediate/interaction_scoring/saintq/NKX25_mergedIntensities.txt"), sep = "\t",
                           header = T, stringsAsFactors = F)
n_intensities$logFC.intensity = log2(n_intensities$Intensity_KO/n_intensities$Intensity_WT)

g = merge(f_black_g, g_intensities, by.x = "hgnc_symbol", by.y = "GeneNames", all.x=T)
t = merge(f_black_t, t_intensities, by.x = "hgnc_symbol", by.y = "GeneNames", all.x=T)

g_tab <- g[,c("hgnc_symbol","logFC.TreatmentGKO", "FDR", "logFC.intensity", "BFDR")]
t_tab <- t[,c("hgnc_symbol","logFC.TreatmentTKO", "FDR", "logFC.intensity", "BFDR")]

########## COMPARE TO REGRESSION INTENSITIES ##########
g_reg_intensities = read.table(paste0(base_path, "/intermediate/interaction_scoring/regression/GATA4_regression.txt"), sep = "\t",
                               header = T, stringsAsFactors = F)
g = merge(f_black_g, g_reg_intensities, by.x = "hgnc_symbol", by.y = "GeneName", all.x=T)


########## WRITE OUT BLACKLIST ############
# Write out all genes except the bait
g_blacklist = rnaseq[which(rnaseq$logFC.TreatmentGKO < -1 & rnaseq$FDR < 0.05 & rnaseq$hgnc_symbol != "GATA4"),]
n_blacklist = rnaseq[which(rnaseq$logFC.TreatmentNKO < -1 & rnaseq$FDR < 0.05 & rnaseq$hgnc_symbol != "NKX25"),]
t_blacklist = rnaseq[which(rnaseq$logFC.TreatmentTKO < -1 & rnaseq$FDR < 0.05 & rnaseq$hgnc_symbol != "TBX5"),]

write.table(g_blacklist$hgnc_symbol, file = paste0(base_path, "/intermediate/rnaseq/GATA4_negativeFC_blacklist.txt"),
            quote = F, row.names = F, col.names = F)
write.table(t_blacklist$hgnc_symbol, file = paste0(base_path, "/intermediate/rnaseq/TBX5_negativeFC_blacklist.txt"),
            quote = F, row.names = F, col.names = F)
write.table(n_blacklist$hgnc_symbol, file = paste0(base_path, "/intermediate/rnaseq/NKX25_negativeFC_blacklist.txt"),
            quote = F, row.names = F, col.names = F)
