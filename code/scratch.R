NK = read.table(paste0(base_path, "/intermediate/interaction_scoring/saintq/NKX25_peptide_norm_true.txt"),
           sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
NK = NK[which(NK$V6 <= 0.1),]
names(NK) = c("Bait", "Prey_proteinname", "Rep", "Count", "Score", "BFDR")
unimap <- read.table(paste0(base_path, "/input/aliases/combined_unimap.txt"),
                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)
NK$GeneSymbol = unimap$GeneSymbol[match(NK$Prey_proteinname, unimap$UniProt)]
write.table(NK, paste0(base_path, "/NKX25_cutoff1.csv"), row.names = FALSE, quote = FALSE, sep = ",")

########################################
# check if there are any outliers in NKX2-5 controls
library(ggfortify)
quant = read.table(file = paste0(base_path,"/input/evidence/expTFs/saintq_inputs/NKX25/saintq_input_peptides.txt"),
                                    sep = "\t", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, comment.char = "",
                   skip = 2)
colMeans(quant[,c(3:8)])
library(vioplot)
vioplot(quant$control.1, quant$control.2, quant$control.3)
boxplot(quant$control.1, quant$control.2, quant$control.3)
autoplot(prcomp(quant[,c(3:8)]), data = quant, color = names(quant)[c(3:8)])

# Run saintq using only 2 as controls
source(paste0(base_path, "/code/functions/score_interactions.R"))
int_type = "NKX25_2ctrl"

# Load saintq norm that only used 2 controls
ctrl2 = read.table(paste0(base_path, "/input/evidence/expTFs/saintq_inputs/NKX25_c2/scores_list__saintq_peptides_test.txt__.tsv"),
                   sep = "\t", comment.char = "", header = TRUE, stringsAsFactors = FALSE)

# Load saintqnorm with 3 controls
ctrl3 = read.table(file = paste0(base_path,"/intermediate/interaction_scoring/saintq/NKX25_peptide_norm_true.txt"),
                          sep = "\t", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, comment.char = "")

ctrl2 <- ctrl2[which(ctrl2$BFDR < 0.1),]
ctrl3 <- ctrl3[which(ctrl3$BFDR < 0.1),]

