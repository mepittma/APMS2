# Lists of variants in cases and controls
LoF_cases <- read.csv(paste0(base_path, "/input/variants/LoF_cases.csv"), skip=1, header = TRUE, stringsAsFactors = FALSE)
LoF_ctrls <- read.csv(paste0(base_path, "/input/variants/LoF_ctrls.csv"), skip=1, header = TRUE, stringsAsFactors = FALSE)

DNV_cases <- read.csv(paste0(base_path, "/input/variants/DNV_cases.csv"), skip=1, header = TRUE, stringsAsFactors = FALSE)
DNV_ctrls <- read.csv(paste0(base_path, "/input/variants/DNV_ctrls.csv"), skip=1, header = TRUE, stringsAsFactors = FALSE)

# Get control dataframe - synonymous DNVs
DNV_case_syn <- DNV_cases[which(DNV_cases$Variant.Class == "syn"),]
DNV_ctrl_syn <- DNV_ctrls[which(DNV_ctrls$Variant.Class == "syn"),]

# Remove DNV entries that weren't predicted to be damaging
DNV_cases <- DNV_cases[which(!(DNV_cases$Variant.Class == "syn")),]
DNV_ctrls <- DNV_ctrls[which(!(DNV_ctrls$Variant.Class == "syn")),]

rec_cases <- read.csv(paste0(base_path, "/input/variants/recessive_cases.csv"), skip=1, 
                      header = TRUE, stringsAsFactors = FALSE)
rec_ctrls <- read.csv(paste0(base_path, "/input/variants/recessive_ctrls.csv"), skip=1, 
                      header = TRUE, stringsAsFactors = FALSE)

combined_cases <- rbind(LoF_cases[,c('Blinded.ID','Gene')], DNV_cases[,c('Blinded.ID','Gene')])
combined_ctrls <- rbind(LoF_ctrls[,c('Blinded.ID','Gene')], DNV_ctrls[,c('Blinded.ID','Gene')])


# Replace any instance of "NKX2-5" with "NKX25"
LoF_cases$Gene[which(as.vector(LoF_cases$Gene) == "NKX2-5")] = "NKX25"
LoF_ctrls$Gene[which(as.vector(LoF_ctrls$Gene) == "NKX2-5")] = "NKX25"
DNV_cases$Gene[which(as.vector(DNV_cases$Gene) == "NKX2-5")] = "NKX25"
DNV_ctrls$Gene[which(as.vector(DNV_ctrls$Gene) == "NKX2-5")] = "NKX25"
rec_cases$Gene[which(as.vector(rec_cases$Gene) == "NKX2-5")] = "NKX25"
rec_ctrls$Gene[which(as.vector(rec_ctrls$Gene) == "NKX2-5")] = "NKX25"
combined_cases$Gene[which(as.vector(combined_cases$Gene) == "NKX2-5")] = "NKX25"
combined_ctrls$Gene[which(as.vector(combined_ctrls$Gene) == "NKX2-5")] = "NKX25"
