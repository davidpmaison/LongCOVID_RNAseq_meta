# === Differential Expression Analysis ===
# Author: David Maison

# Load required packages
library(DESeq2)
library(data.table)

# === Step 1: Load Metadata ===
# Read metadata (sample IDs, filenames, conditions, covariates)
metadata <- fread("metadata.csv")
rownames(metadata) <- metadata$sample_id

# === Step 2: Load Count Data ===
# Collect all STAR ReadsPerGene.out.tab files
count_files <- list.files(pattern = "*ReadsPerGene.out.tab$")
count_data <- fread(count_files[1])[, c(1,4)]  # GeneID and 4th column (unique mapping)

# Merge counts from all files
for (i in 2:length(count_files)) {
  temp <- fread(count_files[i])[, 2]  # 2nd column contains unstranded counts
  count_data <- cbind(count_data, temp)
}

# Remove first 4 summary lines
count_data <- count_data[5:.N, ]

# Set gene names and adjust column names
colnames(count_data) <- c("GeneID", gsub("_paired_trimmed_HUMAN_ReadsPerGene.out.tab", "", count_files))
rownames(count_data) <- count_data$GeneID
count_data <- count_data[, -1, with = FALSE]

# Filter to samples present in metadata
count_data <- count_data[, metadata$sample_id, with = FALSE]

# Add pseudocounts (+1)
count_data <- count_data + 1

# === Step 3: Differential Expression Analysis ===
dds <- DESeqDataSetFromMatrix(countData = as.matrix(count_data),
                              colData = metadata,
                              design = ~ condition + study_id)

dds <- DESeq(dds)

# === Step 4: Save DESeq2 Object (Optional) ===
save(dds, file = "dds_object.RData")

# === Step 5: Extract Results ===
# Get the results ("condition", "Treatment", "Control")
res <- results(dds, contrast = c("condition", "LC", "CR"), alpha = 0.05)

# Save DEG results
write.csv(as.data.frame(res), file = "DEG_results_CR_vs_LC.csv", row.names = TRUE)
