# ---------------------------------------------------------------------------------------------
# Script to prepare taxonomy data for each sample from the GTDBtk table

# R 3.6.3

# Author: Eva Sanmartín Vázquez

# ---------------------------------------------------------------------------------------------

# Load packages
library('data.table')

# Load data
paths <- c(`20IB32` = "~/GTDBTK_20IB32/20IB32_gtdbtk_todo_summary.tsv",
           `23IB32` = "~/GTDBTK_23IB32/23IB32_gtdbtk_todo_summary.tsv",
           CO = "~/GTDBTK_CO/CO_gtdbtk_todo_summary.tsv")

# -----------------------------
# Editing the table to keep the classification only and separate the taxa
# -----------------------------

GTDBTK <- list()
taxonomy <- list()
separated <- list()

for (sample in names(paths)) {
   GTDBTK[[sample]] <- read.table(paths[sample], sep = "\t", header = T) # Load GTDB table
   taxonomy[[sample]] <- as.vector(GTDBTK[[sample]][, "classification"]) # Getting the classification
   separated[[sample]] <- strsplit(taxonomy[[sample]], split = "[;]") # Separating the taxa
   separated[[sample]] <- transpose(as.data.frame(separated[[sample]])) # Transposing the data frame
   rownames(separated[[sample]]) <- GTDBTK[[sample]][,"user_genome"] # Bin names to row names
   colnames(separated[[sample]]) <- c("domain", "phylum", "class", "order",
                                      "family", "genus", "species") # Taxa names to column names
}

# Saving individual sample tables
# write.table(separated$`20IB32`, "20IB32_taxonomy.tsv", sep = "\t")
# write.table(separated$`23IB32`, "23IB32_taxonomy.tsv", sep = "\t")
# write.table(separated$CO, "CO_taxonomy.tsv", sep = "\t")

# -----------------------------
# Adding the three samples together into a single data frame and adding sample names
# -----------------------------

taxonomy <- data.frame()

for (sample in names(paths)) {
    # this.sample <- read.table(paths[sample], sep = "\t")
    bin.names <- rownames(separated[[sample]])
    rownames(separated[[sample]]) <- NULL
    separated[[sample]] <- cbind(sample = rep(sample, nrow(separated[[sample]])), bin.names, separated[[sample]])
    taxonomy <- rbind(taxonomy, separated[[sample]])
}

write.table(taxonomy, "taxonomia_todo.tsv", sep = "\t")