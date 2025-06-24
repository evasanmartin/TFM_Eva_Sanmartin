# -----------------------------------------------------------------------------
# Script for counting the reads that belong to each KEGG in the sample from 
# SQM data, preparing the data for rarefaction plot and plotting the rarefaction 
# curve.

# R 3.6.3 and 4.2.1

# Author: Eva Sanmartín Vázquez and Alicia García Roldán

# -----------------------------------------------------------------------------

# R 3.6.3

library('SQMtools')

# ---------------------------------
# ----- LOAD SQM PROJECT DATA -----
# ---------------------------------

projPaths = c(`20IB32` = "~/SQM_20IB32", `23IB32` = "~/SQM_23IB32", CO = "~/SQM_CO_20_23")
projs=list()
for(sample in names(projPaths)){
	projs[[sample]] = loadSQM(projPaths[sample], engine='data.table')
}

# ---------------------------------
# ----- DIVERSITY DATA FRAME ------
# ---------------------------------

# Modified script to read SQM project data instead of tables.

KEGG <- list()

for (sample in c("20IB32", "23IB32")) {
    KEGG[[sample]] <- projs[[sample]]$functions$KEGG$abund # Table with KEGG ids and abundances (reads)
    # write.table(KEGG[[sample]], file = paste("Diversity_table_KEGG_", sample, ".tsv", sep = ""), quote = F, col.name = F, sep = "\t")
}

# -----------------------------------------------------------------------------

# R 4.2.1

library('vegan')
library('ggplot2')
library('extrafont')

# ---------------------------------------
# ---- MERGING THE SAMPLES TOGETHER -----
# ---------------------------------------

# Reading data from the saved diversity tables

alldiv <- c(`20IB32` = "Rarefaction/Diversity_table_KEGG_20IB32.tsv", 
            `23IB32` = "Rarefaction/Diversity_table_KEGG_23IB32.tsv")

div <- list()

for (sample in names(alldiv)) {
   div[[sample]] <- read.table(alldiv[[sample]], header = F, sep = "\t")
}

# Merging samples into a single data frame

merged_20_23 <- merge(div$`20IB32`, div$`23IB32`, all = T, by = 1)
merged_20_23[is.na(merged_20_23)] <- 0
rownames(merged_20_23) <- merged_20_23[,1] # Turning first column (tax) into row names
merged_20_23[,1] <- NULL # Removing first column
colnames(merged_20_23) <- c("20IB32", "23IB32") # Naming columns with samples

# For the "coassembly" sample, adding the samples 20 and 23 together

CO <- rowSums(merged_20_23)
merged_20_23 <- cbind(merged_20_23, CO)

# Save the table
# write.table(merged_20_23, "Rarefaction/Diversity_table_KEGG_all.tsv", sep = "\t", quote = F)

# ------------------------------
# ---- DATASET FOR THE PLOT ----
# ------------------------------

# Read the saved table
# merged_20_23 <- read.table("Rarefaction/Diversity_table_KEGG_all.tsv", 
#                 sep = "\t", header = T, quote = "", check.names = F)

merged_final <- t(merged_20_23) # Transpose the matrix

# Richness calculations

S <- specnumber(merged_final)
raremax <- min(rowSums(merged_final))
Srare <- rarefy(merged_final, raremax)

# Function rarecurve with tidy = T returns a data frame that can be used
# to plot with ggplot. It has 3 columns: Site (sample), Sample (sample
# size -> x axis) and Species (number of families -> y axis).

# plot.df <- rarecurve(merged_final, step = 20, sample = raremax, tidy = T)

plot.df <- rarecurve(merged_final, step = 2000, sample = raremax, tidy = T)

# --------------
# ---- PLOT ----
# --------------

ggplot(plot.df, aes(x = Sample, y = Species, colour = Site)) +
    geom_vline(xintercept = raremax, color = "#a0a0a0", linewidth = 0.2) +
    geom_hline(yintercept = S, color = "#a0a0a0", linewidth = 0.2) +
    geom_path(linewidth = 1) +
    xlab("Sample size") +
    ylab("KEGG IDs") +
    labs(color = "Sample") +
    scale_color_manual(values = c("#8cabd9", "#f6a7b8", "#dbd55f")) +
    scale_y_continuous(expand = expansion(mult = c(0,0.05))) + # Removing space below 0 in Y axis
    theme(
        text = element_text(family = "Poppins"),
        legend.position = "right",
        legend.key = element_blank(),
        legend.title = element_text(size = 17, face = "bold"),
        legend.text = element_text(size = 15),
        panel.background = element_blank(), # White background
        axis.title.x = element_text(colour="black", margin = margin(r = 10, t=15), 
                      size = 17, face = "bold"), # X axis title
        axis.title.y = element_text(colour="black", margin = margin(r = 10, t=10), 
                      size = 17, face = "bold"), # Y axis title
        axis.text.x = element_text(colour="black", margin = margin(r = 5, t=5), 
                      size= 15), # X axis numbers
        axis.text.y = element_text(colour="black", margin = margin(r = 5, t=5),  
                      size= 15), # Y axis numbers
        axis.line = element_line(size=0.75), # Axis lines
        axis.ticks = element_line(size=0.75) # Little lines from axis to numbers
    )

ggsave("Rarefaction_reads_KEGG.png", width = 25, height = 20, units = "cm")
ggsave("Rarefaction_reads_KEGG.svg", width = 25, height = 20, units = "cm")
