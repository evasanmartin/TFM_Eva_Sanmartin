# -------------------------------------------------------------------
# This scripts creates the tables with isoelectric point data from 
# archaea and bacteria in our metagenomes.

# R 3.6.3

# Author: Eva Sanmartín Vázquez and Alicia García Roldán
# -------------------------------------------------------------------

# ---- LOAD PACKAGES ----

library('stringr')

# ---- LOAD SQM PROJECT DATA ----
# ----
# We will use SQM tables instead of working directly with SMQtools.

# We grab coding sequences (CDS) from SQM table 13, using the
# following commands in BASH:

# grep "CDS" 13.SQM_20IB32.orftable > 20IB32_file1.txt
# grep "CDS" 13.SQM_23IB32.orftable > 23IB32_file1.txt
# grep "CDS" 13.SQM_20IB32.orftable | cut -f1 > 20IB32_file2.txt
# grep "CDS" 13.SQM_23IB32.orftable | cut -f1 > 23IB32_file2.txt

# Then we get the isoelectric points of the CDS in file 2 from
# the isoelectric points table with the following command:

# grep -f 20IB32_file2.txt PI/03.SQM_20IB32.faa.fasta2pi.txt > 20IB32_file3.txt
# grep -f 23IB32_file2.txt PI/03.SQM_23IB32.faa.fasta2pi.txt > 23IB32_file3.txt

# Now we have the tables with SQM data (file 1) and isoelectric points (file 3)
# for all coding sequences. We will read those files and extract info:
# - From SQM data: orf name, type (CDS), taxonomy
# - From PI table: orf name (check they are the same), PI

file1paths <- c(`20IB32` = "~/PI/20IB32_file1.txt", `23IB32` = "~/PI/23IB32_file1.txt")
file3paths <- c(`20IB32` = "~/PI/20IB32_file3.txt", `23IB32` = "~/PI/23IB32_file3.txt")

data_SQM <- list()
data_PI <- list()
data_all <- list()

for (sample in names(file1paths)) {
    # Reading SQM data
    data_SQM[[sample]] <- read.table(file1paths[[sample]], sep = "\t", header = F, quote = "")
    # Dividing column 9 (taxonomy) and splitting it into domain and the rest of levels.
    taxonomia <- str_split_fixed(data_SQM[[sample]]$V9, ";", n = 2)
    samples <- rep(sample, nrow(data_SQM[[sample]])) # Preparing a vector with sample names
    data_SQM[[sample]] <- cbind(samples, data_SQM[[sample]][,c(1,3)], taxonomia) # Grabbing columns of interest
    # Reading PI data
    data_PI[[sample]] <- read.table(file3paths[[sample]], sep = "\t", header = F, quote = "")
    data_PI[[sample]] <- data_PI[[sample]][,1:2]
    # Merging SQM and PI data
    data_all[[sample]] <- cbind(data_SQM[[sample]], data_PI[[sample]])
    colnames(data_all[[sample]]) <- c("Sample", "Orf_name","Type", "Domain", "Tax", "Orf_name2", "PI")
}

# Finally we merge the samples into a single data frame

df_allsamples <- data.frame()
for (sample in names(data_SQM)) {
    df_allsamples <- rbind(df_allsamples, data_all[[sample]])
}

# We only want archaeal and bacterial data, so we filter by domain column
# and save to another data frame

Arch_Bact <- subset(df_allsamples, df_allsamples$Domain == "k_Archaea" | df_allsamples$Domain == "k_Bacteria")

# Saving the tables with IP for every protein and only archaeal and bacterial proteins

# write.table(df_allsamples, file = "Isoelectric_point_total.txt",sep = "\t")
# write.table(Arch_Bact, file = "Isoelectric_point_ArchBact.txt",sep = "\t")
