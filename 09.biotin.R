# -----------------------------------------------------------------------
# This script measures the abundance of genes related to biotin metabolism
# in medium and high quality bins obtained with SqueezeMeta.

# R 3.6.3

# Author: Eva Sanmartín Vázquez and Alicia García Roldán
# -----------------------------------------------------------------------

# Load packages

library('SQMtools')

# Load SQM project data

projPaths <- c(`20IB32` = "~/SQM_20IB32", `23IB32` = "~/SQM_23IB32")
projs <- list()
for(sample in names(projPaths)){
	projs[[sample]] = loadSQM(projPaths[sample], engine='data.table')
}

# Load biotin genes list (names and KO):

genes <- read.table("Genes/biotin.tsv", header = T)

# ----------------------------
# ---- PREPARING THE DATA ----
# ----------------------------

# We create two lists, one with the "real" names of the bins and another
# one with "fake" names (sample and tidy bin name).

# List with real names (bin.fa.contigs)

sample_realname <- list()

for (i in names(projPaths)){
  sample_realname[[i]] = rownames(projs[[i]]$bins$table)
}

# List with fake names (sample_bin)

sample_fakename <- list()

for (sample in names(sample_realname)) {
  realbin = sample_realname[[sample]][0:length(sample_realname[[sample]])]
  sample_fakename[[sample]] = paste(sample, 
                                    unlist(strsplit(realbin, ".fa"))[seq(1,length(unlist(strsplit(realbin, ".fa"))), by = 2)], 
                                    sep = "_")
}

# Now we filter high and medium quality bins according to completeness
# and contamination percentages and we save them in two other lists.

sample2coolrealbin <- list()
sample2coolfakebin <- list()

for (sample in names(sample_realname)){
  bins = sample_realname[[sample]]
  for(b in bins) {
    comp = projs[[sample]]$bin$table[b,'Completeness']
    cont = projs[[sample]]$bin$table[b,'Contamination']
    if(!is.na(comp) & !is.na(cont) & comp >= 50 & cont < 10) { # Selecting high and medium quality
      sample2coolrealbin[[sample]] = c(sample2coolrealbin[[sample]], b) # List of real names
      sample2coolfakebin[[sample]] = c(sample2coolfakebin[[sample]], # Lista de fake names
                                       sample_fakename[[sample]][which(sample_realname[[sample]] == b)])
    }
  }
}

# To be able to access functions data for each bin, we have to
# ~*~subsetBins~*~, turning them into SQM objects in a new list
# (which we access with the fake names).

mis_bincitos <- list()

for (sample in names(sample2coolrealbin)){
  realbins = sample2coolrealbin[[sample]]
  fakebins = sample2coolfakebin[[sample]]
  for (i in 1:length(sample2coolrealbin[[sample]])){
    mis_bincitos[[fakebins[i]]] = subsetBins(projs[[sample]], realbins[i])
  }
}

# ---------------------------------
# ---- SEARCHING FOR THE GENES ----
# ---------------------------------

# Adding together homologous genes in order to see if the reaction occurs:
# fabB + fabF y bioA (K00833|K19563) + bio3-bio1 (K19562)

# First we create the table to fill out

genes.table <- paste(genes$Gene, " (", genes$KO, ")", sep = "") # Name and KO of the genes
genes.table[2] <- paste(genes.table[2], genes.table[3], sep = " + ") # Adding fabB and fabF
genes.table <- genes.table[-3]
genes.table[8] <- paste("bioA (K00833|K19563)", genes.table[10], sep = " + ") # Adding bioA and bio3-bio1
genes.table <- genes.table[c(-9,-10)]

table <- matrix(0, nrow = length(mis_bincitos), ncol = length(genes.table), 
                dimnames = list(names(mis_bincitos), genes.table))

# Now we search for the genes in each bin

for(bin in names(mis_bincitos)){
    here <- mis_bincitos[[bin]]
    for (gene in genes$KO) {
        if (gene %in% rownames(here$functions$KEGG$copy_number) == T) { # Checking that the gene is present
            for (i in grep(gene, colnames(table))) { # Searching for the gene in the columns, loops in case there are several matches
              table[bin,i] <- table[bin,i] +
                              here$functions$KEGG$copy_number[gene,]
            }
        }
    }
}

# Saving the table for later
# write.table(table, "biotin_table.tsv", sep = "\t", quote = F)

# ----------------------------------
# ---- ADDING COMPLETENESS DATA ----
# ----------------------------------

# Getting completeness data from checkm and adding to the
# gene abundance matrix

checkm <- read.table("checkm_todo_copy.tsv", sep = "\t", header = T, stringsAsFactors = F)

table2 <- cbind(table, "Completeness" = NA)

for (sample in names(sample2coolfakebin)) {
   for (i in 1:length(sample2coolfakebin[[sample]])) {
      table2[sample2coolfakebin[[sample]][i],"Completeness"] <- checkm[checkm[,"Sample"] == sample & checkm[,"Bin"] == sample2coolrealbin[[sample]][i], "Completeness"]
   }
}

# Saving the table for later
# write.table(table2, "biotin_comp_table.tsv", sep = "\t", quote = F)
