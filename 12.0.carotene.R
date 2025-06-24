# -----------------------------------------------------------------------
# Script that measures the abundance of genes related to carotenoid
# (beta-carotene and bacterioruberin) methabolism in medium and high 
# quality bins from SQM.

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

# Load carotenoid methabolism genes list (names and KO):

genes <- read.table("Genes/carotene.tsv", header = T)

# ----------------------------
# ---- PREPARING THE DATA ----
# ----------------------------

# Filtering bins according to completeness and contamination and saving
# them in a list under their sample names.

sample2coolrealbin <- list()

for (sample in names(projPaths)) {
    tablabins <- projs[[sample]]$bin$table
    comp <- tablabins[,"Completeness"]
    cont <- tablabins[,"Contamination"]
    sample2coolrealbin[[sample]] <- rownames(tablabins[!is.na(comp) & !is.na(cont) & comp >= 50 & cont < 10,])
}

# Creating "fake names" and saving them in another list

sample2coolfakebin <- list()

for (sample in names(projPaths)) {
    realbin <- sample2coolrealbin[[sample]]
    sample2coolfakebin[[sample]] = paste(sample, 
                                    unlist(strsplit(realbin, ".fa"))[seq(1,length(unlist(strsplit(realbin, ".fa"))), by = 2)], 
                                    sep = "_")
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

# Some genes are homologous, so we will add their copy numbers, 
# because we just want to know if that reaction takes place:
# crtB + al-2, crtI + pds and lcyB + al-2

# First we create the table to fill out

genes.table <- c("crtB (K02291) + al-2 (K17841)", "crtI (K10027) + pds (K02293)",
                "lcyB (K06443) + al-2 (K17841)", "crtY (K22502)", "lyeJ (K20616)",
                "crtD (K20611)", "cruF (K08977)")

table <- matrix(0, nrow = length(mis_bincitos), ncol = length(genes.table), 
                dimnames = list(names(mis_bincitos), genes.table))

# Now we search for the genes in each bin

for(bin in names(mis_bincitos)){
    here <- mis_bincitos[[bin]]
    for (gene in genes$KO) {
        if (gene %in% rownames(here$functions$KEGG$copy_number) == T) { # Comprueba que el gen exista en la muestra
          for (i in grep(gene, colnames(table))) { # Busca el gen en las columnas, el loop es por si está en más de una
            table[bin,i] <- table[bin,i] +
                            here$functions$KEGG$copy_number[gene,]
          }
        }
    }
}

# for(bin in names(mis_bincitos)){
#     aqui <- mis_bincitos[[bin]]
#     for (gene in genes$KO) {
#        if (gene %in% rownames(aqui$functions$KEGG$copy_number) == T) {
#           if (gene == "K02291") {
#              tabla[bin, "crtB (K02291) + al-2 (K17841)"] <- tabla[bin, "crtB (K02291) + al-2 (K17841)"] + 
#                                                             aqui$functions$KEGG$copy_number["K02291",]
#           }
#           if (gene == "K17841") {
#              tabla[bin, "crtB (K02291) + al-2 (K17841)"] <- tabla[bin, "crtB (K02291) + al-2 (K17841)"] +
#                                                             aqui$functions$KEGG$copy_number["K17841",]
#              tabla[bin, "lcyB (K06443) + al-2 (K17841)"] <- tabla[bin, "lcyB (K06443) + al-2 (K17841)"] +
#                                                             aqui$functions$KEGG$copy_number["K17841",]
#           }
#           if (gene == "K10027") {
#              tabla[bin, "crtI (K10027) + pds (K02293)"] <- tabla[bin, "crtI (K10027) + pds (K02293)"] +
#                                                            aqui$functions$KEGG$copy_number["K10027",]
#           }
#           if (gene == "K02293") {
#              tabla[bin, "crtI (K10027) + pds (K02293)"] <- tabla[bin, "crtI (K10027) + pds (K02293)"] +
#                                                            aqui$functions$KEGG$copy_number["K02293",]
#           }
#           if (gene == "K06443") {
#              tabla[bin, "lcyB (K06443) + al-2 (K17841)"] <- tabla[bin, "lcyB (K06443) + al-2 (K17841)"] +
#                                                             aqui$functions$KEGG$copy_number["K06443",]
#           }
#           if (gene %in% c("K22502", "K20616", "K20611", "K08977") == T) {
#              tabla[bin,grep(gene, colnames(tabla))] <- aqui$functions$KEGG$copy_number[gene,]
#           }
#        }
#     }    
# }

# Saving the table for later
# write.table(table, "carotene_table.tsv", sep = "\t", quote = F)

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
# write.table(table2, "carotene_comp_table.tsv", sep = "\t", quote = F)


