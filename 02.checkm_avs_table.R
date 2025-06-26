# -----------------------------------------------------------------------------
# Script that calculates average and median completeness and contamination for:
# - high and medium quality bins (>50/<10)
# - high quality bins (>90/,<5)
# - medium quality bins (>50/<10 - >90/<5)
# from the data frame with quality data (checkm_todo)

# R 4.2.1

# Author: Eva Sanmartín Vázquez

# -----------------------------------------------------------------------------

# Load data: checkm table and sample names
checkm_todo <- read.table("checkm_todo.tsv", sep = "\t")
sample.names <- c("20IB32", "23IB32", "CO")

# ------------
# CÁLCULOS
# ------------

# -----------------------------------------------------------------------------
# HIGH AND MEDIUM QUALITY BINS

# Average and median completeness and contamination in >50/<10 bins from each sample

Complet.av <- c()
Contam.av <- c()
Complet.Q2 <- c()
Contam.Q2 <- c()

for (sample in sample.names) {
     Complet.av <- c(Complet.av, mean(checkm_todo[checkm_todo$Sample==sample,"Completeness"]))
     Contam.av <- c(Contam.av, mean(checkm_todo[checkm_todo$Sample==sample,"Contamination"]))
     Complet.Q2 <- c(Complet.Q2, median(checkm_todo[checkm_todo$Sample==sample,"Completeness"]))
     Contam.Q2 <- c(Contam.Q2, median(checkm_todo[checkm_todo$Sample==sample,"Contamination"]))
}

table <- matrix(c(Complet.av, Complet.Q2, Contam.av, Contam.Q2), nrow = 4, ncol = 3, byrow = T,
                dimnames = list(c("Complet.av", "Complet.Q2", "Contam.av", "Contam.Q2"), sample.names))

# ---------------------------------------------------------------------------------------
# MEDIUM QUALITY BINS

# Selecting medium quality bins

checkm_MQ <- checkm_todo[checkm_todo$Quality=="MQ",]

# Average and median completeness and contamination in MQ bins from each sample

Complet.av.MQ <- c()
Contam.av.MQ <- c()
Complet.Q2.MQ <- c()
Contam.Q2.MQ <- c()

for (sample in sample.names) {
     Complet.av.MQ <- c(Complet.av.MQ, mean(checkm_MQ[checkm_MQ$Sample==sample,"Completeness"]))
     Contam.av.MQ <- c(Contam.av.MQ, mean(checkm_MQ[checkm_MQ$Sample==sample,"Contamination"]))
     Complet.Q2.MQ <- c(Complet.Q2.MQ, median(checkm_MQ[checkm_MQ$Sample==sample,"Completeness"]))
     Contam.Q2.MQ <- c(Contam.Q2.MQ, median(checkm_MQ[checkm_MQ$Sample==sample,"Contamination"]))
}

table <- rbind(table, Complet.av.MQ, Complet.Q2.MQ, Contam.av.MQ, Contam.Q2.MQ)

# ---------------------------------------------------------------------------------
# HIGH QUALITY BINS

# Selecting high quality bins

checkm_HQ <- checkm_todo[checkm_todo$Quality=="HQ",]

# Average and median completeness and contamination in HQ bins from each sample

Complet.av.HQ <- c()
Contam.av.HQ <- c()
Complet.Q2.HQ <- c()
Contam.Q2.HQ <- c()

for (sample in sample.names) {
     Complet.av.HQ <- c(Complet.av.HQ, mean(checkm_HQ[checkm_HQ$Sample==sample,"Completeness"]))
     Contam.av.HQ <- c(Contam.av.HQ, mean(checkm_HQ[checkm_HQ$Sample==sample,"Contamination"]))
     Complet.Q2.HQ <- c(Complet.Q2.HQ, median(checkm_HQ[checkm_HQ$Sample==sample,"Completeness"]))
     Contam.Q2.HQ <- c(Contam.Q2.HQ, median(checkm_HQ[checkm_HQ$Sample==sample,"Contamination"]))
}

table <- rbind(table, Complet.av.HQ, Complet.Q2.HQ, Contam.av.HQ, Contam.Q2.HQ)

# ---------------------------------------------------------------------------------

write.table(table, "checkm_tabla_medias.tsv", sep = "\t", quote = F)

