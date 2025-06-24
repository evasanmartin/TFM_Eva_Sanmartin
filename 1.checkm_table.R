# ------------------------------------------------------------------------------------
# Script to create the table with quality (completeness and contamination) data from 
# medium and high quality bins (checkm_todo.tsv).

# R 3.6.3

# Author: Eva Sanmartín Vázquez

# ------------------------------------------------------------------------------------


# Loading packages
library('SQMtools')

# Load SQM project data into R
projPaths = c(`20IB32` = "~/SQM_20IB32", `23IB32` = "~/SQM_23IB32", CO = "~/SQM_CO_20_23")
projs=list()
for(sample in names(projPaths)){
	projs[[sample]] = loadSQM(projPaths[sample], engine='data.table')
}
Merlin = combineSQMlite(projs)

# --------------------

# Creating a list with a data frame for each sample that contains the names of
# the bins and completeness and contamination percentages.

checkm <- list()

for (sample in names(projPaths)) {
   checkm[[sample]] <- projs[[sample]]$bins$table[, c("Completeness", "Contamination")]
}

# Filtering according to medium quality thresholds (50% completeness and 10% contamination)

for (sample in names(checkm)) {
    goodbins <- checkm[[sample]][,"Completeness"]>50 & checkm[[sample]][,"Contamination"]<10 
    checkm[[sample]] <- checkm[[sample]][goodbins,]
    checkm[[sample]] <- checkm[[sample]][!is.na(checkm[[sample]][,"Contamination"]),] # Remove NAs
}

# --------------------

# Merging samples together into a single data frame

checkm_todo <- data.frame()

for (sample in names(checkm)){
    checkm_todo <- rbind(checkm_todo, as.data.frame(checkm[[sample]]))
}

# --------------------

# Adding a column with sample names

# Bin names (row names) to a new column
checkm_todo <- cbind(rownames(checkm_todo), data.frame(checkm_todo, row.names = NULL)) 
colnames(checkm_todo)[1] <- "Bin"

# Creating a vector with sample names
sample.names <- c()
for (sample in names(checkm)) {
   sample.names <- c(sample.names, rep(sample, nrow(checkm[[sample]])))
}

# Adding the vector as the first column
checkm_todo <- cbind(sample.names, checkm_todo) 
colnames(checkm_todo)[1] <- "Sample"

# --------------------

# Adding a column indicating if the bins are medium or high quality

Quality <- c()

for (i in rownames(checkm_todo)) {
   if ((checkm_todo[i,"Completeness"]>90&checkm_todo[i,"Contamination"]<5) == T) {
      Quality <- c(Quality, "HQ")
   } else {
      if (((checkm_todo[i,"Completeness"]>90&checkm_todo[i,"Contamination"]>5)| # Medium quality due to contamination
          (checkm_todo[i,"Completeness"]<90&checkm_todo[i,"Contamination"]<5)| # Medium quality due to completeness
          (checkm_todo[i,"Completeness"]<90&checkm_todo[i,"Contamination"]>5)) & # "Classic" medium quality
          (checkm_todo[i,"Completeness"]>50&checkm_todo[i,"Contamination"]<10)) {  # Ensuring every bin is >50% completeness and <10% contamination
         Quality <- c(Quality, "MQ")
      } else {
         Quality <- c(Quality, "LQ")
      } 
   }
}

checkm_todo <- cbind(checkm_todo, Quality)

# Saving the table
# write.table(checkm_todo, file = "checkm_todo.tsv", sep = "\t")


