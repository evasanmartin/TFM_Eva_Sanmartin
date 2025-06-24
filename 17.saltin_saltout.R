# -----------------------------------------------------------------------
# This script measures the abundance of genes related to osmoregulation
# strategies (salt-in / salt-out) and creates a heatmap.

# Author: Eva Sanmartín Vázquez and Alicia García Roldán

# R 3.6.3
# -----------------------------------------------------------------------

# Load packages
library('SQMtools')
library('tidyr')
library('tibble')
library('ggplot2')
library('extrafont')
library('paletteer')

# Load SQM project data
projPaths <- c(`20IB32` = "~/SQM_20IB32", `23IB32` = "~/SQM_23IB32")
projs <- list()
for(sample in names(projPaths)){
	projs[[sample]] = loadSQM(projPaths[sample], engine='data.table')
}

# Load list of osmoregulation genes
genes <- read.table("Genes/new_osmorregulation.tsv", sep = "\t", header = T)

# --------------------------
# PREPARING THE DATA
# --------------------------

# Filtering bins by completeness and contamination and saving them
# in a list under the corresponding sample name.

sample2coolrealbin <- list()

for (sample in names(projPaths)) {
    bintable <- projs[[sample]]$bin$table
    comp <- bintable[,"Completeness"]
    cont <- bintable[,"Contamination"]
    sample2coolrealbin[[sample]] <- rownames(bintable[!is.na(comp) & !is.na(cont) & comp >= 50 & cont < 10,])
}

# Creating "fake names" for the bins and saving them in another list

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

# --------------------------
# SEARCHING FOR THE GENES
# --------------------------

# First, we create the table that will be filled out

genes.table <- paste(genes$Name, " (", genes$KO, ")", sep = "") # Name and KEGG orthology code

table <- matrix(0, nrow = length(mis_bincitos), ncol = length(genes.table), 
                dimnames = list(names(mis_bincitos), genes.table))

# Now we search for each gene in each bin

for(bin in names(mis_bincitos)){
    here <- mis_bincitos[[bin]]
    for (gene in genes$KO) {
        if (gene %in% rownames(here$functions$KEGG$copy_number) == T) { # Checks that the gene is present in the bin
            for (i in grep(gene, colnames(table))) { # Searchs for the gene in the columns, loops in case it's in more than one
              table[bin,i] <- table[bin,i] +
                              here$functions$KEGG$copy_number[gene,]
            }
        }
    }
}

# Saving the table for later
# write.table(table, "osmorregulation_table.tsv", sep = "\t")
# Reading the saved table
# table <- read.table("osmorregulation_table.tsv", sep = "\t", check.names = F)

# Finally, we transform the table to long format for the plot

df.long <- as.data.frame(table) %>% 
           rownames_to_column(var = "Bin") %>% 
           gather(KO, Copy_number, 2:(length(colnames(table))+1))


# --------------
# ---- PLOT ----
# --------------

# Coloring bin names according to if they are bacteria or archaea:

taxonomy <- read.table("taxonomia_todo.tsv", sep = "\t", header = T) # Reading taxonomy table
taxonomy <- taxonomy[taxonomy[,"sample"] != "CO",] # Get samples 20IB32 and 23IB32
colors <- ifelse(as.vector(taxonomy[,"domain"]) == "d__Archaea", "#398b74", "#9057c6") # Colors for bacteria and archaea
color.table <- cbind(taxonomy[,c("sample", "bin.names", "domain")], colors) # Creating another table that includes colors
color.table$bin.names <- paste(color.table$sample, color.table$bin.names, sep = "_") # Adding the sample to the bin name
color.table <- color.table[order(color.table$bin.names),] # Ordering alphabetically as it will go in the plot

ggplot(df.long, aes(x = KO, y = Bin, fill = Copy_number)) +
    geom_tile() +
    scale_fill_gradientn(colours = c(paletteer_d("ggsci::orange_material")),
                         values = scales::rescale(c(0.001,0.005,0.2,0.4,0.6,0.9,2,4,6,11,16))) +
    xlab(label = "") +
    ylab (label = "") +
    labs(fill="Abundance\n(copy number)") +
    coord_fixed(ratio = 1) + # Aspect ratio of the tiles, () or 1 for squares
    scale_y_discrete(limits = rev(levels(as.factor(df.long$Bin)))) + # Ordering bins from A to Z
    scale_x_discrete(limits = colnames(table)) +
    theme(
        text = element_text(family = "Poppins"),
        panel.background = element_blank(),
        legend.title = element_text(size = 10, face = "bold", margin = margin(b = 15)),
        legend.position = "right",
        axis.title.x = element_text(colour="black", margin = margin(r = 10, t=5), 
                      size = 13, face = "bold"), # X axis title
        axis.title.y = element_text(colour="black", margin = margin(r = 5, t=5), 
                      size = 13, face = "bold"), # Y axis title
        axis.text.x = element_text(colour="black", margin = margin(t = 0.1, r = 0.2, l = 0.2, b = 0.1, unit = "cm"), 
                      size= 10, angle = 45, hjust = 1), # X axis numbers
        axis.text.y = element_text(colour = rev(as.vector(color.table$colors)), margin = margin(r = 5, t=2),  
                       size= 10), # Y axis numbers
        axis.line = element_blank(), # X and Y axis lines
        axis.ticks = element_blank() # X and Y axis little marks
    )

ggsave("osmorregulation_heatmap.png", width = 50, height = 30, units = "cm")

# ----------------------
# Plotting only the genes that are present in any bin

# Removing from the table columns that only have zeroes

table.pres <- table
for (gene in colnames(table.pres)) {
   if (sum(table.pres[,gene] == rep(0, length(table.pres[,gene]))) == length(table.pres[,gene])) {
      table.pres <- table.pres[,-which(colnames(table.pres) == gene)]
   }
}

df.long.pres <- as.data.frame(table.pres) %>% 
                rownames_to_column(var = "Bin") %>% 
                gather(KO, Copy_number, 2:(length(colnames(table.pres))+1))

ggplot(df.long.pres, aes(x = KO, y = Bin, fill = Copy_number)) +
    geom_tile() +
    scale_fill_gradientn(colours = c(paletteer_d("ggsci::orange_material")),
                         values = scales::rescale(c(0.001,0.005,0.2,0.4,0.6,0.9,2,4,6,11,16))) +
    xlab(label = "") +
    ylab (label = "") +
    labs(fill="Abundance\n(copy number)") +
    coord_fixed(ratio = 1) + 
    scale_y_discrete(limits = rev(levels(as.factor(df.long.pres$Bin)))) + 
    scale_x_discrete(limits = colnames(table.pres)) +
    theme(
        text = element_text(family = "Poppins"),
        panel.background = element_blank(),
        legend.title = element_text(size = 10, face = "bold", margin = margin(b = 15)),
        legend.position = "right",
        axis.title.x = element_text(colour="black", margin = margin(r = 10, t=5), 
                      size = 13, face = "bold"), 
        axis.title.y = element_text(colour="black", margin = margin(r = 5, t=5), 
                      size = 13, face = "bold"), 
        axis.text.x = element_text(colour="black", margin = margin(t = 0.1, r = 0.2, l = 0.2, b = 0.1, unit = "cm"), 
                      size= 10, angle = 45, hjust = 1), 
        axis.text.y = element_text(colour = rev(as.vector(color.table$colors)), margin = margin(r = 5, t=2),  
                       size= 10), 
        axis.line = element_blank(), 
        axis.ticks = element_blank() 
    )

ggsave("osmorregulation_heatmap_pres.png", width = 40, height = 30, units = "cm")
ggsave("osmorregulation_heatmap_pres.svg", width = 40, height = 30, units = "cm")

# ----------------------------------------------
# Barplot counting the times each gene is present in each family

# Adding taxonomy data to the data frame

tax <- read.table("taxonomia_todo.tsv", header = T, sep = "\t", as.is = T)

family <- c()

for (sample in names(sample2coolrealbin)) {
   for (i in 1:length(sample2coolrealbin[[sample]])) {
      this.family <- tax[tax[,"bin.names"] == sample2coolrealbin[[sample]][i] & tax[,"sample"] == sample, "family"]
      family <- c(family, this.family)
   }
}

table2 <- cbind(as.data.frame(table), family)

# Removing genes that aren't in any bin and for those that are present,
# transforming all values to 1 to count only presence or absence.

for (gene in colnames(table2)) {
   if (sum(table2[,gene] == rep(0, length(table2[,gene]))) == length(table2[,gene])) {
      table2 <- table2[,-which(colnames(table2) == gene)]
   } else {
      if (gene != "family") {
         table2[,gene][table2[,gene] != 0] <- 1
      }
   }  
}

df.long2 <- as.data.frame(table2) %>% 
           rownames_to_column(var = "Bin") %>% 
           gather(KO, Present, 2:(length(colnames(table2))))

ggplot(df.long2, aes(x = KO, y = Present, fill = family)) +
    geom_col(position = "stack") +
    coord_flip() +
    xlab("KO") +
    ylab("Number of occurrences") +
    scale_x_discrete(limits = df.long2[,"KO"]) +
    labs(fill = "Family") +
    theme(
        text = element_text(family = "Poppins"),
        panel.background = element_blank(),
        legend.position = "right",
        legend.key = element_rect(fill = "white", colour = "white"),
        axis.title.x = element_text(colour="black", margin = margin(r = 10, t=5), 
                      size = 17, face = "bold"),
        axis.title.y = element_text(colour="black", margin = margin(r = 10, t=10), 
                      size = 17, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1,colour="black", margin = margin(r = 5, t=5), 
                      size= 10),
        axis.text.y = element_text(colour="black", margin = margin(r = 5, t=5),  
                       size= 10),
        axis.line = element_line(size=0.75),
        axis.ticks = element_line(size=0.75)
    )+
    scale_fill_manual(values=c("#8cabd9", "#f6a7b8", "#f1ec7a", "#61bea4", "#f08838", "#9057c6", "#ffe1bd"))

ggsave("osmo_genes_family.png", width = 20, height = 30, units = "cm", bg = "transparent")
