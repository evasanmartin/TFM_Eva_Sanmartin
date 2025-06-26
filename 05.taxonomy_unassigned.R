# -------------------------------------------------------------------------------------
# Script that counts how many bins have an unassigned taxa at a certain level and 
# creates a barplot. Uses GTDBtk taxonomy table (script 4).

# R 4.2.1

# Author: Eva Sanmartín Vázquez

# -------------------------------------------------------------------------------------

# Load packages
library('ggplot2')
library('extrafont')
library('dplyr')
library('tibble')

# Load data
taxonomy <- read.table("taxonomia_todo.tsv", sep = "\t", header = T)

# -------------------------------------------
# Count how many unassigned taxa there are
# -------------------------------------------

NAs <- data.frame(matrix(nrow = 3, ncol = 3))
colnames(NAs) <- c("family", "genus", "species")
rownames(NAs) <- c("20IB32", "23IB32", "CO")

for (sample in rownames(NAs)) {
   NAs[sample, "species"] <- sum(taxonomy$sample == sample & taxonomy$species == "s__")
   NAs[sample, "genus"] <- sum(taxonomy$sample == sample & taxonomy$genus == "g__")
   NAs[sample, "family"] <- sum(taxonomy$sample == sample & taxonomy$family == "f__")
}

# ----------------
# Plot
# ----------------

# We can generate several plots at once with a for loop. To do that, we have 
# to write the variables that change in each step in a matrix:
# variable to plot, axis name and file name to save the plot.

save.names <- c("NA_species.png", "NA_genus.png")
y.labs <- c("NA species", "NA genus")
graph.prop <- cbind(save.names, y.labs)
rownames(graph.prop) <- c("species", "genus")

for (taxon in rownames(graph.prop)) {
    var<-sym(taxon) # This is necessary so it interprets it as a value and not a character string
    ggplot(NAs, aes(x = sample, y = !!var, fill = sample)) + # Write !! so it gets the value of taxon
        geom_bar(stat = "identity") +
        xlab("Sample") +
        ylab(graph.prop[taxon,"y.labs"]) +
        theme(
            text = element_text(family = "Poppins"),
            panel.background = element_blank(), # No background
            legend.position = "none", # No legend
            axis.title.x = element_text(colour="black", margin = margin(r = 10, t=15), 
                        size = 17, face = "bold"), # X axis title
            axis.title.y = element_text(colour="black", margin = margin(r = 10, t=10), 
                        size = 17, face = "bold"), # Y axis title
            axis.text.x = element_text(colour="black", margin = margin(r = 5, t=5), 
                        size= 15), # X axis numbers
            axis.text.y = element_text(colour="black", margin = margin(r = 5, t=5),  
                        size= 15), # Y axis numbers
            axis.line = element_line(size=0.75), # Axis lines
            axis.ticks = element_line(size=0.75) # Lines from the axis lines to the numbers
        )+
        scale_fill_manual(values=c("#8cabd9", "#f6a7b8", "#f1ec7a"))

    ggsave(graph.prop[taxon, "save.names"], width = 10, height = 10, units = "cm")
}


# ---------------------------------------------
# Calculating percentages instead of absolute values
# ---------------------------------------------

NAs.percent <- NAs

for (sample in rownames(NAs)) {
   NAs.percent[sample, "species"] <- NAs[sample, "species"]/sum(taxonomy$sample == sample)
   NAs.percent[sample, "genus"] <- NAs[sample, "genus"]/sum(taxonomy$sample == sample)
   NAs.percent[sample, "family"] <- NAs[sample, "family"]/sum(taxonomy$sample == sample)
}

NAs.percent <- rownames_to_column(NAs.percent, "sample")

ggplot(NAs.percent, aes(x = sample, y = genus, fill = sample)) + 
    geom_bar(stat = "identity") +
    xlab("Sample") +
    ylab("NA genus (% bins)") +
    theme(
        text = element_text(family = "Poppins"),
        panel.background = element_blank(), # Fondo blanco
        legend.position = "none", # Sin leyenda
        axis.title.x = element_text(colour="black", margin = margin(r = 10, t=15), 
                    size = 17, face = "bold"), # Título del eje X
        axis.title.y = element_text(colour="black", margin = margin(r = 10, t=10), 
                    size = 17, face = "bold"), # Título del eje Y
        axis.text.x = element_text(colour="black", margin = margin(r = 5, t=5), 
                    size= 15), # Números del eje X
        axis.text.y = element_text(colour="black", margin = margin(r = 5, t=5),  
                    size= 15), # Números del eje Y
        axis.line = element_line(size=0.75), # Líneas de los ejes
        axis.ticks = element_line(size=0.75), # Rayitas de los números de los ejes
        strip.background = element_blank(),
        strip.text = element_text(size = 15)
    )+
    scale_fill_manual(values=c("#8cabd9", "#f6a7b8", "#f1ec7a"))

ggsave("NAs_perc_genus.png", width = 10, height = 10, units = "cm")
ggsave("NAs_perc_genus.svg", width = 10, height = 10, units = "cm")

ggplot(NAs.percent, aes(x = sample, y = species, fill = sample)) + 
    geom_bar(stat = "identity") +
    xlab("Sample") +
    ylab("NA species (% bins)") +
    labs(fill = "Sample") +
    theme(
        text = element_text(family = "Poppins"),
        panel.background = element_blank(), # Fondo blanco
        legend.position = "none", # Sin leyenda
        axis.title.x = element_text(colour="black", margin = margin(r = 10, t=15), 
                    size = 17, face = "bold"), # Título del eje X
        axis.title.y = element_text(colour="black", margin = margin(r = 10, t=10), 
                    size = 17, face = "bold"), # Título del eje Y
        axis.text.x = element_text(colour="black", margin = margin(r = 5, t=5), 
                    size= 15), # Números del eje X
        axis.text.y = element_text(colour="black", margin = margin(r = 5, t=5),  
                    size= 15), # Números del eje Y
        axis.line = element_line(size=0.75), # Líneas de los ejes
        axis.ticks = element_line(size=0.75), # Rayitas de los números de los ejes
        strip.background = element_blank(),
        strip.text = element_text(size = 15)
    )+
    scale_fill_manual(values=c("#8cabd9", "#f6a7b8", "#f1ec7a"))

ggsave("NAs_perc_species.png", width = 10, height = 10, units = "cm")
ggsave("NAs_perc_species.svg", width = 10, height = 10, units = "cm")

# We can also create both plots at once using facetting, but they are
# saved together in a file instead of several files with each plot.

library('tidyr')

df.long <- NAs.percent %>% gather(taxon, value, family:species)

ggplot(df.long, aes(x = sample, y = value, fill = sample)) + 
    geom_bar(stat = "identity") +
    xlab("Sample") +
    ylab("NAs (% bins)") +
    theme(
        text = element_text(family = "Poppins"),
        panel.background = element_blank(), # Fondo blanco
        legend.position = "none", # Sin leyenda
        axis.title.x = element_text(colour="black", margin = margin(r = 10, t=15), 
                    size = 17, face = "bold"), # Título del eje X
        axis.title.y = element_text(colour="black", margin = margin(r = 10, t=10), 
                    size = 17, face = "bold"), # Título del eje Y
        axis.text.x = element_text(colour="black", margin = margin(r = 5, t=5), 
                    size= 15), # Números del eje X
        axis.text.y = element_text(colour="black", margin = margin(r = 5, t=5),  
                    size= 15), # Números del eje Y
        axis.line = element_line(size=0.75), # Líneas de los ejes
        axis.ticks = element_line(size=0.75), # Rayitas de los números de los ejes
        strip.background = element_blank(),
        strip.text = element_text(size = 15)
    )+
    facet_wrap(~taxon, scales = "free", 
    labeller = as_labeller(c("family" = "Family", "genus" = "Genus", "species" = "Species"))) + # Change grid titles
    scale_fill_manual(values=c("#8cabd9", "#f6a7b8", "#f1ec7a"))

ggsave("NA_percent.png", width = 30, height = 10, units = "cm", bg = "transparent")
