# --------------------------------------------------------------------
# This scripts finds the common taxa between the samples, the 
# intersections and creates a Venn diagram. Uses the taxonomy tables.

# R 4.2.1

# Author: Eva Sanmartín Vázquez

# --------------------------------------------------------------------

# Load packages

library('ggplot2')
library('extrafont')
library('ggvenn')

# Load data (taxonomy tables)

paths <- c(`20IB32` = "~/20IB32_taxonomy.tsv",
           `23IB32` = "~/23IB32_taxonomy.tsv",
            CO = "~/CO_taxonomy.tsv")

taxonomy <- list()

for (sample in names(paths)) {
    taxonomy[[sample]] <- read.table(paths[sample], sep = "\t")
}

# -------------------------
# PREPARING THE DATA
# -------------------------

# Removing unassigned species and genera (s__ and g__)

no.NA <- list()

for (sample in names(paths)) {
    no.NA$especies[[sample]] <- taxonomy[[sample]]$species[taxonomy[[sample]]$species!="s__"]
    no.NA$generos[[sample]] <- taxonomy[[sample]]$genus[taxonomy[[sample]]$genus!="g__"]
    no.NA$familias[[sample]] <- taxonomy[[sample]]$family[taxonomy[[sample]]$family!="f__"]
}

# Transforming into a data frame in order to plot with ggvenn

df <- list()
for (taxon in names(no.NA)) {
   df[[taxon]] <- list_to_data_frame(no.NA[[taxon]])
}

# ---------------------------------------
# SUBSETS
# ---------------------------------------

common <- list()

for (taxon in names(no.NA)) {
    # Common to the three samples
    common[[taxon]][["20.23.CO"]] <- Reduce(intersect, no.NA[[taxon]])
    # Common two by two. First we intersect two sets (intersect) 
    # and then we remove the part that is common to the three (setdiff)
    common[[taxon]][["20.23"]] <- setdiff(intersect(no.NA[[taxon]]$`20IB32`, no.NA[[taxon]]$`23IB32`), 
                                          common[[taxon]]$`20.23.CO`)
    common[[taxon]][["20.CO"]] <- setdiff(intersect(no.NA[[taxon]]$`20IB32`, no.NA[[taxon]]$CO), 
                                          common[[taxon]]$`20.23.CO`)
    common[[taxon]][["23.CO"]] <- setdiff(intersect(no.NA[[taxon]]$`23IB32`, no.NA[[taxon]]$CO), 
                                          common[[taxon]]$`20.23.CO`)
    # Present in only one sample: removing (setdiff) the three intersections with the other groups
    common[[taxon]][["20"]] <- setdiff(setdiff(setdiff(no.NA[[taxon]]$`20IB32`, common[[taxon]][["20.23"]]), 
                                               common[[taxon]][["20.CO"]]), 
                                       common[[taxon]]$`20.23.CO`)
    common[[taxon]][["23"]] <- setdiff(setdiff(setdiff(no.NA[[taxon]]$`23IB32`, common[[taxon]][["20.23"]]), 
                                               common[[taxon]][["23.CO"]]), 
                                       common[[taxon]]$`20.23.CO`)
    common[[taxon]][["CO"]] <- setdiff(setdiff(setdiff(no.NA[[taxon]]$CO, common[[taxon]][["20.CO"]]), 
                                               common[[taxon]][["23.CO"]]),
                                       common[[taxon]]$`20.23.CO`)
}


# ---------
# PLOTS
# ---------

# Creating three plots at three different levels: families, genera and species
save.names <- c("venn_species.png", "venn_genera.png", "venn_families.png")
names(save.names) <- names(no.NA)


for (taxon in names(save.names)) { # Creating the three plots at once with the loop
   ggplot(df[[taxon]], aes(A = `20IB32`, B = `23IB32`, C = CO)) + # Indicate the three sets
    geom_venn(
        fill_color = c("#8cabd9", "#f6a7b8", "#f1ec7a"),
        fill_alpha = 0.75,
        stroke_size = 0.5,
        show_percentage = F,
        set_name_size = 5,
        text_size = 5) +
    theme(
        panel.background = element_blank(), # White background
        legend.position = "none", # No legend
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()
    ) 
    
    ggsave(save.names[taxon], width = 10, height = 10, units = "cm", bg = "white")
}
