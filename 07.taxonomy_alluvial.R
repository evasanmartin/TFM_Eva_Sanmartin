# ---------------------------------------------------------------------
# This script creates a Sankey diagram / alluvial plot representing
# the taxonomy of the samples from the table with GTDBtk taxonomy
# (script 4).

# R 4.2.1

# Author: Eva Sanmartín Vázquez
# ---------------------------------------------------------------------

# ---- Load packages -----

library('ggplot2')
library('ggalluvial')
library('extrafont')

# ---- Load data (taxonomy table) ----

taxonomy <- read.table("taxonomia_todo.tsv", sep = "\t", header = T, stringsAsFactors = F)

# ---- PREPARE DATASET FOR THE PLOT ----

# In order to have the strata (taxa) ordered in each level, 
# we have to order the data frame, turning the taxa columns into
# factor type and assigning the levels in the order we want.

# More info about ordering:
# https://matthewdharris.com/2017/11/11/a-brief-diversion-into-static-alluvial-sankey-diagrams-in-r/

ordered <- taxonomy[order(taxonomy$domain),] # First we sort by domain
order.phylum <- c("p__Halobacteriota", "p__Nanohaloarchaeota", "p__Bacteroidota", 
                  "p__Pseudomonadota", "p__Patescibacteria")
ordered$phylum <- factor(ordered$phylum, levels = order.phylum)

order.class <- c("c__Halobacteria", "c__Nanosalinia", "c__Rhodothermia", "c__Bacteroidia",
                 "c__Gammaproteobacteria", "c__Paceibacteria")
ordered$class <- factor(ordered$class, levels = order.class)

order.order <- c("o__Halobacteriales", "o__Nanosalinales", "o__Rhodothermales", 
                 "o__CAILMK01", "o__Nitrococcales", "o__UBA9973")
ordered$order <- factor(ordered$order, levels = order.order)

order.family <- c("f__Haloarculaceae", "f__Haloferacaceae", "f__Halobacteriaceae",
                  "f__Nanosalinaceae", "f__Salinibacteraceae", "f__JAAYUY01",
                  "f__Nitrococcaceae", "f__SW-6-46-9")
ordered$family <- factor(ordered$family, levels = order.family)

order.genus <- c("g__Natronomonas","g__Halosegnis","g__Halovenus","g__Salinirubellus",
                 "g__Salinirussus","g__Haloarcula","g__Haloglomus",
                 "g__Haloquadratum","g__A07HB70","g__Haloplanus","g__Halobaculum",
                 "g__Halobellus","g__Halorubrum","g__B1-Br10-U2g21", "g__B1-Br10-U2g19",
                 "g__","g__Salinibacter", "g__Spiribacter", "g__SW-6-46-9")
ordered$genus <- factor(ordered$genus, levels = order.genus)

order.species <-c( "s__Natronomonas aquatica","s__Natronomonas salsuginis","s__Halosegnis rubeus",
                  "s__Halovenus sp028275145", "s__Salinirussus sp028275005",
                  "s__Haloglomus irregulare", "s__Haloquadratum walsbyi", 
                  "s__Haloquadratum sp028275165","s__A07HB70 sp937138675" ,"s__",
                  "s__Halorubrum sp937146585", "s__Halorubrum rutilum", "s__Salinibacter pepae")
ordered$species <- factor(ordered$species, levels = order.species)

# ---- PLOT ----

ggplot(ordered, aes(axis1 = domain, axis2 = phylum, axis3 = class, 
    axis4 = order, axis5 = family, axis6 = genus, axis7 = species, fill = sample)) +
    geom_alluvium(aes(fill = sample), aes.bind = "alluvia", decreasing = NA, lode.guidance = "forward") + 
    # aes.bind = "alluvia" groups the same colors (samples) together
    geom_stratum(width = 1/20, fill = "white", alpha = 0, col = "grey", decreasing = NA) + # Draw rectangles for the strata
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2, decreasing = NA) + # Write strata (taxa) names
    scale_x_discrete(limits = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) +
    theme(
    text = element_text(family = "Poppins"),
        panel.background = element_blank(), # White background
        legend.position = "right", # Legend to the right
        legend.key = element_rect(fill = "white", colour = "white"), # White background for the legend
        axis.title.x = element_text(colour="black", margin = margin(r = 10, t=15), 
                      size = 17, face = "bold"), # X axis title
        axis.text.x = element_text(colour="black", margin = margin(r = 5, t=0), 
                      size= 10), # X axis text (taxa levels)
        axis.text.y = element_blank(), # Remove Y axis text
        axis.line = element_blank(), # Remove axis lines
        axis.ticks = element_blank() # Remove axis ticks
    ) +
    scale_fill_manual(values = c("#8cabd9", "#f6a7b8", "#f1ec7a"))

ggsave("alluvial_ordered.png", width = 30, height = 15, units = "cm")
ggsave("alluvial_ordered.svg", width = 30, height = 15, units = "cm")
