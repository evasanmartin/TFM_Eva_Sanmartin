# -----------------------------------------------------------------------
# This script creates a beautiful heatmap from biotin metabolism gene 
# abundances and adds a barplot with completeness data. Uses biotin
# gene table obtained with previous script.

# R 4.2.1

# Author: Eva Sanmartín Vázquez
# -----------------------------------------------------------------------

# ---- Loading packages ----

library('tidyr')
library('tibble')
library('ggplot2')
library('extrafont')
library('paletteer')
library('ggpubr')

# -----------------
# ---- DATASET ----
# -----------------

# Reading data and transforming to long format

df <- read.table("Funciones/biotin_table.tsv", sep = "\t", check.names = F)

df.long <- as.data.frame(df) %>% 
           rownames_to_column(var = "Bin") %>% 
           gather(KO, Copy_number, 2:(length(colnames(df))+1))

df2 <- read.table("Funciones/biotin_comp_table.tsv", sep = "\t", check.names = F)

df.long2 <- as.data.frame(df2) %>% 
            rownames_to_column(var = "Bin") %>% 
            gather(KO, Copy_number, 2:(length(colnames(df))+1))


# --------------
# ---- PLOT ----
# --------------

# ------------------------
# Heatmap (no barplot)

heatmap <- ggplot(df.long, aes(x = KO, y = Bin, fill = Copy_number)) +
    geom_tile(colour = "white") +
    scale_fill_gradientn(colours = paletteer_d("ggsci::orange_material"),
                         values = scales::rescale(c(0,0.005,0.2,0.4,0.6,0.9,2,4,6,11,16))) +
    xlab(label = "") +
    ylab (label = "") +
    labs(fill="Abundance\n(copy number)") +
    coord_fixed(ratio = 1/1.5) + # Proportion of the rectangles
    scale_y_discrete(limits = rev(levels(as.factor(df.long$Bin)))) + # Ordering bins alphabetically
    scale_x_discrete(limits = as.factor(colnames(df))) + # Ordering genes
    theme(
        text = element_text(family = "Poppins"),
        panel.background = element_blank(),
        legend.title = element_text(size = 10, face = "bold", margin = margin(b = 15)), # hjust = 0.5 for centered position
        legend.position = "bottom", # Bottom if we combine it with barplot, or else right
        axis.title.x = element_text(colour="black", margin = margin(r = 10, t=5), 
                      size = 13, face = "bold"), # X axis title
        axis.title.y = element_text(colour="black", margin = margin(r = 5, t=5), 
                      size = 13, face = "bold"), # Y axis title
        axis.text.x = element_text(colour="black", margin = margin(t = 0.1, r = 0.2, l = 0.2, b = 0.1, unit = "cm"), 
                      size= 10, angle = 45, hjust = 1), # X axis text
        axis.text.y = element_text(colour="black", margin = margin(r = 5, t=2),  
                       size= 10), # Y axis text
        axis.line = element_blank(), # Remove axis lines
        axis.ticks = element_blank(), # Remove little lines from the axis to text
        plot.margin = margin(t = 20, r = 40)
    )

ggsave("biotin_heatmap.png", width = 20, height = 25, units = "cm")
ggsave("biotin_heatmap.svg", width = 20, height = 25, units = "cm")

# ------------------------
# Barplot

barplot <- ggplot(df.long2, aes(x = Completeness/17, y = Bin)) +
    geom_col(fill = "#FFCC80FF") +
    xlab(label = "Completeness (%)") +
    ylab(label = "") +
    scale_y_discrete(limits = rev(levels(as.factor(df.long2$Bin)))) + # Ordering bins alphabetically
    scale_x_continuous(limits = c(0,100), breaks = seq(0,100,25)) +
    theme(
        text = element_text(family = "Poppins"),
        panel.background = element_blank(),
        panel.grid.major.x = element_line(colour = "#c9c9c9"),
        legend.position = "none",
        axis.title.x = element_text(colour="black", margin = margin(r = 10, t=5), 
                      size = 12, face = "bold"), # X axis title
        axis.text.x =  element_text(colour="black", margin = margin(t = 0.1, r = 0.2, l = 0.2, b = 0.1, unit = "cm"), 
                      size= 9), # X axis text
        axis.text.y = element_blank(), # Remove Y axis text
        axis.line.x = element_line(), # Print X axis line
        axis.ticks.x = element_line(), # Print X axis little lines
        axis.line.y = element_blank(), # Remove Y axis line
        axis.ticks.y = element_blank(), # Remove little lines from the axis to text
        plot.margin = margin(t = 20, r = 40)
    )

ggarrange(heatmap, barplot, ncol = 2, nrow = 1, widths = c(1,0.2), align = "h")

ggsave("biotin_heatmap_barplot.png", width = 25, height = 30, units = "cm", bg = "white")
ggsave("biotin_heatmap_barplot.svg", width = 25, height = 30, units = "cm")
