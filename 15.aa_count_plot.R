# ---------------------------------------------------------------
# Script to create a radar plot from amino acid frequencies

# R 4.2.1

# Author: Eva Sanmartín Vázquez and Alicia García Roldán
# ---------------------------------------------------------------

# Load package

library(fmsb) # To create the radar plot

# Load data

aaPaths <- c(`20IB32` = "Count_AA/20IB32.faa.clean.fasta.aas.out", 
             `23IB32` = "Count_AA/23IB32.faa.clean.fasta.aas.out")

aaNames <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
aaCount <- data.frame(aaNames)
muestras <- list()

for (sample in names(aaPaths)) {
    muestras[[sample]] <- read.table(aaPaths[[sample]], header = T)
    a <- colnames(aaCount)
    aaCount <- cbind(aaCount, muestras[[sample]][order(muestras[[sample]]$aas),3])
    colnames(aaCount) <- c(a, sample)
}

# --------------------
# Preparing the data frame for plotting
# --------------------

# We have to transpose the data frame
aaCount.t <- t(aaCount[,-1])
colnames(aaCount.t) <- aaNames

# fmsb requires two rows at the start with the minimun and maximun
# values, so the plot is normalized
maxmin <- matrix(ncol = ncol(aaCount.t), nrow = 2, dimnames = list(c("max", "min"), aaNames))
maxmin["max",] <- 12
maxmin["min",] <- 0
df_graf <- as.data.frame(rbind(maxmin, aaCount.t))

# ---------------
# PLOT
# ---------------

colors <- c("#8cabd9", "#f6a7b8")

png("aa_spider_plot.png", width = 15, height = 15, units = "cm", res = 300)
pdf("aa_spider_plot.pdf")

par(mar = c(1, 2, 2, 1)) # Margins
par(family = "Poppins")

radarchartcirc(df_graf,
           axistype = 1, # Draw central axis label
           cglcol = "grey", cglty = 3, cglwd = 0.8, # Radar grid color, type (dashed) and width
           axislabcol = "black",
           vlcex = 1, # Size of the labels with aa names
           caxislabels = c(0, NA, 6, NA, 12), # Show max, min and middle value 
           calcex = 0.8, # Max and min size
           pty = 32, # No point
           plwd = 3, # Line width
           plty = 1, # Line type (normal line)
           pcol = colors # Line colours
)

legend(
  x = "bottomright", 
  legend = rownames(df_graf[-c(1,2),]), 
  horiz = F,
  xjust = 0, # Left justified
  bty = "n", # No box around the legend
  pch = 20, # Code for the point
  col = colors,
  text.col = "black", 
  cex = 0.8, pt.cex = 1.5 # Text and point sizes
)

dev.off()
# Embedding fonts into the PDF files so they display correctly
embed_fonts("aa_spider_plot.pdf", outfile="aa_spider_plot.pdfembed.pdf")

# -----------------------------------------
# Alternative way: using package ggradar (R 4.2.1)

library(ggplot2)
library(ggradar)
library(extrafont)

df <- as.data.frame(aaCount.t)
df <- cbind(rownames(df), df)
rownames(df) <- NULL

ggradar(df,
    font.radar = "Poppins",
    grid.min = 0,
    grid.mid = 6,
    grid.max = 12,
    draw.points = F,
    legend.position = "none",
    axis.label.size = 5,
    background.circle.colour = "white",
    group.line.width = 0.75,
    values.radar = c("0", "6", "12"),
    group.colours = c("#8cabd9", "#f6a7b8"))


ggsave("prueba.png", width = 20, height = 20, units = "cm")

