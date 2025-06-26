# ----------------------------------------------------------
# This script plots quast quality results.
# Every parameter is plotted in a single file with facetting.

# R 4.2.1

# Author: Eva Sanmartín Vázquez

# ----------------------------------------------------------

# ---- LOAD PACKAGES ----
library('ggplot2')
library('extrafont')
library('tidyr')

# ---- LOAD DATA ----
quast <- read.table("quast.tsv", sep = "\t")

# Transforming to "long format" so that we can group by the 
# columns of the orginial data frame (metrics).
df_long <- quast %>% gather(metric, value, n.contigs.0.bp:Ns.per.100.kbp)


# -------------------
# ---- P L O T S ----
# -------------------

## ---- Boxplots with ALL RESULTS AT THE SAME TIME ----

pdf("quast_todo_boxplot.pdf", width = 17.7, height = 17.7)

ggplot(df_long, aes(x = sample, y = value, fill = sample)) +
    stat_boxplot(geom = "errorbar", width = 0.2, lwd = 0.5) + # Adding ends to the lines
    geom_boxplot(width = 0.5, lwd = 0.5, fatten = 1) +
    xlab("Sample") +
    stat_summary(fun = "mean", shape = 4, cex = 0.2) + # Adding the mean
    theme(
        text = element_text(family = "Poppins"),
        panel.background = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(colour="black", margin = margin(r = 10, t=15), 
                      size = 17, face = "bold"),
        axis.title.y = element_text(colour="black", margin = margin(r = 10, t=10), 
                      size = 17, face = "bold"),
        axis.text.x = element_text(colour="black", margin = margin(r = 5, t=5), 
                      size= 15), 
        axis.text.y = element_text(colour="black", margin = margin(r = 5, t=5),  
                      size= 15), 
        axis.line = element_line(size=0.75),
        axis.ticks = element_line(size=0.75) 
    ) +
    facet_wrap(~metric, scales = "free") +
    scale_fill_manual(values=c("#8cabd9", "#f6a7b8", "#f1ec7a")) 

dev.off()

ggsave("quast_todo_boxplot.png", width = 45, height = 45, units = "cm", bg = "transparent")

## ---- Violin plots with ALL RESULTS AT THE SAME TIME ----

ggplot(df_long, aes(x = sample, y = value, fill = sample)) +
    geom_violin() +
    xlab("Sample") +
    theme(
        text = element_text(family = "Poppins"),
        panel.background = element_blank(), 
        legend.position = "none", 
        axis.title.x = element_text(colour="black", margin = margin(r = 10, t=15), 
                      size = 17, face = "bold"),
        axis.title.y = element_text(colour="black", margin = margin(r = 10, t=10), 
                      size = 17, face = "bold"),
        axis.text.x = element_text(colour="black", margin = margin(r = 5, t=5), 
                      size= 15), 
        axis.text.y = element_text(colour="black", margin = margin(r = 5, t=5),  
                      size= 15),
        axis.line = element_line(size=0.75), 
        axis.ticks = element_line(size=0.75) 
    )+
    facet_wrap(~metric, scales = "free") +
    scale_fill_manual(values=c("#8cabd9", "#f6a7b8", "#f1ec7a")) 

ggsave("quast_todo_violinplot.png", width = 45, height = 45, units = "cm", bg = "transparent")
