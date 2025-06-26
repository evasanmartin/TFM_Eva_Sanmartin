# ----------------------------------------------------------------------------------------
# Script to create different plots depicting medium and high quality bin statistical
# data from the different assembly strategies. Data from SQM checkm results (table
# from the first script).

# - barplots with the number of bins obtained with each method and their quality (HQ / MQ).
# - scatterplot with completeness and contamination of all bins.
# - boxplots with completeness and contamination percentages of HQ+MQ and MQ bins.
# - barplots with completeness and contamination percentages of HQ bins.

# R 4.2.1

# Author: Eva Sanmartín Vázquez
# ----------------------------------------------------------------------------------------

# Load packages
library('ggplot2')
library('extrafont')

# Load data
checkm_todo <- read.table("checkm_todo.tsv", sep = "\t")

# ---------------
# ---- PLOTS ----
# ---------------

# ---------------------------------------------------------------------------------------
# Barplots with the number of bins obtained with each method and their quality (HQ / MQ).
# ---------------------------------------------------------------------------------------

ggplot(checkm_todo, aes(x = Sample, fill = Quality)) +
    geom_bar(position = "stack") +
    xlab("Sample") +
    ylab("Bins") +
    theme(
        text = element_text(family = "Poppins"),
        panel.background = element_blank(), # White background
        legend.position = "right", # Legend to the right
        axis.title.x = element_text(colour="black", margin = margin(r = 10, t=15), 
                      size = 17, face = "bold"), # X axis title
        axis.title.y = element_text(colour="black", margin = margin(r = 10, t=10), 
                      size = 17, face = "bold"), # Y axis title
        axis.text.x = element_text(colour="black", margin = margin(r = 5, t=5), 
                      size= 15), # X axis text / numbers
        axis.text.y = element_text(colour="black", margin = margin(r = 5, t=5),  
                       size= 15), # Y axis text / numbers
        axis.line = element_line(size=0.75), # Axis lines
        axis.ticks = element_line(size=0.75) # Little lines from axis lines to numbers
    )+
    scale_fill_manual(values=c("#8cabd9", "#f6a7b8"))

ggsave("bins_HQ_MQ_barplot.png", width = 10.5, height = 10, units = "cm", bg = "transparent")
ggsave("bins_HQ_MQ_barplot.svg", width = 10.5, height = 10, units = "cm", bg = "transparent")

# ----------------------------------------------------------------------
# Scatterplot with completeness and contamination of all bins
# ----------------------------------------------------------------------

ggplot(checkm_todo, aes(x = Completeness, y = Contamination, colour = Sample)) +
    geom_point() +
    theme(
        text = element_text(family = "Poppins"),
        panel.background = element_blank(), # White background
        panel.grid.major = element_line(color ="#dcdbdb", linetype = 5, size = 0.5),
        legend.position = "right", 
        legend.key = element_rect(colour = "white"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title.x = element_text(colour="black", margin = margin(r = 10, t=15), 
                      size = 17, face = "bold"), # X axis title
        axis.title.y = element_text(colour="black", margin = margin(r = 10, t=10), 
                      size = 17, face = "bold"), # Y axis title
        axis.text.x = element_text(colour="black", margin = margin(r = 5, t=5), 
                      size= 15), # X axis text / numbers
        axis.text.y = element_text(colour="black", margin = margin(r = 5, t=5),  
                      size= 15), # Y axis text / numbers
        axis.line = element_line(size=0.75), # Axis lines
        axis.ticks = element_line(size=0.75) # Little lines from axis lines to numbers
    ) +
    scale_x_continuous(limits = c(50,100)) + # Axis limits
    scale_colour_manual(values=c("#8cabd9", "#f6a7b8", "#61bea4")) + # Point colours
    annotate("rect", xmin=90, xmax=100, ymin=0 , ymax=5, alpha = 0.2, fill ="#f6ca94") # Adding a rectangle

ggsave("completeness_contamination_scatterplot.png", width = 13, height = 10, units = "cm", bg = "transparent")
ggsave("completeness_contamination_scatterplot.svg", width = 13, height = 10, units = "cm")


# -------------------------------------------------------------------------------
# Boxplots with completeness and contamination percentages of HQ + MQ bins
# -------------------------------------------------------------------------------

# Contamination

ggplot(checkm_todo, aes(x = Sample, y = Contamination, fill = Sample)) +
    stat_boxplot(geom = "errorbar", width = 0.2, lwd = 0.75) + # Adding "stops" to the end of the bars
    geom_boxplot(width = 0.75, lwd = 0.75, fatten = 1.5) +
    xlab("Sample") +
    stat_summary(fun = "mean", shape = 4) + # Adding the mean
    theme(
        text = element_text(family = "Poppins"),
        panel.background = element_blank(), # White background
        legend.position = "none", # No legend
        axis.title.x = element_text(colour="black", margin = margin(r = 10, t=15), 
                      size = 17, face = "bold"), # X axis title
        axis.title.y = element_text(colour="black", margin = margin(r = 10, t=10), 
                      size = 17, face = "bold"), # Y axis title
        axis.text.x = element_text(colour="black", margin = margin(r = 5, t=5), 
                      size= 15), # X axis text / numbers
        axis.text.y = element_text(colour="black", margin = margin(r = 5, t=5),  
                      size= 15), # Y axis text / numbers
        axis.line = element_line(size=0.75), # Axis lines
        axis.ticks = element_line(size=0.75) # Little lines from axis lines to numbers
    )+
    scale_fill_manual(values=c("#8cabd9", "#f6a7b8", "#f1ec7a")) 

ggsave("contamination_boxplot.png", width = 10, height = 10, units = "cm", bg = "transparent")
ggsave("contamination_boxplot.svg", width = 10, height = 10, units = "cm")

# Completeness

ggplot(checkm_todo, aes(x = Sample, y = Completeness, fill = Sample)) +
    stat_boxplot(geom = "errorbar", width = 0.2, lwd = 0.75) +
    geom_boxplot(width = 0.75, lwd = 0.75, fatten = 1.5) +
    xlab("Sample") +
    stat_summary(fun = "mean", shape = 4) + # Adding the mean
    theme(
        text = element_text(family = "Poppins"),
        panel.background = element_blank(), # White background
        legend.position = "none", # No legend
        axis.title.x = element_text(colour="black", margin = margin(r = 10, t=15), 
                      size = 17, face = "bold"), # X axis title
        axis.title.y = element_text(colour="black", margin = margin(r = 10, t=10), 
                      size = 17, face = "bold"), # Y axis title
        axis.text.x = element_text(colour="black", margin = margin(r = 5, t=5), 
                      size= 15), # X axis text / numbers
        axis.text.y = element_text(colour="black", margin = margin(r = 5, t=5),  
                      size= 15), # Y axis text / numbers
        axis.line = element_line(size=0.75), # Axis lines
        axis.ticks = element_line(size=0.75) # Little lines from axis lines to numbers
    )+
    scale_fill_manual(values=c("#8cabd9", "#f6a7b8", "#f1ec7a"))

ggsave("completeness_boxplot.png", width = 10, height = 10, units = "cm", bg = "transparent")
ggsave("completeness_boxplot.svg", width = 10, height = 10, units = "cm")

# -------------------------------------------------------------------------------
# Boxplots with completeness and contamination percentages of MQ bins
# -------------------------------------------------------------------------------

# Contamination

ggplot(subset(checkm_todo, Quality == "MQ"), aes(x = Sample, y = Contamination, fill = Sample)) +
    stat_boxplot(geom = "errorbar", width = 0.2, lwd = 0.75) +
    geom_boxplot(width = 0.75, lwd = 0.75, fatten = 1.5) +
    xlab("Sample") +
    stat_summary(fun = "mean", shape = 4) + # Adding the mean
    theme(
        text = element_text(family = "Poppins"),
        panel.background = element_blank(), # White background
        legend.position = "none", # No legend
        axis.title.x = element_text(colour="black", margin = margin(r = 10, t=15), 
                      size = 17, face = "bold"), # X axis title
        axis.title.y = element_text(colour="black", margin = margin(r = 10, t=10), 
                      size = 17, face = "bold"), # Y axis title
        axis.text.x = element_text(colour="black", margin = margin(r = 5, t=5), 
                      size= 15), # X axis text / numbers
        axis.text.y = element_text(colour="black", margin = margin(r = 5, t=5),  
                      size= 15), # Y axis text / numbers
        axis.line = element_line(size=0.75), # Axis lines
        axis.ticks = element_line(size=0.75) # Little lines from axis lines to numbers
    )+
    scale_fill_manual(values=c("#8cabd9", "#f6a7b8", "#f1ec7a"))


ggsave("contamination_boxplot_MQ.png", width = 10, height = 10, units = "cm", bg = "transparent")

# Completeness

ggplot(subset(checkm_todo, Quality == "MQ"), aes(x = Sample, y = Completeness, fill = Sample)) +
    stat_boxplot(geom = "errorbar", width = 0.2, lwd = 0.75) +
    geom_boxplot(width = 0.75, lwd = 0.75, fatten = 1.5) +
    xlab("Sample") +
    stat_summary(fun = "mean", shape = 4) + # Adding the mean
    theme(
        text = element_text(family = "Poppins"),
        panel.background = element_blank(), # White background
        legend.position = "none", # No legend
        axis.title.x = element_text(colour="black", margin = margin(r = 10, t=15), 
                      size = 17, face = "bold"), # X axis title
        axis.title.y = element_text(colour="black", margin = margin(r = 10, t=10), 
                      size = 17, face = "bold"), # Y axis title
        axis.text.x = element_text(colour="black", margin = margin(r = 5, t=5), 
                      size= 15), # X axis text / numbers
        axis.text.y = element_text(colour="black", margin = margin(r = 5, t=5),  
                      size= 15), # Y axis text / numbers
        axis.line = element_line(size=0.75), # Axis lines
        axis.ticks = element_line(size=0.75) # Little lines from axis lines to numbers
    )+
    scale_fill_manual(values=c("#8cabd9", "#f6a7b8", "#f1ec7a"))


ggsave("completeness_boxplot_MQ.png", width = 10, height = 10, units = "cm", bg = "transparent")

# ----------------------------------------------------------------------
# Barplots with completeness and contamination percentages of HQ bins
# ----------------------------------------------------------------------

# Contamination

ggplot(subset(checkm_todo, Quality == "HQ"), aes(x = Sample, y = Contamination, fill = Sample)) +
    geom_bar(stat = "summary", fun = "mean") + # Height of the bars corresponds to the mean
    geom_point(data = subset(checkm_todo, (Sample == "23IB32"|Sample == "CO") & Quality == "HQ"),
               position = position_identity()) + # Adding centered points
    geom_point(data = subset(checkm_todo, Sample == "20IB32" & Quality == "HQ"), 
               position = position_jitter(width = 0.2, height = 0, seed = 1)) + # Adding points with jitter
    xlab("Sample") +
    theme(
        text = element_text(family = "Poppins"),
        panel.background = element_blank(), # White background
        legend.position = "none", # No legend
        axis.title.x = element_text(colour="black", margin = margin(r = 10, t=15), 
                      size = 17, face = "bold"), # X axis title
        axis.title.y = element_text(colour="black", margin = margin(r = 10, t=10), 
                      size = 17, face = "bold"), # Y axis title
        axis.text.x = element_text(colour="black", margin = margin(r = 5, t=5), 
                      size= 15), # X axis text / numbers
        axis.text.y = element_text(colour="black", margin = margin(r = 5, t=5),  
                      size= 15), # Y axis text / numbers
        axis.line = element_line(size=0.75), # Axis lines
        axis.ticks = element_line(size=0.75) # Little lines from axis lines to numbers
    )+
    scale_fill_manual(values=c("#8cabd9", "#f6a7b8", "#f1ec7a"))

ggsave("contamination_barplot_HQ.png", width = 10, height = 10, units = "cm", bg = "transparent")

# Completeness

ggplot(subset(checkm_todo, Quality == "HQ"), aes(x = Sample, y = Completeness, fill = Sample)) +
    geom_bar(stat = "summary", fun = "mean") + # Height of the bars corresponds to the mean
    geom_point(data = subset(checkm_todo, (Sample == "23IB32"|Sample == "CO") & Quality == "HQ"),
               position = position_identity()) + # Adding centered points
    geom_point(data = subset(checkm_todo, Sample == "20IB32" & Quality == "HQ"), 
               position = position_jitter(width = 0.2, height = 0, seed = 1)) + # Adding points with jitter
    xlab("Sample") +
    ylim(c(0,100)) +
    theme(
        text = element_text(family = "Poppins"),
        panel.background = element_blank(), # White background
        legend.position = "none", # No legend
        axis.title.x = element_text(colour="black", margin = margin(r = 10, t=15), 
                      size = 17, face = "bold"), # X axis title
        axis.title.y = element_text(colour="black", margin = margin(r = 10, t=10), 
                      size = 17, face = "bold"), # Y axis title
        axis.text.x = element_text(colour="black", margin = margin(r = 5, t=5), 
                      size= 15), # X axis text / numbers
        axis.text.y = element_text(colour="black", margin = margin(r = 5, t=5),  
                      size= 15), # Y axis text / numbers
        axis.line = element_line(size=0.75), # Axis lines
        axis.ticks = element_line(size=0.75) # Little lines from axis lines to numbers
    )+
    scale_fill_manual(values=c("#8cabd9", "#f6a7b8", "#f1ec7a"))

ggsave("completeness_barplot_HQ.png", width = 10, height = 10, units = "cm", bg = "transparent")

