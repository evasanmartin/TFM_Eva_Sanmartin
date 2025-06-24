# This script creates a violin plot from isoelectric point tables
# calculated in the previous script.

# R 4.2.1

# Author: Eva Sanmartín Vázquez and Alicia García Roldán
# ----------------------------------------------------------------

# ---- LOAD PACKAGES ----

library(ggplot2)
library(extrafont)

# ---- LOAD DATA ----

df = read.table("~/PI/Isoelectric_point_ArchBact.txt", header=T, sep="\t")

# ---- PLOT ----

ggplot(df, aes(y = PI, x = Domain, fill = Domain, colour = Domain)) +
  geom_violin (linewidth = 0.5) + 
  scale_y_continuous(breaks = seq(2,14, by = 2)) + 
  scale_x_discrete(limits = levels(as.factor(df$Domain)), labels = c("Archaea", "Bacteria")) +
  ylab("Isoelectric point") + 
  theme(
    text = element_text(family = "Poppins"),
    panel.background = element_blank(), # White background
    legend.position = "none", # No legend 
    axis.title.x = element_text(colour="black", margin = margin(r = 10, t=15), 
                      size = 17, face = "bold"), # X axis title
    axis.title.y = element_text(colour="black", margin = margin(r = 10, t=10), 
                      size = 17, face = "bold"), # Y axis title
    axis.text.x = element_text(colour="black", margin = margin(r = 5, t=5), 
                      size= 15, face = "italic"), # X axis numbers
    axis.text.y = element_text(colour="black", margin = margin(r = 5, t=5),  
                      size= 15), # Y axis numbers
    axis.line = element_line(size=0.75), # Axis lines
    axis.ticks = element_line(size=0.75) # Axis lines to the numbers
  ) + 
  scale_fill_manual(values=c("#61bea4", "#9057c6")) +
  scale_colour_manual(values=c("#398b74", "#632d95"))

ggsave("PI_archbact.png", width = 15, height = 15, units = "cm", bg = "transparent")
ggsave("PI_archbact.svg", width = 15, height = 15, units = "cm", bg = "transparent")
