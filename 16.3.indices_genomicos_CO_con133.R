# ---------------------------------------------------------------------
# This script creates heatmaps of OGRIs results (GGDC, ANI, AAI).
# Based on Manzanilla script by Alicia García Roldán & Cristina Galisteo.

# R 4.2.1

# Authors: Eva Sanmartín Vázquez, Alicia García Roldán & Cristina Galisteo
# ---------------------------------------------------------------------

# ---- LOAD PACKAGES ----
library('ggplot2')
library('paletteer')
library('tidyr')
library('tibble')
library('stringr')
library('extrafont')

# -----------------------------
# ---- DATASET PREPARATION ----
# -----------------------------

## ---- GGDC (DDH) ----

GGDC <- read.table("Indices_genomicos/CO_concoct.133/Resultados/GGDC_formula3.txt", header = T)
samples <- colnames(GGDC) # Getting sample names from GGDC matrix
rownames(GGDC) <- paste(samples, ".", sep = "") # If I don't add something to the end it doesn't understand rows and columns
matrix.og <- GGDC
mybin <- c(1,4) # Important so that names and format are rendered correctly later
names.pers <- data.frame(pos = mybin, name = c("23IB32_concoct.14", "CO_concoct.133")) # New names for my bins

## ---- ANI ----

ANI <- read.table("Indices_genomicos/CO_concoct.133/Resultados/Resultados_orthoANI.txt", header = T, row.names = 1)
names <- read.table("Indices_genomicos/CO_concoct.133/Resultados/Resultados_orthoANI_filenames.txt", header = TRUE)
samples <- unlist(strsplit(names$file, split = ".fsa", fixed = T)) # Removing .fsa from names
colnames(ANI) <- samples
rownames(ANI) <- paste(samples, ".", sep = "")
matrix.og <- ANI
mybin <- c(1,13) # Important so that names and format are rendered correctly later
names.pers <- data.frame(pos = mybin, name = c("CO_concoct.133", "23IB32_concoct.14")) # New names for my bins

## ---- AAI ----
AAI <- read.table("Indices_genomicos/CO_concoct.133/Resultados/Resultados_AAI.txt", header = T)
samples <- str_split_fixed(colnames(AAI), "[.]f", n = 2)[,1] # Getting the names and removing .fsa / .fa
colnames(AAI) <- samples
rownames(AAI) <- paste(samples, ".", sep = "")
matrix.og <- AAI
mybin <- c(1,4) # Important so that names and format are rendered correctly later
names.pers <- data.frame(pos = mybin, name = c("23IB32_concoct.14", "CO_concoct.133")) # New names for my bins

## ---- COMMON FOR ALL OF THEM ----

# Reordering the data frame
my.order <- c(4,1,16,3,2,13,17,5,7,14,10,11,12,9,8,6,15,18) # GGDC & AAI
# my.order <- c(1,13,12,17,16,9,15,18,14,10,6,7,8,5,4,3,11,2) # ANI
matrix.og <- matrix.og[my.order,my.order]
samples <- samples[my.order]
mybin <- match(mybin, my.order) # New position of my bin according to the new order
names.pers$pos <- match(names.pers$pos, my.order) # New position of my bin according to the new order

# Adding a value out of the limits so that half of the matrix isn't plotted.
matrix.og[lower.tri(matrix.og)] <- 101 # DDH & AAI
# matrix.og[upper.tri(matrix.og)] <- 101 # ANI

# Transforming to long format for the plot
matrix_superfinal <- as.data.frame(matrix.og) %>% 
                     rownames_to_column(var = "species1") %>% 
                     gather(species2, value, 2:(length(colnames(matrix.og))+1))

# Converting samples to factors so that the order is respected in the plot
# and making sure numbers are numeric (depends on R version):
matrix_superfinal$value <- as.numeric(matrix_superfinal$value)
matrix_superfinal$species1 <- factor(matrix_superfinal$species1, levels = paste(samples, ".", sep = ""))
matrix_superfinal$species2 <- factor(matrix_superfinal$species2, levels = samples)

# ------------------------
# ---- MY MAGNUM OPUS ----
# ------------------------

# DISPLAYING NAMES WITH CORRECT FORMAT IN THE PLOT

names <- list() # Creating a list that will contain each bin

# Under each bin, a list with its genus, species and strain

for (bin in 1:length(samples)) {
   names[[bin]] <- list(genus = NULL, species = NULL, strain = NULL)
}

# Getting each part of the name from the "samples"
# (that we grabbed from the names of the original matrix)

for (bin in 1:length(samples)) {
   if (bin %in% mybin) { # This is key: bins in "mybin" have a customised name so they aren't split
      names[[bin]]$strain <- names.pers[names.pers$pos == bin, "name"]
   } else {
      names[[bin]]$genus <- unlist(strsplit(samples[bin], split = "[_]"))[1]
      names[[bin]]$species <- unlist(strsplit(samples[bin], split = "[_]"))[2]
      names[[bin]]$strain <- unlist(strsplit(samples[bin], split = "[_]"))[3]
   }
  names[[bin]][is.na(names[[bin]])] <- NULL # Remove NAs
}

# Now we'll create a list with the text that we need to write each name
# In order to have the genus and species in italics, but not the strain,
# we need an expression like:
# expression(paste(italic("Genero "), italic("especie"), "strain"))
# We can also add the superindex T for type strains if we write:
# expression(paste(italic("Genero "), italic("especie"), strain^{'T'}))

names.graph <- list()

for (bin in 1:length(names)) {
  if (bin %in% mybin) { # My bin is not written in italics, but it is bold
    names.graph[bin] <- paste("expression(paste(bold(\"", names[[bin]]$strain, "\")))", sep = "")
      # Since we need to write quotes inside the expression we use \ so it ignores them
  }
  else { 
    if (length(grep("01", names[[bin]]$genus)) != 0) { # Names for alphanumeric genera are not in italics
      names.graph[bin] <- paste("expression(paste(\"", names[[bin]]$genus, " \", \"", 
                                  names[[bin]]$species, " \", \"", names[[bin]]$strain, "\"))", sep = "")
    } else {
    # In the rest of cases, genus and species are in italics and T is added
    names.graph[bin] <- paste("expression(paste(italic(\"", names[[bin]]$genus, 
                              " \"), italic(\"", names[[bin]]$species, 
                              " \"), ", names[[bin]]$strain, "^{'T'}))", sep = "") # No quotes for the strain if we add superindex
    }
  }
}

names.graph <- paste(names.graph, collapse = ", ") # We collapse the list into a single string
names.graph <- paste("c(", names.graph, ")", sep = "") # We add c() so it will be interpreted as a vector

names.graph.final <- eval(parse(text = names.graph)) # Evaluating the expression and saving it for the plot

# --------------
# ---- PLOT ----
# --------------

ggplot(matrix_superfinal, aes(x = species1, y = species2, fill = value))  + 
  geom_tile(color="white") + # Plot heatmap
  coord_fixed()  +   # Squares
  scale_fill_gradientn(colours = c(paletteer_d("ggsci::cyan_material"),"white"),
                       values = scales::rescale(c(10, 20, 25, 30, 35, 50, 70, 80, 90, 99, 100))
                       )+
  # scale_fill_gradientn(colours = c(paletteer_d("ggsci::orange_material"),"white"),
  #                       values = scales::rescale(c(60, 65, 70, 73, 75, 80, 85, 90, 95, 99, 100))
  #                       )+ # ANI
  # scale_fill_gradientn(colours = c(paletteer_d("ggsci::cyan_material"),"white"),
  #                      values = scales::rescale(c(60, 62, 64, 66, 68, 70, 80, 85, 90, 99,99.5))
  #                      )+ # AAI
  scale_y_discrete(limits = rev(levels(matrix_superfinal$species2)), labels = rev(names.graph.final)) + # Italics ARE WORKING
  scale_x_discrete(labels = names.graph.final) + 
  xlab(label = NULL) +
  ylab(label = NULL) +
  labs(fill = "DDH (%)") +
  # labs(fill = "ANI (%)") +
  # labs(fill = "AAI (%)") +
  theme(
    text = element_text(family = "Poppins"),
    legend.position = "right",  # Maybe you want it left?
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.text = element_text(size = 8),
    legend.title = element_text(size=10, face="bold", margin = margin(b = 10)),
    axis.title.x = element_text(size=12, face="bold",  
                                margin = margin(t = 0, r = 0, b = 20, l = 0)), # Sets a little spaces between the text and the coord.
    axis.title.y = element_text(size=12, face="bold",
                                margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.text.x = element_text(angle=45, hjust=1, size = 7, color = "black"),  # It's angled 45 degress. It can be set to '0'.
    axis.text.y = element_text(angle=0, hjust=1, size = 7, color = "black"),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "transparent", colour = NA),  # No background
    plot.background = element_rect(fill = "transparent", colour = NA),   # No background
  ) +
geom_text(aes(label= ifelse(round(value) < 70, round(value, digits = 0), "" )), # The cutoff for the text color was 90 for ANI & AAI
          size = 2.5, colour = "black", family = "Poppins") +                 # (depends on how dark the matrix is)
geom_text(aes(label= ifelse(round(value) >=70, round(value, digits = 0), "" )), 
          size = 2.5, colour = "white", family = "Poppins") 
    # Values above 100 are not plotted. The upper matrix has values above 100 (101), so those numbers won't be labelled.  

ggsave("DDH_CO_concoct.133.png", width = 20, height = 15, units = "cm", bg = "white", dpi = 600)
ggsave("DDH_CO_concoct.133.svg", width = 20, height = 15, units = "cm")
# ggsave("ANI_CO_concoct.133.png", width = 20, height = 15, units = "cm", bg = "white", dpi = 600)
# ggsave("ANI_CO_concoct.133.svg", width = 20, height = 15, units = "cm")
# ggsave("AAI_CO_concoct.133.png", width = 20, height = 15, units = "cm", bg = "white", dpi = 600)
# ggsave("AAI_CO_concoct.133.svg", width = 20, height = 15, units = "cm")
