# --------------------------------------------------------------------

# Script that creates an "upset plot", representing the intersections
# and crossovers between families of osmorregulation genes present.
# Some functions are modified from the "UpSetR" package.

# Uses the osmorregulation gene table obtained with script 17.

# R 4.2.1

# Author: Eva Sanmartín Vázquez

# --------------------------------------------------------------------

# Load libraries

library('UpSetR')
library('tidyr')
library('tibble')
library('dplyr')
library('stringr')
library('ggplot2')
library('extrafont')
library('paletteer')

# Load data

table <- read.table("osmorregulation_table.tsv", sep = "\t", header = T, 
                    row.names = 1, check.names = F)

tax <- read.table("taxonomia_todo.tsv", header = T, sep = "\t", as.is = T)

# -------------------
# PREPARING THE DATA
# -------------------

# Add taxonomy info (family level) to the osmorregulation gene data frame

family <- c()
bins <- c()

for (bin in rownames(table)) {
   sample <- unlist(strsplit(bin, split = "_", fixed = T))[1]
   bin.name <- unlist(strsplit(bin, split = "_", fixed = T))[2]
   bin.rows <- grep(paste(bin.name, ".", sep = ""), tax$bin.names, fixed = T)
   for (row in bin.rows) {
      if (tax[row, "sample"] == sample) {
         family <- c(family, tax[row, "family"])
         bins <- c(bins, tax[row, "bin.names"])
      }
   }
}

table2 <- cbind(as.data.frame(table), family)

# Remove the genes that are not present in any of the bins.
# For genes that are present, transform all values to 1 in 
# order to evaluate only presence or absence.

for (gene in colnames(table2)) {
   if (sum(table2[,gene] == rep(0, length(table2[,gene]))) == length(table2[,gene])) {
      table2 <- table2[,-which(colnames(table2) == gene)]
   } else {
      if (gene != "family") {
         table2[,gene][table2[,gene] != 0] <- 1
      }
   }  
}

# Transform into long format

df.long2 <- as.data.frame(table2) %>% 
           rownames_to_column(var = "Bin") %>% 
           gather(KO, Present, 2:(length(colnames(table2))))

# Create a table with the count of genes present in each family

freqs <- count(df.long2[df.long2[,"Present"] == 1,], family, KO)


# The input to the upset plot is a list of vectors, each one
# is named after the set (family) and contains all the elements
# of that set (genes).

data <- list()

for (fam in levels(as.factor(freqs$family))) {
   fam.good <- unlist(strsplit(fam, "__", fixed = T))[2]
   data[[fam.good]] <- freqs[freqs[,"family"] == fam, "KO"]
}

# Expression for the labels of the plot in italics:

names.graph <- list()

names.ordered <- names(data)[rev(c(4, 3, 5, 6, 1, 2))]

for (i in 1:length(names.ordered)) {
  if (length(grep("01", names.ordered[i])) != 0) {
    names.graph[(i)] <- paste("expression(paste(\"", names.ordered[i], "\"))", sep = "")
  }
  else {
    names.graph[(i)] <- paste("expression(paste(italic(\"", names.ordered[i], "\")))", sep = "")
  }
}

names.graph <- paste(names.graph, collapse = ", ") # Collapsing the list into a single string
names.graph <- paste("c(", names.graph, ")", sep = "") # Adding the "vector"

names.graph.final <- eval(parse(text = names.graph)) # Evaluate and save in a variable for the plot

# ----------
# PLOT
# ----------

#### MODIFIED FUNCTIONS FROM THE PACKAGE UPSETR ####

Make_main_bar <- function(Main_bar_data, Q, show_num, ratios, customQ, number_angles, number.colors,
                          ebar, ylabel, ymax, scale_intersections, text_scale, attribute_plots, plot.title){

  bottom_margin <- (-1)*0.65

  if(is.null(attribute_plots) == FALSE){
    bottom_margin <- (-1)*0.45
  }
  
  if(length(text_scale) > 1 && length(text_scale) <= 6){
    y_axis_title_scale <- text_scale[1]
    y_axis_tick_label_scale <- text_scale[2]
    intersection_size_number_scale <- text_scale[6]
  }
  else{
    y_axis_title_scale <- text_scale
    y_axis_tick_label_scale <- text_scale
    intersection_size_number_scale <- text_scale
  }
  
  if(is.null(Q) == F){
    inter_data <- Q
    if(nrow(inter_data) != 0){
      inter_data <- inter_data[order(inter_data$x), ]
    }
    else{inter_data <- NULL}
  }
  else{inter_data <- NULL}
  
  if(is.null(ebar) == F){
    elem_data <- ebar
    if(nrow(elem_data) != 0){
      elem_data <- elem_data[order(elem_data$x), ]
    }
    else{elem_data <- NULL}
  }
  else{elem_data <- NULL}
  
  #ten_perc creates appropriate space above highest bar so number doesnt get cut off
  if(is.null(ymax) == T){
  ten_perc <- ((max(Main_bar_data$freq)) * 0.1)
  ymax <- max(Main_bar_data$freq) + ten_perc
  }
  
  if(ylabel == "Intersection Size" && scale_intersections != "identity"){
    ylabel <- paste("Intersection Size", paste0("( ", scale_intersections, " )"))
  }
  if(scale_intersections == "log2"){
    Main_bar_data$freq <- round(log2(Main_bar_data$freq), 2)
    ymax <- log2(ymax)
  }
  if(scale_intersections == "log10"){
    Main_bar_data$freq <- round(log10(Main_bar_data$freq), 2)
    ymax <- log10(ymax)
  }
  Main_bar_plot <- (ggplot(data = Main_bar_data, aes_string(x = "x", y = "freq")) 
                    + scale_y_continuous(trans = scale_intersections)
                    + ylim(0, ymax)
                    + geom_bar(stat = "identity", width = 0.6,
                               fill = Main_bar_data$color)
                    + scale_x_continuous(limits = c(0,(nrow(Main_bar_data)+1 )), expand = c(0,0),
                                         breaks = NULL)
                    + xlab(NULL) + ylab(ylabel) +labs(title = NULL)
                    + theme(text = element_text(family = "Poppins"), 
                            panel.background = element_rect(fill = "white"),
                            plot.margin = unit(c(0.5,0.5,bottom_margin,0.5), "lines"), panel.border = element_blank(),
                            axis.title.y = element_text(vjust = -0.8, size = 8.3*y_axis_title_scale), axis.text.y = element_text(vjust=0.3,
                                                                                                            size=7*y_axis_tick_label_scale)))
  if((show_num == "yes") || (show_num == "Yes")){
    if(is.null(number.colors)) {
      Main_bar_plot <- (Main_bar_plot + geom_text(aes_string(label = "freq"), size = 2.2*intersection_size_number_scale, vjust = -1,
                                                  angle = number_angles, colour = Main_bar_data$color,  family = "Poppins"))
    } else {
      Main_bar_plot <- (Main_bar_plot + geom_text(aes_string(label = "freq"), size = 2.2*intersection_size_number_scale, vjust = -1,
                                                  angle = number_angles, colour = number.colors, family = "Poppins"))
    }
  }

  bInterDat <- NULL
  pInterDat <- NULL
  bCustomDat <- NULL
  pCustomDat <- NULL
  bElemDat <- NULL
  pElemDat <- NULL
  if(is.null(elem_data) == F){
    bElemDat <- elem_data[which(elem_data$act == T), ]
    bElemDat <- bElemDat[order(bElemDat$x), ]
    pElemDat <- elem_data[which(elem_data$act == F), ]
  }
  if(is.null(inter_data) == F){
    bInterDat <- inter_data[which(inter_data$act == T), ]
    bInterDat <- bInterDat[order(bInterDat$x), ]
    pInterDat <- inter_data[which(inter_data$act == F), ]
  }
  if(length(customQ) != 0){
    pCustomDat <- customQ[which(customQ$act == F), ]
    bCustomDat <- customQ[which(customQ$act == T), ]
    bCustomDat <- bCustomDat[order(bCustomDat$x), ]
  }
  if(length(bInterDat) != 0){
    Main_bar_plot <- Main_bar_plot + geom_bar(data = bInterDat,
                                              aes_string(x="x", y = "freq"),
                                              fill = bInterDat$color,
                                              stat = "identity", position = "identity", width = 0.6)
  }
  if(length(bElemDat) != 0){
    Main_bar_plot <- Main_bar_plot + geom_bar(data = bElemDat,
                                              aes_string(x="x", y = "freq"),
                                              fill = bElemDat$color,
                                              stat = "identity", position = "identity", width = 0.6)
  }
  if(length(bCustomDat) != 0){
    Main_bar_plot <- (Main_bar_plot + geom_bar(data = bCustomDat, aes_string(x="x", y = "freq2"),
                                               fill = bCustomDat$color2,
                                               stat = "identity", position ="identity", width = 0.6))
  }
  if(length(pCustomDat) != 0){
    Main_bar_plot <- (Main_bar_plot + geom_point(data = pCustomDat, aes_string(x="x", y = "freq2"), colour = pCustomDat$color2,
                                                 size = 2, shape = 17, position = position_jitter(width = 0.2, height = 0.2)))
  }
  if(length(pInterDat) != 0){
    Main_bar_plot <- (Main_bar_plot + geom_point(data = pInterDat, aes_string(x="x", y = "freq"),
                                                 position = position_jitter(width = 0.2, height = 0.2),
                                                 colour = pInterDat$color, size = 2, shape = 17))
  }
  if(length(pElemDat) != 0){
    Main_bar_plot <- (Main_bar_plot + geom_point(data = pElemDat, aes_string(x="x", y = "freq"),
                                                 position = position_jitter(width = 0.2, height = 0.2),
                                                 colour = pElemDat$color, size = 2, shape = 17))
  }
  
  Main_bar_plot <- (Main_bar_plot 
                    + geom_vline(xintercept = 0, color = "gray0")
                    + geom_hline(yintercept = 0, color = "gray0"))
  
  if(!is.na(plot.title)) {
    Main_bar_plot <- c(Main_bar_plot + ggtitle(plot.title))
  }
  Main_bar_plot <- ggplotGrob(Main_bar_plot)
  return(Main_bar_plot)
}

Make_size_plot <- function(Set_size_data, sbar_color, ratios, ylabel, scale_sets, text_scale, set_size_angle, set_size.show, set_size.scale_max,
                           set_size.number_size){
  
  if(length(text_scale) > 1 && length(text_scale) <= 6){
    x_axis_title_scale <- text_scale[3]
    x_axis_tick_label_scale <- text_scale[4]
  }
  else{
    x_axis_title_scale <- text_scale
    x_axis_tick_label_scale <- text_scale
  }
  
  if(ylabel == "Set Size" && scale_sets != "identity"){
    ylabel <- paste("Set Size", paste0("( ", scale_sets, " )"))
    if(scale_sets == "log2"){
      Set_size_data$y <- log2(Set_size_data$y)
    }
    if(scale_sets == "log10"){
      Set_size_data$y <- log10(Set_size_data$y)
    }
  }
  
  if(!is.null(set_size.number_size)) {
    num.size <- (set_size.number_size/2.845276)*x_axis_tick_label_scale
  } else {
    num.size <- (7/2.845276)*x_axis_tick_label_scale
  }
  
  Size_plot <- (ggplot(data = Set_size_data, aes_string(x ="x", y = "y"))
                + geom_bar(stat = "identity",colour = sbar_color, width = 0.4,
                           fill = sbar_color, position = "identity")
                + scale_x_continuous(limits = c(0.5, (nrow(Set_size_data) + 0.5)),
                                     breaks = c(0, max(Set_size_data)),
                                     expand = c(0,0))
                + theme(text = element_text(family = "Poppins"),
                        panel.background = element_rect(fill = "white"),
                        plot.margin=unit(c(-0.11,-1.3,0.5,0.5), "lines"),
                        axis.title.x = element_text(size = 8.3*x_axis_title_scale),
                        axis.text.x = element_text(size = 7*x_axis_tick_label_scale,
                                                   vjust = 1, hjust = 0.5),
                        axis.line = element_line(colour = "gray0"),
                        axis.line.y = element_blank(),
                        axis.line.x = element_line(colour = "gray0", size = 0.3),
                        axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.grid.major = element_blank())
                + xlab(NULL) + ylab(ylabel)
                + coord_flip())
  
  if(set_size.show == TRUE){
    Size_plot <- (Size_plot + geom_text(aes(label=y,vjust=0.5,hjust=1.2, angle = set_size_angle), size=num.size))
  }
    
  if(scale_sets == "log10"){
    if(!is.null(set_size.scale_max)) {
      Size_plot <- (Size_plot + scale_y_continuous(limits = c(set_size.scale_max, 0), 
                                                   trans = log10_reverse_trans()))
    } else {
      Size_plot <- (Size_plot + scale_y_continuous(trans = log10_reverse_trans()))
    }
  }
  else if (scale_sets == "log2"){
    if(!is.null(set_size.scale_max)) {
      Size_plot <- (Size_plot + scale_y_continuous(limits = c(set_size.scale_max, 0), 
                                                   trans = log2_reverse_trans()))
    } else {
      Size_plot <- (Size_plot + scale_y_continuous(trans = log2_reverse_trans()))
    }
  }
  else{
    if(!is.null(set_size.scale_max)) {
      Size_plot <- (Size_plot + scale_y_continuous(limits = c(set_size.scale_max, 0),
                                                   trans = "reverse"))
    } else {
      Size_plot <- (Size_plot + scale_y_continuous(trans = "reverse"))
    }
  }
  
  Size_plot <- ggplot_gtable(ggplot_build(Size_plot))
  return(Size_plot)
}

Make_matrix_plot <- function(Mat_data,Set_size_data, Main_bar_data, point_size, line_size, text_scale, labels,
                             shading_data, shade_alpha){

                              
  if(length(text_scale) == 1){
    name_size_scale <- text_scale
  }


  if(length(text_scale) > 1 && length(text_scale) <= 6){
    name_size_scale <- text_scale[5]
  }


  Matrix_plot <- (ggplot()
                  + theme(text = element_text(family = "Poppins"),
                          panel.background = element_rect(fill = "white"),
                          plot.margin=unit(c(-0.2,0.5,0.5,0.5), "lines"),
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.text.y = element_text(colour = "gray0",
                                                     size = 7*name_size_scale, hjust = 0.4),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())
                  + xlab(NULL) + ylab("   ")
                  + scale_y_continuous(breaks = c(1:nrow(Set_size_data)),
                                       limits = c(0.5,(nrow(Set_size_data) +0.5)),
                                       labels = names.graph.final, expand = c(0,0))
                  + scale_x_continuous(limits = c(0,(nrow(Main_bar_data)+1 )), expand = c(0,0))
                  + geom_rect(data = shading_data, aes_string(xmin = "min", xmax = "max",
                                                              ymin = "y_min", ymax = "y_max"),
                              fill = shading_data$shade_color, alpha = shade_alpha)
                  + geom_point(data= Mat_data, aes_string(x= "x", y= "y"), colour = Mat_data$color,
                               size= point_size, alpha = Mat_data$alpha, shape=16)
                  + geom_line(data= Mat_data, aes_string(group = "Intersection", x="x", y="y",
                                                         colour = "color"), size = line_size)
                  + scale_color_identity())


  Matrix_plot <- ggplot_gtable(ggplot_build(Matrix_plot))


  return(Matrix_plot)


}

upset <- function(data, nsets = 5, nintersects = 40, sets = NULL, keep.order = F, set.metadata = NULL, intersections = NULL,
                  matrix.color = "gray23", main.bar.color = "gray23", mainbar.y.label = "Intersection Size", mainbar.y.max = NULL,
                  sets.bar.color = "gray23", plot.title = NA, sets.x.label = "Set Size", point.size = 2.2, line.size = 0.7,
                  mb.ratio = c(0.70,0.30), expression = NULL, att.pos = NULL, att.color = main.bar.color, order.by = c("freq", "degree"),
                  decreasing = c(T, F), show.numbers = "yes", number.angles = 0, number.colors=NULL, group.by = "degree",cutoff = NULL,
                  queries = NULL, query.legend = "none", shade.color = "gray88", shade.alpha = 0.25, matrix.dot.alpha =0.5,
                  empty.intersections = NULL, color.pal = 1, boxplot.summary = NULL, attribute.plots = NULL, scale.intersections = "identity",
                  scale.sets = "identity", text.scale = 1, set_size.angles = 0 , set_size.show = FALSE, set_size.numbers_size = NULL, set_size.scale_max = NULL){

  startend <-UpSetR:::FindStartEnd(data)
  first.col <- startend[1]
  last.col <- startend[2]

  if(color.pal == 1){
    palette <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2",
                 "#7F7F7F", "#BCBD22", "#17BECF")
  }
  else{
    palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
                 "#CC79A7")
  }

  if(is.null(intersections) == F){
    Set_names <- unique((unlist(intersections)))
    Sets_to_remove <- UpSetR:::Remove(data, first.col, last.col, Set_names)
    New_data <- UpSetR:::Wanted(data, Sets_to_remove)
    Num_of_set <- UpSetR:::Number_of_sets(Set_names)
    if(keep.order == F){
      Set_names <- UpSetR:::order_sets(New_data, Set_names)
    }
    All_Freqs <- UpSetR:::specific_intersections(data, first.col, last.col, intersections, order.by, group.by, decreasing,
                                        cutoff, main.bar.color, Set_names)
  }
  else if(is.null(intersections) == T){
    Set_names <- sets
    if(is.null(Set_names) == T || length(Set_names) == 0 ){
      Set_names <- UpSetR:::FindMostFreq(data, first.col, last.col, nsets)
    }
    Sets_to_remove <- UpSetR:::Remove(data, first.col, last.col, Set_names)
    New_data <- UpSetR:::Wanted(data, Sets_to_remove)
    Num_of_set <- UpSetR:::Number_of_sets(Set_names)
    if(keep.order == F){
    Set_names <- UpSetR:::order_sets(New_data, Set_names)
    }
    All_Freqs <- UpSetR:::Counter(New_data, Num_of_set, first.col, Set_names, nintersects, main.bar.color,
                         order.by, group.by, cutoff, empty.intersections, decreasing)
  }
  Matrix_setup <- UpSetR:::Create_matrix(All_Freqs)
  labels <- UpSetR:::Make_labels(Matrix_setup)
  #Chose NA to represent NULL case as result of NA being inserted when at least one contained both x and y
  #i.e. if one custom plot had both x and y, and others had only x, the y's for the other plots were NA
  #if I decided to make the NULL case (all x and no y, or vice versa), there would have been alot more if/else statements
  #NA can be indexed so that we still get the non NA y aesthetics on correct plot. NULL cant be indexed.
  att.x <- c(); att.y <- c();
  if(is.null(attribute.plots) == F){
    for(i in seq_along(attribute.plots$plots)){
      if(length(attribute.plots$plots[[i]]$x) != 0){
        att.x[i] <- attribute.plots$plots[[i]]$x
      }
      else if(length(attribute.plots$plots[[i]]$x) == 0){
        att.x[i] <- NA
      }
      if(length(attribute.plots$plots[[i]]$y) != 0){
        att.y[i] <- attribute.plots$plots[[i]]$y
      }
      else if(length(attribute.plots$plots[[i]]$y) == 0){
        att.y[i] <- NA
      }
    }
  }

  BoxPlots <- NULL
  if(is.null(boxplot.summary) == F){
    BoxData <- UpSetR:::IntersectionBoxPlot(All_Freqs, New_data, first.col, Set_names)
    BoxPlots <- list()
    for(i in seq_along(boxplot.summary)){
      BoxPlots[[i]] <- UpSetR:::BoxPlotsPlot(BoxData, boxplot.summary[i], att.color)
    }
  }

  customAttDat <- NULL
  customQBar <- NULL
  Intersection <- NULL
  Element <- NULL
  legend <- NULL
  EBar_data <- NULL
  if(is.null(queries) == F){
    custom.queries <- UpSetR:::SeperateQueries(queries, 2, palette)
    customDat <- UpSetR:::customQueries(New_data, custom.queries, Set_names)
    legend <- UpSetR:::GuideGenerator(queries, palette)
    legend <- UpSetR:::Make_legend(legend)
    if(is.null(att.x) == F && is.null(customDat) == F){
      customAttDat <- UpSetR:::CustomAttData(customDat, Set_names)
    }
    customQBar <- UpSetR:::customQueriesBar(customDat, Set_names, All_Freqs, custom.queries)
  }
  if(is.null(queries) == F){
    Intersection <- UpSetR:::SeperateQueries(queries, 1, palette)
    Matrix_col <- UpSetR:::intersects(QuerieInterData, Intersection, New_data, first.col, Num_of_set,
                             All_Freqs, expression, Set_names, palette)
    Element <- UpSetR:::SeperateQueries(queries, 1, palette)
    EBar_data <-UpSetR:::ElemBarDat(Element, New_data, first.col, expression, Set_names,palette, All_Freqs)
  }
  else{
    Matrix_col <- NULL
  }

  Matrix_layout <- UpSetR:::Create_layout(Matrix_setup, matrix.color, Matrix_col, matrix.dot.alpha)
  Set_sizes <- UpSetR:::FindSetFreqs(New_data, first.col, Num_of_set, Set_names, keep.order)
  Bar_Q <- NULL
  if(is.null(queries) == F){
    Bar_Q <- UpSetR:::intersects(QuerieInterBar, Intersection, New_data, first.col, Num_of_set, All_Freqs, expression, Set_names, palette)
  }
  QInter_att_data <- NULL
  QElem_att_data <- NULL
  if((is.null(queries) == F) & (is.null(att.x) == F)){
    QInter_att_data <- UpSetR:::intersects(QuerieInterAtt, Intersection, New_data, first.col, Num_of_set, att.x, att.y,
                                  expression, Set_names, palette)
    QElem_att_data <- UpSetR:::elements(QuerieElemAtt, Element, New_data, first.col, expression, Set_names, att.x, att.y,
                               palette)
  }
  AllQueryData <- UpSetR:::combineQueriesData(QInter_att_data, QElem_att_data, customAttDat, att.x, att.y)

  ShadingData <- NULL

  if(is.null(set.metadata) == F){
    ShadingData <- UpSetR:::get_shade_groups(set.metadata, Set_names, Matrix_layout, shade.alpha)
    output <- UpSetR:::Make_set_metadata_plot(set.metadata, Set_names)
    set.metadata.plots <- output[[1]]
    set.metadata <- output[[2]]

    if(is.null(ShadingData) == FALSE){
    shade.alpha <- unique(ShadingData$alpha)
    }
  } else {
    set.metadata.plots <- NULL
  }
  if(is.null(ShadingData) == TRUE){
  ShadingData <- UpSetR:::MakeShading(Matrix_layout, shade.color)
  }
  Main_bar <- suppressMessages(Make_main_bar(All_Freqs, Bar_Q, show.numbers, mb.ratio, customQBar, number.angles, number.colors, EBar_data, mainbar.y.label,
                            mainbar.y.max, scale.intersections, text.scale, attribute.plots, plot.title))
  Matrix <- Make_matrix_plot(Matrix_layout, Set_sizes, All_Freqs, point.size, line.size,
                             text.scale, labels, ShadingData, shade.alpha)
  Sizes <- Make_size_plot(Set_sizes, sets.bar.color, mb.ratio, sets.x.label, scale.sets, text.scale, set_size.angles,set_size.show,
                          set_size.scale_max, set_size.numbers_size)
  
  structure(class = "upset",
    .Data=list(
      Main_bar = Main_bar,
      Matrix = Matrix,
      Sizes = Sizes,
      labels = labels,
      mb.ratio = mb.ratio,
      att.x = att.x,
      att.y = att.y,
      New_data = New_data,
      expression = expression,
      att.pos = att.pos,
      first.col = first.col,
      att.color = att.color,
      AllQueryData = AllQueryData,
      attribute.plots = attribute.plots,
      legend = legend,
      query.legend = query.legend,
      BoxPlots = BoxPlots,
      Set_names = Set_names,
      set.metadata = set.metadata,
      set.metadata.plots = set.metadata.plots)
  )
}


#### PLOT ####

png("osmo_upset.png", width = 13, height = 13, units = "cm", res = 300)

upset_plot <- upset(
    fromList(data), order.by = "freq", nsets = 6,
    mb.ratio = c(0.7,0.3), # Height ratio between the bar plot and matrix part of the plot
    matrix.color = "#8cabd9",
    main.bar.color = "#8cabd9",
    number.colors = "gray23",
    sets.bar.color = c(rep("#8eddc6",2), rep("#b888e5",3), "#8eddc6"), # Different colors for archaea and bacteria
    shade.color = "#c5daf9",
    sets.x.label = "Number of genes",
    text.scale = c(1.2,1,1.1,1,1.2,1) # Vector with text sizes in the following order:
                                      # c(intersection size title, intersection size tick labels, 
                                      # set size title, set size tick labels, set names, numbers above bars)
    ) 

dev.off()

ggsave("osmo_upset_matrix.svg", plot = upset_plot$Matrix, width = 13, height = 3, units = "cm")
ggsave("osmo_upset_sizes.svg", plot = upset_plot$Sizes, width = 5, height = 5, units = "cm")

# ---- SETS ----

# What genes do each of the sets in the graph correspond to?

sets <- list()

# Haloarculaceae + Haloferacaceae

sets$HaloHalo <- setdiff(
  setdiff(
    setdiff(
      setdiff(
        intersect(data$Haloarculaceae, data$Haloferacaceae), 
        data$Salinibacteraceae), 
    data$Nitrococcaceae), 
  data$JAAYUY01),
data$Nanosalinaceae)

# Haloarculaceae + Haloferacaceae + Salinibacteraceae + Nitrococcaceae

sets$HaloHaloSaliNitro <- setdiff(
  setdiff(
    intersect(
      intersect(
        intersect(data$Haloarculaceae, data$Haloferacaceae), 
        data$Salinibacteraceae), 
    data$Nitrococcaceae), 
  data$JAAYUY01),
data$Nanosalinaceae)

# Haloarculaceae + Haloferacaceae + Nitrococcaceae

sets$HaloHaloNitro <- setdiff(
  setdiff(
    setdiff(
      intersect(
        intersect(data$Haloarculaceae, data$Haloferacaceae), 
        data$Nitrococcaceae), 
    data$Salinibacteraceae), 
  data$JAAYUY01),
data$Nanosalinaceae)

# Haloarculaceae + Haloferacaceae + Salinibacteraceae

sets$HaloHaloSali <- setdiff(
  setdiff(
    setdiff(
      intersect(
        intersect(data$Haloarculaceae, data$Haloferacaceae), 
        data$Salinibacteraceae), 
    data$Nitrococcaceae), 
  data$JAAYUY01),
data$Nanosalinaceae)

# Haloarculaceae + Haloferacaceae + Salinibacteraceae + Nitrococcaceae + JAAYUY01

sets$HaloHaloSaliNitroJA <- setdiff(
  intersect(
    intersect(
      intersect(
        intersect(data$Haloarculaceae, data$Haloferacaceae), 
        data$Salinibacteraceae), 
    data$Nitrococcaceae), 
  data$JAAYUY01),
data$Nanosalinaceae)

# Nitrococcaceae

sets$Nitro <- setdiff(
  setdiff(
    setdiff(
      setdiff(
        setdiff(data$Nitrococcaceae, data$Haloferacaceae), 
        data$Salinibacteraceae), 
    data$Haloarculaceae), 
  data$JAAYUY01),
data$Nanosalinaceae)

# Nanosalinaceae

sets$Nano <- setdiff(
  setdiff(
    setdiff(
      setdiff(
        setdiff(data$Nanosalinaceae, data$Haloferacaceae), 
        data$Salinibacteraceae), 
    data$Haloarculaceae), 
  data$JAAYUY01),
data$Nitrococcaceae)

# Salinibacteraceae

sets$Sali <- setdiff(
  setdiff(
    setdiff(
      setdiff(
        setdiff(data$Salinibacteraceae, data$Haloferacaceae), 
        data$Haloarculaceae), 
    data$Nanosalinaceae), 
  data$JAAYUY01),
data$Nitrococcaceae)

# JAAYUY01 + Salinibacteraceae

sets$JASali <- setdiff(
  setdiff(
    setdiff(
      setdiff(
        intersect(data$Salinibacteraceae, data$JAAYUY01), 
        data$Haloarculaceae), 
    data$Nanosalinaceae), 
  data$Haloferacaceae),
data$Nitrococcaceae)

# JAAYUY01 + Haloarculaceae

sets$JAHaloa <- setdiff(
  setdiff(
    setdiff(
      setdiff(
        intersect(data$Haloarculaceae, data$JAAYUY01), 
        data$Salinibacteraceae), 
    data$Nanosalinaceae), 
  data$Haloferacaceae),
data$Nitrococcaceae)

# Haloferacaceae

sets$Halof <- setdiff(
  setdiff(
    setdiff(
      setdiff(
        setdiff(data$Haloferacaceae, data$JAAYUY01), 
        data$Haloarculaceae), 
    data$Nanosalinaceae), 
  data$Salinibacteraceae),
data$Nitrococcaceae)

# All families

sets$all <- intersect(
  intersect(
    intersect(
      intersect(
        intersect(data$Haloferacaceae, data$JAAYUY01), 
        data$Haloarculaceae), 
    data$Nanosalinaceae), 
  data$Salinibacteraceae),
data$Nitrococcaceae)

# JAAYUY01

sets$JAAYUY01 <- setdiff(
  setdiff(
    setdiff(
      setdiff(
        setdiff(data$JAAYUY01, data$Haloferacaceae), 
        data$Haloarculaceae), 
    data$Nanosalinaceae), 
  data$Salinibacteraceae),
data$Nitrococcaceae)

# Nitrococcaceae + Salinibacteraceae

sets$NitroSali <- setdiff(
  setdiff(
    setdiff(
      setdiff(
        intersect(data$Nitrococcaceae, data$Salinibacteraceae), 
        data$Haloarculaceae), 
    data$Nanosalinaceae), 
  data$JAAYUY01),
data$Haloferacaceae)

# Salinibacteraceae + Haloarculaceae

sets$SaliHaloa <- setdiff(
  setdiff(
    setdiff(
      setdiff(
        intersect(data$Salinibacteraceae, data$Haloarculaceae), 
        data$JAAYUY01), 
    data$Nanosalinaceae), 
  data$Haloferacaceae),
data$Nitrococcaceae)

# Salinibacteraceae + Haloferacaceae

sets$SaliHalof <- setdiff(
  setdiff(
    setdiff(
      setdiff(
        intersect(data$Salinibacteraceae, data$Haloferacaceae), 
        data$Haloarculaceae), 
    data$Nanosalinaceae), 
  data$JAAYUY01),
data$Nitrococcaceae)

# All archaea families

sets$archaea <- setdiff(
  setdiff(
    setdiff(
      intersect(
        intersect(data$Haloarculaceae, data$Haloferacaceae), 
        data$Nanosalinaceae), 
    data$JAAYUY01), 
  data$Salinibacteraceae),
data$Nitrococcaceae)

# Haloarculaceae + Haloferacaceae + Nitrococcaceae + JAAYUY01

sets$HaloaHalofNitroJA <- setdiff(
  setdiff(
    intersect(
      intersect(
        intersect(data$Haloarculaceae, data$Haloferacaceae), 
        data$Nitrococcaceae), 
    data$JAAYUY01), 
  data$Salinibacteraceae),
data$Nanosalinaceae)

# ---- BAR PLOT ----

## --- DATASET ----

# Barplot similar to the one in the upset plot but coloring
# genes by their category

# Reading the data frame with category info
genes <- read.table("Genes/osmorregulation.tsv", sep = "\t", header = T)

# We create a data frame with the genes present in each set
# of the upset plot from the sets list and we assign those genes to 
# their categories.

sets.cat <- data.frame(set.name = NULL, gene = NULL, category = NULL)

for (set in names(sets)) {
  sets.cat.each <- data.frame(set.name = rep(set, length(sets[[set]])), 
    gene = sets[[set]], 
    # Separating the name of the gene from the KO and searching the gene table to find its category:
    category = genes[match(str_split_fixed(sets[[set]], " [(]", 2)[,1], genes$Name), "Category"])
  sets.cat <- rbind(sets.cat, sets.cat.each)
}

# Transforming the set names to factors so that the order is respected in the plot
sets.cat$set.name <- factor(sets.cat$set.name, levels = names(sets))

## --- PLOTTING ----

ggplot(sets.cat, aes(x = set.name, fill = category)) +
    geom_bar(position = "stack", width = 0.5) +
    scale_x_discrete(labels = NULL) +
    scale_y_continuous(limits = c(0,12.5), breaks = seq(0, 12, 2), expand = expansion(mult = c(0,0.02))) +
    ylab("Number of genes") +
    xlab("") +
    labs(fill = "Category") +
    scale_fill_manual(values = c("#DCB0F2", "#B497E7", "#9EB9F3",
    "#91dffe", "#F6CF71", "#F89C74", "#8BE0A4FF"), # From paletteer_d("rcartocolor::Pastel")
    breaks = c("Glycine-betaine synthesis", "Glycine-betaine degradation", 
    "Glycine-betaine transport", "Choline transport", "Ectoine synthesis", 
    "Ectoine degradation", "Salt-in")) + # Changing the order of the legend
    theme(
        text = element_text(family = "Poppins"),
        panel.background = element_blank(), # No background color
        legend.position = "right", 
        axis.title.y = element_text(colour="black", margin = margin(r = 10, t=10), 
                      size = 17, face = "bold"), # Título del eje Y
        axis.text.y = element_text(colour="black", margin = margin(r = 5, t=5),  
                      size= 15), # Números del eje Y
        axis.line = element_line(size=0.75), # Líneas de los ejes
        axis.ticks.y = element_line(size=0.75), # Rayitas de los números de los ejes
        axis.ticks.x = element_blank()
    )

ggsave("osmo_sets_barplot.png", height = 20, width = 20, units = "cm")
ggsave("osmo_sets_barplot.svg", height = 20, width = 20, units = "cm")
