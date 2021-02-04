Tree and tree to map
================

Load Data
---------

First we load the required packages and files.

Note that I use a cladogram which I manually constructed with FigTree. The coordinates for the Outgroups are just fake numbers (coordinates for a lake in Turkmenistan) so they align nicely next to our samples.

``` r
rm(list = ls())

library(phytools)
library(ggplot2)
library(mapdata)
library(viridis)
library(ggtree)
library(extrafont)
library(cowplot)

Nt_Tree = read.tree("../00_Data/Trees/RAxML_Nucleotide_Cladogram.nwk")

Coordinates = as.matrix(read.csv("../00_Data/Miscellaneous/Coordinates.csv", sep = "\t", header = T, row.names = 1))
```

Plot phylogeny projected on a map
=================================

Build PhyloToMap object
-----------------------

The xlim and ylim are roughly the border coordinates for Europe. Then we build a color vector based on the lineage (Turkey, viridis, bilineata, Adriatic, outgroup).

``` r
obj = phylo.to.map(Nt_Tree, 
                   Coordinates, 
                   xlim = c(-10, 59), 
                   ylim = c(30, 71),
                   fill = T, 
                   col = "moccasin",
                   plot=F, 
                   rotate = T, 
                   direction = "downwards", 
                   asp = 1)

colo = rep("x", length(Nt_Tree$tip.label))
for (i in seq(1:length(Nt_Tree$tip.label))){
  if (grepl("Lv_TR", Nt_Tree$tip.label[i])){
    colo[i] = viridis(4, end = 0.9)[4]
  } else if (grepl("Lv", Nt_Tree$tip.label[i])) {
    colo[i] = viridis(4, end = 0.9)[3]
  } else if (grepl("Lb", Nt_Tree$tip.label[i])) {
    colo[i] = viridis(4, end = 0.9)[1]
  } else if (grepl("AL", Nt_Tree$tip.label[i])) {
    colo[i] = viridis(4, end = 0.9)[2]
  } else {
    colo[i] = NA_character_
  }
}
cols = setNames(colo, Nt_Tree$tip.label)
```

Plot Phylo to Map
-----------------

Next we plot the object and save it as a PDF, because in the default plot window it looks quite weird.

``` r
pdf(file="TreetoMap.pdf") 
plot(obj, 
     type = "phylogram", 
     direction = "downwards", 
     colors = cols, 
     fill = T, 
     col = "moccasin", 
     cex.points = c(1.1, 1), 
     lwd = 1.2, 
     bg = "darkblue", 
     xlim = c(-10, 59), 
     ylim = c(30, 71)) 
invisible(dev.off())
```

![](Readme_files/figure-markdown_github/Plot%20PhyloToMap-1.png)

Align labels
------------

Thats already very nice, but let's try aligning the tip labels manually. So first we check which label is the longest, and for the shorter labels we paste "·" as many times as required so that in the end all labels have the same number of characters.

``` r
obj_aligned = obj

# get the longest label
MaxLabelLength = max(nchar(obj$tree$tip.label))

# now loop over the labels and append "·" 
for (i in (seq(1:length(obj$tree$tip.label)))){
  if (nchar(obj$tree$tip.label[i]) != MaxLabelLength){
    dash = paste0(rep("-", each = MaxLabelLength-nchar(obj$tree$tip.label[i])), collapse = "")
    obj_aligned$tree$tip.label[i] = paste0(dash, obj$tree$tip.label[i])
  }
}

for (i in (seq(1:length(rownames(obj$coords))))){
  if (nchar(rownames(obj$coords)[i]) != MaxLabelLength){
    dash = paste0(rep("-", each = MaxLabelLength-nchar(rownames(obj$coords)[i])), collapse = "")
    rownames(obj_aligned$coords)[i] = paste0(dash, rownames(obj$coords)[i])
  }
}

#  now assign the colors to the new names
colo_aligned = rep("x", length(obj_aligned$tree$tip.label))
for (i in seq(1:length(obj_aligned$tree$tip.label))){
  if (grepl("Lv_TR", obj_aligned$tree$tip.label[i])){
    colo_aligned[i] = viridis(4, end = 0.9)[4]
  } else if (grepl("Lv", obj_aligned$tree$tip.label[i])) {
    colo_aligned[i] = viridis(4, end = 0.9)[3]
  } else if (grepl("Lb", obj_aligned$tree$tip.label[i])) {
    colo_aligned[i] = viridis(4, end = 0.9)[1]
  } else if (grepl("AL", obj_aligned$tree$tip.label[i])) {
    colo_aligned[i] = viridis(4, end = 0.9)[2]
  } else {
    colo_aligned[i] = NA_character_
  }
}
cols_aligned = setNames(colo_aligned, obj_aligned$tree$tip.label)
```

Plot with aligned labels
------------------------

``` r
font_import(prompt = F, pattern = "DejaVuSans")
loadfonts()
#pdf(file="TreetoMap_alignedLabels.pdf", family = "DejaVu Sans Mono")
png(file="TreetoMap_alignedLabels.png", family = "DejaVu Sans Mono", res = 300, width = 18, height = 16, units = "cm")
plot(obj_aligned, 
     type = "phylogram", 
     direction = "downwards", 
     colors = cols_aligned, 
     ftype = "b", 
     cex.points = c(1.2, 1.3), 
     cex.text = 1.2,
     lwd = 3)
par(family = "DejaVu Sans")
legend(43.5, 63, legend = c(expression(italic("bilineata")), "Adriatic",
                          expression(italic("viridis")), 
                          expression(italic("viridis")*" (Turkey)"), 
                          "Outgroups"), 
       pt.bg = c(viridis(4, end = 0.9), "white"),
       col = "black", 
       title = expression(bold(bolditalic("Lacerta")*" clades")),
       bg = "grey85", pch = 21)
invisible(dev.off())
```

![](Readme_files/figure-markdown_github/Plot%20with%20aligned%20labels-1.png)

Plot phylogram with `ggtree`
============================

Load tree
---------

Now we use `ggplot` and `ggtree` to plot just the tree, with bootstrap support given by colored node points.

``` r
Nt_Tree_phylo = read.tree("../00_Data/Trees/RAxML_Nucleotide_Phylogram.nwk")

# sort tip labels so they are in the same order as in the phylo.to.map object
Nt_Tree_phylo = ape::rotateConstr(Nt_Tree_phylo, rev(obj$tree$tip.label))


colo[is.na(colo)] = "black"

g = ggtree(Nt_Tree_phylo, size = 1.1) + 
  geom_tiplab(align = T, hjust = 1, linetype = "dotted", 
              fontface = "bold", offset = 0.2, 
              geom = "label", color = colo, fill = "white", 
              family = "DejaVu Sans Mono", size = 4.5) + 
  geom_nodepoint(aes(color = as.numeric(label)), size = 4.5) + 
  scale_color_viridis_c(option = "cividis", limits = c(0, 100), 
                        direction = -1, name = "Bootstrap Support", 
                        na.value = NA) +
  geom_cladelabel(node = 48, extend = 1, 
                  label = expression(bold("Outgroups")), align = T, 
                  geom = "label", fill = "white", 
                  offset = 0.3, barsize = 2, 
                  family = "DejaVu Sans") +
  geom_cladelabel(node = 30, extend = 0.25, 
                  label = expression(bold(atop(paste(bolditalic("Lacerta viridis")), paste("(Turkey)")))), 
                  align = T, 
                  geom = "label", fill = viridis(4, end = 0.9)[4], 
                  color = c(viridis(4, end = 0.9)[4], "white"), 
                  offset = 0.3, parse = T, barsize = 2, 
                  family = "DejaVu Sans") + 
  geom_cladelabel(node = 32, extend = 0.25, 
                  label = 'bolditalic("Lacerta viridis")', align = T, 
                  geom = "label", fill = viridis(4, end = 0.9)[3], 
                  color = c(viridis(4, end = 0.9)[3], "white"), 
                  offset = 0.3, parse = T, barsize = 2, 
                  family = "DejaVu Sans") + 
  geom_cladelabel(node = 12, extend = 1.5, 
                  label = expression(bold("Adriatic clade")), align = T, 
                  geom = "label", fill = viridis(4, end = 0.9)[2], 
                  color = c(viridis(4, end = 0.9)[2], "white"), 
                  offset = 0.3, barsize = 2, 
                  family = "DejaVu Sans") + 
  geom_cladelabel(node = 44, extend = 1, 
                  label = 'bolditalic("Lacerta bilineata")', align = T, 
                  geom = "label", fill = viridis(4, end = 0.9)[1], 
                  color = c(viridis(4, end = 0.9)[1], "white"), 
                  offset = 0.3, parse = T, barsize = 2, 
                  family = "DejaVu Sans") + 
  xlim(c(-0.05, 1.25)) +
  scale_y_reverse()  + 
  theme_tree() +
  theme(legend.position = c(0.2, 0.4), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(face = "bold", size = 18, 
                  family = "DejaVu Sans")) +
  geom_treescale(x = 0, y = -7.5, width = 0.1) + 
  annotate("text", x = 0, y = 23.2, 
           label = expression(bold("Note on sample labels\n(excl. Outgroups):")), 
           parse = T, hjust = 0, family = "DejaVu Sans", size = 5) + 
  annotate("label", x = 0, y = 24.6, 
           label = expression(paste(atop("9004" , 1), atop("_Lv" , 2), atop("_TR" , 3))), 
           parse = T, hjust = 0, family = "DejaVu Sans Mono", size = 5) + 
  annotate("text", x = 0, y = 27.5, 
           label = "1: Genbank/Tissuecollection ID\n2: Clade (Lb bilineata, Lv viridis, AL adriatic)\n3: ISO Country Code", 
           hjust = 0, family = "DejaVu Sans", size = 5)


ggsave(plot = g, filename = "Tree.pdf", device = cairo_pdf, 
       width = 8.5, height = 8.5, dpi = 300, units = "in")
ggsave(plot = g, filename = "Tree.png", device = "png", 
       width = 8.5, height = 8.5, dpi = 300)

g
```

![](Readme_files/figure-markdown_github/Plot%20regular%20tree-1.png)

Plot Bayes Tree
---------------

Now we do the same for the tree calculated with MrBayes.

``` r
Bayes_Tree_phylo = read.tree("../00_Data/Trees/MrBayes_Nucleotide_Phylogram.nwk")
Bayes_Tree_clado = read.tree("../00_Data/Trees/MrBayes_Nucleotide_Cladogram.nwk")

# sort tip labels so they are in the same order as in the phylo.to.map object
Bayes_Tree_phylo = ape::rotateConstr(Bayes_Tree_phylo,
                                     rev(obj$tree$tip.label))
Bayes_Tree_clado = ape::rotateConstr(Bayes_Tree_clado,
                                     rev(obj$tree$tip.label))

colo_Bayes = rep("x", length(Bayes_Tree_phylo$tip.label))
for (i in seq(1:length(Bayes_Tree_phylo$tip.label))){
  if (grepl("Lv_TR", Bayes_Tree_phylo$tip.label[i])){
    colo_Bayes[i] = viridis(4, end = 0.9)[4]
  } else if (grepl("Lv", Bayes_Tree_phylo$tip.label[i])) {
    colo_Bayes[i] = viridis(4, end = 0.9)[3]
  } else if (grepl("Lb", Bayes_Tree_phylo$tip.label[i])) {
    colo_Bayes[i] = viridis(4, end = 0.9)[1]
  } else if (grepl("AL", Bayes_Tree_phylo$tip.label[i])) {
    colo_Bayes[i] = viridis(4, end = 0.9)[2]
  } else {
    colo_Bayes[i] = NA_character_
  }
}
colo_Bayes[is.na(colo_Bayes)] = "black"

g_Bayes = ggtree(Bayes_Tree_phylo, size = 1.1) + 
  geom_tiplab(align = T, hjust = 1, linetype = "dotted", 
              fontface = "bold", offset = 0.2, 
              geom = "label", color = colo_Bayes, fill = "white", 
              family = "DejaVu Sans Mono", size = 4.5) + 
  geom_nodepoint(aes(color = as.numeric(label)), size = 4.5) + 
  scale_color_viridis_c(option = "cividis", limits = c(0, 1), 
                        direction = -1, name = "Posterior\nProbability", 
                        na.value = NA) +
  geom_cladelabel(node = 28, extend = 1, 
                  label = expression(bold("Outgroups")), align = T, 
                  geom = "label", fill = "white", 
                  offset = 0.3, barsize = 2, 
                  family = "DejaVu Sans") +
  geom_cladelabel(node = 46, extend = 0.5, 
                  label = expression(bold(atop(paste(bolditalic("Lacerta viridis")), paste("(Turkey)")))), 
                  align = T, 
                  geom = "label", fill = viridis(4, end = 0.9)[4], 
                  color = c(viridis(4, end = 0.9)[4], "white"), 
                  offset = 0.3, parse = T, barsize = 2, 
                  family = "DejaVu Sans") + 
  geom_cladelabel(node = 33, extend = 0.25, 
                  label = 'bolditalic("Lacerta viridis")', align = T, 
                  geom = "label", fill = viridis(4, end = 0.9)[3], 
                  color = c(viridis(4, end = 0.9)[3], "white"), 
                  offset = 0.3, parse = T, barsize = 2, 
                  family = "DejaVu Sans") + 
  geom_cladelabel(node = 16, extend = 1.5, 
                  label = expression(bold("Adriatic clade")), align = T, 
                  geom = "label", fill = viridis(4, end = 0.9)[2], 
                  color = c(viridis(4, end = 0.9)[2], "white"), 
                  offset = 0.3, barsize = 2, 
                  family = "DejaVu Sans") + 
  geom_cladelabel(node = 42, extend = 1, 
                  label = 'bolditalic("Lacerta bilineata")', align = T, 
                  geom = "label", fill = viridis(4, end = 0.9)[1], 
                  color = c(viridis(4, end = 0.9)[1], "white"), 
                  offset = 0.3, parse = T, barsize = 2, 
                  family = "DejaVu Sans") + 
  xlim(c(-0.05, 1.25)) +
  scale_y_reverse()  + 
  theme_tree() +
  theme(legend.position = c(0.1, 0.4), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(face = "bold", size = 18, 
                  family = "DejaVu Sans")) +
  geom_treescale(x = 0, y = -7.5, width = 0.1)


ggsave(plot = g_Bayes, filename = "BayesTree.pdf", device = cairo_pdf, 
       width = 8.5, height = 8.5, dpi = 300, units = "in")
ggsave(plot = g_Bayes, filename = "BayesTree.png", device = "png", 
       width = 8.5, height = 8.5, dpi = 300)

g_Bayes_clado = ggtree(Bayes_Tree_clado, size = 1.1) + 
  geom_tiplab(align = T, hjust = 1, linetype = "dotted", 
              fontface = "bold", offset = 5, 
              geom = "label", color = colo_Bayes, fill = "white", 
              family = "DejaVu Sans Mono", size = 4.5) + 
  geom_nodepoint(aes(color = as.numeric(label)), 
                 size = 4.5, show.legend = F) + 
  scale_color_viridis_c(option = "cividis", limits = c(0, 1), 
                        direction = -1, name = "Posterior\nProbability", 
                        na.value = NA) +
  geom_cladelabel(node = 28, extend = 1, 
                  label = expression(bold("Outgroups")), align = T, 
                  geom = "label", fill = "white", 
                  offset = 5.3, barsize = 2, 
                  family = "DejaVu Sans") +
  geom_cladelabel(node = 46, extend = 0.5, 
                  label = expression(bold(atop(paste(bolditalic("Lacerta viridis")), paste("(Turkey)")))), 
                  align = T, 
                  geom = "label", fill = viridis(4, end = 0.9)[4], 
                  color = c(viridis(4, end = 0.9)[4], "white"), 
                  offset = 5.3, parse = T, barsize = 2, 
                  family = "DejaVu Sans") + 
  geom_cladelabel(node = 33, extend = 0.25, 
                  label = 'bolditalic("Lacerta viridis")', align = T, 
                  geom = "label", fill = viridis(4, end = 0.9)[3], 
                  color = c(viridis(4, end = 0.9)[3], "white"), 
                  offset = 5.3, parse = T, barsize = 2, 
                  family = "DejaVu Sans") + 
  geom_cladelabel(node = 16, extend = 1.5, 
                  label = expression(bold("Adriatic clade")), align = T, 
                  geom = "label", fill = viridis(4, end = 0.9)[2], 
                  color = c(viridis(4, end = 0.9)[2], "white"), 
                  offset = 5.3, barsize = 2, 
                  family = "DejaVu Sans") + 
  geom_cladelabel(node = 42, extend = 1, 
                  label = 'bolditalic("Lacerta bilineata")', align = T, 
                  geom = "label", fill = viridis(4, end = 0.9)[1], 
                  color = c(viridis(4, end = 0.9)[1], "white"), 
                  offset = 5.3, parse = T, barsize = 2, 
                  family = "DejaVu Sans") + 
  xlim(c(-0.05, 20.25)) +
  scale_y_reverse()  + 
  theme_tree()

ggsave(plot = g_Bayes_clado, filename = "BayesTreeClado.pdf", device = cairo_pdf, 
       width = 8.5, height = 8.5, dpi = 300, units = "in")
ggsave(plot = g_Bayes_clado, filename = "BayesTreeClado.png", device = "png", 
       width = 8.5, height = 8.5, dpi = 300)
```

Combine Tree and Phylo to Map
-----------------------------

We can plot the tree and the phylomap next to each other by assigning a function to the r base plot and then use this within the `plot_grid` function.

``` r
PhyloMap = function(){
par(family = "DejaVu Sans Mono")
plot(obj_aligned, 
     type = "phylogram", 
     direction = "downwards", 
     colors = cols_aligned, 
     ftype = "b", 
     cex.points = c(1.2, 1.3), 
     lwd = 3)
par(family = "DejaVu Sans")
legend(43.5, 63, legend = c(expression(italic("bilineata")), "Adriatic",
                          expression(italic("viridis")), 
                          expression(italic("viridis")*" (Turkey)"), 
                          "Outgroups"), 
       pt.bg = c(viridis(4, end = 0.9), "white"),
       col = "black", 
       title = expression(bold("Lacerta lineages")),
       bg = "grey85", pch = 21)
}

combi = plot_grid(g, PhyloMap, labels = "AUTO", 
                  nrow = 1, ncol = 2, 
                  label_size = 18, label_fontfamily = "DejaVu Sans", 
                  vjust = 3)


ggsave(plot = combi, filename = "Tree_PhyloToMap_Combined.pdf", 
       device = cairo_pdf, width = 16, height = 9, 
       dpi = 600)
ggsave(plot = combi, filename = "Tree_PhyloToMap_Combined.png", 
       device = "png", width = 16, height = 9, 
       dpi = 600)

combi
```

![](Readme_files/figure-markdown_github/combine%20Tree%20and%20PhyloToMap-1.png)

Pairwise Distance Heatmap
=========================

Plot the uncorrected p-distances as a Heatmap. The distance matrix was generated with the program Geneious.

``` r
distmat = as.dist(read.csv("../00_Data/Miscellaneous/pDistanceMatrix.csv", 
                           header = T, row.names = 1), 
                  diag = T, upper = T)

# Distances are here given as %identity, so we need to invert the values
# then we get the uncorrected p-distance in %
distmat = 100-distmat

# get the tip labels and sort by their y value
# this will be used to sort the distance matrix in the same way
Labels = as.data.frame(cbind(g$data[g$data$isTip == T,]$label,
                             as.numeric(g$data[g$data$isTip == T,]$y)), 
                       stringsAsFactors = F)
Labels$V2 = as.numeric(Labels$V2)
Labels = Labels[order(Labels$V2, decreasing = T),]
distmat = as.matrix(distmat)[Labels$V1, Labels$V1]

distmat_melted = reshape2::melt(as.matrix(distmat))
distmat_melted$value = ifelse(distmat_melted$Var1 == distmat_melted$Var2, 
                              NA, distmat_melted$value)
distmat_melted$HeatmapValue = ifelse(distmat_melted$value > 8.5, 
                                     8.5, distmat_melted$value)

# get rectangle positions for each lineage
rectdata = data.frame(Clade = c("Lacerta viridis", 
                                    "Lacerta viridis (Turkey)", 
                                    "Lacerta bilineata", 
                                    "Adriatic clade", 
                                    "Outgroups"), 
                      xminval = 0, xmaxval = 0, 
                      yminval = 0, ymaxval = 0)

rectdata[rectdata$Clade == "Lacerta viridis", "xminval"] = min(which(grepl("Lv_(?!.*TR)", head(distmat_melted, n = 25)$Var1, perl = T)))-0.5
rectdata[rectdata$Clade == "Lacerta viridis", "xmaxval"] = max(which(grepl("Lv_(?!.*TR)", head(distmat_melted, n = 25)$Var1, perl = T)))+0.5

rectdata[rectdata$Clade == "Lacerta viridis (Turkey)", "xminval"] = min(which(grepl("Lv_TR", head(distmat_melted, n = 25)$Var1, perl = T)))-0.5
rectdata[rectdata$Clade == "Lacerta viridis (Turkey)", "xmaxval"] = max(which(grepl("Lv_TR", head(distmat_melted, n = 25)$Var1, perl = T)))+0.5

rectdata[rectdata$Clade == "Lacerta bilineata", "xminval"] = min(which(grepl("Lb_", head(distmat_melted, n = 25)$Var1, perl = T)))-0.5
rectdata[rectdata$Clade == "Lacerta bilineata", "xmaxval"] = max(which(grepl("Lb_", head(distmat_melted, n = 25)$Var1, perl = T)))+0.5

rectdata[rectdata$Clade == "Adriatic clade", "xminval"] = min(which(grepl("AL_", head(distmat_melted, n = 25)$Var1, perl = T)))-0.5
rectdata[rectdata$Clade == "Adriatic clade", "xmaxval"] = max(which(grepl("AL_", head(distmat_melted, n = 25)$Var1, perl = T)))+0.5

rectdata[rectdata$Clade == "Outgroups", "xminval"] = min(which(grepl("[KN]C", head(distmat_melted, n = 25)$Var1, perl = T)))-0.5
rectdata[rectdata$Clade == "Outgroups", "xmaxval"] = max(which(grepl("[KN]C", head(distmat_melted, n = 25)$Var1, perl = T)))+0.5

rectdata$yminval = rectdata$xminval
rectdata$ymaxval = rectdata$xmaxval

d = ggplot(distmat_melted) + 
  geom_raster(aes(x = Var1, y = Var2, fill = HeatmapValue)) + 
  # uncomment these lines if you want to include the distance given in % per tile
  #geom_text(aes(x = Var2, y = Var1, 
  #              label = ifelse(is.na(value), "", 
  #                             ifelse(value <= 10, 
  #                                    sprintf("%0.2f", round(value, digits = 2)), 
  #                                    sprintf("%0.1f", round(value, digits = 1))))), 
  #          color = ifelse(distmat_melted$value >= 1, "white", "black"), 
  #          size = 2.5, na.rm = T) +
  geom_rect(data = rectdata, aes(xmin = xminval, xmax = xmaxval, 
                                 ymin = yminval, ymax = ymaxval, 
                                 color = Clade), 
            fill = NA, size = 2) +
  scale_fill_viridis(option = "cividis", na.value = "white", direction = -1, 
                     breaks = seq(0, 8, by = 2), 
                     labels = c(0, 2, 4, 6, ">8"), 
                     name = "Uncorrected\np distance (%)") +
  scale_color_manual(limits = c("Outgroups", "Lacerta viridis (Turkey)", 
                                "Lacerta viridis", "Adriatic clade", 
                                "Lacerta bilineata"), 
                     values = c("black", 
                                viridis(4, end = 0.9, direction = -1))) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.text = element_text(size = 12, 
                                 family = "DejaVu Sans Mono", 
                                 face = "bold"), 
        legend.text = element_text(size = 12, family = "DejaVu Sans"), 
        legend.title = element_text(size = 14, face = "bold", 
                                    family = "DejaVu Sans"), 
        panel.background = element_blank(), 
        panel.grid = element_blank()) + 
  coord_equal()


ggsave(plot = d, filename = "Uncorrected-p-Distance_Heatmap.pdf", 
       device = cairo_pdf, width = 8.5, height = 6, dpi = 300)
ggsave(plot = d, filename = "Uncorrected-p-Distance_Heatmap.png", 
       device = "png", width = 8.5, height = 6, dpi = 300)

d
```

![](Readme_files/figure-markdown_github/pDistance%20Heatmap-1.png)

Control Region Visualisation
============================

This visualisation will be a stacked bar chart to outline the structure of the mt control region adapted from [Böhme et al. 2007](https://doi.org/10.1016/j.gene.2007.02.006).

``` r
CRdata = as.data.frame(read.csv("../00_Data/Miscellaneous/DummyDataControlRegion.csv", header = T, 
                                stringsAsFactors = F))

CRdata$xmin = 0
CRdata$xmax = 1
CRdata$ymin = 0
CRdata$ymax = 0

for(i in seq(1:nrow(CRdata))){
  if(i == 1){
    CRdata[i, "ymin"] = 0
    CRdata[i, "ymax"] = CRdata[i, "Height"]}
  else {
    CRdata[i, "ymin"] = CRdata[i-1, "ymax"]
    CRdata[i, "ymax"] = CRdata[i-1, "ymax"] + CRdata[i, "Height"]
  }
}

c = ggplot(CRdata, aes(xmin = xmin, ymin = ymin, 
                   xmax = xmax, ymax = ymax, 
                   fill = Group)) + 
  geom_rect() + 
  geom_text(aes(x = (xmin+xmax)/2, y = (ymin+ymax)/2, 
                label = Label), 
            na.rm = T, angle = 90, 
            color = "white", fontface = "bold", family = "DejaVu Sans", 
            size = 6.5) +
  geom_rect(xmin = 0, xmax = 1, ymin = 0, ymax = 100, 
            color = "black", fill = NA, size = 1.5) +
  scale_fill_viridis_d(option = "cividis", na.value = "white", 
                       begin = 0.1, end = 0.8, name = "Feature", 
                       na.translate = F) +
  theme_void() + 
  coord_flip() +
  theme(legend.title = element_text(face = "bold", size = 14, 
                                    family = "DejaVu Sans"), 
        legend.text = element_text(size = 12, family = "DejaVu Sans"), 
        aspect.ratio = 0.3, 
        legend.position = "right")

ggsave(plot = c, filename = "SchematicOverwiev_Controlregion.pdf", 
       device = cairo_pdf, width = 16, height = 4, dpi = 300)
ggsave(plot = c, filename = "SchematicOverwiev_Controlregion.png", 
       device = "png", width = 16, height = 4, dpi = 300)

c
```

![](Readme_files/figure-markdown_github/Control%20Region%20Visualisation-1.png)

Motif Evolution
===============

This section covers the plotting of the sequence motiv evolution on the cladogram.

``` r
Nt_Tree = ape::rotateConstr(Nt_Tree, rev(obj$tree$tip.label))

t = ggtree(Nt_Tree, size = 1.1) + 
  geom_tiplab(align = T, hjust = 1, linetype = "dotted", 
              fontface = "bold", offset = 5, 
              geom = "label", color = colo, fill = "white", 
              family = "DejaVu Sans Mono", size = 4.5) + 
  scale_color_viridis_c(option = "cividis", limits = c(0, 100), 
                        direction = -1, name = "Bootstrap Support", 
                        na.value = NA) +
  geom_cladelabel(node = 48, extend = 1, 
                  label = expression(bold("Outgroups")), align = T, 
                  geom = "label", fill = "white", 
                  offset = 5.25, barsize = 2, 
                  family = "DejaVu Sans") +
  geom_cladelabel(node = 30, extend = 0.25, 
                  label = expression(bold(atop(paste(bolditalic("Lacerta viridis")), paste("(Turkey)")))), 
                  align = T, 
                  geom = "label", fill = viridis(4, end = 0.9)[4], 
                  color = c(viridis(4, end = 0.9)[4], "white"), 
                  offset = 5.25, parse = T, barsize = 2, 
                  family = "DejaVu Sans") + 
  geom_cladelabel(node = 32, extend = 0.25, 
                  label = 'bolditalic("Lacerta viridis")', align = T, 
                  geom = "label", fill = viridis(4, end = 0.9)[3], 
                  color = c(viridis(4, end = 0.9)[3], "white"), 
                  offset = 5.25, parse = T, barsize = 2, 
                  family = "DejaVu Sans") + 
  geom_cladelabel(node = 12, extend = 1.5, 
                  label = expression(bold("Adriatic clade")), align = T, 
                  geom = "label", fill = viridis(4, end = 0.9)[2], 
                  color = c(viridis(4, end = 0.9)[2], "white"), 
                  offset = 5.25, barsize = 2, 
                  family = "DejaVu Sans") + 
  geom_cladelabel(node = 44, extend = 1, 
                  label = 'bolditalic("Lacerta bilineata")', align = T, 
                  geom = "label", fill = viridis(4, end = 0.9)[1], 
                  color = c(viridis(4, end = 0.9)[1], "white"), 
                  offset = 5.25, parse = T, barsize = 2, 
                  family = "DejaVu Sans") + 
  scale_y_reverse()  + 
  theme_tree() 
  

# look at the nodepoints and add the annotation manually
# with the vectors we build now

TAS2a_tip = c(1:2, 12:14, 16:20)
TAS2a_node = 29
TAS2b_tip = 4:11
TAS2b_node = 33
TAS2c_tip = c(3, 15)
Insert_tip = 3:11
Insert_node = 32

tl = t +
  geom_tippoint(aes(subset = (node %in% TAS2a_tip), 
                    color = "TAS2a", 
                    shape = "TAS2a"), 
                size = 4.5)+
  geom_nodepoint(aes(subset = (node %in% TAS2a_node), 
                     color = "TAS2a", 
                     shape = "TAS2a"), 
                 size = 4.5)+
  geom_tippoint(aes(subset = (node %in% TAS2b_tip), 
                    color = "TAS2b", 
                    shape = "TAS2b"), 
                size = 4.5)+
  geom_nodepoint(aes(subset = (node %in% TAS2b_node), 
                     color = "TAS2b", 
                     shape = "TAS2b"), 
                 size = 4.5)+
  geom_tippoint(aes(subset = (node %in% TAS2c_tip), 
                    color = "TAS2c", 
                    shape = "TAS2c"), 
                size = 4.5)+ 
  geom_tippoint(aes(subset = (node %in% Insert_tip), 
                    color = "11bp Insertion", 
                    shape = "11bp Insertion"), 
                size = 4.5, 
                position = position_nudge(x = -0.4)) + 
  geom_nodepoint(aes(subset = (node %in% Insert_node), 
                     color = "11bp Insertion", 
                     shape = "11bp Insertion"), 
                 size = 4.5)+
  scale_color_manual(values = cividis(4, begin = 0.1, 
                                      end = 0.95, direction = 1), 
                     labels = c("TAS2a", "TAS2b", "TAS2c", "11bp Insertion"), 
                     breaks = c("TAS2a", "TAS2b", "TAS2c", "11bp Insertion"),
                     name = "Sequence Motif") +
  scale_shape_manual(values = c(16, 16, 16, 15), 
                     labels = c("TAS2a", "TAS2b", "TAS2c", "11bp Insertion"), 
                     breaks = c("TAS2a", "TAS2b", "TAS2c", "11bp Insertion"),
                     name = "Sequence Motif") +
  theme(legend.position = c(0.15, 0.15), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(face = "bold", size = 18, 
                  family = "DejaVu Sans")) +
  xlim(c(0, 20.5))

ggsave(plot = tl, filename = "SequenceMotifEvolution.pdf", 
       device = cairo_pdf, width = 9, height = 8.5, dpi = 300)
ggsave(plot = tl, filename = "SequenceMotifEvolution.png", 
       device = "png", width = 9, height = 8.5, dpi = 300)

tl  
```

![](Readme_files/figure-markdown_github/MotifEvolution-1.png)
