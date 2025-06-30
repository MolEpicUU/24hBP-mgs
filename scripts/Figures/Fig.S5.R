
rm(list=ls())
#set the seed to make it reproducible
set.seed(123)

##load libraries
library(rio)
library(BiocParallel)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(sandwich)
library(lmtest)

met=res1$Metabolite
res <- res2[res2$Metabolite %in% met, ]
#res=res3[,c("MGS","Metabolite","Estimate","q.value")]


#complete with all the models
res.plot.model=dcast(res,MGS~Metabolite,value.var = "Estimate")
res.plot.model.q=dcast(res, MGS ~ Metabolite,value.var="q.value")
rownames(res.plot.model)=rownames(res.plot.model.q)=res.plot.model[,1]
res.plot.model=res.plot.model[,2:ncol(res.plot.model)]
res.plot.model.q=res.plot.model.q[,2:ncol(res.plot.model.q)]


# to list the x-axis
#res.plot.model <- res.plot.model[,c(7,8,9,1,2,3,10,11,12,4,5,6)]
#res.plot.model.q <- res.plot.model.q[,c(7,8,9,1,2,3,10,11,12,4,5,6)]
#res.plot.model <- as.matrix(res.plot.model)

#res.plot.model <-res.plot.model[,c(7,8,9,1,2,3,10,11,12,4,5,6)]
#res.plot.model.q <- res.plot.model.q[,c(7,8,9,1,2,3,10,11,12,4,5,6)]
res.plot.model <- as.matrix(res.plot.model)
res.plot.model.q <- as.matrix(res.plot.model.q)
res.plot.model.q <- ifelse(res.plot.model.q<0.05,"**","")



#############################
# Match colnames of `a` with `b$metabolite` and get the corresponding `subclass`
matched_subclass <- res1[res1$Metabolite %in% colnames(res.plot.model), c("Metabolite", "Subclass")]

# Preserve the original order of colnames in `a`
#matched_subclass <- matched_subclass[match(colnames(res1), matched_subclass$Metabolite, noMatch = 0), ]
##############################


# Generate a softer, more natural color palette
nature_colors <- c(
  "#8ECAE6",  # Soft blue
  "#95B8D1",  # Muted blue-gray
  "#B8D4E3",  # Light steel blue
  "#219EBC",  # Teal
  "#A8DADC",  # Pale cyan
  "#E9C46A",  # Muted yellow
  "#F4A261",  # Soft orange
  "#CCD5AE",  # Sage green
  "#E9EDC9",  # Pale olive
  "#FEFAE0"   # Cream
)

# Generate the annotation with the nature-inspired colors
num_subclasses <- length(unique(res$Metabolite_subclass))

# If we need more colors than provided, interpolate between existing ones
if (num_subclasses > length(nature_colors)) {
  nature_colors <- colorRampPalette(nature_colors)(num_subclasses)
}

bottom_annotation <- HeatmapAnnotation(
  subclass = matched_subclass$Subclass[c(5,4,6,11,7,12,2,1,9,3)],
  col = list(
    subclass = structure(
      nature_colors[1:num_subclasses],
      names = unique(matched_subclass$Subclass)
    )
  ),
  annotation_name_side = "left",
  annotation_legend_param = list(
    title = "Subclass",
    labels_gp = gpar(fontsize = 8),
    title_gp = gpar(fontsize = 8)
  )
)

h1 = Heatmap(res.plot.model, 
             #top_annotation = ha1,
             column_names_rot =  30,
             column_names_side = "top",
             cluster_columns = FALSE,
             cluster_rows = TRUE,
             #column_title = "Univariate",
             #name="BPs  (mmHg)",
             #column_labels = c("SBP", "DBP", "SBP-SD", "DBP-SD"),
             column_names_gp = gpar(fontsize = 11),
             #column_title_gp = gpar(fontsize = 7,fontface='bold'),
             heatmap_legend_param = list(labels_gp = gpar(fontsize=9),
                                         title_gp = gpar(fontsize=8),
                                         grid_width = unit(5, "mm")),
             col = colorRamp2(
               c(-1, -0.5, -0.2, 0, 0.2, 0.5, 1),
               c("#8BA5D9", "#A8C0E0", "#CAD5E0", "white", "#F7C4BB", "#F4A69B", "#EC887B")
             ),
             #cell_fun = function(j, i, x, y, width, height, fill) 
             #{
             #   grid.text(sprintf("%s", res.plot.model.q[i, j]), x, y, gp = gpar(fontsize = 8))
             #},
             #left_annotation = ha, 
             show_row_dend = FALSE,
             row_names_side = "left",
             row_names_gp = gpar(fontsize = 11,fontface = "italic"), 
             #row_names_gp = gpar(fontsize = 8, fontface = "italic"), 
             row_dend_width = unit(1.3, "cm"), 
             clustering_distance_rows = "euclidean",
             bottom_annotation = bottom_annotation
             #layer_fun = function(j, i, x, y, w, h, fill){
             #   pindex(hm.matrix1,j, i)
             #   grid.rect(gp = gpar(lwd = 2, fill = "transparent"))
             #  }
)   

h1

dev.off()
