#' ComponentPlot
#'
#' Plots two columns from the meta.data within a \code{Seurat} object to visualise how cells compare between two identities
#'
#' @param seurat \code{Seurat} object for information to be drawn from
#' @param identity The identity to be plotted on the x axis
#' @param component The second identity, plotted on the y axis and will show the breakdown of this identity within the first identity
#' @param show.pct Whether to show the total cells (FALSE, default) or the percentage of cells (TRUE)
#'
#' @return The graph
#'
#'#' @importFrom ggplot2 ggplot aes geom_bar coord_flip geom_text position_stack element_text ggtitle rel
#'
#'
#'@export
#'
ComponentBar <- function(seurat.md, identity, component, show.pct = FALSE, show.text = 25) {
  Pct = NULL
  Cells = NULL
  Component = NULL
  Identity = NULL

  df <- as.data.frame(table(seurat.md[[paste0(identity)]],
                            seurat.md[[paste0(component)]]))
  id <- as.data.frame(table(seurat.md[[paste0(identity)]]))
  colnames(id) <- c("Identity", "CellTotal")
  colnames(df) <- c("Identity", "Component", "Cells")
  df <- dplyr::left_join(df, id, "Identity")
  df$Pct <- round(df$Cells / df$CellTotal * 100, digits = 2)
  
  if (show.pct == FALSE) {
    p <- ggplot2::ggplot(df, ggplot2::aes(Identity, Cells, fill = Component))
    id <- "Cells and Percentage"
    
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(Identity, Pct, fill = Component))
    id <- "Percentage"
    
  }
  p +
    ggplot2::geom_bar(stat = "Identity", colour = "black", width = 1, position = position_stack(reverse = TRUE)) +
    ggplot2::geom_text(ggplot2::aes(label=ifelse(Pct >= show.text, levels(Component)[Component], "")), size = 4,
                       position=ggplot2::position_stack(vjust=0.5, reverse = TRUE), colour="black") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = -45, hjust = 0, vjust = 0.5)) +
    ggplot2::theme(legend.title=ggplot2::element_text(size=ggplot2::rel(1.1))) +
    ggplot2::ggtitle(paste0(id, " of ", component, " in ", identity))
}


#' StrDotPlot
#'
#' @param x
#' @param features
#' @param group.by
#' @param assay
#' @param col.min
#' @param col.max
#'
#' @return
#' @export
#'
#' @examples
StrDotPlot <- function(x, features, group.by = "Identity", assay = "RNA", col.min = 0,
                       col.max = 5, dot.min = 0, dot.scale = 6, scale = F, rainbow = F) {
  if (rainbow){
    c <- scale_color_gradientn(colours = rainbow(5))
  } else {
    c <- scale_colour_gradient(low = "lightgrey", high = "red")
  }
  DotPlot(object = x, features = features, col.min = col.min, col.max = col.max,
          cols = c("lightgrey", "red"), dot.scale = dot.scale,
          group.by = group.by, dot.min = dot.min, scale = scale,
          assay = assay) +
    c +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 1),
          panel.grid.major = element_line(colour = "lightgray")) +
    theme(legend.title=element_text(size=rel(0.5)))
}



'%!in%' <- function(x,y)!('%in%'(x,y))



vcol <- viridis::viridis_pal()(10)



#' DimDot_Comparison
#'
#' @param seurat
#' @param identity
#' @param features
#' @param label
#'
#' @return
#' @export
#'
#' @examples
DimDot_Comparison <- function(seurat,
                              identity,
                              features,
                              label=T){
  DimPlot(seurat, group.by = identity, label = label) + (DotPlot(seurat, group.by = identity, features = features) +
                                                           coord_flip() +
                                                           theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5),
                                                                 panel.grid.major = element_line(colour = "lightgray")) +
                                                           theme(legend.title=element_text(size=rel(0.5))))
}


#' MultiDotPlot
#'
#' @param seurat
#' @param cluster.sequence
#' @param features
#'
#' @return
#' @export
#'
#' @examples
MultiDotPlot <- function(seurat,
                         cluster.sequence,
                         features) {
  seq <- cluster.sequence
  clusters <- c(paste0("SCT_snn_res.", seq))
  
  #plot <- map(clusters, function(x) DimPlot(seurat, group.by = x, label = T))
  plot <- map(clusters, function(x) (DotPlot(seurat, group.by = x, features = features) +
                                       coord_flip() +
                                       ggtitle(label = paste0(x)) +
                                       theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5),
                                             panel.grid.major = element_line(colour = "lightgray")) +
                                       theme(legend.title=element_text(size=rel(0.5)))))
  wrap_plots(plot, ncol = round(sqrt(length(cluster.sequence))))
}




#' MultiDimPlot
#'
#' @param seurat
#' @param cluster.sequence
#'
#' @return
#' @export
#'
#' @examples
MultiDimPlot <- function(seurat,
                         cluster.sequence,
                         reduction = "umap") {
  seq <- cluster.sequence
  clusters <- c(paste0("SCT_snn_res.", seq))
  
  plot <- map(clusters, function(x) DimPlot(seurat, group.by = x, label = T, reduction = reduction) +
                ggtitle(label = paste0(x)))
  
  wrap_plots(plot, ncol = round(sqrt(length(cluster.sequence))))
}

##' SankeyPlot
##' @description Plots a Sankey diagram shown how the cells are from a category into the final classification
##'
##' @param data seurat object fot the data to be pulled from
##' @param origin the meta.data column describing the origin identity
##' @param height the height of the plot
##' @param width the width of the plot
##' @param fontsize the fontsize of the plot
##' @param simple defualt is FALSE; TRUE changes to a direct plot going from the origin to final classification without intermediate steps shown
##' @param forced default FALSE IF TRUE changes the outcome to show all cells force called into a cell type instead of being unassigned
##' @param sinksRight default FALSE; graphing parameter from base plotting function that forces final identity to right side of graph
##' @param LinkGroup Changes the way in which linking colours are plotted. Default is "Var3" representing the origin population
##' @param remove.unassigned default FALSE; determines if unassigned values are removed from the plot
##' @param nephron.only default FALSE; TRUE will plot only cell types that are classified into a nephron lineage
##'
##' @importFrom networkD3 sankeyNetwork
##'
##' @return the plot
##' @export
##'
##'
#SankeyPlot <- function(data,
#                       simple = T,
#                       origin = "sample.id",
#                       end = "DKCC",
#                       height = 900,
#                       width = 1200,
#                       fontsize = 12,
#                       sinksRight = TRUE,
#                       LinkGroup = "Var1",
#                       remove.unassigned = F,
#                       nephron.only = F){
#
#  if (remove.unassigned==T){
#    .data$data <- data[, data$DKCC != "unassigned"]
#
#  }
#
#  if (nephron.only == T){
#    .data$data[, data$DKCC %!in% c("Stroma", "unassigned", "Non_Renal")]
#  }
#
#
#
#  if (simple == T){
#    links <- as.data.frame(table(data@meta.data[[origin]], data@meta.data[[end]]))
#  } else {
#    links <- as.data.frame(table(data@meta.data[[origin]], data@meta.data$LineageID)) %>% mutate(Var3 = .data$Var1)
#    links <- bind_rows(links, (as.data.frame(table(data$LineageID, data$DKCC, data@meta.data[[origin]]),
#                                             stringsAsFactors = F) %>% filter(Var1 != Var2))) %>% filter(.data$Freq!=0)
#
#  }
#  links <- filter(links, .data$Freq != 0)
#
#  nodes <- data.frame(name = c(as.character(links$Var1), as.character(links$Var2)),
#                      stringsAsFactors = F) %>% unique()
#
#  links$IDsource <- match(links$Var1, nodes$name)-1
#  links$IDtarget <- match(links$Var2, nodes$name)-1
#
#  p <- sankeyNetwork(Links = links, Nodes = nodes,
#                     Source = "IDsource", Target = "IDtarget",
#                     Value = "Freq", NodeID = "name",
#                     sinksRight=sinksRight, fontSize = fontsize,
#                     height = height, width = width, LinkGroup = LinkGroup)
#  return(p)
#
#}





ComponentHeatMapDF <- function(seurat.md, identity, component, scale = "none", show.pct = FALSE, tidy = T) {
  Pct = NULL
  Cells = NULL
  Component = NULL
  Identity = NULL
  
  df <- as.data.frame(table(as.character(seurat.md[[paste0(identity)]]),
                            as.character(seurat.md[[paste0(component)]])))
  id <- as.data.frame(table(seurat.md[[paste0(identity)]]))
  colnames(id) <- c("Identity", "CellTotal")
  colnames(df) <- c("Identity", "Component", "Cells")
  df <- dplyr::left_join(df, id, "Identity")
  df$Pct <- round(df$Cells / df$CellTotal * 100, digits = 2)
  if (!tidy){
    df <- pivot_wider(df[c("Identity", "Component", "Pct")], names_from = "Component", values_from = "Pct") %>% column_to_rownames("Identity")
  }
  return(df)
}
