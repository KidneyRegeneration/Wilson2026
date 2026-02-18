library(viridis)

subset_fs19 <- function(x = 0.1, by = "capture", times = 1) {
  fs19 <- readr::read_rds(here::here("output/rds/FS19.rds"))
  ds <- caret::createDataPartition(fs19@meta.data[, by], times = times, p = x)
  if (times==1){
    fs19 <- fs19[, ds$Resample1]
  } else {
    fs19 <- map(ds, ~fs19[, .x])
  }
  return(fs19)
}

subset_seu <- function(seu, x = 0.1, by = "capture", times = 1) {
  
  ds <- caret::createDataPartition(seu@meta.data[, by], times = times, p = x)
  if (times==1){
    seu <- seu[, ds$Resample1]
  } else {
    seu <- map(ds, ~seu[, .x])
  }
  return(seu)
}

plotly3d <- function(seu, col) {
  plotly::plot_ly(data.frame(cell = colnames(seu),
                     dim1 = seu@reductions$umap3d@cell.embeddings[,1],
                     dim2 = seu@reductions$umap3d@cell.embeddings[,2],
                     dim3 = seu@reductions$umap3d@cell.embeddings[,3]),
          x = ~dim1, 
          y = ~dim2,
          z = ~dim3, type="scatter3d", mode = 'markers',
          marker = list(opacity = 0.7, size=2),
          
          color = seu@meta.data[,col])
}

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

get.psudobulk <- function(seu, column){
  pb <- data.frame(row.names = rownames(seu@assays$RNA@counts))
  for (i in 1:length(unique(seu@meta.data[,column]))){
    temp.seurat <- seu[, seu@meta.data[,column] == unique(seu@meta.data[,column])[i]]
    #temp.seurat <- subset(matorg, ident = unique(matorg$DKCC)[i])
    temp.counts <- as.data.frame(temp.seurat@assays$RNA@counts)
    temp.bulk <- data.frame(rowSums(temp.counts))
    colnames(temp.bulk) <- c(unique(as.character(seu@meta.data[, column]))[i])
    pb <- cbind(pb, temp.bulk)
  }
  pb
}

cc <- viridis::viridis(3)

pseudobulk <- function(seu, ident) {
  counts <- data.frame(row.names = rownames(seu@assays$RNA@counts))
  for (i in 1:length(unique(seu@meta.data[,ident]))){
    temp.seurat <- seu[, seu@meta.data[,ident] == unique(seu@meta.data[,ident])[i]]
    #temp.seurat <- subset(matorg, ident = unique(matorg$DKCC)[i])
    temp.counts <- as.data.frame(temp.seurat@assays$RNA@counts)
    temp.bulk <- data.frame(rowSums(temp.counts) %>% edgeR::cpm())
    colnames(temp.bulk) <- c(unique(as.character(seu@meta.data[,ident]))[i])
    counts <- cbind(counts, temp.bulk)
  }
  return(counts)
}



allcols <- read_rds(here::here("data/rds/allcols.rds"))


#gcols <- ggplotColors(5)

DFDotPlot <- function (object, assay = NULL, features, cols = c("lightgrey", 
                                                                "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                       idents = NULL, group.by = NULL, split.by = NULL, cluster.idents = FALSE, 
                       scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) 
{
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% 
                                                   rownames(x = brewer.pal.info))
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(X = 1:length(features), 
                                        FUN = function(x) {
                                          return(rep(x = names(x = features)[x], each = length(features[[x]])))
                                        }))
    if (any(is.na(x = feature.groups))) {
      warning("Some feature groups are unnamed.", call. = FALSE, 
              immediate. = TRUE)
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
  data.features <- FetchData(object = object, vars = features, 
                             cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  }
  else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, 
                              sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(what = rbind, args = lapply(X = data.plot, 
                                               FUN = unlist))
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  ngroup <- length(x = levels(x = data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled", 
            call. = FALSE, immediate. = TRUE)
  }
  else if (ngroup < 5 & scale) {
    warning("Scaling data with a low number of groups may produce misleading results", 
            call. = FALSE, immediate. = TRUE)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             if (scale) {
                               data.use <- scale(x = data.use)
                               data.use <- MinMax(data = data.use, min = col.min, 
                                                  max = col.max)
                             }
                             else {
                               data.use <- log1p(x = data.use)
                             }
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(X = as.character(x = data.plot$id), 
                         FUN = gsub, FUN.VALUE = character(length = 1L), 
                         pattern = paste0("^((", paste(sort(x = levels(x = object), 
                                                            decreasing = TRUE), collapse = "|"), ")_)"), 
                         replacement = "", USE.NAMES = FALSE)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = split.colors, yes = "colors", 
                     no = "avg.exp.scaled")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot], 
                                       levels = unique(x = feature.groups))
  }
  
  return(data.plot)
}
