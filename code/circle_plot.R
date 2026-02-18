CoruscantPlot <- function(seurat, column, group, segment, show.pct = T, plot.out = "globe", get.df = F, 
                          ylim=20, lab.dist=6, text.size=3, position = "stack") {
  Pct = NULL
  Cells = NULL
  Component = NULL
  Identity = NULL
  df <- as.data.frame(table(seurat@meta.data[[paste0(column)]],
                            seurat@meta.data[[paste0(group)]],
                            seurat@meta.data[[paste0(segment)]]))
  
  id <- as.data.frame(table(seurat[[paste0(group)]]))
  
  colnames(id) <- c("Group", "CellTotal")
  colnames(df) <- c("Identity", "Group", "Segment", "Cells")
  df <- dplyr::left_join(df, id, "Group")
  if (show.pct){
    df$Pct <- round(df$Cells / df$CellTotal * 100, digits = 2)
  } else {
    df$Pct <- Cells
  }
  
  #df$Identity <- df$Identity %>% as.character()
  #df$Group <- df$Group %>% as.character()
  #df$Segment <- df$Segment %>% as.character()
  df$ig <- paste0(df$Identity, df$Group)
  a <- paste0(seurat@meta.data[[paste0(column)]], seurat@meta.data[[paste0(group)]])
  df <- df %>% filter(ig %in% unique(a))
  data <- df
  
  if (get.df==T){
    data
  } else {
    #data <- data %>%  gather(key = "observation", value="value", -c(1,2)) 
    #data <- data %>% filter(CellTotal>0, Group != "iPSC") %>% dplyr::select(-c("ig"))
    
    data <- data %>% filter(CellTotal>0) %>% dplyr::select(-c("ig"))
    
    target <- data.frame(identity = c(unique(seurat[[column]]))[[column]], stringsAsFactors = F)
    
    idents <- factor(target$identity, levels = target$identity)
    
    #data$Identity <- factor(data$Identity, levels = target$identity)
    #data$Group <- factor(data$Group, levels = unique(data$Group))
    
    # Set a number of 'empty bar' to add at the end of each group
    empty_bar <- 2
    nObsType <- nlevels(as.factor(data$Segment))
    to_add <- data.frame( matrix(NA, empty_bar*length(unique(data$Group))*nObsType, ncol(data)) )
    colnames(to_add) <- colnames(data)
    to_add$Group <- rep(unique(data$Group), each=empty_bar*nObsType )
    data <- rbind(data, to_add)
    data <- data %>% arrange(Group, Identity)
    
    data$id <- rep( seq(1, nrow(data)/nObsType) , each=nObsType)
    
    
    # Get the name and the y position of each label
    label_data <- data %>% group_by(id, Identity) %>% summarize(tot=sum(Pct))
    number_of_bar <- nrow(label_data)
    angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
    label_data$hjust <- ifelse( angle < -90, 1, 0)
    label_data$angle <- ifelse(angle < -90, angle+180, angle)
    
    # prepare a data frame for base lines
    base_data <- data %>% 
      group_by(Group) %>% 
      summarize(start=min(id), end=max(id) - empty_bar) %>% 
      rowwise() %>% 
      mutate(title=mean(c(start, end)))
    base_data$CellTotal <- (id %>% filter(CellTotal>0))$CellTotal
    
    # prepare a data frame for grid (scales)
    grid_data <- base_data
    grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
    grid_data$start <- grid_data$start - 1
    grid_data <- grid_data[-1,]
    
    
    # Make the plot
    flat <- ggplot(data) +      
      
      # Add the stacked bar
      geom_bar(aes(x=as.integer(id), y=Pct, fill=Segment), colour = "black", stat="identity", alpha=0.75, position = position) +
      viridis::scale_fill_viridis(discrete=TRUE) +
      
      # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
      geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
      geom_segment(data=grid_data, aes(x = end, y = 25, xend = start, yend = 25), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
      geom_segment(data=grid_data, aes(x = end, y = 50, xend = start, yend = 50), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
      geom_segment(data=grid_data, aes(x = end, y = 75, xend = start, yend = 75), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
      geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
      
      # Add text showing the value of each 100/75/50/25 lines
      ggplot2::annotate("text", x = rep(-1,5), y = c(0, 25, 50, 75, 100), label = c("0", "25", "50", "75", "100") , color="grey", size=6 , angle=0, fontface="bold", hjust=0.5) +
      
      ylim(-5,max(label_data$tot+25, na.rm=T)) +
      xlim(-4, max(data$id)) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(2,4), "cm") 
      ) +
      # coord_polar() +
      
      # Add labels on top of each bar
      #geom_text(data=label_data, aes(x=id, y=tot+10, label=Identity, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=1, angle= label_data$angle, inherit.aes = FALSE ) +
      ###geom_text(data=label_data, aes(x=id, y=label_data$tot+10, label=Identity), color="black", fontface="bold",alpha=0.6, size=3, angle= 90, inherit.aes = FALSE) +
      # Add base line information
      geom_segment(data=base_data, aes(x = start, y = -2, xend = end, yend = -2), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
      #geom_text(data=base_data, aes(x = title, y = -5, label=Group), hjust=c(1,1,1,1,0,0,0), vjust=c(1,1,0,-1,-1,1,1), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
      geom_text(data=base_data, aes(x = title, y = -5, label=Group), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
    
    
    # Save at png
    #ggsave(p, file="output.png", width=10, height=10)
    globe <- ggplot(data) +      
      
      # Add the stacked bar
      #geom_hline(yintercept = 50, alpha=0.3, colour="red", size=0.3) +
      geom_bar(aes(x=as.integer(id), y=Pct, fill=Segment), colour = "black", stat="identity", alpha=0.75, position = position) +
      #geom_bar(data = base_data, aes(x=as.integer(Group), y=-CellTotal), colour = "black", stat="identity", alpha=0.75) +
      
      viridis::scale_fill_viridis(discrete=TRUE) +
      
      # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
      geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
      geom_segment(data=grid_data, aes(x = end, y = 25, xend = start, yend = 25), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
      geom_segment(data=grid_data, aes(x = end, y = 50, xend = start, yend = 50), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
      geom_segment(data=grid_data, aes(x = end, y = 75, xend = start, yend = 75), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
      geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 101), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
      
      # Add text showing the value of each 100/75/50/25 lines
      ggplot2::annotate("text", x = rep(max(data$id),5), y = c(0, 25, 50, 75, 100), label = c("0", "25", "50", "75", "100") , color="grey", size=4 , angle=0, fontface="bold", hjust=1) +
      
      ylim(-ylim,max(label_data$tot+10, na.rm=T)) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        #plot.margin = unit(rep(2,4), "cm") 
      ) +
      coord_polar() +
      
      # Add labels on top of each bar
      #geom_text(data=label_data, aes(x=id, y=tot+10, label=Identity, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=1, angle= label_data$angle, inherit.aes = FALSE ) +
      geom_text(data=label_data, aes(x=id, y=label_data$tot+4, label=Identity), color="black", fontface="bold",alpha=0.6, size=text.size, angle= label_data$angle, hjust = label_data$hjust, inherit.aes = FALSE) +
      # Add base line information
      geom_segment(data=base_data, aes(x = start, y = -2, xend = end, yend = -2), colour = "black", alpha=0.8, size=1 , inherit.aes = FALSE )  +
      #geom_text(data=base_data, aes(x = title, y = -6, label=Group), hjust=c(0.75,1.25,0,0), vjust=c(1.2,2,-1,0), colour = "black", alpha=0.8, size=3, fontface="bold", inherit.aes = FALSE)
      geom_text(data=base_data, aes(x = title, y = -lab.dist, label=Group), colour = "black", alpha=0.8, size=text.size, fontface="bold", inherit.aes = FALSE) +
      #geom_text(data=base_data, aes(x = title, y = -4, label=Group), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
      
      geom_text(data=base_data, aes(x = title+1, y = -((lab.dist)*2), label=CellTotal), colour = "black", alpha=0.8, size=text.size*0.75, fontface="bold", inherit.aes = FALSE)
    
    
    if (plot.out == "both"){
      flat | globe
    } else if (plot.out == "flat"){
      flat
    } else {
      globe
    }
    
    
  }
}

#CoruscantPlot(full.timecourse,
#               column = "global_high_res_ann", 
#               group = "global_low_res_ann",
#               segment = "chir", show.pct =  T, plot.out = "globe") 



EmpirePlot <- function(seurat, column, group, segment, show.pct = T, get.df = F, lab.dist=50, text.size=4) {
  Pct = NULL
  Cells = NULL
  Component = NULL
  Identity = NULL
  df <- as.data.frame(table(seurat@meta.data[[paste0(column)]],
                            seurat@meta.data[[paste0(group)]],
                            seurat@meta.data[[paste0(segment)]]))
  
  id <- as.data.frame(table(seurat[[paste0(group)]]))
  
  colnames(id) <- c("Group", "CellTotal")
  colnames(df) <- c("Identity", "Group", "Segment", "Cells")
  df <- dplyr::left_join(df, id, "Group")
  if (show.pct){
    df$Pct <- round(df$Cells / df$CellTotal * 100, digits = 2)
  } else {
    df$Pct <- Cells
  }
  
  #df$Identity <- df$Identity %>% as.character()
  #df$Group <- df$Group %>% as.character()
  #df$Segment <- df$Segment %>% as.character()
  df$ig <- paste0(df$Identity, df$Group)
  a <- paste0(seurat@meta.data[[paste0(column)]], seurat@meta.data[[paste0(group)]])
  df <- df %>% filter(ig %in% unique(a))
  data <- df
  
  if (get.df==T){
    data
  } else {
    #data <- data %>%  gather(key = "observation", value="value", -c(1,2)) 
    #data <- data %>% filter(CellTotal>0, Group != "iPSC") %>% dplyr::select(-c("ig"))
    
    data <- data %>% filter(CellTotal>0) %>% dplyr::select(-c("ig"))
    
    target <- data.frame(identity = c(unique(seurat[[column]]))[[column]], stringsAsFactors = F)
    
    idents <- factor(target$identity, levels = target$identity)
    
    #data$Identity <- factor(data$Identity, levels = target$identity)
    #data$Group <- factor(data$Group, levels = unique(data$Group))
    
    # Set a number of 'empty bar' to add at the end of each group
    empty_bar <- 2
    nObsType <- nlevels(as.factor(data$Segment))
    to_add <- data.frame( matrix(NA, empty_bar*length(unique(data$Group))*nObsType, ncol(data)) )
    colnames(to_add) <- colnames(data)
    to_add$Group <- rep(unique(data$Group), each=empty_bar*nObsType )
    data <- rbind(data, to_add)
    data <- data %>% arrange(Group, Identity)
    
    data$id <- rep( seq(1, nrow(data)/nObsType) , each=nObsType)
    
    
    
    # Get the name and the y position of each label
    label_data <- data %>% group_by(id, Identity) %>% summarize(tot=sum(Pct), totcells = sum(Cells))
    number_of_bar <- nrow(label_data)
    angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
    label_data$hjust <- ifelse( angle < -90, 1, 0)
    label_data$angle <- ifelse(angle < -90, angle+180, angle)
    label_data$Label <- ifelse(angle < -90, paste0(label_data$totcells, "  ", label_data$Identity), paste0(label_data$Identity, "  ", label_data$totcells))
    label_data$Label <- gsub("NA  NA", NA, label_data$Label)
    
    # prepare a data frame for base lines
    base_data <- data %>% 
      group_by(Group) %>% 
      summarize(start=min(id), end=max(id) - empty_bar) %>% 
      rowwise() %>% 
      mutate(title=mean(c(start, end)))
    base_data$CellTotal <- (id %>% filter(CellTotal>0))$CellTotal
    
    # prepare a data frame for grid (scales)
    grid_data <- base_data
    grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
    grid_data$start <- grid_data$start - 1
    grid_data <- grid_data[-1,]
    
    # prep the labels
    
    #if (max(data$Pct[!is.na(data$Pct)]) < 25){
    #  lab <- ggplot2::annotate("text", x = rep(max(data$id),2), y = c(0, -25), label = c("0", "25") , color="grey", size=4 , angle=0, fontface="bold", hjust=1)
    #  ylim <- 30
    #  
    #} else if (max(data$Pct[!is.na(data$Pct)]) < 50){
    #  lab <- ggplot2::annotate("text", x = rep(max(data$id),3), y = c(0, -25, -50), label = c("0", "25", "50") , color="grey", size=4 , angle=0, fontface="bold", hjust=1)
    #  ylim <- 60
    #  
    #} else if (max(data$Pct[!is.na(data$Pct)]) < 75){
    #  lab <- ggplot2::annotate("text", x = rep(max(data$id),4), y = c(0, -25, -50, -75), label = c("0", "25", "50", "75") , color="grey", size=4 , angle=0, fontface="bold", hjust=1)
    #  ylim <- 85
    #  
    #} else {
    #  lab <- ggplot2::annotate("text", x = rep(max(data$id),5), y = c(0, -25, -50, -75, -100), label = c("0", "25", "50", "75", "100") , color="grey", size=4 , angle=0, fontface="bold", hjust=1)
    #  ylim <- 115
    #  
    #}
    
    
    max.pct <- max(data$Pct[!is.na(data$Pct)])
    lab <- ggplot2::annotate("text", x = rep(max(data$id),3), y = c(0, -(round((max.pct+5)/2)), -(max.pct+5)), label = c("0", paste0(round((max.pct+5)/2)), paste0(round(max.pct+5))),
                             color="grey", size=4 , angle=0, fontface="bold", hjust=1)
    ylim <- max.pct+10
    
    
    
    globe <- ggplot(data) +      
      
      # Add the stacked bar
      #geom_hline(yintercept = 50, alpha=0.3, colour="red", size=0.3) +
      geom_bar(aes(x=as.integer(id), y=-Pct, fill=Segment), colour = "black", stat="identity", alpha=0.75) +
      #geom_bar(data = base_data, aes(x=as.integer(Group), y=-CellTotal), colour = "black", stat="identity", alpha=0.75) +
      
      viridis::scale_fill_viridis(discrete=TRUE) +
      
      # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
      geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
      geom_segment(data=grid_data, aes(x = end, y = -25, xend = start, yend = -25), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
      geom_segment(data=grid_data, aes(x = end, y = -50, xend = start, yend = -50), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
      geom_segment(data=grid_data, aes(x = end, y = -75, xend = start, yend = -75), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
      geom_segment(data=grid_data, aes(x = end, y = -100, xend = start, yend = -101), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
      
      # Add text showing the value of each 100/75/50/25 lines
      lab + 
      
      
      
      ylim(-ylim,max(label_data$tot+10, na.rm=T)) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        #plot.margin = unit(rep(2,4), "cm") 
      ) +
      coord_polar() +
      
      # Add labels on top of each bar
      geom_text(data=label_data, aes(x=id, y=5, label=Label, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=text.size, angle= label_data$angle, inherit.aes = FALSE ) +
      #geom_text(data=label_data, aes(x=id, y=label_data$tot+4, label=Identity), color="black", fontface="bold",alpha=0.6, size=text.size, angle= label_data$angle, hjust = label_data$hjust, inherit.aes = FALSE) +
      # Add base line information
      geom_segment(data=base_data, aes(x = start, y = 3, xend = end, yend = 3), colour = "black", alpha=0.8, size=1 , inherit.aes = FALSE )  +
      #geom_text(data=base_data, aes(x = title, y = -6, label=Group), hjust=c(0.75,1.25,0,0), vjust=c(1.2,2,-1,0), colour = "black", alpha=0.8, size=3, fontface="bold", inherit.aes = FALSE)
      # label based on lab.dist 
      geom_text(data=base_data, aes(x = title, y = lab.dist, label=Group), colour = "black", alpha=0.8, size=text.size, fontface="bold", inherit.aes = FALSE) 
    # label based on fixed pos
    #geom_text(data=base_data, aes(x = title, y = 25, label=Group), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
    
    #geom_text(data=base_data, aes(x = title+1, y = -((lab.dist)*2), label=CellTotal), colour = "black", alpha=0.8, size=text.size*0.75, fontface="bold", inherit.aes = FALSE)
    #geom_text(data=base_data, aes(x = title+1, y = lab.dist, label=CellTotal), colour = "black", alpha=0.8, size=text.size*0.75, fontface="bold", inherit.aes = FALSE)
    
    
    
    globe
    
    
  }
}
