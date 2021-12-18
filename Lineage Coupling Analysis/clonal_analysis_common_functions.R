# Load all the necessary packages
packages <- c("Seurat", "ggplot2", "tidyr", "readr"
              , "plyr", "stringr", "harmony"
              , "cowplot", "reshape2","ggpubr"
              , "gsubfn", "tibble", "gplots"
              , "Matrix", "dplyr", "pbapply", "schex"
              , "UpSetR", "extraDistr", "ape"
              , "stats"
              )
lapply(packages, library, character.only = TRUE)



append.cellID = function(seurat.obj){
  # Define how we're making the transformations
  
  # I think there must be a better way of doing this
  cellID.column <- colnames(seurat.obj)
  cellID.column <- substr(cellID.column, 1, nchar(cellID.column)-2)
  # Doesn't convince me. I'm manually appending the trailing '_1' when I suspect it actually comes from the dataset name. Must think how to do it dynamically
  cellID.column <- paste(seurat.obj$orig.ident, cellID.column, sep="_")
  seurat.obj@meta.data$cellID <- cellID.column
  return(seurat.obj)
}

# Left Join way
compute.cols.to.keep = function(dt.source, df.cols.discard, except.cols){
  cols.dt.source <- colnames(dt.source)
  cols.df.cols.discard <- colnames(df.cols.discard)
  cols.to.discard <- cols.df.cols.discard[!(cols.df.cols.discard %in% except.cols)]
  cols.to.keep <- cols.dt.source[!(cols.dt.source %in% cols.to.discard)]
  return(cols.to.keep)
}


merge.transcr.ling <- function(seurat.obj, pool, clusters.column = "seurat_clusters"
                               , only.cells.from.multicell.clones = FALSE){
  seurat.obj <- append.cellID(seurat.obj)
  
  library(data.table)
  seur.obj.dt <- data.table(seurat.obj@meta.data)
  pool.dt <- data.table(pool)
  
  # Should I use a data table?
  # There must be a way to join and select only the columns that I'm interested in
  # I think maybe editing removing the columns of pool I'm not interested in beforehand could be the way
  merged.data.dt <- merge(seur.obj.dt, pool, by="cellID", all.x=TRUE)
  
  # Discard rows that didn't have a cloneID, and all columns I'm not interested in
  # This would fail if there are common column names in both -pools and seurat.obj- data frames; need to manually add the exceptions
  columns.to.keep.for.seurat.obj <- compute.cols.to.keep(merged.data.dt, pool, c("cloneID"))
  columns.to.keep.for.pool <- compute.cols.to.keep(merged.data.dt, seurat.obj@meta.data, c(clusters.column, "cellID"))
  
  seurat.obj@meta.data <- setDF(merged.data.dt[, ..columns.to.keep.for.seurat.obj])
  rownames(seurat.obj@meta.data) <- colnames(seurat.obj)
  if(only.cells.from.multicell.clones){
    multicell.clonesID <- seurat.obj@meta.data %>% ungroup() %>% group_by(cloneID) %>% filter(n() != 1)
    multicell.clonesID <- unique(multicell.clonesID$cloneID)
    seurat.obj@meta.data <- seurat.obj@meta.data[(seurat.obj@meta.data$cloneID %in% multicell.clonesID) & (!is.na(seurat.obj@meta.data$cloneID)), ]
  }
  
  pool <- setDF(merged.data.dt[!is.na(cloneID), ..columns.to.keep.for.pool])
  rownames(pool) <- pool$cellID
  
  pool <- pool %>% dplyr::rename(ident = .data[[clusters.column]])
  pool <- pool[,c("dataset", "cell", "cloneID", "cellID", "ident")]
  pool <- pool[complete.cases(pool), ] # REMOVE everything that includes NA 
  pool <- pool %>% ungroup() %>% group_by(cloneID) %>% filter(n() != 1) # Take everything that is not single clone (from the cells in pool, not the cells in the transcriptome dataset)
  
  return(list(seurat.obj, pool))
}


merge.transcr.ling.2 <- function(seurat.obj, pool, only.cells.from.multicell.clones = FALSE){
  # Focused on return pool, and also change the column name of the seurat clusters to use
  seurat.obj <- append.cellID(seurat.obj)
  
  library(data.table)
  seur.obj.dt <- data.table(seurat.obj@meta.data)
  pool.dt <- data.table(pool)
  
  # Should I use a data table?
  # There must be a way to join and select only the columns that I'm interested in
  # I think maybe editing removing the columns of pool I'm not interested in beforehand could be the way
  merged.data.dt <- merge(seur.obj.dt, pool, by="cellID", all.x=TRUE)
  
  # Discard rows that didn't have a cloneID, and all columns I'm not interested in
  # This would fail if there are common column names in both -pools and seurat.obj- data frames; need to manually add the exceptions
  columns.to.keep.for.pool <- compute.cols.to.keep(merged.data.dt, seurat.obj@meta.data, c("seurat_clusters_z_score", "cellID"))
  
  pool <- setDF(merged.data.dt[!is.na(cloneID), ..columns.to.keep.for.pool])
  rownames(pool) <- pool$cellID
  
  pool <- pool %>% dplyr::rename(ident = seurat_clusters_z_score)
  pool <- pool[,c("dataset", "cell", "cloneID", "cellID", "ident")]
  pool <- pool[complete.cases(pool), ] # REMOVE everything that includes NA 
  pool <- pool %>% ungroup() %>% group_by(cloneID) %>% filter(n() != 1) %>% ungroup() # Take everything that is not single clone (from the cells in pool, not the cells in the transcriptome dataset)
  
  return(pool)
}

merge.transcr.ling.proj.ca <- function(seurat.obj, pool, clusters.column="HARMONY_Clusters", only.cells.from.multicell.clones = FALSE){
  library(data.table)
  pool$dataset <- NULL
  pool$cell <- NULL
  seurat.obj@meta.data$cloneID <- NULL
  seurat.obj@meta.data$cellID <- colnames(x = seurat.obj)
  seurat.obj.dt <- data.table(seurat.obj@meta.data)
  pool.dt <- data.table(pool)
  
  # Should I use a data table?
  # There must be a way to join and select only the columns that I'm interested in
  # I think maybe editing removing the columns of pool I'm not interested in beforehand could be the way
  merged.data.dt <- merge(seurat.obj.dt, pool, by="cellID", all.x=TRUE)
  
  # Discard rows that didn't have a cloneID, and all columns I'm not interested in
  # This would fail if there are common column names in both -pools and seurat.obj- data frames; need to manually add the exceptions
  columns.to.keep.for.seurat.obj <- compute.cols.to.keep(merged.data.dt, pool, c("cloneID"))
  columns.to.keep.for.pool <- compute.cols.to.keep(merged.data.dt, seurat.obj@meta.data, c(clusters.column,"cellID","dataset","stage","collected","HARMONY_Assignments","HARMONY_Assignment","HARMONY_Assignment2"))
  
  seurat.obj@meta.data <- setDF(merged.data.dt[, ..columns.to.keep.for.seurat.obj])
  rownames(seurat.obj@meta.data) <- colnames(seurat.obj)
  if(only.cells.from.multicell.clones){
    multicell.clonesID <- seurat.obj@meta.data %>% ungroup() %>% group_by(cloneID) %>% filter(n() != 1)
    multicell.clonesID <- unique(multicell.clonesID$cloneID)
    seurat.obj@meta.data <- seurat.obj@meta.data[(seurat.obj@meta.data$cloneID %in% multicell.clonesID) & (!is.na(seurat.obj@meta.data$cloneID)), ]
  }
  
  pool <- setDF(merged.data.dt[!is.na(cloneID), ..columns.to.keep.for.pool])
  rownames(pool) <- pool$cellID
  
  #pool <- pool %>% dplyr::rename(ident = seurat_clusters)
  #pool <- pool[,c("dataset", "cell", "cloneID", "cellID", "ident")]
  pool <- pool[complete.cases(pool), ] # REMOVE everything that includes NA 
  pool <- pool %>% ungroup() %>% group_by(cloneID) %>% filter(n() != 1) # Take everything that is not single clone (from the cells in pool, not the cells in the transcriptome dataset)
  
  return(list(seurat.obj, pool))
}

compute.and.plot.freq.table.of.freq <- function(cloneIDs, min.clonesize=NULL
                                                , plot.bool=FALSE
                                                , title="Histogram of clonesizes"
                                                , n.breaks=NULL, limits.x=NULL
                                                , limits.y=NULL
                                                , by=NULL, rel=FALSE){
                                                # , breaks=NULL, limits.x=NULL
                                                # , num.ticks=NULL, by=NULL)){
  # Histogram and table of clonesizes
  cloneIDs <- cloneIDs[!is.na(cloneIDs)]
  cloneIDs.table <- table(cloneIDs)
  clonesizes.table <- table(cloneIDs.table)
  if(!is.null(min.clonesize)){
    cloneIDs.table <- cloneIDs.table[cloneIDs.table>=min.clonesize]
    clonesizes.table <- clonesizes.table[as.numeric(rownames(clonesizes.table))>=min.clonesize]
  }
  # if(is.null(limits.x)){
  #   # For some reason the limits are apparently triming 2 values
  #   limits.x <- c(0, max(cloneIDs.table)+(5-(max(cloneIDs.table)%%5))%%5)
  #   print(limits.x)
  # }
  # if(is.null(breaks)){
  #   # by has priority over num.ticks
  #   if(is.null(by)){
  #     if(is.null(num.ticks)){
  #       num.ticks <- 4
  #     }
  #     breaks <- seq.int(from=limits.x[1], to=limits.x[2], length.out=num.ticks)
  #   }
  #   else{
  #     breaks <- seq.int(from=limits.x[1], to=limits.x[2], by=by)
  #   }
  # }
  if(plot.bool){
    g <- ggplot(as.data.frame(cloneIDs.table), aes(x=Freq))
    if(rel){
      g <- (g
            + geom_histogram(aes(y = stat(count) / sum(count)), binwidth=1, color="black", fill="white")
            )
    }
    else{
      g <- g + geom_histogram(binwidth=1, color="black", fill="white")
    }
    g <- (g +  scale_x_continuous(name="clonesize"
                                 #, breaks = breaks
                                 , limits=limits.x)
            + scale_y_continuous(name="Number of clones"
                                 , limits=limits.y)
            + ggtitle(title))
    print(g)
  }
  return(list(clonesizes.table, g))
}


plot.seur.obj.color.all.clones <- function(seurat.obj){
  color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  col=sample(color, length(table(seurat.obj@meta.data$cloneID)), replace = T) # Why table()?
  print(DimPlot(seurat.obj, reduction = "umap", pt.size = 0.8, group.by = "cloneID",
          label = F, label.size = 15) + scale_colour_manual(values = col, na.value = "grey98")+
    NoLegend())
}

plot.seur.obj.color.selected.clones <- function(seurat.obj, clones, same.color=TRUE
                                                , print.bool=TRUE){
  # Append a new column with only the cloneIDs in clones, and make it the new identity class for plotting it
  cells.highlight <- as.list(rownames(seurat.obj@meta.data[seurat.obj@meta.data$cloneID %in% clones, ]))
  cloneIDs.table <- table(seurat.obj@meta.data$cloneID)
  selected.clones.table <- cloneIDs.table[as.character(clones)]
  if(same.color){
    cols.highlight <- "red" # Why table()?
  }
  else{
    color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    cols.highlight <- sample(color, length(clones), replace = F) # Why table()?
    cols.highlight <- as.list(rep.int(cols.highlight, selected.clones.table))
  }
  g <- (DimPlot(seurat.obj, reduction = "umap", pt.size = 0.8,
          label = F, label.size = 15,
          # For some reason I can't get it to plot the different clones with different colors
          cells.highlight = cells.highlight, cols.highlight = cols.highlight,
          sizes.highlight = 2.0) + NoLegend())
  if(print.bool){
    print(g)
  }
  return(g)
}

plot.seur.obj.color.selected.cells.by.clone <- function(seurat.obj, cells, same.color=FALSE
                                                        , print.bool=TRUE){
  # Append a new column with only the cloneIDs in clones, and make it the new identity class for plotting it
  cloneIDs.table <- table(seurat.obj@meta.data[cells, "cloneID"])
  cells <- as.list(cells)
  if(same.color){
    cols.highlight <- "red" # Why table()?
  }
  else{
    color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    cols.highlight <- sample(color, length(cloneIDs.table), replace = F) # Why table()?
    cols.highlight <- as.list(rep.int(cols.highlight, cloneIDs.table))
  }
  g <- (DimPlot(seurat.obj, reduction = "umap", pt.size = 0.8,
          label = F, label.size = 15,
          # For some reason I can't get it to plot the different clones with different colors
          cells.highlight = cells, cols.highlight = cols.highlight,
          sizes.highlight = 2.0) + NoLegend())
  if(print.bool){
    print(g)
  }
  return(g)
}

plot.seur.obj.color.clones.from.clusters <- function(seurat.obj, clusters
                                                     , same.color=FALSE
                                                     , only.multicell.clones=TRUE
                                                     , only.cells.from.clusters=TRUE
                                                     , print.bool=TRUE){
  # Plot and highlight cells from clones in clusters
  clones.of.clusters <- NULL
  idxs <- seurat.obj@meta.data$seurat_clusters %in% clusters
  if(only.multicell.clones){
    multicell.clonesID <- seurat.obj@meta.data %>% filter(!is.na(cloneID)) %>% group_by(cloneID) %>% filter(n() > 1)
    multicell.clonesID <- unique(multicell.clonesID$cloneID)
    idxs <- idxs & (seurat.obj@meta.data$cloneID %in% multicell.clonesID)
  }
  if(only.cells.from.clusters){
    cells.to.highlight <- rownames(seurat.obj@meta.data[idxs, ])
    return(
      plot.seur.obj.color.selected.cells.by.clone(
        seurat.obj, cells.to.highlight, same.color=FALSE)
    )
  }
  clones.of.clusters <- seurat.obj@meta.data[idxs, "cloneID"]
  
  clones.of.clusters <- union(
    NULL, clones.of.clusters
  )
  g <- plot.seur.obj.color.selected.clones(seurat.obj, clones.of.clusters, same.color)
  if(print.bool){
    print(g)
  }
  return(g)
}

plot.seur.obj.color.shared.clones.between.clusters <- function(seurat.obj, clusters
                                                               , same.color=FALSE
                                                               , only.cells.from.clusters=TRUE
                                                               , print.bool=TRUE){
  seurat.obj.meta.data.cells.from.shared.clones <- seurat.obj@meta.data
  seurat.obj.meta.data.cells.from.shared.clones[, "cellID"] <- rownames(
                            seurat.obj.meta.data.cells.from.shared.clones)
  seurat.obj.meta.data.cells.from.shared.clones <- (seurat.obj.meta.data.cells.from.shared.clones
                                                    %>% filter(!is.na(cloneID))
                                                    %>% group_by(cloneID)
                                                    %>% filter(Reduce("&", clusters %in% unique(seurat_clusters)))
                                                    %>% ungroup())
  shared.clonesID <- unique(seurat.obj.meta.data.cells.from.shared.clones$cloneID)
  
  shared.clonesID.table <- table(seurat.obj@meta.data[seurat.obj@meta.data$cloneID %in% shared.clonesID, "cloneID"])
  print("shared.clonesID.table:")
  print(shared.clonesID.table)
  print(paste0("length(shared.clonesID) = ", length(shared.clonesID)))
  if(only.cells.from.clusters){
    seurat.obj.meta.data.cells.from.shared.clones <- (seurat.obj.meta.data.cells.from.shared.clones
                                                      %>% filter(seurat_clusters %in% clusters))
    cells.to.highlight <- seurat.obj.meta.data.cells.from.shared.clones$cellID
    g <- plot.seur.obj.color.selected.cells.by.clone(
              seurat.obj, cells.to.highlight, same.color=FALSE, print.bool=FALSE)
  }
  else {
    g <- plot.seur.obj.color.selected.clones(seurat.obj, shared.clonesID, same.color, print.bool=FALSE)
    
  }
  g <- g + ggtitle(label="Dispersion of shared clones"
                   , subtitle=paste("Clusters", paste(clusters
                                                      , collapse=", ")
                                    , sep=": "))
  if(print.bool){
    print(g)
  }
  return(g)
}

plot.hexbin.all.clones = function(seurat.obj, nbins=150){
  # Sets a new column, which is basically a copy of cloneID but where cells that don't have a clone have a value of 0
  seurat.obj@meta.data[,"cloneID"] <- ifelse(is.na(seurat.obj@meta.data[,"cloneID"]), 0, seurat.obj@meta.data[,"cloneID"])
  seurat.obj <- make_hexbin(seurat.obj, nbins=nbins, dimension_reduction = "UMAP")
  # print(plot_hexbin_density(seurat.obj))
  # readline(prompt = "Press enter to continue...")
  gg <- plot_hexbin_meta(seurat.obj, col = "cloneID", action = "prop_0")
  print(gg + theme_void())
}

plot.hexbin.selected.clones = function(seurat.obj, clones, nbins=150){# Sets a new column, which is basically a copy of cloneID but where cells that are not in clones have a value of 0
  seurat.obj@meta.data[,"cloneID"] <- ifelse((seurat.obj@meta.data[,"cloneID"] %in% clones), seurat.obj@meta.data[,"cloneID"], 0)
  seurat.obj <- make_hexbin(seurat.obj, nbins=nbins, dimension_reduction = "UMAP")
  # plot_hexbin_density(seurat.obj)
  # readline(prompt = "Press enter to continue...")
  gg <- plot_hexbin_meta(seurat.obj, col = "cloneID", action = "prop_0")
  print(gg + theme_void())
}

plot.hexbin.clones.of.sizes = function(seurat.obj, clonesizes, nbins=150){# Sets a new column, which is basically a copy of cloneID but where cells that are not in clones have a value of 0
  clones.tibble <- seurat.obj@meta.data %>% filter(!is.na(cloneID)) %>% group_by(cloneID) %>% summarize(count=n()) %>% filter(count %in% clonesizes)
  clones <- clones.tibble$cloneID
  plot.hexbin.selected.clones(seurat.obj, clones, nbins=nbins)
}

plot.hexbin.clones.of.size.at.least = function(seurat.obj, min.clonesize, nbins=150){# Sets a new column, which is basically a copy of cloneID but where cells that are not in clones have a value of 0
  clones.tibble <- seurat.obj@meta.data %>% filter(!is.na(cloneID)) %>% group_by(cloneID) %>% summarize(count=n()) %>% filter(count >= min.clonesize)
  clones <- clones.tibble$cloneID
  plot.hexbin.selected.clones(seurat.obj, clones, nbins=nbins)
}

plot.hexbin.cells.matched.values = function(seurat.obj, columns.values=list("cloneID"=NULL)
                                            , subset.not.nulls=FALSE, nbins=150){
  # columns.values must be a named list which stores the columns and values
  # a cell must have in order to take it into account
  # , or NULL if only interested in presence/absence of a value
  # for now I only consider the case when the first column is "cloneID", its value NULL
  # , and there's only 1 other variable, only 1 of whose values I'm interested in
  
  columns.null <- names(columns.values[sapply(columns.values, is.null)])
  columns.not.null <- names(columns.values[sapply(columns.values, function(x){!is.null(x)})])
  # cells.fulfilled.conditions.bool <- (!is.na(seurat.obj@meta.data[, columns.null])
  #                     & seurat.obj@meta.data[, columns.not.null]==unlist(columns.values[columns.not.null]))
  rows.matched.column.values.bool <- as.numeric(seurat.obj@meta.data[, columns.not.null]
                                                ==unlist(columns.values[columns.not.null]))
  # Sets a new column, which is basically a copy of cloneID but where cells that don't have a clone have a value of 0
  seurat.obj@meta.data[, "value"] <- ifelse(is.na(seurat.obj@meta.data[, columns.null])
                                            , NA
                                            , rows.matched.column.values.bool)
  if(subset.not.nulls){
    seurat.obj <- subset(x = seurat.obj, subset = (value==0 | value==1))
  }
  seurat.obj <- make_hexbin(seurat.obj, nbins=nbins, dimension_reduction = "UMAP")
  # print(plot_hexbin_density(seurat.obj))
  # readline(prompt = "Press enter to continue...")
  gg <- plot_hexbin_meta(seurat.obj, col = "value", action = "prop_0")
  print(gg + theme_void())
  return(gg)
}

compute.tibble.clones.cells.per.cluster = function(seurat.obj, only.cells.w.cloneID=FALSE){
  if(only.cells.w.cloneID){
    # Tibble with number of different clones and cells (which had a cloneID) per cluster, ordered by the first
    seurat.obj@meta.data %>% filter(!is.na(cloneID)) %>% group_by(seurat_clusters) %>% summarize(num.diff.clones=length(unique(cloneID)), num.cells=n()) %>% arrange(desc(num.diff.clones))
  }
  else{
    # Tibble with number of different clones and cells per cluster, ordered by the first
    seurat.obj@meta.data %>% group_by(seurat_clusters) %>% summarize(num.diff.clones=n_distinct(cloneID, na.rm=TRUE), num.cells=n()) %>% arrange(desc(num.diff.clones))
  }
}

compute.tibble.num.groups.and.cells.per.clone = function(
        seurat.obj, min.clonesize=NULL, max.clonesize=NULL
         , groups.column="seurat_clusters"){
  if(is.null(min.clonesize) || min.clonesize==0){
    # Tibble with number of different groups and cells per clone, ordered by the first
    seur.obj.grouped.cloneID <- seurat.obj@meta.data %>% group_by(cloneID)
    if(!is.null(max.clonesize)){
      seur.obj.grouped.cloneID <- seur.obj.grouped.cloneID %>% mutate(count=n()) %>% filter(count<=max.clonesize)
    }
    seur.obj.grouped.cloneID <- seur.obj.grouped.cloneID %>% summarize(num.diff.groups=n_distinct(.data[[groups.column]], na.rm=TRUE), num.cells=n()) %>% arrange(desc(num.diff.groups))
  }
  else{
    # Tibble with number of different groups and cells (which had a cloneID) per clone, ordered by the first
    seur.obj.grouped.cloneID <- seurat.obj@meta.data %>% filter(!is.na(cloneID)) %>% group_by(cloneID)
    seur.obj.grouped.cloneID <- seur.obj.grouped.cloneID %>% mutate(count=n()) %>% filter(count>=min.clonesize)
    if(!is.null(max.clonesize)){
      seur.obj.grouped.cloneID <- seur.obj.grouped.cloneID %>% filter(count<=max.clonesize)
    }
    seur.obj.grouped.cloneID <- seur.obj.grouped.cloneID %>% summarize(num.diff.groups=length(unique(.data[[groups.column]])), num.cells=n()) %>% arrange(desc(num.diff.groups))
  }
}

# Tibble with empirical average clonesize per cluster, ordered by it
compute.tibble.empiric.avg.clonesize.per.cluster = function(seurat.obj){
  return(seurat.obj@meta.data %>% filter(!is.na(cloneID)) %>% group_by(seurat_clusters, cloneID) %>% summarize(clonesize=n()) %>% summarize(empiric.avg.clonesize=round(mean(clonesize), digits=1)) %>% arrange(desc(empiric.avg.clonesize)))
}

cv = function(x){
  return(sd(x)/mean(x))
}

append.lineage.barcode.bool.col <- function(seur.obj){
  seur.obj@meta.data$lineage.barcode.bool <- !is.na(seur.obj@meta.data$cloneID)
  return(seur.obj)
}

compute.num.shared.clones.matrix <- function(lineage_dataset){
  # Pending...
  lineage_dataset %>% dplyr::group_by(cloneID) %>% summarise(cells.of.clone.in.ident=n())
}

map.clusters.idxs <- function(array.clusters){
  array.clusters.idxs <- list()
  array.clusters.idxs[array.clusters] <- 1:length(array.clusters)
  return(array.clusters.idxs)
}

# transform.dispersion.motif <- function(orig.dispersion.motif, array.clusters.idxs){
#   # orig.dispersion.motif.table <- table(factor(orig.dispersion.motif, levels=as.numeric(array.clusters)))
#   orig.dispersion.motif <- array.clusters.idxs
#   transformed.dispersion.motif <- 
# }





plot.n.bigger.clones <- function(seurat.obj, clones.tibble=NULL, n=5){
  # Select top n clones (in terms of size) and plot them
  if(is.null(clones.tibble)){
    clones.tibble <- seur.obj@meta.data %>% filter(!is.na(cloneID)) %>% group_by(cloneID) %>% summarize(clonesize=n()) %>% arrange(desc(clonesize))
  }
  n.bigger.clones <- clones.tibble$cloneID[1:n]
  plot.seur.obj.color.selected.clones(seur.obj.integr.w.lb, n.bigger.clones, same.color=F)
  # View(seur.obj@meta.data %>% filter(cloneID %in% n.bigger.clones))
}


perform.computations.clones.cells.per.cluster <- function(seur.obj.integr.w.lb
                                                             , only.cells.w.cloneID=TRUE){
  # Compute some metrics regarding the number of cells/clones per cluster
  # Tibble with number of different clones and cells per cluster, ordered by the first
  tibble.clones.cells.per.cluster <- compute.tibble.clones.cells.per.cluster(seur.obj.integr.w.lb, only.cells.w.cloneID=only.cells.w.cloneID)
  
  tibble.clones.cells.per.cluster[,"num.diff.clones.pctg"] <- round(tibble.clones.cells.per.cluster[,"num.diff.clones"] / sum(tibble.clones.cells.per.cluster[,"num.diff.clones"]) * 100, digits=2)
  tibble.clones.cells.per.cluster[,"num.cells.pctg"] <- round(tibble.clones.cells.per.cluster[,"num.cells"] / sum(tibble.clones.cells.per.cluster[,"num.cells"]) * 100, digits=2)
  tibble.clones.cells.per.cluster <- tibble.clones.cells.per.cluster[, c("seurat_clusters", "num.diff.clones", "num.diff.clones.pctg", "num.cells", "num.cells.pctg")]
  tibble.clones.cells.per.cluster[,"avg.clonesize"] <- round(tibble.clones.cells.per.cluster[,"num.cells"] / tibble.clones.cells.per.cluster[,"num.diff.clones"], digits=1)
  
  return(tibble.clones.cells.per.cluster)
}

compute.print.descr.stats.tibble.clones.cell.per.cluster <- function(tibble.clones.cells.per.cluster){
  # Compute and plot some descriptive statistics on tibble.clones.cell.per.cluster
  cor(tibble.clones.cells.per.cluster$avg.clonesize, tibble.clones.cells.per.cluster$num.cells)
  print(plot(
        tibble.clones.cells.per.cluster$num.cells, tibble.clones.cells.per.cluster$num.diff.clones
        , main="Number of different clones in a cluster vs. its number of cells"
        , xlab="Number of cells", ylab="Number of different clones"
      ))
  print(ggplot(tibble.clones.cells.per.cluster, aes(x=num.cells, y=num.diff.clones))
    + geom_point()
    + geom_smooth(method="lm", se=FALSE, fullrange=FALSE, level=0.95)
    + stat_cor(method = "pearson")
    + labs(title="Number of different clones in a cluster vs. its number of cells",
           x="Number of cells", y = "Number of different clones"))
  cor(tibble.clones.cells.per.cluster$num.diff.clones, tibble.clones.cells.per.cluster$num.cells)
  sapply(tibble.clones.cells.per.cluster[,-c(1)], mean)
  sapply(tibble.clones.cells.per.cluster[,-c(1)], stats::sd)
  sapply(tibble.clones.cells.per.cluster[,-c(1)], cv)
}
  
perform.computations.num.groups.and.cells.per.clone <- function(seur.obj.integr.w.lb
                                                          , min.clonesize=NULL
                                                          , max.clonesize=NULL
                                                          , groups.column="seurat_clusters"){
  # Compute some metrics regarding the number of group/cells per clones
  # Tibble with number of different groups and cells (which had a cloneID) per clone, ordered by the first
  tibble.groups.cells.per.clone <- compute.tibble.num.groups.and.cells.per.clone(
        seur.obj.integr.w.lb, min.clonesize=min.clonesize, max.clonesize=max.clonesize
        , groups.column=groups.column)
  
  tibble.groups.cells.per.clone[,"num.diff.groups.pctg"] <- round(tibble.groups.cells.per.clone[,"num.diff.groups"] / sum(tibble.groups.cells.per.clone[,"num.diff.groups"]) * 100, digits=2)
  tibble.groups.cells.per.clone[,"num.cells.pctg"] <- round(tibble.groups.cells.per.clone[,"num.cells"] / sum(tibble.groups.cells.per.clone[,"num.cells"]) * 100, digits=2)
  tibble.groups.cells.per.clone <- tibble.groups.cells.per.clone[, c("cloneID", "num.diff.groups", "num.diff.groups.pctg", "num.cells", "num.cells.pctg")]
  tibble.groups.cells.per.clone[,"avg.cells.per.group"] <- round(tibble.groups.cells.per.clone[,"num.cells"] / tibble.groups.cells.per.clone[,"num.diff.groups"], digits=1)
  
  return(tibble.groups.cells.per.clone)
}

perform.computations.clusters.cells.per.clone <- function(seur.obj.integr.w.lb
                                                          , min.clonesize=NULL
                                                          , max.clonesize=NULL
                                                          , groups.column="seurat_clusters"){
  # Function kept only for backwards compatibility purposes
  # (in some files it's still referenced by this name)
  perform.computations.num.groups.and.cells.per.clone(seur.obj.integr.w.lb
                                                      , min.clonesize=min.clonesize
                                                      , max.clonesize=min.clonesize
                                                      , groups.column=groups.column)
}

compute.print.descr.stats.tibble.clusters.cells.per.clone <- function(tibble.clusters.cells.per.clone){
  # Compute and plot some descriptive statistics on tibble.clusters.cells.per.clone
  print(cor(tibble.clusters.cells.per.clone$avg.cells.per.cluster, tibble.clusters.cells.per.clone$num.cells))
  print(cor(tibble.clusters.cells.per.clone$num.diff.clusters, tibble.clusters.cells.per.clone$num.cells))
  print(sapply(tibble.clusters.cells.per.clone[,-c(1)], mean))
  print(sapply(tibble.clusters.cells.per.clone[,-c(1)], stats::sd))
  print(sapply(tibble.clusters.cells.per.clone[,-c(1)], cv))
}

compute.plot.num.dif.groups.freq <- function(tibble.groups.cells.per.clone
                                             , groups.column="refined_COUP_class"){
  # Frequencies of numbers of different groups
  num.diff.groups.table <- table(tibble.groups.cells.per.clone$num.diff.groups)
  g <- (ggplot(tibble.groups.cells.per.clone, aes(x=num.diff.groups))
        + geom_histogram(binwidth=1, color="black", fill="white")
        +  scale_x_continuous(name = paste0("Number of different '"
                                            , groups.column, "' values"))
                              #, breaks = breaks)
        + scale_y_continuous(name = "Number of clones")
        + ggtitle(paste0("Histogram of different '"
                         , groups.column, "' values within clone"))
        )
  print(g)
  # print(hist(tibble.groups.cells.per.clone$num.diff.groups
  #      , breaks = (0:max(as.numeric(names(num.diff.groups.table))) + 0.5)
  #      , freq = TRUE, labels = FALSE
  #      , main = paste0("Histogram of different '"
  #                        , groups.column, "' values within clone")
  #      , xlab = paste0("Number of different '", groups.column, "' values")
  #      , ylab = "Number of clones"))
  View(num.diff.groups.table)
  return(g)
}

compute.plot.num.dif.groups.rel.freq <- function(tibble.groups.cells.per.clone
                                                 , groups.column="refined_COUP_class"){
  # Frequencies of numbers of different groups
  # rel.num.diff.groups.table <- table(tibble.groups.cells.per.clone$num.diff.groups)
  g <- (ggplot(tibble.groups.cells.per.clone, aes(x=num.diff.groups))
        + geom_histogram(aes(y = stat(count) / sum(count)), binwidth=1, color="black", fill="white")
        +  scale_x_continuous(name = paste0("Number of different '"
                                            , groups.column, "' values"))
                              #, breaks = breaks)
        + scale_y_continuous(name = "Number of clones", labels = scales::percent)
        + ggtitle(paste0("Histogram of different '"
                         , groups.column, "' values within clone"))
        )
  print(g)
  # print(hist(tibble.groups.cells.per.clone$num.diff.groups
  #      , breaks = (0:max(as.numeric(names(rel.num.diff.groups.table))) + 0.5)
  #      , freq = TRUE, labels = FALSE
  #      , main = paste0("Histogram of different '"
  #                        , groups.column, "' values within clone")
  #      , xlab = paste0("Number of different '", groups.column, "' values")
  #      , ylab = "Number of clones"))
  
  # View(rel.num.diff.groups.table)   Replace with rel.freq
  
  return(g)
}

compute.plot.avg.num.cells.per.cluster.within.clone.freq <- function(tibble.clusters.cells.per.clone){
  # Frequencies of average numbers of cells per cluster, for the different clones
  avg.cells.per.cluster.table <- table(tibble.clusters.cells.per.clone$avg.cells.per.cluster)
  print(hist(tibble.clusters.cells.per.clone$avg.cells.per.cluster
       , breaks = (0:ceiling(max(as.numeric(names(avg.cells.per.cluster.table)))) + 0.5)
       , freq = TRUE, labels = FALSE
       , main = "Histogram of average cells per cluster within clones"
       , xlab = "Avg. cells per cluster", ylab = "Number of clones"))
  View(avg.cells.per.cluster.table)
}

compute.dispersion.motifs.freq <- function(seur.obj.integr.w.lb, min.clonesize=2
                                           , groups.column="seurat_clusters"
                                           , multiset=TRUE){
  # Frequencies of groups of clusters among which the clones dispersed
  if(multiset){
    # clusters.of.cells.of.clones <- seur.obj.integr.w.lb@meta.data %>% group_by(cloneID) %>% summarise(diff.clusters=paste(sort(.data[[groups.column]]), collapse=","))
    clusters.of.cells.of.clones <- seur.obj.integr.w.lb@meta.data %>%
            group_by(cloneID) %>% filter(n()>=min.clonesize) %>%
            summarise(diff.clusters=paste(sort(.data[[groups.column]]), collapse=","))
    
  }
  else{
    # seur.obj.integr.w.lb@meta.data %>% group_by(cloneID) %>% summarise(diff.clusters=c(unique(.data[[groups.column]])))
    # Without counting repeated clusters
    clusters.of.cells.of.clones <- seur.obj.integr.w.lb@meta.data %>%
            group_by(cloneID) %>% filter(n()>=min.clonesize) %>%
            summarise(diff.clusters=paste(sort(unique(.data[[groups.column]])), collapse=","))
    
  }
  # clusters.of.cells.of.clones.split.into.columns <- clusters.of.cells.of.clones %>% separate(col = diff.clusters, into = paste("ident", 1:16, sep = "_"))
  clusters.of.cells.of.clones.table <- table(clusters.of.cells.of.clones$diff.clusters)
  # print(barplot2(clusters.of.cells.of.clones.table, main = "Histogram of groups of clusters within clones", xlab = "Cluster groups", ylab = "Number of clones"))
  View(sort(clusters.of.cells.of.clones.table, decreasing = T))
  return(clusters.of.cells.of.clones.table)
}

sum.num.clusters.clones <- function(dispers.motif.clusters, num.clones
                                     , dispers.clones.per.cluster){
  # Count the clone as how many different clusters it dispersed into, into the corresponding clusters
  na.clusters.bool <- is.na(dispers.clones.per.cluster[
              dispers.motif.clusters, as.character(length(dispers.motif.clusters))
              ]
          )
  if(any(na.clusters.bool)){
    na.clusters <- dispers.clones.per.cluster[na.clusters.bool]
    dispers.clones.per.cluster[
        na.clusters] <- 0
  }
  dispers.clones.per.cluster[
      dispers.motif.clusters, as.character(length(dispers.motif.clusters))] <- (
          dispers.clones.per.cluster[
              dispers.motif.clusters, as.character(length(dispers.motif.clusters))]
        + num.clones)
  return(dispers.clones.per.cluster)
}

compute.num.dispers.clones.per.cluster <- function(dispersion.motifs.freq){
  dispers.clones.per.cluster <- data.frame(row.names = levels(clusters.factor))
  dispers.clones.per.cluster$cluster <- rownames(dispers.clones.per.cluster)
  # Por alguna fucking razon parece que num.dispers.clones termina con 1 clon extra
  # Esto es un fix improvisadisimo que no me gusta
  for(dispersion.motif in names(dispersion.motifs.freq)){
    dispers.motif.clusters <- as.array(unlist(strsplit(dispersion.motif, ",")))
    dispers.clones.per.cluster <- sum.num.clusters.clones(
        dispers.motif.clusters
        , dispersion.motifs.freq[dispersion.motif]
        , dispers.clones.per.cluster)
  }
  cols.drop <- c("clusters")
  dispers.clones.per.cluster$num.clones <- apply(
      dispers.clones.per.cluster[
          , !(colnames(dispers.clones.per.cluster) %in% cols.drop), drop=TRUE
          ]
      , 1, sum
      )
  # dispers.clones.per.cluster$pctg.num.restr.clones <- (
  #   dispers.clones.per.cluster$num.restr.clones
  #   / dispers.clones.per.cluster$num.clones)
  # dispers.clones.per.cluster$pctg.num.dispers.clones <- (
  #   dispers.clones.per.cluster$num.dispers.clones
  #   / dispers.clones.per.cluster$num.clones)
  return(dispers.clones.per.cluster)
}

sum.restr.dispers.clones <- function(dispers.motif.clusters, num.clones
                                      , restr.dispers.clones.per.cluster){
  # Count the clone as restricted or dispersed into each corresponding clusters
  if(length(dispers.motif.clusters)==1){
    column.to.add.to <- "num.restr.clones"
  } else{
    column.to.add.to <- "num.dispers.clones"
  }
  restr.dispers.clones.per.cluster[dispers.motif.clusters, column.to.add.to] <- (
          restr.dispers.clones.per.cluster[dispers.motif.clusters, column.to.add.to]
          + num.clones)
  return(restr.dispers.clones.per.cluster)
}

compute.num.restr.dispers.clones.per.cluster <- function(dispersion.motifs.freq){
  restr.dispers.clones.per.cluster <- data.frame(row.names = levels(clusters.factor))
  restr.dispers.clones.per.cluster$cluster <- rownames(restr.dispers.clones.per.cluster)
  restr.dispers.clones.per.cluster$num.restr.clones <- 0
  # Por alguna fucking razon parece que num.dispers.clones termina con 1 clon extra
  # Esto es un fix improvisadisimo que no me gusta
  restr.dispers.clones.per.cluster$num.dispers.clones <- -1
  for(dispersion.motif in names(dispersion.motifs.freq)){
    dispers.motif.clusters <- as.array(unlist(strsplit(dispersion.motif, ",")))
    restr.dispers.clones.per.cluster <- sum.restr.dispers.clones(
                                    dispers.motif.clusters
                                    , dispersion.motifs.freq[dispersion.motif]
                                    , restr.dispers.clones.per.cluster)
  }
  restr.dispers.clones.per.cluster$num.clones <- (
      restr.dispers.clones.per.cluster$num.restr.clones
      + restr.dispers.clones.per.cluster$num.dispers.clones)
  restr.dispers.clones.per.cluster$pctg.num.restr.clones <- (
    restr.dispers.clones.per.cluster$num.restr.clones
    / restr.dispers.clones.per.cluster$num.clones)
  restr.dispers.clones.per.cluster$pctg.num.dispers.clones <- (
    restr.dispers.clones.per.cluster$num.dispers.clones
    / restr.dispers.clones.per.cluster$num.clones)
  return(restr.dispers.clones.per.cluster)
}

scatter.plot.pctg.and.identity.func.line <- function(tibble.clones.cells.per.cluster){
  # Scatter plot of % of Number of cells vs. % of number of clones per cluster, and y=x line
  df.to.plot <- data.frame(tibble.clones.cells.per.cluster$num.diff.clones.pctg
                           , tibble.clones.cells.per.cluster$num.diff.clones.pctg
                           , tibble.clones.cells.per.cluster$num.cells.pctg)
  g <- ggplot(df.to.plot, aes(tibble.clones.cells.per.cluster$num.diff.clones.pctg))
  g <- g + geom_line(aes(y=tibble.clones.cells.per.cluster$num.diff.clones.pctg), colour="red")
  g <- g + geom_point(aes(y=tibble.clones.cells.per.cluster$num.cells.pctg), colour="green")
  g <- (g + ylab("Pctg of cells") + xlab("Pctg of clones")
            + ggtitle("Percentage of cells vs. of clones in each cluster"))
  print(g)
  return(g)
}

upset.plot <- function(seur.obj.integr.w.lb, seurat_clusters_z_score
                       , min.clonesize = 2
                       , nintersects = 40, order.by = c("degree", "freq")
                       , decreasing= c(TRUE, TRUE), filename=NULL){
  # Groups of clones present in each cluster
  # # The next line isn't working for some reason
  # # clones.of.cells.of.clusters <- seur.obj.integr.w.lb@meta.data %>% group_by_at(seurat_clusters_z_score) %>% dplyr::summarise(diff.clones=c(base::unique(cloneID)), .groups = "drop_last")
  # If I don't want to avoid duplicate clones
  # clones.of.cells.of.clusters <- seur.obj.integr.w.lb@meta.data %>%
  #           group_by(cloneID) %>% filter(n()>=min.clonesize) %>% ungroup() %>%
  #           group_by_at(seurat_clusters_z_score)
  #           %>% summarise(diff.clones=paste(sort(cloneID), collapse=","))
  clones.of.cells.of.clusters <- seur.obj.integr.w.lb@meta.data %>%
            group_by(cloneID) %>% filter(n()>=min.clonesize) %>% ungroup() %>%
            group_by_at(seurat_clusters_z_score) %>%
            summarise(diff.clones=paste(sort(unique(cloneID)), collapse=","))
  # View(clones.of.cells.of.clusters)
  # # clones.of.cells.of.clusters <- clones.of.cells.of.clusters %>% mutate()
  # num.clones.of.cells.of.clusters <- clones.of.cells.of.clusters %>% group_by_at(seurat_clusters_z_score) %>% mutate(num.diff.clones=length(unlist(strsplit(diff.clones, ",")))) %>% ungroup()
  # View(num.clones.of.cells.of.clusters)
  # # clones.of.cells.of.clusters.split.into.columns <- clones.of.cells.of.clusters %>%
  # #                                                       separate(col = diff.clones
  # #                                                                , into = paste(
  # #                                                                               "cloneID"
  # #                                                                               , 1:max(num.clones.of.cells.of.clusters$num.diff.clones)
  # #                                                                               , sep = "_")
  # #                                                                )
  clones.of.cells.of.clusters.list <- setNames(strsplit(clones.of.cells.of.clusters$diff.clones, ",")
                                               , paste("ident"
                                                       , clones.of.cells.of.clusters[
                                                            , seurat_clusters_z_score, drop=TRUE]
                                                       , sep="_")
                                              )
  clones.of.cells.of.clusters.list <- lapply(clones.of.cells.of.clusters.list, as.array)
  # Just a verification that a random element on the list is indeed an array, and the correct one
  # identical(lapply(clones.of.cells.of.clusters.list, as.array)[['14']]
  #           , as.array((setNames(strsplit(clones.of.cells.of.clusters$diff.clones, ",")
  #                                , clones.of.cells.of.clusters[, seurat_clusters_z_score, drop=TRUE])[['14']]
  #                       )
  #                      )
  #           )
  clones.of.cells.of.clusters.list
  if(!is.null(filename)){
    pdf(file=filename, title="Clonal intersections of Clusters", paper = "legal")
  }
  print(upset(fromList(clones.of.cells.of.clusters.list), sets = names(clones.of.cells.of.clusters.list)
        , nintersects = nintersects
        , mainbar.y.label = "Cluster Intersections", sets.x.label = "Cluster size" 
        , order.by = order.by
        , decreasing = decreasing))
  if(!is.null(filename)){
    dev.off()
  }
}

replace.df.col.non.na.vals <- function(df, column, val=NULL, column.to.copy=NULL){
  if(!is.null(column.to.copy)){
    df[!is.na(df[, column, drop=TRUE]), column] <- df[!is.na(df[, column, drop=TRUE]), column.to.copy]
  }
  if(!is.null(val)){
    df[!is.na(df[, column, drop=TRUE]), column] <- val
  }
  return(df)
}

create.lb.calling.methods.array <- function(
  root="Ryan_toptier", suffix.min=NULL, suffix.max=NULL
  , suffix.array=NULL){
  if(!is.null(suffix.min) && !is.null(suffix.max)){
    suffix.array <- suffix.min:suffix.max
  }
  return(paste(root, as.character(suffix.array), sep=ifelse(is.null(suffix.array),"","_")))
}

plot.stacked.barchart <- function(df, column, fill, position
                                  , reorder.revers=FALSE, reorder.by.freq=FALSE
                                  , coord.flip.bool=TRUE
                                  , colors.bool=FALSE){
  if(reorder.by.freq){
    df[, column] <- reorder(df[, column, drop=TRUE]
                                         , df[, column, drop=TRUE], FUN=length)
  }
  if(reorder.revers){
    df[, column] <- factor(df[, column, drop=TRUE]
                                      , levels=rev(levels(df[, column, drop=TRUE])))
  }
  p <- (ggplot(df, aes(fill=.data[[fill]], x=.data[[column]]))
        + geom_bar(position=position))
  if(coord.flip.bool){
    p <- p + coord_flip()
  }
  if(colors.bool){
    c25 <- c(
      "dodgerblue2", "#E31A1C", # red
      "green4",
      "#6A3D9A", # purple
      "#FF7F00", # orange
      "black", "gold1",
      "skyblue2", "#FB9A99", # lt pink
      "palegreen2",
      "#CAB2D6", # lt purple
      "#FDBF6F", # lt orange
      "gray70", "khaki2",
      "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
      "darkturquoise", "green1", "yellow4", "yellow3",
      "darkorange4", "brown"
    )
    colors <- c25[1:length(unique(df[, fill, drop=TRUE]))]
    p <- p + scale_fill_manual(values=colors)
  }
  print(p)
  # View(df %>% group_by(seurat_clusters, stage) %>% summarise(count=n()) %>% mutate(pctg=count/sum(count)))
  return(p)
}

plot.stat.count <- function(df, x, fill, position="fill", width=0.9, coord.flip.bool=TRUE){
  # it's the same than plot.stacked.barchart but allows customization of width
  p <- (ggplot(df, aes(x=.data[[x]], fill=.data[[fill]]))
         + stat_count(position=position, width=width))
  if(coord.flip.bool){
    p <- p + coord_flip()
  }
  print(p)
  View(df %>% group_by(seurat_clusters, stage) %>% summarise(count=n()) %>% mutate(pctg=count/sum(count)))
  return(p)
}

## function to return n largest values and position for matrix m
# I have to adapt it
n.largest.matrix <- function(m, n, sym=TRUE) {
  mult <- 1;
  if (sym) mult <- 2;
  res <- order(m, decreasing=TRUE)[seq_len(n) * mult];
  pos <- arrayInd(res, dim(m), useNames=TRUE);
  list(values = m[res],
       position = pos)
}

subset.seurat.obj.col.eq.val <- function(seurat.obj, col, val){
  # Sets a new column, which is basically a copy of cloneID but where cells that don't have a clone have a value of 0
  seurat.obj@meta.data[, "subset.bool"] <- ifelse(
                            seurat.obj@meta.data[, col] %in% val
                            , FALSE, TRUE)
  seurat.obj <- subset(x = seurat.obj, subset = subset.bool)
  seurat.obj@meta.data[, "subset.bool"] <- NULL
  return(seurat.obj)
}

prepare.input.file.for.lineage.coupling.analysis <- function(
        seur.obj.proj.ca, min.clonesize=2
        , lb.calling.method="Ryan_toptier_9"
        , stage="", stage.keep.bool=TRUE
        , clusters.column="refined_COUP_clust"
        , use.default.names=TRUE
        , parent.folder.path="./results/"
        , provided.file.name=""
        , dataset.name="clonal_dataset"
        , write.csv.bool=TRUE
        , col.subset=clusters.column
        , val.subset.exclude="Remove"){
  seur.obj.proj.ca@meta.data$cloneID <- seur.obj.proj.ca@meta.data[, lb.calling.method, drop=TRUE]
  
  pool.w.clusters <- seur.obj.proj.ca@meta.data %>%
              filter(!is.na(cloneID)) %>%
              rownames_to_column(var = "cellID") %>%
              group_by(cloneID) %>% filter(n()>=min.clonesize) %>% ungroup()
  
  # If indicated, subset the tibble so that it only has
  # cells that were [not] injected at a specific stage
  stage.bool <- is.character(stage) && stage!=""
  if(stage.bool){
    if(stage.keep.bool){
      pool.w.clusters <- pool.w.clusters %>% filter(stage==!!stage)
    }
    else{
      pool.w.clusters <- pool.w.clusters %>% filter(stage!=!!stage)
    }
  }
  
  subset.seurat.obj.col.eq.val.bool = (is.character(col.subset)
                                       && is.character(val.subset.exclude)
                                       && col.subset!=""
                                       && val.subset.exclude!="")
  print(table(pool.w.clusters[, clusters.column, drop=TRUE]))
  if(subset.seurat.obj.col.eq.val.bool){
    # No I don't actually subset the seurat object because it doesn't make a difference
    # and is faster this way
    pool.w.clusters <- pool.w.clusters %>% filter(.data[[col.subset]]!=!!val.subset.exclude)
  }
  
  table(pool.w.clusters[, clusters.column, drop=TRUE])
  # write.csv(table(pool.w.clusters$ident),file="results/QC/HARMONY_Clusters_Lineage_only.csv")
  # Filter clusters with less than 10 clones
  # pool.w.clusters <- pool.w.clusters %>%
  #         group_by_at(clusters.column) %>%
  #         filter(n_distinct(cloneID)>=10) %>% ungroup()
  # pool.w.clusters[, clusters.column, drop=TRUE] <- as.numeric(pool.w.clusters[, clusters.column, drop=TRUE])
  print(table(pool.w.clusters[, clusters.column, drop=TRUE]))
  # Number of clones
  nrow(table(pool.w.clusters$cloneID))
  # Number of cells
  sum(table(pool.w.clusters$cloneID))
  lineage_dataset <- pool.w.clusters[,c("cellID", "cloneID", clusters.column)]
  lineage_dataset <- lineage_dataset %>% dplyr::rename(ident:=!!clusters.column)
  # num.shared.clones.matrix <- compute.num.shared.clones.matrix(lineage_dataset)
  if(use.default.names){
    folder.name <- ifelse(stage.bool, paste(dataset.name, stage, sep="_"), dataset.name)
    # trailing.clusterization <- paste0("_", clusters.column)
    trailing.clusterization <- paste0("_", clusters.column)
    trailing.lb.method <- ifelse(lb.calling.method=="May" || lb.calling.method=="cloneID", "", paste0("_", lb.calling.method))
    file.name <- paste0(folder.name, trailing.clusterization, trailing.lb.method)
    sub.folder <- "/lb_pool/"
    file.path.lineage.dataset <- paste0(parent.folder.path, folder.name
                                       , sub.folder, file.name, ".csv")
  }
  else{
    file.path.lineage.dataset <- provided.file.name
  }
  print(lineage_dataset)
  print(file.path.lineage.dataset)
  if(write.csv.bool){
    write.csv(lineage_dataset, file=file.path.lineage.dataset)
  }
}

prepare.input.files.for.lineage.coupling.analysis <- function(
        seur.objs.proj.ca, min.clonesizes=c(2)
        , lb.calling.methods=c("Ryan_toptier_9")
        , stages=c("", "e10", "e12", "e13", "e14")
        , stage.keep.bools=c(TRUE)
        , clusters.columns=c("refined_clust", "refined_class")
        , use.default.names.array=c(TRUE)
        , parent.folder.paths=c("./results/")
        , provided.file.names=c("")
        , dataset.names=c("clonal_dataset")
        , write.csv.bools=c(TRUE)
        , cols.subset=clusters.columns
        , val.subset.exclude=c("Remove")){
  # I didn't find a better way of iterating through the cartesian product
  # of all the possible parameter values than this way
  for(seur.obj.proj.ca in seur.objs.proj.ca){
    # Create data frame with every combination of possible parameters
    if(identical(cols.subset, clusters.columns) && identical(use.default.names.array, c(TRUE))){
      dat <- expand.grid(min.clonesize=min.clonesizes
                        , lb.calling.method=lb.calling.methods
                        , stage=stages, stage.keep.bool=stage.keep.bools
                        , clusters.column=clusters.columns
                        , use.default.names=use.default.names.array
                        , parent.folder.path=parent.folder.paths
                        , dataset.names=dataset.names
                        , write.csv.bool=write.csv.bools
                        , val.subset.exclude=val.subset.exclude
                        , stringsAsFactors=FALSE)
      provided.file.names <- ""
      col.subsets <- dat[[5]]
      vals.subset.exclude <- dat[[10]]
    }
    else if(!identical(cols.subset, clusters.columns) && identical(use.default.names.array, c(TRUE))){
      dat <- expand.grid(min.clonesize=min.clonesizes
                        , lb.calling.method=lb.calling.methods
                        , stage=stages, stage.keep.bool=stage.keep.bools
                        , clusters.column=clusters.columns
                        , use.default.names=use.default.names.array
                        , parent.folder.path=parent.folder.paths
                        , dataset.names=dataset.names
                        , write.csv.bool=write.csv.bools
                        , cols.subset=cols.subset
                        , val.subset.exclude=val.subset.exclude
                        , stringsAsFactors=FALSE)
      provided.file.names <- ""
      col.subsets <- dat[[10]]
      vals.subset.exclude <- dat[[11]]
    }
    else if(identical(cols.subset, clusters.columns) && identical(use.default.names.array, c(FALSE))){
      dat <- expand.grid(min.clonesize=min.clonesizes
                        , lb.calling.method=lb.calling.methods
                        , stage=stages, stage.keep.bool=stage.keep.bools
                        , clusters.column=clusters.columns
                        , use.default.names=use.default.names.array
                        , parent.folder.path=parent.folder.paths
                        , dataset.names=dataset.names
                        , write.csv.bool=write.csv.bools
                        , val.subset.exclude=val.subset.exclude
                        , stringsAsFactors=FALSE)
      col.subsets <- dat[[5]]
      vals.subset.exclude <- dat[[10]]
    }
    else if(!identical(cols.subset, clusters.columns) && identical(use.default.names.array, c(FALSE))){
      dat <- expand.grid(min.clonesize=min.clonesizes
                        , lb.calling.method=lb.calling.methods
                        , stage=stages, stage.keep.bool=stage.keep.bools
                        , clusters.column=clusters.columns
                        , use.default.names=use.default.names.array
                        , parent.folder.path=parent.folder.paths
                        , dataset.names=dataset.names
                        , write.csv.bool=write.csv.bools
                        , val.subset.exclude=val.subset.exclude
                        , stringsAsFactors=FALSE)
      col.subsets <- dat[[10]]
      vals.subset.exclude <- dat[[11]]
    }
    min.clonesize <- dat[[1]]
    lb.calling.methods <- dat[[2]]
    stages <- dat[[3]]
    stage.keep.bools <- dat[[4]]
    clusters.columns <- dat[[5]]
    use.default.names.array <- dat[[6]]
    parent.folder.paths <- dat[[7]]
    dataset.names <- dat[[8]]
    write.csv.bools <- dat[[9]]
    
    # Apply function by row
    function_map <- mapply(prepare.input.file.for.lineage.coupling.analysis
                            , as.array(rep(list(seur.obj.proj.ca), nrow(dat)))
                            , min.clonesize, lb.calling.methods, stages
                            , stage.keep.bools, clusters.columns, use.default.names.array
                            , parent.folder.paths, provided.file.names, dataset.names
                            , write.csv.bools, col.subsets, vals.subset.exclude
                    )
    
  }
}


panel.neuron.glia.composition.two.stages <- function(
              seur.obj.proj.ca, min.clonesize=2
              , clusters.column="refined_COUP_class"
              , neuron.glia.mapping=list(
                                  "Astro"="Glia"
                                  , "OPC"="Glia"
                                  , "Oligo"="Glia"
                                  , "Ependymal"="Glia"
                                  , "Neuron"="Neuron"
                                  , "Neuronal Precursor"="Neuron"
                                )
              , columns.values=list("stage"=c("e10", "e12", "e13", "e14"))
              , columns.na=list("cloneID"=FALSE, "neuron.glia.class"=FALSE)
              , parent.folder.results="./results/clonal_dataset"
              , folder.analysis.name="Quantif_Neurons_Glia_per_clone"
              , write.csv.bool=FALSE, write.pdf.bool=FALSE){
  folder.analysis.results.path <- paste(parent.folder.results, folder.analysis.name, sep="/")
  file.name <- paste(tolower(folder.analysis.name)
                     , "with_line_rr"
                     , paste(columns.values[["stage"]], collapse="_")
                     , ifelse("Neuronal Precursor" %in% names(neuron.glia.mapping), "with_NP_data" , "data")
                     , sep="_")
  file.name <- paste(file.name, "csv", sep=".")
  csv.file.path.clones.pctg.most.freq.class <- paste(folder.analysis.results.path, file.name, sep="/")
  
  columns.values[[clusters.column]] <- names(neuron.glia.mapping)
  
  seur.obj.proj.ca@meta.data[, names(columns.na)[2]] <- sapply(seur.obj.proj.ca@meta.data[, names(columns.values)[2]]
                                                    , function(x) ifelse(is.null(neuron.glia.mapping[[x]])
                                                                         , NA, neuron.glia.mapping[[x]])
                                                    )
  # columns.values must be a named list which stores the columns and values
  # a cell must have in order to take it into account.
  # For now I only consider the case when there are two columns,
  # the 2nd of which corresponds to the cluster/classes column to use.
  # columns.na is analogous but specific for NA values.
  # It should store the result of evaluating the expression is.na(names(columns.na)[idx]).
  # For now I only consider the case when there is just one column
  clones.pctg.most.freq.class <- seur.obj.proj.ca@meta.data %>%
          filter(is.na(.data[[names(columns.na)[1]]]) == columns.na[[1]]) %>%
          group_by(cloneID) %>%
          filter(n()>=min.clonesize
                 , Reduce("&", unique(.data[[names(columns.values)[1]]]) %in% columns.values[[1]])
                 , Reduce("&", unique(.data[[names(columns.values)[2]]]) %in% columns.values[[2]])
                 , Reduce("&", is.na(.data[[names(columns.na)[2]]]) == columns.na[[2]]))
  if(length(columns.values[[1]])>1){
    clones.pctg.most.freq.class <- clones.pctg.most.freq.class %>%
                    dplyr::count(.data[[names(columns.na)[2]]], .data[[names(columns.values)[1]]])
  }
  clones.pctg.most.freq.class <- clones.pctg.most.freq.class %>% ungroup() %>%
          pivot_wider(names_from=.data[[names(columns.na)[2]]], values_from=n, values_fill=list(n=0)) %>%
          mutate(category=ifelse(Glia>0 & Neuron==0, "Glia"
                                 , ifelse(Glia==0 & Neuron>0, "Neuron"
                                          , ifelse(Glia>0 & Neuron>0, "Neuron + Glia", "Other")
                                          )
                                 )
                    , num.cells.clone=Glia+Neuron) %>%
          arrange(desc(num.cells.clone))
  p <- (ggplot(clones.pctg.most.freq.class
               , aes(x=Neuron, y=Glia
                     , color=.data[[names(columns.values)[1]]]))
          + geom_count(position=position_dodge(0.8)) + scale_size_area()
          + labs(title="Quantifying neurons and glia per clone"
                 , x="Number of Neurons per clone"
                 , y="Number of Glia per clone"
                 , size="Number of clones"
                 , color="Stage")
          + geom_smooth(method=lm, se=FALSE, linetype="dashed")
          + stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")))
          + scale_x_continuous(breaks=seq(from=0, to=max(clones.pctg.most.freq.class$Neuron)+1, by=1))
          # + scale_y_continuous(breaks=seq(from=0, to=max(clones.pctg.most.freq.class$Glia)+1, by=5))
        )
  pdf.file.path <- paste(
                    substr(csv.file.path.clones.pctg.most.freq.class
                           , 1, nchar(csv.file.path.clones.pctg.most.freq.class)-nchar("_data.csv"))
                    , "pdf", sep=".")
  # View(clones.pctg.most.freq.class)
  print(pdf.file.path)
  if(write.pdf.bool){
    pdf(pdf.file.path, width = 11, height = 9)
  }
  print(p)
  if(write.pdf.bool){
    dev.off()
  }
  if(write.csv.bool){
    write.csv(clones.pctg.most.freq.class, file=csv.file.path.clones.pctg.most.freq.class)
  }
  return(p)
}


neuron.glia.mappings=list(
                        list(
                              "Astro"="Glia"
                              , "OPC"="Glia"
                              , "Oligo"="Glia"
                              , "Ependymal"="Glia"
                              , "Neuron"="Neuron"
                            ),
                        list(
                              "Astro"="Glia"
                              , "OPC"="Glia"
                              , "Oligo"="Glia"
                              , "Ependymal"="Glia"
                              , "Neuron"="Neuron"
                              , "Neuronal Precursor"="Neuron"
                            )
                        )
min.clonesizes <- c(2)
clusters.columns <- c("refined_COUP_class")
columns.values.list <- list(
                          list("stage"=c("e10"))
                          , list("stage"=c("e12"))
                          , list("stage"=c("e13"))
                          , list("stage"=c("e14"))
                          )
columns.na.list <- list(
                      list(
                          "cloneID"=FALSE
                          , "neuron.glia.class"=FALSE
                          )
                        )
parent.folders.results <- c("./results/clonal_dataset")
folder.analysis.names <- c("Quantif_Neurons_Glia_per_clone")
write.csv.bools <- c(FALSE)
write.pdf.bools <- c(FALSE)

panels.neuron.glia.composition.two.stages <- function(
            seur.objs.proj.ca
            , neuron.glia.mappings=neuron.glia.mappings
            , min.clonesizes=min.clonesizes
            , clusters.columns=clusters.columns
            , columns.values.list=columns.values.list
            , columns.na.list=columns.na.list
            , parent.folders.results=parent.folders.results
            , folder.analysis.names=folder.analysis.names
            , write.csv.bools, write.pdf.bools){
  # I didn't find a better way of iterating through the cartesian product
  # of all the possible parameter values than this way
  i <- 0
  for(seur.obj.proj.ca in seur.objs.proj.ca){
    for(neuron.glia.mapping in neuron.glia.mappings){
      for(columns.values in columns.values.list){
        for(columns.na in columns.na.list){
          # Create data frame with every combination of possible parameters
          dat <- expand.grid(min.clonesize=min.clonesizes
                             , clusters.columns=clusters.columns
                             , parent.folders.results=parent.folders.results
                             , folder.analysis.names=folder.analysis.names
                             , write.csv.bools, write.pdf.bools
                             , stringsAsFactors=FALSE)
          # print(paste("iteration", i))
          # i <- i + 1
          # print(seur.obj.proj.ca)
          # print(neuron.glia.mapping)
          # print(dat[[1]])
          # print(dat[[2]])
          # print(columns.values)
          # print(columns.na)
          # print(dat[[3]])
          # print(dat[[4]])
          # print(dat[[5]])
          # print(dat[[6]])
          panel.neuron.glia.composition.two.stages(seur.obj.proj.ca=seur.obj.proj.ca
                                        , neuron.glia.mapping=neuron.glia.mapping
                                        , min.clonesize=dat[[1]]
                                        , clusters.column=dat[[2]]
                                        , columns.values=columns.values
                                        , columns.na=columns.na
                                        , parent.folder.results=dat[[3]]
                                        , folder.analysis.name=dat[[4]]
                                        , write.pdf.bool=dat[[5]]
                                        , write.csv.bool=dat[[6]])
        }
      }
    }
  }
}

panel.neuron.glia.composition <- function(
              seur.obj.proj.ca, min.clonesize=2
              , clusters.column="refined_COUP_class"
              , neuron.glia.mapping=list(
                                  "Astro"="Glia"
                                  , "OPC"="Glia"
                                  , "Oligo"="Glia"
                                  , "Ependymal"="Glia"
                                  , "Neuron"="Neuron"
                                  , "Neuronal Precursor"="Neuron"
                                )
              , columns.values=list("stage"=c("e10", "e12", "e13", "e14"))
              , columns.na=list("cloneID"=FALSE, "neuron.glia.class"=FALSE)
              , parent.folder.results="./results/clonal_dataset"
              , folder.analysis.name="Neuron_Glia_Pie_Chart"
              , write.csv.bool=FALSE, write.pdf.bool=FALSE){
  folder.analysis.results.path <- paste(parent.folder.results, folder.analysis.name, sep="/")
  file.name <- paste(tolower(folder.analysis.name)
                     , paste(columns.values[["stage"]], collapse="_")
                     , ifelse("Neuronal Precursor" %in% names(neuron.glia.mapping), "with_NP_data" , "data")
                     , sep="_")
  file.name <- paste(file.name, "csv", sep=".")
  csv.file.path.clones.pctg.most.freq.class <- paste(folder.analysis.results.path, file.name, sep="/")
  
  columns.values[[clusters.column]] <- names(neuron.glia.mapping)
  
  seur.obj.proj.ca@meta.data[, names(columns.na)[2]] <- sapply(seur.obj.proj.ca@meta.data[, names(columns.values)[2]]
                                                    , function(x) ifelse(is.null(neuron.glia.mapping[[x]])
                                                                         , NA, neuron.glia.mapping[[x]])
                                                    )
  # columns.values must be a named list which stores the columns and values
  # a cell must have in order to take it into account.
  # For now I only consider the case when there are two columns,
  # the 2nd of which corresponds to the cluster/classes column to use.
  # columns.na is analogous but specific for NA values.
  # It should store the result of evaluating the expression is.na(names(columns.na)[idx]).
  # For now I only consider the case when there is just one column
  clones.pctg.most.freq.class <- seur.obj.proj.ca@meta.data %>%
          filter(is.na(.data[[names(columns.na)[1]]]) == columns.na[[1]]) %>%
          group_by(cloneID) %>%
          filter(n()>=min.clonesize
                 , Reduce("&", unique(.data[[names(columns.values)[1]]]) %in% columns.values[[1]])
                 , Reduce("&", unique(.data[[names(columns.values)[2]]]) %in% columns.values[[2]])
                 , Reduce("&", is.na(.data[[names(columns.na)[2]]]) == columns.na[[2]]))
  if(length(columns.values[[1]])>1){
    clones.pctg.most.freq.class <- clones.pctg.most.freq.class %>%
                    dplyr::count(.data[[names(columns.na)[2]]], .data[[names(columns.values)[1]]])
  }
  clones.pctg.most.freq.class <- clones.pctg.most.freq.class %>% ungroup() %>%
          pivot_wider(names_from=.data[[names(columns.na)[2]]], values_from=n, values_fill=list(n=0)) %>%
          mutate(category=ifelse(Glia>0 & Neuron==0, "Glia"
                                 , ifelse(Glia==0 & Neuron>0, "Neuron"
                                          , ifelse(Glia>0 & Neuron>0, "Neuron + Glia", "Other")
                                          )
                                 )
                    , num.cells.clone=Glia+Neuron) %>%
          arrange(desc(num.cells.clone))
  p <- (ggplot(clones.pctg.most.freq.class
               , aes(x=Neuron, y=Glia
                     , fill=.data[[names(columns.values)[1]]]))
          + geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(0.8),
                         stackratio=0.7, dotsize=0.2)
          + labs(title="Quantifying neurons and glia per clone"
                 , x="Number of Neurons per clone"
                 , y="Number of Glia per clone"
                 , size="Number of clones"
                 , fill="Stage")
          # + geom_smooth(method=lm, se=FALSE, linetype="dashed",
          #               color="darkred")
          + scale_x_continuous(breaks=seq(from=0, to=max(clones.pctg.most.freq.class$Neuron)+1, by=1))
          # + scale_y_continuous(breaks=seq(from=0, to=max(clones.pctg.most.freq.class$Glia)+1, by=5))
        )
  pdf.file.path <- paste(
                    substr(csv.file.path.clones.pctg.most.freq.class
                           , 1, nchar(csv.file.path.clones.pctg.most.freq.class)-nchar("_data.csv"))
                    , "pdf", sep=".")
  # View(clones.pctg.most.freq.class)
  print(pdf.file.path)
  if(write.pdf.bool){
    pdf(pdf.file.path, width = 11, height = 9)
  }
  print(p)
  if(write.pdf.bool){
    dev.off()
  }
  if(write.csv.bool){
    write.csv(clones.pctg.most.freq.class, file=csv.file.path.clones.pctg.most.freq.class)
  }
  return(p)
}


neuron.glia.mappings=list(
                        list(
                              "Astro"="Glia"
                              , "OPC"="Glia"
                              , "Oligo"="Glia"
                              , "Ependymal"="Glia"
                              , "Neuron"="Neuron"
                            ),
                        list(
                              "Astro"="Glia"
                              , "OPC"="Glia"
                              , "Oligo"="Glia"
                              , "Ependymal"="Glia"
                              , "Neuron"="Neuron"
                              , "Neuronal Precursor"="Neuron"
                            )
                        )
min.clonesizes <- c(2)
clusters.columns <- c("refined_COUP_class")
columns.values.list <- list(
                          list("stage"=c("e10"))
                          , list("stage"=c("e12"))
                          , list("stage"=c("e13"))
                          , list("stage"=c("e14"))
                          )
columns.na.list <- list(
                      list(
                          "cloneID"=FALSE
                          , "neuron.glia.class"=FALSE
                          )
                        )
parent.folders.results <- c("./results/clonal_dataset")
folder.analysis.names <- c("Quantif_Neurons_Glia_per_clone")
write.csv.bools <- c(FALSE)
write.pdf.bools <- c(FALSE)

panels.neuron.glia.composition <- function(
            seur.objs.proj.ca
            , neuron.glia.mappings=neuron.glia.mappings
            , min.clonesizes=min.clonesizes
            , clusters.columns=clusters.columns
            , columns.values.list=columns.values.list
            , columns.na.list=columns.na.list
            , parent.folders.results=parent.folders.results
            , folder.analysis.names=folder.analysis.names
            , write.csv.bools, write.pdf.bools){
  # I didn't find a better way of iterating through the cartesian product
  # of all the possible parameter values than this way
  i <- 0
  for(seur.obj.proj.ca in seur.objs.proj.ca){
    for(neuron.glia.mapping in neuron.glia.mappings){
      for(columns.values in columns.values.list){
        for(columns.na in columns.na.list){
          # Create data frame with every combination of possible parameters
          dat <- expand.grid(min.clonesize=min.clonesizes
                             , clusters.columns=clusters.columns
                             , parent.folders.results=parent.folders.results
                             , folder.analysis.names=folder.analysis.names
                             , write.csv.bools, write.pdf.bools
                             , stringsAsFactors=FALSE)
          # print(paste("iteration", i))
          # i <- i + 1
          # print(seur.obj.proj.ca)
          # print(neuron.glia.mapping)
          # print(dat[[1]])
          # print(dat[[2]])
          # print(columns.values)
          # print(columns.na)
          # print(dat[[3]])
          # print(dat[[4]])
          # print(dat[[5]])
          # print(dat[[6]])
          panel.neuron.glia.composition(seur.obj.proj.ca=seur.obj.proj.ca
                                        , neuron.glia.mapping=neuron.glia.mapping
                                        , min.clonesize=dat[[1]]
                                        , clusters.column=dat[[2]]
                                        , columns.values=columns.values
                                        , columns.na=columns.na
                                        , parent.folder.results=dat[[3]]
                                        , folder.analysis.name=dat[[4]]
                                        , write.pdf.bool=dat[[5]]
                                        , write.csv.bool=dat[[6]])
        }
      }
    }
  }
}