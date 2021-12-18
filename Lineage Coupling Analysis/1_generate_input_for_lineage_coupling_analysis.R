source("./clonal_analysis_common_functions.R")

# Parameters

dataset.name <- "clonal_dataset"
# Seurat object to analyze
seur.objs.proj.ca <- c(readRDS(file = "./STICR.seuratobject.RDS"))
# Clones of size less than min.clonesizes will be filtered out
min.clonesizes <- c(2)
# Lineage barcode calling method, which will determine
# the column from where cloneIDs will be retrieved
lb.calling.methods <- c("Ryan_toptier_9")
# If specified, in conjunction with the next parameter,
# defines whether and what filter to apply on cells
# depending on the stage at which they were sampled
stages <- c("", "e10", "e12", "e13", "e14")
# If a stage was specified, indicates whether to
# keep only cells sampled at that stage, or exclude only these
stage.keep.bools <- c(TRUE)
# Clusterization of the cells
clusters.columns <- c("refined_COUP_clust", "refined_COUP_class")
# Output folder path
parent.folder.paths <- c("./results/")
# Name of dataset. Used for output directory/file names
dataset.names <- c(dataset.name)
# Whether to actually output the results into a file
write.csv.bools <- c(FALSE)
# If any of the cells has any column with this value, it will be
val.subset.exclude <- c("Remove")


prepare.input.files.for.lineage.coupling.analysis(
        seur.objs.proj.ca=seur.objs.proj.ca
        , min.clonesizes=min.clonesizes
        , lb.calling.methods=lb.calling.methods
        , stages=stages, stage.keep.bools=stage.keep.bools
        , clusters.columns=clusters.columns
        , parent.folder.paths=parent.folder.paths
        , dataset.names=dataset.names
        , write.csv.bools=write.csv.bools
        , val.subset.exclude=val.subset.exclude)