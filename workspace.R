library(tercen)
library(tidyverse)
library(flowCore)
library(FlowSOM)
# library(clValid)
# library(clusterSim)

options("tercen.workflowId"= "7eee20aa9d6cc4eb9d7f2cc2430313b6")
options("tercen.stepId"= "a3f464fd-cd95-41fa-97e2-6e9b058a6269")
getOption("tercen.workflowId")
getOption("tercen.stepId")


as_flowset <- function(a_matrix) {
  minRange <- matrixStats::colMins(a_matrix)
  maxRange <- matrixStats::colMaxs(a_matrix)
  range <- minRange - maxRange
  
  df_params <-
    data.frame(
      name = colnames(a_matrix),
      desc = colnames(a_matrix),
      range = range,
      minRange = minRange,
      maxRange = maxRange
    )
  a_params <- Biobase::AnnotatedDataFrame()
  Biobase::pData(a_params) <- df_params
  Biobase::varMetadata(a_params) <-
    data.frame(
      labelDescription = c(
        "Name of Parameter",
        "Description of Parameter",
        "Range of Parameter",
        "Minimum Parameter Value after Transformation",
        "Maximum Parameter Value after Transformation"
      )
    )
  
  a_flowset <- flowCore::flowFrame(a_matrix, a_params)
  
}

ctx = tercenCtx()
data = ctx$as.matrix() 

# put row names on matrix
rownames(data) <- unite(ctx$rselect(), "rows_united", 1:ncol(ctx$rselect()))[,1, drop=TRUE]

# transpose for flowCore
data = t(data)

# convert to flowset
f_set <-  as_flowset(data)

# initialize 
# start_clus =  as.integer(ctx$op.value('min_cluster_number'))
# end_clus =  as.integer(ctx$op.value('max_cluster_number'))
start_clus = 3
end_clus = 5

# run_dunn <- as.logical(ctx$op.value('dunn'))
# run_davies_bouldin <- as.logical(ctx$op.value('davies_bouldin'))
# run_pseudo_f <- as.logical(ctx$op.value('pseudo_f'))
# run_silhouette<- as.logical(ctx$op.value('silhouette'))

transform_flag <- FALSE
#transform_flag <- as.logical(ctx$op.value('transform'))
transform_cols <- NULL
if (transform_flag) transform_cols = c(1:ncol(data))

# repeatedly call flowsom and calculate metrics
tuning <- lapply(start_clus:end_clus, function(x) {
  
  som <- FlowSOM(f_set, nClus = x, transform  = transform_flag, toTransform = transform_cols,  colsToUse = c(1:ncol(data)))
  # Extracting vector of meta-clusters for each cell
  cluster_vector_num <- as.integer(som[[2]][som[[1]]$map$mapping[,1]])
  cluster_vector_label <- str_pad(as.character(cluster_vector_num), 2, pad = "0")

  return(cluster_vector_label)
})

names(tuning) <- paste0(start_clus:end_clus, "_clusters")

tuning <- as.data.frame(do.call(cbind, tuning))
tuning$.ci = seq_len(nrow(tuning)) - 1

tuning <- ctx$addNamespace(tuning)
ctx$save(tuning)
