library(tercen)
library(tidyverse)
library(flowCore)
library(FlowSOM)
library(clValid)
library(clusterSim)



as_flowset <- function(a_matrix){ 
  
  minRange<- matrixStats::colMins(a_matrix)
  maxRange<- matrixStats::colMaxs(a_matrix)
  range<- minRange-maxRange
  
  df_params <- data.frame(name=colnames(a_matrix), desc=colnames(a_matrix), range=range, minRange=minRange, maxRange=maxRange)
  a_params <- Biobase::AnnotatedDataFrame()
  Biobase::pData(a_params) <- df_params
  Biobase::varMetadata(a_params) <- data.frame(labelDescription=c("Name of Parameter", "Description of Parameter","Range of Parameter","Minimum Parameter Value after Transformation","Maximum Parameter Value after Transformation"))
  
  a_flowset <- flowCore::flowFrame(a_matrix, a_params)
  
}


# http://127.0.0.1:5402/#ds/e466d41ce38cb8344c394688130067d1/11-3
# options("tercen.workflowId"= "e466d41ce38cb8344c394688130067d1")
# options("tercen.stepId"= "11-3")
ctx = tercenCtx()
data = ctx$as.matrix() 

# put row names on matrix
rownames(data) <- unite(ctx$rselect(), "rows_united", 1:ncol(ctx$rselect()))[,1, drop=TRUE]

# transpose for flowCore
data = t(data)

# convert to flowset
f_set <-  as_flowset(data)

# initialize 
start_clus =  as.integer(ctx$op.value('min_cluster_number'))
end_clus =  as.integer(ctx$op.value('max_cluster_number'))

run_dunn <- as.logical(ctx$op.value('dunn'))
run_davies_bouldin <- as.logical(ctx$op.value('davies-boudin'))
run_pseudo_f <- as.logical(ctx$op.value('pseudo_f'))
run_silhoutte <- as.logical(ctx$op.value('silhoutte'))

transform_flag <- as.logical(ctx$op.value('transform'))
transform_cols <- NULL
if (transform_flag) transform_cols = c(1:ncol(data))

# repeatedly call flowsom and calculate metrics
tuning <- lapply(start_clus:end_clus, function(x) {
  
  som <- FlowSOM(f_set, nClus = x, transform  = transform_flag, toTransform = transform_cols,  colsToUse = c(1:ncol(data)))
  # Extracting vector of meta-clusters for each cell
  cluster_vector_num <- as.integer(som[[2]][som[[1]]$map$mapping[,1]])
  cluster_vector_label <- str_pad(as.character(cluster_vector_num), 2, pad = "0")
  
  # extract the transformed data from the som object
  som_data <- as.data.frame(som$FlowSOM$data)
  
  # Calculate internal clustering metrics
  results_df <- data.frame(.ci = seq(from=0, to=length(cluster_vector_num)-1), 
                           cluster_label = as.character(cluster_vector_label),
                           cluster_setting = as.numeric(x),
                           stringsAsFactors = FALSE)
  
  if (run_dunn)           {
    results_df <- data.frame(
      results_df,
      dunn = dunn(dist(som_data), cluster_vector_num),
      stringsAsFactors = FALSE
    )
  }
  if (run_davies_bouldin) {
    results_df <- data.frame(
      results_df,
      davies_bouldin = index.DB(som_data, cluster_vector_num)$DB,
      stringsAsFactors = FALSE
    )
  }
  if (run_pseudo_f)       {
    results_df <- data.frame(
      results_df,
      pseudo_f = index.G1(som_data, cluster_vector_num),
      stringsAsFactors = FALSE
    )
  }
  if (run_silhoutte)      {
    results_df <- data.frame(
      results_df,
      silhouette = index.S(dist(som_data), cluster_vector_num),
      stringsAsFactors = FALSE
    )
  }
  return(results_df)
})

tuning <- do.call(rbind, tuning)
tuning <- tidyr::gather(tuning, key ="metrics_name", value = "metrics_value", -.ci, -cluster_label, -cluster_setting)

tuning <- ctx$addNamespace(tuning)
ctx$save(tuning)
