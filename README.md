# FlowSOM tuning operator

##### Description

`flowsomtuning` operator performs flowSOM clustering for different numbers of clusters.

##### Usage

Input projection|.
---|---
`row`   | represents the variables (e.g. channels, markers)
`col`   | represents the clusters (e.g. cells) 
`y-axis`| is the value of measurement signal of the channel/marker

Input parameters|.
---|---
`min_cluster_number`   | Minimum number of clusters to make
`max_cluster_number`   | Maximal number of clusters to make
`dunn`   | Compute Dunn Index
`davies_bouldin`   | Compute Davies Bouldin Index
`pseudo_f`| Compute the pseudo-F statistic
`silhouette`| Compute the silhouette value
`transform`| Transform data?
`seed`   | Random seed

Output relations|.
---|---
`cluster`| character, cluster label

##### Details

The operator is a wrapper for the `FlowSOM` function of the `FlowSOM` R/Bioconductor package.

#### References

https://bioconductor.org/packages/FlowSOM/

##### See Also

[flowsom_operator](https://github.com/tercen/flowsom_operator)

[flowsom_mst_shiny_operator](https://github.com/tercen/flowsom_mst_shiny_operator)
