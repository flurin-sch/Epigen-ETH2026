# Assignment

* Download ATAC-seq peak counts in the hippocampus upon stress (subset of the original data, already in SummarizedExperiment format) :
  * https://ethz-ins.org/content/mouse_mm38_hippocampus.peakCounts.SE.rds
* Using this object, run 2 differential analyses, respectively:
  * comparing stressed (denoted ‘FSS’ – forced swim stress) and control animals
  * comparing male and female animals
* (No need to worry about normalization here, you can just use the default TMM)
* For each analysis, report the top 10 most significant regions, and plot a heatmap of their normalized counts across the samples.

Save your assignment in a R markdown named `assignment.Rmd`, render it, and push the html file to this folder in your github repository
