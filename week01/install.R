install.packages(c("gsl","ggplot2","remotes","rmarkdown"))
# If you haven't installed Bioconductor yet:
install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c("GenomicRanges", "rtracklayer", "EnrichedHeatmap", "AnnotationHub", 
                       "ensembldb", "edgeR", "esATAC", "sechm", "motifmatchr","rGREAT",
                       "bsseq","DMRcate","data.table","InteractionSet","limma",
                       "universalmotif", "MotifDb", "TFBSTools", "Biostrings", "PWMEnrich",
                       "Rsubread","Rfastp"), force=TRUE)
BiocManager::install(c("Bioconductor/BiocFileCache", "chromVAR"), force=TRUE)
BiocManager::install("ETHZ-INS/epiwraps", ref="devel", force=TRUE)
