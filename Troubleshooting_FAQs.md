# Troubleshooting - Frequently Asked Questions (FAQs)

This is a list of some technical errors people often encounter, especially when trying to set-up the environment.

## Error: "the package XXX is not available"
Unless specified otherwise, always install packages using `BiocManager::install("packagename")` instead of using `install.packages` .

## Error: Rtools is required to build R packages but is not currently installed
On windows, you'll have to install Rtools to be able to build packages. On https://cran.r-project.org/bin/windows/Rtools/ you can download the Rtools version matching your R version (in R, you can see your R version with `R.version$version.string` )

## Error: some file can’t be read / is not found when knitting the Rmarkdown, but interactively it works
In Rmarkdown document, all filepaths are relative to the document itself. This means that if your markdown is, say, in /whatever/folder/assignment.Rmd, then trying to read the file raw/file.fastq.gz will try to find this file in /whatever/folder/raw/file.fastq.gz .
To avoid discrepancy between the current working directory in your Rstudio environment and that of the markdown, when you have your markdown open in Rstudio it's a good idea to do, from the top menu: Session -> Set Working Directory -> To Source File Location
(That will set the working directory of your R session to the directory containing the markdown, so that paths that work interactively will also work when knitting)

## Error: ... Permission denied (typically windows users): 
1. Most likely this is happening because you said yes ("all") to the question of whether or not to update other packages than the one you were trying to install. An easy way out (though not a solution in the long term) is simply not to do so.
Alternatively:
1. Run Rstudio / IDE as an Administrator
2. If this still doesn't work: prevent Rstudio to load .RData from previous session. https://stackoverflow.com/questions/45173367/how-can-i-stop-a-package-from-loading-on-startup-in-rstudio
=> then restart session and install packages needed

If it's about rlang and it still doesn’t work => install from source in the fresh session: install.packages("rlang", type = "source")

## Error: UnableToOpenBlob when trying to plotEnrichedHeatmaps

If you are on windows and this is accompanied by a warning, just before, saying "Failed to load cairo DLL", then you should try to build the Cairo package from source, i.e. `install.packages("Cairo", force=TRUE, type="source")`.
Otherwise this could be simply a caching issue, and you should try, in the cogweel menu above the markdown, to "Clear all outputs" and restart in a fresh R session.

## With the same code, heatmaps look uglier than in class
You probably are also having a message that you should install the `magick` R package. Just install that package, 
(When plotting heatmaps with more data points than there are pixels to display them, the matrix needs to be rasterized and shrinked. The best way to do this is to use the imagemagick image manipulation library, which will be used if available. If not, some lousy alternative gets used instead.)

