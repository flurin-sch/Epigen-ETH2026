# Assignment

* Choose a transcription factor (e.g. p300), and obtain peaks from ENCODE (ChIP-seq in a human context!)
* Find the ENCODE-rE2G enhancer-gene predictions for the corresponding cell-type (the thresholded element-gene links bed file)
* Subset to those peaks that have predicted rE2G target(s)
* Label peaks with the genes it is in 3D contact with (henceforth called ‘predicted target(s)’)
* Find the nearest TSS for each peak
* Split the peaks (or genes) into those for which the peak is upstream or downstream of the predicted target
* For each of the two sets, report:
  - 1) in what proportion of the cases the predicted target is the closest gene?
  - 2) plot a distribution of the distances between the peak and the (TSS of the) predicted target gene

Hint:
beware not to count, when calculating proportions, peaks that don’t have predicted interactions with any TSS!

Save your assignment in a R markdown named `assignment.Rmd`, render it, and push the html file to this folder in your github repository
