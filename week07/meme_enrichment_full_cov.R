#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(memes)
  library(MotifDb)
  library(universalmotif)
  library(Biostrings)
  library(Rsamtools)
})

bw_path <- "data/full_cov.bw"
out_dir <- "data/task1_meme_fullcov"
fasta_path <- "../week06/data/mm10.fa"
target_chr <- "chr19"
quantile_cutoff <- 0.95
window_width <- 150L
min_region_width <- 20L

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(bw_path)) {
  stop("Missing input bigWig: ", bw_path)
}

if (dir.exists("/opt/meme/bin")) {
  options(meme_bin = "/opt/meme/bin")
}

ame_path_from_option <- file.path(getOption("meme_bin", ""), "ame")
if (Sys.which("ame") == "" && !file.exists(ame_path_from_option)) {
  stop("MEME AME binary not found. Ensure MEME is installed and 'ame' is in PATH or set options(meme_bin).")
}

message("Importing bigWig signal...")
gr <- import(bw_path)
scores <- mcols(gr)$score
valid <- is.finite(scores) & scores > 0
if (!any(valid)) {
  stop("No positive finite signal values found in ", bw_path)
}

threshold <- as.numeric(quantile(scores[valid], probs = quantile_cutoff, names = FALSE))
message("Using coverage threshold at q", quantile_cutoff, ": ", signif(threshold, 4))

high_signal <- gr[is.finite(scores) & scores >= threshold]
suppressWarnings(seqlevelsStyle(high_signal) <- "UCSC")
chr_regions <- high_signal[as.character(seqnames(high_signal)) == target_chr]
if (length(chr_regions) == 0) {
  stop("No high-signal bins found on ", target_chr, " after thresholding.")
}

peak_like_regions <- reduce(chr_regions, min.gapwidth = 1L)
peak_like_regions <- peak_like_regions[width(peak_like_regions) >= min_region_width]
if (length(peak_like_regions) == 0) {
  stop("No regions remained after merging and minimum width filtering.")
}

peak_windows <- resize(peak_like_regions, fix = "center", width = window_width)
peak_windows <- trim(peak_windows)
peak_windows$name <- paste0("peak_", seq_along(peak_windows))

export(peak_windows, file.path(out_dir, "peak_windows.bed"))

if (!file.exists(fasta_path)) {
  stop("Reference FASTA not found at ", fasta_path, ". Run week06 setup or provide mm10 FASTA there.")
}
if (!file.exists(paste0(fasta_path, ".fai"))) {
  indexFa(fasta_path)
}
genome_fa <- FaFile(fasta_path)

message("Extracting peak-centered sequences...")
peak_seqs <- get_sequence(peak_windows, genome_fa)
names(peak_seqs) <- peak_windows$name
writeXStringSet(peak_seqs, file.path(out_dir, "peak_windows.fa"))

message("Loading motifs and running AME...")
mouse_motifs <- query(MotifDb, c("Mmusculus", "HOCOMOCOv13"))
if (length(mouse_motifs) == 0) {
  mouse_motifs <- query(MotifDb, c("Mmusculus", "HOCOMOCOv10"))
}
if (length(mouse_motifs) == 0) {
  mouse_motifs <- query(MotifDb, c("Mus", "HOCOMOCO"))
}
if (length(mouse_motifs) == 0) {
  stop("No mouse HOCOMOCO motifs found in MotifDb.")
}

ame <- runAme(
  input = peak_seqs,
  database = convert_motifs(mouse_motifs),
  silent = TRUE
)

ame_df <- as.data.frame(ame)
ame_df$tf_symbol <- gsub("\\..*", "", ame_df$motif_id)
ame_df <- ame_df[order(ame_df$adj.pvalue, ame_df$pvalue), ]

write.table(
  ame_df,
  file = file.path(out_dir, "ame_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
saveRDS(ame_df, file.path(out_dir, "ame_results.rds"))

tf_shortlist <- ame_df[!duplicated(ame_df$tf_symbol), ]
tf_shortlist <- head(tf_shortlist, 25)
write.table(
  tf_shortlist,
  file = file.path(out_dir, "enriched_tf_shortlist.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Done.")
message("Regions analyzed: ", length(peak_windows))
message("Top enriched TF motifs:")
print(tf_shortlist[seq_len(min(15, nrow(tf_shortlist))), c("motif_id", "tf_symbol", "adj.pvalue", "pvalue")], row.names = FALSE)
