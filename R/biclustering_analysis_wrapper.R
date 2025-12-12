

annotatePeakByOverlappingChIPSeekerCCRE <- function(peak, annotation, categories, featureColname="V10") {
  anno_data <- .get_annotations()
  txdb <- anno_data$txdb
  allFeatures <- anno_data$allFeatures

  chipSeekerAnno <- ChIPseeker::annotatePeak(peak, TxDb=txdb, tssRegion=c(-1000, 1000), verbose=FALSE)
  chipSeekerAnnoFreq <-  c()
  chipSeekerAnnoFreq[allFeatures] <- 0
  chipSeekerAnnoFreq[as.character(chipSeekerAnno@annoStat$Feature)] <- chipSeekerAnno@annoStat$Frequency
  condensedFeatures <- c("5' UTR", "3' UTR", "Exon", "Intron")
  chipSeekerAnnoFreqCondensed <- c()
  chipSeekerAnnoFreqCondensed["5' UTR"] <- chipSeekerAnnoFreq["5' UTR"]
  chipSeekerAnnoFreqCondensed["3' UTR"] <- chipSeekerAnnoFreq["3' UTR"]
  chipSeekerAnnoFreqCondensed["Exon"] <- chipSeekerAnnoFreq["1st Exon"] + chipSeekerAnnoFreq["Other Exon"]
  chipSeekerAnnoFreqCondensed["Intron"] <- chipSeekerAnnoFreq["1st Intron"] + chipSeekerAnnoFreq["Other Intron"]
  chipSeekerAnnoFreqCondensed <- chipSeekerAnnoFreqCondensed * chipSeekerAnno@peakNum / 100

  chipSeekerAnno@anno$annotation <- gsub("\\s*\\([^\\)]+\\)","",chipSeekerAnno@anno$annotation)
  finalAnno <- rep("other", length(peak))
  finalAnno[chipSeekerAnno@anno$annotation %in% condensedFeatures] <- chipSeekerAnno@anno$annotation[chipSeekerAnno@anno$annotation %in% condensedFeatures]
  chipSeekerAnno@anno$index <- 1: length(peak)
  disposeFeatures <- c("Promoter", "Downstream (<=300)", "Distal Intergenic")
  unmappedPeakGR <- chipSeekerAnno@anno[chipSeekerAnno@anno$annotation %in% disposeFeatures]
  if (length(unmappedPeakGR) > 0) {
    annoResult <- annotatePeakByOverlappingClosestFeatureHelper(unmappedPeakGR, annotation, categories, featureColname)
    peakCategories <- annoResult[[featureColname]]
    peakCategoriesTable <- table(peakCategories)
    finalAnno[annoResult$queryIndex] <- annoResult[[featureColname]]
  }

  if (length(unmappedPeakGR)==length(peak)) {
    otherCt <- length(peak) - sum(peakCategoriesTable)
    peakFeature <- c(names(peakCategoriesTable), "other")
    peakFreq <- c(unname(peakCategoriesTable), otherCt) / length(peak) * 100
  } else if (length(unmappedPeakGR)==0) {
    otherCt <- length(peak) - sum(chipSeekerAnnoFreqCondensed)
    peakFeature <- c(names(chipSeekerAnnoFreqCondensed), "other")
    peakFreq <- c(unname(chipSeekerAnnoFreqCondensed), otherCt) / length(peak) * 100
  } else {
    otherCt <- length(peak) - sum(peakCategoriesTable) - sum(chipSeekerAnnoFreqCondensed)
    peakFeature <- c(names(peakCategoriesTable), names(chipSeekerAnnoFreqCondensed), "other")
    peakFreq <- c(unname(peakCategoriesTable), unname(chipSeekerAnnoFreqCondensed), otherCt) / length(peak) * 100
  }
  peak$annotation <- finalAnno
  res <- data.frame(Feature=peakFeature, Frequency = peakFreq)
  x <- new("csCCREAnno", annoStat = res, peakNum=length(peak), anno=peak)
  return(x)
}




# Annotate and Plot CCRE and ChromHMM Composition for Row Cluster
# Post: Analyzes genomic annotation composition of clustered regions by overlapping with CCRE and ChromHMM annotations, then generates comparative bar plots showing regulatory element distribution across clusters.
# Parameter: row_cluster_file_path: Path to row cluster assignment .tsv file with feature positions and cluster labels
#            out_dir: Directory to save annotation plot outputs
# Output: Generates and saves two bar plot PDF files for CCRE and ChromHMM composition analysis

annotation_ccre_hmm <- function(row_cluster_file_path, out_dir, ref_genome = "hg38", txdb = NULL) {
  # load library
  suppressPackageStartupMessages({
    library(ggplot2)
    library(cowplot)
    library(GenomeInfoDb)
    library(glue)
    library(tidyverse)
    library(ChIPseeker)
    library(tools)
    library(stringr)
    library(dplyr)
    library(ComplexHeatmap)
    library(grid)
    library("gridExtra")
  })

  # split pos
  biclustering_result <- read.table(row_cluster_file_path, header = TRUE, sep = "\t")
  # biclustering_result <- biclustering_result[!biclustering_result$label %in% c("Epitope_Specific", "Background"), ]
  pos_df <- do.call(rbind, (strsplit(biclustering_result$feature, "_")))
  colnames(pos_df) <- c("seqnames", "start", "end")
  biclustering_result <- cbind(pos_df, biclustering_result)
  biclustering_gr <- makeGRangesFromDataFrame(
    biclustering_result,
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
  )
  current_style <- seqlevelsStyle(biclustering_gr)[1]

  # To-do: need to check whether use as an input variable
  cluster_all_levels <- c(
    as.character(1:100),
    "A",
    "B",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "J",
    "K",
    "L",
    "M",
    "N",
    "O",
    "CRF_specific",
    "Background"
  )

  if (is.null(txdb)) {
    if (ref_genome == "hg38") {
      txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    } else if (ref_genome == "mm10") {
      txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
    } else {
      stop("Please provide TxDb object or use hg19/mm10 as ref_genome")
    }
  }
  target_style <- seqlevelsStyle(txdb)[1]

  if (current_style != target_style) {
    message(sprintf("Converting chromosome style from %s to %s", current_style, target_style))
    tryCatch({
      seqlevelsStyle(biclustering_gr) <- target_style
    }, error = function(e) {
      warning("Chromosome style conversion failed. Trying manual conversion...")
      if (target_style == "UCSC" && !any(grepl("^chr", seqlevels(biclustering_gr)))) {
        seqlevels(biclustering_gr) <- paste0("chr", seqlevels(biclustering_gr))
      } else if (target_style == "NCBI" && any(grepl("^chr", seqlevels(biclustering_gr)))) {
        seqlevels(biclustering_gr) <- gsub("^chr", "", seqlevels(biclustering_gr))
      }
    })
  }




  unique_clusters <- sort(unique(biclustering_result$label), method = "radix")
  chipseeker_ccre_annotations <- list()
  chipseeker_ccre_celltype_agnostic_annotations <- list()
  repeatmasker_annotations <- list()
  chromhmm_short_annotations <- list()
  chromhmm_full_annotations <- list()
  chromhmm_group_annotations <- list()
  for (cluster_id in unique_clusters) {
    biclustering_result_i <- biclustering_result[biclustering_result$label == cluster_id, ]
    biclustering_grange_i <- makeGRangesFromDataFrame(biclustering_result_i, seqnames.field = "seqnames", start.field = "start", end.field = "end", keep.extra.columns = TRUE)
    #####
    data("annotation", package = "epigenomeR", envir = environment())
    data("annotation_celltype_agnostic", package = "epigenomeR", envir = environment())
    data("annotationChromHMM", package = "epigenomeR", envir = environment())
    data("categories", package = "epigenomeR", envir = environment())
    data("categoriesChromHMM", package = "epigenomeR", envir = environment())
    data("repeatMaskerFeatures", package = "epigenomeR", envir = environment())
    data("ChIPSeekerCCRECategoriesOrder", package = "epigenomeR", envir = environment())

    annotationRepeatMasker <- epigenomeR::getRepeatMaskerAnnotation()
    annotationRepeatMasker <- epigenomeR::getRepeatMaskerAnnotation()
    #####
    chipseeker_ccre_annotation <- annotatePeakByOverlappingChIPSeekerCCRE(peak = biclustering_grange_i, annotation = annotation, categories = categories, txdb = txdb)
    chipseeker_ccre_celltype_agnostic_annotation <- annotatePeakByOverlappingChIPSeekerCCRE(peak = biclustering_grange_i, annotation = annotation_celltype_agnostic, categories = categories, featureColname="V6", txdb = txdb)
    repeatmasker_annotation <- annotatePeakByOverlappingRepeatMasker(biclustering_grange_i,annotationRepeatMasker, repeatMaskerFeatures)
    #####
    chromhmm_short_annotation <- annotatepeakByOverlappingChromHMM(peak = biclustering_grange_i, annotation = annotationChromHMM, categoriesChromHMM = categoriesChromHMM, featureColname = "V4")
    chromhmm_full_annotation <- annotatepeakByOverlappingChromHMM(peak = biclustering_grange_i, annotation = annotationChromHMM, categoriesChromHMM = unique(annotationChromHMM$full_anno), featureColname = "full_anno")
    chromhmm_group_annotation <- annotatepeakByOverlappingChromHMM(peak = biclustering_grange_i, annotation = annotationChromHMM, categoriesChromHMM = unique(annotationChromHMM$group), featureColname = "group")
    chipseeker_ccre_annotations[[glue("{cluster_id}")]] <- chipseeker_ccre_annotation
    repeatmasker_annotations[[glue("{cluster_id}")]] <- repeatmasker_annotation
    chromhmm_short_annotations[[glue("{cluster_id}")]] <- chromhmm_short_annotation
    chromhmm_full_annotations[[glue("{cluster_id}")]] <- chromhmm_full_annotation
    chromhmm_group_annotations[[glue("{cluster_id}")]] <- chromhmm_group_annotation
    chipseeker_ccre_celltype_agnostic_annotations[[glue("{cluster_id}")]] <- chipseeker_ccre_celltype_agnostic_annotation
  }

  # plotting chipseeker + ccre
  anno_ccre <- lapply(chipseeker_ccre_celltype_agnostic_annotations, getAnnoStatCCRE)
  anno_ccre.df <- list_to_dataframe(anno_ccre)
  anno_ccre.df$Feature <- factor(anno_ccre.df$Feature, levels = c(ChIPSeekerCCRECategoriesOrder))
  categoryColumn <- ".id"

  cluster_levels <- cluster_all_levels[cluster_all_levels %in% unique(biclustering_result$label)]
  anno_ccre.df$.id <- factor(anno_ccre.df$.id, levels = rev(cluster_levels))
  print(anno_ccre.df$.id)
  p_ccre_agnostic <- plotAnnoBar.data.frame.one.target(anno_ccre.df,
                                                       categoryColumn = categoryColumn,
                                                       colorOption = 1,
                                                       features=ChIPSeekerCCRECategoriesOrder) + theme(
                                                         legend.position = "right",
                                                         plot.title = element_text(size = 18, ),
                                                         axis.text.x = element_text(size = 12, ),
                                                         axis.text.y = element_text(size = 12, ),
                                                         legend.text = element_text(size = 9),
                                                         legend.key.size = unit(0.5, 'cm'),
                                                         axis.title.x = element_text(size = 20)
                                                       ) +
    ggtitle("Cis-regulatory Elements (cell type agnostic)") +
    guides(fill = guide_legend(title = NULL, ncol = 1, reverse = TRUE))


  # plotting chromhmm
  anno_chromhmm_short <- lapply(chromhmm_short_annotations, getAnnoStatCCRE)
  anno_chromhmm_short.df <- list_to_dataframe(anno_chromhmm_short)
  anno_chromhmm_short.df$Feature <- factor(anno_chromhmm_short.df$Feature, levels = c(unique(annotationChromHMM$V4), "other"))
  categoryColumn <- ".id"
  anno_chromhmm_short.df$.id <- factor(anno_chromhmm_short.df$.id, levels = rev(cluster_levels))
  print(anno_chromhmm_short.df$.id)
  p_chromhmm_short <- plotAnnoBar.data.frame.one.target(anno_chromhmm_short.df,
                                                        categoryColumn = categoryColumn,
                                                        colorOption = 6,
                                                        features=c(unique(annotationChromHMM$V4), "other")) + theme(
                                                          legend.position = "right",
                                                          plot.title = element_text(size = 18, ),
                                                          axis.text.x = element_text(size = 12, ),
                                                          axis.text.y = element_text(size = 12, ),
                                                          legend.text = element_text(size = 9),
                                                          legend.key.size = unit(0.5, 'cm'),
                                                          axis.title.x = element_text(size = 20)
                                                        ) +
    ggtitle("ChromHMM Elements") +
    guides(fill = guide_legend(title = NULL, ncol = 1, reverse = TRUE))

  # save
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  ggsave(file.path(out_dir, "p_ccre_agnostic.pdf"), plot = p_ccre_agnostic, width = 6, height = 5)
  ggsave(file.path(out_dir, "p_chromhmm_short.pdf"), plot = p_chromhmm_short, width = 6, height = 5)

}
