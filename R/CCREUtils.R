
# CCREUtils.R
# ========== refine from CCREUtils.R ==========

setClass("csCCREAnno",
         representation=representation(
           annoStat="data.frame",
           peakNum="numeric",
           anno="GRanges"
         ))

getAnnoStatCCRE <- function(x) {
  if (!is(x, "csCCREAnno"))
    stop("not supported...")
  return(x@annoStat)
}

getAnnotateStatPeakByOverlappingClosestFeatureHelper <- function(peak, annotation, categories, featureColname = "V10") {
  allPeakCtsList <- list()
  for (category in categories) {
    allPeakCtsList[[category]] <- 0
  }
  overlapRes <- findOverlaps(peak, annotation)
  peakHitsUnique <- unique(queryHits(overlapRes))
  peakHits <- queryHits(overlapRes)
  peakHitsGR <- peak[peakHits]
  annoHits <- subjectHits(overlapRes)
  annoHitsGR <- annotation[annoHits]
  peakCenter <- ceiling((end(peakHitsGR)-start(peakHitsGR))/2 + start(peakHitsGR))
  annoStarts <- start(annoHitsGR)
  annoEnds <- end(annoHitsGR)
  startsDistance <- abs(peakCenter - annoStarts)
  endsDistance <- abs(peakCenter - annoEnds)
  minDistance <- pmin(startsDistance, endsDistance)

  annoHitsDT <- as.data.table(annoHitsGR)
  annoHitsDT$peakHits <- peakHits
  annoHitsDT$minDistance <- minDistance
  annoHitsDTGroup <- annoHitsDT[annoHitsDT[ , .I[which.min(minDistance)], by = peakHits]$V1]

  peakCategories <- annoHitsDTGroup[[featureColname]]
  peakCategoriesTable <- table(peakCategories)
  return(peakCategoriesTable)
}

annotatePeakByOverlappingClosestFeature <- function(peak, annotation, categories, featureColname = "V10") {
  peakCategoriesTable <- getAnnotateStatPeakByOverlappingClosestFeatureHelper(peak, annotation, categories, featureColname)
  otherLength <- length(peak)-sum(peakCategoriesTable)
  peakFreq <- c(unname(peakCategoriesTable), otherLength) / length(peak) * 100
  res <- data.frame(Feature=c(names(peakCategoriesTable), "other"), Frequency = peakFreq)
  #x <- new("csCCREAnno", annoStat = res, peakNum=length(peak))
  x <- new("csCCREAnno", annoStat = res, peakNum = length(peak), anno = peak)
  return(x)
}

annotatePeakByOverlappingClosestFeatureHelper <- function(peak, annotation, categories, featureColname = "V10") {
  allPeakCtsList <- list()
  for (category in categories) {
    allPeakCtsList[[category]] <- 0
  }
  overlapRes <- findOverlaps(peak, annotation)
  peakHitsUnique <- unique(queryHits(overlapRes))
  peakHits <- queryHits(overlapRes)
  peakHitsGR <- peak[peakHits]
  annoHits <- subjectHits(overlapRes)
  annoHitsGR <- annotation[annoHits]
  peakCenter <- ceiling((end(peakHitsGR)-start(peakHitsGR))/2 + start(peakHitsGR))
  annoStarts <- start(annoHitsGR)
  annoEnds <- end(annoHitsGR)
  startsDistance <- abs(peakCenter - annoStarts)
  endsDistance <- abs(peakCenter - annoEnds)
  minDistance <- pmin(startsDistance, endsDistance)

  annoHitsDT <- as.data.table(annoHitsGR)
  annoHitsDT$peakHits <- peakHits
  annoHitsDT$minDistance <- minDistance
  annoHitsDT$queryStart <- start(peakHitsGR)
  annoHitsDT$queryEnd <- end(peakHitsGR)
  annoHitsDT$queryIndex <- peakHitsGR$index

  annoHitsDTGroup <- annoHitsDT[annoHitsDT[ , .I[which.min(minDistance)], by = peakHits]$V1]
  peakCategories <- annoHitsDTGroup[[featureColname]]

  peakCategoriesTable <- table(peakCategories)
  return(annoHitsDTGroup)
}

# ChIPSeeker + CCRE annotation
annotatePeakByOverlappingChIPSeekerCCRE <- function(peak, annotation, categories, featureColname="V10", txdb) {
  chipSeekerAnno <- annotatePeak(peak, TxDb=txdb, tssRegion=c(-1000, 1000), verbose=FALSE)
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

# RepeatMasker annotation
annotatepeakByOverlappingRepeatMasker <- function(peak, annotationRepeatMasker, repeatMaskerFeatures) {
  return(annotatePeakByOverlappingClosestFeature(peak, annotationRepeatMasker, repeatMaskerFeatures, featureColname = "X11"))
}

# ChromHMM annotation
annotatepeakByOverlappingChromHMM <- function(peak, annotationChromHMM, categoriesChromHMM, featureColname = "V4") {
  return(annotatePeakByOverlappingClosestFeature(peak, annotationChromHMM, categoriesChromHMM, featureColname = featureColname))
}


plotAnnoBar.data.frame.one.target <- function(anno.df,
                                              xlab="",
                                              ylab="Percentage(%)",
                                              title="Feature Distribution",
                                              categoryColumn = ".id",
                                              colorOption = 3,
                                              features) {
  anno.df$Feature <- factor(anno.df$Feature, levels = rev(levels(anno.df$Feature)))
  missingIndices <- c()
  print(anno.df)
  for (i in seq_along(features)) {
    feature <- features[i]
    if (!feature %in% anno.df$Feature) {
      missingIndices <- append(missingIndices, i)
    }
    else {
      if (sum(anno.df$Frequency[anno.df$Feature == feature]) == 0) {
        # missingIndices <- append(missingIndices, i)
      }
    }
  }
  anno.df <- anno.df[anno.df$Frequency != 0,]
  print(missingIndices)
  p <- ggplot(anno.df, aes_string(x = categoryColumn,
                                  fill = "Feature",
                                  y = "Frequency"))

  p <- p + geom_bar(stat="identity") + coord_flip() + theme_bw() + theme(axis.text=element_text(size=16))
  p <- p + ylab(ylab) + xlab(xlab) + ggtitle(title)

  if (categoryColumn == 1) {
    p <- p + scale_x_continuous(breaks=NULL)

    p <- p+scale_fill_manual(values=rev(getCols(nrow(anno.df), option = colorOption, missingIndices)), guide=guide_legend(reverse=TRUE))
  } else {
    col <- getCols(length(unique(anno.df$Feature)), option = colorOption, missingIndices)
    print(col)
    col <- col[!is.na(col)]
    print(col)
    p <- p+scale_fill_manual(values=rev(col), guide=guide_legend(reverse=TRUE))
  }
  # p <- p + theme(plot.title = element_blank(), plot.margin=unit(c(1,1,-0.5,1),  "cm"), panel.border = element_blank(), panel.background = element_rect(fill='transparent'), axis.text=element_text(size=16),
  #                  plot.background = element_rect(fill='transparent', color=NA),
  #                  panel.grid.major = element_blank(),
  #                  panel.grid.minor = element_blank(),
  #                  legend.background = element_rect(fill='transparent'),
  #                  legend.box.background = element_rect(fill='transparent'))+ labs(x=NULL, y=NULL)
  return(p)
}

getCols <- function(n, option = 3, missingIdx = NULL) {
  col <- c("#8dd3c7", "#ffffb3", "#bebada",
           "#fb8072", "#80b1d3", "#fdb462",
           "#b3de69", "#fccde5", "#d9d9d9",
           "#bc80bd", "#ccebc5", "#ffed6f")

  col2 <- c("#1f78b4", "#ffff33", "#c2a5cf",
            "#ff7f00", "#810f7c", "#a6cee3",
            "#006d2c", "#d73027", "#4d4d4d",
            "#8c510a", "#78c679", "#7f0000",
            "#41b6c4", "#e7298a", "#54278f")

  col3 <- c("#a6cee3", "#1f78b4", "#b2df8a",
            "#33a02c", "#fb9a99", "#e31a1c",
            "#fdbf6f", "#ff7f00", "#cab2d6",
            "#6a3d9a", "#ffff99", "#b15928")
  col4 <- c(col3, col2, col)
  col5 <- rainbow(120)
  col6 <- c(col2, col, col3)

  ## colorRampPalette(brewer.pal(12, "Set3"))(n)
  cols <- list(col, col2, col3, col4, col5, col6)
  col <- cols[[option]]
  # print(col)
  if (!is.null(missingIdx)) {
    col <- col[-missingIdx]
  }

  print(col)
  # print("n:")
  # print(n)
  return(col[1:n])
}
