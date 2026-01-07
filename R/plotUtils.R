library(ggplot2)
library(ChIPseeker)
library("dplyr")
library(ComplexHeatmap)
library(grid)
library(stats)
library(latex2exp)

# allFeatures <- c("Promoter", "5' UTR", "3' UTR", "1st Exon", "Other Exon", "1st Intron", "Other Intron", "Distal Intergenic", "Downstream (<=300)")

allFeaturesNoEnhancer <- c("Distal Intergenic", "Promoter", "5' UTR", "1st Exon", "1st Intron", "Other Exon", "Other Intron", "3' UTR", "Downstream (<=300)")
allFeaturesEnhancer <- c("Enhancer", "Other Distal Intergenic", "Promoter", "5' UTR", "1st Exon", "1st Intron", "Other Exon", "Other Intron", "3' UTR", "Downstream (<=300)")
allFeatures <- c("Distal Intergenic", "Promoter", "5' UTR", "1st Exon", "1st Intron", "Other Exon", "Other Intron", "3' UTR", "Downstream (<=300)")

getAnnoStat <- function(x) {
  if (!is(x, "csAnno")) {
    stop("not supported...")
  }
  return(x@annoStat)
}

list_to_dataframe <- function(dataList) {
  if (is.null(names(dataList))) {
    return(do.call("rbind", dataList))
  }

  cn <- lapply(dataList, colnames) %>%
    unlist() %>%
    unique()
  cn <- c(".id", cn)
  dataList2 <- lapply(seq_along(dataList), function(i) {
    data <- dataList[[i]]
    data$.id <- names(dataList)[i]
    idx <- !cn %in% colnames(data)
    if (sum(idx) > 0) {
      for (i in cn[idx]) {
        data[, i] <- NA
      }
    }
    return(data[, cn])
  })
  res <- do.call("rbind", dataList2)
  res$.id <- factor(res$.id, levels = rev(names(dataList)))
  return(res)
}


getCols <- function(n, option = 3, missingIdx = NULL) {
  col <- c(
    "#8dd3c7", "#ffffb3", "#bebada",
    "#fb8072", "#80b1d3", "#fdb462",
    "#b3de69", "#fccde5", "#d9d9d9",
    "#bc80bd", "#ccebc5", "#ffed6f"
  )

  col2 <- c(
    "#1f78b4", "#ffff33", "#c2a5cf",
    "#ff7f00", "#810f7c", "#a6cee3",
    "#006d2c", "#d73027", "#4d4d4d",
    "#8c510a", "#78c679", "#7f0000",
    "#41b6c4", "#e7298a", "#54278f"
  )

  col3 <- c(
    "#a6cee3", "#1f78b4", "#b2df8a",
    "#33a02c", "#fb9a99", "#e31a1c",
    "#fdbf6f", "#ff7f00", "#cab2d6",
    "#6a3d9a", "#ffff99", "#b15928"
  )
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


plotAnnoBar.data.frame <- function(anno.df,
                                   xlab = "",
                                   ylab = "Percentage(%)",
                                   title = "Feature Distribution",
                                   categoryColumn,
                                   colorOption = 3, width = NULL) {
  anno.df$Feature <- factor(anno.df$Feature, levels = rev(levels(anno.df$Feature)))

  p <- ggplot(anno.df, aes_string(
    x = categoryColumn,
    fill = "Feature",
    y = "Frequency"
  )) +
    scale_x_discrete(labels = TeX)

  p <- p + geom_bar(stat = "identity", width = width) + coord_flip() + theme_bw() + theme(axis.text = element_text(size = 16))
  p <- p + ylab(ylab) + xlab(xlab) + ggtitle(title)

  if (categoryColumn == 1) {
    p <- p + scale_x_continuous(breaks = NULL)
    p <- p + scale_fill_manual(values = rev(getCols(nrow(anno.df), option = colorOption)), guide = guide_legend(reverse = TRUE))
  } else {
    p <- p + scale_fill_manual(values = rev(getCols(length(unique(anno.df$Feature)), option = colorOption)), guide = guide_legend(reverse = TRUE))
  }
  return(p)
}


plotAnnoBar.data.frame.given.color <- function(anno.df,
                                               xlab = "",
                                               ylab = "Percentage(%)",
                                               title = "Feature Distribution",
                                               categoryColumn,
                                               colors, width = NULL) {
  anno.df$Feature <- factor(anno.df$Feature, levels = rev(levels(anno.df$Feature)))

  p <- ggplot(anno.df, aes_string(
    x = categoryColumn,
    fill = "Feature",
    y = "Frequency"
  ))

  p <- p + geom_bar(stat = "identity", width = width) + coord_flip() + theme_bw() + theme(axis.text = element_text(size = 16))
  p <- p + ylab(ylab) + xlab(xlab) + ggtitle(title)

  if (categoryColumn == 1) {
    p <- p + scale_x_continuous(breaks = NULL)
    p <- p + scale_fill_manual(values = rev(colors), guide = guide_legend(reverse = TRUE))
  } else {
    p <- p + scale_fill_manual(values = rev(colors), guide = guide_legend(reverse = TRUE))
  }
  return(p)
}


plotAnnoBar.data.frame.one.target <- function(anno.df,
                                              xlab = "",
                                              ylab = "Percentage(%)",
                                              title = "Feature Distribution",
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
    } else {
      if (sum(anno.df$Frequency[anno.df$Feature == feature]) == 0) {
        # missingIndices <- append(missingIndices, i)
      }
    }
  }
  anno.df <- anno.df[anno.df$Frequency != 0, ]
  print(missingIndices)
  p <- ggplot(anno.df, aes_string(
    x = categoryColumn,
    fill = "Feature",
    y = "Frequency"
  ))

  p <- p + geom_bar(stat = "identity") + coord_flip() + theme_bw() + theme(axis.text = element_text(size = 16))
  p <- p + ylab(ylab) + xlab(xlab) + ggtitle(title)

  if (categoryColumn == 1) {
    p <- p + scale_x_continuous(breaks = NULL)

    p <- p + scale_fill_manual(values = rev(getCols(nrow(anno.df), option = colorOption, missingIndices)), guide = guide_legend(reverse = TRUE))
  } else {
    col <- getCols(length(unique(anno.df$Feature)), option = colorOption, missingIndices)
    print(col)
    col <- col[!is.na(col)]
    print(col)
    p <- p + scale_fill_manual(values = rev(col), guide = guide_legend(reverse = TRUE))
  }
  # p <- p + theme(plot.title = element_blank(), plot.margin=unit(c(1,1,-0.5,1),  "cm"), panel.border = element_blank(), panel.background = element_rect(fill='transparent'), axis.text=element_text(size=16),
  #                  plot.background = element_rect(fill='transparent', color=NA),
  #                  panel.grid.major = element_blank(),
  #                  panel.grid.minor = element_blank(),
  #                  legend.background = element_rect(fill='transparent'),
  #                  legend.box.background = element_rect(fill='transparent'))+ labs(x=NULL, y=NULL)
  return(p)
}

plotAnnoBar.data.frame.red.even.rows <- function(anno.df,
                                                 xlab = "",
                                                 ylab = "Percentage(%)",
                                                 title = "Feature Distribution",
                                                 categoryColumn,
                                                 colorOption = 3) {
  anno.df$Feature <- factor(anno.df$Feature, levels = rev(levels(anno.df$Feature)))

  anno.df$outline <- ifelse(grepl("treated", anno.df[[categoryColumn]]), "red", NA)

  # p <- ggplot(anno.df, aes_string(x = categoryColumn,
  #                                 fill = "Feature",
  #                                 y = "Frequency",
  #                                 col = "outline"))
  p <- ggplot(anno.df)
  p <- p + geom_bar(aes(
    x = .id,
    fill = Feature,
    y = Frequency,
    col = outline
  ), stat = "identity") + coord_flip() + theme_bw() + theme(axis.text = element_text(size = 16))
  p <- p + ylab(ylab) + xlab(xlab) + ggtitle(title) + geom_bar(aes(x = .id, y = Frequency, group = .id),
    stat = "summary",
    fun = sum,
    fill = "transparent"
  ) + scale_color_identity()


  if (categoryColumn == 1) {
    p <- p + scale_x_continuous(breaks = NULL)
    p <- p + scale_fill_manual(values = rev(getCols(nrow(anno.df), option = colorOption)), guide = guide_legend(reverse = TRUE))
  } else {
    p <- p + scale_fill_manual(values = rev(getCols(length(unique(anno.df$Feature)), option = colorOption)), guide = guide_legend(reverse = TRUE))
  }
  p <- p
  return(p)
}


getAnnotationMatrix <- function(anno, allFeatures) {
  resultList <- list()
  for (f in allFeatures) {
    resultList[[f]] <- c()
  }
  for (comb in names(anno)) {
    features <- as.character(anno[[comb]]$Feature)
    frequencies <- as.numeric(anno[[comb]]$Frequency)
    names(frequencies) <- features
    frequencies <- as.list(frequencies)
    missingFeatures <- allFeatures[!allFeatures %in% features]
    for (f in missingFeatures) {
      frequencies[[f]] <- 0
    }
    for (f in allFeatures) {
      resultList[[f]] <- append(resultList[[f]], frequencies[[f]])
    }
  }
  resultDF <- as.data.frame(resultList, check.names = FALSE)
  rownames(resultDF) <- names(anno)
  resultMat <- as.matrix(resultDF)
  return(resultMat)
}

getEnhancerAnnoStats <- function(xT, ehancerDb) {
  annoStats <- list()
  for (tag in names(xT)) {
    distalElem <- xT[[tag]]@anno[xT[[tag]]@anno$annotation == "Distal Intergenic"]
    enhancerPercentage <- length(unique(queryHits(findOverlaps(distalElem, ehancerDb)))) / length(distalElem)
    annoStat <- getAnnoStat(xT[[tag]])
    if (!is.nan(enhancerPercentage)) {
      distalFreq <- annoStat$Frequency[annoStat$Feature == "Distal Intergenic"]
      enhancerFreq <- distalFreq * enhancerPercentage
      otherDistalFreq <- distalFreq * (1 - enhancerPercentage)
      annoStat <- rbind(annoStat, as.data.frame(list(Feature = "Enhancer", Frequency = enhancerFreq)))
      annoStat <- rbind(annoStat, as.data.frame(list(Feature = "Other Distal Intergenic", Frequency = otherDistalFreq)))
      annoStat <- annoStat[(annoStat$Feature != "Distal Intergenic"), ]
      annoStats[[tag]] <- annoStat
    }
  }
  return(annoStats)
}


filterByChisq <- function(resultMatT, resultMatV, xT, xV, cutoff) {
  tags <- c()
  for (tag in rownames(resultMatT)) {
    resultTagT <- resultMatT[tag, ] / 100 * length(xT[[tag]]@anno)
    resultTagV <- resultMatV[tag, ] / 100 * length(xV[[tag]]@anno)
    a <- chisq.test(rbind(resultTagT, resultTagV))
    if (!is.nan(a$p.value)) {
      if (a$p.value < 0.05) {
        tags <- append(tags, tag)
        cat(tag, " ", a$p.value, "\n")
      }
    }
  }
  return(tags)
}


getRevHeterotone <- function(tag) {
  tag <- as.character(tag)
  tags <- strsplit(tag, "-")[[1]]
  tag1 <- tags[1]
  tag2 <- tags[2]
  return(paste0(tag2, "-", tag1))
}


getUniqueHeterotones <- function(anno.df, tag) {
  unqiueHeterotones <- c()
  for (comb in unique(anno.df$.id)) {
    if (substring(comb, 1, nchar(tag)) == tag) {
      unqiueHeterotones <- append(unqiueHeterotones, comb)
    } else {
      unqiueHeterotones <- append(unqiueHeterotones, getRevHeterotone(comb))
    }
  }
  unqiueHeterotones <- unique(unqiueHeterotones)
  return(unqiueHeterotones)
}

getAnnoStatForUniqueHeterotones <- function(x, allFeatures, unqiueHeterotones) {
  newAnnoStatV <- list()
  for (tag in unqiueHeterotones) {
    annoStat <- x[[tag]]@annoStat
    peakNum <- x[[tag]]@peakNum
    annoStat$Frequency <- annoStat$Frequency * peakNum
    revAnnoStat <- x[[getRevHeterotone(tag)]]@annoStat
    revPeakNum <- x[[getRevHeterotone(tag)]]@peakNum
    revAnnoStat$Frequency <- revAnnoStat$Frequency * revPeakNum
    freq <- list()
    for (f in allFeatures) {
      freq[[f]] <- 0
    }
    for (f in annoStat$Feature) {
      freq[[f]] <- freq[[f]] + annoStat$Frequency[annoStat$Feature == f]
    }
    for (f in revAnnoStat$Feature) {
      freq[[f]] <- freq[[f]] + revAnnoStat$Frequency[revAnnoStat$Feature == f]
    }
    newannoStat <- data.frame(Feature = names(freq), Frequency = (as.numeric(freq) / sum(as.numeric(freq)) * 100))
    newAnnoStatV[[tag]] <- newannoStat
  }
  return(newAnnoStatV)
}

extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}
