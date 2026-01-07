biclustering_genomic_distribution <- function(row_cluster_file_path, out_dir = "./", distributions = c("genic", "ccre"), ref_genome = "hg38", ref_source = "knownGene", mode = "nearest", plot = TRUE) {
  # Validate parameters
  if (!ref_genome %in% c("hg38", "mm10")) {
    stop("ref_genome must be 'hg38' or 'mm10'")
  }
  if (!ref_source %in% c("knownGene", "GENCODE")) {
    stop("ref_source must be 'knownGene' or 'GENCODE'")
  }
  if (!mode %in% c("nearest", "weighted")) {
    stop("mode must be 'nearest' or 'weighted'")
  }

  valid_distributions <- c("genic", "ccre", "chromhmm", "repeat")
  invalid_annos <- setdiff(distributions, valid_distributions)
  if (length(invalid_annos) > 0) {
    stop(
      "Invalid annotation types: ", paste(invalid_annos, collapse = ", "),
      "\nValid options are: ", paste(valid_distributions, collapse = ", ")
    )
  }
  if (length(distributions) == 0) {
    stop(
      "distributions parameter must contain at least one annotation type: ",
      paste(valid_distributions, collapse = ", ")
    )
  }

  # Read row cluster file and convert to GRangesList
  row_cluster <- read.table(row_cluster_file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (!all(c("feature", "label") %in% colnames(row_cluster))) {
    stop("Input file must contain 'feature' and 'label' columns")
  }
  pos_df <- do.call(rbind, strsplit(row_cluster$feature, "_"))
  colnames(pos_df) <- c("seqnames", "start", "end")
  row_cluster <- cbind(pos_df, row_cluster)
  row_gr <- makeGRangesFromDataFrame(
    row_cluster,
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
  )
  row_grl <- split(row_gr, row_gr$label)
  message("Loaded ", length(row_gr), " regions across ", length(row_grl), " clusters")
  message("Clusters: ", paste(names(row_grl), collapse = ", "))

  genomic_distribution(query = row_grl, out_dir = out_dir, distributions = distributions, ref_genome = ref_genome, ref_source = ref_source, mode = mode, plot = plot)
  message("Distribution annotation complete")
}
