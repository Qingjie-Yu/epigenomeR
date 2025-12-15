#' ### Parameters
#'
#' | Parameter | Type | Required | Description | Example |
#' |-----------|------|----------|-------------|---------|
#' | `row_cluster_file_path` | character | Yes | Path to TSV file containing cluster assignments. Must have columns: 'feature' (genomic coordinates in format chr_start_end) and 'label' (cluster ID). Typically uses row_table_clean.tsv from biclustering results | `row_cluster_file_path = "./row_table_clean.tsv"` |
#' | `out_dir` | character | No (default: "./") | Directory to save all annotation outputs. Subdirectories will be created for each annotation type | `out_dir = "./annotations"` |
#' | `ref_genome` | character | No (default: "hg38") | Reference genome version. Must be either "hg38" (Human GRCh38) or "mm10" (Mouse GRCm38) | `ref_genome = "hg38"` |
#' | `ref_source` | character | No (default: "knownGene") | Gene annotation source for cCRE annotation. Options: "knownGene" (UCSC knownGene) or "GENCODE" (GENCODE annotations). Only used if "ccre" is in annotations parameter | `ref_source = "knownGene"` |
#' | `mode` | character | No (default: "nearest") | Annotation mode for all annotation types. Options: "nearest" (assigns each region to closest feature) or "weighted" (proportional assignment by overlap length) | `mode = "nearest"` |
#' | `annotations` | character vector | No (default: c("ccre", "chromhmm", "repeat")) | Vector specifying which annotation types to perform. Valid options: "ccre" (cCRE and gene features), "chromhmm" (chromatin states), "repeat" (repeat elements). Can specify any combination | `annotations = c("ccre", "chromhmm")` |
#' | `plot` | logical | No (default: TRUE) | Whether to generate stacked barplot visualizations for each annotation type | `plot = TRUE` |
#'
#' 
biclustering_genomic_distribution <- function(row_cluster_file_path, out_dir = "./", ref_genome = "hg38", ref_source = "knownGene", mode = "nearest", annotations = c("ccre")) {
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

  valid_annotations <- c("ccre", "chromhmm", "repeat")
  invalid_annos <- setdiff(annotations, valid_annotations)
  if (length(invalid_annos) > 0) {
    stop("Invalid annotation types: ", paste(invalid_annos, collapse = ", "), 
         "\nValid options are: ", paste(valid_annotations, collapse = ", "))
  }
  if (length(annotations) == 0) {
    stop("annotations parameter must contain at least one annotation type: ", 
         paste(valid_annotations, collapse = ", "))
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

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

 # Perform selected annotations
  if ("ccre" %in% annotations) {
    message("\n========== Running cCRE annotation ==========")
    annotation_ccre(
      query_grl = row_grl,
      out_dir = out_dir,
      ref_genome = ref_genome,
      ref_source = ref_source,
      mode = mode,
      plot = plot
    )
    message("cCRE annotation complete. Results saved to: ", out_dir)
  }
  
  if ("chromhmm" %in% annotations) {
    message("\n========== Running ChromHMM annotation ==========")
    annotation_chromhmm(
      query_grl = row_grl,
      out_dir = out_dir,
      ref_genome = ref_genome,
      mode = mode,
      plot = plot
    )
    message("ChromHMM annotation complete. Results saved to: ", out_dir)
  }
  
  if ("repeat" %in% annotations) {
    message("\n========== Running Repeat annotation ==========")
    annotation_repeat(
      query_grl = row_grl,
      out_dir = out_dir,
      ref_genome = ref_genome,
      mode = mode,
      plot = plot
    )
    message("Repeat annotation complete. Results saved to: ", out_dir)
  }
  
}