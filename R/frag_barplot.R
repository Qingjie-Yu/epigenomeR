# Post: Generate comprehensive fragment length distribution bar plots showing nucleosome occupancy patterns across samples, with options for top performers and full sample visualization.
# Parameter: frag_decomp_file: Vector of file paths to TSV files containing fragment percentage data (output from df_for_barplot_with_perc)
#            target_pair_list: Vector of sample names to include in analysis
#            target_pair_mapping_df: Optional data frame for mapping sample names to display names (default: NULL)
#            out_dir: Output directory path to save generated plots
# Output: Saves multiple PDF plots including top 20 samples by category, summary plots, and full sample plots; returns list of ggplot objects
frag_barplot <- function(frag_decomp_file, target_pair_list, target_pair_mapping_df = NULL, out_dir) {
  # load library
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(glue)
    library(latex2exp)
    library(cowplot)
  })
  print(out_dir)
  print(dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  print(out_dir)
  print(dir.exists(out_dir))

  # 1. Top20 plots for each size_category
  plts <- list()
  for (i in seq_along(frag_decomp_file)) {
    basenames <- tools::file_path_sans_ext(basename(frag_decomp_file[i]))
    frag_len_percent <- read.table(frag_decomp_file[i], sep = "\t", header = TRUE, row.names = 1)
    frag_len_percent_filtered <- frag_len_percent[target_pair_list, ]
    colnames(frag_len_percent_filtered)[3] <- "dimer+"
    plts[[basenames]] <- list()

    for (size_category in colnames(frag_len_percent_filtered)) {
      frag_len_percent_filtered_ordered <- frag_len_percent_filtered[
        order(frag_len_percent_filtered[[size_category]], decreasing = TRUE),
      ]
      frag_len_percent_filtered_ordered_top50 <- frag_len_percent_filtered_ordered[1:20, ]

      frag_len_list <- list()
      for (tag in rownames(frag_len_percent_filtered_ordered_top50)) {
        percentage <- as.numeric(as.vector(frag_len_percent_filtered_ordered_top50[tag, ]))
        frag_len_list[[tag]] <- data.frame(
          type = colnames(frag_len_percent_filtered_ordered_top50),
          percentage = percentage
        )
      }

      barplot_df <- list_to_dataframe(frag_len_list)
      barplot_df$.id <- map_target_names(barplot_df$.id, target_pair_mapping_df)
      barplot_df$.id <- factor(
        barplot_df$.id,
        levels = rev(map_target_names(rownames(frag_len_percent_filtered_ordered_top50), target_pair_mapping_df))
      )
      barplot_df$type <- factor(barplot_df$type, levels = rev(c("subnucleo", "monomer", "dimer+")))

      p <- ggplot(barplot_df, aes(fill = type, y = percentage, x = .id)) +
        scale_x_discrete(labels = TeX) +
        geom_bar(position = "fill", stat = "identity") +
        theme_gray(base_size = 18) +
        scale_fill_manual(values = c("#00b0be", "#ff8ca1", "#9ac9db", "#FBE7C6")) +
        theme_classic() +
        xlab("type") +
        ylab("Peak Number") +
        theme_bw() +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size = 2),
          axis.title.y = element_blank()
        ) +
        ggtitle(glue("Fragment length distribution with \nmost {size_category}")) +
        labs(color = "") +
        guides(col = guide_legend(override.aes = list(linetype = 2))) +
        labs(
          fill = "",
          y = "Reads percentage",
          x = "condition"
        ) +
        coord_flip()

      plts[[basenames]][[size_category]] <- p
      pdf(glue("{out_dir}/barplot_{basenames}_sort_on_{size_category}.pdf"), height = 5, width = 8)
      print(p)
      dev.off()
    }
  }

  # 2. Summary plot with legend
  legend <- get_legend(p)
  for (i in seq_along(frag_decomp_file)) {
    basenames <- tools::file_path_sans_ext(basename(frag_decomp_file[i]))
    pdf(glue("{out_dir}/barplot_{basenames}_summary.pdf"), height = 6, width = 14)
    print(plot_grid(
      plts[[basenames]][["subnucleo"]] + theme(legend.position = "none"),
      plts[[basenames]][["monomer"]] + theme(legend.position = "none"),
      plts[[basenames]][["dimer+"]] + theme(legend.position = "none"),
      legend,
      align = "h",
      axis = "bt",
      ncol = 4,
      rel_widths = c(1, 1, 1, 0.3)
    ))
    dev.off()
  }

  # ------------------------
  # 3. Full plot for all samples
  # ------------------------
  plts <- list()
  for (i in seq_along(frag_decomp_file)) {
    basenames <- tools::file_path_sans_ext(basename(frag_decomp_file[i]))
    frag_len_percent <- read.table(frag_decomp_file[i], sep = "\t", header = TRUE, row.names = 1)
    frag_len_percent_filtered <- frag_len_percent[target_pair_list, ]
    colnames(frag_len_percent_filtered)[3] <- "dimer+"
    plts[[basenames]] <- list()

    size_category <- colnames(frag_len_percent_filtered)[1]
    frag_len_percent_filtered_ordered <- frag_len_percent_filtered[
      order(frag_len_percent_filtered[[size_category]], decreasing = TRUE),
    ]

    frag_len_list <- list()
    for (tag in rownames(frag_len_percent_filtered_ordered)) {
      percentage <- as.numeric(as.vector(frag_len_percent_filtered_ordered[tag, ]))
      frag_len_list[[tag]] <- data.frame(
        type = colnames(frag_len_percent_filtered_ordered),
        percentage = percentage
      )
    }

    barplot_df <- list_to_dataframe(frag_len_list)
    barplot_df$.id <- map_target_names(barplot_df$.id, target_pair_mapping_df)
    barplot_df$.id <- factor(
      barplot_df$.id,
      levels = rev(map_target_names(rownames(frag_len_percent_filtered_ordered), target_pair_mapping_df))
    )
    barplot_df$type <- factor(barplot_df$type, levels = rev(c("subnucleo", "monomer", "dimer+")))

    write.csv(barplot_df, glue("{out_dir}/barplot_{basenames}_sort_on_all.csv"), quote = FALSE, row.names = FALSE)

    p <- ggplot(barplot_df, aes(fill = type, y = percentage, x = .id)) +
      scale_x_discrete(labels = TeX) +
      geom_bar(position = "fill", stat = "identity") +
      theme_gray(base_size = 18) +
      scale_fill_manual(values = c("#00b0be", "#ff8ca1", "#9ac9db", "#FBE7C6")) +
      theme_classic() +
      xlab("type") +
      ylab("Peak Number") +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size = 2),
        axis.title.y = element_blank()
      ) +
      ggtitle(glue("Fragment length distribution")) +
      labs(color = "") +
      guides(col = guide_legend(override.aes = list(linetype = 2))) +
      labs(
        fill = "",
        y = "Reads percentage",
        x = "condition"
      ) +
      coord_flip()

    plts[[basenames]][[size_category]] <- p
    pdf(glue("{out_dir}/barplot_{basenames}_sort_on_all.pdf"), height = 120, width = 8)
    print(p)
    dev.off()
  }

  invisible(plts)
}
