#' Create a bar plot of TaxSEA results
#'
#' This function takes a TaxSEA result data frame and creates a bar plot
#' of the results
#' using ggplot2. It highlights Enriched and Depleted taxon sets
#' based on their #' normalized enrichment scores (NES) and false discovery
#' rate (FDR).Note this function assumes there are both enriched and
#' depleted sets.
#' @param taxsea_results A data frame containing the results of the TaxSEA
#' function
#' @param threshold Numeric value representing the FDR threshold for
#' displaying
#' taxon sets (default: 0.2)
#' @param custom_colors A character vector of length 2 with colors for
#' Enriched and Depleted taxon sets, respectively (default: NULL)
#' @return A ggplot2 bar plot of the TaxSEA results
#' @examples
#' data("TaxSEA_test_data")
#' taxsea_results <- TaxSEA(TaxSEA_test_data)
#' metabolite_plot <- TaxSEA_barplot(taxsea_results$Metabolite_producers)
#'
#' @import ggplot2
#' @export
TaxSEA_barplot <- function(taxsea_results, threshold = 0.2,
                           custom_colors = NULL) {
  # Check if ggplot2 is installed
  if (!requireNamespace("ggplot2", quietly = TRUE))
    cat("ggplot2 is not installed\n")

  # Check that there are more than four rows with FDR < 0.2
  if (sum(taxsea_results$FDR < threshold) <= 5)
    stop("There are very few taxon sets meeting the plotting threshold.
         I suggest viusalising using an alternative approach.")

  # Check that NES column contains both positive and negative values
  if (!(any(taxsea_results$NES < 0) && any(taxsea_results$NES > 0)))
    stop("There are only changes in one direction. This function assumes
         there are both enriched and depleted taxon sets.")

  # Extract relevant columns and calculate log10FDR
  taxsea_results$log10FDR <-
    -log10(taxsea_results$FDR) * ifelse(taxsea_results$NES < 0, -1, 1)

  # Sort the dataframe by NES sign and PValue
  taxsea_results <- taxsea_results[rev(order((taxsea_results$NES < 0),
                                             taxsea_results$PValue)), ]

  # Convert taxonSetName column to a factor with the current order
  taxsea_results$taxonSetName <-
    factor(taxsea_results$taxonSetName,
                                        levels =
             taxsea_results$taxonSetName)

  # Set default colors if custom_colors is NULL
  if (is.null(custom_colors)) {
    custom_colors <- c("#E74C3C", "#3498DB")
  }

  # Create a bar plot using ggplot2
  p <- ggplot2::ggplot(taxsea_results[taxsea_results$FDR < threshold, ],
                       aes(x = log10FDR,
                           y = taxonSetName, fill = NES < 0)) +
    geom_col(colour="black") +
    geom_vline(xintercept = c(-log10(0.1), log10(0.1)), linetype = 3) +
    scale_fill_manual(values = custom_colors,
                      labels = c("Enriched", "Depleted")) +
    theme_classic() +
    theme(
      legend.title = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      legend.text = element_text(size = 12),
      panel.grid.major = element_line(color = "#D3D3D3"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white")
    ) +
    xlab("-log10 FDR") +
    ggtitle("TaxSEA results")

  # Remove fill legend if there are no items with NES < 0
  if (all(taxsea_results$NES >= 0)) {
    p <- p + theme(legend.position = "none")
  }

  return(p+ scale_x_continuous(labels = function(x) ifelse(x < 0, -x, x)))
}



