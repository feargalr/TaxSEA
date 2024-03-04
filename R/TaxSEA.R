#' Retrieve NCBI Taxonomy IDs for a list of taxon names
#'
#' This function takes a vector of taxon names and returns a vector of 
#' NCBI taxonomy IDs
#' by querying the NCBI Entrez API.
#'
#' @param taxon_names A character vector of taxon names
#' @return A character vector of NCBI taxonomy IDs corresponding to the 
#' input taxon names
#' @examples
#' taxon_names <- c("Escherichia coli", "Staphylococcus aureus", 
#' "Bacillus subtilis")
#' taxon_ids <- get_ncbi_taxon_ids(taxon_names)
#' print(taxon_ids)
#' @import data/NCBI_ids.rds
#' @export
get_ncbi_taxon_ids <- function(taxon_names) {
  get_ncbi_taxon_id <- function(taxon_name) {
    # Define the URL for the API call
    base_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    esearch_url <- paste0(base_url, "esearch.fcgi?db=taxonomy&term=",
                          URLencode(taxon_name, reserved = TRUE), 
                          "&retmode=xml")

    # Retrieve the XML data
    xml_data <- readLines(url(esearch_url, "r"))

    # Find the line containing the taxonomy ID
    id_line <- grep("<Id>", xml_data, value = TRUE)

    # Extract the taxonomy ID
    taxon_id <- gsub("<Id>|</Id>", "", id_line)

    # Return the taxonomy ID
    if (length(taxon_id) > 0) {
      return(taxon_id)
    } else {
      warning(paste("No taxonomy ID found for taxon:", taxon_name))
      return(NULL)
    }
  }

  NCBI_ids <- readRDS(system.file("data/NCBI_ids.rds", package = "TaxSEA"))
  ids2fetch = taxon_names[!taxon_names %in% names(NCBI_ids)]
  taxon_names = taxon_names[taxon_names %in% names(NCBI_ids)]
  local_ids = unlist(NCBI_ids[taxon_names])
  # Retrieve the taxonomy IDs for each taxon name in the vector
  if (length(ids2fetch) > 0){
    print("Fetching some NCBI IDs for input taxa. Please wait")
    taxon_ids <- sapply(ids2fetch, get_ncbi_taxon_id)
    taxon_ids = c(local_ids,unlist(fetched_ids))
  } else {
    taxon_ids = local_ids
  }

  return(taxon_ids)
  close(xml_data)

}



#' TaxSEA: Taxon Set Enrichment Analysis
#'
#' TaxSEA is designed to enable rapid annotation of changes observed in a 
#' microbiome association study
#' by testing for enrichment for producers of particular metabolites, or 
#' associations with marker taxa
#' for particular diseases. It focuses specifically on human gut microbiome 
#' associations and uses
#' a Kolmogorov-Smirnov test to test if a particular set of taxa is changed 
#' relative to a control group.
#' The input taxon_ranks are log2 fold changes between control and test group 
#' (e.g., healthy and IBD).
#'
#' @param taxon_ranks A named vector of log2 fold changes between control and 
#' test group.
#' @param lookup_missing Logical value indicating whether to fetch missing 
#' NCBI IDs. Default is FALSE.
#' @param min_set_size Minimum size of taxon sets to include in the analysis. 
#' Default is 5.
#' @param max_set_size Maximum size of taxon sets to include in the analysis. 
#' Default is 100.
#' @param database A character specifying the database to use for enrichment 
#' analysis.
#'   Options are "All", "GutMGene", "MiMeDB", and "GMRepoV2". Default is "All".
#' @return A data frame with taxon set enrichment results.
#' @seealso
#' \itemize{
#'   \item \url{https://doi.org/10.1093/nar/gkac868} for MiMeDB
#'   \item \url{https://doi.org/10.1093/nar/gkab1019} for GMrepo
#'   \item \url{https://doi.org/10.1093/nar/gkab786} for gutMGene
#' }
#' @examples
#' # Example data
#' taxon_ranks <- runif(10, -3, 3)
#' names(taxon_ranks) <- paste("Taxon", 1:10)
#' taxon_sets <- list(
#'   set1 = c("Taxon 1", "Taxon 2", "Taxon 3"),
#'   set2 = c("Taxon 4", "Taxon 5", "Taxon 6"),
#'   set3 = c("Taxon 7", "Taxon 8", "Taxon 9")
#' )
#' # Run TaxSEA
#' result_df <- TaxSEA(taxon_ranks)
#' @import data/NCBI_ids.rds
#' @export


TaxSEA <- function(taxon_ranks, database = "All", lookup_missing = FALSE,
                   min_set_size = 5, max_set_size = 100) {
  # Function implementation
}


TaxSEA <- function(taxon_ranks, database = "All",lookup_missing = FALSE,
                   min_set_size = 5, max_set_size = 100) {
  NCBI_ids <- readRDS(system.file("data/NCBI_ids.rds", package = "TaxSEA"))

  if(length(taxon_ranks) < 5) {
    stop("Error: Very few taxa provided. Unadvisable to continue. Stopping")
  }
  # Select database

  if (database == "All") {
    taxon_sets <- TaxSEA_db
  } else if (!database %in% c("GutMGene", "MiMeDB", "GMRepoV2")) {
    stop("INCORRECT DATABASE SPECIFIED")
  } else {
    taxon_sets <- TaxSEA_db[grep(database, names(TaxSEA_db))]
  }

  if (lookup_missing == TRUE) {

  #Relabel taxa with NCBI IDs
  ids2fetch = names(taxon_ranks[!(names(taxon_ranks) %in% names(NCBI_ids))])
  if (length(ids2fetch) > 0){
    print("Fetching some NCBI IDs for input taxa. Please wait")
    fetched_ids = suppressWarnings(get_ncbi_taxon_ids(ids2fetch))
    if (length(unlist(fetched_ids)) > 0){
      NCBI_ids = c(NCBI_ids,unlist(fetched_ids))
    }
  }
}

  taxon_ranks = taxon_ranks[names(taxon_ranks) %in% names(NCBI_ids)]
  original_ranks = taxon_ranks
  names(taxon_ranks) = NCBI_ids[names(taxon_ranks)]
  NCBI_ids = NCBI_ids[names(NCBI_ids) %in% names(original_ranks)]
  NCBI_ids = NCBI_ids[!duplicated(NCBI_ids)]
  legible_names = names(NCBI_ids)
  names(legible_names) = NCBI_ids[legible_names]

  if(length(taxon_ranks) < 3) {
    stop("Error: Very few taxa provided. Unadvisable to continue. Stopping")
  }

  # Filter taxon_ranks and taxon_sets
  taxon_ranks = taxon_ranks[names(taxon_ranks) %in% unique(unlist(taxon_sets))]
  taxon_sets <- lapply(taxon_sets, function(taxon_set) {
    return(unique(taxon_set[taxon_set %in% names(taxon_ranks)]))})
  set_sizes <- sapply(taxon_sets, function(taxon_set) { length(taxon_set) })
  taxon_sets = taxon_sets[set_sizes>= min_set_size & set_sizes <= max_set_size]
  taxon_ranks = taxon_ranks[names(taxon_ranks) %in% unique(unlist(taxon_sets))]

  if (!is.vector(taxon_ranks) || !is.list(taxon_sets)) {
    stop("Input error")
  }

  taxon_sets <- lapply(taxon_sets, function(taxon_set) {
    intersect(names(taxon_ranks), taxon_set)
  })

  # Perform KS tests
  ks_results <- lapply(taxon_sets, function(taxon_set) {
    taxon_set_ranks <- taxon_ranks[taxon_set]


    nes <- median(taxon_ranks[taxon_set])

    ks_result <- suppressWarnings(ks.test(taxon_set_ranks, taxon_ranks))

    return(list(nes = nes, ks_result = ks_result))
  })
  taxon_sets = lapply(taxon_sets,function(X){legible_names[X]})
  result_df <- data.frame(
    taxonSetName = names(taxon_sets),
    NES = sapply(ks_results, function(res) res$nes),
    PValue = sapply(ks_results, function(res) res$ks_result$p.value),
    FDR = p.adjust(sapply(ks_results, function(res) res$ks_result$p.value),
                   method = "fdr"),
    TaxonSet = sapply(taxon_sets, function(taxon_set) paste(taxon_set,
                                                            collapse = ", "))
  )
  result_df = result_df[order(result_df$PValue,decreasing = FALSE),]
  metabolites_df = results_df[grepl("producers_of",results_df),]
  disease_df = results_df[!grepl("producers_of",results_df),]
  metabolites_df$FDR = p.adjust(metabolites_df$PValue,method = "fdr")
  disease_df$FDR = p.adjust(disease_df$PValue,method = "fdr")
  
  res_list = list(Metabolite_producers = metabolites_df,
                  Health_associations = disease_df)
  return(res_list)
  rm(TaxSEA_db)
  rm(NCBI_ids)

}


########################################################################
###############################################################################
#' Create a bar plot of TaxSEA results
#'
#' This function takes a TaxSEA result data frame and creates a bar plot of 
#' the results
#' using ggplot2. It highlights upregulated and downregulated taxon sets 
#' based on their #' normalized enrichment scores (NES) and false discovery 
#' rate (FDR).
#'
#' @param taxsea_results A data frame containing the results of the TaxSEA 
#' function
#' @param threshold Numeric value representing the FDR threshold for displaying 
#' taxon sets (default: 0.2)
#' @param custom_colors A character vector of length 2 with colors for 
#' upregulated and downregulated taxon sets, respectively (default: NULL)
#' @return A ggplot2 bar plot of the TaxSEA results
#' @examples
#' # Assuming you have already run the TaxSEA function and have a result 
#' data frame named 'taxsea_results'
#' p <- TaxSEA_barplot(taxsea_results)
#' print(p)
#' @import ggplot2
#' @export
TaxSEA_barplot <- function(taxsea_results, threshold = 0.2, 
                           custom_colors = NULL) {
  # Check if ggplot2 is installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) 
    cat("ggplot2 is not installed\n")
  require(ggplot2)

  # Extract relevant columns and calculate log10FDR
  taxsea_results$log10FDR <- 
    -log10(taxsea_results$FDR) * ifelse(taxsea_results$NES < 0, -1, 1)

  # Sort the dataframe by NES sign and PValue
  taxsea_results <- taxsea_results[rev(order((taxsea_results$NES < 0), 
                                             taxsea_results$PValue)), ]

  # Convert taxonSetName column to a factor with the current order
  taxsea_results$taxonSetName <- factor(taxsea_results$taxonSetName, 
                                        levels = taxsea_results$taxonSetName)

  # Set default colors if custom_colors is NULL
  if (is.null(custom_colors)) {
    custom_colors <- c("#E74C3C", "#3498DB")
  }

  # Create a bar plot using ggplot2
  p <- ggplot(taxsea_results[taxsea_results$FDR < threshold, ], 
              aes(x = log10FDR, y = taxonSetName, fill = NES < 0)) +
    geom_col(colour="black") +
    geom_vline(xintercept = c(-log10(0.1), log10(0.1)), linetype = 3) +
    scale_fill_manual(values = custom_colors, 
                      labels = c("Upregulated", "Downregulated")) +
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

  return(p)
}
###############################################################################
#' Create a barcode plot of TaxSEA results
#'
#' This function takes a TaxSEA result data frame, taxon ranks, and taxon 
#' sets to create a barcode plot
#' using ggplot2. The barcode plot shows the taxon sets with the top 
#' significant P values.
#'
#' @param taxsea_results A data frame containing the results of the 
#' TaxSEA function
#' @param taxon_ranks A named numeric vector with taxon names as names 
#' and their corresponding ranks as values
#' @param taxon_sets A named list of character vectors with taxon set 
#' names as names and taxon names as values
#' @param axis_limits A numeric vector of length 2 specifying the x-axis
#'  limits (default: c(-7, 7))
#' @param n_to_plot Integer value specifying the number of top taxon
#'  sets to display (default: 10)
#' @param boxplot Boolean value specifying whether or not to include a 
#' boxplot (default: FALSE)

#' @return A ggplot2 barcode plot of the TaxSEA results
#' @examples
#' # Assuming you have already run the TaxSEA function and have a result
#'  data frame named 'taxsea_results', taxon ranks 'taxon_ranks', and taxon sets 'taxon_sets'
#' plot <- barcode_plot(taxsea_results, taxon_ranks, taxon_sets)
#' print(plot)
#' @import ggplot2
#' @export
TaxSEA_barcode <- function(taxsea_results, taxon_ranks,
                           axis_limits = c(-7,7),n_to_plot = 10,boxplot=FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) 
    cat("ggplot2 is not installed\n")
  require(ggplot2)

    taxon_sets <- TaxSEA_db

  #Relabel taxa with NCBI IDs
  ids2fetch = names(taxon_ranks[!(names(taxon_ranks) %in% names(NCBI_ids))])
  if (length(ids2fetch) > 0){
    print("Fetching some NCBI IDs for input taxa. Please wait")
    fetched_ids = suppressWarnings(get_ncbi_taxon_ids(ids2fetch))
    if (length(unlist(fetched_ids)) > 0){
      NCBI_ids = c(NCBI_ids,unlist(fetched_ids))
    }
  }


  taxon_ranks = taxon_ranks[names(taxon_ranks) %in% names(NCBI_ids)]
  original_ranks = taxon_ranks
  names(taxon_ranks) = NCBI_ids[names(taxon_ranks)]









  taxon_ranks = taxon_ranks[names(taxon_ranks) %in% unique(unlist(taxon_sets))]
  taxon_sets <- lapply(taxon_sets, function(taxon_set) {
    return(unique(taxon_set[taxon_set %in% names(taxon_ranks)]))})
  set_sizes <- sapply(taxon_sets, function(taxon_set) { length(taxon_set) })
  taxon_sets = taxon_sets[set_sizes>5]
  # Select the top 5 taxon sets based on P value
  top_taxon_sets <- taxsea_results[order(taxsea_results$PValue), 
                                   "taxonSetName"][1:n_to_plot]

  # Create a data frame for ggplot2
  barcode_df <- data.frame(Rank = integer(), TaxonSet = character(), 
                           InSet = logical())

  for (taxon_set_name in top_taxon_sets) {
    taxon_set <- taxon_sets[[taxon_set_name]]
    temp_df <- data.frame(
      Rank = taxon_ranks,
      TaxonSet = taxon_set_name,
      InSet = names(taxon_ranks) %in% taxon_set
    )
    barcode_df <- rbind(barcode_df, temp_df)
  }

  # Generate the barcode plot
  barcode_plot <- ggplot(barcode_df, aes(x = Rank, y = TaxonSet, 
                                         col = InSet)) +
    geom_point(data = barcode_df[barcode_df$InSet,],size = 4,
               shape="|",show.legend = TRUE,color="black",
               fill="grey44",alpha=1) +
    labs(title = "Barcode plot",
         x = "log2 FC",
         y = "Taxon Set") +
    geom_vline(xintercept = 0,linetype=3)+
    xlim(axis_limits)+
    theme_bw() +
    theme(axis.title.y = element_blank(),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          panel.grid.major.y = element_line(color = "gray"))


  # Generate the barcode plot
  boxcode_plot <- ggplot(barcode_df, aes(y = Rank, x = TaxonSet, 
                                         col = InSet)) +
    geom_point(data = barcode_df[barcode_df$InSet,],size = 2,
               shape=21,show.legend = TRUE,color="black",fill="grey44",alpha=1) +
    geom_boxplot(data = barcode_df[barcode_df$InSet,],color="black",
                 alpha=0)+
    labs(title = "Boxcode plot",
         y = "log2 FC",
         x = "Taxon Set") +
    geom_hline(yintercept = 0,linetype=3,size=2)+
    ylim(-7,7)+
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold",
                                hjust = 0.5))+coord_flip()

  if(boxplot==FALSE)
  {return(barcode_plot)}
  else
  {return(boxcode_plot)}
}
