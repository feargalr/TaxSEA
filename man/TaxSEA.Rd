% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TaxSEA.R
\name{TaxSEA}
\alias{TaxSEA}
\title{TaxSEA: Taxon Set Enrichment Analysis}
\usage{
TaxSEA(
  taxon_ranks,
  lookup_missing = FALSE,
  min_set_size = 5,
  max_set_size = 100,
  custom_db = NULL
)
}
\arguments{
\item{taxon_ranks}{A named vector of log2 fold changes between control and test groups.}

\item{lookup_missing}{Logical indicating whether to fetch missing NCBI IDs. Default is FALSE.}

\item{min_set_size}{Minimum size of taxon sets to include in the analysis. Default is 5.}

\item{max_set_size}{Maximum size of taxon sets to include in the analysis. Default is 100.}

\item{custom_db}{A user-provided list of taxon sets. If NULL (default), the built-in database is used.}
}
\value{
A list of data frames with taxon set enrichment results.
}
\description{
TaxSEA enables rapid annotation of changes observed in a microbiome association study
by testing for enrichment of producers of particular metabolites or associations with marker taxa
for specific diseases. It focuses on human gut microbiome associations and uses
a Kolmogorov-Smirnov test to determine if a particular set of taxa is enriched or depleted
relative to a control group.
}
\examples{
data("TaxSEA_test_data")
taxsea_results <- TaxSEA(TaxSEA_test_data)

}
\seealso{
\itemize{
  \item \url{https://doi.org/10.1093/nar/gkac868} for MiMeDB
  \item \url{https://doi.org/10.1093/nar/gkab1019} for GMrepo
  \item \url{https://doi.org/10.1093/nar/gkab786} for gutMGene
  \item \url{https://doi.org/10.1038/s41587-023-01872-y} for BugSigDB
}
}
