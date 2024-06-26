% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{TaxSEA_test_data}
\alias{TaxSEA_test_data}
\title{TaxSEA Test Data}
\format{
A data frame with two columns:
\describe{
  \item{rank}{Character vector representing taxon ranks}
  \item{id}{Character vector representing taxon IDs}
}
}
\source{
See READ ME.
}
\usage{
TaxSEA_test_data
}
\description{
A dataset containing taxon ranks and taxon IDs.
}
\examples{
data(TaxSEA_test_data)
test_results <- TaxSEA(TaxSEA_test_data)
}
\keyword{datasets}
\value{A test dataset}
