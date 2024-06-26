% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{TaxSEA_db}
\alias{TaxSEA_db}
\title{TaxSEA Database
A dataset containing taxon sets. Each item in the list is a taxon set,
and each member within a taxon set is a taxon.}
\format{
A list of vectors. Each vector contains character strings
representing taxa.
}
\source{
See READ ME.
}
\usage{
TaxSEA_db
}
\description{
TaxSEA Database
A dataset containing taxon sets. Each item in the list is a taxon set,
and each member within a taxon set is a taxon.
}
\examples{
data(TaxSEA_db)
all_sets = names(TaxSEA_db)
GABA_producers = TaxSEA_db[["MiMeDB_producers_of_GABA"]]
}
\keyword{datasets}
\value{A list of vectors. Each item is a taxon set}
