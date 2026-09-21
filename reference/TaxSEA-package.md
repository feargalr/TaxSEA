# TaxSEA: Taxon Set Enrichment Analysis

TaxSEA makes metagenomic data easier to interpret. Rather than reading
species one at a time, it tests whether a-priori defined sets of taxa
shift together, drawing on public databases (BugSigDB, MiMeDB, GutMGene,
mBodyMap, BacDive and GMRepoV2) and on sets collated from the
literature. Two approaches are provided. TaxSEA takes taxon names or
NCBI IDs with a rank per taxon, such as a fold change or correlation
coefficient, and tests each set with a Kolmogorov-Smirnov test. ssTaxSEA
scores every set in a single sample as the mean centered log-ratio of
its members, giving a per-sample value to use in any downstream model.
Both work with any taxonomic profiling technology, including 16S rRNA
gene sequencing, shotgun metagenomics and metatranscriptomics.

## See also

Useful links:

- <https://github.com/feargalr/taxsea>

- <https://feargalr.github.io/TaxSEA/>

- Report bugs at <https://github.com/feargalr/taxsea/issues>

## Author

**Maintainer**: Feargal Ryan <feargalr@gmail.com>
([ORCID](https://orcid.org/0000-0002-1565-4598)) (funding: Supported by
NHMRC Investigator Grant) \[funder\]

Authors:

- Feargal Ryan <feargalr@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-1565-4598)) (funding: Supported
  by NHMRC Investigator Grant) \[funder\]
