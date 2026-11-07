\name{Prognoser}
\alias{Prognoser}
\title{Integrative Prognostic Gene Signature Discovery}
\usage{
Prognoser(expression_accession_vector, denovo_threshold = 0.05,
  knowledge_genesets, na_ratio = 0.3, gene_threshold = 0.001,
  min_genes = 5, similarity_threshold = 0.5,
  linkage = "ward.D",
  de_redundant_basis = c("agglomeration", "igraph_component"))
}
\arguments{
  \item{expression_accession_vector}{A character vector of dataset names.
    Each dataset (e.g., \code{"GSE12345"}) must exist as a data frame in the
    global environment, with rows as genes and columns as samples. For each
    dataset, a corresponding survival data frame named
    \code{paste0(dataset, "_survival")} must also exist, containing columns
    \code{time} (numeric) and \code{status} (numeric, 0 = alive, 1 = dead).}
  \item{denovo_threshold}{Numeric. P-value threshold for univariate Cox
    regression to classify a gene as a "protection" (HR < 1) or "risk"
    (HR >= 1) gene within each dataset. Default is 0.05.}
  \item{knowledge_genesets}{A list of pre-defined gene sets (character
    vectors) to be incorporated into the candidate gene set pool.}
  \item{na_ratio}{Numeric. Maximum allowed proportion of missing values (NA)
    for a gene across datasets when building the effect-size matrix. Genes
    with NA proportion exceeding this value are removed. Default is 0.3.}
  \item{gene_threshold}{Numeric. P-value threshold (from robust rank
    aggregation) for selecting genes into the final gene pools. Default is
    0.001.}
  \item{min_genes}{Integer. Minimum number of genes required in a final
    clustered gene set. Sets with fewer genes are discarded. Default is 5.}
  \item{similarity_threshold}{Numeric. Jaccard similarity threshold for
    clustering gene sets. Gene sets with similarity >= this value are merged.
    Default is 0.5.}
  \item{linkage}{Character. Linkage method for the hierarchical-clustering
    step. One of \code{"ward.D"}, \code{"ward.D2"}, \code{"single"},
    \code{"complete"}, \code{"average"}, \code{"mcquitty"},
    \code{"median"} or \code{"centroid"}. Ignored when
    \code{de_redundant_basis = "igraph_component"}. Default is
    \code{"ward.D"}.}
  \item{de_redundant_basis}{Character. How biologically redundant gene sets
    are consolidated. \code{"agglomeration"} (default) cuts the
    hierarchical clustering of the Jaccard distances at
    \code{1 - similarity_threshold} and unions each cluster.
    \code{"igraph_component"} builds an undirected graph whose edges are
    the gene-set pairs with Jaccard similarity at or above
    \code{similarity_threshold}, then unions each connected component.
    With \code{linkage = "single"} the two partitions are identical, so
    \code{"igraph_component"} names the same consolidation explicitly,
    is faster on the sparse graph and returns the network itself.}
}
\value{
  A list with three elements, plus two more when
  \code{de_redundant_basis = "igraph_component"}:
  \item{genesets}{A list of final non-redundant gene signatures, named with
    \code{"LRGS"} (low-risk / protection) and \code{"HRGS"} (high-risk)
    prefixes.}
  \item{combined_matrix_protection}{A matrix of combined effect scores
    (sign * -log10(P)) for protection-direction genes, including a column
    \code{p_protection} with robust rank aggregation scores.}
  \item{combined_matrix_risk}{A matrix of combined effect scores for
    risk-direction genes, including a column \code{p_risk} with robust rank
    aggregation scores.}
  \item{network_protection}{The Jaccard edge list (\code{from},
    \code{to}, \code{jaccard}) consolidated into the protection
    signatures. Only present with
    \code{de_redundant_basis = "igraph_component"}.}
  \item{network_risk}{The same edge list for the risk signatures.}
}
\description{
  This function integrates multiple expression datasets and corresponding
  survival information to discover prognostic gene signatures. It performs
  univariate Cox regression in each dataset, meta-analyzes the results,
  screens genes based on effect direction and significance, and finally
  clusters gene sets to remove redundancy.
}
\examples{
\dontrun{
my_known_genesets <- list(immune = c("CD8A", "GZMB", "PRF1"))
result <- Prognoser(
  expression_accession_vector = c("GSE10000", "GSE20000"),
  knowledge_genesets = my_known_genesets,
  denovo_threshold = 0.05,
  gene_threshold = 0.001,
  similarity_threshold = 0.6
)
final_signatures <- result$genesets
}
}
