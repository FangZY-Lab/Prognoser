# Prognoser

`Prognoser` is an R package for integrative prognostic gene-signature discovery.
It runs univariate Cox regression in each dataset, screens protective
(HR < 1) and risk (HR >= 1) genes across multiple datasets, and clusters gene
sets to remove redundancy, returning non-redundant gene signatures.

## Dependencies

Make sure the following R packages are installed:

- `survival`
- `RobustRankAggreg`
- `meta`

```r
install.packages(c("survival", "RobustRankAggreg", "meta"))
```

## Installation

### Option 1: install from the source tarball

```r
install.packages("~/Desktop/Prognoser_0.1.0.tar.gz", repos = NULL, type = "source")
```

### Option 2: install from GitHub

```r
install.packages("remotes")
remotes::install_github("FangZY-Lab/Prognoser")
```

Then load the package:

```r
library(Prognoser)
```

## Data format

`Prognoser()` reads data from the global environment via `get()`. It requires:

- Each dataset is a `data.frame` with genes as rows and samples as columns.
- Each dataset has a matching survival data frame named `<dataset>_survival`
  with numeric columns `time` and `status` (`0` = alive, `1` = dead). Its row
  names must match the column names of the expression data frame.

## Simulated data

The code below creates two simulated datasets (`SimDataA`, `SimDataB`) and
their survival data. For demonstration purposes, we embed 10 "risk" genes
(`GENE1`-`GENE10`, higher expression -> worse survival) and 10 "protective"
genes (`GENE11`-`GENE20`, higher expression -> better survival).

```r
set.seed(123)

make_dataset <- function(n_samples = 80, n_genes = 100, seed = 1) {
  set.seed(seed)

  # Latent risk score: higher score -> worse survival
  risk_score <- rnorm(n_samples)

  # Survival data: higher risk -> shorter survival time
  time <- rweibull(n_samples, shape = 1, scale = exp(2 - risk_score))
  status <- rbinom(n_samples, 1, 0.7)
  surv <- data.frame(
    time = round(time, 3),
    status = status,
    row.names = paste0("S", seq_len(n_samples))
  )

  # Expression matrix: genes as rows, samples as columns
  mat <- matrix(rnorm(n_samples * n_genes), nrow = n_genes, ncol = n_samples)

  # First 10 genes are risk genes (positively correlated with risk_score)
  mat[1:10, ] <- mat[1:10, ] + matrix(rep(risk_score * 2, each = 10), nrow = 10)

  # Genes 11-20 are protective genes (negatively correlated with risk_score)
  mat[11:20, ] <- mat[11:20, ] - matrix(rep(risk_score * 2, each = 10), nrow = 10)

  expr <- as.data.frame(mat)
  rownames(expr) <- paste0("GENE", seq_len(n_genes))
  colnames(expr) <- paste0("S", seq_len(n_samples))

  list(expr = expr, surv = surv)
}

dA <- make_dataset(seed = 11)
dB <- make_dataset(seed = 22)

# Put data into the global environment (the function reads them via get())
SimDataA <- dA$expr
SimDataA_survival <- dA$surv
SimDataB <- dB$expr
SimDataB_survival <- dB$surv
```

## Run

```r
# Optional prior / knowledge gene sets (added to the candidate gene-set pool)
knowledge <- list(immune = c("GENE1", "GENE2", "GENE3", "GENE11", "GENE12"))

result <- Prognoser(
  expression_accession_vector = c("SimDataA", "SimDataB"),
  knowledge_genesets = knowledge,
  denovo_threshold = 0.05,
  gene_threshold = 0.01,
  min_genes = 3,
  similarity_threshold = 0.5,
  linkage = "ward.D"
)

# Final non-redundant gene sets
str(result$genesets)

# Combined effect matrices for the protective / risk directions
head(result$combined_matrix_protection)
head(result$combined_matrix_risk)
```

The example uses `gene_threshold = 0.01` and `min_genes = 3` to keep the demo
output stable; adjust them as needed. Defaults are listed below.

## Arguments

| Argument | Description | Default |
| --- | --- | --- |
| `expression_accession_vector` | Dataset names; each must exist in the global environment together with `<name>_survival` | none |
| `denovo_threshold` | P-value threshold for univariate Cox screening (protective / risk genes) | `0.05` |
| `knowledge_genesets` | A `list` of prior gene sets (character vectors) | none |
| `na_ratio` | Maximum allowed proportion of missing values per gene across datasets | `0.3` |
| `gene_threshold` | P-value threshold from robust rank aggregation (RRA) | `0.001` |
| `min_genes` | Minimum number of genes in a final clustered gene set | `5` |
| `similarity_threshold` | Jaccard similarity threshold for clustering gene sets | `0.5` |
| `linkage` | Hierarchical clustering method; pass a single value (e.g. `"ward.D"`) | see note below |

## Return value

- `genesets`: final non-redundant gene sets. Low-risk (protective) sets are
  prefixed with `LRGS`, high-risk sets with `HRGS`.
- `combined_matrix_protection`: combined effect matrix (sign x `-log10(P)`) for
  protective-direction genes, including a `p_protection` column.
- `combined_matrix_risk`: combined effect matrix for risk-direction genes,
  including a `p_risk` column.

## Notes

1. **Pass a single value for `linkage`** (e.g. `"ward.D"`, `"average"`, or
   `"complete"`). The function's default is a vector of choices; omitting it
   will cause `hclust()` to error.
2. Data must be in the global environment and named exactly `<dataset>` and
   `<dataset>_survival`.
