# Prognoser 0.2.1

Bug fixes found by auditing the pipeline end to end.

## Robustness

* `pmatrix` and `hrmatrix` were subset without `drop = FALSE`. With a single
  cohort group they collapsed to a vector and `rowSums(is.na(pmatrix))` aborted
  with "x must be an array of at least two dimensions". A one-dataset run now
  completes.
* Loops written as `for (i in 1:length(x))` iterate over `c(1, 0)` when `x` is
  empty, and `paste0("LRGS", 1:length(x))` builds two names for a zero-length
  vector. Every such site now uses `seq_along()` and `sprintf()`, so a run in
  which one direction or one gene pool is empty returns an empty set instead of
  aborting with "subscript out of bounds" or "'names' attribute [2] must be the
  same length as the vector [0]".
* `aux.merge()` returns the input unchanged when fewer than two gene sets
  survive, instead of handing a one-row distance matrix to `hclust()`.
* Rows whose meta-analysis failed (HR or p-value `NA`) are dropped before gene
  screening, so they can no longer leak into the protective / risk pools.

## Linkage

* `median` and `centroid` linkage can produce non-monotonic merge heights, which
  `cutree()` refuses. The previous workaround reordered `hc$merge`, which breaks
  the parent references of the tree: at `h >= 0.3` the returned labels could be
  structurally impossible (for example `1 2 2 3 1 4 5 3`). The heights are now
  made monotone with `cummax()` and the merge matrix is left untouched.

## Package hygiene

* the three `library()` calls were redundant — the packages are imported in
  `NAMESPACE` — and made `R CMD check` warn; removed.
* the `meta` non-standard-evaluation symbols are declared with
  `utils::globalVariables()`, clearing the "no visible binding for global
  variable" note.

`R CMD check` now reports `Status: OK`.

## Behaviour

Multi-cohort results are unchanged. The same synthetic four-cohort dataset,
run before and after the fixes, returns identical signatures for the default,
the empty-pool, the single-cohort and the centroid configurations.
