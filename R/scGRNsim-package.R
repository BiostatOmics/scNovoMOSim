# scNovoMOSim-package.R
# Package-level documentation and imports.

#' scNovoMOSim: Simulate Single-Cell Multiome Data with Ground-Truth GRNs
#'
#' @description
#' `scNovoMOSim` generates synthetic paired scRNA-seq and scATAC-seq count
#' matrices together with a fully specified ground-truth gene regulatory
#' network (GRN).  It is designed for benchmarking GRN inference algorithms
#' and differential network analysis tools applied to single-cell multiomics
#' data.
#'
#' ## Key features
#'
#' - **Two regulatory network layers**: TF → gene (GRN) and region ← TF
#'   (chromatin accessibility).
#' - **Hierarchical network structure**: shared backbone → cell-type-specific
#'   layers → condition-specific rewiring.
#' - **Controlled rewiring**: edge addition, deletion, sign-flip and
#'   weight drift, all recorded in a tidy ground-truth event table.
#' - **Realistic count models**: Negative Binomial for RNA, Poisson
#'   (optionally binarised) for ATAC, with cell-level correlated library
#'   sizes and logistic dropout.
#' - **Batch effects**: replicate-level shifts in expression, accessibility
#'   and sequencing depth.
#'
#' ## Typical workflow
#'
#' ```r
#' library(scNovoMOSim)
#'
#' cfg <- default_config()          # start from defaults
#' cfg$n_genes  <- 500              # customise as needed
#' cfg$seed     <- 2024
#'
#' sim <- simulate_dataset(cfg)     # run simulation
#' quick_qc(sim)                    # sanity check
#' export_simulation(sim, "out/", cfg)  # write to disk
#' ```
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Matrix sparseMatrix writeMM Matrix rowSums colSums
#' @importFrom jsonlite write_json
## usethis namespace: end
NULL
