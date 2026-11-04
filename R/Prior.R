# prior.R
# Functions for loading, validating and sampling from a prior knowledge
# table of TF-gene regulatory interactions.
#
# Expected input format (data.frame):
#   source  : character — TF name
#   target  : character — gene name
#   effect  : numeric   — direction of regulation: 1 (activation),
#                         -1 (repression), or NA (unknown). OPTIONAL column.

# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------

#' Validate and normalise a prior knowledge table
#'
#' Checks that the table has the required columns, coerces types, enforces
#' that `effect` (if present) only contains values in `{-1, 1, NA}`, and
#' removes duplicate `source`–`target` pairs.
#'
#' @param prior A `data.frame` with at least columns `source` and `target`.
#'   An optional `effect` column may contain `1`, `-1`, or `NA`.
#'
#' @return A cleaned `data.frame` with columns `source` (character),
#'   `target` (character), and `effect` (numeric; all `NA` if the column
#'   was absent in the input).
#'
#' @export
#'
#' @examples
#' db <- data.frame(source = c("TF1","TF2"), target = c("gene_1","gene_2"),
#'                  effect = c(1, -1))
#' validate_prior(db)
validate_prior <- function(prior) {
  if (!is.data.frame(prior)) {
    stop("`prior` must be a data.frame.")
  }

  required_cols <- c("source", "target")
  missing_cols  <- setdiff(required_cols, names(prior))
  if (length(missing_cols) > 0) {
    stop("Prior table is missing required column(s): ",
         paste(missing_cols, collapse = ", "), ".")
  }

  prior$source <- as.character(prior$source)
  prior$target <- as.character(prior$target)

  # Add effect column if absent (all NA = direction unknown)
  if (!"effect" %in% names(prior)) {
    message("Prior table has no 'effect' column. ",
            "Regulation direction will be sampled randomly for all edges.")
    prior$effect <- NA_real_
  } else {
    prior$effect <- as.numeric(prior$effect)
    invalid <- !is.na(prior$effect) & !prior$effect %in% c(-1, 1)
    if (any(invalid)) {
      warning(sum(invalid), " row(s) with 'effect' not in {-1, 1, NA} ",
              "will have their effect set to NA.")
      prior$effect[invalid] <- NA_real_
    }
  }

  # Resolve duplicate source-target pairs.
  # A duplicate is any pair (source, target) that appears more than once.
  # Two cases:
  #   - Consistent duplicates (same effect): collapse to one row.
  #   - Contradictory duplicates (conflicting effects): set effect = NA,
  #     meaning direction is unknown and will be sampled randomly later.
  pair_key <- paste(prior$source, prior$target, sep = "|")

  if (anyDuplicated(pair_key)) {
    unique_keys  <- unique(pair_key)
    resolved     <- vector("list", length(unique_keys))

    n_conflicted <- 0L

    for (k in seq_along(unique_keys)) {
      rows <- prior[pair_key == unique_keys[k], , drop = FALSE]

      if (nrow(rows) == 1L) {
        resolved[[k]] <- rows
      } else {
        unique_effects <- unique(rows$effect[!is.na(rows$effect)])

        if (length(unique_effects) > 1L) {
          # Contradictory evidence: more than one distinct non-NA effect
          rows           <- rows[1L, , drop = FALSE]
          rows$effect    <- NA_real_
          n_conflicted   <- n_conflicted + 1L
        } else {
          # Consistent (or all NA): keep one row with the agreed effect
          rows <- rows[1L, , drop = FALSE]
        }
        resolved[[k]] <- rows
      }
    }

    prior <- do.call(rbind, resolved)

    n_removed <- length(pair_key) - nrow(prior)
    message(n_removed, " duplicate source-target pair(s) collapsed; ",
            n_conflicted, " pair(s) had contradictory effects and were ",
            "set to effect = NA (direction will be sampled randomly).")
  }

  prior[, c("source", "target", "effect")]
}

#' Filter a prior knowledge table to the simulated gene/TF universe
#'
#' Retains only rows where `source` is a known TF **and** `target` is a
#' known gene in the simulation.  Rows involving unknown entities are
#' silently dropped, and a summary message is printed.
#'
#' The returned table also contains integer index columns (`tf_idx`,
#' `gene_idx`) for fast matrix look-up.
#'
#' @param prior A validated `data.frame` as returned by [validate_prior()].
#' @param genes Character vector of all simulated gene names.
#' @param tfs   Character vector of TF names (subset of `genes`).
#'
#' @return A `data.frame` with columns `source`, `target`, `effect`,
#'   `tf_idx` (integer index into `tfs`) and `gene_idx` (integer index
#'   into `genes`).  Returns an empty `data.frame` (0 rows) with those
#'   columns if no prior edges match the simulation universe.
#'
#' @export
filter_prior_to_universe <- function(prior, genes, tfs) {
  n_in <- nrow(prior)

  keep <- prior$source %in% tfs & prior$target %in% genes
  if (!any(keep)) {
    warning("No prior edges match the simulated TF/gene universe. ",
            "Falling back to synthetic GRN construction.")
    return(
      data.frame(source   = character(0), target  = character(0),
                 effect   = numeric(0),
                 tf_idx   = integer(0),   gene_idx = integer(0),
                 stringsAsFactors = FALSE)
    )
  }

  out <- prior[keep, , drop = FALSE]
  out$tf_idx   <- match(out$source, tfs)
  out$gene_idx <- match(out$target, genes)

  message(nrow(out), " / ", n_in,
          " prior edges retained after filtering to simulation universe.")
  out
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Build a GRN weight matrix from a set of active prior edges
#'
#' Populates a `genes × TFs` matrix using the active edges sampled from the
#' prior.  Edge weights are sampled from `N(0, weight_sd)` with their sign
#' constrained to the `effect` column when that value is known (`1` or `-1`).
#' For edges with `effect = NA` the sign is drawn uniformly at random.
#'
#' @param active_edges `data.frame` with columns `gene_idx`, `tf_idx`,
#'   `effect` (as returned by [sample_active_edges()]`$active`).
#' @param n_genes  Integer. Number of genes (rows of the output matrix).
#' @param n_tfs    Integer. Number of TFs (columns of the output matrix).
#' @param weight_sd Numeric. SD of `N(0, sd)` for edge magnitudes.
#'
#' @return A numeric matrix of dimensions `n_genes × n_tfs`.
#'
#' @keywords internal
build_W_from_prior_edges <- function(active_edges, n_genes, n_tfs,
                                    weight_sd, weight_mean = 0) {
  # In prior-informed mode the sign is always determined by the prior
  # (or sampled randomly for edges with unknown direction).
  # weight_mean controls the mean magnitude: when > 0 edges have a
  # guaranteed non-trivial effect size; when 0 magnitudes come from
  # |N(0, weight_sd)|, which also avoids exact zeros.
  W <- matrix(0, n_genes, n_tfs)
  if (nrow(active_edges) == 0) return(W)

  magnitudes <- abs(rnorm(nrow(active_edges), weight_mean, weight_sd))
  magnitudes <- pmax(magnitudes, 0.05)   # ensure a non-negligible effect size

  signs <- active_edges$effect
  unknown_sign <- is.na(signs)
  signs[unknown_sign] <- sample(c(-1, 1), sum(unknown_sign), replace = TRUE)

  for (k in seq_len(nrow(active_edges))) {
    i <- active_edges$gene_idx[k]
    j <- active_edges$tf_idx[k]
    W[i, j] <- signs[k] * magnitudes[k]
  }
  W
}