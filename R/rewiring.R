# rewiring.R
# Functions that introduce condition-specific structural changes
# (rewiring) into weight matrices representing GRN or region-TF links.

#' Format a single rewiring event as a one-row data frame
#'
#' Helper called inside [rewire_matrix()] to record every structural change
#' applied to a weight matrix.
#'
#' @param event      Character string: one of `"add"`, `"delete"`, `"flip"`,
#'   `"drift"`.
#' @param i          Integer row index of the modified entry.
#' @param j          Integer column index of the modified entry.
#' @param old        Numeric. Weight value before the change.
#' @param new        Numeric. Weight value after the change.
#' @param row_names  Character vector of row names for the matrix.
#' @param col_names  Character vector of column names for the matrix.
#' @param direction  One of `"col_to_row"` (TF -> gene) or `"row_to_col"`
#'   (region <- TF).
#' @param row_type   Character label for the row entity type.
#' @param col_type   Character label for the column entity type.
#'
#' @return A one-row `data.frame` with columns `event`, `source`, `target`,
#'   `source_type`, `target_type`, `old_effect`, `new_effect`.
#'
#' @keywords internal
make_rewire_event <- function(event, i, j, old, new,
                              row_names, col_names,
                              direction = c("col_to_row", "row_to_col"),
                              row_type, col_type) {
  direction <- match.arg(direction)

  if (direction == "col_to_row") {
    source      <- col_names[j];  target      <- row_names[i]
    source_type <- col_type;      target_type <- row_type
  } else {
    source      <- row_names[i];  target      <- col_names[j]
    source_type <- row_type;      target_type <- col_type
  }

  data.frame(
    event       = event,
    source      = source,
    target      = target,
    source_type = source_type,
    target_type = target_type,
    old_effect  = old,
    new_effect  = new,
    stringsAsFactors = FALSE
  )
}

#' Apply stochastic rewiring to a regulatory weight matrix
#'
#' Applies four types of perturbation in strict sequential order, where each
#' step only affects edges that have not been modified by any previous step:
#'
#' 1. **Deletion** - a fraction `q_del` of existing edges is set to zero.
#' 2. **Addition** - new edges proportional to `q_add` are introduced.
#' 3. **Sign flip** - a fraction `q_flip` of edges that were **not** deleted
#'    or added have their sign reversed.  When `prior_effects` is supplied,
#'    only edges with `NA` documented direction are eligible.
#' 4. **Weight drift** - a fraction `q_drift` of edges that were **not**
#'    deleted, added, or flipped are multiplied by a log-normal factor.
#'
#' This ordering ensures that each edge undergoes at most one type of change,
#' making `rewiring_events` a clean, non-redundant record of Disease differences.
#'
#' @param W             Numeric matrix (rows = target, cols = source).
#' @param row_names     Character vector of row names.
#' @param col_names     Character vector of column names.
#' @param row_type      Character label for row entity type (e.g. `"gene"`).
#' @param col_type      Character label for col entity type (e.g. `"tf"`).
#' @param direction     `"col_to_row"` or `"row_to_col"`.
#' @param q_del         Fraction of existing edges to delete. Default `0.15`.
#' @param q_add         Fraction of existing edges to add. Default `0.10`.
#' @param q_flip        Fraction of unmodified edges to sign-flip. Default
#'   `0.05`. Ignored when `allow_flip = FALSE`.
#' @param q_drift       Fraction of unmodified edges that undergo weight drift.
#'   Default `0.30`.
#' @param scale_sd      SD of the log-normal weight drift. Default `0.25`.
#' @param new_weight_sd SD of `N(0, sd)` for newly added edge weights.
#'   Default `0.8`.
#' @param allow_flip    Logical. Whether sign-flipping is allowed at all.
#'   Default `TRUE`.
#' @param add_pool      Optional two-column integer matrix of candidate
#'   positions for edge addition (prior-informed mode). `NULL` means any
#'   zero entry is a candidate (synthetic mode). An empty matrix (0 rows)
#'   means no additions are allowed (prior pool exhausted).
#' @param add_pool_effects Optional numeric vector (length = `nrow(add_pool)`)
#'   with documented effect directions for each candidate.
#' @param prior_effects Optional numeric matrix (same dimensions as `W`)
#'   with documented effect directions (`1`, `-1`, `NA`) for every position
#'   in the full prior universe. Edges with a non-NA documented direction
#'   are never flipped. Should be built from the complete prior universe
#'   (including inactive edges) so that newly added edges with known
#'   directions are also protected.
#'
#' @return A named list:
#'   \describe{
#'     \item{`matrix`}{The rewired weight matrix (same dimensions as `W`).}
#'     \item{`events`}{`data.frame` of all changes, or empty `data.frame`.}
#'   }
#'
#' @seealso [build_ground_truth()], [make_rewire_event()]
#' @keywords internal
rewire_matrix <- function(W,
                          row_names,
                          col_names,
                          row_type,
                          col_type,
                          direction,
                          q_del            = 0.15,
                          q_add            = 0.10,
                          q_flip           = 0.05,
                          q_drift          = 0.30,
                          scale_sd         = 0.25,
                          new_weight_sd    = 0.8,
                          allow_flip       = TRUE,
                          add_pool         = NULL,
                          add_pool_effects = NULL,
                          prior_effects    = NULL,
                          cfg              = NULL) {
  W2     <- W
  events <- list()
  eidx   <- 1L

  existing <- which(abs(W2) > 0, arr.ind = TRUE)
  n_edges  <- nrow(existing)

  # Track all modified positions to enforce "one change per edge" rule.
  # Keys are "i j" strings for fast membership testing.
  modified_keys <- character(0)

  # ------------------------------------------------------------------
  # 1. Edge deletion
  # ------------------------------------------------------------------
  n_del <- floor(q_del * n_edges)
  if (n_del > 0) {
    sel <- sample(seq_len(n_edges), n_del, replace = FALSE)
    for (s in sel) {
      i <- existing[s, 1]; j <- existing[s, 2]
      old      <- W2[i, j]
      W2[i, j] <- 0
      modified_keys <- c(modified_keys, paste(i, j))
      events[[eidx]] <- make_rewire_event(
        "delete", i, j, old, 0,
        row_names, col_names, direction = direction,
        row_type = row_type, col_type = col_type
      )
      eidx <- eidx + 1L
    }
  }

  # ------------------------------------------------------------------
  # 2. Edge addition
  # Three cases:
  #   add_pool = NULL       → synthetic: any zero entry is a candidate
  #   add_pool has rows     → prior mode: restrict to those positions
  #   add_pool is empty     → prior mode, pool exhausted: no additions
  # ------------------------------------------------------------------
  if (is.null(add_pool)) {
    candidates   <- which(abs(W2) == 0, arr.ind = TRUE)
    cand_effects <- NULL
  } else if (nrow(add_pool) > 0) {
    still_zero   <- W2[cbind(add_pool[, 1], add_pool[, 2])] == 0
    candidates   <- add_pool[still_zero, , drop = FALSE]
    cand_effects <- if (!is.null(add_pool_effects)) add_pool_effects[still_zero] else NULL
  } else {
    candidates   <- matrix(integer(0), nrow = 0, ncol = 2)
    cand_effects <- NULL
  }

  n_add <- floor(q_add * max(n_edges, 1L))

  if (n_add > 0 && nrow(candidates) > 0) {
    sel <- sample(seq_len(nrow(candidates)),
                  min(n_add, nrow(candidates)), replace = FALSE)
    for (s in sel) {
      i <- candidates[s, 1]; j <- candidates[s, 2]

      doc_effect <- if (!is.null(cand_effects)) cand_effects[s] else NA_real_
      # if (!is.na(doc_effect)) {
      #   new <- doc_effect * abs(rnorm(1, 0, new_weight_sd))
      # } else {
      #   new <- rnorm(1, 0, new_weight_sd)
      # }
      if (!is.na(doc_effect)) {
        new <- doc_effect * abs(rnorm(1, 0, new_weight_sd))
      } else {
        new <- rnorm(1, 0, new_weight_sd)
      }

      # DESPUÉS:
      if (!is.null(cfg)) {
        raw <- .sample_grn_weights(1L, cfg)
      } else {
        raw <- rnorm(1L, 0, new_weight_sd)
      }
      if (!is.na(doc_effect)) {
        new <- doc_effect * abs(raw)
      } else {
        new <- raw
      }

      if (abs(new) < 0.05) new <- sign(new + 1e-6) * 0.05

      W2[i, j] <- new
      modified_keys <- c(modified_keys, paste(i, j))
      events[[eidx]] <- make_rewire_event(
        "add", i, j, 0, new,
        row_names, col_names, direction = direction,
        row_type = row_type, col_type = col_type
      )
      eidx <- eidx + 1L
    }
  }

  # ------------------------------------------------------------------
  # 3. Sign flip
  # Only edges that were NOT deleted or added are eligible.
  # When prior_effects is supplied, edges with a known documented
  # direction (1 or -1) are additionally excluded — this covers both
  # active Control edges and add_pool edges (since prior_effects is
  # built from the full prior_universe).
  # ------------------------------------------------------------------
  if (allow_flip && q_flip > 0) {
    all_active <- which(abs(W2) > 0, arr.ind = TRUE)
    all_active_keys <- paste(all_active[, 1], all_active[, 2])

    # Remove edges modified in steps 1-2
    unmodified_mask <- !all_active_keys %in% modified_keys
    flippable       <- all_active[unmodified_mask, , drop = FALSE]

    # Additionally restrict by documented direction
    if (!is.null(prior_effects) && nrow(flippable) > 0) {
      doc_dir   <- prior_effects[cbind(flippable[, 1], flippable[, 2])]
      flippable <- flippable[is.na(doc_dir), , drop = FALSE]
    }

    n_flip <- floor(q_flip * nrow(all_active))
    n_flip <- min(n_flip, nrow(flippable))

    if (n_flip > 0) {
      sel <- sample(seq_len(nrow(flippable)), n_flip, replace = FALSE)
      for (s in sel) {
        i <- flippable[s, 1]; j <- flippable[s, 2]
        old      <- W2[i, j]
        W2[i, j] <- -old
        modified_keys <- c(modified_keys, paste(i, j))
        events[[eidx]] <- make_rewire_event(
          "flip", i, j, old, -old,
          row_names, col_names, direction = direction,
          row_type = row_type, col_type = col_type
        )
        eidx <- eidx + 1L
      }
    }
  }

  # ------------------------------------------------------------------
  # 4. Weight drift
  # Only edges that were NOT deleted, added, or flipped are eligible.
  # Models stable interactions that change only in magnitude.
  # ------------------------------------------------------------------
  if (q_drift > 0) {
    all_active  <- which(abs(W2) > 0, arr.ind = TRUE)
    active_keys <- paste(all_active[, 1], all_active[, 2])

    driftable <- all_active[!active_keys %in% modified_keys, , drop = FALSE]
    n_drift   <- floor(q_drift * nrow(driftable))

    if (n_drift > 0) {
      sel_drift     <- sample(seq_len(nrow(driftable)), n_drift, replace = FALSE)
      drift_idx     <- driftable[sel_drift, , drop = FALSE]
      drift_factors <- exp(rnorm(n_drift, 0, scale_sd))
      old_weights   <- W2[drift_idx]
      W2[drift_idx] <- old_weights * drift_factors

      for (k in seq_len(n_drift)) {
        i <- drift_idx[k, 1]; j <- drift_idx[k, 2]
        events[[eidx]] <- make_rewire_event(
          "drift", i, j, old_weights[k], W2[i, j],
          row_names, col_names, direction = direction,
          row_type = row_type, col_type = col_type
        )
        eidx <- eidx + 1L
      }
    }
  }

  events_df <- if (length(events) > 0) do.call(rbind, events) else data.frame()
  list(matrix = W2, events = events_df)
}
