# metadata.R
# Functions for building the cell-level metadata table.

#' Build cell metadata for a simulation experiment
#'
#' Creates a `data.frame` in which every row represents one cell.  Cells are
#' organised hierarchically: conditions → replicates within each condition →
#' cell types within each replicate.  The number of cells per cell type is
#' taken from `cfg$n_cells_per_type`.
#'
#' In addition to the biological labels the function appends three integer
#' index columns (`condition_idx`, `cell_type_idx`, `sample_idx`) that map
#' each cell to the corresponding row/column of the parameter matrices built
#' later in the simulation pipeline.
#'
#' @param cfg A named list produced by [default_config()] or the user, with at
#'   least the following fields:
#'   \describe{
#'     \item{`conditions`}{Character vector of condition names
#'           (e.g. `c("Control", "Disease")`)`.}
#'     \item{`n_reps_per_condition`}{Integer. Number of biological replicates
#'           per condition.}
#'     \item{`cell_types`}{Character vector of cell-type names.}
#'     \item{`n_cells_per_type`}{Named integer vector; names must match
#'           `cell_types`.}
#'   }
#'
#' @return A `data.frame` with one row per cell and columns:
#'   \describe{
#'     \item{`cell_id`}{Unique zero-padded character identifier
#'           (`"cell_0000000"` …).}
#'     \item{`condition`}{Condition label.}
#'     \item{`replicate`}{Integer replicate index within the condition.}
#'     \item{`sample`}{Sample identifier formed as
#'           `"<condition>_R<replicate>"`.}
#'     \item{`cell_type`}{Cell-type label.}
#'     \item{`condition_idx`}{Integer position of the condition in
#'           `cfg$conditions`.}
#'     \item{`cell_type_idx`}{Integer position of the cell type in
#'           `cfg$cell_types`.}
#'     \item{`sample_idx`}{Integer position of the sample among all unique
#'           samples.}
#'   }
#'
#' @seealso [simulate_dataset()]
#' @keywords internal
make_metadata <- function(cfg) {
  rows <- list()
  idx  <- 1
  cid  <- 0   # running zero-based cell counter

  for (condition in cfg$conditions) {
    for (rep in seq_len(cfg$n_reps_per_condition)) {
      sample_id <- paste0(condition, "_R", rep)

      for (ct in cfg$cell_types) {
        n   <- cfg$n_cells_per_type[[ct]]
        ids <- cid + seq_len(n) - 1

        rows[[idx]] <- data.frame(
          cell_id   = sprintf("cell_%07d", ids),
          condition = condition,
          replicate = rep,
          sample    = sample_id,
          cell_type = ct,
          stringsAsFactors = FALSE
        )

        cid <- cid + n
        idx <- idx + 1
      }
    }
  }

  meta <- do.call(rbind, rows)

  # Integer index columns for fast matrix look-up later
  meta$condition_idx <- match(meta$condition, cfg$conditions)
  meta$cell_type_idx <- match(meta$cell_type,  cfg$cell_types)
  meta$sample_idx    <- match(meta$sample,      unique(meta$sample))

  meta
}
