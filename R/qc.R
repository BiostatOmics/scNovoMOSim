# qc.R
# Post-simulation quality control summaries.

#' Print quick quality-control summaries for a simulation
#'
#' Computes and prints four summary tables intended for a rapid sanity check
#' after running [simulate_dataset()]:
#'
#' 1. **Library-size and feature-detection medians** – median RNA/ATAC
#'    library sizes and number of detected genes/regions, broken down by
#'    condition, replicate and cell type.
#' 2. **GRN edge counts** – number of TF-gene edges per condition and cell
#'    type.
#' 3. **Region-TF edge counts** – number of region-TF edges per condition
#'    and cell type.
#' 4. **Rewiring event table** – count of `add`, `delete` and `flip` events
#'    per network layer and cell type (only printed if rewiring occurred).
#'
#' @param sim A simulation object as returned by [simulate_dataset()].
#'
#' @return Invisibly returns `NULL`.  The function is called for its side
#'   effect of printing to the console.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cfg <- default_config()
#' sim <- simulate_dataset(cfg)
#' quick_qc(sim)
#' }
#'
#' @seealso [simulate_dataset()]
quick_qc <- function(sim) {
  meta <- sim$metadata
  RNA  <- sim$RNA
  ATAC <- sim$ATAC

  qc <- meta
  qc$rna_library_size      <- Matrix::rowSums(RNA)
  qc$atac_library_size     <- Matrix::rowSums(ATAC)
  qc$rna_detected_genes    <- Matrix::rowSums(RNA  > 0)
  qc$atac_detected_regions <- Matrix::rowSums(ATAC > 0)

  message("\n--- Median library sizes and detected features ---")
  print(
    aggregate(
      cbind(rna_library_size, atac_library_size,
            rna_detected_genes, atac_detected_regions) ~
        condition + replicate + cell_type,
      data = qc,
      FUN  = median
    )
  )

  message("\n--- TF-gene GRN edge counts (ground truth) ---")
  print(
    aggregate(effect ~ condition + cell_type,
              data = sim$grn_edges, FUN = length)
  )

  message("\n--- Region-TF edge counts (ground truth) ---")
  print(
    aggregate(effect ~ condition + cell_type,
              data = sim$region_tf_edges, FUN = length)
  )

  if (nrow(sim$rewiring_events) > 0) {
    message("\n--- Rewiring events by network, cell type and event type ---")
    print(
      as.data.frame(
        table(
          sim$rewiring_events$network,
          sim$rewiring_events$cell_type,
          sim$rewiring_events$event
        )
      )
    )
  }

  invisible(NULL)
}
