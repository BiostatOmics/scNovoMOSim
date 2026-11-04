# export.R
# Functions for writing simulation outputs to disk.

#' Export a simulation to a directory
#'
#' Writes all components of a simulation object to `outdir` in standard,
#' interoperable formats:
#'
#' | File | Format | Content |
#' |------|--------|---------|
#' | `RNA_counts_cells_x_genes.mtx` | Matrix Market | RNA count matrix |
#' | `ATAC_counts_cells_x_regions.mtx` | Matrix Market | ATAC count matrix |
#' | `cell_metadata.csv` | CSV | Per-cell annotations |
#' | `rna_features_genes.csv` | CSV | Gene names |
#' | `atac_features_regions.csv` | CSV | Region (peak) names |
#' | `tf_list.csv` | CSV | TF names |
#' | `ground_truth_GRN_TF_gene_edges.csv` | CSV | TF → gene edge list |
#' | `ground_truth_region_TF_edges.csv` | CSV | Region ← TF edge list |
#' | `ground_truth_region_gene_links.csv` | CSV | Region → gene assignment |
#' | `ground_truth_rewiring_events.csv` | CSV | Structural changes per condition |
#' | `simulation_config.json` | JSON | Full configuration used |
#'
#' @param sim    A simulation object as returned by [simulate_dataset()].
#' @param outdir Character scalar. Path to the output directory.  Created
#'               recursively if it does not exist.
#' @param cfg    The configuration list used to produce `sim`.  Saved as
#'               `simulation_config.json` for reproducibility.
#'
#' @return Invisibly returns `outdir`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cfg <- default_config()
#' sim <- simulate_dataset(cfg)
#' export_simulation(sim, "output/my_simulation", cfg)
#' }
#'
#' @seealso [simulate_dataset()]
export_simulation <- function(sim, outdir, cfg) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  Matrix::writeMM(sim$RNA,
                  file.path(outdir, "RNA_counts_cells_x_genes.mtx"))
  Matrix::writeMM(sim$ATAC,
                  file.path(outdir, "ATAC_counts_cells_x_regions.mtx"))

  write.csv(sim$metadata,
            file.path(outdir, "cell_metadata.csv"),
            row.names = FALSE)
  write.csv(data.frame(gene   = sim$genes),
            file.path(outdir, "rna_features_genes.csv"),
            row.names = FALSE)
  write.csv(data.frame(region = sim$regions),
            file.path(outdir, "atac_features_regions.csv"),
            row.names = FALSE)
  write.csv(data.frame(tf     = sim$tfs),
            file.path(outdir, "tf_list.csv"),
            row.names = FALSE)
  write.csv(sim$grn_edges,
            file.path(outdir, "ground_truth_GRN_TF_gene_edges.csv"),
            row.names = FALSE)
  write.csv(sim$region_tf_edges,
            file.path(outdir, "ground_truth_region_TF_edges.csv"),
            row.names = FALSE)
  write.csv(sim$region_gene_links,
            file.path(outdir, "ground_truth_region_gene_links.csv"),
            row.names = FALSE)
  write.csv(sim$rewiring_events,
            file.path(outdir, "ground_truth_rewiring_events.csv"),
            row.names = FALSE)

  jsonlite::write_json(cfg,
                       file.path(outdir, "simulation_config.json"),
                       pretty     = TRUE,
                       auto_unbox = TRUE)

  message("Simulation exported to: ", outdir)
  invisible(outdir)
}
