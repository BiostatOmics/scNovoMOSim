# simulate.R
# Core simulation engine: generates paired scRNA-seq and scATAC-seq
# count matrices from the ground-truth regulatory networks.
#
#' Simulate paired scRNA-seq and scATAC-seq data with ground-truth GRNs
#'
#' This is the main entry point of the package.  It orchestrates all
#' simulation steps:
#'
#' 1. **Cell metadata** – assigns every cell to a condition, replicate and
#'    cell type via [make_metadata()].
#' 2. **Ground-truth networks** – builds cell-type- and condition-specific
#'    TF-gene (GRN) and region-TF weight matrices via [build_ground_truth()].
#' 3. **Biological baselines** – samples gene/region base expression levels,
#'    cell-type shifts, TF activity profiles and replicate batch effects.
#' 4. **Correlated sequencing depths** – RNA and ATAC library sizes per cell
#'    are drawn from a bivariate log-Normal to reflect the within-cell
#'    correlation observed in multiome experiments.
#' 5. **Count generation (chunked)** – for memory efficiency, cells are
#'    processed in chunks of size `cfg$chunk_size`.  Within each chunk:
#'    - Latent TF activity is sampled and used to drive ATAC accessibility
#'      (via `M`) and gene expression (via `W`).
#'    - ATAC counts are Poisson-sampled from the accessibility signal,
#'      optionally binarised.
#'    - RNA counts are Negative-Binomial-sampled from the transcriptional
#'      signal, which integrates GRN drive and a chromatin accessibility
#'      term.
#'    - Dropout is applied to both modalities independently.
#'
# RNA predictor:
#   grn_drive[i,j] = sum_k W[j,k] * softplus(A_regulation[i,k]) * pi[i, r(k,j)]
#   where r(k,j) is the region connecting TF k to gene j, and softplus(x) =
#   log(1+exp(x)) (see .softplus()). The gate pi[i, r(k,j)] is ~0 when the
#   peak is closed, silencing the regulatory contribution of TF k on gene j
#   in that cell. softplus() is applied only inside this gate computation;
#   sim$TFActivity still exports raw, untransformed A.
#
# ATAC predictor:
#   Independent of TF activity (no TF -> accessibility feedback yet).
#   pi is determined by region_base, ct_region_shift, replicate batch,
#   a per-cell noise term (atac_accessibility_noise_sd), and an additive
#   boost (atac_regulator_boost) applied to regions that mediate an
#   active TF-gene edge in that cell's (cell_type, condition) group.
#   ATAC counts are then generated in two steps, not a shared-budget
#   Poisson: n_open ~ Binomial(2, pi) (accessible copies, at most 2, a
#   diploid locus), then ATAC ~ Binomial(n_open, p_capture) where
#   p_capture is a per-cell sequencing-depth-dependent detection
#   probability (atac_capture_intercept/slope). This keeps counts
#   physically bounded and avoids the dilution a shared per-cell budget
#   would cause when many regions are simultaneously accessible.

#' Simulate paired scRNA-seq and scATAC-seq with ground-truth GRNs
#'
#' @param cfg A named list of simulation parameters.  Use [default_config()]
#'   to obtain a valid starting point.  Key fields:
#'   \describe{
#'     \item{`seed`}{Integer random seed for reproducibility.}
#'     \item{`n_genes`, `n_tfs`, `n_regions`}{Feature dimensions.}
#'     \item{`conditions`, `cell_types`, `n_cells_per_type`,
#'           `n_reps_per_condition`}{Experimental design.}
#'     \item{`affected_cell_types`}{Which cell types are rewired in the
#'           non-reference condition.}
#'     \item{`ct_gene_sigma`, `ct_region_sigma`}{SD of cell-type-specific
#'           shifts in gene and region baseline levels.}
#'     \item{`tf_activity_ct_sigma`}{SD of cell-type-level TF activity
#'           baseline.}
#'     \item{`disease_tf_activity_sigma`}{Additional TF activity shift in
#'           affected cell types under the Disease condition.}
#'     \item{`tf_activity_noise_sd`}{Per-cell TF activity noise SD.}
#'     \item{`atac_regulator_boost`, `atac_accessibility_noise_sd`}{See
#'           [default_config()]: additive accessibility boost for
#'           regions mediating an active edge in a given (cell_type,
#'           condition) group, and per-cell accessibility noise SD
#'           (the ATAC analogue of `tf_activity_noise_sd`).}
#'     \item{`debug_export_pi`}{See [default_config()]: when `TRUE`,
#'           returns the latent accessibility matrix as `sim$debug_pi`.}
#'     \item{`replicate_gene_sigma`, `replicate_region_sigma`}{Batch-effect
#'           SD for genes and regions.}
#'     \item{`replicate_tf_sigma`}{SD of the per-sample TF activity batch
#'           effect. Set to 0 (default) to disable. When > 0, each sample
#'           receives a random shift on every TF's activity.}
#'     \item{`tf_batch_propagates_to_rna`}{Logical. When TRUE (Version A),
#'           the TF batch effect is included in the activity used to drive
#'           RNA and ATAC, creating correlated batch effects between
#'           co-regulated genes. When FALSE (Version B), the batch is only
#'           an observation bias in the exported TFActivity matrix; RNA and
#'           ATAC are driven by batch-free activity.}
#'     \item{`replicate_depth_sigma`}{Sample-level sequencing depth
#'           variability (log scale).}
#'     \item{`replicate_dropout_sigma`}{Sample-level noise added to the RNA
#'           dropout logit and, per-cell, to the ATAC capture logit (see
#'           `atac_capture_intercept`/`atac_capture_slope`).}
#'     \item{`rna_depth_mean`, `rna_depth_logsd`}{Mean and log-SD of RNA
#'           library sizes (log-Normal).}
#'     \item{`atac_depth_mean`, `atac_depth_logsd`}{Mean and log-SD of ATAC
#'           library sizes.}
#'     \item{`depth_corr`}{Pearson correlation between log RNA and log ATAC
#'           library sizes within a cell.}
#'     \item{`rna_theta_mean`, `rna_theta_logsd`}{Mean and log-SD of the
#'           per-gene NB dispersion parameter (sampled on log scale).}
#'     \item{`grn_strength`}{Scalar multiplier for the GRN drive term in
#'           the RNA linear predictor.}
#'     \item{`atac_tf_strength`}{Scalar multiplier for the TF drive in the
#'           ATAC logit linear predictor.}
#'     \item{`chromatin_to_rna_strength`}{Scalar multiplier for the
#'           chromatin accessibility term in the RNA linear predictor.}
#'     \item{`rna_dropout_intercept`, `rna_dropout_slope`}{Parameters of
#'           the logistic dropout model for RNA.}
#'     \item{`atac_capture_intercept`, `atac_capture_slope`}{See
#'           [default_config()]: parameters of the per-cell capture
#'           probability that replaces the old ATAC dropout model
#'           (`atac_dropout_intercept`/`atac_dropout_slope`, now unused).}
#'     \item{`binarize_atac`}{Logical. Whether to binarise ATAC counts to
#'           0/1 after sampling.}
#'     \item{`chunk_size`}{Number of cells processed per iteration.}
#'   }
#' @param prior Optional `data.frame` of prior knowledge edges with columns
#'   `source` (TF), `target` (gene), and optionally `effect` (1 / -1 / NA)
#'   and `confidence` ("A" / "B" / "C" / "D").  When supplied, the TF-gene
#'   GRN is built by sampling from this table rather than synthetically.
#'   Region-TF networks are always generated synthetically.  Pass `NULL`
#'   (default) to use the original synthetic construction.
#'
#' @return A named list with:
#'   \describe{
#'     \item{`RNA`}{Sparse Matrix (cells x genes) of integer RNA counts.}
#'     \item{`ATAC`}{Sparse Matrix (cells x regions) of integer (or
#'           binary) ATAC counts.}
#'     \item{`metadata`}{`data.frame` with cell-level annotations.}
#'     \item{`genes`, `regions`, `tfs`}{Character vectors of feature names.}
#'     \item{`grn_edges`}{Tidy edge list of all TF-gene weights.}
#'     \item{`region_tf_edges`}{Tidy edge list of all region-TF weights.}
#'     \item{`region_gene_links`}{`data.frame` mapping each region to its
#'           target gene.}
#'     \item{`rewiring_events`}{`data.frame` of all ground-truth rewiring
#'           events (edge additions, deletions, flips).}
#'     \item{`TFActivity`}{Numeric matrix (cells x TFs) of observed TF
#'           activity. Includes the replicate batch effect when
#'           `replicate_tf_sigma > 0`, regardless of
#'           `tf_batch_propagates_to_rna`.}
#'     \item{`prior_universe`}{Only present in prior-informed mode.
#'           `data.frame` of all prior edges filtered to the simulated
#'           TF/gene universe, with columns `source`, `target`, `effect`,
#'           `tf_idx`, `gene_idx`.  Edges active in `grn_edges` are the
#'           ground-truth positives; the remaining rows are inactive prior
#'           edges that serve as hard-negative false-positive controls.}
#'     \item{`debug_pi`}{Only present when `cfg$debug_export_pi = TRUE`.
#'           Numeric matrix (cells x regions) of the latent accessibility
#'           `pi` that actually gates `grn_drive`, before the Poisson
#'           count / dropout / library-size steps applied to produce
#'           `ATAC`. Useful for diagnosing the TF-accessibility-RNA
#'           mechanism without single-peak measurement noise.}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cfg <- default_config()
#' sim <- simulate_dataset(cfg)
#' dim(sim$RNA)
#' dim(sim$ATAC)
#' head(sim$grn_edges)
#' }
#'
#' @seealso [default_config()], [export_simulation()], [quick_qc()]
simulate_dataset <- function(cfg, prior = NULL) {
  set.seed(cfg$seed)

  simulate_atac <- isTRUE(cfg$simulate_atac)

  meta <- make_metadata(cfg)
  gt   <- build_ground_truth(cfg, prior = prior)

  G  <- length(gt$genes)
  TF <- length(gt$tfs)
  P  <- if (simulate_atac) length(gt$regions) else 0L
  K  <- length(cfg$cell_types)
  D  <- length(cfg$conditions)
  S  <- length(unique(meta$sample))
  N  <- nrow(meta)

  # ------------------------------------------------------------------
  # Biological baselines
  # ------------------------------------------------------------------
  gene_base <- rnorm(G, -2.0, 0.7)

  ct_gene_shift <- matrix(rnorm(K * G, 0, cfg$ct_gene_sigma), K, G)

  if (simulate_atac) {
    region_base     <- rnorm(P, -4.0, 0.9)
    ct_region_shift <- matrix(rnorm(K * P, 0, cfg$ct_region_sigma), K, P)
  }

  tf_base_ct <- matrix(rnorm(K * TF, 0, cfg$tf_activity_ct_sigma), K, TF)

  cond_tf_shift <- array(0, dim = c(K, D, TF))
  for (k in seq_along(cfg$cell_types)) {
    ct <- cfg$cell_types[k]
    if (ct %in% cfg$affected_cell_types) {
      di <- match("Disease", cfg$conditions)
      cond_tf_shift[k, di, ] <- rnorm(TF, 0, cfg$disease_tf_activity_sigma)
    }
  }

  sample_gene_eff <- matrix(rnorm(S * G, 0, cfg$replicate_gene_sigma), S, G)
  if (simulate_atac)
    sample_region_eff <- matrix(rnorm(S * P, 0, cfg$replicate_region_sigma), S, P)

  # Per-sample TF batch effect (S x TF); zero when replicate_tf_sigma = 0
  sample_tf_eff <- matrix(rnorm(S * TF, 0, cfg$replicate_tf_sigma), S, TF)

  theta_gene <- exp(rnorm(G, log(cfg$rna_theta_mean), cfg$rna_theta_logsd))

  # ------------------------------------------------------------------
  # Library sizes
  # When simulate_atac = TRUE: bivariate log-Normal (RNA and ATAC correlated).
  # When simulate_atac = FALSE: univariate log-Normal for RNA only.
  # ------------------------------------------------------------------
  if (simulate_atac) {
    cov_depth <- matrix(c(
      cfg$rna_depth_logsd^2,
      cfg$depth_corr * cfg$rna_depth_logsd * cfg$atac_depth_logsd,
      cfg$depth_corr * cfg$rna_depth_logsd * cfg$atac_depth_logsd,
      cfg$atac_depth_logsd^2
    ), nrow = 2, byrow = TRUE)
    Z            <- matrix(rnorm(N * 2), N, 2)
    depth_latent <- Z %*% chol(cov_depth)
    depth_latent[, 1] <- depth_latent[, 1] + log(cfg$rna_depth_mean)
    depth_latent[, 2] <- depth_latent[, 2] + log(cfg$atac_depth_mean)
    sample_depth_shift <- matrix(rnorm(S * 2, 0, cfg$replicate_depth_sigma), S, 2)
    depth_latent       <- depth_latent + sample_depth_shift[meta$sample_idx, ]
    lib_rna  <- exp(depth_latent[, 1])
    lib_atac <- exp(depth_latent[, 2])
  } else {
    rna_depth_raw  <- rnorm(N, log(cfg$rna_depth_mean), cfg$rna_depth_logsd)
    depth_shift    <- rnorm(S, 0, cfg$replicate_depth_sigma)[meta$sample_idx]
    lib_rna        <- exp(rna_depth_raw + depth_shift)
  }

  # ------------------------------------------------------------------
  # Normalisation factor for grn_drive
  # With ATAC: based on the mean number of active edges per gene in tf_region_gene.
  # Without ATAC: based on grn_shared_edge_prob over the full G x TF matrix.
  # ------------------------------------------------------------------
  if (simulate_atac) {
    ref_key         <- paste(cfg$cell_types[1], "Control", sep = "|")
    ref_edges       <- gt$tf_region_gene[[ref_key]]
    mean_edges_gene <- if (!is.null(ref_edges) && nrow(ref_edges) > 0)
      nrow(ref_edges) / G else 1
    grn_norm <- sqrt(max(mean_edges_gene, 1))
  } else {
    grn_norm <- sqrt(max(TF * cfg$grn_shared_edge_prob, 1))
  }

  # ------------------------------------------------------------------
  # ATAC regulator-boost mask: for each (cell type, condition), which
  # regions mediate at least one active TF-gene edge in tf_region_gene.
  # Precomputed once (static across the chunked loop below); regions in
  # this mask get a higher baseline accessibility via
  # cfg$atac_regulator_boost. Because active edges differ between Control
  # and Disease (via rewiring), the mask differs between conditions too,
  # giving the ATAC ground truth a genuine condition-specific signal.
  # ------------------------------------------------------------------
  if (simulate_atac) {
    regulated_region_mask <- list()
    for (ct in cfg$cell_types) {
      for (condition in cfg$conditions) {
        key   <- paste(ct, condition, sep = "|")
        edges <- gt$tf_region_gene[[key]]
        mask  <- rep(FALSE, P)
        if (!is.null(edges) && nrow(edges) > 0L)
          mask[unique(edges$region_idx)] <- TRUE
        regulated_region_mask[[key]] <- mask
      }
    }
  }

  # ------------------------------------------------------------------
  # Chunked count generation
  # ------------------------------------------------------------------
  RNA_blocks  <- list()
  ATAC_blocks <- list()
  TF_blocks   <- list()
  PI_blocks   <- list()
  bidx <- tidx <- 1L
  chunk <- cfg$chunk_size

  for (start in seq(1L, N, by = chunk)) {
    end       <- min(start + chunk - 1L, N)
    ix_global <- start:end
    n         <- length(ix_global)

    ct_chunk   <- meta$cell_type[ix_global]
    cond_chunk <- meta$condition[ix_global]
    k_chunk    <- meta$cell_type_idx[ix_global]
    d_chunk    <- meta$condition_idx[ix_global]
    s_chunk    <- meta$sample_idx[ix_global]

    # ----------------------------------------------------------------
    # TF activity
    # ----------------------------------------------------------------
    tf_mu <- tf_base_ct[k_chunk, , drop = FALSE]
    for (ii in seq_len(n))
      tf_mu[ii, ] <- tf_mu[ii, ] + cond_tf_shift[k_chunk[ii], d_chunk[ii], ]

    A_base <- matrix(
      rnorm(n * TF, mean = as.vector(tf_mu), sd = cfg$tf_activity_noise_sd),
      nrow = n, ncol = TF
    )

    # Observed TF activity (includes per-sample batch effect)
    A <- A_base + sample_tf_eff[s_chunk, , drop = FALSE]

    # Activity used to drive RNA / ATAC:
    #   Version A (tf_batch_propagates_to_rna = TRUE):
    #     batch propagates into expression, creating correlated batch
    #     effects between co-regulated genes.
    #   Version B (tf_batch_propagates_to_rna = FALSE):
    #     expression uses batch-free A_base; the batch is only an
    #     observation bias visible in TFActivity.
    if (isTRUE(cfg$tf_batch_propagates_to_rna)) {
      A_regulation <- A
    } else {
      A_regulation <- A_base
    }

    TF_blocks[[tidx]] <- A   # always export observed (biased) activity
    tidx <- tidx + 1L

    # ----------------------------------------------------------------
    # ATAC and RNA drive: bifurcated on simulate_atac.
    # ----------------------------------------------------------------
    if (simulate_atac) {
      # ATAC: independent of TF activity, but regions that mediate an
      # active regulatory edge in this cell's (cell type, condition) group
      # get a higher baseline accessibility (cfg$atac_regulator_boost)
      # than non-regulating background regions, plus per-cell noise
      # (cfg$atac_accessibility_noise_sd) so the signal is not a
      # deterministic group-level label recoverable from pseudobulk alone.
      boost_mat <- matrix(0, n, P)
      for (ct in cfg$cell_types) {
        for (condition in cfg$conditions) {
          local_ix <- which(ct_chunk == ct & cond_chunk == condition)
          if (length(local_ix) == 0L) next
          mask <- regulated_region_mask[[paste(ct, condition, sep = "|")]]
          boost_mat[local_ix, ] <-
            matrix(mask, length(local_ix), P, byrow = TRUE) * cfg$atac_regulator_boost
        }
      }

      logit_pi <- matrix(region_base, n, P, byrow = TRUE) +
        ct_region_shift[k_chunk, , drop = FALSE] +
        sample_region_eff[s_chunk, , drop = FALSE] +
        boost_mat +
        matrix(rnorm(n * P, 0, cfg$atac_accessibility_noise_sd), n, P)

      pi       <- sigmoid(logit_pi)
      if (isTRUE(cfg$debug_export_pi))
        PI_blocks[[bidx]] <- pi

      # Biological accessibility state: each region has at most 2
      # accessible copies (diploid alleles), each independently open with
      # probability pi. This replaces the old shared-budget compositional
      # count model (rel_atac = pi / rowSums(pi)), which is appropriate
      # for RNA (relative expression competing for a fixed total) but not
      # for ATAC (a physically bounded, near-binary per-locus state):
      # under the old model, a truly open peak's expected count shrank
      # with however many OTHER regions happened to be active in that
      # cell's group, diluting detection regardless of that peak's own pi.
      n_open <- matrix(rbinom(n * P, size = 2L, prob = as.vector(pi)), n, P)

      # Capture/sequencing step: each open copy is independently captured
      # with a per-cell probability that increases with sequencing depth
      # (cfg$atac_capture_intercept/slope), plus per-cell noise
      # (cfg$replicate_dropout_sigma, same role it had on the old dropout
      # logit). This is where depth affects sensitivity, without letting
      # counts exceed the number of biologically open copies.
      capture_noise <- rnorm(n, 0, cfg$replicate_dropout_sigma)
      p_capture <- sigmoid(
        cfg$atac_capture_intercept +
          cfg$atac_capture_slope * log(lib_atac[ix_global]) +
          capture_noise
      )
      p_capture_mat <- matrix(p_capture, n, P)
      atac_counts   <- matrix(
        rbinom(n * P, size = as.vector(n_open), prob = as.vector(p_capture_mat)),
        n, P
      )
      if (isTRUE(cfg$binarize_atac))
        atac_counts <- 1L * (atac_counts > 0L)

      # RNA drive with gate: contribution = W[j,k] * softplus(A[i,k]) * pi[i, r(k,j)]
      # softplus() is applied only here, inside the gate -- sim$TFActivity
      # (TF_blocks, below) still exports raw A, untransformed. See .softplus().
      grn_drive <- matrix(0, n, G)
      for (ct in cfg$cell_types) {
        for (condition in cfg$conditions) {
          local_ix <- which(ct_chunk == ct & cond_chunk == condition)
          if (length(local_ix) == 0L) next
          edges <- gt$tf_region_gene[[paste(ct, condition, sep = "|")]]
          if (is.null(edges) || nrow(edges) == 0L) next
          for (e in seq_len(nrow(edges))) {
            k <- edges$tf_idx[e]; r <- edges$region_idx[e]
            j <- edges$gene_idx[e]; w <- edges$effect[e]
            grn_drive[local_ix, j] <- grn_drive[local_ix, j] +
              w * .softplus(A_regulation[local_ix, k]) * pi[local_ix, r]
          }
        }
      }

    } else {
      # RNA-only: ungated GRN drive = A * W (no peak accessibility).
      grn_drive <- matrix(0, n, G)
      for (ct in cfg$cell_types) {
        for (condition in cfg$conditions) {
          local_ix <- which(ct_chunk == ct & cond_chunk == condition)
          if (length(local_ix) == 0L) next
          W_ctd <- gt$W[[paste(ct, condition, sep = "|")]]
          grn_drive[local_ix, ] <- A_regulation[local_ix, , drop = FALSE] %*% t(W_ctd)
        }
      }
    }

    # ----------------------------------------------------------------
    # RNA counts (NB)
    # ----------------------------------------------------------------
    eta <- matrix(gene_base, n, G, byrow = TRUE) +
      ct_gene_shift[k_chunk, , drop = FALSE] +
      cfg$grn_strength * grn_drive / grn_norm +
      sample_gene_eff[s_chunk, , drop = FALSE]

    rel_rna <- exp(pmin(pmax(eta, -12), 12))
    rel_rna <- rel_rna / pmax(rowSums(rel_rna), 1e-12)
    mu_rna  <- sweep(rel_rna, 1, lib_rna[ix_global], `*`)

    rna_counts <- nb_counts(mu_rna, theta_gene)
    rna_drop_shift <- rnorm(n, 0, cfg$replicate_dropout_sigma)
    p_drop_rna <- sigmoid(
      cfg$rna_dropout_intercept -
        cfg$rna_dropout_slope * log1p(mu_rna) +
        rna_drop_shift
    )
    rna_counts[matrix(runif(length(rna_counts)), n, G) < p_drop_rna] <- 0L

    RNA_blocks[[bidx]] <- Matrix::Matrix(rna_counts, sparse = TRUE)
    if (simulate_atac)
      ATAC_blocks[[bidx]] <- Matrix::Matrix(atac_counts, sparse = TRUE)
    bidx <- bidx + 1L

    message("Simulated cells ", start, ":", end, " / ", N)
  }

  RNA <- do.call(rbind, RNA_blocks)
  rownames(RNA) <- meta$cell_id
  colnames(RNA) <- gt$genes

  # ------------------------------------------------------------------
  # Ground-truth edge lists
  # ------------------------------------------------------------------
  grn_edges <- list()
  eidx      <- 1L
  for (ct in cfg$cell_types) {
    for (condition in cfg$conditions) {
      grn_edges[[eidx]] <- matrix_to_grn_edges(
        gt$W[[paste(ct, condition, sep = "|")]],
        gt$genes, gt$tfs, condition, ct,
        prior_universe = if (gt$prior_mode) gt$prior_universe else NULL
      )
      eidx <- eidx + 1L
    }
  }
  grn_edges <- do.call(rbind, grn_edges)

  TF_activity           <- do.call(rbind, TF_blocks)
  rownames(TF_activity) <- meta$cell_id
  colnames(TF_activity) <- gt$tfs

  meta_out <- meta[, c("cell_id", "condition", "replicate", "sample", "cell_type")]

  out <- list(
    RNA             = RNA,
    metadata        = meta_out,
    genes           = gt$genes,
    tfs             = gt$tfs,
    grn_edges       = grn_edges,
    rewiring_events = gt$rewiring_events,
    TFActivity      = TF_activity
  )

  if (simulate_atac) {
    ATAC <- do.call(rbind, ATAC_blocks)
    rownames(ATAC) <- meta$cell_id
    colnames(ATAC) <- gt$regions

    region_tf_edges <- list(); ridx <- 1L
    for (ct in cfg$cell_types) {
      for (condition in cfg$conditions) {
        region_tf_edges[[ridx]] <- matrix_to_region_tf_edges(
          gt$M[[paste(ct, condition, sep = "|")]],
          gt$regions, gt$tfs, condition, ct
        )
        ridx <- ridx + 1L
      }
    }

    region_gene_links <- do.call(rbind, lapply(seq_len(P), function(r) {
      js <- which(as.matrix(gt$R2G)[r, ] > 0)
      if (length(js) == 0L) return(NULL)
      data.frame(region = gt$regions[r], gene = gt$genes[js],
                 stringsAsFactors = FALSE)
    }))

    # Convert tf_region_gene from integer indices to feature names for export
    tf_rg_named <- lapply(names(gt$tf_region_gene), function(key) {
      df <- gt$tf_region_gene[[key]]
      data.frame(
        tf        = gt$tfs[df$tf_idx],
        region    = gt$regions[df$region_idx],
        gene      = gt$genes[df$gene_idx],
        effect    = df$effect,
        stringsAsFactors = FALSE
      )
    })
    names(tf_rg_named) <- names(gt$tf_region_gene)

    # Convert B to tidy data frame: all TF-region binding sites per CT
    tf_region_binding <- lapply(names(gt$B), function(key) {
      B_mat <- gt$B[[key]]
      idx   <- which(B_mat > 0, arr.ind = TRUE)
      if (nrow(idx) == 0L) return(data.frame(
        ct_condition = character(0), tf = character(0),
        region = character(0), stringsAsFactors = FALSE))
      data.frame(
        ct_condition = key,
        tf           = gt$tfs[idx[, 2L]],
        region       = gt$regions[idx[, 1L]],
        stringsAsFactors = FALSE
      )
    })
    tf_region_binding <- do.call(rbind, tf_region_binding)

    out$ATAC              <- ATAC
    out$regions           <- gt$regions
    out$region_tf_edges   <- do.call(rbind, region_tf_edges)
    out$region_gene_links <- region_gene_links
    out$tf_region_gene    <- tf_rg_named
    out$tf_region_binding <- tf_region_binding

    if (isTRUE(cfg$debug_export_pi)) {
      debug_pi <- do.call(rbind, PI_blocks)
      rownames(debug_pi) <- meta$cell_id
      colnames(debug_pi) <- gt$regions
      out$debug_pi <- debug_pi
    }
  }

  if (gt$prior_mode)
    out$prior_universe <- gt$prior_universe

  out
}

#' Numerically stable softplus: log(1 + exp(x))
#'
#' Used only inside the RNA regulatory gate (see [simulate_dataset()]) to
#' transform TF activity before it multiplies the accessibility gate `pi`.
#' Unlike `A` itself (symmetric around 0 by design), `softplus(A)` is
#' always positive and its mean is always > 0 regardless of where `A` is
#' centred, which keeps the accessibility term `pi` from having an exactly
#' null marginal correlation with expression (see project discussion).
#' Stays close to the identity for A >> 0 (preserving TF resolution where
#' it matters most) and smoothly approaches 0 for A << 0, unlike a
#' sigmoid, which saturates and compresses resolution on both tails.
#'
#' @keywords internal
.softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))
