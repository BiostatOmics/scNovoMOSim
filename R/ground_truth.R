# ground_truth.R
# Construction of ground-truth regulatory networks (GRN, TF-region binding
# and region-gene links) and helpers to serialise weight matrices as tidy
# edge lists.
#
# The core data structure is tf_region_gene: a tidy table where each row is
# one active TF-region-gene regulatory triplet, with columns tf_idx,
# region_idx, gene_idx and weight. The RNA predictor uses the accessibility
# pi[i, r(k,j)] of the specific region connecting TF k to gene j as a
# multiplicative gate on the TF contribution.
#
# build_ground_truth() supports two modes:
#
#   Synthetic mode  (prior = NULL, default)
#     Gene and TF names are generated synthetically ("gene_1", "gene_2", ...).
#     R2G is built uniformly: each gene receives one promoter region and a
#     random number of enhancers; a fraction p_shared_region of enhancers
#     regulate a second gene.
#     B (TF-region binding, binary) is drawn from Bernoulli(b_shared_prob)
#     for the shared backbone and Bernoulli(b_ct_prob) for the cell-type layer.
#     tf_region_gene is the cross-product of B and R2G; active edges are
#     sampled from it using grn_shared_edge_prob and grn_ct_edge_prob (both
#     interpreted as fractions of G x TF, capped at the universe size).
#     W is derived from the sampled active edges.
#
#   Prior-informed mode  (prior = <data.frame>)
#     Gene and TF names are taken directly from the prior table.
#     cfg$n_tfs TFs are sampled from unique(prior$source).
#     cfg$n_genes genes are: all selected TFs + a random sample from the
#     remaining genes up to cfg$n_genes total.
#     R2G is built with an enhancer bias proportional to the number of
#     regulators each gene has in the prior; p_shared_region applies as in
#     synthetic mode.
#     W edges are sampled from prior_universe using grn_shared_edge_prob and
#     grn_ct_edge_prob (fractions of G x TF, capped at prior_universe size).
#     B is generated the same way as in synthetic mode (independent Bernoulli
#     draws), with binding probability boosted to b_prior_boost_prob at
#     positions that could structurally realise a documented prior pair, so
#     that B has the same generative status and statistical properties in
#     both modes. Each active edge is then assigned one connecting region
#     from this structural universe via preferential attachment and
#     region_reuse_prob (TFs that regulate multiple genes tend to reuse the
#     same region across those genes); pairs with no structural candidate are
#     dropped and reported via a message.

#' Build ground-truth regulatory networks for all cell types and conditions
#'
#' Constructs the latent regulatory architecture used to generate count data.
#' Three network layers are produced:
#'
#' - **TF-gene** weight matrices `W` (dimensions: genes x TFs).
#' - **TF-region-gene** triplet table `tf_region_gene` mapping each active
#'   regulatory edge to the specific region that mediates it. The RNA
#'   predictor scales each TF contribution by the accessibility of that
#'   region: `contribution = W[j,k] * A[i,k] * pi[i, r(k,j)]`.
#' - **Region-TF** weight matrices `M` (dimensions: regions x TFs).
#'   Currently set to zero; reserved for future pioneer-TF functionality.
#'
#' ## Synthetic mode (`prior = NULL`)
#'
#' Gene and TF names are generated as `"gene_1"`, `"gene_2"`, ... with
#' dimensions taken from `cfg$n_genes` and `cfg$n_tfs`. The construction
#' follows three steps:
#'
#' 1. **R2G**: each gene receives one promoter (region j for gene j) and a
#'    random number of enhancers drawn uniformly. A fraction
#'    `p_shared_region` of enhancers regulate a second gene.
#' 2. **B**: binary TF-region binding matrix drawn from
#'    `Bernoulli(b_shared_prob)` (shared backbone) and
#'    `Bernoulli(b_ct_prob)` (cell-type layer). The cross-product of B and
#'    R2G defines the universe of possible TF-region-gene triplets.
#' 3. **W**: `grn_shared_edge_prob` and `grn_ct_edge_prob` are applied as
#'    fractions of G x TF (capped at the universe size) to sample the active
#'    subset of triplets. W is derived from those active triplets.
#'
#' ## Prior-informed mode (`prior = <data.frame>`)
#'
#' Gene and TF names come from the prior table. The simulation universe is
#' built as follows:
#'
#' 1. All unique TFs (`prior$source`) and genes (`prior$source` union
#'    `prior$target`) are extracted from the validated prior.
#' 2. `cfg$n_tfs` TFs are sampled uniformly at random.
#' 3. The gene universe is: selected TFs + a random sample of non-TF genes
#'    to reach `cfg$n_genes` total.
#' 4. The prior is filtered to edges where both source and target are in
#'    this universe, giving `prior_universe`.
#' 5. **R2G**: built as in synthetic mode but with enhancer probability
#'    proportional to the number of regulators each gene has in the prior.
#' 6. **W**: `grn_shared_edge_prob` and `grn_ct_edge_prob` are applied as
#'    fractions of G x TF (capped at `nrow(prior_universe)`) to activate
#'    edges from the prior.
#' 7. **tf_region_gene**: each active edge is assigned one connecting region
#'    using preferential attachment (`b_attachment_alpha`) and optional
#'    region reuse across genes regulated by the same TF
#'    (`region_reuse_prob`).
#'
#' @param cfg   Named list of simulation parameters (see [default_config()]).
#'   All parameters are the same in both modes. In prior-informed mode
#'   `n_genes` and `n_tfs` control the size of the sampled universe, and
#'   `grn_shared_edge_prob` / `grn_ct_edge_prob` are reinterpreted as edge
#'   activation fractions rather than Bernoulli probabilities.
#' @param prior Optional `data.frame` with columns `source` (TF name),
#'   `target` (gene name) and optionally `effect` (`1` / `-1` / `NA`).
#'   Pass `NULL` (default) for synthetic construction.
#'
#' @return A named list with:
#'   \describe{
#'     \item{`genes`, `tfs`, `regions`}{Character vectors of feature names.}
#'     \item{`R2G`}{Sparse binary matrix (regions x genes). Each row may
#'           have more than one non-zero entry when `p_shared_region > 0`.}
#'     \item{`W`}{Named list of TF-gene weight matrices (genes x TFs), one
#'           per `"<cell_type>|<condition>"` key.}
#'     \item{`M`}{Named list of region-TF weight matrices (regions x TFs).
#'           Currently all-zero; reserved for pioneer-TF functionality.}
#'     \item{`tf_region_gene`}{Named list of data frames, one per
#'           `"<cell_type>|<condition>"` key. Each row is one active
#'           regulatory triplet with columns `tf_idx`, `region_idx`,
#'           `gene_idx` and `weight`.}
#'     \item{`rewiring_events`}{Tidy `data.frame` of all structural network
#'           changes introduced during rewiring.}
#'     \item{`prior_mode`}{Logical. `TRUE` if a prior was used.}
#'     \item{`prior_universe`}{Filtered prior edge table (only in prior mode).
#'           Columns: `source`, `target`, `effect`, `tf_idx`, `gene_idx`.}
#'   }
#'
#' @seealso [validate_prior()], [filter_prior_to_universe()],
#'   [rewire_matrix()], [simulate_dataset()]
build_ground_truth <- function(cfg, prior = NULL) {

  use_prior      <- !is.null(prior) && nrow(prior) > 0
  prior_universe <- NULL

  # ------------------------------------------------------------------
  # Universe definition
  # ------------------------------------------------------------------
  if (use_prior) {
    prior     <- validate_prior(prior)
    all_tfs   <- sort(unique(prior$source))
    all_genes <- sort(unique(c(prior$source, prior$target)))

    n_tfs_req <- cfg$n_tfs
    if (n_tfs_req > length(all_tfs)) {
      warning("cfg$n_tfs (", n_tfs_req, ") exceeds unique TFs in prior (",
              length(all_tfs), "). Using all available.")
      n_tfs_req <- length(all_tfs)
    }
    tfs <- sample(all_tfs, n_tfs_req)

    other_genes <- setdiff(all_genes, tfs)
    n_extra     <- cfg$n_genes - length(tfs)
    if (n_extra < 0)
      stop("cfg$n_genes must be >= cfg$n_tfs.")
    if (n_extra > length(other_genes)) {
      warning("Only ", length(other_genes), " non-TF genes available.")
      n_extra <- length(other_genes)
    }
    genes <- c(tfs, sample(other_genes, n_extra))
    G     <- length(genes)
    TF    <- length(tfs)

    prior_universe <- filter_prior_to_universe(prior, genes, tfs)
    if (nrow(prior_universe) == 0) {
      warning("No prior edges in universe. Falling back to synthetic.")
      use_prior <- FALSE
    }
  }

  if (!use_prior) {
    G     <- cfg$n_genes
    TF    <- cfg$n_tfs
    genes <- paste0("gene_", seq_len(G))
    tfs   <- genes[seq_len(TF)]
  }

  # ------------------------------------------------------------------
  # RNA-only mode: skip regions, B and tf_region_gene entirely.
  # W is built with a simple Bernoulli or prior-sampling approach,
  # identical to the pre-ATAC version of the simulator.
  # ------------------------------------------------------------------
  if (!isTRUE(cfg$simulate_atac)) {
    grn_list <- if (use_prior) {
      .build_grn_from_prior_rna_only(cfg, prior_universe, genes, tfs)
    } else {
      .build_grn_synthetic_rna_only(cfg, genes, tfs)
    }
    out <- list(
      genes           = genes,
      tfs             = tfs,
      regions         = character(0),
      W               = grn_list$W,
      M               = list(),
      tf_region_gene  = NULL,
      rewiring_events = grn_list$rewiring_events,
      prior_mode      = use_prior,
      simulate_atac   = FALSE
    )
    if (use_prior) out$prior_universe <- prior_universe
    return(out)
  }

  P       <- cfg$n_regions
  regions <- paste0("peak-", (seq_len(P) - 1L) * 500L + 1L, "-", seq_len(P) * 500L)

  if (P < G)
    stop("cfg$n_regions (", P, ") must be >= cfg$n_genes (", G, ").")

  # ------------------------------------------------------------------
  # R2G: promoters (1 per gene) + enhancers with p_shared_region
  # ------------------------------------------------------------------
  # Promoter of gene j = region j (bijection)
  i_vec <- seq_len(G)
  j_vec <- seq_len(G)

  if (P > G) {
    # In prior-informed mode with enhancer_bias_alpha > 0, genes with more
    # regulators in the prior attract enhancers with higher probability.
    # The sampling weight for gene j is 1 + alpha * n_regulators[j],
    # giving a soft bias without making it strictly proportional.
    # In synthetic mode (or when alpha = 0) enhancers are distributed uniformly.
    if (use_prior && nrow(prior_universe) > 0 &&
        !is.null(cfg$enhancer_bias_alpha) && cfg$enhancer_bias_alpha > 0) {
      n_reg     <- tabulate(match(prior_universe$target, genes), nbins = G)
      prob_enh  <- 1 + cfg$enhancer_bias_alpha * n_reg
      prob_enh  <- prob_enh / sum(prob_enh)
      enh_primary <- sample(seq_len(G), P - G, replace = TRUE, prob = prob_enh)
    } else {
      enh_primary <- sample(seq_len(G), P - G, replace = TRUE)
    }

    i_vec <- c(i_vec, (G + 1L):P)
    j_vec <- c(j_vec, enh_primary)

    # Some enhancers also regulate a second gene
    shared_mask <- runif(P - G) < cfg$p_shared_region
    if (any(shared_mask)) {
      enh_r  <- (G + which(shared_mask))
      enh_g1 <- enh_primary[shared_mask]
      enh_g2 <- vapply(enh_g1, function(g)
        sample(setdiff(seq_len(G), g), 1L), integer(1L))
      i_vec <- c(i_vec, enh_r)
      j_vec <- c(j_vec, enh_g2)
    }
  }

  R2G <- Matrix::sparseMatrix(
    i = i_vec, j = j_vec, x = 1, dims = c(P, G)
  )
  R2G_dense <- as.matrix(R2G)

  # ------------------------------------------------------------------
  # GRN (W) and tf_region_gene
  # ------------------------------------------------------------------
  if (use_prior) {
    grn_list <- .build_grn_from_prior(cfg, prior_universe, genes, tfs, R2G_dense)
  } else {
    grn_list <- .build_grn_synthetic(cfg, genes, tfs, R2G_dense)
  }

  W               <- grn_list$W
  B_per_ct        <- grn_list$B
  tf_region_gene  <- grn_list$tf_region_gene
  rewiring_events <- grn_list$rewiring_events

  # ------------------------------------------------------------------
  # M (region-TF, kept for optional pioneer TF use; currently zeros)
  # ATAC accessibility is independent of TFs in this configuration.
  # ------------------------------------------------------------------
  M <- list()
  zero_M <- matrix(0, P, TF)
  for (ct in cfg$cell_types) {
    M[[paste(ct, "Control",  sep = "|")]] <- zero_M
    M[[paste(ct, "Disease",  sep = "|")]] <- zero_M
  }

  # ------------------------------------------------------------------
  # Output
  # ------------------------------------------------------------------
  out <- list(
    genes          = genes,
    tfs            = tfs,
    regions        = regions,
    R2G            = R2G,
    W              = W,
    M              = M,
    B              = B_per_ct,
    tf_region_gene = tf_region_gene,
    rewiring_events = rewiring_events,
    prior_mode     = use_prior
  )
  if (use_prior) out$prior_universe <- prior_universe
  out
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Build tf_region_gene universe from a binary B matrix and R2G
#'
#' Returns a data.frame with columns tf_idx, region_idx, gene_idx for every
#' (k, r, j) triplet where B[r,k] = 1 AND R2G_dense[r,j] = 1.
#'
#' @keywords internal
.build_tf_region_gene <- function(B, R2G_dense) {
  B_entries <- which(B > 0, arr.ind = TRUE)   # [region, tf]
  if (nrow(B_entries) == 0L)
    return(data.frame(tf_idx     = integer(0),
                      region_idx = integer(0),
                      gene_idx   = integer(0)))

  rows <- vector("list", nrow(B_entries))
  for (e in seq_len(nrow(B_entries))) {
    r        <- B_entries[e, 1L]
    k        <- B_entries[e, 2L]
    genes_r  <- which(R2G_dense[r, ] > 0)
    if (length(genes_r) == 0L) next
    rows[[e]] <- data.frame(tf_idx     = k,
                            region_idx = r,
                            gene_idx   = genes_r)
  }
  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (length(rows) == 0L)
    return(data.frame(tf_idx     = integer(0),
                      region_idx = integer(0),
                      gene_idx   = integer(0)))
  do.call(rbind, rows)
}

#' Build W matrix from a tf_region_gene active table
#' @keywords internal
.build_W_from_tfrg <- function(active, G, TF) {
  W <- matrix(0, G, TF)
  if (nrow(active) == 0L) return(W)
  for (e in seq_len(nrow(active)))
    W[active$gene_idx[e], active$tf_idx[e]] <- active$effect[e]
  W
}

#' Update tf_region_gene active table after rewiring W
#'
#' Handles deletions, additions and weight updates consistently.
#'
#' @param struct_pool Optional `data.frame` (`tf_idx`, `region_idx`,
#'   `gene_idx`) of structurally valid candidates, used as a fallback when
#'   `inactive_pool` has no candidate for an added `(k, j)` pair. When
#'   supplied (prior-informed mode), an added edge with no candidate in
#'   either pool is dropped (zeroed in the returned `dropped` table)
#'   instead of being assigned an arbitrary region via `R2G_dense`. When
#'   `NULL` (synthetic mode, default), behaviour is unchanged.
#'
#' @return A named list with `tfrg` (the updated active table, as returned
#'   before) and `dropped` (`data.frame` with columns `gene_idx`, `tf_idx`
#'   for added edges that could not be assigned a structural region; empty
#'   when `struct_pool` is `NULL` or coverage is complete).
#'
#' @keywords internal
.update_tfrg_after_rewiring <- function(tfrg_ctrl, W_ctrl, W_dis,
                                        R2G_dense, inactive_pool,
                                        struct_pool = NULL) {
  dis_active <- tfrg_ctrl
  dropped    <- data.frame(gene_idx = integer(0), tf_idx = integer(0))

  # Deletions: edges present in Control but absent in Disease
  deleted <- which(abs(W_ctrl) > 0 & abs(W_dis) == 0, arr.ind = TRUE)
  if (nrow(deleted) > 0L) {
    del_key    <- paste(deleted[, 1L], deleted[, 2L])
    active_key <- paste(dis_active$gene_idx, dis_active$tf_idx)
    dis_active <- dis_active[!active_key %in% del_key, , drop = FALSE]
  }

  # Additions: edges absent in Control but present in Disease
  added <- which(abs(W_dis) > 0 & abs(W_ctrl) == 0, arr.ind = TRUE)
  if (nrow(added) > 0L) {
    new_rows <- vector("list", nrow(added))
    for (p in seq_len(nrow(added))) {
      j <- added[p, 1L]; k <- added[p, 2L]
      # Prefer a region from the inactive pool for this (k,j) pair
      if (!is.null(inactive_pool) && nrow(inactive_pool) > 0L) {
        cands <- inactive_pool[
          inactive_pool$tf_idx == k & inactive_pool$gene_idx == j, ,
          drop = FALSE
        ]
      } else {
        cands <- data.frame()
      }

      if (nrow(cands) > 0L) {
        r_new <- cands$region_idx[sample(nrow(cands), 1L)]
      } else if (!is.null(struct_pool)) {
        # Prior-informed mode: restrict the fallback to structurally valid
        # candidates instead of any region that happens to regulate gene j.
        cand2 <- struct_pool[
          struct_pool$tf_idx == k & struct_pool$gene_idx == j, , drop = FALSE
        ]
        if (nrow(cand2) > 0L) {
          r_new <- cand2$region_idx[sample(nrow(cand2), 1L)]
        } else {
          # No structural candidate anywhere: drop this addition rather
          # than forcing an inconsistent (tf, region, gene) triplet.
          dropped <- rbind(dropped, data.frame(gene_idx = j, tf_idx = k))
          next
        }
      } else {
        # Synthetic mode fallback: any region that regulates gene j.
        r_new <- sample(which(R2G_dense[, j] > 0), 1L)
      }
      new_rows[[p]] <- data.frame(tf_idx     = k,
                                  region_idx = r_new,
                                  gene_idx   = j,
                                  effect     = W_dis[j, k])
    }
    new_rows   <- new_rows[!vapply(new_rows, is.null, logical(1L))]
    if (length(new_rows) > 0L)
      dis_active <- rbind(dis_active, do.call(rbind, new_rows))
  }

  # Weight updates for surviving edges (drift / flip)
  surv <- which(abs(W_ctrl) > 0 & abs(W_dis) > 0, arr.ind = TRUE)
  if (nrow(surv) > 0L) {
    for (p in seq_len(nrow(surv))) {
      j   <- surv[p, 1L]; k <- surv[p, 2L]
      idx <- which(dis_active$gene_idx == j & dis_active$tf_idx == k)
      if (length(idx) > 0L)
        dis_active$effect[idx[1L]] <- W_dis[j, k]
    }
  }
  list(tfrg = dis_active, dropped = dropped)
}

#' Pick one element of a vector without `sample()`'s length-1 pitfall
#'
#' `sample(x, 1)` treats a length-1 numeric `x` as `sample(seq_len(x), 1)`,
#' which silently returns the wrong value when `x` is a single candidate
#' index. This helper always treats `x` as a value set to draw from.
#'
#' @keywords internal
.sample_one <- function(x, prob = NULL) {
  if (length(x) == 1L) return(x)
  sample(x, 1L, prob = prob)
}

#' Assign one region per TF-gene edge from the structural candidate universe
#'
#' Used in prior-informed mode where edges come from the prior and a
#' connecting region must be assigned to each. Candidates are restricted to
#' `struct_pool`, the TF-region-gene universe built from the (boosted) B
#' matrix and R2G (see [.build_tf_region_gene()]), so a region is only ever
#' assigned if the TF actually binds it in the simulated universe. Pairs
#' with no structural candidate are dropped rather than forced onto an
#' arbitrary region, and counted in the returned `n_uncovered`.
#'
#' @param struct_pool `data.frame` with columns `tf_idx`, `region_idx`,
#'   `gene_idx`: the structurally valid triplets to choose from (as
#'   returned by [.build_tf_region_gene()]).
#'
#' @return A named list with `assigned` (`data.frame` as before, with an
#'   `effect` column initialised to `0`) and `n_uncovered` (integer count
#'   of edges dropped for lack of a structural candidate).
#'
#' @keywords internal
.assign_regions_from_prior <- function(edges, struct_pool, cfg,
                                       existing_tfrg = NULL) {
  empty <- data.frame(tf_idx     = integer(0),
                      region_idx = integer(0),
                      gene_idx   = integer(0),
                      effect     = numeric(0))
  if (nrow(edges) == 0L)
    return(list(assigned = empty, n_uncovered = 0L))

  result        <- vector("list", nrow(edges))
  assigned_tfrg <- existing_tfrg
  n_uncovered   <- 0L
  pool_key      <- paste(struct_pool$tf_idx, struct_pool$gene_idx)

  for (e in seq_len(nrow(edges))) {
    k          <- edges$tf_idx[e]
    j          <- edges$gene_idx[e]
    cands      <- struct_pool[pool_key == paste(k, j), , drop = FALSE]
    regiones_j <- cands$region_idx
    if (length(regiones_j) == 0L) {
      n_uncovered <- n_uncovered + 1L
      next
    }

    r_chosen <- NULL

    # Region reuse: check if TF k already has regions that also regulate j
    if (!is.null(assigned_tfrg) && nrow(assigned_tfrg) > 0L) {
      regs_k      <- assigned_tfrg$region_idx[assigned_tfrg$tf_idx == k]
      compartidas <- intersect(regs_k, regiones_j)
      if (length(compartidas) > 0L &&
          runif(1L) < cfg$region_reuse_prob)
        r_chosen <- .sample_one(compartidas)
    }

    if (is.null(r_chosen)) {
      # Preferential attachment among structurally valid candidates
      if (!is.null(assigned_tfrg) && nrow(assigned_tfrg) > 0L) {
        counts <- tabulate(
          assigned_tfrg$region_idx[assigned_tfrg$region_idx %in% regiones_j],
          nbins = max(regiones_j)
        )[regiones_j]
      } else {
        counts <- rep(0L, length(regiones_j))
      }
      p_r      <- 1 + cfg$b_attachment_alpha * counts
      r_chosen <- .sample_one(regiones_j, prob = p_r / sum(p_r))
    }

    new_row <- data.frame(tf_idx = k, region_idx = r_chosen,
                          gene_idx = j, effect = 0)
    result[[e]]   <- new_row
    assigned_tfrg <- rbind(assigned_tfrg, new_row)
  }

  result <- result[!vapply(result, is.null, logical(1L))]
  assigned <- if (length(result) == 0L) empty else do.call(rbind, result)
  list(assigned = assigned, n_uncovered = n_uncovered)
}


# ---------------------------------------------------------------------------
# Weight sampling helper
# ---------------------------------------------------------------------------

#' Sample TF-gene edge weights for synthetic mode
#'
#' When `cfg$grn_weight_mean = 0` (default), weights are drawn from
#' `N(0, grn_weight_sd)`, a symmetric distribution centred at zero.
#' When `grn_weight_mean > 0`, each weight has a magnitude drawn from
#' `|N(grn_weight_mean, grn_weight_sd)|` and a sign sampled according to
#' `grn_repression_prob`.  This ensures all active edges have a
#' non-trivial effect size, which improves detection power in benchmarks.
#'
#' @keywords internal
.sample_grn_weights <- function(n, cfg) {
  if (n == 0L) return(numeric(0))
  wm <- if (is.null(cfg$grn_weight_mean)) 0 else cfg$grn_weight_mean
  if (wm > 0) {
    rp  <- if (is.null(cfg$grn_repression_prob)) 0.25 else cfg$grn_repression_prob
    mag <- abs(rnorm(n, wm, cfg$grn_weight_sd))
    sgn <- sample(c(-1L, 1L), n, replace = TRUE,
                  prob = c(rp, 1 - rp))
    mag * sgn
  } else {
    rnorm(n, 0, cfg$grn_weight_sd)
  }
}

# ---------------------------------------------------------------------------
# Synthetic GRN construction
# ---------------------------------------------------------------------------

.build_grn_synthetic <- function(cfg, genes, tfs, R2G_dense) {
  G  <- length(genes)
  TF <- length(tfs)
  P  <- nrow(R2G_dense)

  # Build shared B backbone
  B_shared <- matrix(runif(P * TF) < cfg$b_shared_prob, P, TF) + 0L

  # Build tf_region_gene universe from B_shared x R2G
  tfrg_shared_univ <- .build_tf_region_gene(B_shared, R2G_dense)

  # Sample shared active edges
  n_shared <- floor(G * TF * cfg$grn_shared_edge_prob)
  n_shared <- max(min(n_shared, nrow(tfrg_shared_univ)), 0L)

  if (n_shared > 0L && nrow(tfrg_shared_univ) > 0L) {
    idx           <- sample(nrow(tfrg_shared_univ), n_shared, replace = FALSE)
    shared_active <- tfrg_shared_univ[idx, , drop = FALSE]
    shared_active$effect <- .sample_grn_weights(n_shared, cfg)
    pool_shared   <- tfrg_shared_univ[-idx, , drop = FALSE]
  } else {
    shared_active        <- tfrg_shared_univ[integer(0), , drop = FALSE]
    shared_active$effect <- numeric(0)
    pool_shared          <- tfrg_shared_univ
  }

  W              <- list()
  B_per_ct       <- list()
  tf_region_gene <- list()
  rewiring_events <- data.frame()

  for (ct in cfg$cell_types) {
    # CT-specific B layer
    B_ct <- matrix(runif(P * TF) < cfg$b_ct_prob, P, TF) + 0L
    tfrg_ct_all <- .build_tf_region_gene(B_ct, R2G_dense)

    # Remove entries already in the shared universe
    if (nrow(tfrg_ct_all) > 0L && nrow(tfrg_shared_univ) > 0L) {
      sk <- paste(tfrg_shared_univ$tf_idx,
                  tfrg_shared_univ$region_idx,
                  tfrg_shared_univ$gene_idx)
      ck <- paste(tfrg_ct_all$tf_idx,
                  tfrg_ct_all$region_idx,
                  tfrg_ct_all$gene_idx)
      tfrg_ct <- tfrg_ct_all[!ck %in% sk, , drop = FALSE]
    } else {
      tfrg_ct <- tfrg_ct_all
    }

    # Pool for CT-specific sampling = shared remainder + CT-only entries
    pool_ct <- rbind(pool_shared, tfrg_ct)
    pool_ct <- pool_ct[
      !duplicated(paste(pool_ct$tf_idx, pool_ct$region_idx, pool_ct$gene_idx)),
      , drop = FALSE
    ]

    n_ct <- floor(G * TF * cfg$grn_ct_edge_prob)
    n_ct <- max(min(n_ct, nrow(pool_ct)), 0L)

    if (n_ct > 0L && nrow(pool_ct) > 0L) {
      idx_ct    <- sample(nrow(pool_ct), n_ct, replace = FALSE)
      ct_active <- pool_ct[idx_ct, , drop = FALSE]
      ct_active$effect <- .sample_grn_weights(n_ct, cfg)
    } else {
      ct_active        <- pool_ct[integer(0), , drop = FALSE]
      ct_active$effect <- numeric(0)
    }

    # Merge; shared takes priority on duplicate (tf, gene) pairs
    all_active <- rbind(shared_active, ct_active)
    dup_key    <- paste(all_active$tf_idx, all_active$gene_idx)
    all_active <- all_active[!duplicated(dup_key), , drop = FALSE]

    # Full B for this CT = shared backbone OR CT-specific sites
    B_ctrl_full <- pmax(B_shared, B_ct)
    B_per_ct[[paste(ct, "Control",  sep = "|")]] <- B_ctrl_full
    B_per_ct[[paste(ct, "Disease", sep = "|")]] <- B_ctrl_full  # disease rewiring only on W

    W_ctrl <- .build_W_from_tfrg(all_active, G, TF)
    W[[paste(ct, "Control", sep = "|")]]              <- W_ctrl
    tf_region_gene[[paste(ct, "Control", sep = "|")]] <- all_active

    if (ct %in% cfg$affected_cell_types) {
      # add_pool: inactive entries from the full universe for this CT
      full_univ <- rbind(tfrg_shared_univ, tfrg_ct)
      full_univ <- full_univ[
        !duplicated(paste(full_univ$tf_idx, full_univ$region_idx, full_univ$gene_idx)),
        , drop = FALSE
      ]
      active_key <- paste(all_active$tf_idx, all_active$gene_idx)
      pool_key   <- paste(full_univ$tf_idx, full_univ$gene_idx)
      inactive   <- full_univ[!pool_key %in% active_key, , drop = FALSE]

      if (nrow(inactive) > 0L) {
        add_pool <- cbind(inactive$gene_idx, inactive$tf_idx)
      } else {
        add_pool <- matrix(integer(0), nrow = 0L, ncol = 2L)
      }

      rw <- rewire_matrix(
        W_ctrl,
        row_names     = genes,
        col_names     = tfs,
        row_type      = "gene",
        col_type      = "tf",
        direction     = "col_to_row",
        q_del         = cfg$disease_grn_q_del,
        q_add         = cfg$disease_grn_q_add,
        q_flip        = cfg$disease_grn_q_flip,
        q_drift       = cfg$disease_grn_q_drift,
        scale_sd      = cfg$disease_grn_scale_sd,
        new_weight_sd = cfg$grn_weight_sd,
        allow_flip    = TRUE,
        add_pool      = add_pool,
        cfg           = cfg
      )
      W_dis <- rw$matrix

      upd        <- .update_tfrg_after_rewiring(
        all_active, W_ctrl, W_dis, R2G_dense, inactive
      )
      dis_active <- upd$tfrg
      if (nrow(upd$dropped) > 0L)
        W_dis[cbind(upd$dropped$gene_idx, upd$dropped$tf_idx)] <- 0

      W[[paste(ct, "Disease", sep = "|")]]              <- W_dis
      tf_region_gene[[paste(ct, "Disease", sep = "|")]] <- dis_active

      if (nrow(rw$events) > 0L) {
        rw$events$network        <- "TF_gene"
        rw$events$cell_type      <- ct
        rw$events$condition_from <- "Control"
        rw$events$condition_to   <- "Disease"
        rewiring_events <- rbind(rewiring_events, rw$events)
      }
    } else {
      W[[paste(ct, "Disease", sep = "|")]]              <- W_ctrl
      tf_region_gene[[paste(ct, "Disease", sep = "|")]] <- all_active
    }
  }

  list(W = W, B = B_per_ct, tf_region_gene = tf_region_gene, rewiring_events = rewiring_events)
}


# ---------------------------------------------------------------------------
# Prior-informed GRN construction
# ---------------------------------------------------------------------------

.build_grn_from_prior <- function(cfg, prior_universe, genes, tfs, R2G_dense) {
  G  <- length(genes)
  TF <- length(tfs)
  P  <- nrow(R2G_dense)
  N  <- nrow(prior_universe)

  # ------------------------------------------------------------------
  # TF-region binding universe (B): generated the same way as in
  # synthetic mode (independent Bernoulli draws), with probability
  # boosted to cfg$b_prior_boost_prob at (region, tf) positions that
  # could structurally realise a documented prior pair (the region
  # regulates a gene the TF is reported to regulate). This guarantees a
  # structural home for active prior edges without deriving B from
  # already-active edges, so B has the same generative status and
  # statistical properties in both modes.
  # ------------------------------------------------------------------
  prior_pair_mask <- matrix(FALSE, G, TF)
  prior_pair_mask[cbind(prior_universe$gene_idx, prior_universe$tf_idx)] <- TRUE
  boost_region_tf <- (R2G_dense %*% prior_pair_mask) > 0   # P x TF

  b_boost <- if (is.null(cfg$b_prior_boost_prob)) 0.9 else cfg$b_prior_boost_prob
  B_prob  <- matrix(cfg$b_shared_prob, P, TF)
  B_prob[boost_region_tf] <- b_boost
  B_shared <- (matrix(runif(P * TF), P, TF) < B_prob) + 0L

  tfrg_shared_univ <- .build_tf_region_gene(B_shared, R2G_dense)

  # Shared backbone from prior (pair-level sampling; same semantics as
  # before: grn_shared_edge_prob/grn_ct_edge_prob govern which documented
  # pairs are active, independent of how many structural regions connect
  # them).
  n_shared <- floor(G * TF * cfg$grn_shared_edge_prob)
  n_shared <- max(min(n_shared, N), 1L)

  idx_shared     <- sample(N, n_shared, replace = FALSE)
  shared_prior   <- prior_universe[ idx_shared, , drop = FALSE]
  pool_remaining <- prior_universe[-idx_shared, , drop = FALSE]

  W               <- list()
  tf_region_gene  <- list()
  rewiring_events <- data.frame()
  B_per_ct        <- list()
  n_uncovered     <- 0L

  for (ct in cfg$cell_types) {
    n_ct <- floor(G * TF * cfg$grn_ct_edge_prob)
    n_ct <- max(min(n_ct, nrow(pool_remaining)), 0L)

    if (n_ct > 0L) {
      idx_ct   <- sample(nrow(pool_remaining), n_ct, replace = FALSE)
      ct_prior <- pool_remaining[idx_ct, , drop = FALSE]
    } else {
      ct_prior <- pool_remaining[integer(0), , drop = FALSE]
    }

    all_prior <- rbind(shared_prior, ct_prior)
    all_prior <- all_prior[
      !duplicated(paste(all_prior$gene_idx, all_prior$tf_idx)),
      , drop = FALSE
    ]

    # Cell-type-specific binding layer, mirrors synthetic mode. No extra
    # boost needed here: B_shared already guarantees structural coverage
    # for every documented pair, regardless of whether it ends up in the
    # shared or cell-type-specific active set.
    B_ct        <- matrix(runif(P * TF) < cfg$b_ct_prob, P, TF) + 0L
    tfrg_ct_all <- .build_tf_region_gene(B_ct, R2G_dense)
    if (nrow(tfrg_ct_all) > 0L && nrow(tfrg_shared_univ) > 0L) {
      sk <- paste(tfrg_shared_univ$tf_idx, tfrg_shared_univ$region_idx,
                  tfrg_shared_univ$gene_idx)
      ck <- paste(tfrg_ct_all$tf_idx, tfrg_ct_all$region_idx,
                  tfrg_ct_all$gene_idx)
      tfrg_ct <- tfrg_ct_all[!ck %in% sk, , drop = FALSE]
    } else {
      tfrg_ct <- tfrg_ct_all
    }
    struct_pool <- rbind(tfrg_shared_univ, tfrg_ct)
    struct_pool <- struct_pool[
      !duplicated(paste(struct_pool$tf_idx, struct_pool$region_idx, struct_pool$gene_idx)),
      , drop = FALSE
    ]

    # Assign one structurally valid region per active pair, restricted to
    # struct_pool (the B x R2G universe for this cell type). Pairs with
    # no structural candidate are dropped, never forced, and counted.
    assign_out  <- .assign_regions_from_prior(all_prior, struct_pool, cfg)
    tfrg_ctrl   <- assign_out$assigned
    n_uncovered <- n_uncovered + assign_out$n_uncovered

    all_prior <- all_prior[
      paste(all_prior$gene_idx, all_prior$tf_idx) %in%
        paste(tfrg_ctrl$gene_idx, tfrg_ctrl$tf_idx),
      , drop = FALSE
    ]

    # Build W from prior edges (respects documented signs); restricted to
    # the pairs that received a structural region above.
    W_ctrl <- build_W_from_prior_edges(all_prior, G, TF, cfg$grn_weight_sd,
                                       weight_mean = if (is.null(cfg$grn_weight_mean)) 0 else cfg$grn_weight_mean)

    # Populate weights from W_ctrl
    for (e in seq_len(nrow(tfrg_ctrl))) {
      j <- tfrg_ctrl$gene_idx[e]; k <- tfrg_ctrl$tf_idx[e]
      tfrg_ctrl$effect[e] <- W_ctrl[j, k]
    }

    W[[paste(ct, "Control", sep = "|")]]              <- W_ctrl
    tf_region_gene[[paste(ct, "Control", sep = "|")]] <- tfrg_ctrl
    B_per_ct[[paste(ct, "Control", sep = "|")]] <- pmax(B_shared, B_ct)
    B_per_ct[[paste(ct, "Disease", sep = "|")]] <- pmax(B_shared, B_ct)

    # Inactive pool for rewiring, restricted to prior pairs (unchanged)
    active_key <- paste(all_prior$gene_idx, all_prior$tf_idx)
    pool_key   <- paste(prior_universe$gene_idx, prior_universe$tf_idx)
    inactive   <- prior_universe[!pool_key %in% active_key, , drop = FALSE]

    if (nrow(inactive) > 0L) {
      add_pool         <- cbind(inactive$gene_idx, inactive$tf_idx)
      add_pool_effects <- inactive$effect
    } else {
      add_pool         <- matrix(integer(0), nrow = 0L, ncol = 2L)
      add_pool_effects <- numeric(0)
    }

    prior_effects <- matrix(NA_real_, G, TF)
    prior_effects[cbind(all_prior$gene_idx, all_prior$tf_idx)] <- all_prior$effect

    if (ct %in% cfg$affected_cell_types) {
      rw <- rewire_matrix(
        W_ctrl,
        row_names        = genes,
        col_names        = tfs,
        row_type         = "gene",
        col_type         = "tf",
        direction        = "col_to_row",
        q_del            = cfg$disease_grn_q_del,
        q_add            = cfg$disease_grn_q_add,
        q_flip           = cfg$disease_grn_q_flip,
        q_drift          = cfg$disease_grn_q_drift,
        scale_sd         = cfg$disease_grn_scale_sd,
        new_weight_sd    = cfg$grn_weight_sd,
        allow_flip       = TRUE,
        add_pool         = add_pool,
        add_pool_effects = add_pool_effects,
        prior_effects    = prior_effects,
        cfg           = cfg
      )
      W_dis <- rw$matrix

      # Build inactive_tfrg pool for region assignment of added edges,
      # restricted to the same structural universe (struct_pool) as
      # tfrg_ctrl above.
      inactive_assign <- .assign_regions_from_prior(
        inactive, struct_pool, cfg, existing_tfrg = tfrg_ctrl
      )
      inactive_tfrg <- inactive_assign$assigned
      n_uncovered   <- n_uncovered + inactive_assign$n_uncovered

      upd        <- .update_tfrg_after_rewiring(
        tfrg_ctrl, W_ctrl, W_dis, R2G_dense, inactive_tfrg,
        struct_pool = struct_pool
      )
      dis_active <- upd$tfrg
      if (nrow(upd$dropped) > 0L) {
        W_dis[cbind(upd$dropped$gene_idx, upd$dropped$tf_idx)] <- 0
        n_uncovered <- n_uncovered + nrow(upd$dropped)
      }

      W[[paste(ct, "Disease", sep = "|")]]              <- W_dis
      tf_region_gene[[paste(ct, "Disease", sep = "|")]] <- dis_active

      if (nrow(rw$events) > 0L) {
        rw$events$network        <- "TF_gene"
        rw$events$cell_type      <- ct
        rw$events$condition_from <- "Control"
        rw$events$condition_to   <- "Disease"
        rewiring_events <- rbind(rewiring_events, rw$events)
      }
    } else {
      W[[paste(ct, "Disease", sep = "|")]]              <- W_ctrl
      tf_region_gene[[paste(ct, "Disease", sep = "|")]] <- tfrg_ctrl
    }
  }

  if (n_uncovered > 0L)
    message(
      n_uncovered, " prior edge instance(s) dropped: no structurally valid ",
      "region connected the TF to the gene in the simulated B x R2G ",
      "universe. Increase cfg$b_prior_boost_prob, cfg$b_shared_prob or ",
      "cfg$n_regions to reduce this."
    )

  list(W = W, B = B_per_ct, tf_region_gene = tf_region_gene, rewiring_events = rewiring_events)
}


# ---------------------------------------------------------------------------
# RNA-only GRN builders (simulate_atac = FALSE)
# These are simplified versions of the ATAC-aware builders that produce
# W directly without regions, B or tf_region_gene.
# ---------------------------------------------------------------------------

#' Build W by Bernoulli sampling (no ATAC, no tf_region_gene)
#' @keywords internal
.build_grn_synthetic_rna_only <- function(cfg, genes, tfs) {
  G  <- length(genes)
  TF <- length(tfs)

  shared_mask <- matrix(runif(G * TF) < cfg$grn_shared_edge_prob, G, TF)
  shared_W    <- matrix(0, G, TF)
  shared_W[shared_mask] <- .sample_grn_weights(sum(shared_mask), cfg)

  W               <- list()
  rewiring_events <- data.frame()

  for (ct in cfg$cell_types) {
    ct_mask <- matrix(runif(G * TF) < cfg$grn_ct_edge_prob, G, TF)
    W_ctrl  <- shared_W
    W_ctrl[ct_mask] <- .sample_grn_weights(sum(ct_mask), cfg)
    W[[paste(ct, "Control", sep = "|")]] <- W_ctrl

    if (ct %in% cfg$affected_cell_types) {
      rw <- rewire_matrix(
        W_ctrl,
        row_names     = genes, col_names     = tfs,
        row_type      = "gene", col_type     = "tf",
        direction     = "col_to_row",
        q_del         = cfg$disease_grn_q_del,
        q_add         = cfg$disease_grn_q_add,
        q_flip        = cfg$disease_grn_q_flip,
        q_drift       = cfg$disease_grn_q_drift,
        scale_sd      = cfg$disease_grn_scale_sd,
        new_weight_sd = cfg$grn_weight_sd,
        allow_flip    = TRUE,
        cfg           = cfg
      )
      W[[paste(ct, "Disease", sep = "|")]] <- rw$matrix
      if (nrow(rw$events) > 0L) {
        rw$events$network        <- "TF_gene"
        rw$events$cell_type      <- ct
        rw$events$condition_from <- "Control"
        rw$events$condition_to   <- "Disease"
        rewiring_events <- rbind(rewiring_events, rw$events)
      }
    } else {
      W[[paste(ct, "Disease", sep = "|")]] <- W_ctrl
    }
  }
  list(W = W, rewiring_events = rewiring_events)
}

#' Build W from prior_universe without region assignment (no ATAC, no tf_region_gene)
#' @keywords internal
.build_grn_from_prior_rna_only <- function(cfg, prior_universe, genes, tfs) {
  G  <- length(genes)
  TF <- length(tfs)
  N  <- nrow(prior_universe)

  n_shared   <- max(min(floor(G * TF * cfg$grn_shared_edge_prob), N), 1L)
  idx_shared <- sample(N, n_shared, replace = FALSE)
  shared_prior   <- prior_universe[ idx_shared, , drop = FALSE]
  pool_remaining <- prior_universe[-idx_shared, , drop = FALSE]

  W               <- list()
  rewiring_events <- data.frame()

  for (ct in cfg$cell_types) {
    n_ct <- max(min(floor(G * TF * cfg$grn_ct_edge_prob), nrow(pool_remaining)), 0L)
    ct_prior <- if (n_ct > 0L)
      pool_remaining[sample(nrow(pool_remaining), n_ct, replace = FALSE), , drop = FALSE]
    else
      pool_remaining[integer(0L), , drop = FALSE]

    all_prior <- rbind(shared_prior, ct_prior)
    all_prior <- all_prior[
      !duplicated(paste(all_prior$gene_idx, all_prior$tf_idx)), , drop = FALSE
    ]

    W_ctrl <- build_W_from_prior_edges(all_prior, G, TF, cfg$grn_weight_sd,
                                       weight_mean = if (is.null(cfg$grn_weight_mean)) 0 else cfg$grn_weight_mean)
    W[[paste(ct, "Control", sep = "|")]] <- W_ctrl

    active_key <- paste(all_prior$gene_idx, all_prior$tf_idx)
    pool_key   <- paste(prior_universe$gene_idx, prior_universe$tf_idx)
    inactive   <- prior_universe[!pool_key %in% active_key, , drop = FALSE]

    add_pool         <- if (nrow(inactive) > 0L) cbind(inactive$gene_idx, inactive$tf_idx) else matrix(integer(0L), 0L, 2L)
    add_pool_effects <- if (nrow(inactive) > 0L) inactive$effect else numeric(0L)

    prior_effects <- matrix(NA_real_, G, TF)
    prior_effects[cbind(all_prior$gene_idx, all_prior$tf_idx)] <- all_prior$effect

    if (ct %in% cfg$affected_cell_types) {
      rw <- rewire_matrix(
        W_ctrl,
        row_names        = genes, col_names        = tfs,
        row_type         = "gene", col_type        = "tf",
        direction        = "col_to_row",
        q_del            = cfg$disease_grn_q_del,
        q_add            = cfg$disease_grn_q_add,
        q_flip           = cfg$disease_grn_q_flip,
        q_drift          = cfg$disease_grn_q_drift,
        scale_sd         = cfg$disease_grn_scale_sd,
        new_weight_sd    = cfg$grn_weight_sd,
        allow_flip       = TRUE,
        add_pool         = add_pool,
        add_pool_effects = add_pool_effects,
        prior_effects    = prior_effects,
        cfg           = cfg
      )
      W[[paste(ct, "Disease", sep = "|")]] <- rw$matrix
      if (nrow(rw$events) > 0L) {
        rw$events$network        <- "TF_gene"
        rw$events$cell_type      <- ct
        rw$events$condition_from <- "Control"
        rw$events$condition_to   <- "Disease"
        rewiring_events <- rbind(rewiring_events, rw$events)
      }
    } else {
      W[[paste(ct, "Disease", sep = "|")]] <- W_ctrl
    }
  }
  list(W = W, rewiring_events = rewiring_events)
}

# ---------------------------------------------------------------------------
# Edge list serialisers
# ---------------------------------------------------------------------------

#' Convert a TF-gene weight matrix to a tidy edge list
#'
#' Extracts all non-zero entries of `W` and returns them as a `data.frame`
#' annotated with condition and cell-type labels.  When `prior_universe` is
#' supplied, each edge is also annotated with whether it appears in the prior
#' knowledge base (`in_prior_db`) and its documented effect direction
#' (`db_effect`).
#'
#' @param W             Numeric matrix (genes x TFs) of regulatory weights.
#' @param genes         Character vector of gene names (row names of `W`).
#' @param tfs           Character vector of TF names (column names of `W`).
#' @param condition     Character scalar. Condition label.
#' @param cell_type     Character scalar. Cell-type label.
#' @param prior_universe Optional filtered prior `data.frame` as returned
#'   by [filter_prior_to_universe()].  When supplied, adds columns
#'   `in_prior_db` and `db_effect`.
#'
#' @return A `data.frame` with columns `network`, `condition`, `cell_type`,
#'   `tf`, `gene`, `effect`, and optionally `in_prior_db`, `db_effect`.
#'
#' @keywords internal
matrix_to_grn_edges <- function(W, genes, tfs, condition, cell_type,
                                prior_universe = NULL) {
  idx <- which(abs(W) > 0, arr.ind = TRUE)
  if (nrow(idx) == 0L)
    return(data.frame(network = character(0), condition = character(0),
                      cell_type = character(0), tf = character(0),
                      gene = character(0), effect = numeric(0),
                      stringsAsFactors = FALSE))
  df  <- data.frame(
    network   = "TF_gene",
    condition = condition,
    cell_type = cell_type,
    tf        = tfs[idx[, 2L]],
    gene      = genes[idx[, 1L]],
    effect    = W[idx],
    stringsAsFactors = FALSE
  )
  if (!is.null(prior_universe) && nrow(prior_universe) > 0L) {
    prior_key <- paste(prior_universe$source, prior_universe$target, sep = "|")
    edge_key  <- paste(df$tf, df$gene, sep = "|")
    in_db     <- edge_key %in% prior_key
    db_eff    <- rep(NA_real_, nrow(df))
    m         <- match(edge_key[in_db], prior_key)
    db_eff[in_db]  <- prior_universe$effect[m]
    df$in_prior_db <- in_db
    df$db_effect   <- db_eff
  }
  df
}

#' Convert a region-TF weight matrix to a tidy edge list
#'
#' @param M         Numeric matrix (regions x TFs) of accessibility weights.
#' @param regions   Character vector of region names.
#' @param tfs       Character vector of TF names.
#' @param condition Character scalar. Condition label.
#' @param cell_type Character scalar. Cell-type label.
#'
#' @return A `data.frame` with columns `network`, `condition`, `cell_type`,
#'   `region`, `tf`, `effect`.
#'
#' @keywords internal
matrix_to_region_tf_edges <- function(M, regions, tfs, condition, cell_type) {
  empty <- data.frame(network = character(0), condition = character(0),
                      cell_type = character(0), region = character(0),
                      tf = character(0), effect = numeric(0),
                      stringsAsFactors = FALSE)
  # Guard against NULL M (e.g. pioneers not yet implemented) or all-zero M
  if (is.null(M)) return(empty)
  idx <- which(abs(M) > 0, arr.ind = TRUE)
  if (is.null(nrow(idx)) || nrow(idx) == 0L) return(empty)
  data.frame(
    network   = "region_TF",
    condition = condition,
    cell_type = cell_type,
    region    = regions[idx[, 1L]],
    tf        = tfs[idx[, 2L]],
    effect    = M[idx],
    stringsAsFactors = FALSE
  )
}