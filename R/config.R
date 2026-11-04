# config.R
# Default simulation configuration.

#' Return the default simulation configuration
#'
#' Provides a fully-specified named list that can be passed directly to
#' [simulate_dataset()].  Every field is documented below.  Users typically
#' call this function, modify the fields they need and pass the result to the
#' simulator.
#'
#' @return A named list with the following fields (grouped by concern):
#'
#' **Reproducibility**
#' \describe{
#'   \item{`seed`}{Integer random seed.  Default: `42`.}
#' }
#'
#' **Experimental design**
#' \describe{
#'   \item{`conditions`}{Character vector of condition labels.  The first
#'         element is treated as the reference.  Default:
#'         `c("Control", "Disease")`.}
#'   \item{`cell_types`}{Character vector of cell-type labels.  Default:
#'         `c("TypeA", "TypeB")`.}
#'   \item{`affected_cell_types`}{Subset of `cell_types` that will have
#'         condition-specific GRN rewiring.  Default: `"TypeA"`.}
#'   \item{`n_reps_per_condition`}{Number of biological replicates per
#'         condition.  Default: `2`.}
#'   \item{`n_cells_per_type`}{Named integer vector of cells per cell type
#'         per replicate per condition.  Default: 200 cells for each type.}
#' }
#'
#' **Feature dimensions**
#' \describe{
#'   \item{`n_genes`}{Number of genes.  Default: `200`.}
#'   \item{`n_tfs`}{Number of transcription factors (subset of the first
#'         `n_tfs` genes in synthetic mode; sampled from the prior in
#'         prior-informed mode).  Default: `20`.}
#'   \item{`n_regions`}{Number of chromatin accessibility peaks.  Must be
#'         >= `n_genes` so that every gene receives at least one promoter
#'         region.  Default: `400`.}
#' }
#'
#' **ATAC simulation**
#' \describe{
#'   \item{`simulate_atac`}{Logical.  When `TRUE` (default), paired scATAC-seq
#'         counts are generated and the RNA predictor uses a multiplicative
#'         gate equal to the accessibility of the connecting peak:
#'         `contribution = W[j,k] * A[i,k] * pi[i, r(k,j)]`.  When `FALSE`,
#'         no ATAC matrix is produced, no peaks or B matrix are constructed,
#'         and the RNA predictor reduces to the ungated form
#'         `grn_drive = A * W`, identical to a GRN-only simulation.
#'         Default: `TRUE`.}
#' }
#'
#' **ATAC accessibility model**
#' \describe{
#'   \item{`atac_regulator_boost`}{Additive term on the accessibility logit
#'         (`logit_pi`) for regions that mediate at least one active
#'         TF-gene edge (i.e. appear in `tf_region_gene`) in a given
#'         `(cell_type, condition)` group. Regions with no regulatory role
#'         in that group keep the unboosted baseline. Because active edges
#'         differ between Control and Disease (via rewiring), this also
#'         gives the ATAC ground truth a genuine condition-specific
#'         accessibility signal, tracking whichever edges are active.
#'         Set to `0` to disable (accessibility becomes purely
#'         non-regulatory background, as before). Default: `2.0`.}
#'   \item{`atac_accessibility_noise_sd`}{SD of per-cell Gaussian noise
#'         added to `logit_pi`. Without this, accessibility is identical
#'         for every cell of the same cell type and replicate, so any
#'         group-level signal (including `atac_regulator_boost`) would be
#'         a deterministic label recoverable from pseudobulk averages
#'         alone. This noise creates genuine cell-to-cell overlap between
#'         regulating and non-regulating regions, analogous to
#'         `tf_activity_noise_sd` for TF activity. Default: `0.5`.}
#'   \item{`debug_export_pi`}{Logical. When `TRUE`, the latent accessibility
#'         matrix `pi` (cells x regions, the sigmoid of `logit_pi`, before
#'         the Poisson count / dropout / library-size steps that produce
#'         the exported `ATAC` matrix) is accumulated per chunk and
#'         returned as `sim$debug_pi`. This is the value that actually
#'         gates `grn_drive`, so it is useful for diagnosing the
#'         TF-accessibility-RNA mechanism without the measurement noise of
#'         Poisson counts and dropout on a single peak. Has no effect on
#'         the simulated RNA/ATAC data, only adds this extra output.
#'         Ignored when `simulate_atac = FALSE`. Default: `FALSE`.}
#' }
#'
#' **GRN architecture (W matrix, TF-gene)**
#' \describe{
#'   \item{`grn_shared_edge_prob`}{Controls the density of the shared
#'         backbone (active in all cell types).  In both modes the target
#'         number of active edges is `floor(n_genes * n_tfs * prob)`.
#'         In synthetic mode this is drawn from the tf_region_gene universe
#'         (cross-product of B and R2G) and capped at its size.  In
#'         prior-informed mode it is drawn from `prior_universe` and capped
#'         at `nrow(prior_universe)`.  Default: `0.05`.}
#'   \item{`grn_ct_edge_prob`}{Cell-type-specific edge density, interpreted
#'         the same way as `grn_shared_edge_prob`.  Default: `0.02`.}
#'   \item{`grn_weight_mean`}{Mean magnitude of TF-gene edge weights.
#'         When `0` (default), weights are drawn from `N(0, grn_weight_sd)`,
#'         a symmetric distribution centred at zero.  When `> 0`, each edge
#'         is assigned a magnitude drawn from `|N(grn_weight_mean,
#'         grn_weight_sd)|` and a sign sampled with probability
#'         `grn_repression_prob` for repression and `1 -
#'         grn_repression_prob` for activation.  This guarantees non-trivial
#'         effect sizes and avoids the many near-zero weights that arise
#'         from a mean-zero distribution.  In prior-informed mode the sign
#'         always comes from the prior (or is sampled randomly for edges
#'         with unknown direction), so only the magnitude is affected.
#'         Default: `0`.}
#'   \item{`grn_repression_prob`}{Fraction of synthetic edges that are
#'         repressive (negative weight).  Only used when
#'         `grn_weight_mean > 0`.  In prior-informed mode the sign is
#'         determined by the prior, so this parameter has no effect.
#'         A value around `0.25` reflects the typical activator-dominated
#'         composition of characterised GRNs.  Default: `0.25`.}
#'   \item{`grn_weight_sd`}{SD of the weight magnitude distribution.
#'         When `grn_weight_mean = 0` this is the SD of the full
#'         `N(0, sd)` distribution.  When `grn_weight_mean > 0` this is
#'         the SD around the mean magnitude before the sign is applied.
#'         Default: `0.3`.}
#' }
#'
#' **TF-region binding (B matrix)**
#' \describe{
#'   \item{`b_shared_prob`}{Probability that a TF has a binding motif in a
#'         given region (shared backbone, same for all cell types).  Together
#'         with `n_regions` and `p_shared_region` this determines the
#'         richness of the tf_region_gene universe from which active edges
#'         are sampled.  Default: `0.05`.}
#'   \item{`b_ct_prob`}{Probability of a cell-type-specific TF-region
#'         binding site (layer added on top of the shared backbone).
#'         Default: `0.02`.}
#'   \item{`b_prior_boost_prob`}{Only used in prior-informed mode.
#'         Binding probability applied instead of `b_shared_prob` at
#'         `(region, tf)` positions that could structurally realise a
#'         documented prior pair (i.e. the region regulates a gene that
#'         the TF is reported to regulate).  This guarantees that active
#'         prior edges have a real region to attach to, while keeping `B`
#'         a single generative (Bernoulli) object shared with synthetic
#'         mode, rather than a matrix reverse-engineered from already
#'         active edges.  Higher values reduce the fraction of prior
#'         edges dropped for lack of a structural candidate (see
#'         `simulate_dataset()` messages for the observed drop rate).
#'         Default: `0.9`.}
#' }
#'
#' **Region-gene links (R2G)**
#' \describe{
#'   \item{`p_shared_region`}{Probability that an enhancer region regulates
#'         two genes instead of one.  Set to `0` to enforce a strict one
#'         region to one gene mapping.  When `> 0`, TFs that regulate
#'         multiple genes have a higher chance of acting through the same
#'         peak, reflecting the biology of shared enhancers.
#'         Default: `0.10`.}
#'   \item{`enhancer_bias_alpha`}{Only used in prior-informed mode.
#'         Controls how strongly the number of regulators a gene has in
#'         the prior influences its probability of receiving enhancer
#'         regions.  The sampling probability for gene j is proportional
#'         to `1 + enhancer_bias_alpha * n_regulators[j]`.
#'         Set to `0` for uniform distribution (same as synthetic mode).
#'         Default: `0.25`.}
#' }
#'
#' **Prior-informed mode: region assignment**
#' \describe{
#'   \item{`region_reuse_prob`}{Only used in prior-informed mode.
#'         Probability that when assigning a connecting region to a TF-gene
#'         edge, a region already used by that TF for another gene is
#'         reused (provided it also regulates the target gene).  This
#'         increases the chance that two genes regulated by the same TF
#'         share the same peak.  Default: `0.30`.}
#'   \item{`b_attachment_alpha`}{Only used in prior-informed mode.  Strength
#'         of the preferential attachment when choosing a region for a
#'         TF-gene edge: the probability of choosing region r is
#'         proportional to `1 + alpha * n_tfs_already_in_r`.  Higher values
#'         concentrate binding sites on a smaller number of busy regions
#'         (promoter-like hubs).  `0` gives uniform selection.
#'         Default: `1.0`.}
#' }
#'
#' **Disease rewiring: TF-gene network**
#' \describe{
#'   \item{`disease_grn_q_del`}{Fraction of edges deleted.  Default: `0.15`.}
#'   \item{`disease_grn_q_add`}{Fraction of edges added.  Default: `0.10`.}
#'   \item{`disease_grn_q_flip`}{Fraction of surviving edges whose sign is
#'         flipped (activation to repression or vice versa).  In
#'         prior-informed mode only edges with unknown direction
#'         (`effect = NA`) are eligible.  Default: `0.05`.}
#'   \item{`disease_grn_q_drift`}{Fraction of surviving edges that undergo
#'         continuous weight drift (multiplied by a log-Normal factor).
#'         Edges not selected keep their exact Control weight.  Set to `1`
#'         to drift all edges.  Default: `0.30`.}
#'   \item{`disease_grn_scale_sd`}{SD of the log-Normal weight drift:
#'         each selected edge weight is multiplied by `exp(N(0, scale_sd))`.
#'         Default: `0.25`.}
#' }
#'
#' **Disease rewiring: region-TF network (M matrix)**
#' \describe{
#'   \item{`disease_region_tf_q_del`}{Fraction of edges deleted.
#'         Default: `0.10`.}
#'   \item{`disease_region_tf_q_add`}{Fraction of edges added.
#'         Default: `0.08`.}
#'   \item{`disease_region_tf_q_drift`}{Fraction of surviving edges that
#'         undergo weight drift.  Sign-flipping is disabled for this layer.
#'         Default: `0.30`.}
#'   \item{`disease_region_tf_scale_sd`}{SD of the log-Normal weight drift.
#'         Default: `0.20`.}
#' }
#'
#' **TF activity**
#' \describe{
#'   \item{`tf_activity_ct_sigma`}{SD of the cell-type-level TF activity
#'         baseline.  Default: `0.5`.}
#'   \item{`disease_tf_activity_sigma`}{Additional TF activity shift SD in
#'         affected cell types under Disease.  Default: `0.8`.}
#'   \item{`tf_activity_noise_sd`}{Per-cell TF activity noise SD.
#'         Default: `0.3`.}
#' }
#'
#' **Cell-type expression shifts**
#' \describe{
#'   \item{`ct_gene_sigma`}{SD of the per-cell-type gene expression shift.
#'         Default: `0.5`.}
#'   \item{`ct_region_sigma`}{SD of the per-cell-type region accessibility
#'         shift.  Default: `0.4`.}
#' }
#'
#' **Batch effects**
#' \describe{
#'   \item{`replicate_gene_sigma`}{Gene-level replicate batch effect SD.
#'         Default: `0.1`.}
#'   \item{`replicate_region_sigma`}{Region-level replicate batch effect SD.
#'         Default: `0.1`.}
#'   \item{`replicate_depth_sigma`}{Replicate-level library-size variability
#'         (log scale).  Default: `0.1`.}
#'   \item{`replicate_dropout_sigma`}{Replicate-level variability in the
#'         dropout intercept.  Default: `0.1`.}
#'   \item{`replicate_tf_sigma`}{SD of the per-sample TF activity batch
#'         effect.  Each sample receives a shift drawn from
#'         `N(0, replicate_tf_sigma)` on every TF's activity.  Set to `0`
#'         (default) to disable.  Default: `0`.}
#'   \item{`tf_batch_propagates_to_rna`}{Logical.  Controls whether the TF
#'         batch effect enters the regulatory model.  `TRUE` (Version A):
#'         the batch-shifted activity drives both RNA and ATAC, creating
#'         correlated batch effects between co-regulated genes.  `FALSE`
#'         (Version B): the batch is an observation bias only in the
#'         exported `TFActivity` matrix; RNA and ATAC are driven by
#'         batch-free activity.  Default: `TRUE`.}
#' }
#'
#' **Sequencing depths**
#' \describe{
#'   \item{`rna_depth_mean`}{Mean RNA library size per cell.
#'         Default: `3000`.}
#'   \item{`rna_depth_logsd`}{Log-SD of RNA library size.  Default: `0.4`.}
#'   \item{`atac_depth_mean`}{Mean ATAC library size per cell.
#'         Default: `2000`.}
#'   \item{`atac_depth_logsd`}{Log-SD of ATAC library size.  Default: `0.4`.}
#'   \item{`depth_corr`}{Pearson correlation between log RNA and log ATAC
#'         depths within a cell.  Default: `0.4`.}
#' }
#'
#' **NB dispersion**
#' \describe{
#'   \item{`rna_theta_mean`}{Mean of the per-gene NB dispersion parameter.
#'         Default: `5`.}
#'   \item{`rna_theta_logsd`}{Log-SD of the per-gene NB dispersion.
#'         Default: `0.5`.}
#' }
#'
#' **Regulatory signal strengths**
#' \describe{
#'   \item{`grn_strength`}{Multiplier for the GRN drive term in the RNA
#'         linear predictor.  The drive is automatically normalised by the
#'         square root of the mean number of active edges per gene, so
#'         `grn_strength` acts as a calibration scalar.  Default: `1.0`.}
#'   \item{`atac_tf_strength`}{Multiplier for the pioneer-TF drive in the
#'         ATAC logit predictor (via M matrix).  Currently M is set to zero
#'         so this parameter has no effect; reserved for future use.
#'         Default: `1.0`.}
#' }
#'
#' **Dropout models**
#' \describe{
#'   \item{`rna_dropout_intercept`}{Logistic dropout intercept for RNA.
#'         More positive values give more dropout.  Default: `2`.}
#'   \item{`rna_dropout_slope`}{Slope of the log-mean term in the RNA
#'         dropout model.  Default: `0.5`.}
#'   \item{`atac_dropout_intercept`, `atac_dropout_slope`}{Superseded by
#'         `atac_capture_intercept`/`atac_capture_slope` (see "ATAC count
#'         model" below), which now govern ATAC detection instead of a
#'         separate dropout step. Kept for backward compatibility with
#'         old configs but no longer used by [simulate_dataset()].}
#'   \item{`binarize_atac`}{Logical.  Binarise ATAC counts after sampling?
#'         Default: `FALSE`.}
#' }
#'
#' **ATAC count model**
#' \describe{
#'   \item{`atac_capture_intercept`, `atac_capture_slope`}{ATAC counts are
#'         generated in two steps: (1) each region has up to 2 accessible
#'         copies (alleles), `n_open ~ Binomial(2, pi)`, driven purely by
#'         the accessibility gate `pi`; (2) each open copy is
#'         independently captured/sequenced with per-cell probability
#'         `p_capture = sigmoid(atac_capture_intercept +
#'         atac_capture_slope * log(lib_atac))`, so sequencing depth
#'         affects detection sensitivity without inflating counts beyond
#'         the physically bounded `n_open` (unlike a shared read-budget
#'         Poisson model, which is appropriate for RNA but not for a
#'         near-binary, per-locus signal like chromatin accessibility).
#'         `replicate_dropout_sigma` still adds per-cell noise to the
#'         capture logit, as it previously did to the dropout logit.
#'         Defaults (`atac_capture_intercept = -7.5`,
#'         `atac_capture_slope = 0.8`) target `p_capture ~ 0.5` at
#'         `atac_depth_mean`; recalibrate empirically for your depth and
#'         desired sensitivity, the same way as `atac_regulator_boost`.}
#' }
#'
#' **Performance**
#' \describe{
#'   \item{`chunk_size`}{Number of cells processed per chunk.  Reduce to
#'         lower peak memory usage.  Default: `500`.}
#' }
#'
#' @export
#'
#' @examples
#' cfg <- default_config()
#' cfg$n_genes <- 500
#' cfg$seed    <- 123
#' # sim <- simulate_dataset(cfg)

default_config <- function() {
  list(
    # Reproducibility
    seed                      = 2808L,

    # Experimental design
    conditions                = c("Control", "Disease"),
    cell_types                = c("Tcell", "Bcell", "Mono"),
    affected_cell_types       = c("Tcell", "Mono"),
    n_reps_per_condition      = 2L,
    n_cells_per_type          = c(Tcell = 800L, Bcell = 600L, Mono = 700L),

    # Feature dimensions
    n_genes                   = 2000L,
    n_tfs                     = 120L,
    n_regions                 = 6000L,

    # ATAC simulation
    simulate_atac             = TRUE,

    # ATAC accessibility model
    atac_regulator_boost        = 2.0,
    atac_accessibility_noise_sd = 0.5,
    debug_export_pi              = FALSE,

    # GRN architecture (W matrix, TF-gene)
    # In synthetic mode: active edges are sampled from the tf_region_gene
    # universe (cross-product of B and R2G); grn_shared_edge_prob targets
    # floor(G * TF * prob) edges, capped at the universe size.
    # In prior-informed mode: same target applied to prior_universe.
    grn_shared_edge_prob      = 0.001,
    grn_ct_edge_prob          = 0.001,
    grn_weight_mean           = 0,     # 0 = N(0, sd); >0 enables bimodal sampling
    grn_repression_prob       = 0.25,  # fraction of repressive edges (used when grn_weight_mean > 0)
    grn_weight_sd             = 0.3,

    # TF-region binding (B matrix, binary)
    b_shared_prob             = 0.05,
    b_ct_prob                 = 0.02,
    b_prior_boost_prob        = 0.9,   # prior mode only: coverage boost, see docs above

    # Region-gene links (R2G)
    p_shared_region           = 0.10,
    enhancer_bias_alpha       = 0.25,  # prior mode only: bias enhancers toward highly regulated genes

    # Prior-informed mode: region assignment
    region_reuse_prob         = 0.30,
    b_attachment_alpha        = 1.0,

    # Disease rewiring — TF-gene (W)
    disease_grn_q_del         = 0.25,
    disease_grn_q_add         = 0.25,
    disease_grn_q_flip        = 0.05,
    disease_grn_q_drift       = 0.05,
    disease_grn_scale_sd      = 0.35,

    # Region-TF architecture (M matrix; currently zero, reserved for pioneers)
    region_tf_shared_prob     = 0.004,
    region_tf_ct_prob         = 0.002,
    region_tf_weight_mean     = 0.5,
    region_tf_weight_sd       = 0.25,

    # Disease rewiring — region-TF (M)
    disease_region_tf_q_del   = 0.15,
    disease_region_tf_q_add   = 0.15,
    disease_region_tf_q_drift = 0.05,
    disease_region_tf_scale_sd = 0.25,

    # TF activity
    tf_activity_ct_sigma      = 0.7,
    disease_tf_activity_sigma = 0.6,
    tf_activity_noise_sd      = 0.35,

    # Cell-type expression shifts
    ct_gene_sigma             = 0.35,
    ct_region_sigma           = 0.35,

    # Batch effects
    replicate_gene_sigma      = 0.12,
    replicate_region_sigma    = 0.12,
    replicate_depth_sigma     = 0.25,
    replicate_dropout_sigma   = 0.20,
    replicate_tf_sigma        = 0.0,
    tf_batch_propagates_to_rna = TRUE,
    # Version A (TRUE): TF batch propagates into
    # RNA/ATAC expression, creating correlated
    # batch effects between co-regulated genes.
    # Version B (FALSE): batch is only an
    # observation bias in TFActivity; RNA/ATAC
    # are driven by batch-free activity.

    # Sequencing depths
    rna_depth_mean            = 6000,
    rna_depth_logsd           = 0.45,
    atac_depth_mean           = 12000,
    atac_depth_logsd          = 0.55,
    depth_corr                = 0.45,

    # NB dispersion
    rna_theta_mean            = 5,
    rna_theta_logsd           = 0.8,

    # Signal strengths
    grn_strength              = 0.9,
    atac_tf_strength          = 0.7,

    # Dropout
    rna_dropout_intercept     = 2,
    rna_dropout_slope         = 0.55,
    atac_dropout_intercept    = 2.5,   # superseded, see atac_capture_* below
    atac_dropout_slope        = 0.4,   # superseded, see atac_capture_* below
    binarize_atac             = FALSE,

    # ATAC count model (capture step; see roxygen doc above)
    atac_capture_intercept    = -7.5,
    atac_capture_slope        = 0.8,

    # Performance
    chunk_size                = 500L
  )
}
