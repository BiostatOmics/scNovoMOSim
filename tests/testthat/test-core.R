test_that("sigmoid output is in (0, 1) and is numerically stable", {
  vals <- sigmoid(c(-1000, -30, 0, 30, 1000))
  expect_true(all(vals > 0 & vals < 1))
  expect_equal(sigmoid(0), 0.5)
})

test_that("nb_counts returns a non-negative integer matrix", {
  mu    <- matrix(c(1, 5, 10, 20), nrow = 2, ncol = 2)
  theta <- c(5, 5)
  out   <- nb_counts(mu, theta)
  expect_equal(dim(out), dim(mu))
  expect_true(all(out >= 0))
  expect_true(is.numeric(out))
})

test_that("make_metadata returns correct number of rows", {
  cfg <- default_config()
  meta <- make_metadata(cfg)
  expected_rows <- length(cfg$conditions) *
    cfg$n_reps_per_condition *
    sum(unlist(cfg$n_cells_per_type))
  expect_equal(nrow(meta), expected_rows)
})

test_that("make_metadata index columns are valid", {
  cfg  <- default_config()
  meta <- make_metadata(cfg)
  expect_true(all(meta$condition_idx %in% seq_along(cfg$conditions)))
  expect_true(all(meta$cell_type_idx %in% seq_along(cfg$cell_types)))
})

test_that("default_config returns a list with expected fields", {
  cfg <- default_config()
  expect_type(cfg, "list")
  expect_true("seed"         %in% names(cfg))
  expect_true("n_genes"      %in% names(cfg))
  expect_true("n_tfs"        %in% names(cfg))
  expect_true("n_regions"    %in% names(cfg))
  expect_true("conditions"   %in% names(cfg))
  expect_true("cell_types"   %in% names(cfg))
  expect_true("chunk_size"   %in% names(cfg))
})

test_that("build_ground_truth returns correct structure", {
  cfg <- default_config()
  cfg$n_genes   <- 30
  cfg$n_tfs     <- 5
  cfg$n_regions <- 20
  gt <- build_ground_truth(cfg)
  expect_equal(length(gt$genes),   30)
  expect_equal(length(gt$tfs),      5)
  expect_equal(length(gt$regions), 20)
  # One W and M matrix per cell-type x condition combination
  expected_keys <- length(cfg$cell_types) * length(cfg$conditions)
  expect_equal(length(gt$W), expected_keys)
  expect_equal(length(gt$M), expected_keys)
})

test_that("rewire_matrix changes edge structure and records events", {
  set.seed(1)
  W <- matrix(0, 10, 5)
  W[sample(50, 10)] <- rnorm(10)
  rw <- rewire_matrix(
    W,
    row_names = paste0("g", 1:10),
    col_names = paste0("tf", 1:5),
    row_type  = "gene",
    col_type  = "tf",
    direction = "col_to_row",
    q_del = 0.2, q_add = 0.1, q_flip = 0.1
  )
  expect_true(is.matrix(rw$matrix))
  expect_equal(dim(rw$matrix), dim(W))
  expect_true(is.data.frame(rw$events))
  if (nrow(rw$events) > 0) {
    expect_true("event" %in% names(rw$events))
  }
})

test_that("simulate_dataset returns expected output structure", {
  cfg <- default_config()
  cfg$n_genes              <- 30
  cfg$n_tfs                <- 5
  cfg$n_regions            <- 20
  cfg$n_cells_per_type     <- c(TypeA = 20L, TypeB = 20L)
  cfg$n_reps_per_condition <- 1L
  cfg$chunk_size           <- 50L

  sim <- simulate_dataset(cfg)

  total_cells <- sum(cfg$n_cells_per_type) *
    cfg$n_reps_per_condition *
    length(cfg$conditions)

  expect_equal(nrow(sim$RNA),      total_cells)
  expect_equal(ncol(sim$RNA),      cfg$n_genes)
  expect_equal(nrow(sim$ATAC),     total_cells)
  expect_equal(ncol(sim$ATAC),     cfg$n_regions)
  expect_equal(nrow(sim$metadata), total_cells)
  expect_true(all(sim$ATAC@x %in% c(0, 1)))  # binarized
})
