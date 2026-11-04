# scNovoMOSim News

## scNovoMOSim 0.1.0

* Initial release.
* Paired scRNA-seq (NB) and scATAC-seq (Poisson / binary) simulation.
* TF-gene GRN and region-TF accessibility networks with shared backbone,
  cell-type-specific layers and condition-specific rewiring.
* `default_config()` for parameter management.
* `simulate_dataset()` as main entry point with chunked cell generation.
* `export_simulation()` writing Matrix Market + CSV + JSON.
* `quick_qc()` for post-simulation sanity checks.
* Ground-truth rewiring event table for benchmarking.
