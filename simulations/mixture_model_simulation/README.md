# PA mixture model inference

`run_pa_mixture_model_inference.py` fits a two-component Bayesian mixture model to per-test eQTL summary statistics, computed separately in two environmental (E) bins, in order to distinguish tests exhibiting proportional amplification (PA: genotype effect scales with the bin-specific expression scale, beta_hat_b ~ N(eta * S_b, se_b^2)) from tests exhibiting a standard interaction (INT: beta_hat_b ~ N(beta + e_b * gamma, se_b^2)). Inference is via collapsed Gibbs sampling with hierarchical Gaussian priors on the per-test effects (eta, beta, gamma), inverse-gamma hyperpriors on their variances, and a Beta hyperprior on the mixture proportion pi (the fraction of PA tests).

## Pipeline overview

Simulated and real data go through the same two steps:

1. `make_pa_mixture_model_input_from_data.py`: expression bed + plink2 genotype + binary E file (+ covariates), i.e. the same inputs as `PA_h2.py`, plus optional variant-gene pairs -> mixture model input files
2. `run_pa_mixture_model_inference.py` (or `run_pa_mixture_model_inference_standard_gibbs.py`): mixture model input file -> posterior summaries

| file | description |
|---|---|
| `make_pa_mixture_model_input_from_data.py` | generates mixture model input files from data (optionally running a cis interaction-QTL scan to choose variant-gene pairs) |
| `mixture_model_input_utils.py` | per-bin summary statistic computation and input file format (shared code) |
| `run_pa_mixture_model_inference.py` | mixture model inference (collapsed Gibbs) |
| `run_pa_mixture_model_inference_standard_gibbs.py` | mixture model inference (standard Gibbs; cross-check) |
| `simulate_mixture_model_data.py` | simulates data in `PA_h2.py` input format |
| `run_single_simulation.sh` | one simulation iteration: simulate -> make input -> inference, for each simulated PA fraction |
| `driver_key.sh` | submits simulation iterations and runs visualization |
| `visualize_simulation_results.R` | simulation result plots (estimated pi, calibration, posterior distributions) |

## Mixture model input format

A tab-delimited file with a header row and one row per test:

| column | description |
|---|---|
| `test_idx` | integer test identifier (0, 1, 2, ...) |
| `S1`, `S2` | within-bin SD of expression in E bin 1 / bin 2 (raw, or with genotype effect residualized out separately in each bin) |
| `beta_1`, `beta_2` | within-bin OLS effect of (globally standardized) genotype on expression, estimated separately in each bin |
| `beta_1_se`, `beta_2_se` | standard errors of those effects |
| `e1_val`, `e2_val` | standardized E values of bin 1 / bin 2 (constant across tests) |

Notes:
- All quantities are computed from within-bin regressions that include an intercept.
- Genotype has mean 0 and variance 1 globally (across both bins).
- Expression has mean 0 and is scaled so that its residual variance from a global (across both bins) fit of `y ~ g + E + g*E` is 1. Scaling expression to total variance 1 instead makes each test's scale depend on its own genetic effect size, which biases the estimated PA fraction (in simulations: downward by ~0.02-0.05 with genotype-residualized S).
- `S1`/`S2` should be the *genotype-residualized* within-bin SDs (the residual SD from the within-bin regression), not the raw within-bin expression SDs. Raw SDs absorb bin-specific genetic variance and upwardly bias the estimated PA fraction (in simulations: by ~0.03-0.09).

## Generating input from data

```bash
python make_pa_mixture_model_input_from_data.py \
    --expression-bed <expression_bed> \
    --plink2-per-chrom-stem <plink2_per_chrom_stem> \
    --binary-E-interaction-covariate-file <E_var_file> \
    --covariate-file <covariate_file or none> \
    --variant-gene-pairs-file <pairs_file or none> \
    --output-stem <output_stem>
```

`--expression-bed`, `--plink2-per-chrom-stem`, `--binary-E-interaction-covariate-file`, and `--covariate-file` have the same format as in `PA_h2.py`. Expression and E file samples must match 1-to-1 (same order); the plink2 psam may contain additional samples and in a different order. Missing genotypes are mean imputed.

Variant-gene pairs to test:
- `--variant-gene-pairs-file <pairs_file>`: tab-delimited with a header; the first two columns are gene id and variant id (one test per line). A tensorQTL `cis_qtl_top_assoc` file (filtered to significant genes) can be passed directly.
- `--variant-gene-pairs-file none` (default): pairs are identified with a cis interaction-QTL scan. For every autosomal gene, each cis variant (within `--cis-radius` of the TSS, default 500kb, as in `PA_h2.py`; MAF >= `--min-maf-per-bin` in each E bin, default 0.05) is tested for the genotype x E term in `y ~ g + E + g*E + covariates`. The top variant per gene is kept, its p-value is Bonferroni corrected across the gene's cis variants, and genes with Benjamini-Hochberg q < `--interaction-fdr` (default 0.05) are used as pairs. All tested genes are written to `<output_stem>_interaction_qtl_top_per_gene.txt` (its first two columns can be reused as a pairs file).
  - Pairs are selected and then estimated on the same samples, and pi is then the PA fraction among detectable interaction QTLs. PA interaction effects (eta * (S2 - S1)) tend to be smaller than INT interaction effects, so this can be much lower than the PA fraction among all tests (in one simulation with 43% PA tests, only 6% of significant genes were PA).

For each pair, covariates (if provided) are regressed out of expression and genotype across all samples, genotype is standardized across all samples, expression is scaled to unit residual variance (see notes above), and per-bin statistics are computed (bin 1: E == 0, bin 2: E == 1).

Outputs:
- `<output_stem>_pa_mixture_model_input.txt`: S = raw within-bin SD
- `<output_stem>_pa_mixture_model_geno_resid_S_input.txt`: S = within-bin residual SD (recommended)
- `<output_stem>_pa_mixture_model_test_key.txt`: maps `test_idx` to `gene_id` / `variant_id`
- `<output_stem>_interaction_qtl_top_per_gene.txt`: only when running the interaction-QTL scan

## Running inference

```bash
python run_pa_mixture_model_inference.py \
    --mixture-model-input-file <input_file> \
    --output-stem <output_stem>
```

Optional arguments (defaults in parentheses):
- `--seed` (1): random seed for the Gibbs sampler
- `--ig-a`, `--ig-b` (1e-5, 1e-5): shape/rate of the inverse-gamma hyperprior on the effect variances (weak/near-flat)
- `--pi-prior-alpha`, `--pi-prior-beta` (0.1, 0.1): Beta hyperprior on the mixture proportion pi

The sampler runs 5000 iterations with a 4000-iteration burn-in (set in the function defaults).

Outputs:
- `<output_stem>_global_parameter_posterior_summary.txt`: posterior mean and 95% credible interval for `pi`, `sig_sq_beta`, `sig_sq_gamma`, and `sig_sq_eta`.
- `<output_stem>_per_test_posterior_summary.txt`: per test, the posterior probability of being a PA test (`PA_posterior_prob`, Rao-Blackwellized average over post-burn-in sweeps) and the conditional posterior mean effects (`posterior_mean_eta` given PA; `posterior_mean_beta`, `posterior_mean_gamma` given INT). Join to gene / variant ids with the test key file.

`run_pa_mixture_model_inference_standard_gibbs.py` is an independent implementation of the same model using a standard (non-collapsed, all-univariate-updates) Gibbs sampler, with identical inputs and outputs; it serves as a correctness cross-check and mixes more slowly.

## Simulation

```bash
sbatch run_single_simulation.sh <simulation_iter> <eqtl_sample_size> <n_tests> <e_variable_ratio> <simulated_data_dir> <simulation_results_dir>
```

For each simulated PA fraction (0.2, 0.4, 0.6, 0.8), `run_single_simulation.sh`:
1. Runs `simulate_mixture_model_data.py`, which simulates one gene per test (each paired with one variant on chr1) and writes data in the `PA_h2.py` input format (expression bed, plink2 genotype, binary E file; no covariates), a variant-gene pairs file, and an oracle file (true test type and S1/S2 per test).
2. Runs `make_pa_mixture_model_input_from_data.py` on the simulated data (with the simulated pairs file), then deletes the raw simulated data to save disk space (it is re-generatable from the `simulation_iter` seed).
3. Runs both inference implementations on both input files (raw S and genotype-residualized S).

`driver_key.sh` submits one job per simulation iteration and then runs `visualize_simulation_results.R`. Requires the `PA-h2` conda environment (numpy, scipy, pgenlib).
