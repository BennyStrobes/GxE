# PA mixture model

Fits a two-component Bayesian mixture model to eQTL effects estimated separately in two environmental (E) bins, in order to distinguish tests exhibiting proportional amplification (PA: genotype effect scales with the bin-specific expression scale, beta_hat_b ~ N(eta * S_b, se_b^2)) from tests exhibiting a standard interaction (INT: beta_hat_b ~ N(beta + e_b * gamma, se_b^2)). Inference is via collapsed Gibbs sampling with hierarchical Gaussian priors on the per-test effects (eta, beta, gamma), inverse-gamma hyperpriors on their variances, and a Beta hyperprior on the mixture proportion pi (the fraction of PA tests).

Running on your own data is a two-step process:
1. Generate mixture model input from expression, genotype, and a binary E variable (`make_pa_mixture_model_input_from_data.py`)
2. Run mixture model inference (`run_pa_mixture_model_inference.py`)

Requires the `PA-h2` conda environment (numpy, scipy, pgenlib).

## Step 1: Generate mixture model input

```bash
python make_pa_mixture_model_input_from_data.py \
    --expression-bed <expression_bed> \
    --plink2-per-chrom-stem <plink2_per_chrom_stem> \
    --binary-E-interaction-covariate-file <E_var_file> \
    --covariate-file <covariate_file or none> \
    --variant-gene-pairs-file <pairs_file or none> \
    --output-stem <output_stem>
```

### Inputs

These are the same inputs as `PA_h2.py`, plus an optional variant-gene pairs file.

- `--expression-bed`: tab-delimited (optionally gzipped) with a header line. Columns are chromosome (`chr1`, ...), gene TSS, end, gene id, then one column per sample. The header gives sample names for the sample columns.
- `--plink2-per-chrom-stem`: plink2 genotype stem; files `<stem><chrom_num>.pgen/.pvar/.psam` for chrom_num 1-22. The psam may contain additional samples, in any order. Missing genotypes are mean imputed.
- `--binary-E-interaction-covariate-file`: tab-delimited. First line is a label followed by sample names; second line is a label followed by the E variable encoded as 0 / 1. Samples must match the expression bed 1-to-1 (same order). E bin 1 is E == 0; E bin 2 is E == 1.
- `--covariate-file` (default `none`): same format as the E file, with one line per covariate (e.g. genotype PCs, expression PCs/PEER factors). Covariates are regressed out of expression and genotype across all samples.
- `--variant-gene-pairs-file`: which variant-gene pairs to test (one test per pair).
  - A file: tab-delimited with a header; the first two columns are gene id and variant id (matching the pvar ID column). A tensorQTL `cis_qtl_top_assoc` file (filtered to significant genes) can be passed directly.
  - `none` (default): pairs are identified with a cis interaction-QTL scan. For every autosomal gene, each cis variant (within `--cis-radius` of the TSS, default 500kb, as in `PA_h2.py`; MAF >= `--min-maf-per-bin` in each E bin, default 0.05) is tested for the genotype x E term in `y ~ g + E + g*E + covariates`. The top variant per gene is kept, its p-value is Bonferroni corrected across the gene's cis variants, and genes with Benjamini-Hochberg q < `--interaction-fdr` (default 0.05) are tested.

### Outputs

- `<output_stem>_pa_mixture_model_geno_resid_S_input.txt`: mixture model input with S = within-bin residual SD (**recommended**)
- `<output_stem>_pa_mixture_model_input.txt`: mixture model input with S = raw within-bin SD
- `<output_stem>_pa_mixture_model_test_key.txt`: maps `test_idx` to `gene_id` / `variant_id`
- `<output_stem>_interaction_qtl_top_per_gene.txt` (interaction-QTL scan only): top variant per tested gene with its interaction effect, p-value, gene-level Bonferroni p-value, BH q-value, and whether it was selected. Its first two columns can be reused as a pairs file.

Each mixture model input file is tab-delimited with one row per test:

| column | description |
|---|---|
| `test_idx` | integer test identifier (0, 1, 2, ...) |
| `S1`, `S2` | within-bin SD of expression in E bin 1 / bin 2 (residual SD from the within-bin regression, or raw SD) |
| `beta_1`, `beta_2` | within-bin OLS effect of genotype on expression, estimated separately in each bin (with intercept) |
| `beta_1_se`, `beta_2_se` | standard errors of those effects |
| `e1_val`, `e2_val` | standardized E values of bin 1 / bin 2 (constant across tests) |

### Notes

- Genotype is standardized across all samples (both bins). Expression is centered and scaled so that its residual variance from a fit of `y ~ g + E + g*E` across all samples is 1. Scaling expression to total variance 1 instead makes each test's scale depend on its own genetic effect size, which biases the estimated PA fraction (in simulations: downward by ~0.02-0.05 with genotype-residualized S).
- Use the genotype-residualized S input. Raw within-bin SDs absorb bin-specific genetic variance and upwardly bias the estimated PA fraction (in simulations: by ~0.03-0.09).
- Expression normalization matters: PA is a difference in expression scale between bins, so expression should be normalized across all samples (never separately within bins). Transformations that force a common distribution (e.g. inverse normal transform) can compress differences in scale between bins.
- When pairs are selected on interaction significance (interaction-QTL scan, or a significance-filtered pairs file), selection and estimation use the same samples, and pi is the PA fraction among *detectable* interaction QTLs. PA interaction effects (eta * (S2 - S1)) tend to be smaller than INT interaction effects, so this can be much lower than the PA fraction among all tests (in one simulation with 43% PA tests, only 6% of significant genes were PA).

## Step 2: Run inference

```bash
python run_pa_mixture_model_inference.py \
    --mixture-model-input-file <output_stem>_pa_mixture_model_geno_resid_S_input.txt \
    --output-stem <results_stem>
```

Optional arguments (defaults in parentheses):
- `--seed` (1): random seed for the Gibbs sampler
- `--ig-a`, `--ig-b` (1e-5, 1e-5): shape/rate of the inverse-gamma hyperprior on the effect variances (weak/near-flat)
- `--pi-prior-alpha`, `--pi-prior-beta` (0.1, 0.1): Beta hyperprior on the mixture proportion pi

The sampler runs 5000 iterations with a 4000-iteration burn-in (set in the function defaults).

### Outputs

- `<results_stem>_global_parameter_posterior_summary.txt`: posterior mean and 95% credible interval for `pi` (fraction of PA tests), `sig_sq_beta`, `sig_sq_gamma`, and `sig_sq_eta`.
- `<results_stem>_per_test_posterior_summary.txt`: per test, the posterior probability of being a PA test (`PA_posterior_prob`, Rao-Blackwellized average over post-burn-in sweeps) and the conditional posterior mean effects (`posterior_mean_eta` given PA; `posterior_mean_beta`, `posterior_mean_gamma` given INT). Join on `test_idx` with the step 1 test key file to get gene / variant ids.

`run_pa_mixture_model_inference_standard_gibbs.py` is an independent implementation of the same model using a standard (non-collapsed, all-univariate-updates) Gibbs sampler, with identical inputs and outputs; it serves as a correctness cross-check and mixes more slowly.

## Simulation

`run_single_simulation.sh` runs one simulation iteration; `driver_key.sh` submits one job per iteration and then runs `visualize_simulation_results.R` (estimated pi, calibration, and posterior distribution plots).

```bash
sbatch run_single_simulation.sh <simulation_iter> <eqtl_sample_size> <n_tests> <e_variable_ratio> <simulated_data_dir> <simulation_results_dir>
```

For each simulated PA fraction (0.2, 0.4, 0.6, 0.8), it:
1. Simulates data with `simulate_mixture_model_data.py` (one gene per test, each paired with one variant on chr1; no covariates), written in the step 1 input format along with a variant-gene pairs file and an oracle file (true test type and S1/S2 per test).
2. Runs step 1 on the simulated data, then deletes the raw simulated data (re-generatable from the `simulation_iter` seed).
3. Runs step 2 (both inference implementations, on both raw-S and genotype-residualized-S input).

Simulated data therefore passes through exactly the same input-generation code as real data. The per-bin summary statistics are computed in `mixture_model_input_utils.py`.
