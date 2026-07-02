from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import pdb






#####################
# Command line args
#####################
simulation_type = sys.argv[1]
simulated_data_dir = sys.argv[2]
simulation_results_dir = sys.argv[3]
organized_simulation_results_dir = sys.argv[4]

sim_values = np.arange(1,101)
eqtl_sss = ['100', '300']

t = open(organized_simulation_results_dir + 'simulation_summary.txt','w' )
t.write('eqtl_ss\tsim_type\test\test_se\ttruth\n')

for eqtl_ss in eqtl_sss:
	sim_pa_interaction_variance_arr = []
	est_pa_interaction_variance_arr = []
	est_perm_pa_interaction_variance_arr = []
	for sim_iter in sim_values:
		# Extract simulated PA interaction variance for this sim
		sim_data_file = simulated_data_dir + 'simulation_' + str(sim_iter) + '_sim_type_total_PA_default_N_' + eqtl_ss + '.simulation_summary.txt'
		sim_pa_interaction_variance = np.mean((np.loadtxt(sim_data_file,dtype=str,delimiter='\t')[1:,-3]).astype(float))


		# Extract estimated PA interaction variance for this sim
		sim_summary_file = simulation_results_dir + 'simulation_' + str(sim_iter) + '_sim_type_total_PA_default_N_' + eqtl_ss + '_PA_h2_summary.txt'
		est_pa_interaction_variance = float(np.loadtxt(sim_summary_file,dtype=str,delimiter='\t')[3,1])

		# Extract estimated perm PA interaction variance for this sim
		sim_summary_file = simulation_results_dir + 'simulation_' + str(sim_iter) + '_sim_type_total_PA_default_N_' + eqtl_ss + '_perm_PA_h2_summary.txt'
		est_perm_pa_interaction_variance = float(np.loadtxt(sim_summary_file,dtype=str,delimiter='\t')[3,1])

		sim_pa_interaction_variance_arr.append(sim_pa_interaction_variance)
		est_pa_interaction_variance_arr.append(est_pa_interaction_variance)
		est_perm_pa_interaction_variance_arr.append(est_perm_pa_interaction_variance)
	sim_pa_interaction_variance_arr = np.asarray(sim_pa_interaction_variance_arr)
	est_pa_interaction_variance_arr = np.asarray(est_pa_interaction_variance_arr)
	est_perm_pa_interaction_variance_arr = np.asarray(est_perm_pa_interaction_variance_arr)

	t.write(eqtl_ss + '\t' + 'real' + '\t' + str(np.mean(est_pa_interaction_variance_arr)) + '\t' + str(np.std(est_pa_interaction_variance_arr)/np.sqrt(len(est_pa_interaction_variance_arr))) + '\t' + str(np.mean(sim_pa_interaction_variance_arr)) + '\n')
	t.write(eqtl_ss + '\t' + 'perm' + '\t' + str(np.mean(est_perm_pa_interaction_variance_arr)) + '\t' + str(np.std(est_perm_pa_interaction_variance_arr)/np.sqrt(len(est_perm_pa_interaction_variance_arr))) + '\t' + str(0.0) + '\n')

t.close()

print(organized_simulation_results_dir + 'simulation_summary.txt')








def load_simulation_summary(path: str) -> pd.DataFrame:
	"""
	Load simulation_summary.txt (tab-delimited) into a DataFrame with typed columns.
	Expected columns: eqtl_ss, sim_type, est, est_se, truth
	"""
	df = pd.read_csv(path, sep="\t")
	required = {"eqtl_ss", "sim_type", "est", "est_se", "truth"}
	missing = required - set(df.columns)
	if missing:
		raise ValueError(f"Missing columns {missing} in file: {path}")

	# Type coercion (robust to accidental string parsing)
	df["eqtl_ss"] = pd.to_numeric(df["eqtl_ss"], errors="raise").astype(int)
	df["sim_type"] = df["sim_type"].astype(str)
	df["est"] = pd.to_numeric(df["est"], errors="raise")
	df["est_se"] = pd.to_numeric(df["est_se"], errors="raise")
	df["truth"] = pd.to_numeric(df["truth"], errors="raise")
	return df


def plot_estimates_with_ci(
	df: pd.DataFrame,
	truth_line: float = 0.0016818645436700386,
	ci_z: float = 1.96,
	title: str | None = None,
):
	"""
	Bar plot of est by (eqtl_ss, sim_type) with 95% CI from est_se.
	Adds a horizontal dashed line at truth_line.
	"""
	# Ensure deterministic order on x and within-group bars
	ss_vals = sorted(df["eqtl_ss"].unique())
	sim_order = ["real", "perm"]
	# Keep only those present (but preserve preferred order)
	sim_vals = [s for s in sim_order if s in set(df["sim_type"]) ] + \
			   [s for s in sorted(set(df["sim_type"])) if s not in sim_order]

	# Build a full grid so missing combos (if any) don't break plotting
	grid = (
		pd.MultiIndex.from_product([ss_vals, sim_vals], names=["eqtl_ss", "sim_type"])
		.to_frame(index=False)
		.merge(df, on=["eqtl_ss", "sim_type"], how="left")
	)

	# Compute 95% CI
	grid["ci_low"] = grid["est"] - ci_z * grid["est_se"]
	grid["ci_high"] = grid["est"] + ci_z * grid["est_se"]
	grid["yerr"] = ci_z * grid["est_se"]

	# Plot: grouped bars
	n_groups = len(ss_vals)
	n_sims = len(sim_vals)
	x = np.arange(n_groups)

	width = 0.8 / max(n_sims, 1)  # total group width ~0.8
	offsets = (np.arange(n_sims) - (n_sims - 1) / 2.0) * width

	fig, ax = plt.subplots(figsize=(7, 4.5))

	for j, sim in enumerate(sim_vals):
		sub = grid[grid["sim_type"] == sim].copy()
		# Align to ss order
		sub = sub.set_index("eqtl_ss").reindex(ss_vals).reset_index()

		heights = sub["est"].to_numpy()
		yerr = sub["yerr"].to_numpy()

		ax.bar(
			x + offsets[j],
			heights,
			width=width,
			yerr=yerr,
			capsize=4,
			label=sim,
			edgecolor="black",
			linewidth=0.7,
		)

	# Truth line
	ax.axhline(truth_line, linestyle="--", linewidth=1.5, color='0.6')

	ax.set_xticks(x)
	ax.set_xticklabels([str(v) for v in ss_vals])
	ax.set_xlabel("eqtl_ss")
	ax.set_ylabel("est (± 95% CI)")
	if title is None:
		title = "Simulation estimates with 95% CI"
	ax.set_title(title)
	ax.legend(title="sim_type", frameon=False)

	# Nice padding
	ax.margins(x=0.05)

	return fig, ax


# ---- Example usage ----
df = load_simulation_summary(organized_simulation_results_dir + 'simulation_summary.txt')
fig, ax = plot_estimates_with_ci(
    df,
    truth_line=0.0016818645436700386,
    title="Simulation estimates with 95% CI",
)

outpath = organized_simulation_results_dir + 'visualize_simulation_summary.png'
fig.savefig(
    outpath,
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)

print(outpath)










