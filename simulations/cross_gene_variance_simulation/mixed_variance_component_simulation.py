import numpy as np
import os
import sys
import pdb
from pandas_plink import read_plink1_bin




def get_gene_tss(protein_coding_gene_file, gene_tss_file, chrom_num):
	# First extract list of protein coding genes
	pc_genes = {}
	f = open(protein_coding_gene_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		gene_id = data[0]
		pc_genes[data[0]] = 1
	f.close()

	# Second get array of tss's of all PC gene on desired chromosome
	tss_arr = []
	f = open(gene_tss_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		gene_id = data[1]
		if gene_id not in pc_genes:
			continue
		line_chrom = data[7]
		if line_chrom != 'chr' + chrom_num:
			continue
		tss = float(data[8])
		tss_arr.append(tss)
	f.close()

	return np.asarray(tss_arr)

def get_standardized_genotype_matrix(gene_geno):
	new_mat = np.copy(gene_geno)
	NN,KK = gene_geno.shape
	for kk in range(KK):
		new_mat[:,kk] = (new_mat[:,kk] - np.mean(new_mat[:,kk]))/np.std(new_mat[:,kk])
	return new_mat

def simulate_mixed_variance_component_expression(gene_geno, e_variable_proportion, proportional_amplification_scaling_factor, n_proportional_amplication_causal_snps_per_gene, orig_ge_h2, n_constant_causal_snps_per_gene, constant_ge_h2, n_interaction_causal_snps_per_gene, interaction_ge_h2):
	# Remove snps that have no variance across indivudals
	polymorphic_snp_indices = np.var(gene_geno,axis=0) != 0.0
	gene_geno = gene_geno[:, polymorphic_snp_indices]

	# Get sample size and number of snps
	NN,KK = gene_geno.shape

	# Simulate E variable
	EE = np.random.binomial(n=1, p=e_variable_proportion, size=NN).astype(float)

	# Generate original causal effect sizes
	orig_causal_effect_sizes = np.zeros(KK)
	causal_indices = np.random.choice(np.arange(KK), size=n_proportional_amplication_causal_snps_per_gene, replace=False)
	orig_causal_effect_sizes[causal_indices] = np.random.normal(loc=0, scale=np.sqrt(orig_ge_h2/n_proportional_amplication_causal_snps_per_gene), size=n_proportional_amplication_causal_snps_per_gene)

	# Get standardized genotype matrix
	standardized_geno = get_standardized_genotype_matrix(gene_geno)

	# Get variance scaling factor
	SS = np.ones(len(EE))
	SS[EE==1.0] = proportional_amplification_scaling_factor
	SS = SS/np.mean(SS)
	#SS = np.square(np.sqrt(SS)/np.mean(np.sqrt(SS)))

	scaled_geno = standardized_geno * np.sqrt(SS[:, np.newaxis])
	mean_effects = np.dot(scaled_geno, orig_causal_effect_sizes)

	only_genetic_effects = np.dot(standardized_geno, orig_causal_effect_sizes)

	sqrtS = np.sqrt(SS)
	aa = np.mean(sqrtS)
	bb = sqrtS - np.mean(sqrtS)

	G_const = np.mean(sqrtS) * only_genetic_effects

	G_int = (sqrtS - np.mean(sqrtS)) * only_genetic_effects
	# Note that G_const + G_int == mean_effects


	# Constant effects
	constant_causal_effect_sizes = np.zeros(KK)
	constant_causal_indices = np.random.choice(np.arange(KK), size=n_constant_causal_snps_per_gene, replace=False)
	constant_causal_effect_sizes[constant_causal_indices] = np.random.normal(loc=0, scale=np.sqrt(constant_ge_h2/n_constant_causal_snps_per_gene), size=n_constant_causal_snps_per_gene)
	constant_mean_effects = np.dot(standardized_geno, constant_causal_effect_sizes)

	# Interaction effects
	interaction_causal_effect_sizes = np.zeros(KK)
	interaction_causal_indices = np.random.choice(np.arange(KK), size=n_interaction_causal_snps_per_gene, replace=False)
	interaction_causal_effect_sizes[interaction_causal_indices] = np.random.normal(loc=0, scale=np.sqrt(interaction_ge_h2/n_interaction_causal_snps_per_gene), size=n_interaction_causal_snps_per_gene)
	standardized_E = (EE - np.mean(EE))/np.std(EE)
	gxe = standardized_geno * standardized_E[:, np.newaxis]
	interaction_mean_effects = np.dot(gxe, interaction_causal_effect_sizes)



	# Residual variance
	scaled_resid_var = np.copy(SS)*(1.0 - orig_ge_h2 - constant_ge_h2 - interaction_ge_h2)
	residiuals = np.random.normal(loc=0, scale=np.sqrt(scaled_resid_var))


	YY = mean_effects + constant_mean_effects + interaction_mean_effects + residiuals

	return YY, standardized_geno, EE,SS, np.var(G_const,ddof=1), np.var(G_int, ddof=1), np.var(constant_mean_effects, ddof=1), np.var(interaction_mean_effects, ddof=1)

def get_reml_resid_var(YY_sub, GG_sub):
	Y_t_Y = np.dot(YY_sub.reshape(-1,1), YY_sub.reshape(1,-1))
	G_T_G = np.dot(GG_sub, np.transpose(GG_sub))
	rsid_var_cov = np.eye(len(YY_sub))
	i_idx, j_idx = np.triu_indices(len(YY_sub), k=0)

	Y_t_Y_vec = Y_t_Y[i_idx, j_idx]
	G_T_G_vec = G_T_G[i_idx, j_idx]
	resid_vec = rsid_var_cov[i_idx, j_idx]
	input_mat = np.transpose(np.vstack((resid_vec, G_T_G_vec)))
	beta, *_ = np.linalg.lstsq(input_mat, Y_t_Y_vec, rcond=None)

	est = np.max([0.01, beta[0]])

	return est



def variance_partitioning(YY, GG, EE,SS, prop_A=True, scale_estimate_w_reml=True):
	# Standardize EE
	stand_EE = (EE - np.mean(EE))/np.std(EE)
	EE_scaled_geno = GG * stand_EE[:, np.newaxis]
	#EE_scaled_geno = (EE_scaled_geno - np.mean(EE_scaled_geno,axis=0))/np.std(EE_scaled_geno,axis=0)   ##LKNELKRJELKEJRL:EJ



	Y_t_Y = np.dot(YY.reshape(-1,1), YY.reshape(1,-1))
	E_t_E = np.dot(stand_EE.reshape(-1,1), stand_EE.reshape(1,-1))
	G_T_G = np.dot(GG, np.transpose(GG))
	std_interaction_cov = np.dot(EE_scaled_geno, np.transpose(EE_scaled_geno))
	rsid_var_cov = np.eye(len(YY))



	pred_SS = np.ones(len(EE))
	if scale_estimate_w_reml:
		pred_SS[EE==1.0] = get_reml_resid_var(YY[EE==1.0], GG[EE==1.0,:])
		pred_SS[EE==0.0] = get_reml_resid_var(YY[EE==0.0], GG[EE==0.0,:])
	else:
		pred_SS[EE==1.0] = np.var(YY[EE==1.0],ddof=1)
		pred_SS[EE==0.0] = np.var(YY[EE==0.0],ddof=1)		
	pred_SS = pred_SS/np.mean(pred_SS)
	sqrt_S_scaled_geno = GG * np.sqrt(pred_SS[:, np.newaxis])
	#pred_SS = SS/np.mean(SS)
	#sqrt_S_scaled_geno = GG * np.sqrt(SS[:, np.newaxis])

	prop_A_cov =np.dot(sqrt_S_scaled_geno, np.transpose(sqrt_S_scaled_geno))


	i_idx, j_idx = np.triu_indices(len(YY), k=1)

	Y_t_Y_vec = Y_t_Y[i_idx, j_idx]
	E_t_E_vec = E_t_E[i_idx, j_idx]
	G_T_G_vec = G_T_G[i_idx, j_idx]
	std_interaction_cov_vec = std_interaction_cov[i_idx, j_idx]
	prop_A_cov_vec = prop_A_cov[i_idx, j_idx]

	if prop_A:
		input_mat = np.transpose(np.vstack((E_t_E_vec, G_T_G_vec, std_interaction_cov_vec, prop_A_cov_vec)))
		beta, *_ = np.linalg.lstsq(input_mat, Y_t_Y_vec, rcond=None)
		variance_parameters = beta[1:]*GG.shape[1]
		h2_est = variance_parameters[0]
		gxe_est = variance_parameters[1]
		pa_est = variance_parameters[2]
		ratio = np.square(np.mean(np.sqrt(pred_SS)))
		pa_mean_est = pa_est*ratio
		pa_interaction_est = pa_est*(1.0-ratio)
	else:
		input_mat = np.transpose(np.vstack((E_t_E_vec, G_T_G_vec, std_interaction_cov_vec)))
		beta, *_ = np.linalg.lstsq(input_mat, Y_t_Y_vec, rcond=None)
		variance_parameters = beta[1:]*GG.shape[1]
		h2_est = variance_parameters[0]
		gxe_est = variance_parameters[1]
		pa_est = np.nan
		#ratio = np.square(np.mean(np.sqrt(pred_SS)))
		pa_mean_est = np.nan
		pa_interaction_est = np.nan

	return h2_est, gxe_est, pa_est, pa_mean_est, pa_interaction_est



def partition_variance_constant_vs_interaction(EE, SS, only_genetic_mean_effects):
	# Center environment
	E_c = EE - np.mean(EE)

	# Fit sqrt(S) = a + b * E_c by least squares
	sqrtS = np.sqrt(SS)
	X = np.vstack([np.ones_like(E_c), E_c]).T  # [1, E_c]
	a, b = np.linalg.lstsq(X, sqrtS, rcond=None)[0]

	Z = only_genetic_mean_effects

	G_const = a * Z
	G_int   = b * E_c * Z

	var_const = np.var(G_const, ddof=1)
	var_int   = np.var(G_int, ddof=1)
	cov_ci    = np.cov(G_const, G_int, ddof=1)[0, 1]
	var_total_genetic = np.var(G_const + G_int, ddof=1)

	return var_const, var_int


def organize_simulation_results_across_simulation_runs(raw_output_file, organized_output_file, inference_methods, variance_components):
	# open up new output file hanlde
	t = open(organized_output_file,'w')
	t.write('inference_method\tvariance_component\testimated_variance\testimated_variance_lb\testimated_variance_ub\ttrue_simulated_variance\ttrue_simulated_variance_lb\ttrue_simulated_variance_ub\n')

	for inference_method in inference_methods:
		for variance_component in variance_components:
			estimated_variances = []
			simulated_variances = []

			f = open(raw_output_file)
			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue

				# Only conside lines corresponding to the current inference_method-variance_component pair
				if data[1] != inference_method or data[2] != variance_component:
					continue

				# Extract estimates and truth for this line
				estimated_variance = float(data[3])
				simulated_variance = float(data[4])

				# Add to global array
				estimated_variances.append(estimated_variance)
				simulated_variances.append(simulated_variance)

			# Convert to numpy arrays
			estimated_variances = np.asarray(estimated_variances)
			simulated_variances = np.asarray(simulated_variances)


			# Stats for estimated variance components
			estimated_mean = np.mean(estimated_variances)
			estimated_mean_se = np.std(estimated_variances)/np.sqrt(len(estimated_variances))
			estimated_mean_lb = estimated_mean - 1.96*estimated_mean_se
			estimated_mean_ub = estimated_mean + 1.96*estimated_mean_se

			# Stats for simulated variance components
			simulated_mean = np.mean(simulated_variances)
			simulated_mean_se = np.std(simulated_variances)/np.sqrt(len(simulated_variances))
			simulated_mean_lb = simulated_mean - 1.96*estimated_mean_se
			simulated_mean_ub = simulated_mean + 1.96*estimated_mean_se

			# Print to output file
			t.write(inference_method + '\t' + variance_component + '\t' + str(estimated_mean) + '\t' + str(estimated_mean_lb) + '\t' + str(estimated_mean_ub) + '\t' + str(simulated_mean) + '\t' + str(simulated_mean_lb) + '\t' + str(simulated_mean_ub) + '\n')

	t.close()
	return

	
#################
# Command line args
#################
plink_genotype_stem = sys.argv[1]
chrom_num = sys.argv[2]
e_variable_proportion = float(sys.argv[3])
proportional_amplification_scaling_factor = float(sys.argv[4])
n_proportional_amplication_causal_snps_per_gene = int(sys.argv[5])
orig_ge_h2 = float(sys.argv[6])
sample_size = int(sys.argv[7])
cis_radius = int(sys.argv[8])
simulation_results_dir = sys.argv[9]
gene_tss_file = sys.argv[10]
protein_coding_gene_file = sys.argv[11]
n_constant_causal_snps_per_gene = int(sys.argv[12])
constant_ge_h2 = float(sys.argv[13])
n_interaction_causal_snps_per_gene = int(sys.argv[14])
interaction_ge_h2 = float(sys.argv[15])


#######################
# Simulation settings
#######################
n_sims = 50000
min_snps_per_gene = 50

#######################
# Get list of gene's tss to sample over
#######################
gene_tss_arr = get_gene_tss(protein_coding_gene_file, gene_tss_file, chrom_num)

#######################
# Load in genotype data
#######################
G_obj = read_plink1_bin(plink_genotype_stem + chrom_num + ".bed", plink_genotype_stem + chrom_num +".bim", plink_genotype_stem + chrom_num +".fam", verbose=False)
G_obj_geno = G_obj.values # Numpy 2d array of dimension num samples X num snps
G_obj_chrom = np.asarray(G_obj.chrom)
G_obj_pos = np.asarray(G_obj.pos)

#######################
# Open output file handle
#######################

raw_output_file = simulation_results_dir + 'mixed_variance_component_simulation_propAmpScale_' + str(proportional_amplification_scaling_factor) + '_sampleSize_' + str(sample_size) + '_raw_simulation_results2.txt'

t = open(raw_output_file,'w')
t.write('simulation_number\tinference_method\tvariance_component\testimated_variance\ttrue_simulated_variance\n')

#######################
# Loop through simulations
#######################
for sim_iter in range(n_sims):
	# Make sure to select a gene with at least min_snps_per_gene snps
	snp_count = 0
	while snp_count < min_snps_per_gene:
		# Select gene for this simulation
		gene_tss = np.random.choice(gene_tss_arr)
		# Filter genotype data for this simulation
		snp_indices = (G_obj_pos >= gene_tss - cis_radius) & (G_obj_pos < gene_tss + cis_radius)
		# Get number of snps passing filter
		snp_count = np.sum(snp_indices)


	# Filter genotype matrix to snps for this gene
	gene_geno = G_obj_geno[:, snp_indices]

	# Filter genotype matrix to correct number of individuals
	individual_indices = np.sort(np.random.choice(np.arange(gene_geno.shape[0]), size=sample_size, replace=False))
	gene_geno = gene_geno[individual_indices, :]

	# Simulate data
	YY, GG, EE,SS, pa_constant_var, pa_interaction_var, constant_var, interaction_var = simulate_mixed_variance_component_expression(gene_geno, e_variable_proportion, proportional_amplification_scaling_factor, n_proportional_amplication_causal_snps_per_gene, orig_ge_h2, n_constant_causal_snps_per_gene, constant_ge_h2, n_interaction_causal_snps_per_gene, interaction_ge_h2)

	# Run inference
	h2_est, gxe_est, pa_est, pa_mean_est, pa_interaction_est = variance_partitioning(YY, GG, EE,SS)

	# Print results to output file
	t.write(str(sim_iter) + '\t' + 'PA+interaction-VCA' + '\t' + 'constant_genetic_variance' + '\t' + str(h2_est + pa_mean_est) + '\t' + str(pa_constant_var + constant_var) + '\n')
	t.write(str(sim_iter) + '\t' + 'PA+interaction-VCA' + '\t' + 'interaction_genetic_variance' + '\t' + str(gxe_est) + '\t' + str(interaction_var) + '\n')
	t.write(str(sim_iter) + '\t' + 'PA+interaction-VCA' + '\t' + 'PA_interaction_genetic_variance' + '\t' + str(pa_interaction_est) + '\t' + str(pa_interaction_var) + '\n')
	t.write(str(sim_iter) + '\t' + 'PA+interaction-VCA' + '\t' + 'total_genetic_variance' + '\t' + str(h2_est + gxe_est + pa_mean_est + pa_interaction_est) + '\t' + str(constant_var + interaction_var + pa_constant_var + pa_interaction_var) + '\n')

	# Run inference without proportional amplification variance components
	no_pa_mod_h2_est, no_pa_mod_gxe_est, no_pa_mod_pa_est, no_pa_mod_pa_mean_est, no_pa_mod_pa_interaction_est = variance_partitioning(YY, GG, EE,SS, False)

	# Print results to output file
	t.write(str(sim_iter) + '\t' + 'interaction-VCA' + '\t' + 'constant_genetic_variance' + '\t' + str(no_pa_mod_h2_est) + '\t' + str(pa_constant_var + constant_var) + '\n')
	t.write(str(sim_iter) + '\t' + 'interaction-VCA' + '\t' + 'interaction_genetic_variance' + '\t' + str(no_pa_mod_gxe_est) + '\t' + str(interaction_var) + '\n')
	t.write(str(sim_iter) + '\t' + 'interaction-VCA' + '\t' + 'PA_interaction_genetic_variance' + '\t' + str(0.0) + '\t' + str(pa_interaction_var) + '\n')
	t.write(str(sim_iter) + '\t' + 'interaction-VCA' + '\t' + 'total_genetic_variance' + '\t' + str(no_pa_mod_h2_est + no_pa_mod_gxe_est) + '\t' + str(constant_var + interaction_var + pa_constant_var + pa_interaction_var) + '\n')
t.close()

################################
# Organize simulation results
# Aggregate across simulation iters
################################
organized_output_file = simulation_results_dir + 'mixed_variance_component_simulation_propAmpScale_' + str(proportional_amplification_scaling_factor) + '_sampleSize_' + str(sample_size) + '_organized_simulation_results2.txt'
inference_methods = np.asarray(['interaction-VCA', 'PA+interaction-VCA'])
variance_components = np.asarray(['constant_genetic_variance', 'interaction_genetic_variance', 'PA_interaction_genetic_variance', 'total_genetic_variance'])
organize_simulation_results_across_simulation_runs(raw_output_file, organized_output_file, inference_methods, variance_components)
print(organized_output_file)




