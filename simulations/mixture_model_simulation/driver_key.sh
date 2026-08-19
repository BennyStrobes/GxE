

#######################
# output directories
#######################
# Output root directory
output_stem="/lab-share/CHIP-Strober-e2/Public/ben/gxe/simulations/mixture_model_simulation/"

# Output directory containing simulated data
simulated_data_dir=${output_stem}"simulation_data/"

# Output directory containing simulation results
simulation_results_dir=${output_stem}"simulation_results/"

# Output directory containing organized simulation results
organized_simulation_results_dir=${output_stem}"organized_simulation_results/"

# Output directory containing visualizations of simulated data
visualization_results_dir=${output_stem}"visualization_results/"





#######################
# Code
#######################

# Simulation parameters
n_simulations=50
eqtl_sample_size=500
n_tests=2000
e_variable_ratio=0.4

# Loop through simulation iterations (one cluster job per iteration)
for simulation_iter in $(seq 1 ${n_simulations}); do
	sbatch run_single_simulation.sh ${simulation_iter} ${eqtl_sample_size} ${n_tests} ${e_variable_ratio} ${simulated_data_dir} ${simulation_results_dir}
done
