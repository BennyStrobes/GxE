#!/bin/bash
#SBATCH -t 0-1:00                         # Runtime in D-HH:MM format
#SBATCH -p bch-compute                        # Partition to run in
#SBATCH --mem=5GB


tissue_read_count_file="${1}"
downsampling_percentage="${2}"
downsampling_output_root="${3}"

source ~/.bashrc
conda activate plink_env


echo $downsampling_percentage

python run_downsampling_expression_quantification.py ${tissue_read_count_file} ${downsampling_percentage} ${downsampling_output_root}

