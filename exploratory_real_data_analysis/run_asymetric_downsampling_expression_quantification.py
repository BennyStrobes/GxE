import numpy as np
import sys
import pdb
import os
import gzip
import preprocess_expression


def downsample_read_counts(tissue_read_count_file, downsampled_read_count_file, downsampling_percentage, indi_to_E_var, seed=0):
    # Fixed seed so the same downsampling percentage reproduces the same counts across runs
    np.random.seed(seed)
    head_count = 0
    f = gzip.open(tissue_read_count_file, 'rt')
    t = gzip.open(downsampled_read_count_file, 'wt')
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if head_count == 0:
            head_count = head_count + 1
            t.write(line + '\n')
            indi_names = np.asarray(data[4:])
            downsampling_boolean = []
            for indi in indi_names:
                if indi not in indi_to_E_var:
                    raise ValueError('Individual ' + indi + ' not found in E_var file')
                E_var = indi_to_E_var[indi]
                if E_var == 1.0:
                    downsampling_boolean.append(True)
                elif E_var == 0.0:
                    downsampling_boolean.append(False)
                else:
                    raise ValueError('E_var value for individual ' + indi + ' is not 0 or 1: ' + str(E_var))
            downsampling_boolean = np.asarray(downsampling_boolean)
            continue
        if downsampling_percentage == 1.0:
            t.write(line + '\n')
        else:
            header = np.asarray(data[:4])
            counts = np.asarray(data[4:]).astype(int)
            # Thin each (gene, individual) count by keeping every read independently
            # with probability downsampling_percentage
            counts[downsampling_boolean] = np.random.binomial(counts[downsampling_boolean], downsampling_percentage)
            t.write('\t'.join(header) + '\t' + '\t'.join(counts.astype(str)) + '\n')
    f.close()
    t.close()
    return


def remove_lowly_expressed_genes(read_count_file, gene_filtered_read_count_file, min_reads=6, min_fraction_of_samples=0.2):
    head_count = 0
    f = gzip.open(read_count_file, 'rt')
    t = gzip.open(gene_filtered_read_count_file, 'wt')
    n_input_genes = 0
    n_output_genes = 0
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if head_count == 0:
            head_count = head_count + 1
            t.write(line + '\n')
            continue
        n_input_genes = n_input_genes + 1
        counts = np.asarray(data[4:]).astype(int)
        # Keep the gene only if at least min_fraction_of_samples of the individuals
        # have min_reads or more; the line is written through unchanged so the output
        # keeps the exact format of the input read count file
        fraction_of_samples_expressed = np.sum(counts >= min_reads)/len(counts)
        if fraction_of_samples_expressed >= min_fraction_of_samples:
            t.write(line + '\n')
            n_output_genes = n_output_genes + 1
    f.close()
    t.close()
    print('Gene expression filter: kept ' + str(n_output_genes) + ' of ' + str(n_input_genes) + ' genes')
    return

def create_indi_id_to_E_var_mapping(E_var_file):
    indi_id_to_E_var = {}
    f = open(E_var_file, 'r')
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if head_count == 0:
            head_count = head_count + 1
            indis = np.asarray(data[1:])
            continue
        elif head_count == 1:
            head_count = head_count + 1
            E_var_values = np.asarray(data[1:]).astype(float)
            for i in range(len(indis)):
                if indis[i] in indi_id_to_E_var:
                    raise ValueError('Duplicate individual ID in E_var file: ' + indis[i])
                indi_id_to_E_var[indis[i]] = E_var_values[i]
            continue
        else:
            raise ValueError('E_var file has more than 2 lines')
    f.close()
    return indi_id_to_E_var

#########################
# Command line args
tissue_read_count_file = sys.argv[1]
downsampling_percentage = float(sys.argv[2])
downsampling_output_root = sys.argv[3]
E_var_file = sys.argv[4]

# First, create mapping from individual to E_var value
indi_id_to_E_var = create_indi_id_to_E_var_mapping(E_var_file)

# All downsampled output files share this stem, which mirrors the full-depth naming
# in preprocess_expression.py with the downsampling percentage inserted
downsampling_output_stem = downsampling_output_root + '.downsample_' + str(downsampling_percentage)

# First randomly downsample read counts per (gene,individual)
downsampled_read_count_file = downsampling_output_stem + '.gene_reads.filtered.txt.gz'
downsample_read_counts(tissue_read_count_file, downsampled_read_count_file, downsampling_percentage, indi_id_to_E_var)

# Second, apply the same normalizations used on the full-depth read counts
downsampled_log_tmm_file = downsampling_output_stem + '.log_tmm.txt.gz'
downsampled_log_tmm_unstandardized_file = downsampling_output_stem + '.log_tmm.unstandardized.txt.gz'
downsampled_int_file = downsampling_output_stem + '.inverse_normal_transform.txt.gz'
preprocess_expression.normalize_expression(downsampled_read_count_file, downsampled_log_tmm_file, downsampled_log_tmm_unstandardized_file, downsampled_int_file)


# First randomly downsample read counts per (gene,individual)
downsampled_gene_filtered_read_count_file = downsampling_output_stem + '.gene_filtered_reads.filtered.txt.gz'
remove_lowly_expressed_genes(downsampled_read_count_file, downsampled_gene_filtered_read_count_file)

# Second, apply the same normalizations used on the full-depth read counts
downsampled_filt_log_tmm_file = downsampling_output_stem + '.filtered_genes.log_tmm.txt.gz'
downsampled_filt_log_tmm_unstandardized_file = downsampling_output_stem + '.filtered_genes.log_tmm.unstandardized.txt.gz'
downsampled_filt_int_file = downsampling_output_stem + '.filtered_genes.inverse_normal_transform.txt.gz'
preprocess_expression.normalize_expression(downsampled_gene_filtered_read_count_file, downsampled_filt_log_tmm_file, downsampled_filt_log_tmm_unstandardized_file, downsampled_filt_int_file)
