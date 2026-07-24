import numpy as np
import os
import sys
import pdb





####################
# Command line args
####################
input_genotype_dir = sys.argv[1]
output_genotype_dir = sys.argv[2]


# Loop through chromosomes
for chrom_num in range(1,23):
    orig_pgen_file = input_genotype_dir + 'gtex_v9_eqtl_chr' + str(chrom_num) + '.pgen'
    new_pgen_file = output_genotype_dir + 'gtex_v9_eqtl_chr' + str(chrom_num) + '.pgen'
    os.system('cp ' + orig_pgen_file + ' ' + new_pgen_file)

    orig_pvar_file = input_genotype_dir + 'gtex_v9_eqtl_chr' + str(chrom_num) + '.pvar'
    new_pvar_file = output_genotype_dir + 'gtex_v9_eqtl_chr' + str(chrom_num) + '.pvar'
    os.system('cp ' + orig_pvar_file + ' ' + new_pvar_file)

    orig_psam_file = input_genotype_dir + 'gtex_v9_eqtl_chr' + str(chrom_num) + '.psam'
    new_psam_file = output_genotype_dir + 'gtex_v9_eqtl_chr' + str(chrom_num) + '.psam'
    f = open(orig_psam_file)
    t = open(new_psam_file,'w')
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split()
        if head_count == 0:
            head_count =head_count + 1
            t.write(line + '\n')
            continue
        words = data[0].split('-')
        new_name = words[0] + '-' + words[1]
        t.write(new_name + '\t' + new_name + '\n')
    f.close()
    t.close()