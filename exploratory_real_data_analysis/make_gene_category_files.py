import numpy as np
import os
import sys
import pdb



def create_gene_categories(per_gene_variance_file, output_stem, num_categories=5, statistic='mean'):
    f = open(per_gene_variance_file)
    arr = []
    head_count = 0
    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        if head_count == 0:
            head_count = head_count + 1
            continue
        # Columns: gene_name, var_E0, var_E1, mean_E0, mean_E1
        if statistic == 'mean':
            stat_value = (float(data[3]) + float(data[4]))/2.0
        elif statistic == 'abs_var_diff':
            stat_value = np.abs(np.log2(float(data[1])/float(data[2])))
        else:
            print('assumption error: unsupported statistic ' + str(statistic))
            pdb.set_trace()
        arr.append((data[0], stat_value))
    f.close()

    # Sort genes by the statistic (ascending) and split into num_categories
    # equal-sized chunks (as equal as possible when the count is not divisible)
    arr.sort(key=lambda pair: pair[1])
    gene_names = np.asarray([pair[0] for pair in arr])
    stat_values = np.asarray([pair[1] for pair in arr])
    category_gene_chunks = np.array_split(gene_names, num_categories)
    category_stat_chunks = np.array_split(stat_values, num_categories)

    # Write one file per category listing its gene names (one per line), lowest
    # statistic in category 0 through highest in category num_categories - 1
    summary = open(output_stem + 'summary.txt', 'w')
    summary.write('category\tn_genes\tmin_statistic\tmax_statistic\n')
    for category_index in range(num_categories):
        category_genes = category_gene_chunks[category_index]
        category_stats = category_stat_chunks[category_index]

        t = open(output_stem + 'category_' + str(category_index) + '.txt', 'w')
        for gene_name in category_genes:
            t.write(gene_name + '\n')
        t.close()

        summary.write(str(category_index) + '\t' + str(len(category_genes)) + '\t' +
                      str(category_stats[0]) + '\t' + str(category_stats[-1]) + '\n')
    summary.close()



######################
# command line args
######################
per_gene_variance_file = sys.argv[1]
gene_categories_output_stem = sys.argv[2]

statistic='mean'
create_gene_categories(per_gene_variance_file, gene_categories_output_stem + '_' + statistic + '_based_categories_', statistic=statistic)

statistic='abs_var_diff'
create_gene_categories(per_gene_variance_file, gene_categories_output_stem + '_' + statistic + '_based_categories_', statistic=statistic)
