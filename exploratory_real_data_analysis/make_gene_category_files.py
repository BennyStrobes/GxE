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
            # Look columns up by name so that appending new columns upstream cannot
            # silently shift what a statistic reads
            column_to_index = {}
            for column_position in range(len(data)):
                column_to_index[data[column_position]] = column_position
            continue
        def column(column_name):
            if column_name not in column_to_index:
                print('assumption error: column ' + column_name + ' missing from ' + per_gene_variance_file)
                pdb.set_trace()
            return float(data[column_to_index[column_name]])

        if statistic == 'mean':
            stat_value = (column('mean_E0') + column('mean_E1'))/2.0
        elif statistic == 'abs_var_diff':
            stat_value = np.abs(np.log2(column('var_E0')/column('var_E1')))
        elif statistic == 'abs_de_t':
            # Unsigned Welch t-statistic for a difference in mean expression between the
            # two E strata; unequal variances, so each stratum contributes its own var/n
            standard_error = np.sqrt(column('var_E0')/column('n_E0') + column('var_E1')/column('n_E1'))
            stat_value = np.abs(column('mean_E1') - column('mean_E0'))/standard_error
        elif statistic == 'pred_abs_var_diff':
            # Same contrast as abs_var_diff, but on each stratum's fitted variance from the
            # mean-variance model rather than its observed variance. The fitted values are
            # already log2, so the log2 ratio is just their difference.
            stat_value = np.abs(column('pred_log2_var_E0') - column('pred_log2_var_E1'))
        else:
            print('assumption error: unsupported statistic ' + str(statistic))
            pdb.set_trace()
        # A gene with a non-finite statistic (eg. a zero variance, or no fitted value) cannot be
        # ranked, so it is dropped rather than allowed to sort arbitrarily into a category
        if np.isfinite(stat_value) == False:
            continue
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
create_gene_categories(per_gene_variance_file, gene_categories_output_stem + '_' + statistic + '_based_categories_', statistic=statistic, num_categories=10)

statistic='abs_var_diff'
create_gene_categories(per_gene_variance_file, gene_categories_output_stem + '_' + statistic + '_based_categories_', statistic=statistic, num_categories=10)

statistic='abs_de_t'
create_gene_categories(per_gene_variance_file, gene_categories_output_stem + '_' + statistic + '_based_categories_', statistic=statistic, num_categories=10)

statistic='pred_abs_var_diff'
create_gene_categories(per_gene_variance_file, gene_categories_output_stem + '_' + statistic + '_based_categories_', statistic=statistic, num_categories=10)
