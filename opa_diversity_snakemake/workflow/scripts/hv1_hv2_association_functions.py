import pandas as pd
import numpy as np
from munkres import Munkres, print_matrix

def parse_opa_name(name):
    split = name.split('_opa_')
    strain = split[0]
    opa_gene_num = split[1]
    return strain, opa_gene_num

def maximize_trace(pivoted):
    # Maximizes the trace using the munkres algorithm (also works for non-square matrices)
    m = Munkres()
    matrix = (-1)*pivoted.values
    matrix = matrix.tolist()
    indexes = m.compute(matrix)
    total = 0
    column_order = []
    for row, column in indexes:
        value = matrix[row][column]
        total += value
        column_order.append(column)
    total = (-1)*total
    
    new_cols = pivoted.columns[column_order]
    cols_to_add = np.array(list(set(pivoted.columns) - set(new_cols)))
    new_cols_all = np.concatenate([new_cols, cols_to_add])
    pivoted_rearranged = pivoted[new_cols_all]
    
    return pivoted_rearranged, total

def clusters_to_counts(clusters, locus_1_name, locus_2_name):
    # Get counts of each HV1/HV2 combination from the clustering output
    hv1_cluster = []
    hv2_cluster = []
    counts = []
    for i, df in clusters.groupby([locus_1_name, locus_2_name]):
        hv1_cluster.append(i[0])
        hv2_cluster.append(i[1])
        counts.append(len(df))
    counts_df = pd.DataFrame({locus_1_name:hv1_cluster, locus_2_name:hv2_cluster, 'counts':counts})
    return counts_df

def counts_to_matrix(counts_df, locus_1_name, locus_2_name):
    # Convert counts of each HV1/HV2 combination into a pivoted dataframe
    pivoted = counts_df.pivot(index=locus_1_name, columns=locus_2_name, values="counts")
    
    # Note that the munkres algorithm doesn't deal with NaN values, so filling with 0.
    pivoted.fillna(0, inplace = True)
    return pivoted

def clusters_to_trace(clusters, locus_1_name, locus_2_name):
    # Wrapper function that gets the trace from the clusters information dataframe
    counts_df = clusters_to_counts(clusters, locus_1_name, locus_2_name)
    pivoted = counts_to_matrix(counts_df, locus_1_name, locus_2_name)
    pivoted_rearranged, trace = maximize_trace(pivoted)
    return pivoted_rearranged, trace

def get_randomized_matrix_traces(clusters, locus_1_name, locus_2_name, num_trials = 100):
    # Randomize HV2 clusters for each opa gene and get the traces
    clusters_randomized = clusters.copy()

    trace_randomized_all_trials = []
    for i in range(num_trials):
        clusters_randomized[locus_2_name] = clusters[locus_2_name].sample(frac = 1, random_state = i).values
        pivoted_rearranged_randomized, trace_randomized = clusters_to_trace(clusters_randomized, locus_1_name, locus_2_name)
        trace_randomized_all_trials.append(trace_randomized)
    trace_randomized_all_trials = np.array(trace_randomized_all_trials)
    return trace_randomized_all_trials

def get_strains_random_order(baps_clusters):
    # Choose one strain per baps cluster, randomize order
    strains_random_order = []
    
    # Choose one strain per baps cluster
    for cluster, df in baps_clusters.groupby('fastbaps'):
        strains_random_order.append(np.random.choice(df['id']))
        
    # Randomize order
    strains_random_order = np.random.choice(strains_random_order, size=len(strains_random_order), replace=False)
    return strains_random_order