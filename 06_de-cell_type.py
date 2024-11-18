#!/usr/bin/env python3

# Import libraries

import os
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
import diffxpy.api as de
import numpy as np

# Set output directory
save_dir = "<DE_RESULTS_DIRECTORY>"

###############################################################################

# Load and prepare data
data_dir = "<DATA_DIRECTORY>"
adata = sc.read_h5ad(os.path.join(data_dir, '2-tbi-annotated-full.h5ad.h5ad'))
adata.X = adata.layers['counts'].copy()  # Use raw counts for DE analysis

###############################################################################

# DE Analysis Function

def run_de_analysis(adata, cell_type, reference_group, comparison_group, min_cells_per_group=40, max_imbalance_ratio=0.9):
    subset_data = adata[(adata.obs['cell_type'] == cell_type) & 
                        (adata.obs['group'].isin([reference_group, comparison_group]))].copy()
    
    # Check cell counts and balance between groups
    group_counts = subset_data.obs['group'].value_counts()
    if any(group_counts < min_cells_per_group) or (max(group_counts) / sum(group_counts) > max_imbalance_ratio):
        return None  # Skip if not meeting criteria
    
    subset_data.X = subset_data.layers['counts'].copy()
    sc.pp.filter_genes(subset_data, min_cells=3)

    if issparse(subset_data.X):
        subset_data.X = subset_data.X.toarray()

    ## Calculate size factors
    size_factors = subset_data.X.sum(axis=1) / np.mean(subset_data.X.sum(axis=1))
    subset_data.obs['size_factors'] = size_factors
    
    ## Run Wald test
    test_result = de.test.wald(
        data=subset_data,
        formula_loc="~ 1 + group",
        factor_loc_totest="group",
        size_factors=size_factors
    )
    result_df = test_result.summary()
    result_df['cell_type'] = cell_type
    result_df['comparison'] = f"{reference_group}{comparison_group}"
    
    return result_df

###############################################################################

# Run DE Analysis

cell_types = ['Excitatory neuron', 'Inhibitory neuron', 'Astrocyte', 'Oligodendrocyte', 'OPC', 'Microglia', 'Choroid-plexus epithelial', 'Fibroblast']
comparisons = [('A', 'C'), ('C', 'E'), ('A', 'B')]

all_results = []
for cell_type in cell_types:
    for ref_group, comp_group in comparisons:
        result = run_de_analysis(adata, cell_type, ref_group, comp_group)
        if result is not None:
            all_results.append(result)
        else:
            print(f"Skipping {cell_type} {ref_group} vs {comp_group} due to low counts or imbalance")

## Combine results
final_results = pd.concat(all_results, ignore_index=True)

###############################################################################

# Filter and Export Results
def filter_de_results(df, max_log2fc=100):
    return df[abs(df['log2fc']) <= max_log2fc]

filtered_results = filter_de_results(final_results)

final_results.to_csv(os.path.join(save_dir, 'wald-test-cell_type-full.csv'), index=False)
filtered_results.to_csv(os.path.join(save_dir, 'wald-test-cell_type-filtered.csv'), index=False)

###############################################################################
