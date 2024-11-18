#!/usr/bin/env python3

# Import libraries

import os
import scanpy as sc
from scipy.sparse import csr_matrix
import scrublet as scr

###############################################################################

# Define data and save paths

data_dir = "<DATA_DIRECTORY_PATH>"
save_dir = "<SAVE_DIRECTORY_PATH>"

###############################################################################

# Load Data

sample_A = sc.read_10x_mtx(os.path.join(data_dir, "1/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)
sample_B = sc.read_10x_mtx(os.path.join(data_dir, "3/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)
sample_C = sc.read_10x_mtx(os.path.join(data_dir, "5/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)
sample_D = sc.read_10x_mtx(os.path.join(data_dir, "6/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)
sample_E = sc.read_10x_mtx(os.path.join(data_dir, "7/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)
sample_F = sc.read_10x_mtx(os.path.join(data_dir, "8/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)

adata_list = [sample_A, sample_B, sample_C, sample_D, sample_E, sample_F]

###############################################################################

# Quality Control (QC) - Calculate Mitochondrial and Ribosomal Gene Metrics

for adata in adata_list:
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    adata.var['ribosomal'] = adata.var_names.str.match('^(Rpl|Rps)\\d+')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribosomal'], percent_top=None, log1p=False, inplace=True)


for i, adata in enumerate(adata_list):
    adata = adata[adata.obs['total_counts'] > 1000, :]
    adata = adata[adata.obs['total_counts'] < adata.obs['total_counts'].quantile(0.95), :]
    adata = adata[adata.obs['n_genes_by_counts'] > 500, :]
    adata = adata[adata.obs['n_genes_by_counts'] < adata.obs['n_genes_by_counts'].quantile(0.97), :]
    adata = adata[adata.obs['pct_counts_ribosomal'] < adata.obs['pct_counts_ribosomal'].quantile(0.99), :]
    adata = adata[adata.obs['pct_counts_mt'] < 1, :]
    adata_list[i] = adata  # Update filtered sample in list


###############################################################################

# Doublet Scoring

def run_scrublet(adata_list):
    for i, adata in enumerate(adata_list):
        scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.10, random_state=0)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=50)
        adata.obs['doublet_scores'] = doublet_scores
        adata.obs['predicted_doublets'] = predicted_doublets
        adata_list[i] = adata[adata.obs['predicted_doublets'] == False, :] 

run_scrublet(adata_list)

###############################################################################

# Concatenate Filtered Data

adata_concat = sample_A.concatenate(sample_B, sample_C, sample_D, sample_E, sample_F, batch_key='group', batch_categories=['A', 'B', 'C', 'D', 'E', 'F'])

###############################################################################

# Final Pre-Processing

sc.pp.filter_genes(adata_concat, min_cells=3)

adata_concat = adata_concat[adata_concat.obs['doublet_scores'] < 0.1].copy()  # Custom doublet score filter

## Save raw counts, normalize, and log-transform
adata_concat.layers['counts'] = adata_concat.X.copy()
sc.pp.normalize_total(adata_concat)
adata_concat.layers['normalized'] = adata_concat.X.copy()
sc.pp.log1p(adata_concat)
adata_concat.layers['log1p'] = adata_concat.X.copy()
adata_concat.raw = adata_concat

## Copy full gene set data
bdata = adata_concat.copy()

## Select highly variable genes (HVGs) and scale
sc.pp.highly_variable_genes(adata_concat, n_top_genes=3000, subset=True, layer='counts', flavor="seurat_v3")
sc.pp.scale(adata_concat)
adata_concat.layers['scaled'] = adata_concat.X.copy()

###############################################################################

## Export Processed Data

adata_concat.X = csr_matrix(adata_concat.X)
bdata.X = csr_matrix(bdata.X)

adata_concat.write_h5ad(os.path.join(save_dir, '1-tbi-seq-hvg.h5ad'))
bdata.write_h5ad(os.path.join(save_dir, '1-tbi-seq-full.h5ad'))

###############################################################################