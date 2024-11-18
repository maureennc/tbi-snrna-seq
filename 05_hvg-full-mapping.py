#!/usr/bin/env python3

# Import  libraries

import os
import scanpy as sc
from scipy.sparse import csr_matrix

###############################################################################

# Define Paths

data_dir = "<DATA_DIRECTORY_PATH>"
save_dir = "<SAVE_DIRECTORY_PATH>"
os.makedirs(save_dir, exist_ok=True)

###############################################################################

# Load Annotated HVG and full genome data

adata = sc.read_h5ad(os.path.join(data_dir, '2-tbi-annotated-hvg.h5ad'))
bdata = sc.read_h5ad(os.path.join(data_dir, '1-tbi-seq-full.h5ad'))

###############################################################################

# Find Common Barcodes Between HVG and Full Datasets

common_barcodes = adata.obs_names.intersection(bdata.obs_names)
bdata = bdata[common_barcodes].copy()

###############################################################################

# Transfer metadata from HVG to full dataset

unique_columns = [col for col in adata.obs.columns if col not in bdata.obs.columns]
bdata.obs = bdata.obs.join(adata.obs[unique_columns], how='left')

for key in adata.uns:
    if key not in bdata.uns:
        bdata.uns[key] = adata.uns[key]

for key in adata.obsm:
    if key not in bdata.obsm:
        bdata.obsm[key] = adata.obsm[key]

for key in adata.obsp:
    if key not in bdata.obsp:
        bdata.obsp[key] = adata.obsp[key]

###############################################################################

# Export Full Annotated Dataset

bdata.X = csr_matrix(bdata.X)
bdata.write_h5ad(os.path.join(save_dir, '2-tbi-annotated-full.h5ad'))

###############################################################################
