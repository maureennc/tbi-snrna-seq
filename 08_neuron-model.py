#!/usr/bin/env python3

import os
import pandas as pd
import scanpy as sc
import scvi
from scipy.sparse import csr_matrix

###############################################################################
# Set Paths

data_dir = "<DATA_DIRECTORY_PATH>"
scvi_dir = "<SCVI_MODEL_DIRECTORY>"
save_dir = "<SAVE_DIRECTORY_PATH>"

###############################################################################
# Load Data

adata = sc.read_h5ad(os.path.join(data_dir, '2-tbi-annotated-full.h5ad'))

###############################################################################
# Prepare Neuronal Subset

groups = ['Excitatory neuron', 'Inhibitory neuron']
bdata = adata[adata.obs['cell_type'].isin(groups)].copy()
sc.pp.highly_variable_genes(bdata, n_top_genes=3000, subset=True)

###############################################################################
# Set Up and Train Neuronal Model

scvi.model.SCVI.setup_anndata(
    bdata,
    layer='counts',
    categorical_covariate_keys=['group'],
    continuous_covariate_keys=['total_counts', 'pct_counts_mt', 'pct_counts_ribosomal'],
)

neuronal_model = scvi.model.SCVI(bdata)
neuronal_model.train()

# Save trained model
neuronal_model_dir = os.path.join(scvi_dir, 'neuronal_model')
neuronal_model.save(neuronal_model_dir)

###############################################################################
# Latent Representation and Clustering

latent = neuronal_model.get_latent_representation()
bdata.obsm["X_scVI"] = latent

# Neighbors, UMAP, and Leiden clustering based on scVI latent representation
sc.pp.neighbors(bdata, use_rep="X_scVI", random_state=0)
sc.tl.umap(bdata)
sc.tl.leiden(bdata, key_added='leiden_scVI2', resolution=1)

# Define cluster annotations for neuron types
neuron_cluster = {
    '0': 'Excitatory 1',
    '1': 'Excitatory 2',
    '2': 'Excitatory 3',
    '3': 'Inhibitory 1',
    '4': 'Excitatory 4',
    '5': 'Mixed Ex/Inhib',
    '6': 'Inhibitory 2',
    '7': 'Inhibitory 3',
    '8': 'Excitatory 5',
    '9': 'Excitatory 6',
    '10': 'Inhibitory 4',
    '11': 'Excitatory 7',
    '12': 'Excitatory 8',
    '13': 'Excitatory 9',
    '14': 'Inhibitory 5',
    '15': 'Excitatory 10',
    '16': 'Inhibitory 6',
    '17': 'Inhibitory 7',
    '18': 'Inhibitory 8',
    '19': 'Excitatory 11',
    '20': 'Excitatory 12',
    '21': 'Excitatory 13',
    '22': 'Excitatory 14',
    '23': 'Inhibitory 9',
    '24': 'Excitatory 15',
    '25': 'Excitatory 16',
    '26': 'Excitatory 17',
    '27': 'Excitatory 18',
    '28': 'Excitatory 19',
    '29': 'Inhibitory 10'
}

bdata.obs['neuron_cluster'] = bdata.obs['leiden_scVI2'].map(neuron_cluster)

###############################################################################
# Metadata Mapping to Full Dataset

common_barcodes = bdata.obs_names.intersection(adata.obs_names)
adata = adata[common_barcodes].copy()

# Transfer metadata
unique_columns = [col for col in bdata.obs.columns if col not in adata.obs.columns]
adata.obs = adata.obs.join(bdata.obs[unique_columns], how='left')

# Transfer other data entries
for key in ['uns', 'obsm', 'obsp']:
    for k in getattr(bdata, key).keys():
        setattr(adata, key, {**getattr(adata, key), k: getattr(bdata, key)[k]})

###############################################################################

# Export Processed Data

# Convert to sparse matrix format
bdata.X = csr_matrix(bdata.X)
adata.X = csr_matrix(adata.X)

# Save processed data
bdata.write_h5ad(os.path.join(save_dir, '3-neuron-recluster-hvg.h5ad'))
adata.write_h5ad(os.path.join(save_dir, '3-neuron-recluster-full.h5ad'))

###############################################################################