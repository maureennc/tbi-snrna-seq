#!/usr/bin/env python3

# Import libraries

import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
import torch
import random
import scvi

###############################################################################

# Settings

## Random seeds 
random.seed(0)
torch.manual_seed(0)
np.random.seed(0)
scvi.settings.seed = 0

# Define data and save paths
data_dir = "<DATA_DIRECTORY_PATH>"
scvi_dir = "<MODEL_DIRECTORY_PATH>"
save_dir = "<SAVE_DIRECTORY_PATH>"
os.makedirs(save_dir, exist_ok=True)

###############################################################################

# Load Data

adata = sc.read_h5ad(os.path.join(data_dir, '1-tbi-seq-hvg.h5ad'))

model = scvi.model.SCVI.load(scvi_dir, adata=adata)

###############################################################################

# Clustering and Annotation

## Leiden clustering using pre-trained scVI latent rep
sc.pp.neighbors(adata, use_rep='X_scVI', random_state=0)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added='leiden_scVI', resolution=1)


## Broad cell class annotations
cell_class = { 
    "0": "Astrocyte", 
    "1": "Neuron", 
    "2": "Neuron", 
    "3": "Oligodendrocyte", 
    "4": "Neuron",
    "5": "Astrocyte", 
    "6": "Choroid-plexus epithelial", 
    "7": "Neuron", 
    "8": "Neuron", 
    "9": "Neuron",
    "10": "OPC", 
    "11": "Neuron", 
    "12": "Neuron", 
    "13": "Neuron", 
    "14": "Neuron",
    "15": "Neuron", 
    "16": "Neuron", 
    "17": "Neuron", 
    "18": "Neuron", 
    "19": "Fibroblast",
    "20": "Microglia", 
    "21": "Neuron", 
    "22": "Neuron", 
    "23": "Neuron", 
    "24": "Neuron",
    "25": "Neuron", 
    "26": "Astrocyte", 
    "27": "Neuron", 
    "28": "Artifact",
    "29": "Neuron", 
    "30": "Artifact", 
    "31": "Neuron", 
    "32": "Unknown",
    "33": "Neuron", 
    "34": "Neuron",
    "35": "Artifact", 
    "36": "Artifact"
}
adata.obs['cell_class'] = adata.obs.leiden_scVI.map(cell_class)

## Specific cell type annotations
cell_type = { 
    "0": "Astrocyte", 
    "1": "Excitatory neuron", 
    "2": "Excitatory neuron", 
    "3": "Oligodendrocyte",
    "4": "Excitatory neuron", 
    "5": "Astrocyte", 
    "6": "Choroid-plexus epithelial",  
    "7": "Inhibitory neuron",
    "8": "Inhibitory neuron", 
    "9": "Excitatory neuron", 
    "10": "OPC", 
    "11": "Inhibitory neuron",
    "12": "Mixed Ex/In", 
    "13": "Excitatory neuron", 
    "14": "Excitatory neuron",
    "15": "Excitatory neuron", 
    "16": "Excitatory neuron", 
    "17": "Excitatory neuron",
    "18": "Inhibitory neuron", 
    "19": "Fibroblast", 
    "20": "Microglia", 
    "21": "Excitatory neuron",
    "22": "Inhibitory neuron", 
    "23": "Inhibitory neuron", 
    "24": "Excitatory neuron",
    "25": "Excitatory neuron", 
    "26": "Astrocyte", 
    "27": "Excitatory neuron", 
    "28": "Artifact",  
    "29": "Inhibitory neuron", 
    "30": "Artifact",  
    "31": "Excitatory neuron", 
    "32": "Unknown", 
    "33": "Excitatory neuron",
    "34": "Excitatory neuron", 
    "35": "Artifact",  
    "36": "Artifact" 
}
adata.obs['cell_type'] = adata.obs.leiden_scVI.map(cell_type)

## Detailed cluster annotations
cluster = { 
    "0": "Astrocyte 1", 
    "1": "Excitatory neuron 1", 
    "2": "Excitatory neuron 2", 
    "3": "Oligodendrocyte",
    "4": "Excitatory neuron 3", 
    "5": "Astrocyte 2", 
    "6": "Choroid-plexus epithelial",  
    "7": "Inhibitory neuron 1", 
    "8": "Inhibitory neuron 2",
    "9": "Excitatory neuron 4", 
    "10": "OPC", 
    "11": "Inhibitory neuron 3",
    "12": "Mixed Ex/In", 
    "13": "Excitatory neuron 5", 
    "14": "Excitatory neuron 6",
    "15": "Excitatory neuron 7", 
    "16": "Excitatory neuron 8", 
    "17": "Excitatory neuron 9",
    "18": "Inhibitory neuron 5", 
    "19": "Fibroblast", 
    "20": "Microglia", 
    "21": "Excitatory neuron 10",
    "22": "Inhibitory neuron 6", 
    "23": "Inhibitory neuron 7", 
    "24": "Excitatory neuron 11",
    "25": "Excitatory neuron 12", 
    "26": "Astrocyte 3", 
    "27": "Excitatory neuron 13",
    "28": "Artifact_1_HighRibo",  
    "29": "Inhibitory neuron 8", 
    "30": "Artifact_2_HighDoublets",  
    "31": "Excitatory neuron 14", 
    "32": "Unknown",
    "33": "Excitatory neuron 15", 
    "34": "Excitatory neuron 16",
    "35": "Artifact_3_MixedType",  
    "36": "Artifact_4_MixedType"
}
adata.obs['cluster'] = adata.obs.leiden_scVI.map(cluster)


## Exclude artifacts for export
mask = ['Artifact', 'Unknown']
adata = adata[~adata.obs['cell_class'].isin(mask)].copy()

###############################################################################

# Export Annotated Data

adata.X = csr_matrix(adata.X)
adata.write_h5ad(os.path.join(save_dir, '2-tbi-annotated-hvg.h5ad'))

###############################################################################
