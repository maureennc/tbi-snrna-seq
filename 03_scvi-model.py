#!/usr/bin/env python3

# Import libraries

import os
import random
import numpy as np
import scanpy as sc
import torch
import scvi

###############################################################################

# Settings

## Set random seeds for reproducibility
random.seed(0)
torch.manual_seed(0)
np.random.seed(0)
scvi.settings.seed = 0

## Define data and model paths
data_dir = "<DATA_DIRECTORY_PATH>"
scvi_dir = "<SAVE_DIRECTORY_PATH_FOR_MODEL>"
os.makedirs(scvi_dir, exist_ok=True)

###############################################################################

# Load Data

adata = sc.read_h5ad(os.path.join(data_dir, '1-tbi-seq-hvg.h5ad'))

###############################################################################

# Setup and Train Model

## Configure the AnnData for scVI with covariates
scvi.model.SCVI.setup_anndata(
    adata,
    layer='counts',
    categorical_covariate_keys=['group'],
    continuous_covariate_keys=['total_counts', 'pct_counts_mt', 'pct_counts_ribosomal', 'doublet_scores']
)

## Initialize and train the model
model = scvi.model.SCVI(adata)
trainer = scvi.train.Trainer(accelerator='cpu', devices=1)
model.train()

## Save trained model
model_dir = os.path.join(scvi_dir, 'model.model')
os.makedirs(model_dir, exist_ok=True)
model.save(model_dir)

###############################################################################
