#!/usr/bin/env python3

import os
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce

# Set scanpy and matplotlib configurations
sc.set_figure_params(scanpy=True)
plt.rcParams['figure.dpi'] = 600
plt.rcParams['font.family'] = 'Arial'

###############################################################################

# Import data

data_dir = "<DATA_DIRECTORY>"
save_dir = "<SAVE_DIRECTORY>"

adata = sc.read_h5ad(os.path.join(data_dir, '2-tbi-annotated-full.h5ad'))
bdata = sc.read_h5ad(os.path.join(data_dir, '3-neuron-recluster-full.h5ad'))

###############################################################################
###############################################################################

# FIGURE 4

## UMAPs
sc.pl.umap(adata, color=['cell_type'], legend_fontsize=14, title='Cell types', save='_cell_type.pdf')
sc.pl.umap(bdata, color=['neuron_cluster'], legend_fontsize=10, title='Neuronal clusters', save='_neuron_cluster.pdf')

###############################################################################

# Cell type differential expression heatmaps

final_results = pd.read_csv(os.path.join(save_dir, 'wald-test-cell_type.csv'))
filtered_results = pd.read_csv(os.path.join(save_dir, 'wald-test-cell_type-filtered.csv'))

## Create heatmaps of upregulated and downregulated genes by cell type
comparisons = ['AC', 'CE', 'AB']
order = ['Astrocyte', 'Excitatory neuron', 'Inhibitory neuron', 'Microglia', 'OPC', 'Oligodendrocyte', 'Fibroblast', 'Choroid-plexus epithelial']
cell_type_reference = pd.DataFrame(order, columns=['cell_type'])

## Prepare lists for plotting
data_list_up = []
data_list_down = []


## Heatmap visualization
mean_expression_threshold = 0.1
for comp in comparisons:
    df_comp = filtered_results[(filtered_results['comparison'] == comp) & (filtered_results['mean'] > mean_expression_threshold)]
    deg_counts_up = df_comp[(df_comp['qval'] < 0.05) & (df_comp['log2fc'] > 0.5)].groupby('cell_type').size().reset_index(name=comp)
    deg_counts_up = pd.merge(cell_type_reference, deg_counts_up, on='cell_type', how='left').fillna(0)
    data_list_up.append(deg_counts_up)

    deg_counts_down = df_comp[(df_comp['qval'] < 0.05) & (df_comp['log2fc'] < -0.5)].groupby('cell_type').size().reset_index(name=comp)
    deg_counts_down = pd.merge(cell_type_reference, deg_counts_down, on='cell_type', how='left').fillna(0)
    data_list_down.append(deg_counts_down)

combined_df_up = reduce(lambda left, right: pd.merge(left, right, on='cell_type', how='outer'), data_list_up).fillna(0)
combined_df_down = reduce(lambda left, right: pd.merge(left, right, on='cell_type', how='outer'), data_list_down).fillna(0)

vmin = min(combined_df_up.drop(columns='cell_type').min().min(), combined_df_down.drop(columns='cell_type').min().min())
vmax = max(combined_df_up.drop(columns='cell_type').max().max(), combined_df_down.drop(columns='cell_type').max().max())

plt.figure(figsize=(6, 3))
sns.heatmap(combined_df_up.set_index('cell_type'), annot=True, cmap='viridis', fmt="g", vmin=vmin, vmax=vmax)
plt.title('Upregulated Genes')
plt.xticks(ticks=range(len(comparisons)), labels=comparisons)
plt.tight_layout()
plt.savefig(os.path.join(save_dir, 'upregulated-genes-summary.pdf'))
plt.show()

plt.figure(figsize=(6, 3))
sns.heatmap(combined_df_down.set_index('cell_type'), annot=True, cmap='viridis', fmt="g", vmin=vmin, vmax=vmax)
plt.title('Downregulated Genes')
plt.xticks(ticks=range(len(comparisons)), labels=comparisons)
plt.tight_layout()
plt.savefig(os.path.join(save_dir, 'downregulated-genes-summary.pdf'))
plt.show()

###############################################################################

# Neuron cluster composition bar chart

groups = ['A', 'B', 'C', 'E']
cdata = bdata[bdata.obs['group'].isin(groups)].copy()

condition = {
    'A': r'Sham + AAV$^{\mathrm{GFP}}$',
    'B': r'Sham + AAV$^{\mathrm{VEGFC}}$',
    'C': r'TBI + AAV$^{\mathrm{GFP}}$',
    'E': r'TBI + AAV$^{\mathrm{VEGFC}}$'
}
cdata.obs['condition'] = cdata.obs['group'].map(condition)

cluster_counts = cdata.obs.groupby(['condition', 'neuron_cluster']).size().reset_index(name='Count')
total_neurons_per_group = cdata.obs.groupby('condition').size().reset_index(name='Total_Neurons_group')
cluster_counts = pd.merge(cluster_counts, total_neurons_per_group, on='condition')
cluster_counts['Frequency'] = cluster_counts['Count'] / cluster_counts['Total_Neurons_group']

plt.figure(figsize=(3, 5))
g = sns.catplot(
    x='neuron_cluster', 
    y='Frequency', 
    hue='condition', 
    data=cluster_counts, 
    kind='bar', 
    height=3, 
    aspect=2, 
    palette='viridis'
)
g.set_axis_labels(" ", "% Neuronal Population\nper Group")
g._legend.set_title('Group')
g._legend.set_bbox_to_anchor((0.95, .9))
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(os.path.join(save_dir, "cluster_frequencies.pdf"), format='pdf')
plt.show()

###############################################################################

# Pathology-associated genes heatmaps (Holtzmann, with Zbtb20)

plt.rcParams['font.size'] = 16

genes_pathology_with_zbtb20 = ['Arpp21', 'R3hdm1', 'Rorb', 'Cux1', 'Cux2', 'Brinp3', 'Mef2c', 'Zbtb20']
sc.pl.matrixplot(cdata, var_names=genes_pathology_with_zbtb20, groupby='neuron_cluster', standard_scale='var', dendrogram=True, swap_axes=True, save='pathology-genes-cluster.pdf')

sc.pl.matrixplot(cdata, var_names=genes_pathology_with_zbtb20, groupby='condition', standard_scale='var', dendrogram=True, swap_axes=True, save='pathology-genes-condition.pdf')

###############################################################################
###############################################################################

# SUPPLEMENTAL FIGURE 4

###############################################################################

# UMAPs for all clusters

sc.pl.umap(adata, color='cluster', title='Cluster')

###############################################################################

# Cell count by cell type

bdata = adata[adata.obs['group'].isin(['A', 'B', 'C', 'E'])].copy()

condition = {
    'A': r'Sham + AAV$^{\mathrm{GFP}}$',
    'B': r'Sham + AAV$^{\mathrm{VEGFC}}$',
    'C': r'TBI + AAV$^{\mathrm{GFP}}$',
    'E': r'TBI + AAV$^{\mathrm{VEGFC}}$'
}
bdata.obs['condition'] = bdata.obs['group'].map(condition)

## Plot cell type counts by condition
df = pd.DataFrame(bdata.obs)
count_df = df.groupby(['condition', 'cell_type']).size().reset_index(name='Count')

g = sns.catplot(
    x='cell_type',
    y='Count',
    hue='condition',  
    data=count_df,
    kind='bar',
    height=5,
    aspect=0.75,
    palette='viridis',
    legend=True
)
g._legend.set_title('Group')
g._legend.set_bbox_to_anchor((.97, 0.85))
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(os.path.join(save_dir, "cell_type-counts-by-group.pdf"), format='pdf')
plt.show()

###############################################################################

# Dotplot of top cluster markers

sc.tl.rank_genes_groups(adata, groupby='cluster', method='wilcoxon')
cluster_markers = sc.get.rank_genes_groups_df(adata, group=None)
gene_list = cluster_markers.groupby('group').head(3)
genes = gene_list['names'].tolist()

sc.pl.dotplot(adata, var_names=genes, groupby='cluster')

###############################################################################

# Neuronal health-associated genes heatmaps (Holtzmann)

genes_neuronal_health = ['Prox1', 'Synpr', 'C1ql2', 'C1ql3', 'Camk2a', 'Camk2b', 'Tmem108', 'Ppfia2', 'Rfx3', 'Lrrtm4', 'Btbd9', 'Cntnap5a', 'Erc2']

sc.pl.matrixplot(cdata, var_names=genes_neuronal_health, groupby='neuron_cluster', standard_scale='var', dendrogram=True, swap_axes=True, save='neuronal-health-genes-cluster.pdf')

sc.pl.matrixplot(cdata, var_names=genes_neuronal_health, groupby='condition', standard_scale='var', dendrogram=True, swap_axes=True, save='neuronal-health-genes-condition.pdf')

###############################################################################
