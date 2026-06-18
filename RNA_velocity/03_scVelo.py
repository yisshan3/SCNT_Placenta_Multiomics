"""
==============================================================================
Description: Combine velocyto loom files and run scVelo RNA velocity analysis.
==============================================================================
"""


# Setup Environment ------------------------------------------------------------
import warnings
warnings.filterwarnings("ignore")

import os
import loompy
import pandas as pd
import scanpy as sc
import scvelo as scv
import sys

group = "SCNT"  # WT or SCNT

base_dir = "/home4/ssyi/Mouse_Placenta"
meta_dir = os.path.join(base_dir, "codes/4.RNA_velocity")
loom_dir = os.path.join(base_dir, "codes/4.RNA_velocity/loom")
out_dir = os.path.join(base_dir, "codes/4.RNA_velocity/scVelo")
fig_dir = os.path.join(base_dir, "Figures")

os.makedirs(loom_dir, exist_ok=True)
os.makedirs(out_dir, exist_ok=True)
os.makedirs(fig_dir, exist_ok=True)

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo')


# Combine Loom Files -----------------------------------------------------------
combined_loom_path = os.path.join(loom_dir, "combined_loom", f"{group}_Tropho_combined.loom")
os.makedirs(os.path.dirname(combined_loom_path), exist_ok=True)

if group == "WT":
    loom_files = [
        os.path.join(loom_dir, "WT_batch1_125.loom"),
        os.path.join(loom_dir, "WT_batch2_145.loom"),
        os.path.join(loom_dir, "WT_batch4_185.loom")
    ]
elif group == "SCNT":
    loom_files = [
        os.path.join(loom_dir, "SCNT_batch3_125.loom"),
        os.path.join(loom_dir, "SCNT_batch5_145.loom"),
        os.path.join(loom_dir, "SCNT_batch6_185.loom")
    ]
else:
    raise ValueError("Invalid group specified. Must be 'WT' or 'SCNT'.")

if not os.path.exists(combined_loom_path):
    loompy.combine(loom_files, combined_loom_path, key="Accession")


# Load and Format Loom Data ----------------------------------------------------
loom_data = sc.read_loom(combined_loom_path)
loom_data.obs.index = loom_data.obs.index.str.replace(':', '_', regex=False).str.replace('x', '-1', regex=False)

for i in range(1, 7):
    loom_data.obs.index = loom_data.obs.index.str.replace(f'_batch{i}_', '', regex=False)


# Load and Align with Metadata -------------------------------------------------
obs_file = f"{group}_cellID_obs.csv"
umap_file = f"{group}_cell_embeddings.csv"
celltype_file = f"{group}_cell_celltype.csv"

sample_obs = pd.read_csv(os.path.join(meta_dir, obs_file))
cell_umap = pd.read_csv(os.path.join(meta_dir, umap_file), header=0, names=["CellID", "UMAP_1", "UMAP_2"])
cell_celltype = pd.read_csv(os.path.join(meta_dir, celltype_file), header=0, names=["CellID", "celltype"])

# Subset loom data to keep only the cells of interest
valid_cells = sample_obs.iloc[:, 0].values
adata = loom_data[loom_data.obs_names.isin(valid_cells)].copy()

cell_umap.set_index("CellID", inplace=True)
cell_celltype.set_index("CellID", inplace=True)

adata.obsm['X_umap'] = cell_umap.loc[adata.obs_names].values
adata.obs['celltype'] = cell_celltype.loc[adata.obs_names, 'celltype'].values
adata.var_names_make_unique()

# Data validation
print(f"Data Summary for {group}:")
print(f"  Number of cells: {adata.n_obs}")
print(f"  Number of genes: {adata.n_vars}")
print(f"  Cell types: {list(adata.obs['celltype'].unique())}")
print(f"  Spliced counts range: {adata.layers['spliced'].sum(axis=1).min():.0f} - {adata.layers['spliced'].sum(axis=1).max():.0f}")
print(f"  Unspliced counts range: {adata.layers['unspliced'].sum(axis=1).min():.0f} - {adata.layers['unspliced'].sum(axis=1).max():.0f}")

# Save h5ad
h5ad_out = os.path.join(out_dir, f"{group}_Tropho_dynamicModel.h5ad")
adata.write(h5ad_out, compression='gzip')


# RNA Velocity Computation -----------------------------------------------------
scv.pl.proportions(adata)

# Filter low-count genes and compute moments
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# Compute velocity (Stochastic mode)
scv.tl.velocity(adata, mode="stochastic")
scv.tl.velocity_graph(adata)


# Visualization ----------------------------------------------------------------
color_dict = {
    "SynTI Precursor": "#1f78b4", "SynTI": "#a6cee3", 
    "S-TGC Precursor": "#33a02c", "S-TGC": "#b2df8a", 
    "SpT Precursor": "#e31a1c", "SpT": "#fb9a99",
    "SynTII Precursor": "#ff7f00", "SynTII": "#fdbf6f", 
    "JZP1": "#6a3d9a", "JZP2": "#cab2d6", 
    "LaTP": "#ffd92f", "LaTP2": "#ffff99", "GlyT": "#8da0cb"
}

n_cells = adata.n_obs

arrow_filename = f"{group}_scvelo_Trophoblast_UMAP_arrow.svg"
stream_filename = f"{group}_scvelo_Trophoblast_UMAP_stream.svg"

# Velocity arrow plot
scv.pl.velocity_embedding(
    adata, basis='X_umap', color="celltype", arrow_length=3, arrow_size=4, 
    alpha=0.5, dpi=120, palette=color_dict, 
    save=os.path.join(fig_dir, arrow_filename),  
    figsize=(8,7), legend_fontsize=9, show=False, 
    title=f'{group} (N={n_cells})', dpi=120
)

# Velocity stream plot
scv.pl.velocity_embedding_stream(
    adata, basis='X_umap', color="celltype", palette=color_dict, alpha=0.5,
    save=os.path.join(fig_dir, stream_filename),
    figsize=(8,7), legend_fontsize=0, show=False, 
    title=f'{group} (N={n_cells})', dpi=120
)

