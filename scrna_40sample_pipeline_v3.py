#!/usr/bin/env python3
"""
scRNA-seq Pipeline: 40 Samples (20 Parental vs 20 Erlotinib Resistant)
Workflow: Load → QC → Doublet Removal → Batch Correction → Annotation → DEG → GSEA
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scrublet as scr
import celltypist
from celltypist import models
import warnings
import os
warnings.filterwarnings('ignore')

sc.settings.verbosity = 2
sc.settings.n_jobs = 8
sc.settings.set_figure_params(dpi=150, frameon=False)

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_DIR = "/data/ankita/scRNA"

with open(f"{BASE_DIR}/parental_f", 'r') as f:
    PARENTAL_SAMPLES = [line.strip() for line in f.readlines() if line.strip()]
with open(f"{BASE_DIR}/resistant_f", 'r') as f:
    RESISTANT_SAMPLES = [line.strip() for line in f.readlines() if line.strip()]

print(f"Parental samples ({len(PARENTAL_SAMPLES)}): {PARENTAL_SAMPLES}")
print(f"Resistant samples ({len(RESISTANT_SAMPLES)}): {RESISTANT_SAMPLES}")

sample_info = []
for sample in PARENTAL_SAMPLES:
    sample_info.append({
        'sample_id': sample,
        'condition': 'parental',
        'h5_path': f"{BASE_DIR}/normal/{sample}_out_filtered_feature_bc_matrix.h5"
    })
for sample in RESISTANT_SAMPLES:
    sample_info.append({
        'sample_id': sample,
        'condition': 'resistant',
        'h5_path': f"{BASE_DIR}/resistant/{sample}_out_filtered_feature_bc_matrix.h5"
    })

sample_df = pd.DataFrame(sample_info)
print(f"\nTotal samples: {len(sample_df)}")
print(sample_df.groupby('condition').size())

# =============================================================================
# STEP 1: LOAD DATA
# =============================================================================
print("\n" + "="*70)
print("LOADING SAMPLES")
print("="*70)

adatas = []
for idx, row in sample_df.iterrows():
    try:
        adata = sc.read_10x_h5(row['h5_path'])
        adata.var_names_make_unique()
        adata.obs['sample_id'] = row['sample_id']
        adata.obs['condition'] = row['condition']
        adatas.append(adata)
        print(f"{row['sample_id']}: {adata.n_obs} cells")
    except FileNotFoundError:
        print(f"ERROR: File not found for {row['sample_id']}: {row['h5_path']}")

adata = sc.concat(adatas, join='outer')
adata.obs_names_make_unique()
print(f"\nTotal: {adata.n_obs} cells, {adata.n_vars} genes")

# =============================================================================
# STEP 2: QC
# =============================================================================
print("\n" + "="*70)
print("QUALITY CONTROL")
print("="*70)

adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             groupby='condition', multi_panel=True, save='_qc.png')

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, max_genes=8000)
adata = adata[adata.obs.pct_counts_mt < 20, :]
sc.pp.filter_genes(adata, min_cells=10)

print(f"After QC: {adata.n_obs} cells, {adata.n_vars} genes")

# =============================================================================
# STEP 3: DOUBLET REMOVAL
# =============================================================================
print("\n" + "="*70)
print("DOUBLET REMOVAL")
print("="*70)

doublet_scores = []
doublet_predictions = []

for sample in adata.obs['sample_id'].unique():
    sample_data = adata[adata.obs['sample_id'] == sample].copy()
    scrub = scr.Scrublet(sample_data.X, expected_doublet_rate=0.06)
    scores, predicted = scrub.scrub_doublets(min_counts=2, min_cells=3,
                                              min_gene_variability_pctl=85, n_prin_comps=30)
    doublet_scores.extend(scores)
    doublet_predictions.extend(predicted)

adata.obs['doublet_score'] = doublet_scores
adata.obs['predicted_doublet'] = doublet_predictions

print(f"Doublets: {adata.obs['predicted_doublet'].sum()} ({(adata.obs['predicted_doublet'].sum()/adata.n_obs)*100:.1f}%)")

adata = adata[~adata.obs['predicted_doublet']].copy()
print(f"After removal: {adata.n_obs} cells")

# =============================================================================
# STEP 4: NORMALIZATION & HVG
# =============================================================================
print("\n" + "="*70)
print("NORMALIZATION")
print("="*70)

adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['lognorm'] = adata.X.copy()

sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor='seurat_v3',
                            batch_key='sample_id', subset=False)
print(f"HVG: {adata.var.highly_variable.sum()}")

# =============================================================================
# STEP 5: BATCH CORRECTION
# =============================================================================
print("\n" + "="*70)
print("BATCH CORRECTION")
print("="*70)

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'], n_jobs=8)
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['condition', 'sample_id'], save='_before_batch.png')

try:
    import harmonypy as hm
    ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, vars_use=['sample_id'], max_iter_harmony=20)
    adata.obsm['X_pca_harmony'] = ho.Z_corr.T
    print("Harmony batch correction applied")
except ImportError:
    print("Harmony not installed, using ComBat")
    sc.pp.combat(adata, key='sample_id')
    adata.obsm['X_pca_harmony'] = adata.obsm['X_pca'].copy()

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, use_rep='X_pca_harmony')
sc.tl.umap(adata)
sc.pl.umap(adata, color=['condition', 'sample_id'], save='_after_batch.png')

# =============================================================================
# STEP 6: CLUSTERING
# =============================================================================
print("\n" + "="*70)
print("CLUSTERING")
print("="*70)

sc.tl.leiden(adata, resolution=0.8)
print(f"Clusters: {adata.obs['leiden'].nunique()}")

sc.pl.umap(adata, color=['leiden', 'condition'], save='_clusters.png')

# =============================================================================
# STEP 7: CELL TYPE ANNOTATION
# =============================================================================
print("\n" + "="*70)
print("CELL TYPE ANNOTATION")
print("="*70)

adata_pred = adata.copy()
adata_pred.X = adata_pred.layers['lognorm'].copy()

print("Running CellTypist annotation with Immune_All_Low model...")
try:
    predictions = celltypist.annotate(adata_pred, model='Immune_All_Low.pkl', majority_voting=True)
    
    adata.obs['cell_type'] = predictions.predicted_labels.majority_voting.values
    
    if 'conf_score' in predictions.predicted_labels.columns:
        adata.obs['cell_type_conf'] = predictions.predicted_labels.conf_score.values
    
    print("CellTypist annotation complete!")
    
except Exception as e:
    print(f"Error: {e}")
    raise

print(f"\nCell types:\n{adata.obs['cell_type'].value_counts()}")

sc.pl.umap(adata, color=['cell_type', 'condition'], save='_celltypes.png')

celltype_by_condition = pd.crosstab(adata.obs['cell_type'], adata.obs['condition'])
celltype_by_condition.to_csv('celltype_by_condition.csv')

# =============================================================================
# STEP 8: SAMPLE-WISE GENE EXPRESSION DATASET
# =============================================================================
print("\n" + "="*70)
print("GENERATING SAMPLE-WISE EXPRESSION DATASET")
print("="*70)

genes_list = adata.var_names.tolist()
samples_list = sorted(adata.obs['sample_id'].unique())

print(f"Total genes: {len(genes_list)}")
print(f"Total samples: {len(samples_list)}")
print("Computing mean expression per sample...")

expr_data = {}

for i, sample in enumerate(samples_list):
    print(f"  Processing sample {i+1}/{len(samples_list)}: {sample}")
    sample_subset = adata[adata.obs['sample_id'] == sample]
    
    if hasattr(sample_subset.X, 'toarray'):
        expr_array = sample_subset.X.toarray()
    else:
        expr_array = np.array(sample_subset.X)
    
    sample_means = expr_array.mean(axis=0)
    expr_data[sample] = sample_means

expr_matrix = pd.DataFrame(expr_data, index=genes_list)

expr_matrix.to_csv('sample_wise_gene_expression.csv')
print(f"Saved: sample_wise_gene_expression.csv")
print(f"Matrix shape: {expr_matrix.shape} (genes x samples)")

# =============================================================================
# STEP 9: GENERATE GSEA INPUT FILES
# =============================================================================
print("\n" + "="*70)
print("STEP 9: GENERATING GSEA INPUT FILES")
print("="*70)

sample_conditions = {}
for sample in samples_list:
    sample_conditions[sample] = adata[adata.obs['sample_id'] == sample].obs['condition'].iloc[0]

resistant_samples = [s for s in samples_list if sample_conditions[s] == 'resistant']
parental_samples = [s for s in samples_list if sample_conditions[s] == 'parental']

print(f"Parental samples: {len(parental_samples)}")
print(f"Resistant samples: {len(resistant_samples)}")

# FILE 1: Combined_All_Samples.gct (ALL SAMPLES)
print("\n>>> Creating Combined_All_Samples.gct...")
with open('Combined_All_Samples.gct', 'w') as f:
    f.write("#1.2\n")
    f.write(f"{len(expr_matrix)}\t{len(samples_list)}\n")
    f.write("NAME\tDescription\t" + "\t".join(samples_list) + "\n")
    for gene in expr_matrix.index:
        f.write(gene + "\t\t" + "\t".join([str(x) for x in expr_matrix.loc[gene].values]) + "\n")
print(">>> SAVED: Combined_All_Samples.gct")

# FILE 2: Parental.gct
print("\n>>> Creating Parental.gct...")
parental_expr = expr_matrix[parental_samples]
with open('Parental.gct', 'w') as f:
    f.write("#1.2\n")
    f.write(f"{len(parental_expr)}\t{len(parental_expr.columns)}\n")
    f.write("NAME\tDescription\t" + "\t".join(parental_expr.columns) + "\n")
    for gene in parental_expr.index:
        f.write(gene + "\t\t" + "\t".join([str(x) for x in parental_expr.loc[gene].values]) + "\n")
print(">>> SAVED: Parental.gct")

# FILE 3: Resistant.gct
print(">>> Creating Resistant.gct...")
resistant_expr = expr_matrix[resistant_samples]
with open('Resistant.gct', 'w') as f:
    f.write("#1.2\n")
    f.write(f"{len(resistant_expr)}\t{len(resistant_expr.columns)}\n")
    f.write("NAME\tDescription\t" + "\t".join(resistant_expr.columns) + "\n")
    for gene in resistant_expr.index:
        f.write(gene + "\t\t" + "\t".join([str(x) for x in resistant_expr.loc[gene].values]) + "\n")
print(">>> SAVED: Resistant.gct")

# FILE 4: Combined.cls
print(">>> Creating Combined.cls...")
with open('Combined.cls', 'w') as f:
    f.write(f"{len(samples_list)} 2 1\n")
    f.write("# parental resistant\n")
    class_assignments = ["0" if sample_conditions[s] == 'parental' else "1" for s in samples_list]
    f.write(" ".join(class_assignments) + "\n")
print(">>> SAVED: Combined.cls")

print("\n" + "="*70)
print("GSEA FILES SAVED")
print("="*70)
print("Files created in current directory:")
print("  - Combined_All_Samples.gct (ALL samples)")
print("  - Parental.gct")
print("  - Resistant.gct")
print("  - Combined.cls")
# =============================================================================
# STEP 10: DEG ANALYSIS
# =============================================================================
print("\n" + "="*70)
print("DIFFERENTIAL EXPRESSION")
print("="*70)

sc.tl.rank_genes_groups(adata, groupby='condition', method='wilcoxon',
                        key_added='deg_condition', use_raw=False, layer='lognorm', pts=True)

deg_all = sc.get.rank_genes_groups_df(adata, group='resistant', key='deg_condition')
deg_all_sig = deg_all[(deg_all['pvals_adj'] < 0.05) & (abs(deg_all['logfoldchanges']) > 0.5)]
deg_all_sig.to_csv('deg_resistant_vs_parental_all.csv', index=False)
print(f"Total DEGs: {len(deg_all_sig)}")

# DEG by cell type
all_ct_degs = []
for ct in adata.obs['cell_type'].unique():
    adata_sub = adata[adata.obs['cell_type'] == ct].copy()

    n_p = (adata_sub.obs['condition'] == 'parental').sum()
    n_r = (adata_sub.obs['condition'] == 'resistant').sum()

    if n_p < 20 or n_r < 20:
        print(f"Skipping {ct}: insufficient cells (P:{n_p}, R:{n_r})")
        continue

    sc.tl.rank_genes_groups(adata_sub, groupby='condition', method='wilcoxon',
                           use_raw=False, layer='lognorm', pts=True)

    deg = sc.get.rank_genes_groups_df(adata_sub, group='resistant')
    deg['cell_type'] = ct
    deg_sig = deg[(deg['pvals_adj'] < 0.05) & (abs(deg['logfoldchanges']) > 0.5)]
    all_ct_degs.append(deg_sig)
    print(f"{ct}: {len(deg_sig)} DEGs")

if all_ct_degs:
    pd.concat(all_ct_degs).to_csv('deg_resistant_vs_parental_by_celltype.csv', index=False)

if len(deg_all_sig) > 0:
    top_genes = deg_all_sig.nlargest(15, 'logfoldchanges')['names'].tolist()
    top_genes += deg_all_sig.nsmallest(15, 'logfoldchanges')['names'].tolist()
    top_genes = [g for g in top_genes if g in adata.var_names][:30]

    sc.pl.dotplot(adata, var_names=top_genes, groupby='condition',
                  standard_scale='var', save='_degs.png')

# =============================================================================
# STEP 11: SAVE RESULTS
# =============================================================================
print("\n" + "="*70)
print("SAVING RESULTS")
print("="*70)

adata.write('processed_data.h5ad')

cell_metadata = pd.DataFrame({
    'cell_id': adata.obs_names,
    'sample_id': adata.obs['sample_id'],
    'condition': adata.obs['condition'],
    'cluster': adata.obs['leiden'],
    'cell_type': adata.obs['cell_type'],
    'n_genes': adata.obs['n_genes_by_counts'],
    'pct_mt': adata.obs['pct_counts_mt'],
    'UMAP1': adata.obsm['X_umap'][:, 0],
    'UMAP2': adata.obsm['X_umap'][:, 1]
})
cell_metadata.to_csv('cell_annotations.csv', index=False)

print("\nANALYSIS COMPLETE!")
print(f"Cells: {adata.n_obs} | Clusters: {adata.obs['leiden'].nunique()} | Cell types: {adata.obs['cell_type'].nunique()}")
n_degs_by_ct = len(pd.concat(all_ct_degs)) if all_ct_degs else 0
print(f"DEGs (all): {len(deg_all_sig)} | DEGs (by cell type): {n_degs_by_ct}")
print("\nKey outputs:")
print("  - processed_data.h5ad")
print("  - cell_annotations.csv")
print("  - sample_wise_gene_expression.csv")
print("  - Parental.gct")
print("  - Resistant.gct")
print("  - Combined.cls")
print("  - deg_resistant_vs_parental_all.csv")
print("  - deg_resistant_vs_parental_by_celltype.csv")
print("  - celltype_by_condition.csv")
