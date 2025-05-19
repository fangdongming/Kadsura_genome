
work_dir = '/hwfssz1/ST_EARTH/P20Z10200N0035/USER/renyongjuan/KCOC/single_cell/05.hotspot/mingene40_2/'

import hotspot
import warnings
warnings.filterwarnings('ignore')
import scanpy as sc
import pandas as pd
from scipy import sparse
import anndata as ad
import numpy as np
import seaborn as sns 
import matplotlib.pyplot as plt
import mplscience
import pickle
import os
import time
from scipy.sparse import csr_matrix
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=300)
# sc.logging.print_versions()
from scanpy.plotting import palettes
my_palettes = ['#1f77b4', '#ff7fee', '#279e68', '#d62728', '#aa40fc', '#8c564b', '#e377c2', '#b5bd61', '#17becf', '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5bed5', '#c49c94']
default_20 = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5']

os.chdir(work_dir)
time.sleep(15)

data = pd.read_csv('transposed_matrix.csv')

#adata1 = ad.read_h5ad('/hwfssz1/ST_EARTH/P20Z10200N0035/USER/nianting/C4_scRNA/Rice/2023.02/harmony_try2.soupx.r1_gene500_8000_lambda1/data.cc.h5ad')
#adata1.X[adata1.X == None] = 0
#adata1.write('noNaN.h5ad')

adata = sc.AnnData(data.iloc[:, 1:].values)
adata.obs_names = data.iloc[:, 0]
adata.var_names = data.columns[1:]

#adata = ad.read_h5ad('noNaN.h5ad')
adata.X  = csr_matrix(adata.X)
adata.obs["total_counts"] = np.asarray(adata.X.sum(1)).ravel()
adata.layers["csc_counts"] = adata.X.tocsc()

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# #renormalize the data for expression viz on plots
# sc.pp.normalize_total(adata)
# sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=2000)
sc.pl.highly_variable_genes(adata)

interested_genes = []
gene_list = pd.read_csv("genelist.csv")["Gene"].tolist()
interested_genes.extend(gene_list)

genes_not_found = [gene for gene in interested_genes if gene not in adata.var_names]
df = pd.DataFrame({'Non_existent_Gene': genes_not_found})
df.to_csv('Non_existent_gene.csv', index=False)
    
interested_genes = [gene for gene in interested_genes if gene in adata.var_names]
adata.var.loc[interested_genes, 'highly_variable'] = True
adata = adata[:, adata.var.highly_variable]

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.layers["log_normalized"] = adata.X.copy()
sc.pp.scale(adata)
sc.tl.pca(adata)

with mplscience.style_context():
    sc.pl.pca_variance_ratio(adata)

hs = hotspot.Hotspot(
    adata,
    layer_key="csc_counts",
    model='danb',
    latent_obsm_key="X_pca",
    umi_counts_obs_key="total_counts",
)

hs.create_knn_graph(
    weighted_graph=False, n_neighbors=300,
)

hs_results = hs.compute_autocorrelations(jobs=8)
hs_results.head()
df1 = pd.DataFrame(hs_results)
df1.to_csv("hs_results.csv")

hs_genes = hs_results.index[hs_results.index.isin(gene_list) | (hs_results.FDR < 0.00001)]
lcz = hs.compute_local_correlations(hs_genes, jobs=8)

df2 = pd.DataFrame(hs.local_correlation_z)
df2.to_csv("hs.local_correlation_z.csv")


with open("hs", "wb") as file:
    pickle.dump(hs, file)

with open("lcz", "wb") as file:
    pickle.dump(lcz, file)

modules = hs.create_modules(min_gene_threshold=40, core_only=False, fdr_threshold=0.05)
modules.value_counts()
df3 = pd.DataFrame(modules)
df3.to_csv("modules.csv")

hs.plot_local_correlations()
# plt.gcf().set_size_inches(12, 10) 
plt.savefig("local_correlations.pdf")

module_scores = hs.calculate_module_scores()
ms = pd.DataFrame(module_scores)
ms.to_csv("module_scores.csv")

module_cols = []
for c in module_scores.columns:
    key = "Module {}".format(c)
    adata.obs[key] = module_scores[c]
    module_cols.append(key)
