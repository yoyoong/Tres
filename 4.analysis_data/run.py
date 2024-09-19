import numpy as np
import pandas as pd
import scanpy as sc

# 读取10x数据
rawdata_path = "/sibcb1/bioinformatics/hongyuyang/dataset/Tres/4.analysis_data/sc/GSE236581/rawdata"
adata = sc.read_10x_mtx(path=rawdata_path, cache_compression='gzip')

# filter blood
tissue = ["N", "T", "T1", "T2", "LN", "TN"]
cell_names = adata.obs_names
tissue_parts = [name.split('-')[1] for name in cell_names]
mask = [part in tissue for part in tissue_parts]
new_adata = adata[mask, :]

sc.pp.normalize_total(new_adata, target_sum=1e4)
sc.pp.log1p(new_adata)
# sc.pp.scale(adata)
sc.pp.highly_variable_genes(new_adata, n_top_genes=2000)
sc.tl.pca(new_adata, svd_solver='arpack')
sc.pp.neighbors(new_adata, n_neighbors=10, n_pcs=40)
sc.tl.leiden(new_adata)
print("labels_pred len = ", len(list(new_adata.obs.leiden)))

new_adata.write('/sibcb1/bioinformatics/hongyuyang/dataset/Tres/4.analysis_data/sc/GSE236581/adata.h5ad')

