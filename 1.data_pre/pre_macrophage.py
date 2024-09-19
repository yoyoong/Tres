import argparse
import os
import numpy as np
import pandas as pd
import h5py
import scanpy as sc
from scipy.sparse import coo_matrix
from tqdm.autonotebook import tqdm
from sklearn.cluster import KMeans

parser = argparse.ArgumentParser()
parser.add_argument('-I', "--input", type=str, required=False, help="input h5 files path",
                    default='/sibcb1/bioinformatics/hongyuyang/dataset/Tres/4.macrophage/GSE168710')
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb1/bioinformatics/hongyuyang/dataset/Tres/4.macrophage/GSE168710')
parser.add_argument('-O', "--output_tag", type=str, required=False, help="Prefix for output files.",
                    default='GSE168710')
args = parser.parse_args()

input_path = '/sibcb1/bioinformatics/hongyuyang/dataset/Tres/4.macrophage/GSE168710'
output_file_directory = '/sibcb1/bioinformatics/hongyuyang/dataset/Tres/4.macrophage/GSE168710'
output_tag = 'GSE168710'

# file_list = os.listdir(input_path)
# h5file_list = [file for file in file_list if file.endswith('.h5')]
# adata_list = []
# for file in h5file_list:
#     adata = sc.read_10x_h5(os.path.join(input_path, file))
#     adata.var_names_make_unique()
#     adata_list.append(adata)
# adata = sc.concat(adata_list, join='inner', index_unique='.')

adata = sc.read_h5ad('/sibcb1/bioinformatics/hongyuyang/dataset/Tres/4.macrophage/GSE168710/GSE168710_QC_normalized_exp.h5ad')

# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata, base=2)
gem_array = adata.T.X
row_medians = np.median(gem_array, axis=1)
gem_centralized = gem_array - row_medians[:, np.newaxis]

gem_df = pd.DataFrame(gem_centralized)
gem_df.index = adata.var.index
gem_df.columns = adata.obs.index

gem_df.to_csv(os.path.join(output_file_directory, f'{output_tag}_normalized.csv'))
print(f"{output_tag} process end")
