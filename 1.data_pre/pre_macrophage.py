import argparse
import os
import numpy as np
import pandas as pd
import h5py
import scanpy as sc
from scipy.sparse import coo_matrix
from tqdm.autonotebook import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('-I', "--input", type=str, required=False, help="input h5 files path",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/4.macrophage/1.raw_data')
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/4.macrophage/1.raw_data')
parser.add_argument('-O', "--output_tag", type=str, required=False, help="Prefix for output files.",
                    default='GSE168710')
args = parser.parse_args()

input_path = args.input
output_file_directory = args.output_file_directory
output_tag = args.output_tag

file_list = os.listdir(input_path)
h5file_list = [file for file in file_list if file.endswith('.h5')]
adata_list = [sc.read_10x_h5(os.path.join(input_path, file)) for file in h5file_list]
for adata in adata_list:
    adata.var_names_make_unique()
merge_adata = sc.concat(adata_list, join='inner', index_unique='.')

sc.pp.normalize_total(merge_adata, target_sum=1e4)
sc.pp.log1p(merge_adata, base=2)
gem_array = merge_adata.T.X.toarray()
row_medians = np.median(gem_array, axis=1)
gem_centralized = gem_array - row_medians[:, np.newaxis]

gem_df = pd.DataFrame(gem_centralized)
gem_df.index = merge_adata.var.index
gem_df.columns = merge_adata.obs.index

gem_df.to_csv(os.path.join(output_file_directory, f'{output_tag}.csv'))
print(f"{output_tag} process end")
