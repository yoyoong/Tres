import argparse
import os
import numpy as np
import pandas as pd
import h5py
import scanpy as sc
from scipy.sparse import coo_matrix
from tqdm.autonotebook import tqdm


parser = argparse.ArgumentParser()
parser.add_argument('-D', "--dataset", type=str, required=False, help="Dataset name", default='Yost2019')
args = parser.parse_args()
                    
dataset = args.dataset
if dataset == 'Zhang2021':
    input_path = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/3.clinical_data/Zhang2021/raw_data'
elif dataset == 'SadeFeldman2018':
    input_path = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/3.clinical_data/SadeFeldman2018/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz'
elif dataset == 'Yost2019':
    input_path = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/3.clinical_data/Yost2019/GSE123813_bcc_scRNA_counts.txt.gz'
elif dataset == 'Fraietta2018':
    input_path = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/3.clinical_data/Fraietta2018/raw_bulk_profile.csv'
sample_annotation_path = f'/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/3.clinical_data/{dataset}/{dataset}.sample_annotation.txt'
output_file_directory = f'/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/3.clinical_data/{dataset}'
output_tag = f'{dataset}.bulk_profile'

# get sample annotation
sample_annotation = pd.read_csv(sample_annotation_path, sep='\t', index_col=0)
sample_list = sample_annotation.index.tolist()

dataset = input_path.split('/')[-2]
gem_bulk = pd.DataFrame()
if dataset == 'Zhang2021':
    gem_adata = sc.read_10x_mtx(input_path, gex_only=False)
    cellname_list = list(gem_adata.obs_names)
    genename_list = list(gem_adata.var_names)
    gem_df = pd.DataFrame(gem_adata.X.T.toarray(), index=genename_list, columns=cellname_list)

    # filter the gem
    columms_filtered = [col for col in gem_df.columns if col.split('.')[1] in sample_list]
    gem_df = gem_df[columms_filtered]
    flag_group = [col.split('.')[1] for col in gem_df.columns]

elif dataset == 'SadeFeldman2018':
    gem_df = pd.read_csv(input_path, compression="gzip", header=[0, 1], sep='\t')
    gem_df.columns = gem_df.columns.map('.'.join)
    columms_filtered = [col for col in gem_df.columns if col.split('.')[1] in sample_list]
    gem_df = gem_df[columms_filtered]  # filter the gem
    flag_group = [col.split('.')[1] for col in gem_df.columns]

elif dataset == 'Yost2019':
    gem_df = pd.read_csv(input_path, compression="gzip", header=0, sep='\t')
    flag_group = [col.split('_')[0] for col in gem_df.columns]

elif dataset == 'Fraietta2018':
    gem_df = pd.read_csv(input_path, header=0, sep='\t', index_col=0)
    flag_group = gem_df.columns

# group the gem and calculate average
gem_bulk = gem_df.groupby(flag_group, axis=1).mean()
gem_bulk.to_csv(os.path.join(output_file_directory, f'{output_tag}.raw.csv'))

gem_log = (gem_df + 1).apply(np.log2)
gem_bulk = gem_log.groupby(flag_group, axis=1).mean()
gem_bulk.to_csv(os.path.join(output_file_directory, f'{output_tag}.log.csv'))

# # normalize
gem_df *= 1E5 / gem_df.sum()
gem_df = (gem_df + 1).apply(np.log2)

# group the gem and calculate average
gem_bulk = gem_df.groupby(flag_group, axis=1).mean()
gem_bulk.to_csv(os.path.join(output_file_directory, f'{output_tag}.normalize_log.csv'))

# # centralize
gem_df = gem_df.sub(gem_df.mean(axis=1), axis=0)

# group the gem and calculate average
gem_bulk = gem_df.groupby(flag_group, axis=1).mean()
gem_bulk.to_csv(os.path.join(output_file_directory, f'{output_tag}.normalize_log_centralize.csv'))
print(f"{output_tag} process end")
