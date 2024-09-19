import argparse
import os
import sys
import numpy as np
import pandas as pd
import h5py
import scanpy as sc
from scipy.sparse import coo_matrix
from tqdm.autonotebook import tqdm
from numpy import float16

parser = argparse.ArgumentParser()
parser.add_argument('-E', "--expression_path", type=str, required=False, help="Gene expression file.",
                    default='/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_data/1.neutrophil_data/Gao2024/gem.csv')
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_data/1.neutrophil_data/Gao2024')
parser.add_argument('-O', "--output_tag", type=str, required=False, help="Prefix for output files.",
                    default='Gao2024.Neutrophils')
args = parser.parse_args()

expression_path = args.expression_path
output_file_directory = args.output_file_directory
output_tag = args.output_tag

if os.path.isdir(expression_path):
    print(f"This dataset don't need process.")
    sys.exit(1)
gem_df = pd.read_csv(expression_path, index_col=0, header=0)
print(f"Get gem success.")

sample_list = []
sample_flag = []
for column in gem_df.columns:
    column_list = column.split('_')
    if 'OA902' in column:
        sample_flag.append(False)
        continue
    sample_flag.append(True)
    if len(column_list) == 2:
        if 'OA089' in column:
            sample_list.append("_".join([column_list[0], column.split('-')[-1]]))
        else:
            sample_list.append(column_list[0])
    elif len(column_list) == 3:
        sample_list.append("_".join(column_list[0:-1]))
    elif len(column_list) == 4:
        sample_list.append("_".join(column_list[0:-1]))
    elif len(column_list) == 5:
        if 'OA091_A012_ICC' in column or 'OA091_A048_ICC' in column:
            sample_list.append("_".join(column_list[0:-1]))
        else:
            sample_list.append("_".join(column_list[0:3]))
    elif len(column_list) == 6:
        sample_list.append("_".join(column_list[0:3]))
    elif len(column_list) == 8:
        sample_list.append("_".join(column_list[0:4]))
    elif len(column_list) == 9:
        sample_list.append("_".join(column_list[0:1] + column_list[-2:]))
    elif len(column_list) == 10:
        sample_list.append("_".join(column_list[0:1]))
gem_filtered = gem_df.loc[:, sample_flag]
gem_filtered.columns = ['.0'.join(['Neutrophils', sample_name, cell_name]) for (sample_name, cell_name) in
                        zip(sample_list, list(gem_filtered.columns)) ]

# normalize
col_sums = gem_filtered.sum() # 计算每列的总和
non_zero_cols = col_sums[col_sums != 0].index # 找到非零列
for col in non_zero_cols: # 对非零列进行归一化
    gem_filtered[col] = gem_filtered[col] / gem_filtered[col].sum() * 10000
# log
gem_filtered = np.log2(gem_filtered + 1)
# centralize
gem_filtered = gem_filtered.subtract(gem_filtered.mean(axis=1), axis=0)
gem_filtered = gem_filtered.astype(float16)

group_flag = ['.'.join(v.split('.')[1]) for v in gem_filtered.columns]
gem_filtered_groupeby_sample = gem_filtered.groupby(group_flag)
for sample, gem_sub in tqdm(gem_filtered_groupeby_sample, desc="Processing sample"):
    gem_sub.to_csv(os.path.join(output_file_directory, f'{output_tag}.{sample}.csv'))

print(f"{output_tag} process end")
