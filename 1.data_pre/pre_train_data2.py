import argparse
import os
import sys
import numpy as np
import pandas as pd
import h5py
import scanpy as sc
from scipy.sparse import coo_matrix
from tqdm.autonotebook import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('-E', "--expression_path", type=str, required=False, help="Gene expression file.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/1.gem_data/AEL_GSE142213.csv')
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb2/bioinformatics2/hongyuyang/code/Tres')
parser.add_argument('-O', "--output_tag", type=str, required=False, help="Prefix for output files.",
                    default='test')
args = parser.parse_args()

expression_path = args.expression_path
output_file_directory = args.output_file_directory
output_tag = args.output_tag

if os.path.isdir(expression_path):
    print(f"This dataset don't need process.")
    sys.exit(1)
gem_df = pd.read_csv(expression_path, index_col=0, header=0)

gem_dir = os.path.join(output_file_directory, output_tag)
if not os.path.exists(gem_dir):
    os.makedirs(gem_dir)
else:
    file_list = os.listdir(gem_dir)
    for file in file_list:
        file_path = os.path.join(gem_dir, file)
        os.remove(file_path)

flag_group = [v.split('.')[0] for v in gem_df.columns]
gem_group = gem_df.groupby(flag_group, axis=1)
for celltype, gem_sub in tqdm(gem_group, desc="Celltype data_process"):
    celltype = celltype.replace("/", "_")
    gem_sub_path = os.path.join(gem_dir, f'{output_tag}.{celltype}.csv')
    gem_sub.to_csv(gem_sub_path)

print(f"{output_tag} process end")
