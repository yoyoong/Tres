import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from tqdm.autonotebook import tqdm
sys.path.append('/sibcb2/bioinformatics2/hongyuyang/code/Tres')
from Util import read_expression
from csn_torch import upperlower, getCSNMatrix
import torch


# Argument parser
parser = argparse.ArgumentParser(description='get_csndm')
parser.add_argument('-E', "--expression_file", type=str, required=False, help="Gene expression file.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.Tres_data/sc_cohorts/Breast.GSE156728.10x.pickle.gz')
parser.add_argument('-B', '--boxsize', type=float, required=False, help='boxsize', default=0.1)
parser.add_argument('-A', '--alpha', type=float, required=False, help='alpha', default=0.01)
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/Breast/GSE156728')
parser.add_argument('-O', "--output_tag", type=str, required=False, help="Prefix for output files.", default='Breast.GSE156728')
args = parser.parse_args()

expression_file = args.expression_file
boxsize = args.boxsize
alpha = args.alpha
output_tag = args.output_tag if args.output_tag is not None else os.path.splitext(os.path.basename(expression_file))[0]
ndm_filename = os.path.join(args.output_file_directory, f'{output_tag}.NDM.csv')

# get gem data
gem_df = read_expression(expression_file)
gem = np.array(gem_df, dtype=np.float32)
gene_num = gem.shape[0]
gene_list = gem_df.index.tolist()
cell_num = gem.shape[1]
cell_list = gem_df.columns.tolist()
print(f'The dataset has {gene_num} genes, {cell_num} cells')

# Build CSN data
boxsize_str = str("{:.0e}".format(boxsize))
alpha_str = str("{:.0e}".format(alpha))
csn_info = f'boxsize{boxsize_str}_alpha{alpha_str}'

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f'Using device: {device}')
if device == 'cuda':
    print("当前系统上的GPU设备数量为：", torch.cuda.device_count())

(upper, lower) = upperlower(gem, boxsize=boxsize)
gem_tensor = torch.as_tensor(gem, dtype=torch.float32).to(device)

upper_tensor = torch.as_tensor(upper, dtype=torch.float32).to(device)
lower_tensor = torch.as_tensor(lower, dtype=torch.float32).to(device)
ndm = np.zeros((gene_num, cell_num), dtype=int)
total_edge_num = 0
for cell in tqdm(range(cell_num), desc="Get CSNDM"):
    adj_sparse = getCSNMatrix(gem_tensor, upper_tensor, lower_tensor, cell, is_weight=False, has_zero=False, alpha=alpha, device=device)
    adj_array = adj_sparse.toarray()
    ndm[:, cell] = np.sum(adj_array, axis=0)
    total_edge_num += np.sum(ndm[:, cell])

ndm_df = pd.DataFrame(ndm, index=gene_list, columns=cell_list)
ndm_df.to_csv(ndm_filename, sep='\t')
print(f"Sample {output_tag} get csn end!")

