import argparse
import numpy
import os

import numpy as np
import pandas as pd
import sys
import warnings
from tqdm.autonotebook import tqdm
warnings.filterwarnings("ignore")
import pickle

parser = argparse.ArgumentParser()
parser.add_argument('-CT', "--celltype", type=str, default='CD8T', required=False, help="cell type")
parser.add_argument('-CIV', "--cytokine_info_version", type=int, default='2', required=False, help="cytokine info version")
parser.add_argument('-CSV', "--cytokine_signature_version", type=int, default='1', required=False, help="cytokine signature version")
parser.add_argument('-SFV', "--sample_filter_version", type=int, default='1', required=False, help="sample filter version")
args = parser.parse_args()

celltype = args.celltype
cytokine_info_version = args.cytokine_info_version
cytokine_signature_version = args.cytokine_signature_version
sample_filter_version = args.sample_filter_version
if cytokine_signature_version == 2 and sample_filter_version >= 1:
    sys.exit("Unvalid parameters.")
if celltype == "CD8T" :
    interaction_path = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-1.CD8T_Interaction/dataset_interaction'
    output_file_directory = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-1.CD8T_Interaction'
elif celltype == "Macrophage" :
    interaction_path = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-2.Macrophage_Interaction/dataset_interaction'
    output_file_directory = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-2.Macrophage_Interaction'
elif celltype == "Neutrophils" :
    interaction_path = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-3.Neutrophils_Interaction/dataset_interaction'
    output_file_directory = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-3.Neutrophils_Interaction'
elif celltype == "NK" :
    interaction_path = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-4.NK_Interaction/dataset_interaction'
    output_file_directory = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-4.NK_Interaction'
elif celltype == "NK_act" :
    interaction_path = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-4-0.NK_act_Interaction/dataset_interaction'
    output_file_directory = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-4-0.NK_act_Interaction'

real_celltype = celltype
if "_" in real_celltype:
    real_celltype = real_celltype.split("_")[0]

cytokine_info_file = f'/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.model_file/cytokine_info_{cytokine_info_version}.{real_celltype}.txt'
cytokine_info_df = pd.read_csv(cytokine_info_file, index_col=0)
cytokine_list = cytokine_info_df.index.values.tolist()
cytokine_signature_dir = os.path.join(output_file_directory, f'cytokine_signature')

all_positive_signature = []
all_negative_signature = []
for cytokine in tqdm(cytokine_list, desc='Processing cytokine list'):
    cytokine_signature_filename = os.path.join(cytokine_signature_dir, f'{cytokine}_{cytokine_signature_version}.signature.csv')
    cytokine_signature = pd.read_csv(cytokine_signature_filename, index_col=0, header=0)

    if cytokine_info_df.loc[cytokine]['flag'] == '+':
        all_positive_signature.append(cytokine_signature)
    else:
        all_negative_signature.append(cytokine_signature)

def get_tres_signature(all_signature):
    # 求正/负cytokine的每个cytokine都有的样本
    common_col = None
    for signature in all_signature:
        if common_col is None:
            common_col = signature.columns
        else:
            common_col = common_col.intersection(signature.columns)
    for i, result in enumerate(all_signature):
        all_signature[i] = result.loc[:, common_col]

    # 求样本在cytokine的中位数
    median_signature = pd.concat(all_signature, axis=1, join='inner') # 合并正/负cytokine的signature
    median_signature = median_signature.groupby(median_signature.columns, axis=1).median() # 求每个样本在所有正/负cytokine中的中值
    # output to file
    # median_signature_filename = os.path.join(output_file_directory, f'Median_signature_{filter_flag}.{pn_flag}.csv')
    # median_signature.to_csv(median_signature_filename)

    # 过滤样本数，然后求样本的中位数
    if sample_filter_version == 0: # 不过滤，只要这个基因有有效值样本即求中位数
        tres_signature = median_signature.dropna(how='all').median(axis=1)
    elif sample_filter_version == 1: # 过滤有效值样本数占比小于50%的基因
        tres_signature = median_signature.loc[median_signature.isnull().mean(axis=1) < 0.5].median(axis=1)
    elif sample_filter_version == 2: # 只保留所有样本都有有效值的基因
        tres_signature = median_signature.dropna(how='any').median(axis=1)

    tres_signature.name = 'Tres'
    return tres_signature

Tres_signature_dir = os.path.join(output_file_directory, f'Tres_signature')
if not os.path.exists(Tres_signature_dir):
    os.mkdir(Tres_signature_dir)
version = f'{cytokine_info_version}_{cytokine_signature_version}_{sample_filter_version}'
positive_tres_signature = get_tres_signature(all_positive_signature)
positive_tres_signature.to_csv(os.path.join(Tres_signature_dir, f'Tres_signature_{version}.positive.csv'))
negative_tres_signature = get_tres_signature(all_negative_signature)
negative_tres_signature.to_csv(os.path.join(Tres_signature_dir, f'Tres_signature_{version}.negative.csv'))
print("Process end!")