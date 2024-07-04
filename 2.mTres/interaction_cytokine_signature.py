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
parser.add_argument('-CSV', "--cytokine_signature_version", type=int, default='2', required=False, help="cytokine signature version")
args = parser.parse_args()

celltype = args.celltype
cytokine_signature_version = args.cytokine_signature_version
model_matrix_file = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.model_file/signature.centroid.expand'
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

model_matrix_file = pd.read_csv(model_matrix_file, sep='\t', index_col=0)
all_cytokine_list = model_matrix_file.columns.values.tolist()

interaction_list = []
if os.path.isdir(interaction_path):
    interaction_list = sorted(os.listdir(interaction_path))
else:
    interaction_list.append(interaction_path)

def get_dataset_cytokine_signature(dataset_name, interaction_data, cytokine):
    # filter the cytokine
    cytokine_flag = [(v.split('.')[-1] == cytokine) for v in interaction_data.columns]
    assert sum(cytokine_flag) > 0
    interaction_data = interaction_data.loc[:, cytokine_flag]
    interaction_data = interaction_data.loc[interaction_data.isnull().mean(axis=1) < 1]

    # extract t-values
    t_flag = [(v.find('t.') == 0) for v in interaction_data.columns]
    assert sum(t_flag) > 0
    interaction_t = interaction_data.loc[:, t_flag]
    interaction_t.columns = ["%s.%s" % (dataset_name, '.'.join(v.split('.')[1:-1])) for v in interaction_t.columns]
    assert interaction_t.columns.value_counts().max() == 1  # if sample ID extraction is correct, must have no redundancy

    # extract q-values
    q_flag = [(v.find('q.') == 0) for v in interaction_data.columns]
    assert sum(q_flag) > 0
    interaction_q = interaction_data.loc[:, q_flag]
    interaction_q.columns = ["%s.%s" % (dataset_name, '.'.join(v.split('.')[1:-1])) for v in interaction_q.columns]
    assert interaction_q.columns.value_counts().max() == 1  # if sample ID extraction is correct, must have no redundancy

    dataset_signature = pd.DataFrame()
    qthres = 0.05
    if cytokine_signature_version == 0: # 不过滤
        dataset_signature = interaction_t
    elif cytokine_signature_version == 1: # 过滤q值小于0.05的比例（在所有基因中）小于0.001的样本
        frac_thres = 0.001
        # filter the rate of q-value < 0.05 smaller than frac_thres
        qvalue_flag = (interaction_q < qthres).mean() > frac_thres
        if qvalue_flag.sum() == 0:  # nothing to include
            return None
        dataset_signature = interaction_t.loc[:, qvalue_flag]
    elif cytokine_signature_version == 2: # 只保留q值小于0.05的样本
        # filter the q-value < 0.05
        interaction_t[interaction_q > qthres] = pd.NA
        dataset_signature = interaction_t.dropna(how='all')

    return dataset_signature

cytokine_signature_dict = {}
for interaction_filename in tqdm(interaction_list, desc='Processing interaction list'):
    interaction_filepath = os.path.join(interaction_path, interaction_filename)
    interaction_data = pd.read_csv(interaction_filepath, index_col=0)
    dataset_name = os.path.basename(interaction_filepath).split('.')[0]

    for cytokine in all_cytokine_list:
        dataset_cytokine_signature = get_dataset_cytokine_signature(dataset_name, interaction_data, cytokine)
        if dataset_cytokine_signature is None or len(dataset_cytokine_signature) < 1: continue
        if cytokine not in cytokine_signature_dict.keys():
            cytokine_signature_dict[cytokine] = [dataset_cytokine_signature]
        else:
            cytokine_signature_dict[cytokine].append(dataset_cytokine_signature)

cytokine_signature_dir = os.path.join(output_file_directory, f'cytokine_signature')
if not os.path.exists(cytokine_signature_dir):
    os.mkdir(cytokine_signature_dir)
for cytokine in tqdm(all_cytokine_list, desc='Processing cytokine list'):
    if cytokine not in cytokine_signature_dict.keys():
        continue
    cytokine_signature = pd.concat(cytokine_signature_dict[cytokine], axis=1, join='outer', sort=False)
    assert cytokine_signature.columns.value_counts().max() == 1

    cytokine_signature_filename = os.path.join(cytokine_signature_dir, f'{cytokine}_{cytokine_signature_version}.signature.csv')
    cytokine_signature.to_csv(cytokine_signature_filename)
print("Process end!")