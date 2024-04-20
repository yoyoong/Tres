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
parser.add_argument('-I', "--interaction_path", type=str, required=False, help="Interaction result path.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-2.Macrophage_Interaction/dataset_interaction')
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-2.Macrophage_Interaction')
parser.add_argument('-O', "--output_tag", type=str, required=False, help="Prefix for output files.", default='Tres_signature')
parser.add_argument('-C', "--cytokine_info", type=str, required=False, help="Name of signaling, str or .txt file"
                    , default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-2.Macrophage_Interaction/cytokine_info.txt')
args = parser.parse_args()

interaction_path = args.interaction_path
output_file_directory = args.output_file_directory
output_tag = args.output_tag
cytokine_info = args.cytokine_info

def get_dataset_cytokine_signature(dataset_name, interaction_data, cytokine):
    qthres = 0.05
    frac_thres = 0.001

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

    # filter the q-value < 0.05
    # interaction_t[interaction_q > qthres] = pd.NA
    # dataset_signature = interaction_t.dropna(how='all')

    # filter the rate of q-value < 0.05 smaller than frac_thres
    qvalue_flag = (interaction_q < qthres).mean() > frac_thres
    if qvalue_flag.sum() == 0:  # nothing to include
        return None
    dataset_signature = interaction_t.loc[:, qvalue_flag]

    return dataset_signature

def get_median_signature(all_signature, pn_flag, filter_flag):
    # filter the sample that are profiled by all signals
    common_col = None
    for signature in all_signature:
        if common_col is None:
            common_col = signature.columns
        else:
            common_col = common_col.intersection(signature.columns)
    for i, result in enumerate(all_signature):
        all_signature[i] = result.loc[:, common_col]

    # create median signature across immuno-suppressive signals
    median_signature = pd.concat(all_signature, axis=1, join='inner')
    median_signature = median_signature.groupby(median_signature.columns, axis=1).median()
    # output to file
    median_signature_filename = os.path.join(output_file_directory, f'Median_signature_{filter_flag}.{pn_flag}.csv')
    median_signature.to_csv(median_signature_filename)

    # filter the rate of valid median sample num < 0.5
    if filter_flag == 1:
        tres_signature = median_signature.dropna(how='all').median(axis=1)
    elif filter_flag == 2:
        tres_signature = median_signature.loc[median_signature.isnull().mean(axis=1) < 0.5].median(axis=1)
    elif filter_flag == 3:
        tres_signature = median_signature.dropna().median(axis=1)

    tres_signature.name = 'Tres'

    return tres_signature

cytokine_info_df = pd.read_csv(cytokine_info, index_col=0)
cytokine_list = cytokine_info_df.index.values.tolist()

cytokine_signature_dict = {}
interaction_list = sorted(os.listdir(interaction_path))
for interaction_filename in tqdm(interaction_list, desc='Processing interaction list'):
    interaction_filepath = os.path.join(interaction_path, interaction_filename)
    interaction_data = pd.read_csv(interaction_filepath, index_col=0)
    dataset_name = os.path.basename(interaction_filepath).split('.')[0]

    for cytokine in cytokine_list:
        dataset_cytokine_signature = get_dataset_cytokine_signature(dataset_name, interaction_data, cytokine)
        if dataset_cytokine_signature is None or len(dataset_cytokine_signature) < 1: continue
        if cytokine not in cytokine_signature_dict.keys():
            cytokine_signature_dict[cytokine] = [dataset_cytokine_signature]
        else:
            cytokine_signature_dict[cytokine].append(dataset_cytokine_signature)

with open(os.path.join(output_file_directory, f'cytokine_signature_dict.pkl'), 'wb') as file:
    pickle.dump(cytokine_signature_dict, file)

# with open(os.path.join(output_file_directory, f'cytokine_signature_dict.pkl'), 'rb') as file:
#     cytokine_signature_dict = pd.read_pickle(file)

all_positive_signature = []
all_negative_signature = []
for cytokine in tqdm(cytokine_list, desc='Processing cytokine list'):
    if cytokine not in cytokine_signature_dict.keys():
        continue
    cytokine_signature = pd.concat(cytokine_signature_dict[cytokine], axis=1, join='outer', sort=False)
    assert cytokine_signature.columns.value_counts().max() == 1

    cytokine_signature_filename = os.path.join(output_file_directory, 'cytokine_signature', f'{cytokine}.signature.csv')
    cytokine_signature.to_csv(cytokine_signature_filename)

    if cytokine_info_df.loc[cytokine]['flag'] == '+':
        all_positive_signature.append(cytokine_signature)
    else:
        all_negative_signature.append(cytokine_signature)

for filter_flag in [1, 2, 3]:
    positive_tres_signature = get_median_signature(all_positive_signature, 'positive', filter_flag)
    positive_tres_signature.to_csv(os.path.join(output_file_directory, f'Tres_signature_{filter_flag}.positive.csv'))
    negative_tres_signature = get_median_signature(all_negative_signature, 'negative', filter_flag)
    negative_tres_signature.to_csv(os.path.join(output_file_directory, f'Tres_signature_{filter_flag}.negative.csv'))
print("Process end!")