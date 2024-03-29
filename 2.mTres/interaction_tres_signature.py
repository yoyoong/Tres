import argparse
import numpy
import os

import numpy as np
import pandas as pd
import sys
import warnings
from tqdm.autonotebook import tqdm
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument('-I', "--interaction_path", type=str, required=False, help="Interaction result path.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/4.Interaction/new_dataset_interaction')
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/4.Interaction')
parser.add_argument('-O', "--output_tag", type=str, required=False, help="Prefix for output files.", default='Tres_signature')
parser.add_argument('-C', "--cytokine_info", type=str, required=False, help="Name of signaling, str or .txt file"
                    , default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/4.Interaction/cytokine_info.txt')
args = parser.parse_args()

interaction_path = args.interaction_path
output_file_directory = args.output_file_directory
output_tag = args.output_tag
cytokine_info = args.cytokine_info

def get_cytokine_signature(interaction_path, cytokine):
    qthres = 0.05
    frac_thres = 1e-3
    cytokine_signature = []
    interaction_list = sorted(os.listdir(interaction_path))
    for interaction_filename in tqdm(interaction_list, desc='Processing interaction list'):
        interaction_filepath = os.path.join(interaction_path, interaction_filename)
        interaction_data = pd.read_csv(interaction_filepath, index_col=0)

        # filter the cytokine
        cytokine_flag = [(v.split('.')[-1] == cytokine) for v in interaction_data.columns]
        assert sum(cytokine_flag) > 0
        interaction_data = interaction_data.loc[:, cytokine_flag]
        interaction_data = interaction_data.loc[interaction_data.isnull().mean(axis=1) < 1]

        # extract t-values
        t_flag = [(v.find('t.') == 0) for v in interaction_data.columns]
        assert sum(t_flag) > 0
        interaction_t = interaction_data.loc[:, t_flag]
        interaction_t.columns = ['.'.join(v.split('.')[1:-1]) for v in interaction_t.columns]
        assert interaction_t.columns.value_counts().max() == 1  # if sample ID extraction is correct, must have no redundancy

        # extract q-values
        q_flag = [(v.find('q.') == 0) for v in interaction_data.columns]
        assert sum(q_flag) > 0
        interaction_q = interaction_data.loc[:, q_flag]
        interaction_q.columns = ['.'.join(v.split('.')[1:-1]) for v in interaction_q.columns]
        assert interaction_q.columns.value_counts().max() == 1  # if sample ID extraction is correct, must have no redundancy

        # filter the rate of q-value < 0.05 smaller than frac_thres
        qvalue_flag = (interaction_q < qthres).mean() > frac_thres
        if qvalue_flag.sum() == 0:  # nothing to include
            continue
        dataset_signature = interaction_t.loc[:, qvalue_flag]

        # append to cytokine signature
        cytokine_signature.append(dataset_signature)

    cytokine_signature = pd.concat(cytokine_signature, axis=1, join='outer', sort=False)
    assert cytokine_signature.columns.value_counts().max() == 1
    cytokine_signature_filename = os.path.join(output_file_directory, 'cytokine_signature', f'{cytokine}.signature.csv')
    cytokine_signature.to_csv(cytokine_signature_filename)

    return cytokine_signature

def get_median_signature(all_signature, cytokine_list):
    # filter the cell that are profiled by all signals
    common_col = None
    for signature in all_signature:
        if common_col is None:
            common_col = signature.columns
        else:
            common_col = common_col.intersection(signature.columns)
    for i, result in enumerate(all_signature):
        all_signature[i] = result.loc[:, common_col]

    # create median signature across three immuno-suppressive signals
    median_signature = pd.concat(all_signature, axis=1, join='inner')
    cnt_map = median_signature.columns.value_counts()
    assert sum(cnt_map != len(cytokine_list)) == 0
    median_signature = median_signature.groupby(median_signature.columns, axis=1).median()
    # output to file
    median_signature_filename = os.path.join(output_file_directory, f'Median.signature.csv')
    median_signature.to_csv(median_signature_filename)

    # filter the rate of valid median sample num < 0.5
    tres_signature = median_signature.loc[median_signature.isnull().mean(axis=1) < 0.5]
    tres_signature = tres_signature.dropna().median(axis=1)
    # output to file
    tres_signature_filename = os.path.join(output_file_directory, f'Tres.signature.csv')
    tres_signature.to_csv(tres_signature_filename)

cytokine_info_df = pd.read_csv(cytokine_info, index_col=0)
cytokine_positive = cytokine_info_df[cytokine_info_df['flag'] == '+'].index.values.tolist()
cytokine_negative = cytokine_info_df[cytokine_info_df['flag'] == '-'].index.values.tolist()

all_positive_signature = []
for cytokine in tqdm(cytokine_positive, desc='Processing positive cytokine'):
    cytokine_signature = get_cytokine_signature(interaction_path, cytokine)
    all_positive_signature.append(cytokine_signature)
all_negative_signature = []
for cytokine in tqdm(cytokine_negative, desc='Processing negative cytokine'):
    cytokine_signature = get_cytokine_signature(interaction_path, cytokine)
    all_negative_signature.append(cytokine_signature)


get_median_signature(all_positive_signature, cytokine_positive)
get_median_signature(all_negative_signature, cytokine_negative)

#
# interaction_positive_filename = os.path.join(output_file_directory, f'interaction.positive.csv')
# interaction_negative_filename = os.path.join(output_file_directory, f'interaction.negative.csv')
#
# interaction_positive = pd.DataFrame()
# interaction_negative = pd.DataFrame()
# interaction_list = sorted(os.listdir(interaction_path))
# for interaction_filename in tqdm(interaction_list, desc='Processing dataset'):
#     interaction_filepath = os.path.join(interaction_path, interaction_filename)
#     interaction_data = pd.read_csv(interaction_filepath, index_col=0)
#
#     # filter cytokine
#     interaction_positive_filtered = interaction_data[interaction_data.index.isin(cytokine_positive)]
#     interaction_negative_filtered = interaction_data[interaction_data.index.isin(cytokine_negative)]
#
#     # filter q < 0.05
#     interaction_positive_filtered = interaction_positive_filtered[interaction_positive_filtered['q'] < 0.05]
#     interaction_negative_filtered = interaction_negative_filtered[interaction_negative_filtered['q'] < 0.05]
#
#     interaction_positive = pd.concat([interaction_positive, interaction_positive_filtered])
#     interaction_negative = pd.concat([interaction_negative, interaction_negative_filtered])
# interaction_positive.to_csv(interaction_positive_filename)
# interaction_negative.to_csv(interaction_negative_filename)
#
# interaction_positive = pd.read_csv(interaction_positive_filename, index_col=0, header=0)
# interaction_negative = pd.read_csv(interaction_negative_filename, index_col=0, header=0)
# all_median_positive = pd.DataFrame(columns=['GeneID', 'Tres'])
# all_median_negative = pd.DataFrame(columns=['GeneID', 'Tres'])
# positive_sample_num = len(np.unique(interaction_positive['SampleID']))
# negative_sample_num = len(np.unique(interaction_negative['SampleID']))
# interaction_positive_group1 = interaction_positive.groupby(['GeneID'])
# interaction_negative_group1 = interaction_negative.groupby(['GeneID'])
# median_num_cut = 10
#
# for group_key, group_data1 in tqdm(interaction_positive_group1, desc='Processing positive'):
#     interaction_positive_group2 = group_data1.groupby(['SampleID'])
#     sample_median = []
#     for SampleID, group_data2 in interaction_positive_group2:
#         sample_median.append(group_data2['t'].median())
#     if len(sample_median) >= median_num_cut:
#         all_median = np.median(sample_median)
#         new_row = pd.Series({'GeneID': group_key[0], 'Tres': all_median})
#         all_median_positive = pd.concat([all_median_positive, new_row.to_frame().T])
# for group_key, group_data1 in tqdm(interaction_negative_group1, desc='Processing negative'):
#     interaction_negative_group2 = group_data1.groupby(['SampleID'])
#     sample_median = []
#     for SampleID, group_data2 in interaction_negative_group2:
#         sample_median.append(group_data2['t'].median())
#     if len(sample_median) >= median_num_cut:
#         all_median = np.median(sample_median)
#         new_row = pd.Series({'GeneID': group_key[0], 'Tres': all_median})
#         all_median_negative = pd.concat([all_median_negative, new_row.to_frame().T])
#
# all_median_positive.set_index('GeneID', inplace=True)
# all_median_negative.set_index('GeneID', inplace=True)
# tres_signature_positive_filename = os.path.join(output_file_directory, f'tres_signature.positive.csv')
# tres_signature_negative_filename = os.path.join(output_file_directory, f'tres_signature.negative.csv')
# all_median_positive.to_csv(tres_signature_positive_filename)
# all_median_negative.to_csv(tres_signature_negative_filename)
print("Process end!")