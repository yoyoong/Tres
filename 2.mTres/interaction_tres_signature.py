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
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/1.paper_data/4.Interaction/dataset_interaction')
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/1.paper_data/4.Interaction')
parser.add_argument('-O', "--output_tag", type=str, required=False, help="Prefix for output files.", default='Tres_signature')
parser.add_argument('-C', "--cytokine_info", type=str, required=False, help="Name of signaling, str or .txt file"
                    , default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/1.paper_data/4.Interaction/cytokine_info.txt')
args = parser.parse_args()

interaction_path = args.interaction_path
output_file_directory = args.output_file_directory
output_tag = args.output_tag
cytokine_info = args.cytokine_info

cytokine_info_df = pd.read_csv(cytokine_info, index_col=0)
cytokine_positive = cytokine_info_df[cytokine_info_df['flag'] == '+'].index.values.tolist()
cytokine_negative = cytokine_info_df[cytokine_info_df['flag'] == '-'].index.values.tolist()

interaction_positive_filename = os.path.join(output_file_directory, f'interaction.positive.csv')
interaction_negative_filename = os.path.join(output_file_directory, f'interaction.negative.csv')

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

interaction_positive = pd.read_csv(interaction_positive_filename, index_col=0, header=0)
interaction_negative = pd.read_csv(interaction_negative_filename, index_col=0, header=0)
all_median_positive = pd.DataFrame(columns=['GeneID', 'Tres'])
all_median_negative = pd.DataFrame(columns=['GeneID', 'Tres'])
positive_sample_num = len(np.unique(interaction_positive['SampleID']))
negative_sample_num = len(np.unique(interaction_negative['SampleID']))
interaction_positive_group1 = interaction_positive.groupby(['GeneID'])
interaction_negative_group1 = interaction_negative.groupby(['GeneID'])
median_num_cut = 10

for group_key, group_data1 in tqdm(interaction_positive_group1, desc='Processing positive'):
    interaction_positive_group2 = group_data1.groupby(['SampleID'])
    sample_median = []
    for SampleID, group_data2 in interaction_positive_group2:
        sample_median.append(group_data2['t'].median())
    if len(sample_median) >= median_num_cut:
        all_median = np.median(sample_median)
        new_row = pd.Series({'GeneID': group_key[0], 'Tres': all_median})
        all_median_positive = pd.concat([all_median_positive, new_row.to_frame().T])
for group_key, group_data1 in tqdm(interaction_negative_group1, desc='Processing negative'):
    interaction_negative_group2 = group_data1.groupby(['SampleID'])
    sample_median = []
    for SampleID, group_data2 in interaction_negative_group2:
        sample_median.append(group_data2['t'].median())
    if len(sample_median) >= median_num_cut:
        all_median = np.median(sample_median)
        new_row = pd.Series({'GeneID': group_key[0], 'Tres': all_median})
        all_median_negative = pd.concat([all_median_negative, new_row.to_frame().T])

all_median_positive.set_index('GeneID', inplace=True)
all_median_negative.set_index('GeneID', inplace=True)
tres_signature_positive_filename = os.path.join(output_file_directory, f'tres_signature.positive.csv')
tres_signature_negative_filename = os.path.join(output_file_directory, f'tres_signature.negative.csv')
all_median_positive.to_csv(tres_signature_positive_filename)
all_median_negative.to_csv(tres_signature_negative_filename)
print("Process end!")