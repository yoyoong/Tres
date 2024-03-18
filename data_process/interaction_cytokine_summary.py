import argparse
import numpy
import os
import pandas as pd
import sys
import warnings
from tqdm.autonotebook import tqdm
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument('-I', "--interaction_path", type=str, required=False, help="Interaction result path.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/4.Interaction/dataset_interaction')
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/4.Interaction/cytokine_summary')
parser.add_argument('-O', "--output_tag", type=str, required=False, help="Prefix for output files.", default='')
parser.add_argument('-C', "--cytokine_key", type=str, required=False, help="Name of signaling, str or .txt file"
                    , default='IL6') # if null, calculate all cytokine
args = parser.parse_args()

interaction_path = args.interaction_path
output_file_directory = args.output_file_directory
output_tag = args.output_tag
cytokine = args.cytokine_key


interaction_list = sorted(os.listdir(interaction_path))
cytokine_summary_df = pd.DataFrame(columns=['GeneID', 'SampleNum', 'Num(t>0,q<0.05)', 'Num(t<0,q<0.05)'])
for interaction_filename in tqdm(interaction_list, desc='Processing dataset'):
    interaction_filepath = os.path.join(interaction_path, interaction_filename)
    interaction_data = pd.read_csv(interaction_filepath, index_col=0)
    # geneid = interaction_data['GeneID'].unique()
    interaction_data_filtered = interaction_data[interaction_data.index == cytokine]
    interaction_data_grouped = interaction_data_filtered.groupby('GeneID')

    cytokine_summary_sub_df = pd.DataFrame(columns=['GeneID', 'SampleNum', 'Num(t>0,q<0.05)', 'Num(t<0,q<0.05)'])
    for interaction_group in interaction_data_grouped:
        GeneID = interaction_group[0]
        interaction_group_data = interaction_group[1]
        SampleNum = interaction_group_data.shape[0]
        Num1 = interaction_group_data[(interaction_group_data['t'] > 0) & (interaction_group_data['q'] < 0.05)].shape[0]
        Num2 = interaction_group_data[(interaction_group_data['t'] < 0) & (interaction_group_data['q'] < 0.05)].shape[0]

        new_row = pd.Series({'GeneID': GeneID, 'SampleNum': SampleNum, 'Num(t>0,q<0.05)': Num1, 'Num(t<0,q<0.05)': Num2})
        cytokine_summary_sub_df = pd.concat([cytokine_summary_sub_df, new_row.to_frame().T])

    cytokine_summary_df = pd.concat([cytokine_summary_df, cytokine_summary_sub_df])

cytokine_summary_df = cytokine_summary_df.groupby(by=['GeneID'])[['SampleNum', 'Num(t>0,q<0.05)', 'Num(t<0,q<0.05)']].sum()
cytokine_summary_filename = os.path.join(output_file_directory, f'{cytokine}.summary.csv')
cytokine_summary_df.to_csv(cytokine_summary_filename)
print("Process end!")