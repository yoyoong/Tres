import argparse
import numpy
import os
import pandas as pd
import sys
import warnings
from tqdm.autonotebook import tqdm
warnings.filterwarnings("ignore")
import scipy

parser = argparse.ArgumentParser()
parser.add_argument('-CT', "--celltype", type=str, default='CD8T', required=False, help="cell type")
args = parser.parse_args()

celltype = args.celltype
model_matrix_file = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.model_file/signature.centroid.expand'
gene_annotation = '/sibcb2/bioinformatics/iGenome/STAR/GENCODE/human_hg38/ID/tx2g.txt'

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

model_matrix_file = pd.read_csv(model_matrix_file, sep='\t', index_col=0)
all_cytokine_list = model_matrix_file.columns.values.tolist()

interaction_list = []
if os.path.isdir(interaction_path):
    interaction_list = sorted(os.listdir(interaction_path))
else:
    interaction_list.append(interaction_path)

cytokine_summary_dict = {}
for interaction_filename in tqdm(interaction_list, desc='Processing dataset'):
    interaction_filepath = os.path.join(interaction_path, interaction_filename)
    interaction_data = pd.read_csv(interaction_filepath, index_col=0)
    dataset_name = os.path.basename(interaction_filepath).split('.')[0]

    for cytokine in all_cytokine_list:
        # filter the cytokine
        cytokine_flag = [col for col in interaction_data.columns if col.split('.')[-1] == cytokine]
        interaction_cytokine_data = interaction_data[cytokine_flag].dropna(how='all')

        # extract t-values
        t_flag = [(v.find('t.') == 0) for v in interaction_cytokine_data.columns]
        interaction_t = interaction_cytokine_data.loc[:, t_flag]
        interaction_t.columns = ["%s.%s" % (dataset_name, '.'.join(v.split('.')[1:-1])) for v in interaction_t.columns]

        # extract q-values
        q_flag = [(v.find('q.') == 0) for v in interaction_cytokine_data.columns]
        interaction_q = interaction_cytokine_data.loc[:, q_flag]
        interaction_q.columns = ["%s.%s" % (dataset_name, '.'.join(v.split('.')[1:-1])) for v in interaction_q.columns]

        # summary the interaction
        cytokine_summary_df = pd.DataFrame(columns=['GeneID', 'SampleNum', 'Num(t>0,q<0.05)', 'Num(t<0,q<0.05)'])
        cytokine_summary_df['GeneID'] = interaction_t.index
        cytokine_summary_df.set_index('GeneID', inplace=True)
        cytokine_summary_df['SampleNum'] = interaction_t.apply(lambda row: row.count(), axis=1)
        interaction_t[interaction_q > 0.05] = pd.NA # mask the q value > 0.05
        cytokine_summary_df['Num(t>0,q<0.05)'] = interaction_t.apply(lambda row: sum(row > 0), axis=1)
        cytokine_summary_df['Num(t<0,q<0.05)'] = interaction_t.apply(lambda row: sum(row < 0), axis=1)
        if cytokine_summary_df.empty:
            print(f'{interaction_filename}-{cytokine} empty')
            continue

        cytokine_summary_df.sort_index(inplace=True)
        if cytokine not in cytokine_summary_dict.keys():
            cytokine_summary_dict[cytokine] = cytokine_summary_df
        else:
            new_df = cytokine_summary_dict[cytokine].add(cytokine_summary_df, fill_value=0)
            cytokine_summary_dict[cytokine] = new_df

# save the cytokine summary
cytokine_summary_dir = os.path.join(output_file_directory, 'cytokine_summary')
if not os.path.exists(cytokine_summary_dir):
    os.mkdir(cytokine_summary_dir)
for cytokine, cytokine_summary in cytokine_summary_dict.items():
    cytokine_summary_filename = os.path.join(cytokine_summary_dir, f'{cytokine}.summary.csv')
    cytokine_summary.to_csv(cytokine_summary_filename)