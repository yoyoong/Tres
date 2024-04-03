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
parser.add_argument('-I', "--interaction_path", type=str, required=False, help="Interaction result path.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/4.Interaction/new_dataset_interaction')
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/4.Interaction')
parser.add_argument('-O', "--output_tag", type=str, required=False, help="Prefix for output files.", default='')
parser.add_argument('-C', "--cytokine_info", type=str, required=False, help="Name of signaling, str or .txt file"
                    , default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/4.Interaction/cytokine_info.txt')
parser.add_argument('-G', "--gene_annotation", type=str, required=False, help="Name of signaling, str or .txt file"
                    , default='/sibcb2/bioinformatics/iGenome/STAR/GENCODE/human_hg38/ID/tx2g.txt')
args = parser.parse_args()

interaction_path = args.interaction_path
output_file_directory = args.output_file_directory
output_tag = args.output_tag
cytokine_info = args.cytokine_info
gene_annotation = args.gene_annotation

cytokine_info_df = pd.read_csv(cytokine_info, index_col=0)
cytokine_list = cytokine_info_df.index.values.tolist()
gene_annotation_df = pd.read_csv(gene_annotation, index_col=0, header=0, delimiter='\t')

cytokine_summary_dict = {}
interaction_list = sorted(os.listdir(interaction_path))
for interaction_filename in tqdm(interaction_list, desc='Processing dataset'):
    interaction_filepath = os.path.join(interaction_path, interaction_filename)
    interaction_data = pd.read_csv(interaction_filepath, index_col=0)
    dataset_name = os.path.basename(interaction_filepath).split('.')[0]

    for cytokine in cytokine_list:
        # filter the cytokine
        cytokine_flag = [col for col in interaction_data.columns if col.split('.')[-1] == cytokine]
        interaction_data = interaction_data[cytokine_flag]
        interaction_data = interaction_data.dropna(how='all')

        # extract t-values
        t_flag = [(v.find('t.') == 0) for v in interaction_data.columns]
        interaction_t = interaction_data.loc[:, t_flag]
        interaction_t.columns = ["%s.%s" % (dataset_name, '.'.join(v.split('.')[1:-1])) for v in interaction_t.columns]

        # extract q-values
        q_flag = [(v.find('q.') == 0) for v in interaction_data.columns]
        interaction_q = interaction_data.loc[:, q_flag]
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
            continue

        if cytokine not in cytokine_summary_dict.keys():
            cytokine_summary_dict[cytokine] = cytokine_summary_df
        else:
            # new_df = (pd.concat([cytokine_summary_dict[cytokine], cytokine_summary_df], axis=0, join='outer')
            #           .groupby('GeneID').sum())
            new_df = cytokine_summary_dict[cytokine].add(cytokine_summary_df, fill_value=0)
            cytokine_summary_dict[cytokine] = new_df

positive_gene_rank_df = pd.DataFrame(columns=['GeneID', 'SampleNum', 'Num(t>0,q<0.05)', 'Num(t<0,q<0.05)',
                                     'Rate(t>0,q<0.05)', 'Rate(t<0,q<0.05)', 'Rank(t>0,q<0.05)', 'Rank(t<0,q<0.05)', 'mType'])
negative_gene_rank_df = pd.DataFrame(columns=['GeneID', 'SampleNum', 'Num(t>0,q<0.05)', 'Num(t<0,q<0.05)',
                                     'Rate(t>0,q<0.05)', 'Rate(t<0,q<0.05)', 'Rank(t>0,q<0.05)', 'Rank(t<0,q<0.05)', 'mType'])
for cytokine in tqdm(cytokine_list, desc='Processing cytokine list'):
    cytokine_summary = cytokine_summary_dict[cytokine]
    cytokine_summary_filename = os.path.join(output_file_directory, 'new_cytokine_summary', f'{cytokine}.summary.csv')
    cytokine_summary.to_csv(cytokine_summary_filename)

    if cytokine_flag == "+":
        positive_gene_rank_df = pd.concat([positive_gene_rank_df, cytokine_summary], ignore_index=True)
    else:
        negative_gene_rank_df = pd.concat([negative_gene_rank_df, cytokine_summary], ignore_index=True)

# merge the num
positive_gene_rank_df = positive_gene_rank_df.groupby(by=['GeneID'])[['SampleNum', 'Num(t>0,q<0.05)', 'Num(t<0,q<0.05)']].sum()
negative_gene_rank_df = negative_gene_rank_df.groupby(by=['GeneID'])[['SampleNum', 'Num(t>0,q<0.05)', 'Num(t<0,q<0.05)']].sum()

# filter sample num
positive_gene_rank_df = positive_gene_rank_df[positive_gene_rank_df['SampleNum'] >= 100]
negative_gene_rank_df = negative_gene_rank_df[negative_gene_rank_df['SampleNum'] >= 100]

# get rate
positive_gene_rank_df["Rate(t>0,q<0.05)"] = positive_gene_rank_df["Num(t>0,q<0.05)"] / positive_gene_rank_df["SampleNum"]
positive_gene_rank_df["Rate(t<0,q<0.05)"] = positive_gene_rank_df["Num(t<0,q<0.05)"] / positive_gene_rank_df["SampleNum"]
negative_gene_rank_df["Rate(t>0,q<0.05)"] = negative_gene_rank_df["Num(t>0,q<0.05)"] / negative_gene_rank_df["SampleNum"]
negative_gene_rank_df["Rate(t<0,q<0.05)"] = negative_gene_rank_df["Num(t<0,q<0.05)"] / negative_gene_rank_df["SampleNum"]

# get rank
positive_gene_rank_df["Rank(t>0,q<0.05)"] = scipy.stats.rankdata(positive_gene_rank_df["Rate(t>0,q<0.05)"] * -1, method='min')
positive_gene_rank_df["Rank(t<0,q<0.05)"] = scipy.stats.rankdata(positive_gene_rank_df["Rate(t<0,q<0.05)"] * -1, method='min')
negative_gene_rank_df["Rank(t>0,q<0.05)"] = scipy.stats.rankdata(negative_gene_rank_df["Rate(t>0,q<0.05)"] * -1, method='min')
negative_gene_rank_df["Rank(t<0,q<0.05)"] = scipy.stats.rankdata(negative_gene_rank_df["Rate(t<0,q<0.05)"] * -1, method='min')

# get type
positive_gene_list = list(set(gene_annotation_df['Symbol']).intersection(positive_gene_rank_df.index))
negative_gene_list = list(set(gene_annotation_df['Symbol']).intersection(negative_gene_rank_df.index))
def merge_set_to_string(row):
    return '/'.join(map(str, row))
gene_mType = gene_annotation_df.groupby('Symbol')['mType'].agg(set).apply(merge_set_to_string)
positive_gene_rank_df.loc[positive_gene_list, 'mType'] = gene_mType.loc[positive_gene_list]
negative_gene_rank_df.loc[negative_gene_list, 'mType'] = gene_mType.loc[negative_gene_list]

negative_gene_rank_filename = os.path.join(output_file_directory, f'New_Gene_rank.negative.csv')
positive_gene_rank_filename = os.path.join(output_file_directory, f'New_Gene_rank.positive.csv')
negative_gene_rank_df.to_csv(negative_gene_rank_filename)
positive_gene_rank_df.to_csv(positive_gene_rank_filename)
print("Process end!")