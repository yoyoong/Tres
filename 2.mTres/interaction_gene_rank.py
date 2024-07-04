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
parser.add_argument('-CT', "--celltype", type=str, default='NK', required=False, help="cell type")
parser.add_argument('-CIV', "--cytokine_info_version", type=int, default='0', required=False, help="cytokine info version")
args = parser.parse_args()

celltype = args.celltype
cytokine_info_version = args.cytokine_info_version
cytokine_info_file = f'/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.model_file/cytokine_info_{cytokine_info_version}.{celltype}.txt'
gene_annotation = '/sibcb2/bioinformatics/iGenome/STAR/GENCODE/human_hg38/ID/tx2g.txt'
if celltype == "CD8T" :
    output_file_directory = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-1.CD8T_Interaction'
elif celltype == "Macrophage" :
    output_file_directory = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-2.Macrophage_Interaction'
elif celltype == "Neutrophils" :
    output_file_directory = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-3.Neutrophils_Interaction'
elif celltype == "NK" :
    output_file_directory = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-4.NK_Interaction'
elif celltype == "NK_act" :
    output_file_directory = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-4-0.NK_act_Interaction'

cytokine_info_df = pd.read_csv(cytokine_info_file, index_col=0)
cytokine_list = cytokine_info_df.index.values.tolist()
gene_annotation_df = pd.read_csv(gene_annotation, index_col=0, header=0, delimiter='\t')

positive_gene_rank_df = pd.DataFrame(columns=['GeneID', 'SampleNum', 'Num(t>0,q<0.05)', 'Num(t<0,q<0.05)',
                                     'Rate(t>0,q<0.05)', 'Rate(t<0,q<0.05)', 'Rank(t>0,q<0.05)', 'Rank(t<0,q<0.05)', 'mType'])
negative_gene_rank_df = pd.DataFrame(columns=['GeneID', 'SampleNum', 'Num(t>0,q<0.05)', 'Num(t<0,q<0.05)',
                                     'Rate(t>0,q<0.05)', 'Rate(t<0,q<0.05)', 'Rank(t>0,q<0.05)', 'Rank(t<0,q<0.05)', 'mType'])
for cytokine in tqdm(cytokine_list, desc='Processing cytokine list'):
    cytokine_summary_filename = os.path.join(output_file_directory, 'cytokine_summary', f'{cytokine}.summary.csv')
    cytokine_summary = pd.read_csv(cytokine_summary_filename)

    cytokine_flag = cytokine_info_df.loc[cytokine]['flag']
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

gene_rank_dir = os.path.join(output_file_directory, 'gene_rank')
if not os.path.exists(gene_rank_dir):
    os.mkdir(gene_rank_dir)
negative_gene_rank_filename = os.path.join(gene_rank_dir, f'Gene_rank_{cytokine_info_version}.negative.csv')
positive_gene_rank_filename = os.path.join(gene_rank_dir, f'Gene_rank_{cytokine_info_version}.positive.csv')
negative_gene_rank_df.to_csv(negative_gene_rank_filename)
positive_gene_rank_df.to_csv(positive_gene_rank_filename)
print("Process end!")