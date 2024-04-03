import argparse
import numpy
import os
import pandas as pd
import sys
import warnings
from tqdm.autonotebook import tqdm
import scipy
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument('-I', "--summary_path", type=str, required=False, help="Summary data path.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/4.Interaction/cytokine_summary')
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/4.Interaction')
parser.add_argument('-O', "--output_tag", type=str, required=False, help="Prefix for output files.", default='Gene_rank')
parser.add_argument('-C', "--cytokine_info", type=str, required=False, help="Name of signaling, str or .txt file"
                    , default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/4.Interaction/cytokine_info.txt')
parser.add_argument('-G', "--gene_annotation", type=str, required=False, help="Name of signaling, str or .txt file"
                    , default='/sibcb2/bioinformatics/iGenome/STAR/GENCODE/human_hg38/ID/tx2g.txt')
args = parser.parse_args()

summary_path = args.summary_path
output_file_directory = args.output_file_directory
output_tag = args.output_tag
cytokine_info = args.cytokine_info
gene_annotation = args.gene_annotation

cytokine_info_df = pd.read_csv(cytokine_info, index_col=0)
gene_annotation_df = pd.read_csv(gene_annotation, index_col=0, header=0, delimiter='\t')

summary_list = sorted(os.listdir(summary_path))
positive_gene_rank_df = pd.DataFrame(columns=['GeneID', 'SampleNum', 'Num(t>0,q<0.05)', 'Num(t<0,q<0.05)',
                                     'Rate(t>0,q<0.05)', 'Rate(t<0,q<0.05)', 'Rank(t>0,q<0.05)', 'Rank(t<0,q<0.05)', 'mType'])
negative_gene_rank_df = pd.DataFrame(columns=['GeneID', 'SampleNum', 'Num(t>0,q<0.05)', 'Num(t<0,q<0.05)',
                                     'Rate(t>0,q<0.05)', 'Rate(t<0,q<0.05)', 'Rank(t>0,q<0.05)', 'Rank(t<0,q<0.05)', 'mType'])
for summary_filename in summary_list:
    cytokine = summary_filename.split(".")[0]
    cytokine_flag = cytokine_info_df.loc[cytokine]['flag']
    cytokine_summary_path = os.path.join(summary_path, summary_filename)
    cytokine_summary_data = pd.read_csv(cytokine_summary_path)
    
    if cytokine_flag == "+":
        positive_gene_rank_df = pd.concat([positive_gene_rank_df, cytokine_summary_data], ignore_index=True)
    else:
        negative_gene_rank_df = pd.concat([negative_gene_rank_df, cytokine_summary_data], ignore_index=True)

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

negative_gene_rank_filename = os.path.join(output_file_directory, f'Gene_rank.negative.csv')
positive_gene_rank_filename = os.path.join(output_file_directory, f'Gene_rank.positive.csv')
negative_gene_rank_df.to_csv(negative_gene_rank_filename)
positive_gene_rank_df.to_csv(positive_gene_rank_filename)
print("Process end!")