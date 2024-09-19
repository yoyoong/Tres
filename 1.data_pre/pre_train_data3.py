import argparse
import os
import sys
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from tqdm.autonotebook import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('-E', "--expression_path", type=str, required=False, help="Gene expression file.",
                    default='/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_count_data/PAAD_CRA001160.csv')
parser.add_argument('-CTR', "--celltype_mapping_rules_file", type=str, required=False, help="Celltype mapping rules file, .txt format",
                    default='/sibcb1/bioinformatics/hongyuyang/dataset/Tres/0.model_file/tisch_celltype_mapping_rule.txt')
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_count_data2')
parser.add_argument('-O', "--output_tag", type=str, required=False, help="Prefix for output files.",
                    default='PAAD_CRA001160')
args = parser.parse_args()

expression_path = args.expression_path
celltype_mapping_rules_file = args.celltype_mapping_rules_file
output_file_directory = args.output_file_directory
output_tag = args.output_tag

if os.path.isdir(expression_path):
    print(f"This dataset don't need process.")
    sys.exit(1)
gem_df = pd.read_csv(expression_path, index_col=0, header=0)

# get cell name list
CellMetainfo_file = os.path.join(f'/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_data/0.raw_data/{output_tag}_CellMetainfo_table.tsv')
CellMetainfo_df = pd.read_csv(CellMetainfo_file, header=0, sep='\t')

mapping_rules_df = pd.read_csv(celltype_mapping_rules_file, sep='\t')
mapping_rules_dict = {key: value for key, value in zip(mapping_rules_df.iloc[:, 0], mapping_rules_df.iloc[:, 1])}

def generate_cellname(cell, celltype, sample=None):
    if celltype in mapping_rules_dict.keys():
        celltype = mapping_rules_dict[celltype]
    if sample is None:
        if "@" in cell:
            sample = str(cell).split("@")[0]
            cell = str(cell).split("@")[1]
        elif "." in cell:
            sample = str(cell).split(".")[1]
            cell = str(cell).split(".")[0]
        elif "_" in cell:
            sample = str(cell).split("_")[1]
            cell = str(cell).split("_")[0]
        return celltype + '.' + str(sample) + '.' + str(cell)
    else:
        return celltype + '.' + str(sample) + '.' + str(cell)


if 'Patient' in CellMetainfo_df:
    CellMetainfo_df['cellname'] = CellMetainfo_df.apply(
        lambda x: generate_cellname(x['Cell'], x['Celltype (major-lineage)'], x['Patient']), axis=1)
elif 'Sample' in CellMetainfo_df:
    CellMetainfo_df['cellname'] = CellMetainfo_df.apply(
        lambda x: generate_cellname(x['Cell'], x['Celltype (major-lineage)'], x['Sample']), axis=1)
else:
    CellMetainfo_df['cellname'] = CellMetainfo_df.apply(
        lambda x: generate_cellname(x['Cell'], x['Celltype (major-lineage)']), axis=1)
cellname_list = list(CellMetainfo_df['cellname'])

gem_df.columns = cellname_list
gem_dir = os.path.join(output_file_directory, output_tag)
if not os.path.exists(gem_dir):
    os.makedirs(gem_dir)
else:
    file_list = os.listdir(gem_dir)
    for file in file_list:
        file_path = os.path.join(gem_dir, file)
        os.remove(file_path)

flag_group = [v.split('.')[0] for v in gem_df.columns]
gem_group = gem_df.groupby(flag_group, axis=1)
for celltype, gem_sub in tqdm(gem_group, desc="Celltype data_process"):
    celltype = celltype.replace("/", "_")
    gem_sub_path = os.path.join(gem_dir, f'{output_tag}.{celltype}.csv')
    gem_sub.to_csv(gem_sub_path)

gem_path = os.path.join(output_file_directory, f'{output_tag}.csv')
gem_df.to_csv(gem_path)
print(f"{output_tag} process end")
