import argparse
import os
import numpy as np
import pandas as pd
import h5py
import scanpy as sc
from scipy.sparse import coo_matrix
from tqdm.autonotebook import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('-I', "--input", type=str, required=False, help="input file path",
                    default='/sibcb1/bioinformatics/fengyan/Project/Mus/scRNAseq/Rawdata/R23405568-NC1')
parser.add_argument('-A', "--annotation", type=str, required=False, help="annotation file path",
                    default='/sibcb1/bioinformatics/fengyan/Project/Mus/Tres/annotation/Prepare/NC1_annotation.txt')
parser.add_argument('-CTR', "--celltype_mapping_rules_file", type=str, required=False, help="Celltype mapping rules file, .txt format",
                    default='/sibcb1/bioinformatics/hongyuyang/dataset/Tres/0.model_file/tisch_celltype_mapping_rule.txt')
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb1/bioinformatics/hongyuyang/code/Tres')
parser.add_argument('-O', "--output_tag", type=str, required=False, help="Prefix for output files.",
                    default='test')
args = parser.parse_args()

input_path = args.input
annotation_path = args.annotation
celltype_mapping_rules_file = args.celltype_mapping_rules_file
output_file_directory = args.output_file_directory
output_tag = args.output_tag

genename_list = []
cellname_list = []

if os.path.isfile(input_path) and input_path.endswith('.h5'):
    h5file = h5py.File(input_path, "r")
    matrix = h5file['matrix']
    data_raw = np.array(matrix['data'])
    shape = np.array(matrix['shape'])
    indices = np.array(matrix['indices'])
    indptr = np.array(matrix['indptr'])

    # get gene name list
    features = matrix['features']
    gene_names = list(features['name'])
    genename_list = [gene.decode('utf-8') for gene in gene_names]

    # get cell name list
    CellMetainfo_file = input_path.replace('expression', 'CellMetainfo_table').replace('.h5', '.tsv')
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

    # convert to dense matrix
    row = indices
    col = np.zeros_like(indices)
    for index in range(1, len(indptr)):
        if index + 1 < len(indptr):
            col[indptr[index]:indptr[index + 1]] = index
        else:
            col[indptr[index]:] = index
    data = coo_matrix((data_raw, (row, col)), shape=tuple(shape)).toarray()

elif os.path.isdir(input_path): # 10x format
    mtx_10x = sc.read_10x_mtx(input_path)
    data = mtx_10x.X.toarray()

    def generate_cellname(celltype, cell):
        return celltype + '.' + str(cell).replace("_", ".")

    annotation_df = pd.read_csv(annotation_path, delimiter="\t", header=None)
    annotation_df['cellname'] = annotation_df.apply(lambda x: generate_cellname(x.iloc[0], x.iloc[1]), axis=1)

    cellname_list = list(annotation_df['cellname'])
    genename_list = list(mtx_10x.var_names)

# normalize
data *= 1E5 / data.sum()
data = np.log2(data + 1)
# centralize
row_means = np.mean(data, axis=1)
gem_data = data - row_means[:, np.newaxis]

gem_df = pd.DataFrame(data=gem_data, index=genename_list, columns=cellname_list)
MAX_CELL_PER_FILE = 100000
if len(cellname_list) <= MAX_CELL_PER_FILE:
    gem_path = os.path.join(output_file_directory, f'{output_tag}.csv')
    gem_df.to_csv(gem_path)
else:  # one celltype a file
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

print(f"{output_tag} convert end")
