import argparse
import os
import numpy as np
import pandas as pd
import h5py
import scanpy as sc
from scipy.sparse import coo_matrix

parser = argparse.ArgumentParser()
parser.add_argument('-I', "--input", type=str, required=False, help="input file path",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.raw_data/tisch2/ALL_GSE153697_expression.h5')
parser.add_argument('-R', "--celltype_mapping_rules_file", type=str, required=False, help="Celltype mapping rules file, .txt format",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.model_data/tisch2_celltype_mapping_rule.txt')
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.gem_data/tisch2')
parser.add_argument('-O', "--output_tag", type=str, required=False, help="Prefix for output files.",
                    default='ALL_GSE153697')
args = parser.parse_args()

input_path = args.input
celltype_mapping_rules_file = args.celltype_mapping_rules_file
output_file_directory = args.output_file_directory
output_tag = args.output_tag

genename_list = []
cellname_list = []
if input_path.endswith('.h5'):
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

    def generate_cellname(Cell, Celltype, Patient):
        if Celltype in mapping_rules_dict.keys():
            Celltype = mapping_rules_dict[Celltype]
        return Celltype + '.' + Patient + '.' + Cell

    CellMetainfo_df['cellname'] = CellMetainfo_df.apply(
        lambda x: generate_cellname(x['Cell'], x['Celltype (major-lineage)'], x['Patient']), axis=1)
    cellname_list = list(CellMetainfo_df['cellname'])

    # convert to dense matrix
    col = np.zeros_like(indices)
    for index in range(1, len(indptr)):
        col[indptr[index]:] = index
    row = indices
    data = coo_matrix((data_raw, (row, col)), shape=tuple(shape)).toarray()

    # centralize
    row_means = np.mean(data, axis=1)
    gem_data = data - row_means[:, np.newaxis]

gem_path = os.path.join(output_file_directory, f'{output_tag}.csv')
gem_df = pd.DataFrame(data=gem_data, index=genename_list, columns=cellname_list)
gem_df.to_csv(gem_path)
print(f"{output_tag} convert end")
