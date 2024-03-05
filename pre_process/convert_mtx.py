import argparse
import os
import numpy as np
import pandas as pd
import h5py
import scanpy as sc
from scipy.sparse import coo_matrix

parser = argparse.ArgumentParser()
parser.add_argument('-I', "--input", type=str, required=False, help="input file path",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/0.raw_data/NSCLC_GSE127471_expression.h5')
parser.add_argument('-CTR', "--celltype_mapping_rules_file", type=str, required=False, help="Celltype mapping rules file, .txt format",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.model_file/tisch2_celltype_mapping_rule.txt')
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/1.gem_data')
parser.add_argument('-O', "--output_tag", type=str, required=False, help="Prefix for output files.",
                    default='NSCLC_GSE127471')
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

    MAX_CELL_PER_FILE = 100
    if len(indptr) <= MAX_CELL_PER_FILE:
        # convert to dense matrix
        row = indices
        col = np.zeros_like(indices)
        for index in range(1, len(indptr)):
            col[indptr[index]:indptr[index + 1]] = index
        data = coo_matrix((data_raw, (row, col)), shape=tuple(shape)).toarray()

        # centralize
        row_means = np.mean(data, axis=1)
        gem_data = data - row_means[:, np.newaxis]

        gem_path = os.path.join(output_file_directory, f'{output_tag}.csv')
        gem_df = pd.DataFrame(data=gem_data, index=genename_list, columns=cellname_list)
        gem_df.to_csv(gem_path)
    else: # one file per 100,000 lines
        gem_dir = os.path.join(output_file_directory, output_tag)
        if not os.path.exists(gem_dir):
            os.makedirs(gem_dir)
        file_num = len(indptr) // MAX_CELL_PER_FILE if len(indptr) % MAX_CELL_PER_FILE == 0 else (len(indptr) // MAX_CELL_PER_FILE + 1)
        for i in range(file_num):
            start_index = i * MAX_CELL_PER_FILE
            end_index = min(len(indptr) - 1, (i + 1) * MAX_CELL_PER_FILE)

            # convert to dense matrix
            row = indices[indptr[start_index]:indptr[end_index]]
            col = np.zeros(indptr[end_index] - indptr[start_index])
            for index in range(start_index + 1, end_index):
                col[indptr[index]:indptr[index + 1]] = index
                # if index < end_index:
                #     col[indptr[index]:indptr[index + 1]] = index
                # else:
                #     col[indptr[index]:] = index
            data_raw_split = data_raw[indptr[start_index]:indptr[end_index]]
            data = coo_matrix((data_raw_split, (row, col)), shape=(shape[0], end_index - start_index)).toarray()

            # centralize
            row_means = np.mean(data, axis=1)
            gem_data = data - row_means[:, np.newaxis]

            gem_path = os.path.join(gem_dir, f'{output_tag}_{str(i + 1)}.csv')
            cellname_list_split = cellname_list[start_index:end_index]
            gem_df = pd.DataFrame(data=gem_data, index=genename_list, columns=cellname_list_split)
            gem_df.to_csv(gem_path)
            print(f"{output_tag}_{str(i + 1)} convert end.")

print(f"{output_tag} convert end")
