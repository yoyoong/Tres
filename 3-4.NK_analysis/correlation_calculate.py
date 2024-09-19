import argparse
import os
import numpy as np
import pandas as pd
import h5py
import scanpy as sc
from scipy.sparse import coo_matrix
from tqdm.autonotebook import tqdm
from scipy.stats import pearsonr
from scipy.stats import normaltest,shapiro

celltype_list = ['CD8T', 'Macrophage', 'Neutrophils', 'NK', 'NK_act']
dataset_list = ["TCGA-BRCA", "TCGA-CESC", "TCGA-SKCM"]

for celltype in celltype_list:
    if celltype == "CD8T":
        tres_signature_dir = '/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_data/5-1.CD8T_Interaction/Tres_signature'
    elif celltype == "Macrophage":
        tres_signature_dir = '/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_data/5-2.Macrophage_Interaction/Tres_signature'
    elif celltype == "Neutrophils":
        tres_signature_dir = '/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_data/5-3.Neutrophils_Interaction/Tres_signature'
    elif celltype == "NK":
        tres_signature_dir = '/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_data/5-4.NK_Interaction/Tres_signature'
    elif celltype == "NK_act":
        tres_signature_dir = '/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_data/5-4-0.NK_act_Interaction/Tres_signature'

    for dataset in dataset_list:
        bulk_profile_path = f'/sibcb1/bioinformatics/hongyuyang/dataset/Tres/3-4.NK_analysis/{dataset}/gem.csv'
        cytokine_info_version = 0
        cytokine_signature_version = 1
        sample_filter_version = 1
        version = f'{cytokine_info_version}_{cytokine_signature_version}_{sample_filter_version}'
        tres_signature_positive_path = os.path.join(tres_signature_dir, f'Tres_signature_{version}.positive.csv')
        tres_signature_negative_path = os.path.join(tres_signature_dir, f'Tres_signature_{version}.negative.csv')
        tres_signature_positive = pd.read_csv(tres_signature_positive_path, index_col=0)
        tres_signature_negative = pd.read_csv(tres_signature_negative_path, index_col=0)

        def calculate_correlation(bulk_profile_df, tres_signature_df):
            gene_list = bulk_profile_df.index.intersection(tres_signature_df.index)  # common gene list

            tres_signature_filtered = tres_signature_df.loc[gene_list]
            bulk_profile_filtered = bulk_profile_df.loc[gene_list]
            tres_signature = tres_signature_filtered['Tres']

            correlation_df = pd.DataFrame(columns=['sample', 'correlation'])
            for sample in bulk_profile_df.columns.values:
                bulk_profile_data = bulk_profile_filtered[sample]
                correlation, _ = pearsonr(np.array(bulk_profile_data), np.array(tres_signature))

                new_row = pd.Series({'sample': sample, 'correlation': correlation})
                correlation_df = pd.concat([correlation_df, new_row.to_frame().T])

            correlation_df.set_index('sample', inplace=True)
            return correlation_df

        bulk_profile_df = pd.read_csv(bulk_profile_path, index_col=0)
        correlation_positive_df = calculate_correlation(bulk_profile_df, tres_signature_positive)
        correlation_negative_df = calculate_correlation(bulk_profile_df, tres_signature_negative)

        output_file_directory = f'/sibcb1/bioinformatics/hongyuyang/dataset/Tres/3-4.NK_analysis/{dataset}'
        output_tag = f'{celltype}.correlation'
        correlation_positive_filename = os.path.join(output_file_directory, f'{output_tag}.positive.csv')
        os.remove(correlation_positive_filename)
        correlation_negative_filename = os.path.join(output_file_directory, f'{output_tag}.negative.csv')
        os.remove(correlation_negative_filename)

        output_file_directory = f'/sibcb1/bioinformatics/hongyuyang/dataset/Tres/3-4.NK_analysis/{dataset}/correlation'
        if not os.path.exists(output_file_directory):
            os.makedirs(output_file_directory)
        output_tag = f'{celltype}.correlation'
        correlation_positive_filename = os.path.join(output_file_directory, f'{output_tag}.positive.csv')
        correlation_positive_df.to_csv(correlation_positive_filename)
        correlation_negative_filename = os.path.join(output_file_directory, f'{output_tag}.negative.csv')
        correlation_negative_df.to_csv(correlation_negative_filename)
    print(f"{celltype} process end")
