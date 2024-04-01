import argparse
import os
import numpy as np
import pandas as pd
import h5py
import scanpy as sc
from scipy.sparse import coo_matrix
from tqdm.autonotebook import tqdm
from scipy.stats import pearsonr

dataset_list = ["Zhang2021", "SadeFeldman2018", "Yost2019", "Fraietta2018"]
for dataset in dataset_list:
    tres_signature_path = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/4.Interaction/Tres_signature.negative.csv'
    sample_annotation_path = f'/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/3.clinical_data/{dataset}/{dataset}.sample_annotation.txt'
    output_file_directory = f'/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/3.clinical_data/{dataset}'
    output_tag = f'{dataset}.correlation'

    # bulk_profile_path = f'/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/3.clinical_data/{dataset}/{dataset}.bulk_profile.log_centralize.csv'
    if dataset == 'Zhang2021':
        bulk_profile_path = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/test/data/Atezolizumab+Paclitaxel_Pre_Zhang2021_TNBC.csv'
    elif dataset == 'SadeFeldman2018':
        bulk_profile_path = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/test/data/ICB_Pre_SadeFeldman2018_Melanoma.csv'
    elif dataset == 'Yost2019':
        bulk_profile_path = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/test/data/anti-PD1_Pre_Yost2019_BCC.csv'
    elif dataset == 'Fraietta2018':
        bulk_profile_path = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/test/data/CD19CAR_Infusion.CAR_Fraietta2018_CLL.csv'

    tres_signature_df = pd.read_csv(tres_signature_path, index_col=0)
    # filter top and down 5000 genes
    top_rows = tres_signature_df[tres_signature_df['Tres'].isin(tres_signature_df['Tres'].nlargest(2000))]
    down_rows = tres_signature_df[tres_signature_df['Tres'].isin(tres_signature_df['Tres'].nsmallest(2000))]
    tres_signature_df = pd.concat([top_rows, down_rows], axis=0)

    bulk_profile_df = pd.read_csv(bulk_profile_path, index_col=0)
    gene_list = bulk_profile_df.index.intersection(tres_signature_df.index) # common gene list

    tres_signature_filtered = tres_signature_df.loc[gene_list]
    bulk_profile_filtered = bulk_profile_df.loc[gene_list]
    tres_signature = tres_signature_filtered['Tres']

    # sample_annotation_df = pd.read_csv(sample_annotation_path, delimiter='\t', index_col=0)
    # sample_list = sample_annotation_df.index.values.tolist()

    correlation_df = pd.DataFrame(columns=['sample', 'correlation'])
    for sample in bulk_profile_df.columns.values:
        bulk_profile_data = bulk_profile_filtered[sample]
        correlation, _ = pearsonr(np.array(bulk_profile_data), np.array(tres_signature))

        new_row = pd.Series({'sample': sample, 'correlation': correlation})
        correlation_df = pd.concat([correlation_df, new_row.to_frame().T])

    correlation_df.set_index('sample', inplace=True)
    correlation_filename = os.path.join(output_file_directory, f'{output_tag}.csv')
    correlation_df.to_csv(correlation_filename)
    print(f"{dataset} process end")
