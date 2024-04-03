import pandas as pd
import numpy as np
import os
from sklearn.preprocessing import StandardScaler


dataset_list = ["Zhang2021", "SadeFeldman2018", "Yost2019", "Fraietta2018"]
for dataset in dataset_list:
    output_tag = f'{dataset}.bulk_profile'
    output_file_directory = f'/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/3.clinical_data/{dataset}'
    correlation_filename = os.path.join(output_file_directory, f'{output_tag}.csv')

    gem_bulk = pd.read_csv(correlation_filename, index_col=0, header=0)

    scaler = StandardScaler()
    gem_bulk_normalized = pd.DataFrame(scaler.fit_transform(gem_bulk), columns=gem_bulk.columns)

    gem_bulk.to_csv(os.path.join(output_file_directory, f'{output_tag}.normalized.csv'))

