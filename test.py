import numpy as np
import pandas as pd
import os

dataset_markers_dir = "/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_data/6.Marker/DatasetMarker"
celltype_markers_dir = "/sibcb1/bioinformatics/hongyuyang/dataset/Tres/2.tisch_data/6.Marker/CelltypeMarker"

all_markers = {}
celltype_filenum = {}
dataset_markers_files = sorted(os.listdir(dataset_markers_dir))
for dataset_markers_file in dataset_markers_files:
    dataset_markers_path = os.path.join(dataset_markers_dir, dataset_markers_file)
    dataset_markers = pd.read_csv(dataset_markers_path, index_col=0, header=0)
    dataset_celltype = dataset_markers.index.unique().tolist()
    for celltype in dataset_celltype: # 统计每个celltype的
        celltype_filenum[celltype] = celltype_filenum.get(celltype, 0) + 1

    # for celltype, row in dataset_markers.iterrows():
    #     if celltype not in all_markers.keys():
    #         new_dataset_markers = pd.DataFrame(columns=['gene', 'include_num'])
    #         new_row = {'gene': row['gene'], 'include_num': 1}
    #         new_dataset_markers.loc[len(new_dataset_markers)] = new_row
    #         all_markers[celltype] = new_dataset_markers
    #     else:
    #         old_dataset_markers = all_markers[celltype]
    #         if row['gene'] not in list(old_dataset_markers['gene']):
    #             new_row = {'gene': row['gene'], 'include_num': 1}
    #             old_dataset_markers.loc[len(old_dataset_markers)] = new_row
    #         else:
    #             index = old_dataset_markers[old_dataset_markers['gene'] == row['gene']].index.values[0]
    #             old_dataset_markers.loc[index, 'include_num'] += 1
    print(f'{dataset_markers_file} process end.')

for celltype, markers in dataset_markers.iterrows():
    markers['dataset_num'] = markers.index.map(lambda x: celltype_filenum[x])
    markers['rate'] = markers['include_num'] / markers['dataset_num']
    markers.to_csv(os.path.join(celltype_markers_dir, f'{celltype}.csv'))