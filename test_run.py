import pandas as pd
import numpy as np

gem_df = pd.read_csv('/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/3.clinical_data/SadeFeldman2018/test.csv', header=[0, 1], sep='\t')
gem_df.columns = gem_df.columns.map('.'.join)# filter the gem

# gem_df *= 1E5 / gem_df.sum()
# gem_df = gem_df.add(1).apply(np.log2)
# gem_df = gem_df.sub(gem_df.mean(axis=1), axis=0)


flag_group = [col.split('.')[1] for col in gem_df.columns]
gem_bulk = gem_df.groupby(flag_group, axis=1).mean()
gem_bulk = gem_bulk.sub(gem_bulk.mean(axis=1), axis=0)
print(gem_bulk)
