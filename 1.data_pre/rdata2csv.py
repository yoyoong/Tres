import os
import numpy as np
import pandas as pd
import pyreadr
import rpy2.robjects as robjects

rdata_path = "/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/4.Macrophage_analysis/TCGA/raw_RData"
output_path = "/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/4.Macrophage_analysis/TCGA/gem"
annatation_path = "/sibcb2/bioinformatics/KnowledgeBase/Firehose_Methylation/RnBeads_450K_hg19_Probes_GR.RData"
# annatation_rdata = pyreadr.read_r(annatation_path)
# for key in annatation_rdata.keys():
#     annatation_df = annatation_rdata[key]
#     print(annatation_df.shape)

dddd = robjects.r['load'](annatation_path)

rdata_list = os.listdir(rdata_path)
for rdata_filename in rdata_list:
    rdata = pyreadr.read_r(os.path.join(rdata_path, rdata_filename))
    for key in rdata.keys():
        df = rdata[key]
        df.to_csv(os.path.join(output_path, f"{key}.csv"), index=True)
