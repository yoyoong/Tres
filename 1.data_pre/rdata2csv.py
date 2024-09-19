import os
import numpy as np
import pandas as pd
import pyreadr
import rpy2.robjects as robjects

parser = argparse.ArgumentParser()
parser.add_argument('-CT', "--celltype", type=str, default='CD8T', required=False, help="cell type")
parser.add_argument('-CIV', "--cytokine_info_version", type=int, default='2', required=False, help="cytokine info version")
parser.add_argument('-CSV', "--cytokine_signature_version", type=int, default='1', required=False, help="cytokine signature version")
parser.add_argument('-SFV', "--sample_filter_version", type=int, default='1', required=False, help="sample filter version")
args = parser.parse_args()


rdata_path = "/sibcb1/bioinformatics/hongyuyang/dataset/Tres/3-2.Macrophage_analysis/TCGA/raw_RData"
output_path = "/sibcb1/bioinformatics/hongyuyang/dataset/Tres/3-2.Macrophage_analysis/TCGA/gem"
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
