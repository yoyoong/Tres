import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument('-CT', "--celltype", type=str, default='CD8', required=False, help="cell type")
args = parser.parse_args()

celltype = args.celltype
qc_result_dir = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.qc_data/qc_result'
qc_result_filename = os.path.join(qc_result_dir, f'{celltype}.qc_result.csv')
qc_result_data = pd.read_csv(qc_result_filename, sep="\t", header=0)

def get_group(x):
    if x.t >= 0 and x.p <= 0.05:
        return "t > 0, p < 0.05"
    elif x.t <= 0 and x.p <= 0.05:
        return "t < 0, p < 0.05"
    elif x.p > 0.05:
        return "p > 0.05"
qc_result_data['cut'] = qc_result_data.apply(lambda x: get_group(x), axis=1)

fig, axes = plt.subplots(2, 1, figsize=(40, 20))
# axes[0].set_title(celltype)

histplot_fig = plt.figure(figsize=(40, 10))
sns.histplot(data=qc_result_data, x="Cytokine", hue="cut",
    multiple="stack", color='blue', ax=axes[0])

axes[1].axhline(0, color='grey', linestyle='--', linewidth=2)
sns.violinplot(x="Cytokine", y="t", data=qc_result_data, ax=axes[1])
sns.stripplot(x="Cytokine", y="t", data=qc_result_data, color='b', size=1, ax=axes[1])

figure = fig.get_figure()
figure_filename = os.path.join(qc_result_dir, f'{celltype}.pdf')
figure.savefig(figure_filename, dpi=100)
figure.show()