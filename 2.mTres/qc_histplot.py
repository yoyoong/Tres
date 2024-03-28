import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser()
parser.add_argument('-CT', "--celltype", type=str, default='CD8T', required=False, help="cell type")
parser.add_argument('-R', "--qc_result_file", type=str, required=False, help="qc result file path",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/3.qc_result/CD8T.qc_result.csv')
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/3.qc_result')
args = parser.parse_args()

celltype = args.celltype
qc_result_path = args.qc_result_file
output_file_directory = args.output_file_directory
qc_result_data = pd.read_csv(qc_result_path, sep="\t", header=0)

def get_group(x):
    if x.t >= 0 and x.p <= 0.05:
        return "t > 0, p < 0.05"
    elif x.t <= 0 and x.p <= 0.05:
        return "t < 0, p < 0.05"
    elif x.p > 0.05:
        return "p > 0.05"
qc_result_data['cut'] = qc_result_data.apply(lambda x: get_group(x), axis=1)

fig = plt.figure(figsize=(50, 20))
histplot = sns.histplot(data=qc_result_data, x="Cytokine", hue="cut", multiple="stack", color='blue')
histplot.xaxis.label.set_size(30)
histplot.yaxis.label.set_size(30)
histplot.set_yticklabels(histplot.get_yticks(), size = 15)

# fig, axes = plt.subplots(2, 1, figsize=(50, 20))

# sns.histplot(data=qc_result_data, x="Cytokine", hue="cut", multiple="stack", color='blue', ax=axes[0])

# axes[1].axhline(0, color='grey', linestyle='--', linewidth=2)
# sns.violinplot(x="Cytokine", y="t", data=qc_result_data, ax=axes[1])
# sns.stripplot(x="Cytokine", y="t", data=qc_result_data, color='b', size=1, ax=axes[1])

figure = fig.get_figure()
figure_filename = os.path.join(output_file_directory, f'{celltype}.hisplot.pdf')
figure.savefig(figure_filename, dpi=100)
figure.show()