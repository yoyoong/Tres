import gseapy as gp
import pandas as pd

gem_path = f'/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/5.Analysis_data/Braun2020/Braun2020.Express.tsv'
gem = pd.read_csv(gem_path, index_col=0, header=0, delimiter='\t')
gem.set_index('GENE_SYMBOL', inplace=True)
gene_sets = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.model_file/Tres_kegg.Neutrophils.gmt'

sample_annotation_path = f'/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/5.Analysis_data/Braun2020/Braun2020.Response.tsv'
sample_annotation_df = pd.read_csv(sample_annotation_path, delimiter='\t', header=0, index_col=0)
sample_annotation_df.set_index('sample_id', inplace=True)

non_response_flag = sample_annotation_df['response_NR'] == 'N'
response_sample = sample_annotation_df[~non_response_flag].index
nonresponse_sample = sample_annotation_df[non_response_flag].index
cls = ['ALL'] * len(response_sample) + ['AML'] * len(nonresponse_sample)

ss_res = gp.ssgsea(data=gem,
               gene_sets=gene_sets,
               outdir=None,
               sample_norm_method='rank', # choose 'custom' will only use the raw value of `data`
               no_plot=True)
nes = ss_res.res2d.pivot(index='Term', columns='Name', values='NES')
neut_nes = nes.loc['Neut_IFN-15']
print(neut_nes)